#include <numeric>      // std::accumulate
//#include <inttypes.h>

#include "TFile.h"
#include "TMath.h"
#include "TDatime.h"
#include "TH2.h"
#include "TStyle.h"
#include "TString.h"
#include "TVirtualFitter.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "THStack.h"
#include "TVirtualPad.h"
#include "TLatex.h"
#include "TVector2.h"
#include "TEllipse.h"
#include "TMinuit.h"




//#include "/ds50/app/user/xxiang89/work/xy4slad_ms/XYReconstructor/PMTGeom.h"
//#include "/ds50/app/user/xxiang89/work/xy4slad_ms/XYReconstructor/XYReconstructor.h"

#include "PMTGeom.h"
#include "XYReconstructor.h"


using namespace std;

Double_t XYReconstructor::S2overS1(Double_t *x, Double_t *par) {
  Double_t offset = par[0]; //offset
  Double_t mag = par[1]; //magnitude
  Double_t scale = par[2]; //magnitude

//  Double_t r = scale*x[0];
  Double_t r = x[0]*scale;
//  Double_t r = x[0];

  if(r>PMTGeom::tpcRadi) r = PMTGeom::tpcRadi;
//  if(r>PMTGeom::tpcRadi) r = scale*PMTGeom::tpcRadi;
//  return par[0]+ par[1]*TMath::Power(scale*PMTGeom::tpcRadi*scale*PMTGeom::tpcRadi-r*r, 2);
  return offset+ mag*TMath::Power(PMTGeom::tpcRadi*PMTGeom::tpcRadi-r*r, 2);
}

Bool_t XYReconstructor::CheckTH2Bounds(TH2* h, Double_t x, Double_t y) const {
  Int_t bin_x = h->GetXaxis()->FindBin(x);
  Int_t bin_y = h->GetYaxis()->FindBin(y);
  if(bin_x<1 || bin_x>h->GetNbinsX() || bin_y<1 || bin_y>h->GetNbinsY()) {
     return false;
  } else return true;
}

//#define S2OVERS1PROFILE //switch to use S2/S1 vs R profile
//#define S2BOTOVERS2PROFILE //switch to use S2_bottom/S2 vs R profile
#define S2BOTOVERS2TOPPROFILE //switch to use S2_bottom/S2_top vs xy profile
#define PROFILE2D //switch to use S2_sh/S2 vs XY profile, not 1D but 2D profile
//#define DENSITY //switch to use density info
#define LOGS2 //switch to use Log(S2ch/S2tot)
#define LOGS2ERR //switch to use error of Log(S2ch/S2tot)
//#define ROTATIONSYM //use 2d profile from pmt categories
//#define USEALEMETHOD //use Alessandro's 1d profile. This need to use with ROTATIONSYM


Int_t XYReconstructor::SetMinimizer(){
  // Initialize TMinuit via Minimizer
  // create minimizer giving a name and a name (optionally) for the specific
  // algorithm
  // possible choices are:
  //     minName                  algoName
  // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
  //  Minuit2                     Fumili2
  //  Fumili
  //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS,
  //                              BFGS2, SteepestDescent
  //  GSLMultiFit
  //   GSLSimAn
  //   Genetic
  const std::string minimizer = "Minuit2";// "GSLMultiMin"
  const std::string algoType  = ""; // "ConjugateFR", "ConjugatePR", "BFGS", "BFGS2"
//  const std::string minimizer = "GSLMultiMin";
//  const std::string algoType  = "ConjugateFR"; // "ConjugateFR", "ConjugatePR", "BFGS", "BFGS2"
  gErrorIgnoreLevel=kWarning; // suppress message

  fMinimizer = ROOT::Math::Factory::CreateMinimizer(minimizer, algoType);
  std::cout << "Using minimizer " << minimizer << "  " << algoType << std::endl;


  if (fMinimizer == 0) {
     std::cerr << "Error using minimizer " << minimizer << "  " << algoType << std::endl;
     return -1;
  }
//  if (Debug) fMinimizer->SetPrintLevel(3);
  fMinimizer->SetMaxFunctionCalls(5000);
  fMinimizer->SetMaxIterations(5000);
//  fMinimizer->SetTolerance(0.001);
//  fMinimizer->SetTolerance(0.1);
//  fMinimizer->SetTolerance(10);
  fMinimizer->SetTolerance(0.1);
//  fMinimizer->SetTolerance(1e10);
  cout<<"Minimizer Precision is "<<fMinimizer->Precision()<<endl;
  
  //to set n-sigma, set n^2
//  cout<<"Minimizer Error def: sigma= "<<2<<endl;
//  fMinimizer->SetErrorDef(4);

  
  return kTRUE;
}

Int_t XYReconstructor::Minimization(SignalProfile &sigProf){
//  Info("XYReconstructor::Minimization()","Minimization ...");
  ndf_inchi2 = 0; //initialize every time
  this->SetSignalProfile(&sigProf);
  fMinimizer->SetFunction(*this);

  Double_t initx(PMTGeom::pmtUnit*PMTGeom::pmt_x[sigProf.max_s2_chan]), inity(PMTGeom::pmtUnit*PMTGeom::pmt_y[sigProf.max_s2_chan]);
//  Double_t initx(sigProf.PosiX), inity(sigProf.PosiY);
//  cout<<"initx: "<<initx<<" inity: "<<inity<<endl;
//  cout<<"    sigProf.PosiX: "<<sigProf.PosiX<<" sigProf.PosiY: "<<sigProf.PosiY<<endl;
  Minimize2DFunc(initx, inity);
//  cout<<"initx: "<<initx<<" inity: "<<inity<<endl;

//  fMinimizer->SetVariable(0, "X posi.", initx, 0.05);//PMTGeom::pmtUnit);//0.8); //Int_t SetParameter(Int_t ipar, const char* parname, Double_t value, Double_t verr, Double_t vlow, Double_t vhigh)
//  fMinimizer->SetVariable(1, "Y posi.", inity, 0.05);//PMTGeom::pmtUnit);//0.8);
//  fMinimizer->SetLimitedVariable(0, "X posi.", initx, 8.e-1, -PMTGeom::tpcRadi, PMTGeom::tpcRadi); //virtual bool SetLimitedVariable(unsigned int /* ivar */ , const std::string & /* name */ , double /*val */ , double /* step  */, double /* lower */, double /* upper */)
//  fMinimizer->SetLimitedVariable(1, "Y posi.", inity, 8.e-1, -PMTGeom::tpcRadi, PMTGeom::tpcRadi);

  fMinimizer->SetVariable(0, "X posi.", initx, 2);
  fMinimizer->SetVariable(1, "Y posi.", inity, 2);


  //set minimizer per event
//  fMinimizer->SetTolerance(0.1);
  bool ret= fMinimizer->Minimize();
//  minval = fMinimizer->MinValue();
//  edm = fMinimizer->Edm();

  if (!ret) {// Minimization Fail
      if (fMinimizer->Status() ==3) {
          fMinimizer->SetTolerance(10);
          fMinimizer->SetVariableValue(0, initx);
          fMinimizer->SetVariableValue(1, inity);
          ret= fMinimizer->Minimize();
          if (!ret){
              if (fMinimizer->Status() ==3) {
                fMinimizer->SetTolerance(1000);
                fMinimizer->SetVariableValue(0, initx);
                fMinimizer->SetVariableValue(1, inity);
                ret= fMinimizer->Minimize();
                if (!ret)return -1;
              }
          }
      }
      else
        return -1;
  }


//  //set minimizer per event
//  fMinimizer->SetTolerance(100000);
//  bool ret= fMinimizer->Minimize();
//  fMinimizer->SetTolerance(1000);
//  fMinimizer->SetVariableStepSize(0, 0.5)
//  fMinimizer->SetVariableStepSize(1, 0.5)
//  bool ret= fMinimizer->Minimize();
//  fMinimizer->SetTolerance(10);
//  fMinimizer->SetVariableStepSize(0, 0.05)
//  fMinimizer->SetVariableStepSize(1, 0.05)
//  bool ret= fMinimizer->Minimize();


// Find the max chi2 as well by brutal grid search
  Double_t x1 = 0; Double_t x2=0;
  Double_t maxchi2 = FindMaxChi2(*gf2D, x1, x2);
  sigProf.maxchi2_x = x1;
  sigProf.maxchi2_y = x2;
  ndf_inchi2 = (ndf_inchi2==0)? CountDF(): ndf_inchi2;
  Int_t nparam = fMinimizer->NFree();
  sigProf.maxchi2 = maxchi2/(ndf_inchi2-nparam);


//  // test Minos (use the default up of 1)
//  if (minos) {
//
//     double el,eu;
//     for (unsigned int i = 0; i < func.NDim(); ++i) {
//        ret  = fMinimizer->GetMinosError(i,el,eu);
//        if (ret) std::cout << "MINOS error  for " << i  << " = " << el << "   " << eu << std::endl;
//        else     std::cout << "MINOS failed for " << i << std::endl;
//     }
//  }

  //get result
  const double * xmin = fMinimizer->X();
  sigProf.PosiX = xmin[0];
  sigProf.PosiY = xmin[1];

  const double * xerror = fMinimizer->Errors();
  sigProf.PosiXErr = xerror[0];
  sigProf.PosiYErr = xerror[1];

  sigProf.chi2XY = ReducedChi2();//fMinimizer->MinValue();

  return 1;
}

Int_t XYReconstructor::MinimizationStatus(){
// status = 1    : Covariance was made pos defined
// status = 2    : Hesse is invalid
// status = 3    : Edm is above max
// status = 4    : Reached call limit
// status = 5    : Any other failure

  return fMinimizer->Status();
}

Double_t XYReconstructor::ReducedChi2(){
	Double_t chi2 = fMinimizer->MinValue();
	Int_t nparam = fMinimizer->NFree();
	ndf_inchi2 = (ndf_inchi2==0)? CountDF(): ndf_inchi2;
	Double_t reducedchi2 = chi2/(ndf_inchi2-nparam);
//  cout<<"chi2="<<chi2<<" ndf_inchi2="<<ndf_inchi2<<" nparam="<<nparam<<"\n";
	return reducedchi2;
}

void XYReconstructor::Minimize2DFunc(Double_t& x1, Double_t& x2){
#if 0
  Double_t delta_x1 = 10.;//10.; // arrow move 18cm
  Double_t min_x1 = x1-delta_x1;
  Double_t max_x1 = x1+delta_x1;
  Double_t delta_x2 = 10.;//10.;
  Double_t min_x2 = x2-delta_x2;
  Double_t max_x2 = x2+delta_x2;

//  TF2 f2D("f2D", this, min_x1, max_x1, min_x2, max_x2, 0, "XYReconstructor");//TF2(const char* name, void* ptr, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Int_t npar, const char* className)
//  cout<<"min_x1: "<<min_x1<<" max_x1: "<<max_x1<<endl;
//  cout<<"f2D.GetXmin(): "<<f2D.GetXmin()<<" f2D.GetXmax(): "<<f2D.GetXmax()<<endl;
//  if(gf2D) delete gf2D;
//  gf2D = new TF2("f2D", this, min_x1, max_x1, min_x2, max_x2, 0, "XYReconstructor");//TF2(const char* name, void* ptr, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Int_t npar, const char* className)
  gf2D->SetRange(min_x1, min_x2, max_x1, max_x2);
  gf2D->GetMinimumXY(x1, x2);
  gMinuit->Command("SET NOW");
  gMinuit->SetPrintLevel(-2);
#else
  Scan2DFunc(*gf2D, x1, x2);
//  Scan2DFunc(f2D, x1, x2);
//  Double_t r(TMath::Sqrt(x1*x1+x2*x2)), theta(TMath::ATan2(x2, x1));
//  cout<<"X: "<<x1<<" Y: "<<x2<<endl; // convert xy
//  cout<<"R: "<<r<<" Theta: "<<theta<<endl;
#endif
}

void XYReconstructor::Scan2DFunc(TF2 &f2D, Double_t& x1, Double_t& x2){

/*  Double_t fXmin = f2D.GetXmin();
  Double_t fXmax = f2D.GetXmax();
  Double_t fYmin = f2D.GetYmin();
  Double_t fYmax = f2D.GetYmax();
  Double_t fNpx = f2D.GetNpx();
  Double_t fNpy = f2D.GetNpy();*/
//  Double_t delt1(1.5*PMTGeom::pmtUnit), delt2(1.5*PMTGeom::pmtUnit);
  Double_t delt1(2.5*PMTGeom::pmtUnit), delt2(2.5*PMTGeom::pmtUnit);
//  Double_t delt1(10.), delt2(10.);
//  cout<<"x1: "<<x1<<" x2: "<<x2<<" delt1: "<<delt1<< endl;
//  cout<<"x1-delt1: "<<x1-delt1<<" f2D.GetXmin(): "<<f2D.GetXmin()<<endl;
  Double_t fXmin = TMath::Max(x1-delt1,f2D.GetXmin());
  Double_t fXmax = TMath::Min(x1+delt1,f2D.GetXmax());
  Double_t fYmin = TMath::Max(x2-delt2,f2D.GetYmin());
  Double_t fYmax = TMath::Min(x2+delt2,f2D.GetYmax());
  Double_t fNpx = 50.;//f2D.GetNpx();
  Double_t fNpy = 50.;//f2D.GetNpy();
//  cout<<"Xmin: "<<fXmin<<" Xmax: "<<fXmax<<" Ymin: "<<fYmin<<" Ymax: "<<fYmax<<endl;
  //First do a grid search with step size fNpx and fNpy
     Double_t xx, yy, zz;
     Double_t dx = (fXmax - fXmin)/fNpx;
     Double_t dy = (fYmax - fYmin)/fNpy;
     Double_t xxmin = fXmin;
     Double_t yymin = fYmin;
     Double_t zzmin = f2D.Eval(xxmin, yymin+dy);
     for (Int_t i=0; i<fNpx; i++){
        xx=fXmin + (i+0.5)*dx;
        for (Int_t j=0; j<fNpy; j++){
           yy=fYmin+(j+0.5)*dy;
           zz = f2D.Eval(xx, yy);
           if (zz<zzmin) {xxmin = xx, yymin = yy; zzmin = zz;}
        }
     }

     x1 = TMath::Min(fXmax, xxmin);
     x2 = TMath::Min(fYmax, yymin);
}


//Xin - to find the max chi2 by brutal grid search
Double_t XYReconstructor::FindMaxChi2(TF2 &f2D, Double_t& x1, Double_t& x2){
  
  Double_t radi = PMTGeom::tpcRadi;
  
  Double_t fXmin = -radi; Double_t fXmax = radi;
  Double_t fYmin = -radi; Double_t fYmax = radi;
  
  Int_t fNpx = 36; Int_t fNpy = 36;
  
  Double_t x, y, z;
  Double_t dx = (fXmax - fXmin)/fNpx;
  Double_t dy = (fYmax - fYmin)/fNpy;
  Double_t xmax = 0;
  Double_t ymax = 0;
  Double_t zmax = 0;
  for (Int_t i=0; i<fNpx; i++){
    x=fXmin + (i+0.5)*dx;
    for (Int_t j=0; j<fNpy; j++){
      y=fYmin+(j+0.5)*dy;
      if (x*x+y*y<=radi*radi) {
        z = f2D.Eval(x, y);
        if (z>zmax) {xmax = x, ymax = y; zmax = z;}
      }
    }
  }
  
  x1 = xmax;
  x2 = ymax;
  return zmax;
}


void XYReconstructor::PrintComponents (const double * x) const {
  Double_t par[]={x[0], x[1]};
  cout<<endl;
  for(Int_t ich=0; ich<PMTGeom::N_CHANNELS; ich++){
      if (ich < PMTGeom::N_CHANNELS/2) continue; // use only top PMTs
//      if (ich >= PMTGeom::N_CHANNELS/2) continue; // use only bottom PMTs
//      if(PMTGeom::IsOuterCornerPmt(ich)) continue;

      Double_t expS2overS2total = XYReconstructor::GetExpectedS2overS2total(par, ich); // Should we make histogram of solid angle for lookup?
      Double_t expS2overS2totalErr = XYReconstructor::GetExpectedS2overS2totalErr(par, ich); // Should we make histogram of solid angle for lookup?
      // solid angle is used only in the first iteration. After 1st, response function should be used.
      Double_t s2 =fSigProf->s2[ich];// -channel.pulses[1].param.fixed_int2 / channel.pmt.spe_mean;
      Double_t chi = (s2/fSigProf->total_s2 - expS2overS2total)/expS2overS2totalErr;
      cout<<Form("ch%2d  s2/s2tot: %4.3f  exps2/s2tot: %4.3f  exps2/s2toterr: %4.3f  chi: %4.3f", ich, s2/fSigProf->total_s2, expS2overS2total, expS2overS2totalErr, chi)<<endl;
  }
  Double_t expS2overS1total = XYReconstructor::GetExpectedS2overS1total(par);
  Double_t expS2overS1totalErr = XYReconstructor::GetExpectedS2overS1totalErr(par);
  Double_t chi = (fSigProf->s2_over_s1_corr - expS2overS1total)/expS2overS1totalErr;
  cout<<Form("sum   s2/s1_corr: %4.3f  exps2/s1_corr: %4.3f  exps2/s1_correrr: %4.3f  chi: %4.3f", fSigProf->s2_over_s1_corr, expS2overS1total, expS2overS1totalErr, chi)<<endl;
}

Bool_t XYReconstructor::IsChanIncludedInChi2 (Int_t ich) const {
  if (ich < PMTGeom::N_CHANNELS/2) return false; // use only top PMTs
//      if (ich >= PMTGeom::N_CHANNELS/2) return false; // use only bottom PMTs
//      if(PMTGeom::IsOuterCornerPmt(ich)) return false;
//      else if(PMTGeom::IsOuterSidePmt(ich)) return false;
//  if(fSigProf->s2[ich]/fSigProf->total_s2< 1./(1.5*PMTGeom::N_CHANNELS)) return false;
//  if(TMath::Log10(fSigProf->s2[ich]/fSigProf->total_s2)< -2.) return false;
//  if(ich==fSigProf->max_s2_chan) return false;
  if(fSigProf->s2[ich] <= 2.e-3/*0.*/) return false; // if s2_ch is neggative, skip the channel since we use log10(s2_ch)
  if(fSigProf->s2[ich]/fSigProf->total_s2 <= 2.e-3/*0.*/) return false; // if s2_ch is neggative, skip the channel since we use log10(s2_ch)

#if 0 // nearby PMT only
  Double_t r2 = (PMTGeom::pmt_x[ich]-PMTGeom::pmt_x[fSigProf->max_s2_chan])*(PMTGeom::pmt_x[ich]-PMTGeom::pmt_x[fSigProf->max_s2_chan]) + (PMTGeom::pmt_y[ich]-PMTGeom::pmt_y[fSigProf->max_s2_chan])*(PMTGeom::pmt_y[ich]-PMTGeom::pmt_y[fSigProf->max_s2_chan]);
  if(r2>100.) return false;
#endif

  return true;
}

Int_t XYReconstructor::CountDF() const {
  Int_t ndf(0);
  for(Int_t ich=0; ich<PMTGeom::N_CHANNELS; ich++){
	 if (!IsChanIncludedInChi2 (ich)) continue;
	 ndf++;
  }
#ifdef S2OVERS1PROFILE //add S2/S1 Likelihood
  ndf++;
#endif

#ifdef S2BOTOVERS2PROFILE //add S2_bottom/S2 Likelihood
  ndf++;
#endif

#ifdef S2BOTOVERS2TOPPROFILE //add S2_bottom/S2_top Likelihood
  ndf++;
#endif

#ifdef DENSITY
  ndf++;
#endif
  return ndf;
}

Double_t XYReconstructor::DoEval (const double * x) const {
  Double_t par[]={x[0], x[1]};
  if(par[0]!=par[0] || par[1]!=par[1]) return 9999.;

#ifndef USEALEMETHOD
  Double_t r = TMath::Sqrt(par[0]*par[0] + par[1]*par[1]);
  if(r>PMTGeom::tpcRadi/*+2.*/) return 99999.;
#endif

  Double_t sum2(0.), sumweight(0.);
  for(Int_t ich=0; ich<PMTGeom::N_CHANNELS; ich++){
      if (!IsChanIncludedInChi2 (ich)) continue;

      Double_t chi2(0.), weight(0.);
      GetChi2Channel(par, ich, chi2, weight);
      if(chi2<0.) continue;
      sum2 += chi2;
      sumweight += weight;
  }
  Double_t chi2(0.), weight(0.);
#ifdef S2OVERS1PROFILE //add S2/S1 Likelihood
  GetChi2S2overS1(par, chi2, weight);
  sum2 += chi2;
  sumweight += weight;
#endif

#ifdef S2BOTOVERS2PROFILE //add S2_bottom/S2 Likelihood
  GetChi2S2BotoverS2(par, chi2, weight);
  sum2 += chi2;
//  cout<<"chi2: "<<chi2<<endl;
  sumweight += weight;
#endif

#ifdef S2BOTOVERS2TOPPROFILE //add S2_bottom/S2_top Likelihood
  GetChi2S2BotoverS2Top(par, chi2, weight);
  sum2 += chi2;
//  cout<<"chi2: "<<chi2<<endl;
  sumweight += weight;
#endif

#ifdef DENSITY
  GetDensityContribution(par, chi2, weight);
  sum2 += chi2;
  sumweight += weight;
//cout<<"chi2: "<<chi2<<endl;
#endif


//  sum2/=sumweight*sumweight;
  return sum2;
}


void XYReconstructor::GetDensityContribution(const double *par, Double_t &chi2, Double_t &weight) const {
  //par are supposed to be x and y
    Double_t density = XYReconstructor::GetDensity(par); // Should we make histogram of solid angle for lookup?
//    weight = 1.5e+4;
//    weight = 1./2.e-4;//1./5.e-4;
//    weight = 2./5.e-4;
    weight = 10.;//20.;//1.7e+3;
//    weight *= (1.01-TMath::Sqrt(par[0]*par[0]+par[1]*par[1])/PMTGeom::tpcRadi);
//    weight = 1.;
    chi2 = weight*density;
}

void XYReconstructor::GetChi2Channel(const double *par, Int_t ich, Double_t &chi2, Double_t &weight) const {
  //par are supposed to be x and y
#ifndef  USEALEMETHOD
//#if 1
  Double_t expS2overS2total = XYReconstructor::GetExpectedS2overS2total(par, ich); // Should we make histogram of solid angle for lookup?
//  Int_t ctg = PMTGeom::GetPmtPositionCategory(ich);
//  Double_t expS2overS2total = (ctg != CENTER && ctg != MIDDLE)? XYReconstructor::GetExpectedS2overS2total(par, ich):
//      Double_t expS2overS2total = (ctg != MIDDLE)? XYReconstructor::GetExpectedS2overS2total(par, ich):
//  Double_t expS2overS2total = (ich != 30)? XYReconstructor::GetExpectedS2overS2total(par, ich):
//  Double_t expS2overS2total = (ich != 30 || fSigProf->max_s1_chan!=30)? XYReconstructor::GetExpectedS2overS2total(par, ich):
//                                           XYReconstructor::GetExpectedS2overS2totalFromGraph(par, ich); // Should we make histogram of solid angle for lookup?
#else
    Double_t expS2overS2total = XYReconstructor::GetExpectedS2overS2totalFromGraph(par, ich);
#endif

#ifdef LOGS2
    weight = 1./0.04;
#else
//    weight = 1./3.5e-3;
    weight = 1.;
#endif
#if 0
    Double_t expS2overS2totalErr = XYReconstructor::GetExpectedS2overS2totalErr(par, ich); // Should we make histogram of solid angle for lookup?
//    weight = (expS2overS2totalErr<0.05)? 1./expS2overS2totalErr: 0.05;
    weight = 1./expS2overS2totalErr;
#endif
//    cout<<expS2overS2totalErr<<endl;
//    if(ich == 30){ // for center channel use actual error
//        Double_t expS2overS2totalErr = XYReconstructor::GetExpectedS2overS2totalErr(par, ich); // Should we make histogram of solid angle for lookup?
//        weight = 1./expS2overS2totalErr;
//    }
    Double_t s2 =fSigProf->s2[ich];// -channel.pulses[1].param.fixed_int2 / channel.pmt.spe_mean;
    if(s2<0.) {chi2=-1.; return;}
//      Double_t weight = 1./TMath::Sqrt(s2/fSigProf->total_s2);
//      Double_t weight = 1./TMath::Sqrt(1./s2+1./fSigProf->total_s2)/(s2/fSigProf->total_s2);
//      Double_t weight = 1./TMath::Log10(expS2overS2totalErr);
#ifdef LOGS2ERR
//    weight = 1./fLogS2FracErr->Eval(TMath::Log10(s2/fSigProf->total_s2)); // from fit to deviation
    weight = (s2>0.)? 2.3026/TMath::Sqrt(1./fSigProf->total_s2 + 1./s2): 0.; //convert s2/total_s2 error to log10(s2/total_s2). 2.3026 is ln10
#endif

#ifdef LOGS2
    Double_t chi = weight*(TMath::Log10(s2/fSigProf->total_s2) - expS2overS2total);
#else
    Double_t chi = weight*(s2/fSigProf->total_s2 - expS2overS2total);
#endif
    chi2 = chi*chi;
}

void XYReconstructor::GetChi2S2overS1(const double *par, Double_t &chi2, Double_t &weight) const {
  //par are supposed to be x and y
#if 1
  Double_t expS2overS1total = XYReconstructor::GetExpectedS2overS1total(par);
    //  Double_t expS2overS1totalErr = XYReconstructor::GetExpectedS2overS1totalErr(par);
    //  Double_t weight = 1./expS2overS1totalErr;
#else
  Double_t s2_x(par[0]), s2_y(par[1]); //x-y position of S2
      Double_t r = TMath::Sqrt(s2_x*s2_x + s2_y*s2_y);
//      Double_t param[] = {14.85, 4.56e-04, 1.}; // for 400<s1_600PE
      Double_t param[] = {1.15436e+01, 2.22915e-04, 1.}; // for 60<s1<300PE
      Double_t expS2overS1total = XYReconstructor::S2overS1(&r, param);
#endif
//      weight = 1./3.;//1.5;//10.;
      weight = 0.5;//1./2.;//1.5;//10.;
      Double_t chi = weight*(fSigProf->s2_over_s1_corr - expS2overS1total);
    chi2 = chi*chi;
}

void XYReconstructor::GetChi2S2BotoverS2(const double *par, Double_t &chi2, Double_t &weight) const {
  //par are supposed to be x and y
  Double_t expS2BotoverS2 = XYReconstructor::GetExpectedS2BotoverS2(par);
//  Double_t expS2BotoverS2Err = 3.*XYReconstructor::GetExpectedS2BotoverS2Err(par);
//  Double_t weight = 1./expS2BotoverS2Err;
//  weight = 1./0.2; // for normal S2_bot/S2_total
//  weight = 1./1.6e-2;
  weight = 1./1.e-3;
//  Double_t weight = 1./TMath::Sqrt(fSigProf->s2_bottom/fSigProf->total_s2);
//  Double_t chi = weight*(fSigProf->s2_bottom/fSigProf->total_s2 - expS2BotoverS2);
  Double_t chi = weight*(TMath::Log10(fSigProf->s2_bottom/fSigProf->total_s2) - expS2BotoverS2);
  chi2 = chi*chi;
}

void XYReconstructor::GetChi2S2BotoverS2Top(const double *par, Double_t &chi2, Double_t &weight) const {
  //par are supposed to be x and y
  Double_t expS2BotoverS2Top = XYReconstructor::GetExpectedS2BotoverS2Top(par);
//  Double_t expS2BotoverS2Err = 3.*XYReconstructor::GetExpectedS2BotoverS2Err(par);
//  Double_t weight = 1./expS2BotoverS2Err;
//  weight = 1./0.2; // for normal S2_bot/S2_total
//  weight = 1./1.6e-2;
//  weight = 0.5e+2; // used in slad v2_2
  Float_t s2_top = fSigProf->total_s2 - fSigProf->s2_bottom;
  weight = (fSigProf->s2_bottom>0. && s2_top>0.)? 2.3026/TMath::Sqrt(1./fSigProf->s2_bottom + 1./s2_top): 0.; //convert s2bot/s2top error to log10(s2bot/s2top). 2.3026 is ln10
//  Double_t weight = 1./TMath::Sqrt(fSigProf->s2_bottom/fSigProf->total_s2);
//  Double_t chi = weight*(fSigProf->s2_bottom/fSigProf->total_s2 - expS2BotoverS2);
  Double_t chi = weight*(TMath::Log10(fSigProf->s2_bottom/s2_top) - expS2BotoverS2Top);
  chi2 = chi*chi;
}

#define SIGNED_R

Double_t XYReconstructor::GetExpectedS2overS2total(const double *par, Int_t ich) const {
  Double_t s2_x(par[0]), s2_y(par[1]); //x-y position of S2
/*  if( ich < PMTGeom::N_CHANNELS/2
//      || PMTGeom::IsOuterCornerPmt(ich)
      ){ // For bottom pmts, use 2D profile
      Double_t val = fProfile2D[ich]->Interpolate(s2_x, s2_y);
      if(val == 0.) val = fProfile2D[ich]->Interpolate(s2_x-0.1, s2_y);
      return val;
  }*/

#ifdef PROFILE2D
#ifndef ROTATIONSYM
  if(!XYReconstructor::CheckTH2Bounds(hLRF[ich], s2_x, s2_y)) return 9999.; //histogram bound check
  return hLRF[ich]->Interpolate(s2_x, s2_y);
#else //#ifndef ROTATIONSYM
  Double_t pmt_phi = TMath::ATan2(PMTGeom::pmt_y[ich], PMTGeom::pmt_x[ich]);
  TVector2 vEvent(s2_x, s2_y);
  TVector2 vEvent_rotate = vEvent.Rotate(-pmt_phi);
  Int_t ctg = PMTGeom::GetPmtPositionCategory(ich);
//  if(ich==37 || ich==34 || ich==36 || ich==32 || ich==21)
  if(ich==37 || ich==34 || ich==36 || ich==32 || ich==33 || ich==22)
    return hLRF[ich]->Interpolate(s2_x, s2_y);
  else
    return hLRF_Rotate[ctg]->Interpolate(vEvent_rotate.X(), vEvent_rotate.Y());
#endif //#ifndef ROTATIONSYM

#else
  //For top pmts, use spline
  Double_t parxy[] = {s2_x, s2_y};
#ifndef SIGNED_R
  Double_t r2 = PMTGeom::GetDist2fromPMT(ich, parxy);
  return fSpline[ich]->Eval(TMath::Sqrt(r2));
#else
//  Double_t r = PMTGeom::GetSignedDistfromPMT(ich, parxy);
  Double_t r(0.), phi(0.);
  PMTGeom::GetPhiFromPmttoCenter(ich, parxy, r, phi);
  Double_t frac = TMath::Abs(phi)/TMath::Pi();// phi is angle between Pmt vector from center and event vector from pmt
  Double_t sign = (frac>0.5)? 1. : -1.;
#if 1
  return (TMath::Abs(r)<5.)? (1.-frac)*fSpline[ich]->Eval(-r)+ frac*fSpline[ich]->Eval(r): fSpline[ich]->Eval(sign*r);
//  return (1.-frac)*fSpline[ich]->Eval(-r)+ frac*fSpline[ich]->Eval(r);
#else
  return fSpline[ich]->Eval(sign*r);
#endif
#endif //SIGNED_R
#endif //PROFILE2D
}

Double_t XYReconstructor::GetExpectedS2overS2totalFromSpline(const double *par, Int_t ich) const {
  Double_t s2_x(par[0]), s2_y(par[1]); //x-y position of S2
  //For top pmts, use spline
  Double_t parxy[] = {s2_x, s2_y};
  Double_t r = PMTGeom::GetSignedDistfromPMT(ich, parxy);
  return fSpline[ich]->Eval(r);

}

Double_t XYReconstructor::GetExpectedS2overS2totalFromGraph(const double *par, Int_t ich) const {
  Double_t s2_x(par[0]), s2_y(par[1]); //x-y position of S2
  //For top pmts, use spline
  Double_t parxy[] = {s2_x, s2_y};

  Double_t r2 = PMTGeom::GetDist2fromPMT(ich, parxy);
  Int_t ctg = PMTGeom::GetPmtPositionCategory(ich);
//  return gLRF_Rotate[ctg]->Eval(TMath::Sqrt(r2));
  return avoid_tgraph_insane_extrapolation(TMath::Sqrt(r2), *gLRF_Rotate[ctg]);
}

Double_t XYReconstructor::avoid_tgraph_insane_extrapolation (float x, const TGraph& g) const { // imported from Alessandro's code
  if (x <= g.GetX ()[0]) return g.GetY ()[0];
  int n = g.GetN ();
  if (x >= g.GetX ()[n - 1]) return g.GetY ()[n - 1];

  return g.Eval (x);
}


Double_t XYReconstructor::GetExpectedS2overS2totalErr(const double *par, Int_t ich) const {
  Double_t s2_x(par[0]), s2_y(par[1]); //x-y position of S2
  if( ich < PMTGeom::N_CHANNELS/2
//      || PMTGeom::IsOuterCornerPmt(ich)
      ){ // For bottom pmts, use 2D profile
     Double_t err = 5.*fProfile2D[ich]->GetBinError(fProfile2D[ich]->FindBin(s2_x, s2_y));// 30. is factor to make error comparable to the one from spline
//     Double_t err = fProfile2D[ich]->GetBinError(fProfile2D[ich]->FindBin(s2_x, s2_y));// 30. is factor to make error comparable to the one from spline
      if(err < 1.e-12) err = 1.e-12;
      return err;// 30. is factor to make error comparable to the one from spline
  }
    //  else
//    return 3.*fProfile2D[ich]->GetBinError(fProfile2D[ich]->FindBin(s2_x, s2_y));// 30. is factor to make error comparable to the one from spline

  Double_t fac = (PMTGeom::IsOuterCornerPmt(ich))? 1. : 1.;
//#ifdef PROFILE2D
#if 0
  if(!XYReconstructor::CheckTH2Bounds(hLRF_Err[ich], s2_x, s2_y)) return 9999.; //histogram bound check
  return fac*hLRF_Err[ich]->Interpolate(s2_x, s2_y);
#else

  //For top pmts, use spline
  Double_t parxy[] = {s2_x, s2_y};
//#ifndef SIGNED_R
#if 1
  Double_t r2 = PMTGeom::GetDist2fromPMT(ich, parxy);
//  return fac*fSpline_sigma[ich]->Eval(TMath::Sqrt(r2));
  return fac*fSpline_sigma[30]->Eval(TMath::Sqrt(r2)); //use error from ch30
//  return 10.*fSpline_sigma[ich]->Eval(TMath::Sqrt(r2));
#else
  Double_t r = PMTGeom::GetSignedDistfromPMT(ich, parxy);
//  return fac*fSpline_sigma[ich]->Eval(r);
  return fac*fSpline_sigma[30]->Eval(r); //use error from ch30
#endif //SIGNED_R
#endif //PROFILE2D
}

Double_t XYReconstructor::GetExpectedS2overS1total(const double *par) const {
  Double_t s2_x(par[0]), s2_y(par[1]); //x-y position of S2
  Double_t r2 = s2_x*s2_x + s2_y*s2_y;
  if(TMath::Sqrt(r2)>20.) return 9999.; // 20 cm is limit of profile
#if 0
  return pS2overS1VsR->Interpolate(TMath::Sqrt(r2));
#else
  return fS2overS1VsR->Eval(TMath::Sqrt(r2));
#endif
}

Double_t XYReconstructor::GetExpectedS2overS1totalErr(const double *par) const {
  Double_t s2_x(par[0]), s2_y(par[1]); //x-y position of S2
  Double_t r2 = s2_x*s2_x + s2_y*s2_y;
  if(TMath::Sqrt(r2)>20.) return -1.; // 20 cm is limit of profile

//  return 1.5*pS2overS1VsR->GetBinError(pS2overS1VsR->FindBin(TMath::Sqrt(r2)));
  return 0.5*pS2overS1VsR->GetBinError(pS2overS1VsR->FindBin(TMath::Sqrt(r2)));
//  return pS2overS1VsR->GetBinError(pS2overS1VsR->FindBin(TMath::Sqrt(r2)));
}

Double_t XYReconstructor::GetExpectedS2BotoverS2(const double *par) const {
  Double_t s2_x(par[0]), s2_y(par[1]); //x-y position of S2
//  cout<<"x: "<<s2_x<<" y: "<<s2_y<<endl;
  if(!XYReconstructor::CheckTH2Bounds(hS2BotOverS2VsXY, s2_x, s2_y)) return 9999.; //histogram bound check
  return hS2BotOverS2VsXY->Interpolate(s2_x, s2_y);
}

Double_t XYReconstructor::GetExpectedS2BotoverS2Err(const double *par) const {
  Double_t s2_x(par[0]), s2_y(par[1]); //x-y position of S2
  if(!XYReconstructor::CheckTH2Bounds(hS2BotOverS2VsXYErr, s2_x, s2_y)) return 9999.; //histogram bound check
   return hS2BotOverS2VsXYErr->Interpolate(s2_x, s2_y);
}

Double_t XYReconstructor::GetExpectedS2BotoverS2Top(const double *par) const {
  Double_t s2_x(par[0]), s2_y(par[1]); //x-y position of S2
//  cout<<"x: "<<s2_x<<" y: "<<s2_y<<endl;
  if(!XYReconstructor::CheckTH2Bounds(hS2BotOverS2TopVsXY, s2_x, s2_y)) return 9999.; //histogram bound check
  return hS2BotOverS2TopVsXY->Interpolate(s2_x, s2_y);
}

Double_t XYReconstructor::GetDensity(const double *par) const {
  Double_t s2_x(par[0]), s2_y(par[1]); //x-y position of S2
  if(!XYReconstructor::CheckTH2Bounds(hRecPosiXY, s2_x, s2_y)) return 9999.; //histogram bound check
   return hRecPosiXY->Interpolate(s2_x, s2_y);
}

void XYReconstructor::CorrectQEforS2(SignalProfile &sigProf){
   //Loop over channels
  Double_t s2total(0.), s2top(0.), s2bot(0.);
  sigProf.max_s2 = 0.;// initialize
  for (Int_t ch = 0; ch < PMTGeom::N_CHANNELS; ch++)
  {
     sigProf.s2[ch] *= PMTGeom::pmt_qe_norm[ch];
     s2total+=sigProf.s2[ch];
     if (ch < PMTGeom::N_CHANNELS/2) s2bot += sigProf.s2[ch];
     else s2top += sigProf.s2[ch];

     if (sigProf.s2[ch] > sigProf.max_s2)
     {
         sigProf.max_s2 = sigProf.s2[ch];
         sigProf.max_s2_chan = ch;
     }

  }
  sigProf.s2_bottom = s2bot;
  sigProf.total_s2 = s2total;
}

void XYReconstructor::ChargeCorrectforS2(SignalProfile &sigProf){
   //Loop over channels
  Double_t s2total(0.), s2top(0.), s2bot(0.);
  sigProf.max_s2 = 0.;// initialize
  for (Int_t ch = 0; ch < PMTGeom::N_CHANNELS; ch++)
  {
     if(ch >= PMTGeom::N_CHANNELS/2) sigProf.s2[ch] /= PMTGeom::charge_normalization[ch%(PMTGeom::N_CHANNELS/2)];// for top pmt
     else sigProf.s2[ch] *= PMTGeom::pmt_qe_norm[ch];
     s2total+=sigProf.s2[ch];
     if (ch < PMTGeom::N_CHANNELS/2) s2bot += sigProf.s2[ch];
     else s2top += sigProf.s2[ch];

     if (sigProf.s2[ch] > sigProf.max_s2)
     {
         sigProf.max_s2 = sigProf.s2[ch];
         sigProf.max_s2_chan = ch;
     }

  }
  sigProf.s2_bottom = s2bot;
  sigProf.total_s2 = s2total;
}

void XYReconstructor::InverseChargeCorrectforS2(SignalProfile &sigProf){
   //Loop over channels
  Double_t s2total(0.), s2top(0.), s2bot(0.);
  sigProf.max_s2 = 0.;// initialize
  for (Int_t ch = 0; ch < PMTGeom::N_CHANNELS; ch++)
  {
     if(ch >= PMTGeom::N_CHANNELS/2) sigProf.s2[ch] *= PMTGeom::charge_normalization[ch%(PMTGeom::N_CHANNELS/2)];// for top pmt
     else sigProf.s2[ch] /= PMTGeom::pmt_qe_norm[ch];
     s2total+=sigProf.s2[ch];
     if (ch < PMTGeom::N_CHANNELS/2) s2bot += sigProf.s2[ch];
     else s2top += sigProf.s2[ch];

     if (sigProf.s2[ch] > sigProf.max_s2)
     {
         sigProf.max_s2 = sigProf.s2[ch];
         sigProf.max_s2_chan = ch;
     }

  }
  sigProf.s2_bottom = s2bot;
  sigProf.total_s2 = s2total;
}

Double_t XYReconstructor::XYCorrection4S2(Double_t x, Double_t y) const {
  if(!XYReconstructor::CheckTH2Bounds(hS2OverS1VsXY_norm_center, x, y)) return 9999.; //histogram bound check
//  return 1./hS2OverS1VsXY_norm_center->Interpolate(x, y);
  return 1./hS2OverS1VsXY_norm_center->GetBinContent(hS2OverS1VsXY_norm_center->FindBin(x, y));
}


void XYReconstructor::LoadProfile(TString filename) {
  TFile *file = new TFile(filename.Data());
  if(!file) cout<<"file: "<<filename.Data()<<" is not found."<<endl;
  else cout<<"profile:  "<<filename.Data()<<" is used."<<endl;
  fProfile2D = new TProfile2D*[PMTGeom::N_CHANNELS];
  fSpline = new TSpline5*[PMTGeom::N_CHANNELS];
  fSpline_sigma = new TSpline5*[PMTGeom::N_CHANNELS];
  for(Int_t ich=0; ich<PMTGeom::N_CHANNELS; ich++){
#ifdef LOGS2
      TString hname = Form("pLogS2OverS2TotVsXY_ch%d", ich);
#else
      TString hname = Form("pS2OverS2TotVsXY_ch%d", ich);
#endif
      fProfile2D[ich] = (TProfile2D*) file->Get(hname.Data());
      if(!fProfile2D[ich]) {
          cout<<Form("profile: %s is not found.", hname.Data())<<endl;
      }
      fProfile2D[ich]->SetName(Form("fProfile2D_ch%d", ich));
      fProfile2D[ich]->SetDirectory(0);
  }

  hLRF = new TH2D*[PMTGeom::N_CHANNELS];
  hLRF_Err = new TH2D*[PMTGeom::N_CHANNELS];
  for(Int_t ich=PMTGeom::N_CHANNELS/2; ich<PMTGeom::N_CHANNELS; ich++){ //only top pmts
//   for(Int_t ich=0; ich<PMTGeom::N_CHANNELS; ich++){ //only top pmts
       TString hname;
       if(ich>=PMTGeom::N_CHANNELS/2){
           hname = Form("spS2OverS2TotVsR_ch%d", ich);
           fSpline[ich] = (TSpline5*) file->Get(hname.Data());
           if(!fSpline[ich]) cout<<Form("spline: %s is not found.", hname.Data())<<endl;
           fSpline[ich]->SetName(Form("fSpline_ch%d", ich));
     //      fSpline[ich]->SetDirectory(0);
           hname = Form("spS2OverS2TotVsR_sigma_ch%d", ich);
           fSpline_sigma[ich] = (TSpline5*) file->Get(hname.Data());
           if(!fSpline_sigma[ich]) cout<<Form("spline_sigma: %s is not found.", hname.Data())<<endl;
           fSpline_sigma[ich]->SetName(Form("fSpline_sigma_ch%d", ich));
       }

#ifdef PROFILE2D
      hname = Form("hS2OverS2TotVsXY_ch%d_Interp", ich);
      hLRF[ich] = (TH2D*) file->Get(hname.Data());
      if(!hLRF[ich]){
          cout<<Form("hLRF: %s is not found.", hname.Data())<<endl;
          hname = Form("hS2OverS2TotVsXY_ch%d", ich);
//          hname = Form("hLogS2OverS2TotVsXY_ch%d", ich);
          hLRF[ich] = (TH2D*) file->Get(hname.Data());
          if(!hLRF[ich]) cout<<Form("hLRF: %s is not found.", hname.Data())<<endl;
          else cout<<Form("hLRF: %s is found.", hname.Data())<<endl;
      }
      hLRF[ich]->SetName(Form("hLRF_ch%d", ich));
      hLRF[ich]->SetDirectory(0);

      hname = Form("hS2OverS2TotVsXYErr_ch%d_Interp", ich);
      hLRF_Err[ich] = (TH2D*) file->Get(hname.Data());
      if(!hLRF_Err[ich]) {
          cout<<Form("hLRF_Err: %s is not found.", hname.Data())<<endl;
          hname = Form("hS2OverS2TotVsXYErr_ch%d", ich);
//          hname = Form("hLogS2OverS2TotVsXYErr_ch%d", ich);
          hLRF_Err[ich] = (TH2D*) file->Get(hname.Data());
          if(!hLRF_Err[ich]) cout<<Form("hLRF_Err: %s is not found.", hname.Data())<<endl;
          else cout<<Form("hLRF_Err: %s is found.", hname.Data())<<endl;
      }
      hLRF_Err[ich]->SetName(Form("hLRF_Err_ch%d", ich));
      hLRF_Err[ich]->SetDirectory(0);
#endif
  }

#ifdef S2OVERS1PROFILE
  pS2overS1VsR = (TProfile*) file->Get("hS2overS1VsR_1");
  if(!pS2overS1VsR) {
      cout<<"profile hS2overS1VsR_1 is not found."<<endl;
  } else{
      pS2overS1VsR->SetName("pS2overS1VsR");
      pS2overS1VsR->SetDirectory(0);
  }

  fS2overS1VsR = (TF1*) file->Get("fS2overS1_R");
  if(!fS2overS1VsR) {
      cout<<"profile fS2overS1_R is not found."<<endl;
  } else{
      fS2overS1VsR->SetName("fS2overS1VsR");
  }
#endif

#ifdef S2BOTOVERS2PROFILE
  hS2BotOverS2VsXY = (TH2D*) file->Get("hS2BotOverS2VsXY_Interp");
  if(!hS2BotOverS2VsXY) {
      cout<<"profile hS2BotOverS2VsXY_Interp is not found."<<endl;
      hS2BotOverS2VsXY = (TH2D*) file->Get("hLogS2BotOverS2VsXY_Interp");
      if(!hS2BotOverS2VsXY) cout<<"profile hLogS2BotOverS2VsXY_Interp is not found."<<endl;
      else cout<<"profile hLogS2BotOverS2VsXY_Interp is found."<<endl;
  }
  hS2BotOverS2VsXY->SetName("hLRF_S2BotOverS2");
  hS2BotOverS2VsXY->SetDirectory(0);


  hS2BotOverS2VsXYErr = (TH2D*) file->Get("hS2BotOverS2VsXYErr_Interp");
  if(!hS2BotOverS2VsXYErr) {
      cout<<"profile hS2BotOverS2VsXYErr_Interp is not found."<<endl;
      hS2BotOverS2VsXYErr = (TH2D*) file->Get("hLogS2BotOverS2VsXYErr_Interp");
      if(!hS2BotOverS2VsXYErr) cout<<"profile hLogS2BotOverS2VsXYErr_Interp is not found."<<endl;
      else cout<<"profile hLogS2BotOverS2VsXYErr_Interp is found."<<endl;
  }
  hS2BotOverS2VsXYErr->SetName("hLRF_S2BotOverS2Err");
  hS2BotOverS2VsXYErr->SetDirectory(0);

#endif

#ifdef S2BOTOVERS2TOPPROFILE
  hS2BotOverS2TopVsXY = (TH2D*) file->Get("hS2BotOverS2TopVsXY_Interp");
  if(!hS2BotOverS2TopVsXY) {
      cout<<"profile hS2BotOverS2TopVsXY_Interp is not found."<<endl;
//      hS2BotOverS2TopVsXY = (TH2D*) file->Get("hLogS2BotOverS2TopVsXY");
      hS2BotOverS2TopVsXY = (TH2D*) file->Get("hLogS2BotOverS2TopVsXY_Interp");
//      hS2BotOverS2TopVsXY = (TH2D*) file->Get("hLogS2BotOverS2TopVsXY_Interp_Forced");
      if(!hS2BotOverS2TopVsXY) cout<<"profile hLogS2BotOverS2TopVsXY_Interp is not found."<<endl;
      else cout<<"profile hLogS2BotOverS2TopVsXY_Interp is found."<<endl;
  }
  hS2BotOverS2TopVsXY->SetName("hLRF_S2BotOverS2Top");
  hS2BotOverS2TopVsXY->SetDirectory(0);

#endif

#ifdef DENSITY
  hRecPosiXY = (TH2D*) file->Get("hRecPosiXY");
  if(!hRecPosiXY) {
      cout<<"profile hRecPosiXY is not found."<<endl;
  }
  hRecPosiXY->SetName("XYDensity");
  hRecPosiXY->SetDirectory(0);
  if(hRecPosiXY->GetNbinsX()<25) hRecPosiXY->Rebin2D(4,4); // too small bins
#endif

#ifdef LOGS2ERR
  fLogS2FracErr = (TF1*) file->Get("fLogS2FracErr");
  if(!fLogS2FracErr) {
      cout<<"TF1 fLogS2FracErr is not found."<<endl;
  }
#endif

//#ifdef ROTATIONSYM
#if 0
  hLRF_Rotate = new TH2D*[PMTGeom::nPmtCategories];
  gLRF_Rotate = new TGraph*[PMTGeom::nPmtCategories];
  for(Int_t iCatg=0; iCatg<PMTGeom::nPmtCategories; iCatg++){
      TString hname = Form("hLogS2OverS2TotVsXY_Rotate_%s_Interp", PMTGeom::pmtlabels[iCatg].Data());
      hLRF_Rotate[iCatg] = (TH2D*) file->Get(hname.Data());
      if(!hLRF_Rotate[iCatg]) cout<<"profile: "<<hname.Data()<<" is not found."<<endl;
      hLRF_Rotate[iCatg]->SetName(Form("hLRF_Rotate_%s", PMTGeom::pmtlabels[iCatg].Data()));
      hLRF_Rotate[iCatg]->SetDirectory(0);

      hname = Form("gLogS2OverS2Tot_Vs_R_%s", PMTGeom::pmtlabels[iCatg].Data());
      gLRF_Rotate[iCatg] = (TGraph*) file->Get(hname.Data());
      if(!gLRF_Rotate[iCatg]) cout<<"profile: "<<hname.Data()<<" is not found."<<endl;
      gLRF_Rotate[iCatg]->SetName(Form("gLRF_Rotate_%s", PMTGeom::pmtlabels[iCatg].Data()));

  }
#endif


  hS2OverS1VsXY_norm_center = (TH2D*) file->Get("hS2OverS1VsXY_norm_center");
  if(!hS2OverS1VsXY_norm_center) cout<<"profile: hS2OverS1VsXY_norm_center is not found."<<endl;
  hS2OverS1VsXY_norm_center->SetDirectory(0);

  file->Close();
  delete file;
}


Bool_t XYReconstructor::Is_S2OVERS1PROFILE_Used() const{
#ifdef S2OVERS1PROFILE //switch to use S2/S1 vs R profile
  return true;
#endif
  return false;

}
Bool_t XYReconstructor::Is_S2BOTOVERS2PROFILE_Used() const{
#ifdef S2BOTOVERS2PROFILE //switch to use S2_bottom/S2 vs R profile
  return true;
#endif
  return false;
}
Bool_t XYReconstructor::Is_S2BOTOVERS2TOPPROFILE_Used() const {
#ifdef S2BOTOVERS2TOPPROFILE //switch to use S2_bottom/S2_top vs xy profile
  return true;
#endif
  return false;

}
Bool_t XYReconstructor::Is_PROFILE2D_Used() const {
#ifdef PROFILE2D //switch to use S2_sh/S2 vs XY profile, not 1D but 2D profile
  return true;
#endif
  return false;

}
Bool_t XYReconstructor::Is_DENSITYPROFILE_Used() const {
#ifdef DENSITY //switch to use density info
  return true;
#endif
  return false;

}
Bool_t XYReconstructor::Is_LOGS2_Used() const {
#ifdef LOGS2 //switch to use Log(S2ch/S2tot)
  return true;
#endif
  return false;

}
Bool_t XYReconstructor::Is_LOGS2ERR_Used() const {
#ifdef LOGS2ERR //switch to use error of Log(S2ch/S2tot)
  return true;
#endif
  return false;

}
Bool_t XYReconstructor::Is_ROTATIONSYM_Used() const {
#ifdef ROTATIONSYM //use 2d profile from pmt categories
  return true;
#endif
  return false;

}
Bool_t XYReconstructor::Is_USEALEMETHOD_Used() const {
#ifdef USEALEMETHOD //use Alessandro's 1d profile. This need to use with ROTATIONSYM
  return true;
#endif
  return false;

}

