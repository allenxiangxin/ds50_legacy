#ifndef XYReconstructor_h
#define XYReconstructor_h

#include <iostream>

#include "TROOT.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TSpline.h"
#include "TF1.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/IFunction.h"

#include <vector>
#include <map>

//
//#include "/ds50/app/user/xxiang89/work/xy4slad_ms/XYReconstructor/PMTGeom.h"
//#include "/ds50/app/user/xxiang89/work/xy4slad_ms/XYReconstructor/SignalProfile.h"

#include "PMTGeom.h"
#include "SignalProfile.h"

using namespace std;

class XYReconstructor : public ROOT::Math::IMultiGenFunction {

public :

  XYReconstructor(): fSigProf(NULL), ndf_inchi2(0),
  fProfile2D(0), fSpline(0), fSpline_sigma(0), pS2overS1VsR(0), fS2overS1VsR(0), fLogS2FracErr(0),
  hLRF(0), hLRF_Err(0), hLRF_Rotate(0), gLRF_Rotate(0),
  hS2BotOverS2VsXY(0), hS2BotOverS2VsXYErr(0), hS2BotOverS2TopVsXY(0), hRecPosiXY(0), hS2OverS1VsXY_norm_center(0)
  {
    //TF2(const char* name, void* ptr, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Int_t npar, const char* className)
    gf2D = new TF2("gf2D", this, -1.2*PMTGeom::tpcRadi, 1.2*PMTGeom::tpcRadi, -1.2*PMTGeom::tpcRadi, 1.2*PMTGeom::tpcRadi, 0, "XYReconstructor");
    gf2D->SetNpx(180);
    gf2D->SetNpy(180);

  }

  ROOT::Math::Minimizer *fMinimizer;

  double operator() (double *x, double *) const { return DoEval(x); }
  void SetSignalProfile(SignalProfile* sigProf){ fSigProf = sigProf; }
  void CorrectQEforS2(SignalProfile &sigProf);
  void ChargeCorrectforS2(SignalProfile &sigProf);
  void InverseChargeCorrectforS2(SignalProfile &sigProf);
  Double_t XYCorrection4S2(Double_t x, Double_t y) const;

  void LoadProfile(TString filename);
  unsigned int NDim() const { return 2; }
  unsigned int NDF() const { return ndf_inchi2; }

  virtual Int_t Minimization(SignalProfile &sigProf);
  virtual Int_t SetMinimizer();
  Int_t MinimizationStatus();
  Double_t ReducedChi2();
  void Minimize2DFunc(Double_t& x1, Double_t& x2);
  void Scan2DFunc(TF2 &f2D, Double_t& x1, Double_t& x2);
  Double_t FindMaxChi2(TF2 &f2D, Double_t& x1, Double_t& x2);
  void PrintComponents (const double * x) const ;
  Bool_t IsChanIncludedInChi2 (Int_t ich) const ;
  Int_t CountDF() const;

  void GetDensityContribution(const double *par, Double_t &chi2, Double_t &weight) const ;
  void GetChi2Channel(const double *par, Int_t ich, Double_t &chi2, Double_t &weight) const ;
  void GetChi2S2overS1(const double *par, Double_t &chi2, Double_t &weight) const ;
  void GetChi2S2BotoverS2(const double *par, Double_t &chi2, Double_t &weight) const ;
  void GetChi2S2BotoverS2Top(const double *par, Double_t &chi2, Double_t &weight) const ;
  Double_t GetExpectedS2overS2total(const double *par, Int_t ich) const ;
  Double_t GetExpectedS2overS2totalErr(const double *par, Int_t ich) const ;
  Double_t GetExpectedS2overS2totalFromSpline(const double *par, Int_t ich) const ;
  Double_t GetExpectedS2overS2totalFromGraph(const double *par, Int_t ich) const ;
  Double_t avoid_tgraph_insane_extrapolation (float x, const TGraph& g) const ;

  Double_t GetExpectedS2overS1total(const double *par) const ;
  Double_t GetExpectedS2overS1totalErr(const double *par) const ;

  Double_t GetExpectedS2BotoverS2(const double *par) const ;
  Double_t GetExpectedS2BotoverS2Err(const double *par) const ;

  Double_t GetExpectedS2BotoverS2Top(const double *par) const ;

  Double_t GetDensity(const double *par) const ;

   ROOT::Math::IMultiGenFunction * Clone() const {
      return new XYReconstructor();
   }

   Bool_t Is_S2OVERS1PROFILE_Used() const;
   Bool_t Is_S2BOTOVERS2PROFILE_Used() const;
   Bool_t Is_S2BOTOVERS2TOPPROFILE_Used() const;
   Bool_t Is_PROFILE2D_Used() const;
   Bool_t Is_DENSITYPROFILE_Used() const;
   Bool_t Is_LOGS2_Used() const;
   Bool_t Is_LOGS2ERR_Used() const;
   Bool_t Is_ROTATIONSYM_Used() const;
   Bool_t Is_USEALEMETHOD_Used() const;

private:

   Double_t DoEval (const double * x) const ;
   Double_t GetSolidAngle(const double *x, const double *par) const ;
   Double_t S2overS1(Double_t *x, Double_t *par) ;
   Bool_t CheckTH2Bounds(TH2* h, Double_t x, Double_t y) const ;

   SignalProfile* fSigProf; //!
   TF2           *gf2D; //!

   Int_t ndf_inchi2; // number of degree of freedom in chi2 (pmt channel + S2bot/S2top + etc....)

   TProfile2D** fProfile2D; //! s2_ch/s2_total profile
   TSpline5**   fSpline; //! s2_ch/s2_total spline
   TSpline5**   fSpline_sigma; //! s2_ch/s2_total spline
   TProfile*    pS2overS1VsR; //! total_s2/total_s1 profile
   TF1*         fS2overS1VsR; //! total_s2/total_s1 profile
   TF1*         fLogS2FracErr; //! error of Log(S2/S2_total) function of Log(S2/S2_total)

   TH2D** hLRF; //! s2_ch/s2_total profile
   TH2D** hLRF_Err; //! s2_ch/s2_total profile
   TH2D** hLRF_Rotate; //! s2_ch/s2_total profile for each pmt category
   TGraph** gLRF_Rotate; //!
   TH2D *hS2BotOverS2VsXY, *hS2BotOverS2VsXYErr; //! s2_bottom/s2_total profile
   TH2D *hS2BotOverS2TopVsXY; //! s2_bottom/s2_top profile
   TH2D *hRecPosiXY; //! density
   TH2D *hS2OverS1VsXY_norm_center; //! S2 correction histogram

};

ClassImp(XYReconstructor)

#endif




