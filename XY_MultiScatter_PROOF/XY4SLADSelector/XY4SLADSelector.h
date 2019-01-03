

#ifndef XY4SLADSelector_h
#define XY4SLADSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TRint.h>

#include "TStopwatch.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include <iostream>

#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include <TVector2.h>

//#include "/ds50/app/user/xxiang89/work/xy4slad_ms/XYReconstructor/XYRecHelper.h"
//#include "/ds50/app/user/xxiang89/work/xy4slad_ms/XYReconstructor/SignalProfile.h"

#include "../XYReconstructor/XYRecHelper.h"
#include "../XYReconstructor/SignalProfile.h"

using namespace std;


// Header file for the classes stored in the TTree if any.
const Int_t MAX_N_PULSES = 47;
const float def_s2_dz =33;
const float minchi2_cx =2;
const int s1_min=20;
const int s1_max = 820;
const int s1_width =20;
const int n_slice = (s1_max-s1_min)/s1_width;
const int N_Evt=100;
const int n_chi2cx = 2;




class XY4SLADSelector : public TSelector {
  public :
  TTree*            fChain;                   //!pointer to the analyzed TTree or TChain
  TFile*            foutfile;                 //!pointer to output file which save histograms
  TString           fOutName;
  XYReconstructor*  fMyFcn;                   //!pointer to Masa's Recon class
  TFile*            fs1mfFile;                //!pointer to file which saved S1 max frac def
  TH2F*             h_c95_s1mf_thresholds;    //!pointer to 95 percentile s1 max frac def
  TH2F*             h_c99_s1mf_thresholds;    //!pointer to 99 percentile s1 max frac


  
  //helper variables
  int N_Entries;
  int counter;
  bool Debug;
  
  
  // variables from main file
  int    run_id;
  int    event_id;
  int    nchannels;
  short  baseline_not_found;
  double live_time;
  double inhibit_time;
  int    npulses;
  int    has_s3;
  short  hasV1724;
  float  s1_start_time;
  float  total_s1;
  float  total_f90;
  float  total_s2;
  float  total_s1_corr;
  float  total_s2_corr;
  float  total_s1_top;
  float  total_s1_bottom;
  float  total_s2_top;
  float  total_s2_bottom;
  double  tdrift;
  
  int    max_s1_frac_cut_exceeds99;
  int    max_s1_frac_cut_exceeds95;
  float  c95_s1mf_thresholds;
  float  c99_s1mf_thresholds;
  int    s1_max_chan;
  float  s1_max_frac;
  
  // variables from s2 file
  float   s2_chan[N_CHAN];
  float   s2_max_frac;
  
  // variables from masa's xy file
  float masas_x;
  float masas_y;
  float masas_chi2;
  float allpulses_x[MAX_N_PULSES];
  float allpulses_y[MAX_N_PULSES];
  float allpulses_chi2[MAX_N_PULSES];
  float allpulses_xycorr_factor[MAX_N_PULSES];
  
  // variables from jason's xy file
  float xyl_best_x;
  float xyl_best_y;
  float xyl_best_chi2;
  float allpulses_xyl_x[MAX_N_PULSES];
  float allpulses_xyl_y[MAX_N_PULSES];
  float allpulses_xyl_chi2[MAX_N_PULSES];
  float allpulses_xyl_xycorr_factor[MAX_N_PULSES];
  
  
  // variables from allpulses file
  int   pulse_npulses;
  float pulse_start_time[MAX_N_PULSES];
  float pulse_total_npe[MAX_N_PULSES];
  float pulse_f90[MAX_N_PULSES];
  int   pulse_max_chan[MAX_N_PULSES];
  float pulse_fixed_int1[MAX_N_PULSES];
  float pulse_fixed_int2[MAX_N_PULSES];
  int   pulse_saturated[MAX_N_PULSES];
  int   pulse_ch_light_size;
  float pulse_ch_light[N_CHAN*MAX_N_PULSES];
  
  //variable calculated on the fly
  float s2_xycorr[MAX_N_PULSES];  //xy correction to all pulses;
  float x[MAX_N_PULSES];
  float y[MAX_N_PULSES];
  float s2_dx;
  float s2_dy;
  float s2_dz;      //spatial separation of two s2, default to be 0
  float s2_dr;      //xy separation of two s2
  Double_t TBChi2;  //Masa's Bot/Top Chi2: from void XYReconstructor::GetChi2S2BotoverS2Top
  
  
  //init variables
  XY4SLADSelector(): fChain(0), fOutName(""), fMyFcn(0), counter(0), Debug(0), tdrift(0.), s2_dx(0.), s2_dy(0.),  s2_dz(0.),  s2_dr(0.), TBChi2(0.), c99_s1mf_thresholds(1.), c95_s1mf_thresholds(1.),

  //init histogram pointer to NULL
  h_tdrift(0),  h_s2_ch_occ(0),  h_xy_sing(0),  h_xy_comb(0), h_xy_sing_chi2cx(0), h_xy_comb_chi2cx(0), h_sladxy_sing_uar(0), h_jasonxy_sing(0),
//  h_xy[0] = 0, h_xy[1](0),
  h_p0_start(0), h_chi2_dr(0),  h_logchi2_dr(0),  h_TBchi2_dr(0),  h_rerr_dr(0),  h_rerr_r_comb(0),  h_rerr_r_sing(0),  h_maxerr_dr(0),
    h_maxerr_r_comb(0),h_maxerr_r_sing(0),  h_top_var_dr(0),  h_top_var_r_sing(0),  h_top_var_r_comb(0),
//  h_top_var_r[0](0), h_top_var_r[1](0),
  h_bot_var_dr(0),  h_bot_var_r_sing(0),  h_bot_var_r_comb(0),
//  h_bot_var_r[0](0), h_bot_var_r[0](0),
  h_topbot_var_dr(0),  h_topbot_var_r_sing(0),  h_topbot_var_r_comb(0),
//  h_topbot_var_r[0](0), h_topbot_var_r[1](0),
  h_chi2_comb(0),  h_chi2_sing(0),  h_sladchi2_sing_uar(0),  h_TBchi2_comb(0),  h_TBchi2_sing(0),  h_TB_comb(0),  h_TB_sing(0),
//  h_chi2[0](0), h_chi2[1](0),
  
  h_rchi2_sing(0),  h_rchi2_sing_aar(0),  h_rchi2_comb(0),  h_rchi2_chi2_sing(0),  h_rchi2_chi2_sing_aar(0),  h_rchi2_chi2_comb(0),
    h_fracchi2_chi2_sing(0), h_fracchi2_chi2_sing_aar(0),  h_fracchi2_chi2_comb(0),  h_fracchi2_logchi2_sing(0),
  h_fracchi2_logchi2_sing_aar(0),  h_fracchi2_logchi2_comb(0),  h_maxchi2_chi2_sing(0),  h_maxchi2_chi2_sing_aar(0),  h_maxchi2_chi2_comb(0),  h_s1_comb(0),  h_s1_comb_goodxy(0),  h_s1_comb_badxy(0),  h_s1_sing(0),  h_s1_sing_goodxy(0),  h_s1_sing_badxy(0),  h_dexpxy_dr(0),  h_R(0),  h_R_masas(0),  h_R_jasons(0),  h_R_s1cx(0),  h_R_masas_s1cx(0),  h_R_jasons_s1cx(0),
  
  h_chi2_sing_aar(0),  h_sladchi2_sing_aar(0),  h_xy_sing_aar(0),  h_xy_sing_aar_chi2cx(0), h_sladxy_sing_aar(0),  h_jasonxy_sing_aar(0),  h_top_var_r_sing_aar(0),  h_bot_var_r_sing_aar(0),  h_topbot_var_r_sing_aar(0),  h_rerr_r_sing_aar(0),  h_maxerr_r_sing_aar(0),  h_TBchi2_sing_aar(0),  h_TB_sing_aar(0),
    h_maxvar_chi2_sing(0), h_maxvar_chi2_comb(0), h_maxvar_chi2_sing_aar(0), h_rvar_chi2_sing(0), h_rvar_chi2_comb(0), h_rvar_chi2_sing_aar(0),  
  h_f90_s1_uar(0)

  {
    h_chi2[0]=0; h_chi2[1]=0;
    h_s2s1_vs_s1[0]=0; h_s2s1_vs_s1[1]=0;
  
    for (int i=0; i<n_chi2cx; i++) {
      h_xy[i] = 0;
      h_top_var_r[i]=0;
      h_bot_var_r[i]=0;
      h_topbot_var_r[i]=0;
      h_f90_s1_chi2cx_uar[i]=0;
    }
  
    for (int i=0; i<n_slice; i++) {
      h_f90_uar[i] = 0;
      for (int j=0; j<n_chi2cx; j++) {
        h_f90_uar_chi2cx[j][i]=0;
      }
    }
  
    for (int i=0; i<N_Evt; i++) {
      h_xy_sing1_ev[i] = 0;
      h_xy_sing2_ev[i] = 0;
      h_xy_comb_ev[i] = 0;
    }
  };
  
  
  virtual ~XY4SLADSelector() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();
  
  void SetDebug(Bool_t debug){Debug = debug;}
  void SetOutName(TString str) {fOutName = str;}
  void SetS1MaxFracFileName(TString str) {fs1mfName = str;}
  void SetProfileName(TString fProfile) { fProfileName = fProfile; }

  
  
private:
  
  Bool_t fMinimizeFlag;
  TString fProfileName;
  TString fs1mfName;
  
  Bool_t          IsValidEvent();
  Bool_t          IsEventOK();
  void            BookHistograms();
  void            FillHistograms(SignalProfile* sigProf);
  Double_t        s1_corr_factor(Double_t t_drift_max, Double_t t_drift);
  Double_t        s1_TBAcorr_factor(Double_t s1_top, Double_t s1_bottom, Double_t s1);
  Float_t         GetS2LogTopBot(SignalProfile const& sp);
  Float_t         GetTopPMTVariance(SignalProfile sigProf);
  Float_t         GetTopPMTVariance(int pls_index);
  Float_t         GetBotPMTVariance(SignalProfile sigProf);
  Float_t         GetBotPMTVariance(int pls_index);
  Float_t         GetPMT1DMaxVariance(SignalProfile);
  Float_t         GetPMT1DMaxVariance(int pls_index);  
  vector<float>   GetPMTMaxMinVariance(SignalProfile);

  bool CX_UAR(); bool CX_AAR();
  bool CX1();  bool CX2();  bool CX3(); bool CX4();
  bool CX_2PLS(); bool CX_3PLS();
  bool CX_TDrift_STD(); bool CX_TDrift();
  bool CX_TRIG();
  bool CX_Sat();  bool CX_S2Sat();
  
  bool CX_C99S1MF(); bool CX_C95S1MF();
  bool CX_S2True();  bool CX_S2CorrTrue();
  bool CX_S1F90();
  bool CX_c99S1MF();
  bool CX_S2dZ();
  
  bool CX_R();
  bool CX_Surface();
  bool CX_GoodXY(); bool CX_GoodXY_MASA(SignalProfile const& sp);
  bool CX_S1Range();
  
  ////All Veto cuts
  //bool CX20() { return (veto_present==true && veto_roi_lsv_charge_vec->size()>0); }
  //bool CX21() { return veto_roi_lsv_charge_vec->at(0)<=1; }
  
  
  
  //Histograms
  TH1F* h_tdrift;
  TH2F* h_s2_ch_occ;
  TH2F* h_xy_sing;
  TH2F* h_xy_comb;
  TH2F* h_xy_sing_chi2cx;
  TH2F* h_xy_comb_chi2cx;
  TH2F* h_sladxy_sing_uar;
  TH2F* h_jasonxy_sing;
  
  TH2F* h_xy[n_chi2cx];
  TH1F* h_p0_start;
  TH2F* h_chi2_dr;
  TH2F* h_logchi2_dr;
  TH2F* h_TBchi2_dr;
  TH2F* h_rerr_dr;
  TH2F* h_rerr_r_comb;
  TH2F* h_rerr_r_sing;
  TH2F* h_maxerr_dr;
  TH2F* h_maxerr_r_comb;
  TH2F* h_maxerr_r_sing;
  
  TH2F* h_top_var_dr;
  TH2F* h_top_var_r_sing;
  TH2F* h_top_var_r_comb;
  TH2F* h_top_var_r[n_chi2cx];
  TH2F* h_bot_var_dr;
  TH2F* h_bot_var_r_sing;
  TH2F* h_bot_var_r_comb;
  TH2F* h_bot_var_r[n_chi2cx];
  TH2F* h_topbot_var_dr;
  TH2F* h_topbot_var_r_sing;
  TH2F* h_topbot_var_r_comb;
  TH2F* h_topbot_var_r[n_chi2cx];
  TH1F* h_chi2_comb;
  TH1F* h_chi2_sing;
  TH1F* h_sladchi2_sing_uar;
  TH1F* h_TBchi2_comb;
  TH1F* h_TBchi2_sing;
  TH1F* h_TB_comb;
  TH1F* h_TB_sing;
  TH1F* h_chi2[2];
  
  TH1F* h_rchi2_sing;
  TH1F* h_rchi2_sing_aar;
  TH1F* h_rchi2_comb;
  TH2F* h_rchi2_chi2_sing;
  TH2F* h_rchi2_chi2_sing_aar;
  TH2F* h_rchi2_chi2_comb;
  TH2F* h_fracchi2_chi2_sing;
  TH2F* h_fracchi2_chi2_sing_aar;
  TH2F* h_fracchi2_chi2_comb;
  TH2F* h_fracchi2_logchi2_sing;
  TH2F* h_fracchi2_logchi2_sing_aar;
  TH2F* h_fracchi2_logchi2_comb;
  
  TH2F* h_maxchi2_chi2_sing;
  TH2F* h_maxchi2_chi2_sing_aar;
  TH2F* h_maxchi2_chi2_comb;
  
  TH1F* h_s1_comb;
  TH1F* h_s1_comb_goodxy;
  TH1F* h_s1_comb_badxy;
  TH1F* h_s1_sing;
  TH1F* h_s1_sing_goodxy;
  TH1F* h_s1_sing_badxy;
  
  TH2F* h_dexpxy_dr;
  TH1F* h_R;
  TH1F* h_R_masas;
  TH1F* h_R_jasons;
  TH1F* h_R_s1cx;
  TH1F* h_R_masas_s1cx;
  TH1F* h_R_jasons_s1cx;
  
  //aar
  TH1F* h_chi2_sing_aar;
  TH1F* h_sladchi2_sing_aar;
  TH2F* h_xy_sing_aar;
  TH2F* h_xy_sing_aar_chi2cx;
  TH2F* h_sladxy_sing_aar;
  TH2F* h_jasonxy_sing_aar;
  TH2F* h_top_var_r_sing_aar;
  TH2F* h_bot_var_r_sing_aar;
  TH2F* h_topbot_var_r_sing_aar;
  TH2F* h_rerr_r_sing_aar;
  TH2F* h_maxerr_r_sing_aar;
  TH1F* h_TBchi2_sing_aar;
  TH1F* h_TB_sing_aar;

  TH2F* h_maxvar_chi2_sing;
  TH2F* h_maxvar_chi2_sing_aar;
  TH2F* h_maxvar_chi2_comb;
  TH2F* h_rvar_chi2_sing;
  TH2F* h_rvar_chi2_sing_aar;
  TH2F* h_rvar_chi2_comb;  


  //AAr UAr after experimental f90
  TH2F* h_f90_s1_uar;
  
  TH2F* h_f90_s1_chi2cx_uar[n_chi2cx];
  TH1F* h_f90_uar[n_slice];
  TH1F* h_f90_uar_chi2cx[n_chi2cx][n_slice];
  
  TH2F* h_s2s1_vs_s1[2]; // 1 for not corr, 1 for corr
  TH2F* h_xy_comb_ev[N_Evt];   TH2F* h_xy_sing1_ev[N_Evt];  TH2F* h_xy_sing2_ev[N_Evt];
  
  
  ClassDef(XY4SLADSelector,0);
};

#endif


#ifdef XY4SLADSelector_cxx
//------------------------------------------------------------------------------
// Define cuts. These definitions are taken from
bool XY4SLADSelector::CX_UAR() {return (run_id >= 11856 && run_id <= 13181);}
bool XY4SLADSelector::CX_AAR() {
  //return (run_id >= 5372 && run_id <= 8433);
  return (run_id >= 7950 && run_id <= 7970);
}
bool XY4SLADSelector::CX1() { return nchannels==38 ; }
bool XY4SLADSelector::CX2() { return baseline_not_found == false; }
bool XY4SLADSelector::CX3() { return (live_time+inhibit_time)>=1.35e-3; }
bool XY4SLADSelector::CX4() { return live_time < 1.; }
bool XY4SLADSelector::CX_2PLS() { return npulses==2 || (npulses==3 && has_s3); }
bool XY4SLADSelector::CX_3PLS() { return npulses==3 && !has_s3;}
bool XY4SLADSelector::CX_TDrift_STD() { return tdrift>40. && tdrift<336.; }
bool XY4SLADSelector::CX_TDrift() {
  if (npulses<1) return true;
  else if(CX_3PLS()) return tdrift>10. && pulse_start_time[2]-pulse_start_time[0]<366;
  else if(CX_2PLS()) return tdrift>10. && tdrift<336.;
  else return tdrift>10. && tdrift<336.;
}

bool XY4SLADSelector::CX_TRIG() { return
  ((run_id >= -999 && run_id < 7344   && s1_start_time >= -0.25 && s1_start_time <= -0.15) ||
   (run_id >= 7344 && run_id < 7641   && s1_start_time >= -4.10 && s1_start_time <= -4.00) ||
   (run_id >= 7641 && run_id < 999999 && s1_start_time >= -6.10 && s1_start_time <= -6.00)); }

bool XY4SLADSelector::CX_Sat() {
  if (npulses<1) return true;
  if (npulses==1) return pulse_saturated[0]==0;
  bool notSat = pulse_saturated[0]==0 && pulse_saturated[1]==0;
  if (XY4SLADSelector::CX_3PLS())
    return notSat && pulse_saturated[2]==0;
  else
    return notSat;
}

bool XY4SLADSelector::CX_S2Sat() {
  if (npulses<=1) return true;
  else if (XY4SLADSelector::CX_3PLS())
    return pulse_saturated[1]==0 && pulse_saturated[2]==0;
  else
    return pulse_saturated[1]==0;
}

bool XY4SLADSelector::CX_C99S1MF() { return s1_max_frac < c99_s1mf_thresholds; }
bool XY4SLADSelector::CX_C95S1MF() { return s1_max_frac < c95_s1mf_thresholds; }
bool XY4SLADSelector::CX_S2True() { //combine s2_f90 cut and s2_size cut
  if (npulses<=1) return false;
  bool isS2 = (pulse_f90[1] < 0.2 && pulse_total_npe[1] > 100.);
  if (npulses==3 && !has_s3)
    return (isS2 && pulse_f90[2] < 0.2 && pulse_total_npe[2] > 100.);
  else
    return isS2;
}
//cut on corr S2
bool XY4SLADSelector::CX_S2CorrTrue() { //combine s2_f90 cut and s2_size cut
  if (npulses<=1) return false;
  bool isS2 = (pulse_f90[1] < 0.2 && pulse_total_npe[1]* allpulses_xycorr_factor[1]> 100.);
  if (npulses==3 && !has_s3)
    return (isS2 && pulse_f90[2] < 0.2 && pulse_total_npe[2]* allpulses_xycorr_factor[2]> 100.);
  else
    return isS2;
}
bool XY4SLADSelector::CX_S1F90() { return total_f90 >= 0.05; }
bool XY4SLADSelector::CX_c99S1MF()    { return max_s1_frac_cut_exceeds99; }
bool XY4SLADSelector::CX_S2dZ()   {return s2_dz>def_s2_dz;}

//XY Fid CUT
bool XY4SLADSelector::CX_R() {return (x[1]*x[1]+y[1]*y[1]) <200;}
bool XY4SLADSelector::CX_Surface() {return ((x[1]*x[1]+y[1]*y[1]) >200 && (x[1]*x[1]+y[1]*y[1]) <400);}
bool XY4SLADSelector::CX_GoodXY() {
  if (npulses<2) return false;
  else if (CX_3PLS()) return (x[1]*x[1]+y[1]*y[1]) <400 && (x[2]*x[2]+y[2]*y[2]) <400;
  else if (CX_2PLS() || npulses>3) return (x[1]*x[1]+y[1]*y[1]) <400;
  else return false;
}

bool XY4SLADSelector::CX_GoodXY_MASA(SignalProfile const& sp) {
  return (sp.PosiX*sp.PosiX+sp.PosiY*sp.PosiY) <400;
}

bool XY4SLADSelector::CX_S1Range(){return (total_s1>20 && total_s1<450);}

////All Veto cuts
//bool XY4SLADSelector::CX20() { return (veto_present==true && veto_roi_lsv_charge_vec->size()>0); }
//bool XY4SLADSelector::CX21() { return veto_roi_lsv_charge_vec->at(0)<=1; }



void XY4SLADSelector::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  //   fChain->SetMakeClass(1);
  
  fChain->SetBranchStatus("*",0); //disable all
  
  //Main file
  fChain->SetBranchStatus("events.run_id", 1);
  fChain->SetBranchAddress("run_id", &run_id);
  
  fChain->SetBranchStatus("events.event_id", 1);
  fChain->SetBranchAddress("event_id", &event_id);
  
  fChain->SetBranchStatus("events.hasV1724",1);
  fChain->SetBranchAddress("hasV1724", &hasV1724);
  
  fChain->SetBranchStatus("nchannel.nchannel", 1);
  fChain->SetBranchAddress("nchannel.nchannel", &nchannels);
  
  fChain->SetBranchStatus("baseline.SumChannelHasNoBaseline",1);
  fChain->SetBranchAddress("baseline.SumChannelHasNoBaseline", & baseline_not_found);
  
  fChain->SetBranchStatus("long_lifetime.lifetime",1);
  fChain->SetBranchAddress("long_lifetime.lifetime", & live_time);
  
  fChain->SetBranchStatus("long_lifetime.inhibittime",1);
  fChain->SetBranchAddress("long_lifetime.inhibittime", & inhibit_time);
  
  fChain->SetBranchStatus("npulses.n_phys_pulses",1);
  fChain->SetBranchAddress("npulses.n_phys_pulses", & npulses);
  
  fChain->SetBranchStatus("npulses.has_s3",1);
  fChain->SetBranchAddress("npulses.has_s3", & has_s3);
  
  fChain->SetBranchStatus("s1_time.s1_start_time", 1);
  fChain->SetBranchAddress("s1_time.s1_start_time", & s1_start_time);
  
  fChain->SetBranchStatus("s1.total_s1", 1);
  fChain->SetBranchAddress("s1.total_s1", & total_s1);
  
  fChain->SetBranchStatus("s1_f90.total_f90", 1);
  fChain->SetBranchAddress("s1_f90.total_f90", & total_f90);
  
  fChain->SetBranchStatus("s1.total_s1_corr", 1);
  fChain->SetBranchAddress("s1.total_s1_corr", & total_s1_corr);
  
  fChain->SetBranchStatus("s1.total_s1_top", 1);
  fChain->SetBranchAddress("s1.total_s1_top", & total_s1_top);
  
  fChain->SetBranchStatus("s1.total_s1_bottom", 1);
  fChain->SetBranchAddress("s1.total_s1_bottom", & total_s1_bottom);
  
  fChain->SetBranchStatus("s1_fraction.s1_max_frac", 1);
  fChain->SetBranchAddress("s1_fraction.s1_max_frac", & s1_max_frac);
  
  fChain->SetBranchStatus("s1_fraction.s1_max_chan", 1);
  fChain->SetBranchAddress("s1_fraction.s1_max_chan", & s1_max_chan);
  
  fChain->SetBranchStatus("max_s1_frac_cut.max_s1_frac_cut_exceeds99", 1);
  fChain->SetBranchAddress("max_s1_frac_cut.max_s1_frac_cut_exceeds99", & max_s1_frac_cut_exceeds99);
  
  fChain->SetBranchStatus("s2.total_s2", 1);
  fChain->SetBranchAddress("s2.total_s2", & total_s2);
  
  fChain->SetBranchStatus("s2.total_s2_corr", 1);
  fChain->SetBranchAddress("s2.total_s2_corr", & total_s2_corr);
  
  fChain->SetBranchStatus("s2_fraction.s2_chan", 1);
  fChain->SetBranchAddress("s2_fraction.s2_chan", s2_chan);
  
  fChain->SetBranchStatus("s2.total_s2_top", 1);
  fChain->SetBranchAddress("s2.total_s2_top", & total_s2_top);
  
  fChain->SetBranchStatus("s2.total_s2_bottom", 1);
  fChain->SetBranchAddress("s2.total_s2_bottom", & total_s2_bottom);
  
  fChain->SetBranchStatus("tdrift.tdrift", 1);
  fChain->SetBranchAddress("tdrift.tdrift", & tdrift);
  
  
  //masas xy file
  fChain->SetBranchStatus("masas_xy.masas_x", 1);
  fChain->SetBranchAddress("masas_xy.masas_x", &masas_x);
  
  fChain->SetBranchStatus("masas_xy.masas_y", 1);
  fChain->SetBranchAddress("masas_xy.masas_y", &masas_y);
  
  fChain->SetBranchStatus("masas_xy.masas_y", 1);
  fChain->SetBranchAddress("masas_xy.masas_chi2", &masas_chi2);
  
  fChain->SetBranchStatus("allpulses_xy.allpulses_x", 1);
  fChain->SetBranchAddress("allpulses_xy.allpulses_x", allpulses_x);
  
  fChain->SetBranchStatus("allpulses_xy.allpulses_y", 1);
  fChain->SetBranchAddress("allpulses_xy.allpulses_y", allpulses_y);
  
  fChain->SetBranchStatus("allpulses_xy.allpulses_chi2", 1);
  fChain->SetBranchAddress("allpulses_xy.allpulses_chi2", allpulses_chi2);
  
  fChain->SetBranchStatus("allpulses_xy.allpulses_xycorr_factor", 1);
  fChain->SetBranchAddress("allpulses_xy.allpulses_xycorr_factor", allpulses_xycorr_factor);
  
  
  //jason's xy file
  fChain->SetBranchStatus("xylocator_xy.xyl_best_x", 1);
  fChain->SetBranchAddress("xylocator_xy.xyl_best_x", &xyl_best_x);
  
  fChain->SetBranchStatus("xylocator_xy.xyl_best_y", 1);
  fChain->SetBranchAddress("xylocator_xy.xyl_best_y", &xyl_best_y);
  
  fChain->SetBranchStatus("xylocator_xy.xyl_best_chi2", 1);
  fChain->SetBranchAddress("xylocator_xy.xyl_best_chi2", &xyl_best_chi2);
  
  fChain->SetBranchStatus("allpulses_xyl_xy.allpulses_xyl_x", 1);
  fChain->SetBranchAddress("allpulses_xyl_xy.allpulses_xyl_x", allpulses_xyl_x);
  
  fChain->SetBranchStatus("allpulses_xyl_xy.allpulses_xyl_y", 1);
  fChain->SetBranchAddress("allpulses_xyl_xy.allpulses_xyl_y", allpulses_xyl_y);
  
  fChain->SetBranchStatus("allpulses_xyl_xy.allpulses_xyl_chi2", 1);
  fChain->SetBranchAddress("allpulses_xyl_xy.allpulses_xyl_chi2", allpulses_xyl_chi2);
  
  fChain->SetBranchStatus("allpulses_xyl_xy.allpulses_xyl_xycorr_factor", 1);
  fChain->SetBranchAddress("allpulses_xyl_xy.allpulses_xyl_xycorr_factor", allpulses_xyl_xycorr_factor);
  
  
  // all pulse file
  fChain->SetBranchStatus("pulse_info.pulse_info_npulses", 1);
  fChain->SetBranchAddress("pulse_info.pulse_info_npulses", &pulse_npulses);
  
  fChain->SetBranchStatus("pulse_info.pulse_info_start_time", 1);
  fChain->SetBranchAddress("pulse_info.pulse_info_start_time", pulse_start_time);
  
  fChain->SetBranchStatus("pulse_info.pulse_info_total_npe", 1);
  fChain->SetBranchAddress("pulse_info.pulse_info_total_npe", pulse_total_npe);
  
  fChain->SetBranchStatus("pulse_info.pulse_info_max_chan", 1);
  fChain->SetBranchAddress("pulse_info.pulse_info_max_chan", pulse_max_chan);
  
  
  fChain->SetBranchStatus("pulse_info.pulse_info_fixed_int1", 1);
  fChain->SetBranchAddress("pulse_info.pulse_info_fixed_int1", pulse_fixed_int1);
  
  fChain->SetBranchStatus("pulse_info.pulse_info_fixed_int2", 1);
  fChain->SetBranchAddress("pulse_info.pulse_info_fixed_int2", pulse_fixed_int2);
  
  fChain->SetBranchStatus("pulse_info.pulse_info_f90", 1);
  fChain->SetBranchAddress("pulse_info.pulse_info_f90", pulse_f90);
  
  fChain->SetBranchStatus("pulse_info.pulse_info_saturated", 1);
  fChain->SetBranchAddress("pulse_info.pulse_info_saturated", pulse_saturated);
  
  fChain->SetBranchStatus("pulse_info.pulse_info_ch_light_size", 1);
  fChain->SetBranchAddress("pulse_info.pulse_info_ch_light_size", &pulse_ch_light_size);
  
  fChain->SetBranchStatus("pulse_info.pulse_info_ch_light", 1);
  fChain->SetBranchAddress("pulse_info.pulse_info_ch_light", pulse_ch_light);
  
  
  
}




Bool_t XY4SLADSelector::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}


//====================================
// The rest are helper functions
// They are called to compute variables
//====================================

Float_t XY4SLADSelector::GetS2LogTopBot(SignalProfile const& sp){
  Float_t s2_bot=0; Float_t s2_top=0;
  for (int i=0; i<38; i++) {
    if (i<19) s2_bot += sp.s2[i];
    else s2_top += sp.s2[i];
  }
  return log10(s2_top/s2_bot);
}

Double_t XY4SLADSelector::s1_corr_factor(Double_t t_drift_max, Double_t t_drift)
{
  Double_t z = t_drift/(0.5*t_drift_max); // note normalization is to 0.5*t_drift_max
                                          // looked at Kr peaks in 15us t_drift windows (Run5330+5340), and fit these to [0]*z^5 + [1]*z^4 + [2]*z^3+[3]*z^2+[4]*z+[5].
  Double_t fit_par0 = 0.0407;
  Double_t fit_par1 = -0.206;
  Double_t fit_par2 = 0.407;
  Double_t fit_par3 = -0.389;
  Double_t fit_par4 = 0.247;
  Double_t fit_par5 = 0.898;
  // normalizing all points on fitted curve to expected Kr peak at t_drift_max/2
  Double_t exp_Kr_peak_at_half_t_drift_max = fit_par0 + fit_par1 + fit_par2 + fit_par3 + fit_par4 + fit_par5;
  Double_t exp_Kr_peak_at_t_drift = fit_par0*z*z*z*z*z + fit_par1*z*z*z*z + fit_par2*z*z*z + fit_par3*z*z + fit_par4*z + fit_par5;
  return exp_Kr_peak_at_half_t_drift_max/exp_Kr_peak_at_t_drift; // s1 correction factor
}

Double_t XY4SLADSelector::s1_TBAcorr_factor(Double_t s1_top, Double_t s1_bottom, Double_t s1)
{
  //TBAsymmetry correction
  Double_t par2[] = {-0.0397956, -0.27216, 0.794036, 1.70427, -3.98323, -8.50783, -2.66051};
  Double_t x = (s1_top-s1_bottom)/ s1;// TBAsym;
  Double_t diff_total_s1 = par2[0]+(par2[1]+(par2[2]+(par2[3]+(par2[4]+par2[5]*x)*x)*x)*x)*x; //(total_s1-total_s1_corr)/total_s1
                                                                                              //   Double_t total_s1_TBAcorr =  s1*(1.-diff_total_s1);
  return (1.-diff_total_s1);
}

//Calculate Top PMT Variance from XY light distribution
// by default, the frist s2 is used
// pls_index=2 for the second s2
// pls_index=1 for the first s2
Float_t XY4SLADSelector::GetTopPMTVariance(SignalProfile sigProf){
  Double_t bary_x(0.), bary_y(0.), bary_x_sq(0.), bary_y_sq(0.);
  Double_t top_s2(0.), top_s2_sq(0.);
  for (Int_t ch = PMTGeom::N_CHANNELS/2; ch < PMTGeom::N_CHANNELS; ch++) {
    Double_t s2 = sigProf.s2[ch];
    
    top_s2 += s2;
    top_s2_sq += s2*s2;
    bary_x += PMTGeom::pmtUnit*PMTGeom::pmt_x[ch] * s2;
    bary_y += PMTGeom::pmtUnit*PMTGeom::pmt_y[ch] * s2;
    bary_x_sq += PMTGeom::pmtUnit*PMTGeom::pmtUnit*PMTGeom::pmt_x[ch]*PMTGeom::pmt_x[ch] * s2;
    bary_y_sq += PMTGeom::pmtUnit*PMTGeom::pmtUnit*PMTGeom::pmt_y[ch]*PMTGeom::pmt_y[ch] * s2;
  }// end loop over channels
  
  //calculate unbiased weighted variance of xy
  Double_t denom = top_s2*top_s2 - top_s2_sq;
  Double_t variance_x = (top_s2*bary_x_sq - bary_x*bary_x)/denom;
  Double_t variance_y = (top_s2*bary_y_sq - bary_y*bary_y)/denom;
  
  return variance_x + variance_y;
}


Float_t XY4SLADSelector::GetPMT1DMaxVariance(SignalProfile sigProf){

  int N_I=90; Double_t d_theta = TMath::Pi()/N_I; 
  Double_t max_var(0.);
  for (int i=0; i<N_I; i++){
    Double_t theta = i*d_theta; 
    Double_t bary_x(0.), bary_y(0.), bary_x_sq(0.), bary_y_sq(0.); 
    Double_t top_s2(0.), top_s2_sq(0.);

   for (Int_t ch = PMTGeom::N_CHANNELS/2; ch < PMTGeom::N_CHANNELS; ch++) {
      Double_t s2 = sigProf.s2[ch];
      top_s2 += s2;
      top_s2_sq += s2*s2;


      Double_t rot_x = PMTGeom::pmtUnit*(PMTGeom::pmt_x[ch]*TMath::Cos(theta)-PMTGeom::pmt_y[ch]*TMath::Sin(theta));
      bary_x += rot_x * s2;
      bary_x_sq += rot_x*rot_x * s2;

    }// end loop over channels    
  
     Double_t denom = top_s2*top_s2 - top_s2_sq;
     Double_t variance_x = (top_s2*bary_x_sq - bary_x*bary_x)/denom;
     if(variance_x>max_var)
       max_var = variance_x;
  }
  return max_var;
  
}

Float_t XY4SLADSelector::GetPMT1DMaxVariance(int index){

  int N_I=90; Double_t d_theta = TMath::Pi()/N_I;
  Double_t max_var(0.);
  for (int i=0; i<N_I; i++){
    Double_t theta = i*d_theta;
    Double_t bary_x(0.), bary_y(0.), bary_x_sq(0.), bary_y_sq(0.);
    Double_t top_s2(0.), top_s2_sq(0.);

    for (Int_t ch = PMTGeom::N_CHANNELS/2; ch < PMTGeom::N_CHANNELS; ch++) {
      Double_t s2 = s2_chan[ch];
      int pls_3rd = 2*38+ch;
      if (index ==2) s2=pulse_ch_light[pls_3rd];
      top_s2 += s2;
      top_s2_sq += s2*s2;

      Double_t rot_x = PMTGeom::pmtUnit*(PMTGeom::pmt_x[ch]*TMath::Cos(theta)-PMTGeom::pmt_y[ch]*TMath::Sin(theta));
      bary_x += rot_x * s2;
      bary_x_sq += rot_x*rot_x * s2;

    }// end loop over channels                                                                                                                        
    Double_t denom = top_s2*top_s2 - top_s2_sq;
    Double_t variance_x = (top_s2*bary_x_sq - bary_x*bary_x)/denom;
    if(variance_x>max_var)
      max_var = variance_x;
  }
  return max_var;

}


//the first element store max var and the second store min var
std::vector<float> XY4SLADSelector::GetPMTMaxMinVariance(SignalProfile sigProf){

  int N_I=90; Double_t d_theta = TMath::Pi()/N_I;
  Double_t min_var(9999.); 
  Double_t max_var(0.);
  for (int i=0; i<N_I; i++){
    Double_t theta = i*d_theta;
    Double_t bary_x(0.), bary_y(0.), bary_x_sq(0.), bary_y_sq(0.);
    Double_t top_s2(0.), top_s2_sq(0.);

    for (Int_t ch = PMTGeom::N_CHANNELS/2; ch < PMTGeom::N_CHANNELS; ch++) {
      Double_t s2 = sigProf.s2[ch];
      top_s2 += s2;
      top_s2_sq += s2*s2;


      Double_t rot_x = PMTGeom::pmtUnit*(PMTGeom::pmt_x[ch]*TMath::Cos(theta)-PMTGeom::pmt_y[ch]*TMath::Sin(theta));
      bary_x += rot_x * s2;
      bary_x_sq += rot_x*rot_x * s2;

    }// end loop over channels                                                                                                                                                                                                                                                                                    

    Double_t denom = top_s2*top_s2 - top_s2_sq;
    Double_t variance_x = (top_s2*bary_x_sq - bary_x*bary_x)/denom;
    if(variance_x>max_var) max_var = variance_x;
    if(variance_x<min_var) min_var = variance_x;
  }

  std::vector<float> var;
  var.push_back(max_var);
  var.push_back(min_var);
  return var;

}


Float_t XY4SLADSelector::GetTopPMTVariance(int pls_index){
  Double_t bary_x(0.), bary_y(0.), bary_x_sq(0.), bary_y_sq(0.);
  Double_t top_s2(0.), top_s2_sq(0.);
  for (Int_t ch = PMTGeom::N_CHANNELS/2; ch < PMTGeom::N_CHANNELS; ch++) {
    Double_t s2 = s2_chan[ch];
    int pls_3rd = 2*38+ch;
    if (pls_index==2) s2=pulse_ch_light[pls_3rd];
    
    
    top_s2 += s2;
    top_s2_sq += s2*s2;
    bary_x += PMTGeom::pmtUnit*PMTGeom::pmt_x[ch] * s2;
    bary_y += PMTGeom::pmtUnit*PMTGeom::pmt_y[ch] * s2;
    bary_x_sq += PMTGeom::pmtUnit*PMTGeom::pmtUnit*PMTGeom::pmt_x[ch]*PMTGeom::pmt_x[ch] * s2;
    bary_y_sq += PMTGeom::pmtUnit*PMTGeom::pmtUnit*PMTGeom::pmt_y[ch]*PMTGeom::pmt_y[ch] * s2;
  }// end loop over channels
  
  //calculate unbiased weighted variance of xy
  Double_t denom = top_s2*top_s2 - top_s2_sq;
  Double_t variance_x = (top_s2*bary_x_sq - bary_x*bary_x)/denom;
  Double_t variance_y = (top_s2*bary_y_sq - bary_y*bary_y)/denom;
  
  
  return variance_x + variance_y;
}


//Calculate Bottom PMT Variance from XY light distribution
// by default, the frist s2 is used
// pls_index=2 for the second s2
// pls_index=1 for the first s2
Float_t XY4SLADSelector::GetBotPMTVariance(SignalProfile sigProf){
  Double_t bary_x(0.), bary_y(0.), bary_x_sq(0.), bary_y_sq(0.);
  Double_t top_s2(0.), top_s2_sq(0.);
  for (Int_t ch = 0; ch < PMTGeom::N_CHANNELS/2; ch++) {
    Double_t s2 = sigProf.s2[ch];
    
    top_s2 += s2;
    top_s2_sq += s2*s2;
    bary_x += PMTGeom::pmtUnit*PMTGeom::pmt_x[ch] * s2;
    bary_y += PMTGeom::pmtUnit*PMTGeom::pmt_y[ch] * s2;
    bary_x_sq += PMTGeom::pmtUnit*PMTGeom::pmtUnit*PMTGeom::pmt_x[ch]*PMTGeom::pmt_x[ch] * s2;
    bary_y_sq += PMTGeom::pmtUnit*PMTGeom::pmtUnit*PMTGeom::pmt_y[ch]*PMTGeom::pmt_y[ch] * s2;
  }// end loop over channels
  
  //calculate unbiased weighted variance of xy
  Double_t denom = top_s2*top_s2 - top_s2_sq;
  Double_t variance_x = (top_s2*bary_x_sq - bary_x*bary_x)/denom;
  Double_t variance_y = (top_s2*bary_y_sq - bary_y*bary_y)/denom;
  
  
  return variance_x + variance_y;
}

Float_t XY4SLADSelector::GetBotPMTVariance(int pls_index){
  Double_t bary_x(0.), bary_y(0.), bary_x_sq(0.), bary_y_sq(0.);
  Double_t top_s2(0.), top_s2_sq(0.);
  for (Int_t ch = 0; ch < PMTGeom::N_CHANNELS/2; ch++) {
    Double_t s2 = s2_chan[ch];
    int pls_3rd = 2*38+ch;
    if (pls_index==2) s2=pulse_ch_light[pls_3rd];
    
    
    top_s2 += s2;
    top_s2_sq += s2*s2;
    bary_x += PMTGeom::pmtUnit*PMTGeom::pmt_x[ch] * s2;
    bary_y += PMTGeom::pmtUnit*PMTGeom::pmt_y[ch] * s2;
    bary_x_sq += PMTGeom::pmtUnit*PMTGeom::pmtUnit*PMTGeom::pmt_x[ch]*PMTGeom::pmt_x[ch] * s2;
    bary_y_sq += PMTGeom::pmtUnit*PMTGeom::pmtUnit*PMTGeom::pmt_y[ch]*PMTGeom::pmt_y[ch] * s2;
  }// end loop over channels
  
  //calculate unbiased weighted variance of xy
  Double_t denom = top_s2*top_s2 - top_s2_sq;
  Double_t variance_x = (top_s2*bary_x_sq - bary_x*bary_x)/denom;
  Double_t variance_y = (top_s2*bary_y_sq - bary_y*bary_y)/denom;
  
  return variance_x + variance_y;
}




#endif // #ifdef XY4SLADSelector_cxx
