/*

  Definition of the SignalProfile struct

  Masayuki Wada 2013-11-14

 */


#ifndef SIGNALPROFILE_h
#define SIGNALPROFILE_h

const Int_t N_CHAN = 38;

// Holds signal information for each event
struct SignalProfile
{
  // Constructor
  SignalProfile() :
    run_id(-1), event_id(-1), nPulses(-1), nChannels(-1), s1_start_time(-999.), s2_start_time(-999.),
    s1_end_time(-999.), s2_end_time(-999.), t_drift(-999.),
    total_s1(-999.), total_s1_corr(-999.),
    total_f90(-999.), total_f20(-999.), total_s2_f90(-999.), total_s2(-999.), total_s2_corr(-999.),
    s2_over_s1(-999.), s2_over_s1_corr(-999.),
    Eee(-999.), LightYield(-999.),
    baryCentX(-999.), baryCentY(-999.), varianceXY(-999.),
    modified_baryCentX(-999.), modified_baryCentY(-999.), modified_varianceXY(-999.),
    PosiX(-999.), PosiY(-999.), PosiXErr(-999.), PosiYErr(-999.), chi2XY(9999.),
    maxchi2(-999), maxchi2_x(-999), maxchi2_y(-999),
    max_s1(-999.), max_s2(-999.), max_s1_chan(-999), max_s2_chan(-999)
  {
  }

  Int_t run_id;
  Int_t event_id;
  Int_t nPulses;
  Int_t nChannels;

  Double_t s1_start_time;
  Double_t s2_start_time;
  Double_t s1_end_time;
  Double_t s2_end_time;
  Double_t t_drift;

  Double_t total_s1;
  Double_t total_s1_corr; //Correction for position-dependence goes here
  Double_t total_f90; //total_s1_90/total_s1;
  Double_t total_f20; //total_s1_20/total_s1;
  Double_t total_s2_f90; //total_s2_90/total_s2;

  Double_t total_s2; //fEvent->pulses[pulseId[1]].param.fixed_int2;
  Double_t total_s2_corr; //sigProf.total_s2/TMath::Exp(-sigProf.t_drift/electron_lifetime);
  Double_t s2_bottom; //s2 from bottom only
  Double_t s2_bottom_corr; //

  Double_t s2_over_s1; //sigProf.total_s2/sigProf.total_s1;
  Double_t s2_over_s1_corr; //sigProf.total_s2_corr/sigProf.total_s1_corr;

  Double_t Eee; //recoiling energy estimated from S1 and S2
  Double_t LightYield; //S1/Eee

  Double_t baryCentX; //x position of bary center
  Double_t baryCentY; //y position of bary center
  Double_t varianceXY; //unbiased weighted variance of xy

  Double_t modified_baryCentX; //x position of bary center weighted differently
  Double_t modified_baryCentY; //y position of bary center weighted differently
  Double_t modified_varianceXY; //unbiased weighted variance of xy weighted differently

  Double_t PosiX; //x position of reconstructed point
  Double_t PosiY; //y position of reconstructed point
  Double_t PosiXErr; //error of x position of reconstructed point
  Double_t PosiYErr; //error of y position of reconstructed point
  Double_t chi2XY; //expected distance reached from the minimum for position reconstruction
  Double_t maxchi2; //max chi2
  Double_t maxchi2_x; //x position of max chi2
  Double_t maxchi2_y; //y position of max chi2


  Double_t max_s1;
  Double_t max_s2;
  Int_t max_s1_chan;
  Int_t max_s2_chan;

  Double_t s1[N_CHAN]; //[N_CHAN]
  Double_t s2[N_CHAN]; //[N_CHAN]

  inline void clear(){
    run_id = -1; event_id = -1; nPulses = -1; nChannels = -1; s1_start_time = -999.; s2_start_time = -999.;
    s1_end_time = -999.; s2_end_time = -999.; t_drift = -999.;
    total_s1 = -999.; total_s1_corr = -999.;
    total_f90 = -999.; total_f20 = -999.; total_s2_f90 = -999.; total_s2 = -999.; total_s2_corr = -999.;
    s2_over_s1 = -999.; s2_over_s1_corr = -999.;
    Eee = -999.; LightYield = -999.;
    baryCentX = -999.; baryCentY = -999.; varianceXY = -999.;
    modified_baryCentX = -999.; modified_baryCentY = -999.; modified_varianceXY = -999.;
    PosiX = -999.; PosiY = -999.; PosiXErr = -999.; PosiYErr = -999.; chi2XY = 9999.;
    maxchi2 = -999; maxchi2_x=-999; maxchi2_y = -999;
    max_s1 = -999.; max_s2 = -999.; max_s1_chan = -999; max_s2_chan = -999;
    for(Int_t ch =0; ch<N_CHAN; ch++){
        s1[ch] = -999.; s2[ch] = -999.;
    }
  }

};

#endif
