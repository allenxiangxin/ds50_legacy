/*
 X. Xiang <xxiang@princeton.edu>
 
 Overview:
 Macro delicated to study XY recon for multi-scatters
 
 Note:
 0) I/O are specified in main (Friend files in the same dir)
 1) To run it, make sure:
    a) source correct root (try v5.34_enable_all)
    b) LRF is in work dir
    c) S1 Max Frac cut def is in work dir
    d) LD_LIBRARY_PATH include work dir
 2) this code works currently for SLAD v2.2

 
 Options:
  - add #define USE_AAR to include AAr data
 
 Usage:
 $ make
 $ ./XY4SLAD_MS
 
 */

#include "TROOT.h"
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

#include "TProof.h"
#include "TSystem.h"
#include "TDSet.h"
#include "TProofLog.h"
#include "TProofMgr.h"
#include "TEnv.h"

//using masa's xy recon
#include "./XY4SLADSelector/XY4SLADSelector.h"




//------------------------------------------------------------------------------
// Main event loop is contained here.
void Process(TChain* events, TString outName)
{

  TString workdir(gSystem->Getenv("PWD"));
  gEnv->SetValue("ProofLite.Sandbox", workdir);
  std::cout<<"work directory="<<workdir.Data()<<"\n";
  cout<<"Tree entries: "<<events->GetEntries()<<endl;
  
  //Create PROOF, load library, and initialize tchain in Proof framework
  TProof* pr = TProof::Open("workers=8");
  pr->Exec(Form("gSystem->Load(\"%s/xy4sladDict.so\")", workdir.Data()));
  events->SetProof();
  
  //Create selector and set path to LRF and path to s1mf cut def
  XY4SLADSelector *selector = new XY4SLADSelector();
  selector->SetInputList(new TList());
  selector->GetInputList()->Add(new TNamed("profilename",Form("%s/xy_profile_329_Iter.root", workdir.Data())));
  selector->GetInputList()->Add(new TNamed("s1mfName",Form("%s/max_s1_frac_cut_fixed_acceptance_full_stats.root", workdir.Data())));

  
  Bool_t withfriends = kTRUE;
  TDSet *dataset = new TDSet(*events, withfriends);
  dataset->Process(selector, outName.Data(), 4800, 0);
  //  dataset->Process(selector, outName.Data(), 48000000, 0);
  //TDSet(const TChain& chain, Bool_t withfriends = kTRUE)
  //dataset->Process("./XY4SLADSelector/XY4SLADSelector.C+", outName.Data(), 100, 0);

  // save log, this need to be done after Process
  TProofMgr* mgr = pr->GetManager();
  mgr->GetSessionLogs()->Save("*", "log_all-workers.txt") ;
  
}



//------------------------------------------------------------------------------
//Load UAr SLAD data
TChain* LoadSLAD(TString mainfile, TString mainfile2){
  
  // message
  std::cout<<"Loading SLAD file "<<mainfile.Data()<<"\n";
  
  // Consruct the name of the friend tree files
  TString s2file = mainfile;
  s2file.Remove(s2file.Length()-5);
  s2file+="_s2.root";
  
  TString xyfile = mainfile;
  xyfile.Remove(xyfile.Length()-5);
  xyfile+="_masas_xy.root";
  
  TString pulsefile = mainfile;
  pulsefile.Remove(pulsefile.Length()-5);
  pulsefile+="_allpulses.root";
  
  TString jasonxyfile = mainfile;
  jasonxyfile.Remove(jasonxyfile.Length()-5);
  jasonxyfile+="_xylocator_xy.root";
  
  
  // Load the TTrees
  TChain* events = new TChain("events");
  events->Add(mainfile);
  
  TChain* s2_fraction = new TChain("s2_fraction");
  s2_fraction->Add(s2file);
  
  TChain* masas_xy = new TChain("masas_xy");
  masas_xy->Add(xyfile);
  
  TChain* masas_xy_allpulses = new TChain("allpulses_xy");
  masas_xy_allpulses->Add(xyfile);
  
  TChain* jasons_xy = new TChain("xylocator_xy");
  jasons_xy->Add(jasonxyfile);
  
  TChain* jasons_xy_allpulse = new TChain("allpulses_xyl_xy");
  jasons_xy_allpulse->Add(jasonxyfile);
  
  TChain* pulse_info = new TChain("pulse_info");
  pulse_info->Add(pulsefile);
  
  
  //using AAr data is optional. Check if AAr path is specified
  if(mainfile2 != ""){
    
    // message
    std::cout<<"Loading SLAD file "<<mainfile2.Data()<<"\n";
    
    // Construct the name of the friend tree files
    s2file = mainfile2;
    s2file.Remove(s2file.Length()-5);
    s2file+="_s2.root";
    
    xyfile = mainfile2;
    xyfile.Remove(xyfile.Length()-5);
    xyfile+="_masas_xy.root";
    
    pulsefile = mainfile2;
    pulsefile.Remove(pulsefile.Length()-5);
    pulsefile+="_allpulses.root";
    
    jasonxyfile = mainfile2;
    jasonxyfile.Remove(jasonxyfile.Length()-5);
    jasonxyfile+="_xylocator_xy.root";
    
    // Load the TTrees
    events->Add(mainfile2);
    s2_fraction->Add(s2file);
    masas_xy->Add(xyfile);
    masas_xy_allpulses->Add(xyfile);
    jasons_xy->Add(jasonxyfile);
    jasons_xy_allpulse->Add(jasonxyfile);
    pulse_info->Add(pulsefile);
    
  }
  
  //Add friends
  events->AddFriend(s2_fraction);
  events->AddFriend(masas_xy);
  events->AddFriend(masas_xy_allpulses);
  events->AddFriend(jasons_xy);
  events->AddFriend(jasons_xy_allpulse);
  events->AddFriend(pulse_info);
  
  return events;
}



int main()
{
		
  std::cout << "\n==========> XY4SLAD_MS <=============" << std::endl;
  
  // Prevent canvases from being drawn.
  gROOT->SetBatch(kTRUE);
  
  TStopwatch* clock = new TStopwatch();
  clock->Start();
  
  
  // The main SLAD file containing the data we want.
  TString mainfile_uar = "/darkside/data/SLAD_v2_2_Merged/UAr_70d_SLAD_v2_2_merged_runlist_v1.root";
  TString mainfile_aar = "/darkside/data/SLAD_v2_2_Merged/AAr_50d_SLAD_v2_2_merged_no_r7131_r7146.root";
  
  
  //output file to save all histograms
  TString outName = "XY4SLAD_MS_Hist_1.root";
  
  
  
#define USE_AAR
#ifdef USE_AAR
  TChain* events = LoadSLAD(mainfile_uar, mainfile_aar);
#else
  TChain* events = LoadSLAD(mainfile_uar, "");
#endif
  
  Process(events, outName);
  std::cout << "Done!"<<" "<<clock->RealTime()<<" s."<<std::endl;
  
  return 1;
}


