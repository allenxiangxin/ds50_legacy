#define XY4SLADSelector_cxx
// The class definition in XY4SLADSelector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("XY4SLADSelector.C")
// root> T->Process("XY4SLADSelector.C","some options")
// root> T->Process("XY4SLADSelector.C+")
//

#include "XY4SLADSelector.h"

void XY4SLADSelector::Begin(TTree * tree)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  
  cout << "GetOption() = " << option.Data() << endl;
  TObjArray* optarr = option.Tokenize(":");
  TString outdir = ((TObjString*)optarr->At(0))->GetString();
  SetOutName(outdir);
  
  
}


void XY4SLADSelector::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  Info("XY4SLADSelector::SlaveBegin()","begining ...");
  
  TString option = GetOption();
  
  //Load Light Responce Function
  TString profilename(fInput->FindObject("profilename")->GetTitle());
  this->SetProfileName(profilename);
  fMyFcn = new XYReconstructor();
  fMyFcn->LoadProfile(fProfileName.Data());
  fMyFcn->SetMinimizer();
  
  //Load S1 Max Frac Cut
  TString s1mfName(fInput->FindObject("s1mfName")->GetTitle());
  this->SetS1MaxFracFileName(s1mfName);
  fs1mfFile = new TFile(fs1mfName.Data());
  if (fs1mfFile==NULL) {
    std::cout<<"WARNING: "<<fs1mfName.Data()<<"is not found!!\n";
    std::cout<<"===> DO NOT USE S1 MAX FRAC CUT <===\n";
  }
  h_c95_s1mf_thresholds = (TH2F*) fs1mfFile->Get("c95");
  h_c99_s1mf_thresholds = (TH2F*) fs1mfFile->Get("c99");
  if (h_c95_s1mf_thresholds==NULL) {
    std::cout<<"95 percentile s1 max frac cut is not found in"<<fs1mfName.Data()<<"\n";
  }
  if (h_c99_s1mf_thresholds==NULL) {
    std::cout<<"99 percentile s1 max frac cut is not found in"<<fs1mfName.Data()<<"\n";
  }
  //  fs1mf_file = new TFile("../Cuts/CXS1MF/s1mf_thresholds.root");
  //  TH2F* h_s1mf_thresholds = (TH2F*) fs1mf_file->Get("c90");
  
  
  //Declare histograms
  BookHistograms();
  

  
}


Bool_t XY4SLADSelector::Process(Long64_t n)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either XY4SLADSelector::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
  
  
  //get the entry
  Int_t chainentry = fChain->GetChainEntryNumber(n);
  if(chainentry%(fChain->GetEntries()/10000)==0) printf("Processing Entry number %ld [%ld %% of %lld]\n", (long int)chainentry, (long int)(chainentry/(fChain->GetEntries()/100)), fChain->GetEntries());

  fChain->GetTree()->GetEntry(n);
  
//  if(chainentry==0 ) {
//    Info("XY4SLADSelector::Process()","chainentry==0: booking histograms");
//    BookHistograms();
//  }
  
  
  //-------------------------------------
  // Calculate variables on the fly
  Int_t c95_s1mf_xbin   = h_c95_s1mf_thresholds->GetXaxis()->FindBin(total_s1_corr);
  Int_t c95_s1mf_ybin   = h_c95_s1mf_thresholds->GetYaxis()->FindBin(tdrift);
  c95_s1mf_thresholds = h_c95_s1mf_thresholds->GetBinContent(c95_s1mf_xbin, c95_s1mf_ybin);
  Int_t c99_s1mf_xbin   = h_c99_s1mf_thresholds->GetXaxis()->FindBin(total_s1_corr);
  Int_t c99_s1mf_ybin   = h_c99_s1mf_thresholds->GetYaxis()->FindBin(tdrift);
  c99_s1mf_thresholds = h_c99_s1mf_thresholds->GetBinContent(c99_s1mf_xbin, c99_s1mf_ybin);
  
  //all pulses
  for (int i=0; i<npulses; i++) {
    x[i] = ( allpulses_x[i] > -99 ? allpulses_x[i] : (allpulses_xyl_x[i] > -99 ? allpulses_xyl_x[i] : -998) );
    y[i] = ( allpulses_y[i] > -99 ? allpulses_y[i] : (allpulses_xyl_y[i] > -99 ? allpulses_xyl_y[i] : -998) );
    s2_xycorr[i] = allpulses_xycorr_factor[i]*pulse_total_npe[i];
  }
  
  //3 pls
  if (npulses>2) {
    s2_dz = pulse_start_time[2]-pulse_start_time[1];
    s2_dx = abs(allpulses_x[2]-allpulses_x[1]);
    s2_dy = abs(allpulses_y[2]-allpulses_y[1]);
    s2_dr = sqrt(pow(s2_dx,2) + pow(s2_dy, 2));
  }
  
  
  
  // Generate cuts.
  bool cx_uar = CX_UAR();
  bool cx_aar = CX_AAR();
  bool basic_cuts_aar = CX1() && CX2() && CX3() && CX4();
  bool basic_cuts_uar = CX1() && CX2() && CX3();
  bool cx_s1f90 = CX_S1F90();
  bool cx_2pls = CX_2PLS();
  bool cx_3pls = CX_3PLS();
  bool cx_s2true = CX_S2True();
  bool cx_s2corrtrue = CX_S2CorrTrue();
  bool cx_tdrift_std = CX_TDrift_STD();
  bool cx_tdrift = CX_TDrift();
  bool cx_goodxy = CX_GoodXY(); //radius is within TPC, use both masa and jason
  bool cx_s2dz = CX_S2dZ();
  bool cx_sat = CX_Sat();
  bool cx_trig = CX_TRIG();
  bool cx_95s1mf = CX_C95S1MF();
  bool cx_99s1mf = CX_C99S1MF();
  
  
#ifdef use_S1_Range
  bool cx_s1range = CX_S1Range();
#else
  bool cx_s1range = 1;
#endif
  
  
  //-------------------------------------------------
  //initialize Singal Profile for Masa's xy
  //I put cuts there to make it faster, since doing xy for every events waste time
  //Please make sure those cuts are more loose than the cuts you applied when filling histogram
  
  SignalProfile sigProf[3];

  if (cx_uar && basic_cuts_uar && cx_3pls && cx_trig && cx_s1f90 && cx_s2true && cx_sat && cx_s2dz) {
    for(int ch=0; ch<N_CHAN; ch++){
      int pls_3rd = 2*38+ch;
      sigProf[0].s2[ch] =s2_chan[ch]+pulse_ch_light[pls_3rd];
      //        cout<<"3rd pls ["<<ch<<"] = "<<pulse_ch_light[pls_3rd]<<"\n";
    }
    fMyFcn->ChargeCorrectforS2(sigProf[0]);
    fMyFcn->Minimization(sigProf[0]);
    
    //This need to be done after minimization
    double par[] = {sigProf[0].PosiX, sigProf[0].PosiY};
    Double_t weight = 0.;
    fMyFcn->GetChi2S2BotoverS2Top(par, TBChi2, weight );
  }

  if (cx_uar && basic_cuts_uar && cx_2pls && cx_trig && cx_s1f90 && cx_s2true && cx_sat ) {
    for(int ch=0; ch<N_CHAN; ch++){
      sigProf[1].s2[ch] =s2_chan[ch];
    }
    fMyFcn->ChargeCorrectforS2(sigProf[1]);
    fMyFcn->Minimization(sigProf[1]);
    
    double par[] = {sigProf[1].PosiX, sigProf[1].PosiY};
    Double_t weight = 0.;
    fMyFcn->GetChi2S2BotoverS2Top(par, TBChi2, weight);
  }
  
  if (cx_aar && basic_cuts_aar && cx_2pls && cx_trig && cx_s1f90 && cx_s2true && cx_sat) {
    for(int ch=0; ch<N_CHAN; ch++){
      sigProf[2].s2[ch] =s2_chan[ch];
    }
    fMyFcn->ChargeCorrectforS2(sigProf[2]);
    fMyFcn->Minimization(sigProf[2]);
    
    double par[] = {sigProf[2].PosiX, sigProf[2].PosiY};
    Double_t weight = 0.;
    fMyFcn->GetChi2S2BotoverS2Top(par, TBChi2, weight );
  }
  
  FillHistograms(sigProf);
  
  return kTRUE;
}

void XY4SLADSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  Info("XY4SLADSelector::SlaveTerminate()","closining ...");
  if (!fMyFcn) delete fMyFcn;
  
}

void XY4SLADSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  Info("XY4SLADSelector::Terminate()","closining ...");
  
  
  
  
  Int_t nhists = 0;
  TList *list = GetOutputList();
  TIter next(list);
//  if (foutfile==NULL) {
//    std::cout<<"XY4SLADSelector::Terminate(): foutfile is a NULL \n";
//    return ;
//  }
  // Open output file
  foutfile = new TFile(fOutName, "recreate");
  std::cout << "XY4SLADSelector::Begin(), Writing output to "<<foutfile->GetName()<<std::endl;
  foutfile->cd();
  
  TH1 *h;
  TObject* obj;
  while ((obj = (TObject*) next()))
  {
    if(!obj) continue;
    if (obj->InheritsFrom(TH1::Class()))
    {
      h = (TH1*) obj;
      h->Write();
      //          if (Debug)cout << "Write " << h->GetName() << " into the file" << endl;
      //          delete h;
      nhists++;
    }
  }
  
  

  
//  foutfile->Write();
//  foutfile->cd("events");
//  for (int i=0; i<N_Evt; i++) {
//    TCanvas c1(Form("c_%d", i),Form("c_%d",i), 1); c1.cd();
//    h_xy_comb_ev[i]->Draw();
//    h_xy_sing1_ev[i]->Draw("same");
//    h_xy_sing2_ev[i]->Draw("same");
//    
//    TLegend leg(0.9, 0.1, 1, 0.25);
//    leg.AddEntry(h_xy_comb_ev[i], "comb", "p");
//    leg.AddEntry(h_xy_sing1_ev[i], "1st s2", "p");
//    leg.AddEntry(h_xy_sing2_ev[i], "2nd s2", "p");
//    leg.Draw();
//    
//    c1.Write();
//  }
//  foutfile->cd();
//  foutfile->Close();
  
  Info("XY4SLADSelector::Terminate()","terminating ...");
  
}



Bool_t XY4SLADSelector::IsValidEvent(){
  return true;
}


Bool_t XY4SLADSelector::IsEventOK(){
  
  return true;
}








void XY4SLADSelector::BookHistograms()
{
  //  if(Debug)
  std::cout<<"XY4SLADSelector::BookHistograms()"<<std::endl;
  
  TList  *list = GetOutputList();

  //-------------------------------------------------
  // Define histograms
  h_tdrift          = new TH1F("h_tdrift", "tdrift; tdrift [us]; events", 400, 0, 400);
  h_s2_ch_occ       = new TH2F("h_s2_ch_occ", "; ch ID; S2_{i}/S2_{tot}", 40, 0, 40, 1000, 0, 1);
  h_p0_start        = new TH1F("h_p0_start", "; pulse start time [#mus]", 1000, -10, 10);
  h_dexpxy_dr       = new TH2F("h_dexpxy_dr", "; r (xy) separation of two S2s [cm]; Distance to expected XY [cm]", 100, 0, 20, 200, 0, 40);
  h_s2s1_vs_s1[0]   = new TH2F("h_s2s1_vs_s1_0", "2nd pls as S2 ; S1_{corr} [PE]; Log10(S2_{xycorr}/S1_{corr}) [PE]", 100, 0, 3000, 100, -0.5, 2.5);
  h_s2s1_vs_s1[1]   = new TH2F("h_s2s1_vs_s1_1", "3rd pls as S2 ; S1_{corr} [PE]; Log10(S2_{xycorr}/S1_{corr}) [PE]", 100, 0, 3000, 100, -0.5, 2.5);
  list->Add(h_tdrift); list->Add(h_s2_ch_occ);list->Add(h_p0_start);list->Add(h_dexpxy_dr);list->Add(h_s2s1_vs_s1[0]); list->Add(h_s2s1_vs_s1[1]);
  
  
  h_xy_sing         = new TH2F("h_xy_sing", "; x [cm]; y [cm]", 80, -20, 20, 80, -20, 20);
  h_xy_comb         = new TH2F("h_xy_comb", "; x [cm]; y [cm]", 80, -20, 20, 80, -20, 20);
  h_xy_sing_chi2cx     = new TH2F("h_xy_sing_chi2cx", "; x [cm]; y [cm]", 80, -20, 20, 80, -20, 20);
  h_xy_comb_chi2cx        = new TH2F("h_xy_comb_chi2cx", "; x [cm]; y [cm]", 80, -20, 20, 80, -20, 20);
  h_sladxy_sing_uar = new TH2F("h_sladxy_sing_uar", "; x [cm]; y [cm]", 80, -20, 20, 80, -20, 20);
  h_jasonxy_sing         = new TH2F("h_jasonxy_sing", "; x [cm]; y [cm]", 80, -20, 20, 80, -20, 20);
  for(int i=0; i<2;i++){
    h_xy[i] = new TH2F(Form("h_xy_%d",i+1), "; x [cm]; y [cm]", 80, -20, 20, 80, -20, 20);
    list->Add(h_xy[i]);
  }
  list->Add(h_xy_sing); list->Add(h_xy_comb);list->Add(h_sladxy_sing_uar);
  list->Add(h_xy_sing_chi2cx); list->Add(h_xy_comb_chi2cx);list->Add(h_jasonxy_sing);
  
  
  h_chi2_dr       = new TH2F("h_chi2_dr", "; XY separation of two S2s [cm]; Masa's Chi2", 100, 0, 20, 6400, 0, 3200);
  h_logchi2_dr    = new TH2F("h_logchi2_dr", ";XY separation of two S2s [cm]; Log10(Chi2)", 100, 0, 20, 480, -0.8, 4);
  h_TBchi2_dr     = new TH2F("h_TBchi2_dr", ";XY separation of two S2s [cm]; Masa's Bot/Top Chi2", 100, 0, 20, 2000, 0, 1000);
  h_rerr_dr       = new TH2F("h_rerr_dr", ";XY separation of two S2s [cm]; Error in r [cm]", 100, 0, 20, 100, 0, 20);
  h_rerr_r_comb   = new TH2F("h_rerr_r_comb", "; r [cm]; Error in r [cm]", 400, 0, 20, 400, 0, 20);
  h_rerr_r_sing   = new TH2F("h_rerr_r_sing", "; r [cm]; Error in r [cm]", 400, 0, 20, 400, 0, 20);
  h_maxerr_dr     = new TH2F("h_maxerr_dr", "; r separation of two S2s [cm]; Max(xerr, yerr) [cm]", 100, 0, 20, 100, 0, 20);
  h_maxerr_r_comb = new TH2F("h_maxerr_r_comb", "; r [cm]; Max(xerr, yerr) [cm]", 400, 0, 20, 400, 0, 20);
  h_maxerr_r_sing = new TH2F("h_maxerr_r_sing", "; r [cm]; Max(xerr, yerr) [cm]", 400, 0, 20, 400, 0, 20);
  list->Add(h_chi2_dr); list->Add(h_logchi2_dr);list->Add(h_TBchi2_dr);
  list->Add(h_rerr_dr); list->Add(h_rerr_r_comb);list->Add(h_rerr_r_sing);
  list->Add(h_maxerr_dr); list->Add(h_maxerr_r_comb);list->Add(h_maxerr_r_sing);


  h_top_var_dr      = new TH2F("h_top_var_dr", ";XY separation of two S2s [cm]; Top PMT Variance", 100, 0, 20, 1000, 0, 1000);
  h_top_var_r_sing  = new TH2F("h_top_var_r_sing", "; R [cm]; Top PMT XY Variance [cm^{2}]", 400, 0, 20, 1000, 0, 1000);
  h_top_var_r_comb  = new TH2F("h_top_var_r_comb", "; R [cm]; Top PMT XY Variance [cm^{2}]", 400, 0, 20, 1000, 0, 1000);
  h_maxvar_chi2_sing  = new TH2F("h_maxvar_chi2_sing", "; Min Chi2; Max 1D Variance [cm^{2}]", 6400, 0, 3200, 1000, 0, 1000);
  h_rvar_chi2_sing  = new TH2F("h_rvar_chi2_sing", "; Min Chi2; Max/Min Variance", 6400, 0, 3200, 1000, 0.8, 50.8);
  h_maxvar_chi2_comb  = new TH2F("h_maxvar_chi2_comb", "; Min Chi2; Max 1D Variance [cm^{2}]", 6400, 0, 3200, 1000, 0, 100);
  h_rvar_chi2_comb  = new TH2F("h_rvar_chi2_comb", "; Min Chi2; Max/Min Variance", 6400, 0, 3200, 1000, 0.8, 50.8);
  for(int i=0; i<2;i++){
    h_top_var_r[i]= new TH2F(Form("h_top_var_r_%d",i+1), "; R [cm]; Top PMT light Variance", 400, 0, 20, 1000, 0, 1000);
    list->Add(h_top_var_r[i]);
  }
  list->Add(h_top_var_dr); list->Add(h_top_var_r_sing);list->Add(h_top_var_r_comb); 
  list->Add(h_maxvar_chi2_sing); list->Add(h_maxvar_chi2_comb); list->Add(h_rvar_chi2_sing); list->Add(h_rvar_chi2_comb);
  
  
  h_bot_var_dr      = new TH2F("h_bot_var_dr", ";XY separation of two S2s [cm]; Bottom PMT XY Variance [cm^{2}]", 100, 0, 20, 1000, 0, 1000);
  h_bot_var_r_sing  = new TH2F("h_bot_var_r_sing", "; R [cm]; Bottom PMT XY Variance [cm^{2}]", 400, 0, 20, 1000, 0, 1000);
  h_bot_var_r_comb  = new TH2F("h_bot_var_r_comb", "; R [cm]; Bottom PMT XY Variance [cm^{2}]", 400, 0, 20, 1000, 0, 1000);
  for(int i=0; i<2;i++){
    h_bot_var_r[i]= new TH2F(Form("h_bot_var_r_%d",i+1), "; R [cm]; Bottom PMT XY Variance [cm^{2}]", 400, 0, 20, 1000, 0, 1000);
    list->Add(h_bot_var_r[i]);
  }
  list->Add(h_bot_var_dr); list->Add(h_bot_var_r_sing);list->Add(h_bot_var_r_comb);
  
  
  h_topbot_var_dr     = new TH2F("h_topbot_var_dr", ";XY separation of two S2s [cm]; Top Variance/Bot Variance", 100, 0, 20, 1000, 0, 1000);
  h_topbot_var_r_sing = new TH2F("h_topbot_var_r_sing", "; R [cm]; Top Variance/Bot Variance", 400, 0, 20, 1000, 0, 1000);
  h_topbot_var_r_comb = new TH2F("h_topbot_var_r_comb", "; R [cm]; Top Variance/Bot Variance", 400, 0, 20, 1000, 0, 1000);
  for(int i=0; i<2;i++){
  h_topbot_var_r[i]= new TH2F(Form("h_topbot_var_r_%d",i+1), "; R [cm]; Top Variance/Bot Variance", 400, 0, 20, 1000, 0, 1000);
  list->Add(h_topbot_var_r[i]);
  }
  list->Add(h_topbot_var_dr); list->Add(h_topbot_var_r_sing);list->Add(h_topbot_var_r_comb);
  
  
  h_chi2_comb         = new TH1F("h_chi2_comb", "UAr 3-pulse events (Combined S2); Minimum Chi2; Events ", 4000, 0, 2000);
  h_chi2_sing         = new TH1F("h_chi2_sing", "UAr 2-pulse events; Minimum Chi2; Events", 4000, 0, 2000);
  h_sladchi2_sing_uar = new TH1F("h_sladchi2_sing_uar", "Single-Scatter;  Minimum Chi2; Events", 4000, 0, 2000);
  h_TBchi2_comb       = new TH1F("h_TBchi2_comb", "Multi-scatter;  Bot/Top Chi2; Events ", 1000, 0, 1000);
  h_TBchi2_sing       = new TH1F("h_TBchi2_sing", "Single-Scatter;  Bot/Top Chi2; Events", 1000, 0, 1000);
  h_TB_comb           = new TH1F("h_TB_comb", "Multi-scatter; Masa's Bot/Top Chi2; Events ", 1000, -2, 3);
  h_TB_sing           = new TH1F("h_TB_sing", "Single-Scatter; Masa's Bot/Top Chi2; Events", 1000, -2, 3);
  for(int i=0; i<2;i++){
    h_chi2[i] = new TH1F(Form("h_chi2_%d",i+1), "; Minimum Chi2; Events", 6400, 0, 6400);
    list->Add(h_chi2[i]);
  }
  list->Add(h_chi2_comb); list->Add(h_chi2_sing);list->Add(h_sladchi2_sing_uar);
  list->Add(h_TBchi2_comb); list->Add(h_TBchi2_sing);list->Add(h_TB_comb);list->Add(h_TB_sing);


  h_rchi2_sing          = new TH1F("h_rchi2_sing", "UAr 2-pulse events; Max_Chi2/Min_Chi2; Events", 5000, 0, 100000 );
  h_rchi2_sing_aar      = new TH1F("h_rchi2_sing_aar", "AAr 2-pulse events; Max_Chi2/Min_Chi2; Events", 5000, 0, 100000 );
  h_rchi2_comb          = new TH1F("h_rchi2_comb", "UAr 3-pulse events (Combined S2); Max_Chi2/Min_Chi2; Events", 5000, 0, 100000 );
  h_rchi2_chi2_sing     = new TH2F("h_rchi2_chi2_sing", "UAr 2-pulse events; Minimum Chi2; Max_Chi2/Min_Chi2;", 4000, 0, 2000, 5000, 0, 100000 );
  h_rchi2_chi2_sing_aar = new TH2F("h_rchi2_chi2_sing_aar", "AAr 2-pulse events; Minimum Chi2; Max_Chi2/Min_Chi2;", 4000, 0, 2000, 5000, 0, 100000 );
  h_rchi2_chi2_comb     = new TH2F("h_rchi2_chi2_comb", "UAr 3-pulse events (Combined S2); Minimum Chi2; Max_Chi2/Min_Chi2;", 4000, 0, 2000, 5000, 0, 100000 );
  list->Add(h_rchi2_sing); list->Add(h_rchi2_sing_aar);list->Add(h_rchi2_comb);
  list->Add(h_rchi2_chi2_sing); list->Add(h_rchi2_chi2_sing_aar);list->Add(h_rchi2_chi2_comb);
  
  
  h_fracchi2_chi2_sing        = new TH2F("h_fracchi2_chi2_sing", "UAr 2-pulse events; Minimum Chi2; (log(max_chi2)-log(min_chi2))/log(min_chi2);", 1000, 0, 200, 2000, -1, 10 );
  h_fracchi2_chi2_sing_aar    = new TH2F("h_fracchi2_chi2_sing_aar", "AAr 2-pulse events; Minimum Chi2; (log(max_chi2)-log(min_chi2))/log(min_chi2);", 1000, 0, 200, 2000, -1, 10 );
  h_fracchi2_chi2_comb        = new TH2F("h_fracchi2_chi2_comb", "UAr 3-pulse events (Combined S2); Minimum Chi2; (log(max_chi2)-log(min_chi2))/log(min_chi2);", 1000, 0, 200, 2000, -1, 10 );
  h_fracchi2_logchi2_sing     = new TH2F("h_fracchi2_logchi2_sing", "UAr 2-pulse events; log(min_chi2); (log(max_chi2)-log(min_chi2))/log(min_chi2);", 1000, -0.2, 2.5, 1000, -0.2, 10 );
  h_fracchi2_logchi2_sing_aar = new TH2F("h_fracchi2_logchi2_sing_aar", "AAr 2-pulse events; log(min_chi2); (log(max_chi2)-log(min_chi2))/log(min_chi2);", 1000, -0.2, 2.2, 1000, -0.2, 10 );
  h_fracchi2_logchi2_comb     = new TH2F("h_fracchi2_logchi2_comb", "UAr 3-pulse events (Combined S2); log(min_chi2); (log(max_chi2)-log(min_chi2))/log(min_chi2);", 1000, -0.2, 2.5, 1000, -0.2, 10 );
  
  h_maxchi2_chi2_sing     = new TH2F("h_maxchi2_chi2_sing", "UAr 2-pulse events; Log10(Min Chi2); Log10(Max Chi2);", 1000, 0, 3.5, 1000, 0, 8 );
  h_maxchi2_chi2_sing_aar = new TH2F("h_maxchi2_chi2_sing_aar", "AAr 2-pulse events; Log10(Min Chi2); Log10(Max Chi2);", 1000, 0, 3.5, 1000, 0, 8 );
  h_maxchi2_chi2_comb     = new TH2F("h_maxchi2_chi2_comb", "UAr 3-pulse events (Combined S2); Log10(Min Chi2); Log10(Max Chi2);", 1000, 0, 3.5, 1000, 0, 8 );
  list->Add(h_fracchi2_chi2_sing); list->Add(h_fracchi2_chi2_sing_aar);list->Add(h_fracchi2_chi2_comb);
  list->Add(h_fracchi2_logchi2_sing); list->Add(h_fracchi2_logchi2_sing_aar);list->Add(h_fracchi2_logchi2_comb);
  list->Add(h_maxchi2_chi2_sing); list->Add(h_maxchi2_chi2_sing_aar);list->Add(h_maxchi2_chi2_comb);
  
  h_s1_comb         = new TH1F("h_s1_comb","h_s1_comb; S1 (no corr) [pe]", 6000, 0, 6000);
  h_s1_comb_goodxy  = new TH1F("h_s1_comb_goodxy","h_s1_comb_goodxy; S1 (no corr) [pe]", 6000, 0, 6000);
  h_s1_comb_badxy   = new TH1F("h_s1_comb_badxy","h_s1_comb_badxy; S1 (no corr) [pe]", 6000, 0, 6000);
  h_s1_sing         = new TH1F("h_s1_sing","h_s1_sing; S1 (no corr) [pe]", 6000, 0, 6000);
  h_s1_sing_goodxy  = new TH1F("h_s1_sing_goodxy","h_s1_sing_goodxy; S1 (no corr) [pe]", 6000, 0, 6000);
  h_s1_sing_badxy   = new TH1F("h_s1_sing_badxy","h_s1_sing_badxy; S1 (no corr) [pe]", 6000, 0, 6000);
  list->Add(h_s1_comb); list->Add(h_s1_comb_goodxy);list->Add(h_s1_comb_badxy);
  list->Add(h_s1_sing); list->Add(h_s1_sing_goodxy);list->Add(h_s1_sing_badxy);
  
  
  h_R             = new TH1F("h_R", "h_R; sqrt(x^2+y^2+z^2); events", 200, 0, 45);
  h_R_masas       = new TH1F("h_R_masas", "h_R; sqrt(x^2+y^2+z^2); events", 200, 0, 45);
  h_R_jasons      = new TH1F("h_R_jasons", "h_R; sqrt(x^2+y^2+z^2); events", 200, 0, 45);
  h_R_s1cx        = new TH1F("h_R_s1cx", "h_R; sqrt(x^2+y^2+z^2); events", 200, 0, 45);
  h_R_masas_s1cx  = new TH1F("h_R_masas_s1cx", "h_R; sqrt(x^2+y^2+z^2); events", 200, 0, 45);
  h_R_jasons_s1cx = new TH1F("h_R_jasons_s1cx", "h_R; sqrt(x^2+y^2+z^2); events", 200, 0, 45);
  list->Add(h_R); list->Add(h_R_masas);list->Add(h_R_jasons);
  list->Add(h_R_s1cx); list->Add(h_R_masas_s1cx);list->Add(h_R_jasons_s1cx);
  

  //aar
  h_chi2_sing_aar         = new TH1F("h_chi2_sing_aar", "Single-Scatter (AAr); Masa's Chi2; Events", 4000, 0, 2000);
  h_sladchi2_sing_aar     = new TH1F("h_sladchi2_sing_aar", "Single-Scatter; Masa's Chi2; Events", 4000, 0, 2000);
  h_xy_sing_aar           = new TH2F("h_xy_sing_aar", "; x [cm]; y [cm]", 80, -20, 20, 80, -20, 20);
  h_xy_sing_aar_chi2cx   = new TH2F("h_xy_sing_aar_chi2cx", "; x [cm]; y [cm]", 80, -20, 20, 80, -20, 20);
  h_sladxy_sing_aar       = new TH2F("h_sladxy_sing_aar", "; x [cm]; y [cm]", 80, -20, 20, 80, -20, 20);
  h_jasonxy_sing_aar           = new TH2F("h_jasonxy_sing_aar", "; x [cm]; y [cm]", 80, -20, 20, 80, -20, 20);
  h_top_var_r_sing_aar    = new TH2F("h_top_var_r_sing_aar", "; R [cm]; Top PMT light Variance", 400, 0, 20, 1000, 0, 1000);
  h_bot_var_r_sing_aar    = new TH2F("h_bot_var_r_sing_aar", "; R [cm]; Bottom PMT light Variance", 400, 0, 20, 1000, 0, 1000);
  h_topbot_var_r_sing_aar = new TH2F("h_topbot_var_r_sing_aar", "; R [cm]; Top Variance / Bottom Variance", 400, 0, 20, 1000, 0, 1000);
  h_maxvar_chi2_sing_aar  = new TH2F("h_maxvar_chi2_sing_aar", "; Min Chi2; Max 1D Variance [cm^{2}]", 6400, 0, 3200, 1000, 0, 1000);
  h_rvar_chi2_sing_aar  = new TH2F("h_rvar_chi2_sing_aar", "; Min Chi2; Max/Min Variance", 6400, 0, 3200, 1000, 0.8, 50.8);
  list->Add(h_chi2_sing_aar); list->Add(h_sladchi2_sing_aar);list->Add(h_xy_sing_aar);list->Add(h_xy_sing_aar_chi2cx);
  list->Add(h_sladxy_sing_aar); list->Add(h_jasonxy_sing_aar); list->Add(h_top_var_r_sing_aar);list->Add(h_bot_var_r_sing_aar);list->Add(h_topbot_var_r_sing_aar);
  list->Add(h_maxvar_chi2_sing_aar); list->Add(h_rvar_chi2_sing_aar);
  
  h_rerr_r_sing_aar   = new TH2F("h_rerr_r_sing_aar", "; r [cm]; Error in r [cm]", 400, 0, 20, 400, 0, 20);
  h_maxerr_r_sing_aar = new TH2F("h_maxerr_r_sing_aar", "; r [cm]; Max(xerr, yerr) [cm]", 400, 0, 20, 400, 0, 20);
  h_TBchi2_sing_aar   = new TH1F("h_TBchi2_sing_aar", "AAr 1-Scatter; Masa's Bot/Top Chi2; Events", 1000, 0, 1000);
  h_TB_sing_aar       = new TH1F("h_TB_sing_aar", "Single-Scatter; Masa's Bot/Top Chi2; Events", 1000, -2, 3);
  list->Add(h_rerr_r_sing_aar); list->Add(h_maxerr_r_sing_aar);list->Add(h_TBchi2_sing_aar);list->Add(h_TB_sing_aar);

  
  //AAr UAr after experimental f90
  int s1_min = 20; int s1_max = 820;
  h_f90_s1_uar = new TH2F("h_f90_s1_uar", "h_f90_s1_uar; S1_corr [pe]; f90", s1_max-s1_min, s1_min, s1_max, 400, 0, 1);
  h_f90_s1_chi2cx_uar[0] = new TH2F("h_f90_s1_chi2cx_uar_0", "h_f90_s1_uar; S1_corr [pe]; f90", s1_max-s1_min, s1_min, s1_max, 400, 0, 1);
  h_f90_s1_chi2cx_uar[1] = new TH2F("h_f90_s1_chi2cx_uar_1", "h_f90_s1_uar; S1_corr [pe]; f90", s1_max-s1_min, s1_min, s1_max, 400, 0, 1);
  list->Add(h_f90_s1_uar); list->Add(h_f90_s1_chi2cx_uar[0]);list->Add(h_f90_s1_chi2cx_uar[1]);

//
//  
//  //slice s1 and make f90 profile
//  foutfile->mkdir("slices")->cd();
//  for (int i=0; i<n_slice; i++) {
//    int s1_l = s1_min + i*s1_width;
//    int s1_h = s1_l+s1_width;
//    int s1_m = s1_l+s1_width/2;
//    h_f90_uar[i] = new TH1F(Form("h_f90_uar_%dPE", s1_m), Form("%dPE-%dPE",s1_l, s1_h), 400, 0, 1);
//    h_f90_uar_chi2cx[0][i] = new TH1F(Form("h_f90_uar_chi2cx_0_%dPE", s1_m), Form("%dPE-%dPE",s1_l, s1_h), 400, 0, 1);
//    h_f90_uar_chi2cx[1][i] = new TH1F(Form("h_f90_uar_chi2cx_1_%dPE", s1_m), Form("%dPE-%dPE",s1_l, s1_h), 400, 0, 1);
//  }
//  foutfile->cd();
//  
//  
//  //save individual event
//  foutfile->mkdir("events")->cd();
//  for(int i=0; i<N_Evt;i++){
//    h_xy_comb_ev[i] = new TH2F(Form("h_xy_comb_ev_%d",i), "; x [cm]; y [cm]", 80, -20, 20, 80, -20, 20);
//    h_xy_comb_ev[i]->SetMarkerColor(kRed); h_xy_comb_ev[i]->SetMarkerStyle(20); h_xy_comb_ev[i]->SetMarkerSize(0.7);
//    h_xy_sing1_ev[i] = new TH2F(Form("h_xy_sing1_%d",i), "; x [cm]; y [cm]", 80, -20, 20, 80, -20, 20);
//    h_xy_sing1_ev[i]->SetMarkerColor(kBlue);h_xy_sing1_ev[i]->SetMarkerStyle(20); h_xy_sing1_ev[i]->SetMarkerSize(0.7);
//    h_xy_sing2_ev[i] = new TH2F(Form("h_xy_sing2_%d",i), "; x [cm]; y [cm]", 80, -20, 20, 80, -20, 20);
//    h_xy_sing2_ev[i]->SetMarkerColor(kGreen); h_xy_sing2_ev[i]->SetMarkerStyle(20); h_xy_sing2_ev[i]->SetMarkerSize(0.7);
//  }
//  foutfile->cd();
  
}



void XY4SLADSelector::FillHistograms(SignalProfile* sigProf){
  
  //  if(Debug) Info("XY4SLADSelector::FillHistograms()","filling histogram...");
  
  
  
  // Generate cuts.
  const bool cx_uar = CX_UAR();
  const bool cx_aar = CX_AAR();
  const bool basic_cuts_aar = CX1() && CX2() && CX3() && CX4();
  const bool basic_cuts_uar = CX1() && CX2() && CX3();
  const bool cx_s1f90 = CX_S1F90();
  const bool cx_2pls = CX_2PLS();
  const bool cx_3pls = CX_3PLS();
  const bool cx_s2true = CX_S2True();
  const bool cx_s2corrtrue = CX_S2CorrTrue();
  const bool cx_tdrift_std = CX_TDrift_STD();
  const bool cx_tdrift = CX_TDrift();
  const bool cx_goodxy = CX_GoodXY(); //radius is within TPC, use both masa and jason
  const bool cx_s2dz = CX_S2dZ();
  const bool cx_sat = CX_Sat();
  const bool cx_trig = CX_TRIG();
  const bool cx_95s1mf = CX_C95S1MF();
  const bool cx_99s1mf = CX_C99S1MF();
  
  
#ifdef use_S1_Range
  const bool cx_s1range = CX_S1Range();
#else
  const bool cx_s1range = 1;
#endif
  


//  //test filling
//  std::cout<<"tdrift: "<<tdrift<<"\n";
//  if (cx_s1range && cx_uar && basic_cuts_uar && cx_3pls && cx_trig && cx_s1f90 && cx_s2true && cx_sat) {
//    std::cout<<"tdrift(pass cut): "<<tdrift<<"\n";
//
//    h_tdrift->Fill(tdrift);
//  }


  
//  -------------------------------------
//   Apply cuts and fill histograms
//  -------------------------------------

  //3 pls events with clear z separation between two s2
  if (cx_s1range && cx_uar && basic_cuts_uar && cx_3pls && cx_trig && cx_s1f90 && cx_s2true && cx_sat && cx_tdrift_std && cx_s2dz) {
    
    //      cout<<"###\n";
    //      cout<<"run_id="<<run_id<<"event_id="<<event_id<<"\n";
    //      cout<<s2_dx<<" "<<s2_dy<<" "<<total_s1_corr<<" "<<total_s2<<" "<<total_s2_top/total_s1_bottom<<"\n";
    //      cout<<s2_dr<<" "<<sigProf[0].chi2XY <<" "<<GetTopPMTVariance(sigProf)<<" "<<sigProf.PosiX<<" "<<sigProf.PosiY<<"\n";
    
    h_s1_comb->Fill(total_s1);
    if (!CX_GoodXY_MASA(sigProf[0])) {
      h_s1_comb_badxy->Fill(total_s1);
      return;
    }
    
    h_s1_comb_goodxy->Fill(total_s1);
    h_chi2_dr->Fill( s2_dr, sigProf[0].chi2XY );
    h_logchi2_dr->Fill( s2_dr, log10(sigProf[0].chi2XY) );
    h_TBchi2_dr->Fill(s2_dr, TBChi2);
    
    Float_t top_var = GetTopPMTVariance(sigProf[0]);
    Float_t bot_var = GetBotPMTVariance(sigProf[0]);
    h_top_var_dr->Fill( s2_dr, top_var);
    h_bot_var_dr->Fill( s2_dr, bot_var);
    h_topbot_var_dr->Fill( s2_dr, top_var/bot_var);
    
    
    float PosiR = sqrt(sigProf[0].PosiX*sigProf[0].PosiX+sigProf[0].PosiY*sigProf[0].PosiY);
    h_top_var_r_comb->Fill(PosiR, top_var);
    h_bot_var_r_comb->Fill(PosiR, bot_var);
    h_topbot_var_r_comb->Fill(PosiR, top_var/bot_var);

    vector<float> var = GetPMTMaxMinVariance(sigProf[0]);
    Float_t max_var = var.at(0); Float_t min_var = var.at(1);
    //Float_t max_var = GetPMT1DMaxVariance(sigProf[0]);    Float_t min_var = GetPMT1DMinVariance(sigProf[0]);
    h_maxvar_chi2_comb->Fill(sigProf[0].chi2XY, max_var);
    h_rvar_chi2_comb->Fill(sigProf[0].chi2XY, max_var/min_var);

    float first_s2_r = sqrt(allpulses_x[1]*allpulses_x[1]+allpulses_y[1]*allpulses_y[1]);
    h_top_var_r[0]->Fill(first_s2_r, GetTopPMTVariance(1));
    h_bot_var_r[0]->Fill(first_s2_r, GetBotPMTVariance(1));
    h_topbot_var_r[0]->Fill(first_s2_r, GetTopPMTVariance(1)/GetBotPMTVariance(1));
    
    float sec_s2_r = sqrt(allpulses_x[2]*allpulses_x[2]+allpulses_y[2]*allpulses_y[2]);
    h_top_var_r[1]->Fill(sec_s2_r, GetTopPMTVariance(2));
    h_bot_var_r[1]->Fill(sec_s2_r, GetBotPMTVariance(2));
    h_topbot_var_r[1]->Fill(first_s2_r, GetTopPMTVariance(2)/GetBotPMTVariance(2));
    
    h_s2s1_vs_s1[0]->Fill(total_s1_corr, log10(total_s2_corr/total_s1_corr));
    h_s2s1_vs_s1[1]->Fill(total_s1_corr, log10(pulse_total_npe[2]*allpulses_xycorr_factor[2]/total_s1_corr));
    h_xy_comb->Fill(sigProf[0].PosiX, sigProf[0].PosiY);
    h_xy[0]->Fill(allpulses_x[1], allpulses_y[1]);
    h_xy[1]->Fill(allpulses_x[2], allpulses_y[2]);
    h_chi2_comb->Fill(sigProf[0].chi2XY);
    h_TBchi2_comb->Fill(TBChi2);
    h_chi2[0]->Fill(allpulses_chi2[1]);
    h_chi2[1]->Fill(allpulses_chi2[2]);
    h_TB_comb->Fill(GetS2LogTopBot(sigProf[0]));
    
    
    //      cout<<run_id<<" "<< event_id<<" "<<sigProf[0].maxchi2<<" "<<sigProf[0].chi2XY<<"\n";
    h_rchi2_comb->Fill(sigProf[0].maxchi2);
    
    float rchi2 = sigProf[0].maxchi2/sigProf[0].chi2XY;
    h_rchi2_chi2_comb->Fill(sigProf[0].chi2XY, rchi2);
    h_fracchi2_chi2_comb->Fill(sigProf[0].chi2XY, log10(rchi2)/log10(sigProf[0].chi2XY));
    h_fracchi2_logchi2_comb->Fill(log10(sigProf[0].chi2XY), log10(rchi2)/log10(sigProf[0].chi2XY));
    h_maxchi2_chi2_comb->Fill(log10(sigProf[0].chi2XY), log10(sigProf[0].maxchi2));
    
    float dr_err = sqrt(sigProf[0].PosiYErr*sigProf[0].PosiYErr+sigProf[0].PosiXErr*sigProf[0].PosiXErr);
    h_rerr_dr->Fill(s2_dr, dr_err);
    h_rerr_r_comb->Fill(PosiR, dr_err);
    float max_err = TMath::Max(sigProf[0].PosiYErr, sigProf[0].PosiXErr);
    h_maxerr_dr->Fill(s2_dr, max_err);
    h_maxerr_r_comb->Fill(PosiR, max_err);

    
    float comb_s2 = pulse_total_npe[1]+pulse_total_npe[2];
    float exp_x = 1./comb_s2 * (pulse_total_npe[1]*allpulses_x[1]+pulse_total_npe[2]*allpulses_x[2]);
    float exp_y = 1./comb_s2 * (pulse_total_npe[1]*allpulses_y[1]+pulse_total_npe[2]*allpulses_y[2]);
    float dd = sqrt(pow((sigProf[0].PosiX-exp_x),2)+pow((sigProf[0].PosiY-exp_y),2));
    h_dexpxy_dr->Fill(s2_dr, dd);
    
    //experimenting cutting on chi2
    if (sigProf[0].chi2XY<minchi2_cx) {
      h_xy_comb_chi2cx->Fill(sigProf[0].PosiX, sigProf[0].PosiY);
    }
    
//    //#define random_save_events
//#ifdef random_save_events
//    //randomly choose N_Evt
//    if (n%20==0 && counter<N_Evt) {
//      cout<<run_id<<" "<<event_id<<" "<<s2_dr<<" "<<sigProf[0].chi2XY<<"\n";
//    }
//#else
//    if (s2_dr>14.9 && s2_dr<15 && counter<N_Evt) {
//      cout<<run_id<<" "<<event_id<<" "<<s2_dr<<" "<<sigProf[0].chi2XY<<"\n";
//#endif
//      h_xy_comb_ev[counter]->Fill(sigProf[0].PosiX, sigProf[0].PosiY);
//      h_xy_sing1_ev[counter]->Fill(allpulses_x[1], allpulses_y[1]);
//      h_xy_sing2_ev[counter]->Fill(allpulses_x[2], allpulses_y[2]);
//      h_xy_comb_ev[counter]->SetTitle(Form("run_id=%d event_id=%d", run_id, event_id));
//      counter++;
//    }
  }
  

    //  3 pls events, UAr, no s2_dz cut
    //  in 50d paper, we use 0.93 mm/us drift velocity
    if (cx_s1range && cx_uar && basic_cuts_uar && cx_3pls && cx_trig && cx_s1f90 && cx_s2true && cx_sat && cx_tdrift) {
  
      //masa+jason xy
      float dR_sq = pow((x[2]-x[1]),2) + pow((y[2]-y[1]),2) + pow((pulse_start_time[2]*0.093-pulse_start_time[1]*0.093),2);
      h_R->Fill(sqrt(dR_sq));
      if (total_s1>20 && total_s1<450) {
        h_R_s1cx->Fill(sqrt(dR_sq));
      }
  
      //masa's xy
      dR_sq = pow((allpulses_x[2]-allpulses_x[1]),2) + pow((allpulses_y[2]-allpulses_y[1]),2) + pow((pulse_start_time[2]*0.093-pulse_start_time[1]*0.093),2);
      h_R_masas->Fill(sqrt(dR_sq));
      if (total_s1>20 && total_s1<450) {
        h_R_masas_s1cx->Fill(sqrt(dR_sq));
      }
  
      //jason's xy
      dR_sq = pow((allpulses_xyl_x[2]-allpulses_xyl_x[1]),2) + pow((allpulses_xyl_y[2]-allpulses_xyl_y[1]),2) + pow((pulse_start_time[2]*0.093-pulse_start_time[1]*0.093),2);
      h_R_jasons->Fill(sqrt(dR_sq));
      if (total_s1>20 && total_s1<450) {
        h_R_jasons_s1cx->Fill(sqrt(dR_sq));
      }
  
    }

  
    //  2pls events, UAr
    if (cx_s1range && cx_uar && basic_cuts_uar && cx_2pls && cx_trig && cx_s1f90 && cx_s2true && cx_sat && cx_tdrift_std) {
  
      h_s1_sing->Fill(total_s1);
      if (!CX_GoodXY_MASA(sigProf[1])) {
        h_s1_sing_badxy->Fill(total_s1);
        return;
      }
  
      h_s1_sing_goodxy->Fill(total_s1);
      h_chi2_sing->Fill(sigProf[1].chi2XY);
      h_sladchi2_sing_uar->Fill(allpulses_chi2[1]);
      h_TBchi2_sing->Fill(TBChi2);
      h_rchi2_sing->Fill(sigProf[1].maxchi2);
      h_xy_sing->Fill(sigProf[1].PosiX, sigProf[1].PosiY);
      h_sladxy_sing_uar->Fill(allpulses_x[1], allpulses_y[1]);
      h_jasonxy_sing->Fill(allpulses_xyl_x[1], allpulses_xyl_y[1]);
      h_TB_sing->Fill(total_s2_top/total_s2_bottom);
  
      float rchi2 = sigProf[1].maxchi2/sigProf[1].chi2XY;
      h_rchi2_chi2_sing->Fill(sigProf[1].chi2XY, rchi2);
      h_fracchi2_chi2_sing->Fill(sigProf[1].chi2XY, log10(rchi2)/log10(sigProf[1].chi2XY));
      h_fracchi2_logchi2_sing->Fill(log10(sigProf[1].chi2XY), log10(rchi2)/log10(sigProf[1].chi2XY));
      h_maxchi2_chi2_sing->Fill(log10(sigProf[1].chi2XY), log10(sigProf[1].maxchi2));
  
      float PosiR = sqrt(sigProf[1].PosiX*sigProf[1].PosiX+sigProf[1].PosiY*sigProf[1].PosiY);
      h_top_var_r_sing->Fill(PosiR, GetTopPMTVariance(1));
      h_bot_var_r_sing->Fill(PosiR, GetBotPMTVariance(1));
      h_topbot_var_r_sing->Fill(PosiR, GetTopPMTVariance(1)/GetBotPMTVariance(1));
  
      vector<float> var = GetPMTMaxMinVariance(sigProf[1]);
      Float_t max_var = var.at(0);  Float_t min_var =var.at(1);
      //      Float_t max_var =  GetPMT1DMaxVariance(sigProf[1]);      Float_t min_var =  GetPMT1DMinVariance(sigProf[1]);
      h_maxvar_chi2_sing->Fill(sigProf[1].chi2XY, GetPMT1DMaxVariance(sigProf[1]));
      h_rvar_chi2_sing->Fill(sigProf[1].chi2XY, max_var/min_var);

      float dr_err = sqrt(sigProf[1].PosiYErr*sigProf[1].PosiYErr+sigProf[1].PosiXErr*sigProf[1].PosiXErr);
      h_rerr_r_sing->Fill( PosiR, dr_err);
  
      float max_err = TMath::Max(sigProf[1].PosiYErr, sigProf[1].PosiXErr);
      h_maxerr_r_sing->Fill( PosiR, max_err);
      
      //experimenting cuts
      if (sigProf[2].chi2XY<minchi2_cx) {
        h_xy_sing_aar_chi2cx->Fill(sigProf[2].PosiX, sigProf[2].PosiY);
      }
  
    }

//  2pls events, AAr
    if (cx_s1range && cx_aar && basic_cuts_aar && cx_2pls && cx_trig && cx_s1f90 && cx_s2true && cx_sat && cx_tdrift_std) {
  
      h_chi2_sing_aar->Fill(sigProf[2].chi2XY);
      h_sladchi2_sing_aar->Fill(allpulses_chi2[1]);
      h_xy_sing_aar->Fill(sigProf[2].PosiX, sigProf[2].PosiY);
      h_sladxy_sing_aar->Fill(allpulses_x[1], allpulses_y[1]);
      h_jasonxy_sing_aar->Fill(allpulses_xyl_x[1],allpulses_xyl_y[1]);
      h_TBchi2_sing_aar->Fill(TBChi2);
      h_TB_sing_aar->Fill(total_s2_top/total_s2_bottom);
      h_rchi2_sing_aar->Fill(sigProf[2].maxchi2);
  
      //      cout<<run_id<<" "<< event_id<<" "<<sigProf[2].maxchi2<<" "<<sigProf[2].chi2XY<<"\n";
  
  
      float rchi2 = sigProf[2].maxchi2/sigProf[2].chi2XY;
      h_rchi2_chi2_sing_aar->Fill(sigProf[2].chi2XY, rchi2);
      h_fracchi2_chi2_sing_aar->Fill(sigProf[2].chi2XY, log10(rchi2)/log10(sigProf[2].chi2XY));
      h_fracchi2_logchi2_sing_aar->Fill(log10(sigProf[2].chi2XY), log10(rchi2)/log10(sigProf[2].chi2XY));
      h_maxchi2_chi2_sing_aar->Fill(log10(sigProf[2].chi2XY),log10(sigProf[2].maxchi2));
  
      float PosiR = sqrt(sigProf[2].PosiX*sigProf[2].PosiX+sigProf[2].PosiY*sigProf[2].PosiY);
      h_top_var_r_sing_aar->Fill(PosiR, GetTopPMTVariance(1));
      h_bot_var_r_sing_aar->Fill(PosiR, GetBotPMTVariance(1));
      h_topbot_var_r_sing_aar->Fill(PosiR, GetTopPMTVariance(1)/GetBotPMTVariance(1));
      
      vector<float> var = GetPMTMaxMinVariance(sigProf[2]);
      // Float_t max_var = GetPMT1DMaxVariance(sigProf[2]);  Float_t min_var = GetPMT1DMinVariance(sigProf[2]);
      Float_t max_var = var.at(0); Float_t min_var = var.at(1);
      h_maxvar_chi2_sing_aar->Fill(sigProf[2].chi2XY, max_var);
      h_rvar_chi2_sing_aar->Fill(sigProf[2].chi2XY, max_var/min_var);

      float dr_err = sqrt(sigProf[2].PosiYErr*sigProf[2].PosiYErr+sigProf[2].PosiXErr*sigProf[2].PosiXErr);
      h_rerr_r_sing_aar->Fill( PosiR, dr_err);
  
      float max_err = TMath::Max(sigProf[2].PosiYErr, sigProf[2].PosiXErr);
      h_maxerr_r_sing_aar->Fill( PosiR, max_err);
  
    }
  
  
  //  standard wimp search cut, UAr
  if (cx_uar && basic_cuts_uar && cx_2pls && cx_trig && cx_s2corrtrue && cx_sat && cx_tdrift_std && cx_95s1mf) {
    h_f90_s1_uar->Fill(total_s1_corr, total_f90);
    if (sigProf[1].chi2XY<6) h_f90_s1_chi2cx_uar[0]->Fill(total_s1_corr, total_f90);
    if (sigProf[1].chi2XY<2) h_f90_s1_chi2cx_uar[1]->Fill(total_s1_corr, total_f90);
    
//    for (int i=0; i<n_slice; i++) {
//      float s1_l = s1_min + i*s1_width;
//      float s1_h = s1_l + s1_width;
//      if(s1_l<total_s1_corr && total_s1_corr<s1_h){
//        h_f90_uar[i]->Fill(total_f90);
//        if (sigProf[1].chi2XY<6) h_f90_uar_chi2cx[0][i]->Fill(total_f90);
//        if (sigProf[1].chi2XY<2) h_f90_uar_chi2cx[1][i]->Fill(total_f90);
//        break;
//      }
//    }
  }
  
} //end FillHistograms






