

//  Usage:
//     root -b -q plotter.C+




#include "../uarlimit2015//Macros/rootstart.h"

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


void plotter(){

  //Glabal style
  SetMyStyle();

  TFile* f = new TFile("XY_MultiScatter_test_3.root");
  TCanvas* c[20];
  
  
  //--------------------------------
  //compare xy in slad to newly tuned xy
  //UAr
  TString str = Form("h_sladxy_xy_uar");
  c[0] = new TCanvas(str,str, 1200, 500);
  c[0]->Divide(2,1); c[0]->cd(1); gPad->SetLogz();
  
  TH1F* h_sladxy_sing_uar = (TH1F*) f->Get("h_sladxy_sing_uar");
  h_sladxy_sing_uar->SetTitle("UAr - xy in SLAD");
  h_sladxy_sing_uar->Draw("colz");

  c[0]->cd(2); gPad->SetLogz();
  TH1F* h_xy_sing = (TH1F*) f->Get("h_xy_sing");
  h_xy_sing->SetTitle("UAr - xy after tuned");
  h_xy_sing->Draw("colz");
  
  c[0]->SaveAs(TString("plots/"+str+".pdf"));
  
  //AAr
  str = Form("h_sladxy_xy_aar");
  c[1] = new TCanvas(str,str, 1200, 500);
  c[1]->Divide(2,1); c[1]->cd(1); gPad->SetLogz();
  
  TH1F* h_sladxy_sing_aar = (TH1F*) f->Get("h_sladxy_sing_aar");
  h_sladxy_sing_aar->SetTitle("AAr - xy in SLAD");
  h_sladxy_sing_aar->Draw("colz");
  
  c[1]->cd(2); gPad->SetLogz();
  TH1F* h_xy_sing_aar = (TH1F*) f->Get("h_xy_sing_aar");
  h_xy_sing_aar->SetTitle("AAr - xy after tuned");
  h_xy_sing_aar->Draw("colz");
  
  c[1]->SaveAs(TString("plots/"+str+".pdf"));
  
  
  
  //----------------------------------
  //Chi2 radial Separation
  str = Form("h_logchi2_dr");
  c[2] = new TCanvas(str,str, 700, 500);
  c[2]->cd(); c[2]->SetLogz(); c[2]->SetGrid();
  
  TH2F* h_logchi2_dr = (TH2F*) f->Get("h_logchi2_dr");
  h_logchi2_dr->GetYaxis()->SetRangeUser(-0.5, 3.8);
  h_logchi2_dr->GetXaxis()->SetTitleOffset(1.1);
  h_logchi2_dr->GetXaxis()->SetTitle("xy separation of two S2s [cm]");
  h_logchi2_dr->Draw("colz");
  
  c[2]->SaveAs(TString("plots/"+str+".pdf"));
  
  
  
  //-------------------------------------
  //Compare SLAD Chi2 and after re-definition
  str = Form("h_compare_sladchi2_chi2");
  c[3] = new TCanvas(str,str, 700, 500);
  c[3]->cd(); c[3]->SetLogz(); c[3]->SetGrid();
  
  TH1F* h_sladchi2_sing_aar = (TH1F*) f->Get("h_sladchi2_sing_aar");
  h_sladchi2_sing_aar->Scale(1./h_sladchi2_sing_aar->Integral());
  h_sladchi2_sing_aar->SetLineColor(kBlack);
  h_sladchi2_sing_aar->SetLineStyle(2);
  h_sladchi2_sing_aar->GetYaxis()->SetTitle("Normalized by total events");
  h_sladchi2_sing_aar->Draw();
  
  TH1F* h_sladchi2_sing_uar = (TH1F*) f->Get("h_sladchi2_sing_uar");
  h_sladchi2_sing_uar->Scale(1./h_sladchi2_sing_uar->Integral());
  
  TH1F* h_chi2_sing_aar = (TH1F*) f->Get("h_chi2_sing_aar");
  h_chi2_sing_aar->Scale(1./h_chi2_sing_aar->Integral());
  
  TH1F* h_chi2_sing = (TH1F*) f->Get("h_chi2_sing");
  h_chi2_sing->Scale(1./h_chi2_sing->Integral());

  

  
  
  


  
  
  
  
  
  
  

}
