#include <numeric>      // std::accumulate

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
//#include "/ds50/app/user/xxiang89/work/xy4slad_ms/XYReconstructor/XYRecHelper.h"

#include "PMTGeom.h"
#include "XYRecHelper.h"

using namespace std;

ClassImp(XYRecHelper)

//#define NOSTACK
TString xtitle("X [cm]"), ytitle("Y [cm]"), ztitle(/*"#chi^{2}"*/"Chi2");

void XYRecHelper::Draw2DFunc(Double_t x1, Double_t x2, SignalProfile sigProf){
  TString title(Form("Run: %d, Event Id: %d", sigProf.run_id, sigProf.event_id));
  TCanvas *c1 = new TCanvas("c1",title.Data(),1000, 1000);//750);
  c1->cd();
  c1->Draw();
  c1->Divide(2,2);
  c1->SetLogz();

  Double_t r(TMath::Sqrt(x1*x1+x2*x2)), theta(TMath::ATan2(x2, x1));
  cout<<"X: "<<x1<<" Y: "<<x2<<endl; // convert xy
  cout<<"R: "<<r<<" Theta: "<<theta<<endl;
//  Double_t min_x1 = -PMTGeom::tpcRadi;
//  Double_t max_x1 = PMTGeom::tpcRadi;
//  Double_t min_x2 = -PMTGeom::tpcRadi;
//  Double_t max_x2 = PMTGeom::tpcRadi;
//  TF2 f2D("f2D", fMyFcn, min_x1, max_x1, min_x2, max_x2, 0, "LikelihoodFunction");//TF2(const char* name, void* ptr, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Int_t npar, const char* className)
//  f2D.SetNpx(180);
//  f2D.SetNpy(180);

  TVirtualPad* pad = c1->cd(1);
  pad->SetRightMargin(0.12);
  pad->SetLogz();
  gf2D->Draw("COLZ");
  TH2* htmp = (TH2*) gf2D->GetHistogram();
  htmp->SetXTitle(xtitle.Data());
  htmp->SetYTitle(ytitle.Data());
  htmp->SetZTitle(ztitle.Data());
  htmp->SetTitleOffset(0.5, "Z");
  htmp->SetTitle(title.Data());
  TEllipse el(0., 0., PMTGeom::tpcRadi);
  el.SetFillStyle(0);
  el.Draw("same");

  TEllipse posi(x1, x2, 0.5);//event position
  posi.SetFillColor(kWhite);//kBlack);
  posi.Draw("same");

  TLatex L;
//  L.SetNDC();
  L.SetTextSize(0.03);

  TGraph2D g2D(PMTGeom::N_CHANNELS/2);
  g2D.SetNpy(100);
  g2D.SetNpy(100);
  g2D.SetMarkerStyle(20);
//  g2D.GetXaxis()->SetTitle(xtitle.Data());
//  g2D.GetYaxis()->SetTitle(ytitle.Data());
  Double_t x_1[PMTGeom::N_CHANNELS/2], x_2[PMTGeom::N_CHANNELS/2], z[PMTGeom::N_CHANNELS/2];
//  TEllipse** pmt = new TEllipse*[PMTGeom::N_CHANNELS/2];
  TEllipse* pmt;
  Int_t count(0);
  for (Int_t ch = PMTGeom::N_CHANNELS/2; ch < PMTGeom::N_CHANNELS; ch++)
  {
      x_1[count] = PMTGeom::pmtUnit*PMTGeom::pmt_x[ch];
      x_2[count] = PMTGeom::pmtUnit*PMTGeom::pmt_y[ch];
      z[count]   = sigProf.s2[ch]/sigProf.total_s2;

      g2D.SetPoint(count, x_1[count], x_2[count], z[count]);
      L.DrawLatex(x_1[count], x_2[count], Form("%d",ch));

      pmt = new TEllipse(PMTGeom::pmtUnit*PMTGeom::pmt_x[ch], PMTGeom::pmtUnit*PMTGeom::pmt_y[ch], 3.81);//position of pmts // 3.81 cm is pmt radius
//      TEllipse pmt(PMTGeom::pmtUnit*PMTGeom::pmt_x[ch], PMTGeom::pmtUnit*PMTGeom::pmt_y[ch], 3.81);//position of pmts // 3.81 cm is pmt radius
      pmt->SetFillStyle(0);
      pmt->SetLineStyle(2);
      pmt->Draw("same");
//      pmt.Draw("same");
//      pmt[ch]->Draw("same");

      count++;
  }
//  TGraph2D g2D(PMTGeom::N_CHANNELS/2, x_1, x_2, z);


  pad = c1->cd(2);
  pad->SetRightMargin(0.15);
  pad->SetLogz();
  g2D.Draw("pcolz");//("surf1");//("colz");//("surf1");
  g2D.GetXaxis()->SetTitle(xtitle.Data());
  g2D.GetYaxis()->SetTitle(ytitle.Data());
  g2D.GetYaxis()->SetTitleOffset(1.3);
  g2D.GetXaxis()->SetTitleOffset(1.3);
  g2D.SetTitle("S2_{ch}/S2_{total}");
//  f2D.Draw("pol");

  pad = c1->cd(3);
  pad->SetLogy();
  THStack* hs_Chi2 = DrawChi2alongR(theta, sigProf);
#ifdef NOSTACK
  hs_Chi2->Draw("nostack");
#else
  hs_Chi2->Draw();
#endif
  hs_Chi2->GetYaxis()->SetTitleOffset(1.3);
  TLegend* leg = gPad->BuildLegend(0.9,0.1,1,0.9);
//  leg->SetNColumns(2);

  pad = c1->cd(4);
  pad->SetLogy();

  THStack* hs_Chi2_theta = DrawChi2alongTheta(r, sigProf);
#ifdef NOSTACK
  hs_Chi2_theta->Draw("nostack");
#else
  hs_Chi2_theta->Draw();
#endif
  hs_Chi2_theta->GetYaxis()->SetTitleOffset(1.3);
  leg = gPad->BuildLegend(0.9,0.1,1,0.9);
//  leg->SetNColumns(2);


  c1->Update();
  TString fullname = "/Users/xinxiang/ds50analysis/xy4slad/single_s2/" + TString(Form("r%de%d.pdf",sigProf.run_id, sigProf.event_id));
  c1->SaveAs(fullname);
  fRint->Run();

//  f2D.Draw("surf1 pol");
//  f2D.SetMaximum(10);
  TList *tmplist = hs_Chi2->GetHists();
  tmplist->Delete();
  tmplist = hs_Chi2_theta->GetHists();
  tmplist->Delete();
//  for (Int_t ch = PMTGeom::N_CHANNELS/2; ch < PMTGeom::N_CHANNELS; ch++) delete pmt[ch];
  delete hs_Chi2;
  delete hs_Chi2_theta;
 delete c1;
}

THStack* XYRecHelper::DrawChi2alongR(Double_t theta, SignalProfile sigProf){
  THStack* hs_Chi2 = new THStack("hs_Chi2", Form("Chi2 along R (#theta = %3.2f);R [cm];%s",theta, ztitle.Data()));

  Int_t nRbins(200);
  Double_t rmin(0.), rmax(PMTGeom::tpcRadi);
  Double_t sumweight(0.);
  for(Int_t ich=0; ich<PMTGeom::N_CHANNELS; ich++){
      if (!fMyFcn->IsChanIncludedInChi2 (ich)) continue;
//#ifdef NOSTACK
//      if(sigProf.s2[ich]/sigProf.total_s2< 1./(1.5*PMTGeom::N_CHANNELS)) continue;
//#endif

      TH1F* htmp = new TH1F(Form("htmp_%d",ich), Form("chan %d",ich), nRbins, rmin, rmax);
      TAxis *xaxis = htmp->GetXaxis();
      for(Int_t rbin=1; rbin <nRbins+1; rbin++){
          Double_t r = xaxis->GetBinCenter(rbin);
          Double_t par[] = {r*TMath::Cos(theta), r*TMath::Sin(theta) }; // convert xy
          Double_t chi2(0.), weight(0.);
          fMyFcn->GetChi2Channel(par, ich, chi2, weight);
          htmp->SetBinContent(rbin, chi2);
          sumweight += weight;
      }
      htmp->SetLineColor(PMTGeom::colors[ich%19]);
#ifndef NOSTACK
      htmp->SetFillColor(PMTGeom::colors[ich%19]);
#endif
      hs_Chi2->Add(htmp);
//      delete htmp;
  }

  TH1F* htmp(NULL);
  TAxis *xaxis(NULL);

  if(fMyFcn->Is_S2OVERS1PROFILE_Used()){
      htmp = new TH1F("htmp_s2overs1", "S2/S1", nRbins, rmin, rmax);
      xaxis = htmp->GetXaxis();
      for (Int_t rbin = 1; rbin < nRbins + 1; rbin++)
        {
          Double_t r = xaxis->GetBinCenter(rbin);
          Double_t par[] =
            { r * TMath::Cos(theta), r * TMath::Sin(theta) }; // convert xy
          Double_t chi2(0.), weight(0.);
          fMyFcn->GetChi2S2overS1(par, chi2, weight);
          htmp->SetBinContent(rbin, chi2);
        }
      htmp->SetLineColor(kGray);
#ifndef NOSTACK
      htmp->SetFillColor(kGray);
#endif
      hs_Chi2->Add(htmp);
  }

  if(fMyFcn->Is_S2BOTOVERS2PROFILE_Used()){
      htmp = new TH1F("htmp_s2Botovers2", "S2_{bottom}/S2", nRbins, rmin, rmax);
      xaxis = htmp->GetXaxis();
      for (Int_t rbin = 1; rbin < nRbins + 1; rbin++)
        {
          Double_t r = xaxis->GetBinCenter(rbin);
          Double_t par[] =
            { r * TMath::Cos(theta), r * TMath::Sin(theta) }; // convert xy
          Double_t chi2(0.), weight(0.);
          fMyFcn->GetChi2S2BotoverS2(par, chi2, weight);
          htmp->SetBinContent(rbin, chi2);
        }
      htmp->SetLineColor(kGray + 3);
#ifndef NOSTACK
      htmp->SetFillColor(kGray + 3);
#endif
      hs_Chi2->Add(htmp);
    }

  if(fMyFcn->Is_S2BOTOVERS2TOPPROFILE_Used()){
      htmp = new TH1F("htmp_s2Botovers2Top", "S2_{bottom}/S2_{top}", nRbins,
          rmin, rmax);
      xaxis = htmp->GetXaxis();
      for (Int_t rbin = 1; rbin < nRbins + 1; rbin++)
        {
          Double_t r = xaxis->GetBinCenter(rbin);
          Double_t par[] =
            { r * TMath::Cos(theta), r * TMath::Sin(theta) }; // convert xy
          Double_t chi2(0.), weight(0.);
          fMyFcn->GetChi2S2BotoverS2Top(par, chi2, weight);
          htmp->SetBinContent(rbin, chi2);
        }
      htmp->SetLineColor(kGray + 3);
#ifndef NOSTACK
      htmp->SetFillColor(kGray + 3);
#endif
      hs_Chi2->Add(htmp);
  }
  if(fMyFcn->Is_DENSITYPROFILE_Used()){
      htmp = new TH1F("htmp_density", "Density", nRbins, rmin, rmax);
      xaxis = htmp->GetXaxis();
      for (Int_t rbin = 1; rbin < nRbins + 1; rbin++)
        {
          Double_t r = xaxis->GetBinCenter(rbin);
          Double_t par[] =
            { r * TMath::Cos(theta), r * TMath::Sin(theta) }; // convert xy
          Double_t chi2(0.), weight(0.);
          fMyFcn->GetDensityContribution(par, chi2, weight);
          htmp->SetBinContent(rbin, chi2);
        }
      htmp->SetLineColor(kPink - 9);
#ifndef NOSTACK
      htmp->SetFillColor(kPink - 9);
#endif
      hs_Chi2->Add(htmp);
  }

  hs_Chi2->GetStack()->SetOwner();

  return hs_Chi2;
}

THStack* XYRecHelper::DrawChi2alongTheta(Double_t r, SignalProfile sigProf){
  THStack* hs_Chi2_theta = new THStack("hs_Chi2_theta", Form("Chi2 along #theta (r = %3.2f);#theta [rad];%s",r, ztitle.Data()));

  Int_t nbins(600);
  Double_t min(-TMath::Pi()), max(TMath::Pi());
  for(Int_t ich=0; ich<PMTGeom::N_CHANNELS; ich++){
      if (!fMyFcn->IsChanIncludedInChi2 (ich)) continue;
//#ifdef NOSTACK
//      if(sigProf.s2[ich]/sigProf.total_s2< 1./(1.5*PMTGeom::N_CHANNELS)) continue;
//#endif

      TH1F* htmp = new TH1F(Form("htmp_theta_%d",ich), Form("chan %d",ich), nbins, min, max);
      TAxis *xaxis = htmp->GetXaxis();
      for(Int_t ibin=1; ibin <nbins+1; ibin++){
          Double_t theta = xaxis->GetBinCenter(ibin);
          Double_t par[] = {r*TMath::Cos(theta), r*TMath::Sin(theta) }; // convert xy
          Double_t chi2(0.), weight(0.);
          fMyFcn->GetChi2Channel(par, ich, chi2, weight);
//          if(ich==26) cout<<"theta: "<<theta<<" chi2: "<<chi2<<endl;
          htmp->SetBinContent(ibin, chi2);
      }
      htmp->SetLineColor(PMTGeom::colors[ich%19]);
#ifndef NOSTACK
      htmp->SetFillColor(PMTGeom::colors[ich%19]);
#endif
      hs_Chi2_theta->Add(htmp);
//      delete htmp;
  }

  TH1F* htmp(NULL);
  TAxis *xaxis(NULL);
  if(fMyFcn->Is_S2OVERS1PROFILE_Used()){
      htmp = new TH1F("htmp_theta_s2overs1", "S2/S1", nbins, min, max);
      xaxis = htmp->GetXaxis();
      for (Int_t ibin = 1; ibin < nbins + 1; ibin++)
        {
          Double_t theta = xaxis->GetBinCenter(ibin);
          Double_t par[] =
            { r * TMath::Cos(theta), r * TMath::Sin(theta) }; // convert xy
          Double_t chi2(0.), weight(0.);
          fMyFcn->GetChi2S2overS1(par, chi2, weight);
          htmp->SetBinContent(ibin, chi2);
        }
      htmp->SetLineColor(kGray);
#ifndef NOSTACK
      htmp->SetFillColor(kGray);
#endif
      hs_Chi2_theta->Add(htmp);
  }

  if(fMyFcn->Is_S2BOTOVERS2PROFILE_Used()){
      htmp = new TH1F("htmp_theta_s2Botovers2", "S2_{bottom}/S2", nbins, min,
          max);
      xaxis = htmp->GetXaxis();
      for (Int_t ibin = 1; ibin < nbins + 1; ibin++)
        {
          Double_t theta = xaxis->GetBinCenter(ibin);
          Double_t par[] =
            { r * TMath::Cos(theta), r * TMath::Sin(theta) }; // convert xy
          Double_t chi2(0.), weight(0.);
          fMyFcn->GetChi2S2BotoverS2(par, chi2, weight);
          htmp->SetBinContent(ibin, chi2);
        }
      htmp->SetLineColor(kGray + 3);
#ifndef NOSTACK
      htmp->SetFillColor(kGray + 3);
#endif
      hs_Chi2_theta->Add(htmp);
  }

  if(fMyFcn->Is_S2BOTOVERS2TOPPROFILE_Used()){
      htmp = new TH1F("htmp_theta_s2Botovers2Top", "S2_{bottom}/S2_{top}",
          nbins, min, max);
      xaxis = htmp->GetXaxis();
      for (Int_t ibin = 1; ibin < nbins + 1; ibin++)
        {
          Double_t theta = xaxis->GetBinCenter(ibin);
          Double_t par[] =
            { r * TMath::Cos(theta), r * TMath::Sin(theta) }; // convert xy
          Double_t chi2(0.), weight(0.);
          fMyFcn->GetChi2S2BotoverS2Top(par, chi2, weight);
          htmp->SetBinContent(ibin, chi2);
        }
      htmp->SetLineColor(kGray + 3);
#ifndef NOSTACK
      htmp->SetFillColor(kGray + 3);
#endif
      hs_Chi2_theta->Add(htmp);
  }

  if(fMyFcn->Is_DENSITYPROFILE_Used()){
      htmp = new TH1F("htmp_theta_density", "Density", nbins, min, max);
      xaxis = htmp->GetXaxis();
      for (Int_t ibin = 1; ibin < nbins + 1; ibin++)
        {
          Double_t theta = xaxis->GetBinCenter(ibin);
          Double_t par[] =
            { r * TMath::Cos(theta), r * TMath::Sin(theta) }; // convert xy
          Double_t chi2(0.), weight(0.);
          fMyFcn->GetDensityContribution(par, chi2, weight);
          htmp->SetBinContent(ibin, chi2);
//      cout<<"density: "<<chi2<<endl;
        }
      htmp->SetLineColor(kPink - 9);
#ifndef NOSTACK
      htmp->SetFillColor(kPink - 9);
#endif
      hs_Chi2_theta->Add(htmp);
  }

  return hs_Chi2_theta;
}

void XYRecHelper::PrintS2(SignalProfile sigProf){
  cout<<"Event ID: "<<sigProf.event_id<<endl;
  for (Int_t ch = 0; ch < PMTGeom::N_CHANNELS; ch++)
  {
      if(sigProf.s2[ch]/sigProf.total_s2< 1./PMTGeom::N_CHANNELS) continue;
     cout<<"ch: "<<ch<<" s2[ch]: "<<sigProf.s2[ch];
     cout<<" s2[ch]/total_s2: "<<sigProf.s2[ch]/sigProf.total_s2<<endl;
//
//     cout<<Form("ch: %2d  s2[ch]: %4.1f   s2[ch]/total_s2: %3.2f", ch, sigProf.s2[ch], sigProf.s2[ch]/sigProf.total_s2)<<endl;
  }
  cout<<Form("s2[bottom]/total_s2: %3.2f", sigProf.s2_bottom/sigProf.total_s2)<<endl;
  cout<<Form("Log10(s2[bottom]/s2[top]): %3.2f", TMath::Log10(sigProf.s2_bottom/(sigProf.total_s2 - sigProf.s2_bottom)))<<endl;
  cout<<Form("total_s2/total_s1: %3.2f", sigProf.total_s2/sigProf.total_s1)<<endl;
}

/*
void XYRecHelper::PrintComponents (const double * x) const {
  Double_t par[]={x[0], x[1]};
  cout<<endl;
  for(Int_t ich=0; ich<PMTGeom::N_CHANNELS; ich++){
      if (ich < PMTGeom::N_CHANNELS/2) continue; // use only top PMTs
//      if (ich >= PMTGeom::N_CHANNELS/2) continue; // use only bottom PMTs
//      if(XYRecHelper::IsOuterCornerPmt(ich)) continue;

      Double_t expS2overS2total = LikelihoodFunction::GetExpectedS2overS2total(par, ich); // Should we make histogram of solid angle for lookup?
      Double_t expS2overS2totalErr = LikelihoodFunction::GetExpectedS2overS2totalErr(par, ich); // Should we make histogram of solid angle for lookup?
      // solid angle is used only in the first iteration. After 1st, response function should be used.
      Double_t s2 =fSigProf->s2[ich];// -channel.pulses[1].param.fixed_int2 / channel.pmt.spe_mean;
      Double_t chi = (s2/fSigProf->total_s2 - expS2overS2total)/expS2overS2totalErr;
      cout<<Form("ch%2d  s2/s2tot: %4.3f  exps2/s2tot: %4.3f  exps2/s2toterr: %4.3f  chi: %4.3f", ich, s2/fSigProf->total_s2, expS2overS2total, expS2overS2totalErr, chi)<<endl;
  }
  Double_t expS2overS1total = LikelihoodFunction::GetExpectedS2overS1total(par);
  Double_t expS2overS1totalErr = LikelihoodFunction::GetExpectedS2overS1totalErr(par);
  Double_t chi = (fSigProf->s2_over_s1_corr - expS2overS1total)/expS2overS1totalErr;
  cout<<Form("sum   s2/s1_corr: %4.3f  exps2/s1_corr: %4.3f  exps2/s1_correrr: %4.3f  chi: %4.3f", fSigProf->s2_over_s1_corr, expS2overS1total, expS2overS1totalErr, chi)<<endl;
}
*/
