#ifndef PMTGeom_h
#define PMTGeom_h

#include "TVector2.h"
#include "TString.h"
#include "TMath.h"

//struct PMTGeom {
namespace PMTGeom {

//   PMTGeom():
//     N_CHANNELS(38), outerPmtCoverage(0.71), pmtUnit(4.1), tpcRadi(17.945)
//     { }
//   virtual ~PMTGeom() { }

   static Double_t GetDist2fromPMT(Int_t ch, Double_t *par);
   static Double_t GetSignedDistfromPMT(Int_t ch, Double_t *par);
   static void GetPhiFromPmttoCenter(Int_t ch, Double_t *par, Double_t &r, Double_t &phi);

   static Bool_t IsOuterCornerPmt(Int_t ch);
   static Bool_t IsOuterSidePmt(Int_t ch);
   static Int_t GetPmtPositionCategory(Int_t ich);

   const Int_t N_CHANNELS(38);

   const Double_t outerPmtCoverage(0.71);// = 0.71;//update considering shrinking of Tefron by 2%.  0.74; //in cm^2 inefficiency due to geometry
   const Double_t pmtUnit(4.1);// = 4.1;//4.018;//4.1; //in cm   center to center distance 8.2cm /2
   const Double_t tpcRadi(17.945);// = 17.945;//35.89/2.; //in cm   center to center distance 8.2cm /2, 35.17 is 0.98% of 35.89
//   const Int_t N_CHANNELS;
//
//   const Double_t outerPmtCoverage;// = 0.71;//update considering shrinking of Tefron by 2%.  0.74; //in cm^2 inefficiency due to geometry
//   const Double_t pmtUnit;// = 4.1;//4.018;//4.1; //in cm   center to center distance 8.2cm /2
//   const Double_t tpcRadi;// = 17.945;//35.89/2.; //in cm   center to center distance 8.2cm /2, 35.17 is 0.98% of 35.89
   //const Double_t tpcRadi = 17.585;//35.17/2.;//35.89/2.; //in cm   center to center distance 8.2cm /2, 35.17 is 0.98% of 35.89
  //                                     0,    1,    2,    3,    4,    5,    6,  7,   8,   9,  10, 11, 12,   13, 14,  15,  16,  17,  18,   19,   20, 21,  22,  23,  24, 25,   26,   27,   28,   29,30,  31,  32,  33,  34,35,   36, 37
   const Double_t pmt_x[N_CHANNELS] = {  0.,-1.73,-3.46,-3.46,-3.46,-1.73,-1.73, 0.,1.73,3.46,1.73, 0., 0.,-1.73,  0,1.73,1.73,3.46,3.46,-3.46,-1.73, 0.,1.73,3.46,1.73, 0.,-1.73,-3.46,-3.46,-1.73,0.,1.73,3.46,3.46,1.73,0.,-1.73, 0. };
   const Double_t pmt_y[N_CHANNELS] = { -4.,  -3.,  -2.,   0.,   2.,   1.,  -1.,-2., -3., -2., -1., 0., 2.,   3., 4.,  3.,  1.,  0.,  2.,  -2.,  -3.,-4., -3., -2., -1.,-2.,  -1.,   0.,   2.,   1.,0.,  1.,  0.,  2.,  3.,2.,   3., 4. };
   //const Double_t pmt_z[N_CHANNELS] = { -1.,  -1.,  -1.,  -1.,  -1.,  -1.,  -1.,-1., -1., -1., -1.,-1.,-1.,  -1.,-1., -1., -1., -1., -1.,   1.,   1., 1.,  1.,  1.,  1., 1.,   1.,   1.,   1.,   1.,1.,  1.,  1.,  1.,  1.,1.,   1., 1. };

   //                                     0,    1,    2,    3,    4,    5,    6,    7,   8,   9,  10,   11,   12,   13,   14,   15,   16,   17,   18,   19,   20,   21,   22,   23,   24,   25,   26,   27,   28,   29,   30,   31,   32,   33,   34,   35,   36,   37
   //const Double_t pmt_qe[N_CHANNELS] = {  32.7,33.4,35.6,34.5, 32.5, 33.7, 34.2, 35.2,32.5,33.6,35.7, 36.8, 34.1, 33.6, 34.7, 31.4, 36.2, 35.5, 33.9, 37.5, 35.5, 34.8, 33.3, 32.9, 33.3, 34.8, 32.2, 35.4, 35.7, 34.7, 35.3, 34.5, 33.9, 31.3, 33.3, 33.4, 32.8, 33.6 };
   //quantum efficiency normalized to ch30
   const Double_t pmt_qe_norm[N_CHANNELS] = {  1.07951, 1.05689, 0.991573, 1.02319, 1.08615, 1.04748, 1.03216, 1.00284, 1.08615, 1.0506, 0.988796, 0.959239, 1.03519, 1.0506, 1.01729, 1.1242, 0.975138, 0.994366, 1.0413,
                                              0.941333, 0.994366, 1.01437, 1.06006, 1.07295, 1.06006, 1.01437, 1.09627, 0.997175, 0.988796, 1.01729, 1, 1.02319, 1.0413, 1.1278, 1.06006, 1.05689, 1.07622, 1.0506};
   const Double_t charge_normalization[] = { 1.21913, 1.18904, 1.40944, 1.16448, 0.978805, 0.973365, 0.965905, 0.964829, 1.21574, 1.15364, 0.965899, 1, 1.05215, 1.01652, 1.0389, 1.0565, 1.0354, 1.00863, 0.934155 }; // 3rd iteration, obtain by myself with normalization
//   const Double_t charge_normalization[] = { 1.06991, 1.21869, 1.14428, 1.21375, 0.949125, 0.990893, 0.9562, 1.01302, 1.27111, 1.1088, 1.00783, 1, 1.08385, 1.00449, 0.839494, 1.06493, 1.05717, 0.991368, 0.810384 }; // 3rd iteration, obtain by myself with normalization in defferent region

   const Int_t nPmtCategories = 4;
   const TString pmtlabels[nPmtCategories] = {"CENTER", "MIDDLE", "SIDE", "CORNER"};
   enum PmtPosiCategory_t {CENTER, MIDDLE, SIDE, CORNER};

   const int colors[19] = {kBlack,
                     kRed, kRed+3, kPink+1,
                     kMagenta, kViolet+1, kViolet+2,
                     kBlue, kBlue-9, kAzure+10,
                     kAzure+4, kCyan+1, kCyan+3,
                     kTeal,kGreen,kGreen+2,
                     kYellow,kOrange-3,kOrange+9};
};



Double_t PMTGeom::GetDist2fromPMT(Int_t ch, Double_t *par){  //par are supposed to be x and y
  if(par[0]!=par[0] || par[1]!=par[1]) return 999.;
  Double_t s2_x(par[0]), s2_y(par[1]); //x-y position of S2
  Double_t pmtx(pmtUnit*pmt_x[ch]), pmty(pmtUnit*pmt_y[ch]); //x-y position of pmt center
  Double_t rho2 = (s2_x-pmtx)*(s2_x-pmtx) + (s2_y-pmty)*(s2_y-pmty); //square of distance between pmt center to S2 signal position in x-y
  return rho2;
}

Double_t PMTGeom::GetSignedDistfromPMT(Int_t ch, Double_t *par){  //par are supposed to be x and y
  if(par[0]!=par[0] || par[1]!=par[1]) return 999.;
  TVector2 vPmt(pmtUnit*pmt_x[ch], pmtUnit*pmt_y[ch]);
  TVector2 vEvtFromPmt(par[0]-pmtUnit*pmt_x[ch], par[1]-pmtUnit*pmt_y[ch]);
  Double_t phi = vPmt.DeltaPhi(vEvtFromPmt);
  Double_t sign = (TMath::Abs(phi)>TMath::Pi()/2.)? 1. : -1.;

  Double_t rho = sign*vEvtFromPmt.Mod(); //square of distance between pmt center to S2 signal position in x-y
  return rho;
}

void PMTGeom::GetPhiFromPmttoCenter(Int_t ch, Double_t *par, Double_t &r, Double_t &phi){  //par are supposed to be x and y
  if(par[0]!=par[0] || par[1]!=par[1]) return;
  TVector2 vPmt(pmtUnit*pmt_x[ch], pmtUnit*pmt_y[ch]);
  TVector2 vEvtFromPmt(par[0]-pmtUnit*pmt_x[ch], par[1]-pmtUnit*pmt_y[ch]);
  r = vEvtFromPmt.Mod();
  phi = vPmt.DeltaPhi(vEvtFromPmt);
}

Bool_t PMTGeom::IsOuterCornerPmt(Int_t ch){
  Bool_t topPmt = ch==19 || ch==21 || ch==23 || ch==28 || ch==33 || ch==37 ;
//  Boot_t bottomPmt = ch==0 || ch==2 || ch==4 || ch==9 || ch==14 || ch==18 ;
  return topPmt;
}

Bool_t PMTGeom::IsOuterSidePmt(Int_t ch){
  Bool_t topPmt = ch==20 || ch==22 || ch==27 || ch==32 || ch==34 || ch==36 ;
  return topPmt;
}

Int_t PMTGeom::GetPmtPositionCategory(Int_t ich){
//  enum PmtPosiCategory_t {CENTER, MIDDLE, SIDE, CORNER};
  if(ich == 30) return CENTER;
  else if(ich == 24 || ich == 25 || ich == 26 || ich == 29 || ich == 31 || ich == 35) return MIDDLE;
  else if(ich == 20 || ich == 22 || ich == 27 || ich == 32 || ich == 34 || ich == 36) return SIDE;
  else return CORNER;
}
#endif



