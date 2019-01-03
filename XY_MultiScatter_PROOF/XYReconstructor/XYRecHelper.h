#ifndef XYRecHelper_h
#define XYRecHelper_h

#include <iostream>

#include "TROOT.h"
#include "TChain.h"
#include "TSelector.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TSpline.h"
#include "TApplication.h"
#include "TRint.h"
#include "TTree.h"
#include "TF1.h"


//#include "/ds50/app/user/xxiang89/work/xy4slad_ms/XYReconstructor/SignalProfile.h"
//#include "/ds50/app/user/xxiang89/work/xy4slad_ms/XYReconstructor/XYReconstructor.h"
//#include "/ds50/app/user/xxiang89/work/xy4slad_ms/XYReconstructor/PMTGeom.h"
#include "SignalProfile.h"
#include "XYReconstructor.h"
#include "PMTGeom.h"

using namespace std;

class THStack;

class XYRecHelper {
protected:
        Bool_t Debug;
private :
        XYReconstructor *fMyFcn;
	TRint* fRint;    //!pointer to the TRint
	TF2 *gf2D;

public:

   XYRecHelper(XYReconstructor *func, TRint *rint): Debug(false), fMyFcn(func), fRint(rint) {
     //TF2(const char* name, void* ptr, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Int_t npar, const char* className)
     gf2D = new TF2("gf2D", fMyFcn, -1.2*PMTGeom::tpcRadi, 1.2*PMTGeom::tpcRadi, -1.2*PMTGeom::tpcRadi, 1.2*PMTGeom::tpcRadi, 0, "XYReconstructor");
     gf2D->SetNpx(180);
     gf2D->SetNpy(180);
   }

   virtual ~XYRecHelper() { }

   void SetDebug(Bool_t debug){Debug = debug;}
   void SetTRint(TRint *rint){fRint = rint;}
   void SetFunc(XYReconstructor *func){fMyFcn = func;}

   void Draw2DFunc(Double_t x1, Double_t x2, SignalProfile sigProf);
   THStack* DrawChi2alongR(Double_t theta, SignalProfile sigProf);
   THStack* DrawChi2alongTheta(Double_t r, SignalProfile sigProf);
   void PrintS2(SignalProfile sigProf);

   ClassDef(XYRecHelper,0);
};
#endif


