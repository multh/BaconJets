//TSetTree.h -class interface
/* #ifndef TDATAGET_H_ */
/* #define TDATAGET_H_ */

#pragma once

#include <TROOT.h>
#include <vector>
#include <TH1F.h>
using namespace std;
//namespace uhh2 {
class TSetTree {  //class declaration
    public:

        TSetTree(); //Default constructor
        int Runnr, Eventnr;
        float Pt_ave, Eta_probjt, Eta_barrjt, Phi_probjt,  Phi_barrjt, Pt_probjt,  Pt_barrjt, Pt_jt1, Pt_jt2, Phi_jt1, Phi_jt2, Eta_jt1, Eta_jt2;
        float PtRaw_probjt,  PtRaw_barrjt, PtRaw_jt1, PtRaw_jt2;
	float alpha;
        TTree * t;
        TFile * f;
        void fillTree(Int_t ev, Int_t rn, Float_t pt_ave, Float_t eta_jt1, Float_t eta_jt2, Float_t phi_jt1, Float_t phi_jt2, Float_t pt_jt1, Float_t pt_jt2, Float_t ptjt1, Float_t ptjt2, Float_t phijt1, Float_t phijt2, Float_t etajt1, Float_t etajt2, Float_t ptraw_jt1, Float_t ptraw_jt2, Float_t ptraw_probjt, Float_t ptraw_barrjt,  Float_t alpha);
        void file();
        void general();

};
//}
//#endif
