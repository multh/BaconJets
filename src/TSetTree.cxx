// TDataGet class implementation
// system headers
#include<iostream>
#include<fstream>
#include <sys/stat.h>
using namespace std;

// ROOT headers
#include<TCollection.h>
#include<TMath.h>
#include<TFile.h>
#include "TF1.h"
#include "TH1.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TString.h"
#include<TRandom3.h>
#include<TSystem.h>
#include<TVector3.h>
#include "TLine.h"
#include "TCanvas.h"
#include "TMatrixDSym.h"
#include "TMatrixFSym.h"
#include "TTree.h"
// my includes
#include "UHH2/BaconJets/include/TSetTree.h"
TSetTree::TSetTree()
{
file();
}
void TSetTree::file() {

    f = new TFile("new_tree.root", "recreate");
    f -> cd();
    t = new TTree("AnalysisTree", "AnalysisTree");
    t ->SetDirectory(f);


    t->Branch("Runnr",   &Runnr,   "Runnr/I");
    t->Branch("Eventnr", &Eventnr, "Eventnr/I");
    t->Branch("Pt_ave",  &Pt_ave,  "Pt_ave/F");
    t->Branch("Eta_probjt", &Eta_probjt, "Eta_probjt/F");
    t->Branch("Eta_barrjt", &Eta_barrjt, "Eta_barrjt/F");
    t->Branch("Phi_probjt", &Phi_probjt, "Phi_probjt/F");
    t->Branch("Phi_barrjt", &Phi_barrjt, "Phi_barrjt/F");
    t->Branch("Pt_probjt",  &Pt_probjt,  "Pt_probjt/F");
    t->Branch("Pt_barrjt",  &Pt_barrjt,  "Pt_barrjt/F");
    t->Branch("Pt_jt1",  &Pt_jt1,  "Pt_jt1/F");
    t->Branch("Pt_jt2",  &Pt_jt2,  "Pt_jt2/F");
    t->Branch("Phi_jt1",  &Phi_jt1,  "Phi_jt1/F");
    t->Branch("Phi_jt2",  &Phi_jt2,  "Phi_jt2/F");
    t->Branch("Eta_jt1",  &Eta_jt1,  "Eta_jt1/F");
    t->Branch("Eta_jt2",  &Eta_jt2,  "Eta_jt2/F");

    t->Branch("PtRaw_probjt",  &PtRaw_probjt,  "PtRaw_probjt/F");
    t->Branch("PtRaw_barrjt",  &PtRaw_barrjt,  "PtRaw_barrjt/F");
    t->Branch("PtRaw_jt1",  &PtRaw_jt1,  "PtRaw_jt1/F");
    t->Branch("PtRaw_jt2",  &PtRaw_jt2,  "PtRaw_jt2/F");
}

void TSetTree::fillTree(Int_t ev, Int_t rn, Float_t pt_ave, Float_t eta_jt1, Float_t eta_jt2, Float_t phi_jt1, Float_t phi_jt2, Float_t pt_jt1, Float_t pt_jt2, Float_t ptjt1, Float_t ptjt2, Float_t phijt1, Float_t phijt2, Float_t etajt1, Float_t etajt2, Float_t ptraw_jt1, Float_t ptraw_jt2, Float_t ptraw_probjt, Float_t ptraw_barrjt, Float_t alpha_t) {
//     file();

  //  for(event=1; event<=1000; ++event){
        Runnr = rn;
        Eventnr = ev;
        Pt_ave = pt_ave;
        Eta_probjt = eta_jt1;
        Eta_barrjt = eta_jt2;
        Phi_probjt = phi_jt1;
        Phi_barrjt = phi_jt2;
        Pt_probjt = pt_jt1;
        Pt_barrjt = pt_jt2;
        Pt_jt1 = ptjt1;
        Pt_jt2 = ptjt2;
        Phi_jt1 = phijt1;
        Phi_jt2 = phijt2;
        Eta_jt1 = etajt1;
        Eta_jt2 = etajt2;
        PtRaw_probjt = ptraw_probjt;
        PtRaw_barrjt = ptraw_barrjt;
        PtRaw_jt1 = ptraw_jt1;
        PtRaw_jt2 = ptraw_jt2;
	alpha = alpha_t;
        t->Fill();
   // }

//     return 0;
}
void TSetTree::general() {
    f->Write();
}
