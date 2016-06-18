// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// this macro creates pile-up histogram for MC (to be used within UHH2 analysis modul)
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "header.h"
#include <iostream>
#include <cmath>
#include "TString.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TMath.h"
#include "TPaveStats.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"
using namespace std;

void PileUp_MC(TString path="/nfs/dust/cms/user/karavdia/JEC_80X_standalone/TEST/", TString file=""){
  //  path="/pnfs/desy.de/cms/tier2/store/user/akaravdi/RunII_80X_25ns_MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/crab_QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/160615_100441/0000/";//pythia8, FLAT
  //  file="Ntuple_41.root";//pythia8
  //Inout files: root files from UHH2 Ntuplewriter
  file="Sum_Pythia_Flat.root";//pythia8
  //  file="Sum_Herwigpp_Flat.root";//herwigpp
  TFile* MCfile = new TFile(path+file,"READ");
  TTreeReader myReader_MC("AnalysisTree", MCfile);
  //  MCfile->Print();
  TTreeReaderValue<Float_t> nPU(myReader_MC, "m_pileup_TrueNumInteractions");
  TTreeReaderValue<std::vector<float>> weight_mc(myReader_MC, "m_weights");
  TFile* outputfile = new TFile("PileUP_MC.root","RECREATE");
  outputfile->mkdir("input_Event"); 
  TH1D *h_pileup = new TH1D("N_TrueInteractions",";TrueNumInteractions",60,0,60);
   while (myReader_MC.Next()) {
     float sum_weight=0;
     for(unsigned int i=0;i<weight_mc->size();i++)
       sum_weight += weight_mc->at(i);
     h_pileup->Fill(*nPU,sum_weight);
     //     std::cout<<" nPU = "<<*nPU<<std::endl;
   }
   outputfile->cd("input_Event"); 
  h_pileup->Write();
  outputfile->Write(); 
  outputfile->Close(); 
}
