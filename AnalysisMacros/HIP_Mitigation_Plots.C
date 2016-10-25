// plot control plots of important parameters for L2Res determination

#include "header.h"


TString ToStringC(int num) {
  ostringstream start;
  start<<num;
  TString start1=start.str();
  return start1;
}

void HIP_Mitigation_Plots(TString path, TString MCFile_HIP, TString MCFile_mtOff, TString DATA_HIP, TString DATA_mtOff, TString jettag, TString AssumedLumis[3], double alpha_cut){

  //Get Files
  TFile* f_HIP_Lumi1 = new TFile(path+MCFile_HIP+AssumedLumis[0]+".root","READ");
  TFile* f_HIP_Lumi2 = new TFile(path+MCFile_HIP+AssumedLumis[1]+".root","READ");
  TFile* f_HIP_Lumi3 = new TFile(path+MCFile_HIP+AssumedLumis[2]+".root","READ");
  TFile* f_mtOff_Lumi1 = new TFile(path+MCFile_mtOff+AssumedLumis[0]+".root","READ");
  TFile* f_mtOff_Lumi2 = new TFile(path+MCFile_mtOff+AssumedLumis[1]+".root","READ");
  TFile* f_mtOff_Lumi3 = new TFile(path+MCFile_mtOff+AssumedLumis[2]+".root","READ");
  TFile* f_HIP_data = new TFile(path+DATA_HIP+".root","READ");
  TFile* f_mtOff_data = new TFile(path+DATA_mtOff+".root","READ");
  TString dirName = "Selection";
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);


  
  //Get histos, that are already stored
  TH1F *eta_jets_HIP1,*eta_jets_HIP2, *eta_jets_HIP3, *eta_jet_mtOff1, *eta_jets_mtOff2, *eta_jets_mtOff3;
  TH1F *pt_jets_HIP1,*pt_jets_HIP2, *pt_jets_HIP3, *pt_jets_mtOff1, *pt_jets_mtOff2, *pt_jets_mtOff3;
  /* TH2F *mpf_eta_HIP1,*mpf_eta_HIP2, *mpf_eta_HIP3, *mpf_eta_mtOff1, *mpf_eta_mtOff2, *mpf_eta_mtOff3;
     TH2F *rrel_eta_HIP1,*rrel_eta_HIP2, *rrel_eta_HIP3, *rrel_eta_mtOff1, *rrel_eta_mtOff2, *rrel_eta_mtOff3;*/
  
  eta_jets_HIP1 = (TH1F*)f_HIP_Lumi1->Get(dirName+"/eta");
  eta_jets_HIP2 = (TH1F*)f_HIP_Lumi2->Get(dirName+"/eta");
  eta_jets_HIP3 = (TH1F*)f_HIP_Lumi3->Get(dirName+"/eta");
  eta_jets_mtOff1 = (TH1F*)f_mtOff_Lumi1->Get(dirName+"/eta");
  eta_jets_mtOff2 = (TH1F*)f_mtOff_Lumi2->Get(dirName+"/eta");
  eta_jets_mtOff3 = (TH1F*)f_mtOff_Lumi3->Get(dirName+"/eta");

  pt_jets_HIP1 = (TH1F*)f_HIP_Lumi1->Get(dirName+"/pt");
  pt_jets_HIP2 = (TH1F*)f_HIP_Lumi2->Get(dirName+"/pt");
  pt_jets_HIP3 = (TH1F*)f_HIP_Lumi3->Get(dirName+"/pt");
  pt_jets_mtOff1 = (TH1F*)f_mtOff_Lumi1->Get(dirName+"/pt");
  pt_jets_mtOff2 = (TH1F*)f_mtOff_Lumi2->Get(dirName+"/pt");
  pt_jets_mtOff3 = (TH1F*)f_mtOff_Lumi3->Get(dirName+"/pt");

  /*
  mpf_eta_HIP1 = (TH2F*)f_HIP_Lumi1->Get(dirName+"/mpf_vs_etaProbe");
  mpf_eta_HIP2 = (TH2F*)f_HIP_Lumi2->Get(dirName+"/mpf_vs_etaProbe");
  mpf_eta_HIP3 = (TH2F*)f_HIP_Lumi3->Get(dirName+"/mpf_vs_etaProbe");
  mpf_eta_mtOff1 = (TH2F*)f_mtOff_Lumi1->Get(dirName+"/mpf_vs_etaProbe");
  mpf_eta_mtOff2 = (TH2F*)f_mtOff_Lumi2->Get(dirName+"/mpf_vs_etaProbe");
  mpf_eta_mtOff3 = (TH2F*)f_mtOff_Lumi3->Get(dirName+"/mpf_vs_etaProbe");

  rrel_eta_HIP1 = (TH2F*)f_HIP_Lumi1->Get(dirName+"/r_rel_vs_etaProbe");
  rrel_eta_HIP2 = (TH2F*)f_HIP_Lumi2->Get(dirName+"/r_rel_vs_etaProbe");
  rrel_eta_HIP3 = (TH2F*)f_HIP_Lumi3->Get(dirName+"/r_rel_vs_etaProbe");
  rrel_eta_mtOff1 = (TH2F*)f_mtOff_Lumi1->Get(dirName+"/r_rel_vs_etaProbe");
  rrel_eta_mtOff2 = (TH2F*)f_mtOff_Lumi2->Get(dirName+"/r_rel_vs_etaProbe");
  rrel_eta_mtOff3 = (TH2F*)f_mtOff_Lumi3->Get(dirName+"/r_rel_vs_etaProbe");
*/

  //Set up new histos
  TH2F *mpf_eta_HIP1 = new TH2F("mpf_eta_HIP1","mtOn, L_{inst} = 0.6;|#eta| probe;MPF response",10,0,5,80,0,4);
  TH2F *mpf_eta_HIP2 = new TH2F("mpf_eta_HIP2","mtOn, L_{inst} = 0.8;|#eta| probe;MPF response",10,0,5,80,0,4);
  TH2F *mpf_eta_HIP3 = new TH2F("mpf_eta_HIP3","mtOn, L_{inst} = 1.2;|#eta| probe;MPF response",10,0,5,80,0,4);
  TH2F *mpf_eta_mtOff1 = new TH2F("mpf_eta_mtOff1","mtOff, L_{inst} = 0.6;|#eta| probe;MPF response",10,0,5,80,0,4);
  TH2F *mpf_eta_mtOff2 = new TH2F("mpf_eta_mtOff2","mtOff, L_{inst} = 0.8;|#eta| probe;MPF response",10,0,5,80,0,4);
  TH2F *mpf_eta_mtOff3 = new TH2F("mpf_eta_mtOff3","mtOff, L_{inst} = 1.2;|#eta| probe;MPF response",10,0,5,80,0,4);
  TH2F *mpf_eta_mtOff_data = new TH2F("mpf_eta_mtOff_data","mtOff;|#eta| probe;MPF response",10,0,5,80,0,4);
  TH2F *mpf_eta_HIP_data = new TH2F("mpf_eta_HIP_data","mtOn;|#eta| probe;MPF response",10,0,5,80,0,4);

  TH2F *rrel_eta_HIP1 = new TH2F("rrel_eta_HIP1","mtOn, L_{inst} = 0.6;|#eta| probe;p_{T} balance response",10,0,5,80,0,4);
  TH2F *rrel_eta_HIP2 = new TH2F("rrel_eta_HIP2","mtOn, L_{inst} = 0.8;|#eta| probe;p_{T} balance response",10,0,5,80,0,4);
  TH2F *rrel_eta_HIP3 = new TH2F("rrel_eta_HIP3","mtOn, L_{inst} = 1.2;|#eta| probe;p_{T} balance response",10,0,5,80,0,4);
  TH2F *rrel_eta_mtOff1 = new TH2F("rrel_eta_mtOff1","mtOff, L_{inst} = 0.6;|#eta| probe;p_{T} balance response",10,0,5,80,0,4);
  TH2F *rrel_eta_mtOff2 = new TH2F("rrel_eta_mtOff2","mtOff, L_{inst} = 0.8;|#eta| probe;p_{T} balance response",10,0,5,80,0,4);
  TH2F *rrel_eta_mtOff3 = new TH2F("rrel_eta_mtOff3","mtOff, L_{inst} = 1.2;|#eta| probe;p_{T} balance response",10,0,5,80,0,4);
  TH2F *rrel_eta_mtOff_data = new TH2F("rrel_eta_mtOff_data","mtOff;|#eta| probe;p_{T} balance response",10,0,5,80,0,4);
  TH2F *rrel_eta_HIP_data = new TH2F("rrel_eta_HIP_data","mtOn;|#eta| probe;p_{T} balance response",10,0,5,80,0,4);

  TH2F *mpf_nvert_HIP1 = new TH2F("mpf_nvert_HIP1","mtOn, L_{inst} = 0.6;N primary vertices;MPF response",61,-0.5,60.5,80,0,4);
  TH2F *mpf_nvert_HIP2 = new TH2F("mpf_nvert_HIP2","mtOn, L_{inst} = 0.8;N primary vertices;MPF response",61,-0.5,60.5,80,0,4);
  TH2F *mpf_nvert_HIP3 = new TH2F("mpf_nvert_HIP3","mtOn, L_{inst} = 1.2;N primary vertices;MPF response",61,-0.5,60.5,80,0,4);
  TH2F *mpf_nvert_mtOff1 = new TH2F("mpf_nvert_mtOff1","mtOff, L_{inst} = 0.6;N primary vertices;MPF response",61,-0.5,60.5,80,0,4);
  TH2F *mpf_nvert_mtOff2 = new TH2F("mpf_nvert_mtOff2","mtOff, L_{inst} = 0.8;N primary vertices;MPF response",61,-0.5,60.5,80,0,4);
  TH2F *mpf_nvert_mtOff3 = new TH2F("mpf_nvert_mtOff3","mtOff, L_{inst} = 1.2;N primary vertices;MPF response",61,-0.5,60.5,80,0,4);
  TH2F *mpf_nvert_mtOff_data = new TH2F("mpf_nvert_mtOff_data","mtOff;N primary vertices;MPF response",61,-0.5,60.5,80,0,4);
  TH2F *mpf_nvert_HIP_data = new TH2F("mpf_nvert_HIP_data","mtOn;N primary vertices;MPF response",61,-0.5,60.5,80,0,4);

  TH2F *rrel_nvert_HIP1 = new TH2F("rrel_nvert_HIP1","mtOn, L_{inst} = 0.6;N primary vertices;p_{T} balance response",61,-0.5,60.5,80,0,4);
  TH2F *rrel_nvert_HIP2 = new TH2F("rrel_nvert_HIP2","mtOn, L_{inst} = 0.8;N primary vertices;p_{T} balance response",61,-0.5,60.5,80,0,4);
  TH2F *rrel_nvert_HIP3 = new TH2F("rrel_nvert_HIP3","mtOn, L_{inst} = 1.2;N primary vertices;p_{T} balance response",61,-0.5,60.5,80,0,4);
  TH2F *rrel_nvert_mtOff1 = new TH2F("rrel_nvert_mtOff1","mtOff, L_{inst} = 0.6;N primary vertices;p_{T} balance response",61,-0.5,60.5,80,0,4);
  TH2F *rrel_nvert_mtOff2 = new TH2F("rrel_nvert_mtOff2","mtOff, L_{inst} = 0.8;N primary vertices;p_{T} balance response",61,-0.5,60.5,80,0,4);
  TH2F *rrel_nvert_mtOff3 = new TH2F("rrel_nvert_mtOff3","mtOff, L_{inst} = 1.2;N primary vertices;p_{T} balance response",61,-0.5,60.5,80,0,4);
  TH2F *rrel_nvert_mtOff_data = new TH2F("rrel_nvert_mtOff_data","mtOff;N primary vertices;p_{T} balance response",61,-0.5,60.5,80,0,4);
  TH2F *rrel_nvert_HIP_data = new TH2F("rrel_nvert_HIP_data","mtOn;N primary vertices;p_{T} balance response",61,-0.5,60.5,80,0,4);

  TH2F *njets_eta_HIP1 = new TH2F("njets_eta_HIP1","mtOn, L_{inst} = 0.6;|#eta| probe;N_{jets} per event",10,0,5,40,-0.5,39.5);
  TH2F *njets_eta_HIP2 = new TH2F("njets_eta_HIP2","mtOn, L_{inst} = 0.8;|#eta| probe;N_{jets} per event",10,0,5,40,-0.5,39.5);
  TH2F *njets_eta_HIP3 = new TH2F("njets_eta_HIP3","mtOn, L_{inst} = 1.2;|#eta| probe;N_{jets} per event",10,0,5,40,-0.5,39.5);
  TH2F *njets_eta_mtOff1 = new TH2F("njets_eta_mtOff1","mtOff, L_{inst} = 0.6;|#eta| probe;N_{jets} per event",10,0,5,40,-0.5,39.5);
  TH2F *njets_eta_mtOff2 = new TH2F("njets_eta_mtOff2","mtOff, L_{inst} = 0.8;|#eta| probe;N_{jets} per event",10,0,5,40,-0.5,39.5);
  TH2F *njets_eta_mtOff3 = new TH2F("njets_eta_mtOff3","mtOff, L_{inst} = 1.2;|#eta| probe;N_{jets} per event",10,0,5,40,-0.5,39.5);
  TH2F *njets_eta_mtOff_data = new TH2F("njets_eta_mtOff_data","mtOff;|#eta| probe;N_{jets} per event",10,0,5,40,-0.5,39.5);
  TH2F *njets_eta_HIP_data = new TH2F("njets_eta_HIPdata","mtOn;|#eta| probe;N_{jets} per event",10,0,5,40,-0.5,39.5);

  TH2F *njets_pt_HIP1 = new TH2F("njets_pt_HIP1","mtOn, L_{inst} = 0.6;p_{T} probe (GeV);N_{jets} per event",20,0,1500,40,-0.5,39.5);
  TH2F *njets_pt_HIP2 = new TH2F("njets_pt_HIP2","mtOn, L_{inst} = 0.8;p_{T} probe (GeV);N_{jets} per event",20,0,1500,40,-0.5,39.5);
  TH2F *njets_pt_HIP3 = new TH2F("njets_pt_HIP3","mtOn, L_{inst} = 1.2;p_{T} probe (GeV);N_{jets} per event",20,0,1500,40,-0.5,39.5);
  TH2F *njets_pt_mtOff1 = new TH2F("njets_pt_mtOff1","mtOff, L_{inst} = 0.6;p_{T} probe (GeV);N_{jets} per event",20,0,1500,40,-0.5,39.5);
  TH2F *njets_pt_mtOff2 = new TH2F("njets_pt_mtOff2","mtOff, L_{inst} = 0.8;p_{T} probe (GeV);N_{jets} per event",20,0,1500,40,-0.5,39.5);
  TH2F *njets_pt_mtOff3 = new TH2F("njets_pt_mtOff3","mtOff, L_{inst} = 1.2;p_{T} probe (GeV);N_{jets} per event",20,0,1500,40,-0.5,39.5);
  TH2F *njets_pt_mtOff_data = new TH2F("njets_pt_mtOff_data","mtOff;p_{T} probe (GeV);N_{jets} per event",20,0,1500,40,-0.5,39.5);
  TH2F *njets_pt_HIP_data = new TH2F("njets_pt_HIP_data","mtOn;p_{T} probe (GeV);N_{jets} per event",20,0,1500,40,-0.5,39.5);



  //Read out TTrees and apply the alpha-cut before filling histos
 //Get values from TTree
    TTreeReader myReader_HIP_Lumi1("AnalysisTree", f_HIP_Lumi1);
    TTreeReader myReader_HIP_Lumi2("AnalysisTree", f_HIP_Lumi2);
    TTreeReader myReader_HIP_Lumi3("AnalysisTree", f_HIP_Lumi3);
    TTreeReader myReader_mtOff_Lumi1("AnalysisTree", f_mtOff_Lumi1);
    TTreeReader myReader_mtOff_Lumi2("AnalysisTree", f_mtOff_Lumi2);
    TTreeReader myReader_mtOff_Lumi3("AnalysisTree", f_mtOff_Lumi3);
    TTreeReader myReader_mtOff_data("AnalysisTree", f_mtOff_data);
    TTreeReader myReader_HIP_data("AnalysisTree", f_HIP_data);

    //MC
    TTreeReaderValue<Float_t> probejet_eta_HIP_Lumi1(myReader_HIP_Lumi1, "probejet_eta");
    TTreeReaderValue<Float_t> probejet_pt_HIP_Lumi1(myReader_HIP_Lumi1, "probejet_pt");
    TTreeReaderValue<Float_t> alpha_HIP_Lumi1(myReader_HIP_Lumi1, "alpha");
    TTreeReaderValue<Float_t> rel_r_HIP_Lumi1(myReader_HIP_Lumi1, "rel_r");
    TTreeReaderValue<Float_t> mpf_r_HIP_Lumi1(myReader_HIP_Lumi1, "mpf_r");
    TTreeReaderValue<Int_t> n_jets_HIP_Lumi1(myReader_HIP_Lumi1, "Njet");
    TTreeReaderValue<Int_t> n_vert_HIP_Lumi1(myReader_HIP_Lumi1, "nGoodvertices");

    TTreeReaderValue<Float_t> probejet_eta_HIP_Lumi2(myReader_HIP_Lumi2, "probejet_eta");
    TTreeReaderValue<Float_t> probejet_pt_HIP_Lumi2(myReader_HIP_Lumi2, "probejet_pt");
    TTreeReaderValue<Float_t> alpha_HIP_Lumi2(myReader_HIP_Lumi2, "alpha");
    TTreeReaderValue<Float_t> rel_r_HIP_Lumi2(myReader_HIP_Lumi2, "rel_r");
    TTreeReaderValue<Float_t> mpf_r_HIP_Lumi2(myReader_HIP_Lumi2, "mpf_r");
    TTreeReaderValue<Int_t> n_jets_HIP_Lumi2(myReader_HIP_Lumi2, "Njet");
    TTreeReaderValue<Int_t> n_vert_HIP_Lumi2(myReader_HIP_Lumi2, "nGoodvertices");

    TTreeReaderValue<Float_t> probejet_eta_HIP_Lumi3(myReader_HIP_Lumi3, "probejet_eta");
    TTreeReaderValue<Float_t> probejet_pt_HIP_Lumi3(myReader_HIP_Lumi3, "probejet_pt");
    TTreeReaderValue<Float_t> alpha_HIP_Lumi3(myReader_HIP_Lumi3, "alpha");
    TTreeReaderValue<Float_t> rel_r_HIP_Lumi3(myReader_HIP_Lumi3, "rel_r");
    TTreeReaderValue<Float_t> mpf_r_HIP_Lumi3(myReader_HIP_Lumi3, "mpf_r");
    TTreeReaderValue<Int_t> n_jets_HIP_Lumi3(myReader_HIP_Lumi3, "Njet");
    TTreeReaderValue<Int_t> n_vert_HIP_Lumi3(myReader_HIP_Lumi3, "nGoodvertices");
   
    TTreeReaderValue<Float_t> probejet_eta_mtOff_Lumi1(myReader_mtOff_Lumi1, "probejet_eta");
    TTreeReaderValue<Float_t> probejet_pt_mtOff_Lumi1(myReader_mtOff_Lumi1, "probejet_pt");
    TTreeReaderValue<Float_t> alpha_mtOff_Lumi1(myReader_mtOff_Lumi1, "alpha");
    TTreeReaderValue<Float_t> rel_r_mtOff_Lumi1(myReader_mtOff_Lumi1, "rel_r");
    TTreeReaderValue<Float_t> mpf_r_mtOff_Lumi1(myReader_mtOff_Lumi1, "mpf_r");
    TTreeReaderValue<Int_t> n_jets_mtOff_Lumi1(myReader_mtOff_Lumi1, "Njet");
    TTreeReaderValue<Int_t> n_vert_mtOff_Lumi1(myReader_mtOff_Lumi1, "nGoodvertices");

    TTreeReaderValue<Float_t> probejet_eta_mtOff_Lumi2(myReader_mtOff_Lumi2, "probejet_eta");
    TTreeReaderValue<Float_t> probejet_pt_mtOff_Lumi2(myReader_mtOff_Lumi2, "probejet_pt");
    TTreeReaderValue<Float_t> alpha_mtOff_Lumi2(myReader_mtOff_Lumi2, "alpha");
    TTreeReaderValue<Float_t> rel_r_mtOff_Lumi2(myReader_mtOff_Lumi2, "rel_r");
    TTreeReaderValue<Float_t> mpf_r_mtOff_Lumi2(myReader_mtOff_Lumi2, "mpf_r");
    TTreeReaderValue<Int_t> n_jets_mtOff_Lumi2(myReader_mtOff_Lumi2, "Njet");
    TTreeReaderValue<Int_t> n_vert_mtOff_Lumi2(myReader_mtOff_Lumi2, "nGoodvertices");
   
    TTreeReaderValue<Float_t> probejet_eta_mtOff_Lumi3(myReader_mtOff_Lumi3, "probejet_eta");
    TTreeReaderValue<Float_t> probejet_pt_mtOff_Lumi3(myReader_mtOff_Lumi3, "probejet_pt");
    TTreeReaderValue<Float_t> alpha_mtOff_Lumi3(myReader_mtOff_Lumi3, "alpha");
    TTreeReaderValue<Float_t> rel_r_mtOff_Lumi3(myReader_mtOff_Lumi3, "rel_r");
    TTreeReaderValue<Float_t> mpf_r_mtOff_Lumi3(myReader_mtOff_Lumi3, "mpf_r");
    TTreeReaderValue<Int_t> n_jets_mtOff_Lumi3(myReader_mtOff_Lumi3, "Njet");
    TTreeReaderValue<Int_t> n_vert_mtOff_Lumi3(myReader_mtOff_Lumi3, "nGoodvertices");

    //data
    TTreeReaderValue<Float_t> probejet_eta_mtOff_data(myReader_mtOff_data, "probejet_eta");
    TTreeReaderValue<Float_t> probejet_pt_mtOff_data(myReader_mtOff_data, "probejet_pt");
    TTreeReaderValue<Float_t> alpha_mtOff_data(myReader_mtOff_data, "alpha");
    TTreeReaderValue<Float_t> rel_r_mtOff_data(myReader_mtOff_data, "rel_r");
    TTreeReaderValue<Float_t> mpf_r_mtOff_data(myReader_mtOff_data, "mpf_r");
    TTreeReaderValue<Int_t> n_jets_mtOff_data(myReader_mtOff_data, "Njet");
    TTreeReaderValue<Int_t> n_vert_mtOff_data(myReader_mtOff_data, "nGoodvertices");
       
    TTreeReaderValue<Float_t> probejet_eta_HIP_data(myReader_HIP_data, "probejet_eta");
    TTreeReaderValue<Float_t> probejet_pt_HIP_data(myReader_HIP_data, "probejet_pt");
    TTreeReaderValue<Float_t> alpha_HIP_data(myReader_HIP_data, "alpha");
    TTreeReaderValue<Float_t> rel_r_HIP_data(myReader_HIP_data, "rel_r");
    TTreeReaderValue<Float_t> mpf_r_HIP_data(myReader_HIP_data, "mpf_r");
    TTreeReaderValue<Int_t> n_jets_HIP_data(myReader_HIP_data, "Njet");
    TTreeReaderValue<Int_t> n_vert_HIP_data(myReader_HIP_data, "nGoodvertices");


    //loop over all events to fill histograms
    //MC
    while (myReader_HIP_Lumi1.Next()) {
      if(*alpha_HIP_Lumi1<=0.3) continue; //only fill for alpha > 0.3
      float abs_eta = fabs(*probejet_eta_HIP_Lumi1);
      mpf_eta_HIP1->Fill(abs_eta,*mpf_r_HIP_Lumi1); //MPF - eta
      rrel_eta_HIP1->Fill(abs_eta,*rel_r_HIP_Lumi1); // ptbalance - eta
      mpf_nvert_HIP1->Fill(*n_vert_HIP_Lumi1,*mpf_r_HIP_Lumi1); //MPF - nvert
      rrel_nvert_HIP1->Fill(*n_vert_HIP_Lumi1,*rel_r_HIP_Lumi1); // ptbalance - nvert
      njets_eta_HIP1->Fill(abs_eta,*n_jets_HIP_Lumi1); //Njets - eta
      njets_pt_HIP1->Fill(*probejet_pt_HIP_Lumi1,*n_jets_HIP_Lumi1); //Njets - pt
    }

    while (myReader_HIP_Lumi2.Next()) {
      if(*alpha_HIP_Lumi2<=0.3) continue; //only fill for alpha > 0.3
      float abs_eta = fabs(*probejet_eta_HIP_Lumi2);
      mpf_eta_HIP2->Fill(abs_eta,*mpf_r_HIP_Lumi2); //MPF - eta
      rrel_eta_HIP2->Fill(abs_eta,*rel_r_HIP_Lumi2); // ptbalance - eta
      mpf_nvert_HIP2->Fill(*n_vert_HIP_Lumi2,*mpf_r_HIP_Lumi2); //MPF - nvert
      rrel_nvert_HIP2->Fill(*n_vert_HIP_Lumi2,*rel_r_HIP_Lumi2); // ptbalance - nvert
      njets_eta_HIP2->Fill(abs_eta,*n_jets_HIP_Lumi2); //Njets - eta
      njets_pt_HIP2->Fill(*probejet_pt_HIP_Lumi2,*n_jets_HIP_Lumi2); //Njets - pt
    }

    while (myReader_HIP_Lumi3.Next()) {
      if(*alpha_HIP_Lumi3<=0.3) continue; //only fill for alpha > 0.3
      float abs_eta = fabs(*probejet_eta_HIP_Lumi3);
      mpf_eta_HIP3->Fill(abs_eta,*mpf_r_HIP_Lumi3); //MPF - eta
      rrel_eta_HIP3->Fill(abs_eta,*rel_r_HIP_Lumi3); // ptbalance - eta
      mpf_nvert_HIP3->Fill(*n_vert_HIP_Lumi3,*mpf_r_HIP_Lumi3); //MPF - nvert
      rrel_nvert_HIP3->Fill(*n_vert_HIP_Lumi3,*rel_r_HIP_Lumi3); // ptbalance - nvert
      njets_eta_HIP3->Fill(abs_eta,*n_jets_HIP_Lumi3); //Njets - eta
      njets_pt_HIP3->Fill(*probejet_pt_HIP_Lumi3,*n_jets_HIP_Lumi3); //Njets - pt
    }

    while (myReader_mtOff_Lumi1.Next()) {
      if(*alpha_mtOff_Lumi1<=0.3) continue; //only fill for alpha > 0.3
      float abs_eta = fabs(*probejet_eta_mtOff_Lumi1);
      mpf_eta_mtOff1->Fill(abs_eta,*mpf_r_mtOff_Lumi1); //MPF - eta
      rrel_eta_mtOff1->Fill(abs_eta,*rel_r_mtOff_Lumi1); // ptbalance - eta
      mpf_nvert_mtOff1->Fill(*n_vert_mtOff_Lumi1,*mpf_r_mtOff_Lumi1); //MPF - nvert
      rrel_nvert_mtOff1->Fill(*n_vert_mtOff_Lumi1,*rel_r_mtOff_Lumi1); // ptbalance - nvert
      njets_eta_mtOff1->Fill(abs_eta,*n_jets_mtOff_Lumi1); //Njets - eta
      njets_pt_mtOff1->Fill(*probejet_pt_mtOff_Lumi1,*n_jets_mtOff_Lumi1); //Njets - pt
    }

    while (myReader_mtOff_Lumi2.Next()) {
      if(*alpha_mtOff_Lumi2<=0.3) continue; //only fill for alpha > 0.3
      float abs_eta = fabs(*probejet_eta_mtOff_Lumi2);
      mpf_eta_mtOff2->Fill(abs_eta,*mpf_r_mtOff_Lumi2); //MPF - eta
      rrel_eta_mtOff2->Fill(abs_eta,*rel_r_mtOff_Lumi2); // ptbalance - eta
      mpf_nvert_mtOff2->Fill(*n_vert_mtOff_Lumi2,*mpf_r_mtOff_Lumi2); //MPF - nvert
      rrel_nvert_mtOff2->Fill(*n_vert_mtOff_Lumi2,*rel_r_mtOff_Lumi2); // ptbalance - nvert
      njets_eta_mtOff2->Fill(abs_eta,*n_jets_mtOff_Lumi2); //Njets - eta
      njets_pt_mtOff2->Fill(*probejet_pt_mtOff_Lumi2,*n_jets_mtOff_Lumi2); //Njets - pt
    }

    while (myReader_mtOff_Lumi3.Next()) {
      if(*alpha_mtOff_Lumi3<=0.3) continue; //only fill for alpha > 0.3
      float abs_eta = fabs(*probejet_eta_mtOff_Lumi3);
      mpf_eta_mtOff3->Fill(abs_eta,*mpf_r_mtOff_Lumi3); //MPF - eta
      rrel_eta_mtOff3->Fill(abs_eta,*rel_r_mtOff_Lumi3); // ptbalance - eta
      mpf_nvert_mtOff3->Fill(*n_vert_mtOff_Lumi3,*mpf_r_mtOff_Lumi3); //MPF - nvert
      rrel_nvert_mtOff3->Fill(*n_vert_mtOff_Lumi3,*rel_r_mtOff_Lumi3); // ptbalance - nvert
      njets_eta_mtOff3->Fill(abs_eta,*n_jets_mtOff_Lumi3); //Njets - eta
      njets_pt_mtOff3->Fill(*probejet_pt_mtOff_Lumi3,*n_jets_mtOff_Lumi3); //Njets - pt
    }

    //data
    while (myReader_mtOff_data.Next()) {
      if(*alpha_mtOff_data<=0.3) continue; //only fill for alpha > 0.3
      float abs_eta = fabs(*probejet_eta_mtOff_data);
      mpf_eta_mtOff_data->Fill(abs_eta,*mpf_r_mtOff_data); //MPF - eta
      rrel_eta_mtOff_data->Fill(abs_eta,*rel_r_mtOff_data); // ptbalance - eta
      mpf_nvert_mtOff_data->Fill(*n_vert_mtOff_data,*mpf_r_mtOff_data); //MPF - nvert
      rrel_nvert_mtOff_data->Fill(*n_vert_mtOff_data,*rel_r_mtOff_data); // ptbalance - nvert
      njets_eta_mtOff_data->Fill(abs_eta,*n_jets_mtOff_data); //Njets - eta
      njets_pt_mtOff_data->Fill(*probejet_pt_mtOff_data,*n_jets_mtOff_data); //Njets - pt
    }

    while (myReader_HIP_data.Next()) {
      if(*alpha_HIP_data<=0.3) continue; //only fill for alpha > 0.3
      float abs_eta = fabs(*probejet_eta_HIP_data);
      mpf_eta_HIP_data->Fill(abs_eta,*mpf_r_HIP_data); //MPF - eta
      rrel_eta_HIP_data->Fill(abs_eta,*rel_r_HIP_data); // ptbalance - eta
      mpf_nvert_HIP_data->Fill(*n_vert_HIP_data,*mpf_r_HIP_data); //MPF - nvert
      rrel_nvert_HIP_data->Fill(*n_vert_HIP_data,*rel_r_HIP_data); // ptbalance - nvert
      njets_eta_HIP_data->Fill(abs_eta,*n_jets_HIP_data); //Njets - eta
      njets_pt_HIP_data->Fill(*probejet_pt_HIP_data,*n_jets_HIP_data); //Njets - pt
    }






  //Create beatiful plots
  
  TCanvas* c_eta = new TCanvas();

  //Marker styles and colors
  eta_jets_mtOff1->SetMarkerStyle(24);
  eta_jets_mtOff2->SetMarkerStyle(24);
  eta_jets_mtOff3->SetMarkerStyle(24);
  eta_jets_HIP1->SetMarkerStyle(20);
  eta_jets_HIP2->SetMarkerStyle(20);
  eta_jets_HIP3->SetMarkerStyle(20);

  eta_jets_mtOff1->SetMarkerSize(0.75);
  eta_jets_mtOff2->SetMarkerSize(0.75);
  eta_jets_mtOff3->SetMarkerSize(0.75);
  eta_jets_HIP1->SetMarkerSize(0.75);
  eta_jets_HIP2->SetMarkerSize(0.75);
  eta_jets_HIP3->SetMarkerSize(0.75);

  eta_jets_mtOff1->SetMarkerColor(2);
  eta_jets_HIP1->SetMarkerColor(2);
  eta_jets_mtOff2->SetMarkerColor(3);
  eta_jets_HIP2->SetMarkerColor(3);
  eta_jets_mtOff3->SetMarkerColor(4);
  eta_jets_HIP3->SetMarkerColor(4);

  eta_jets_mtOff1->SetLineColor(2);
  eta_jets_HIP1->SetLineColor(2);
  eta_jets_mtOff2->SetLineColor(3);
  eta_jets_HIP2->SetLineColor(3);
  eta_jets_mtOff3->SetLineColor(4);
  eta_jets_HIP3->SetLineColor(4);

  //Set axis titles
  eta_jets_HIP1->GetXaxis()->SetTitle("#eta_{jet}");
  eta_jets_HIP1->GetYaxis()->SetTitle("N_{jets}");

  //Draw
  eta_jets_HIP1->Draw();
  eta_jets_HIP2->Draw("SAME");
  eta_jets_HIP3->Draw("SAME");
  eta_jets_mtOff1->Draw("SAME");
  eta_jets_mtOff2->Draw("SAME");
  eta_jets_mtOff3->Draw("SAME");

  TLegend *leg1;
  leg1 = new TLegend(0.63,0.63,0.9,0.85,"","brNDC");
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.045);
  leg1->SetFillColor(10);
  leg1->SetLineColor(1);
  leg1->SetTextFont(42);
  leg1->AddEntry(eta_jets_HIP1,"mtOn, L_{inst} = 0.6","lp");
  leg1->AddEntry(eta_jets_mtOff1,"mtOff, L_{inst} = 0.6","lp");
  leg1->AddEntry(eta_jets_HIP2,"mtOn, L_{inst} = 0.8","lp");
  leg1->AddEntry(eta_jets_mtOff2,"mtOff, L_{inst} = 0.8","lp");
  leg1->AddEntry(eta_jets_HIP3,"mtOn, L_{inst} = 1.2","lp");
  leg1->AddEntry(eta_jets_mtOff3,"mtOff, L_{inst} = 1.2","lp");
  leg1->Draw();

  //Only with mitigation
  TCanvas* c_eta_HIP = new TCanvas();
  eta_jets_HIP1->Draw();
  eta_jets_HIP2->Draw("SAME");
  eta_jets_HIP3->Draw("SAME");

  //Legend
  TLegend *leg3;
  leg3 = new TLegend(0.63,0.63,0.9,0.85,"","brNDC");
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.045);
  leg3->SetFillColor(10);
  leg3->SetLineColor(1);
  leg3->SetTextFont(42);
  leg3->AddEntry(eta_jets_HIP1,"mtOn, L_{inst} = 0.6","lp");
  leg3->AddEntry(eta_jets_HIP2,"mtOn, L_{inst} = 0.8","lp");
  leg3->AddEntry(eta_jets_HIP3,"mtOn, L_{inst} = 1.2","lp");
  leg3->Draw();

  //Only without mitigation
  TCanvas* c_eta_mtOff = new TCanvas();
  eta_jets_mtOff1->Draw();
  eta_jets_mtOff2->Draw("SAME");
  eta_jets_mtOff3->Draw("SAME");

  //Legend
  TLegend *leg4;
  leg4 = new TLegend(0.63,0.63,0.9,0.85,"","brNDC");
  leg4->SetBorderSize(0);
  leg4->SetTextSize(0.045);
  leg4->SetFillColor(10);
  leg4->SetLineColor(1);
  leg4->SetTextFont(42);
  leg4->AddEntry(eta_jets_mtOff1,"mtOff, L_{inst} = 0.6","lp");
  leg4->AddEntry(eta_jets_mtOff2,"mtOff, L_{inst} = 0.8","lp");
  leg4->AddEntry(eta_jets_mtOff3,"mtOff, L_{inst} = 1.2","lp");
  leg4->Draw();




  TCanvas* c_eta_norm = new TCanvas();

  //Clone histos
  TH1F* eta_jets_mtOff1_norm = (TH1F*)eta_jets_mtOff1->Clone();
  TH1F* eta_jets_mtOff2_norm = (TH1F*)eta_jets_mtOff2->Clone();
  TH1F* eta_jets_mtOff3_norm = (TH1F*)eta_jets_mtOff3->Clone();
  TH1F* eta_jets_HIP1_norm = (TH1F*)eta_jets_HIP1->Clone();
  TH1F* eta_jets_HIP2_norm = (TH1F*)eta_jets_HIP2->Clone();
  TH1F* eta_jets_HIP3_norm = (TH1F*)eta_jets_HIP3->Clone();

  eta_jets_HIP1_norm->GetYaxis()->SetTitle("N_{jets} normalized");

  //Normalise histos
  eta_jets_mtOff1_norm->Scale(1./eta_jets_mtOff1_norm->Integral());
  eta_jets_mtOff2_norm->Scale(1./eta_jets_mtOff2_norm->Integral());
  eta_jets_mtOff3_norm->Scale(1./eta_jets_mtOff3_norm->Integral());
  eta_jets_HIP1_norm->Scale(1./eta_jets_HIP1_norm->Integral());
  eta_jets_HIP2_norm->Scale(1./eta_jets_HIP2_norm->Integral());
  eta_jets_HIP3_norm->Scale(1./eta_jets_HIP3_norm->Integral());

  eta_jets_HIP1_norm->Draw();
  eta_jets_HIP2_norm->Draw("SAME");
  eta_jets_HIP3_norm->Draw("SAME");
  eta_jets_mtOff1_norm->Draw("SAME");
  eta_jets_mtOff2_norm->Draw("SAME");
  eta_jets_mtOff3_norm->Draw("SAME");

  //Draw legend
  leg1->Draw();

  //Only with mitigation
  TCanvas* c_eta_norm_HIP = new TCanvas();
  eta_jets_HIP1_norm->Draw();
  eta_jets_HIP2_norm->Draw("SAME");
  eta_jets_HIP3_norm->Draw("SAME");

  //Draw legend
  leg3->Draw();

  //Only without mitigation
  TCanvas* c_eta_norm_mtOff = new TCanvas();
  eta_jets_mtOff1_norm->Draw();
  eta_jets_mtOff2_norm->Draw("SAME");
  eta_jets_mtOff3_norm->Draw("SAME");

  //Draw legend
  leg4->Draw();



  TCanvas* c_pt = new TCanvas();
  gPad->SetLogy();

  //Marker styles and colors
  pt_jets_mtOff1->SetMarkerStyle(24);
  pt_jets_mtOff2->SetMarkerStyle(24);
  pt_jets_mtOff3->SetMarkerStyle(24);
  pt_jets_HIP1->SetMarkerStyle(20);
  pt_jets_HIP2->SetMarkerStyle(20);
  pt_jets_HIP3->SetMarkerStyle(20);

  pt_jets_mtOff1->SetMarkerSize(0.75);
  pt_jets_mtOff2->SetMarkerSize(0.75);
  pt_jets_mtOff3->SetMarkerSize(0.75);
  pt_jets_HIP1->SetMarkerSize(0.75);
  pt_jets_HIP2->SetMarkerSize(0.75);
  pt_jets_HIP3->SetMarkerSize(0.75);

  pt_jets_mtOff1->SetMarkerColor(2);
  pt_jets_HIP1->SetMarkerColor(2);
  pt_jets_mtOff2->SetMarkerColor(3);
  pt_jets_HIP2->SetMarkerColor(3);
  pt_jets_mtOff3->SetMarkerColor(4);
  pt_jets_HIP3->SetMarkerColor(4);

  pt_jets_mtOff1->SetLineColor(2);
  pt_jets_HIP1->SetLineColor(2);
  pt_jets_mtOff2->SetLineColor(3);
  pt_jets_HIP2->SetLineColor(3);
  pt_jets_mtOff3->SetLineColor(4);
  pt_jets_HIP3->SetLineColor(4);

  //Set axis titles
  pt_jets_HIP1->GetXaxis()->SetTitle("p_{T,jet}");
  pt_jets_HIP1->GetYaxis()->SetTitle("N_{jets}");

  //Draw
  pt_jets_HIP1->Draw();
  pt_jets_HIP2->Draw("SAME");
  pt_jets_HIP3->Draw("SAME");
  pt_jets_mtOff1->Draw("SAME");
  pt_jets_mtOff2->Draw("SAME");
  pt_jets_mtOff3->Draw("SAME");

 //Legend
  TLegend *leg2;
  leg2 = new TLegend(0.63,0.63,0.9,0.85,"","brNDC");
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.045);
  leg2->SetFillColor(10);
  leg2->SetLineColor(1);
  leg2->SetTextFont(42);
  leg2->AddEntry(pt_jets_HIP1,"mtOn, L_{inst} = 0.6","lp");
  leg2->AddEntry(pt_jets_mtOff1,"mtOff, L_{inst} = 0.6","lp");
  leg2->AddEntry(pt_jets_HIP2,"mtOn, L_{inst} = 0.8","lp");
  leg2->AddEntry(pt_jets_mtOff2,"mtOff, L_{inst} = 0.8","lp");
  leg2->AddEntry(pt_jets_HIP3,"mtOn, L_{inst} = 1.2","lp");
  leg2->AddEntry(pt_jets_mtOff3,"mtOff, L_{inst} = 1.2","lp");
  leg2->Draw();


 //Only with mitigation
  TCanvas* c_pt_HIP = new TCanvas();
  gPad->SetLogy();
  pt_jets_HIP1->Draw();
  pt_jets_HIP2->Draw("SAME");
  pt_jets_HIP3->Draw("SAME");

  //Legend
  TLegend *leg5;
  leg5 = new TLegend(0.63,0.63,0.9,0.85,"","brNDC");
  leg5->SetBorderSize(0);
  leg5->SetTextSize(0.045);
  leg5->SetFillColor(10);
  leg5->SetLineColor(1);
  leg5->SetTextFont(42);
  leg5->AddEntry(pt_jets_HIP1,"mtOn, L_{inst} = 0.6","lp");
  leg5->AddEntry(pt_jets_HIP2,"mtOn, L_{inst} = 0.8","lp");
  leg5->AddEntry(pt_jets_HIP3,"mtOn, L_{inst} = 1.2","lp");
  leg5->Draw();

  //Only without mitigation
  TCanvas* c_pt_mtOff = new TCanvas();
  gPad->SetLogy();
  pt_jets_mtOff1->Draw();
  pt_jets_mtOff2->Draw("SAME");
  pt_jets_mtOff3->Draw("SAME");

  //Legend
  TLegend *leg6;
  leg6 = new TLegend(0.63,0.63,0.9,0.85,"","brNDC");
  leg6->SetBorderSize(0);
  leg6->SetTextSize(0.045);
  leg6->SetFillColor(10);
  leg6->SetLineColor(1);
  leg6->SetTextFont(42);
  leg6->AddEntry(pt_jets_mtOff1,"mtOff, L_{inst} = 0.6","lp");
  leg6->AddEntry(pt_jets_mtOff2,"mtOff, L_{inst} = 0.8","lp");
  leg6->AddEntry(pt_jets_mtOff3,"mtOff, L_{inst} = 1.2","lp");
  leg6->Draw();





  TCanvas* c_pt_norm = new TCanvas();
  gPad->SetLogy();

  //Clone histos
  TH1F* pt_jets_mtOff1_norm = (TH1F*)pt_jets_mtOff1->Clone();
  TH1F* pt_jets_mtOff2_norm = (TH1F*)pt_jets_mtOff2->Clone();
  TH1F* pt_jets_mtOff3_norm = (TH1F*)pt_jets_mtOff3->Clone();
  TH1F* pt_jets_HIP1_norm = (TH1F*)pt_jets_HIP1->Clone();
  TH1F* pt_jets_HIP2_norm = (TH1F*)pt_jets_HIP2->Clone();
  TH1F* pt_jets_HIP3_norm = (TH1F*)pt_jets_HIP3->Clone();

  pt_jets_HIP1_norm->GetYaxis()->SetTitle("N_{jets} normalized");

  //Normalise histos
  pt_jets_mtOff1_norm->Scale(1./pt_jets_mtOff1_norm->Integral());
  pt_jets_mtOff2_norm->Scale(1./pt_jets_mtOff2_norm->Integral());
  pt_jets_mtOff3_norm->Scale(1./pt_jets_mtOff3_norm->Integral());
  pt_jets_HIP1_norm->Scale(1./pt_jets_HIP1_norm->Integral());
  pt_jets_HIP2_norm->Scale(1./pt_jets_HIP2_norm->Integral());
  pt_jets_HIP3_norm->Scale(1./pt_jets_HIP3_norm->Integral());

  pt_jets_HIP1_norm->Draw();
  pt_jets_HIP2_norm->Draw("SAME");
  pt_jets_HIP3_norm->Draw("SAME");
  pt_jets_mtOff1_norm->Draw("SAME");
  pt_jets_mtOff2_norm->Draw("SAME");
  pt_jets_mtOff3_norm->Draw("SAME");

  //Draw legend
  leg2->Draw();

  //Only with mitigation
  TCanvas* c_pt_norm_HIP = new TCanvas();
  gPad->SetLogy();
  pt_jets_HIP1_norm->Draw();
  pt_jets_HIP2_norm->Draw("SAME");
  pt_jets_HIP3_norm->Draw("SAME");

  //Draw legend
  leg5->Draw();

  //Only without mitigation
  TCanvas* c_pt_norm_mtOff = new TCanvas();
  gPad->SetLogy();
  pt_jets_mtOff1_norm->Draw();
  pt_jets_mtOff2_norm->Draw("SAME");
  pt_jets_mtOff3_norm->Draw("SAME");

  //Draw legend
  leg6->Draw();

 


  TCanvas* c_mpf = new TCanvas("c_mpf","c_mpf",1000,600);
  c_mpf->Divide(4,2);


  //Draw
  c_mpf->cd(1);
  mpf_eta_HIP1->SetTitle("mtOn, L_{inst} = 0.6");
  mpf_eta_HIP1->Draw("COLZ");
  //gPad->SaveAs(path+"HIP_Mitigation_On_0p6_MPF_vs_Eta_Probe.pdf");
  c_mpf->cd(2);
  mpf_eta_HIP2->SetTitle("mtOn, L_{inst} = 0.8");
  mpf_eta_HIP2->Draw("COLZ");
  //gPad->SaveAs(path+"HIP_Mitigation_On_0p8_MPF_vs_Eta_Probe.pdf");
  c_mpf->cd(3);
  mpf_eta_HIP3->SetTitle("mtOn, L_{inst} = 1.2");
  mpf_eta_HIP3->Draw("COLZ");
  //gPad->SaveAs(path+"HIP_Mitigation_On_1p2_MPF_vs_Eta_Probe.pdf");
  c_mpf->cd(4);
  mpf_eta_HIP_data->SetTitle("DATA with mitigation");
  mpf_eta_HIP_data->Draw("COLZ");
  //gPad->SaveAs(path+"HIP_Mitigation_On_DATA_MPF_vs_Eta_Probe.pdf");
  c_mpf->cd(5);
  mpf_eta_mtOff1->SetTitle("mtOff, L_{inst} = 0.6");
  mpf_eta_mtOff1->Draw("COLZ");
  //gPad->SaveAs(path+"HIP_Mitigation_Off_0p6_MPF_vs_Eta_Probe.pdf");
  c_mpf->cd(6);
  mpf_eta_mtOff2->SetTitle("mtOff, L_{inst} = 0.8");
  mpf_eta_mtOff2->Draw("COLZ");
  //gPad->SaveAs(path+"HIP_Mitigation_Off_0p6_MPF_vs_Eta_Probe.pdf");
  c_mpf->cd(7);
  mpf_eta_mtOff3->SetTitle("mtOff, L_{inst} = 1.2");
  mpf_eta_mtOff3->Draw("COLZ");
  //gPad->SaveAs(path+"HIP_Mitigation_Off_1p2_MPF_vs_Eta_Probe.pdf");
  c_mpf->cd(8);
  mpf_eta_mtOff_data->SetTitle("DATA RunBCD");
  mpf_eta_mtOff_data->Draw("COLZ");
  //gPad->SaveAs(path+"HIP_Mitigation_Off_DATA_MPF_vs_Eta_Probe.pdf");


  TCanvas* c_rrel = new TCanvas("c_rrel","c_rrel",1000,600);
  c_rrel->Divide(4,2);


  //Draw
  c_rrel->cd(1);
  rrel_eta_HIP1->SetTitle("mtOn, L_{inst} = 0.6");
  rrel_eta_HIP1->Draw("COLZ");
  c_rrel->cd(2);
  rrel_eta_HIP2->SetTitle("mtOn, L_{inst} = 0.8");
  rrel_eta_HIP2->Draw("COLZ");
  c_rrel->cd(3);
  rrel_eta_HIP3->SetTitle("mtOn, L_{inst} = 1.2");
  rrel_eta_HIP3->Draw("COLZ");
  c_rrel->cd(4);
  rrel_eta_HIP_data->SetTitle("DATA with mitigation");
  rrel_eta_HIP_data->Draw("COLZ");
  c_rrel->cd(5);
  rrel_eta_mtOff1->SetTitle("mtOff, L_{inst} = 0.6");
  rrel_eta_mtOff1->Draw("COLZ");
  c_rrel->cd(6);
  rrel_eta_mtOff2->SetTitle("mtOff, L_{inst} = 0.8");
  rrel_eta_mtOff2->Draw("COLZ");
  c_rrel->cd(7);
  rrel_eta_mtOff3->SetTitle("mtOff, L_{inst} = 1.2");
  rrel_eta_mtOff3->Draw("COLZ");
  c_rrel->cd(8);
  rrel_eta_mtOff_data->SetTitle("DATA RunBCD");
  rrel_eta_mtOff_data->Draw("COLZ");




  TCanvas* c_mpf_nvert = new TCanvas("c_mpf_nvert","c_mpf_nvert",1000,600);
  c_mpf_nvert->Divide(4,2);


  //Draw
  c_mpf_nvert->cd(1);
  mpf_nvert_HIP1->SetTitle("mtOn, L_{inst} = 0.6");
  mpf_nvert_HIP1->Draw("COLZ");
  c_mpf_nvert->cd(2);
  mpf_nvert_HIP2->SetTitle("mtOn, L_{inst} = 0.8");
  mpf_nvert_HIP2->Draw("COLZ");
  c_mpf_nvert->cd(3);
  mpf_nvert_HIP3->SetTitle("mtOn, L_{inst} = 1.2");
  mpf_nvert_HIP3->Draw("COLZ");
  c_mpf_nvert->cd(4);
  mpf_nvert_HIP_data->SetTitle("DATA with mitigation");
  mpf_nvert_HIP_data->Draw("COLZ");
  c_mpf_nvert->cd(5);
  mpf_nvert_mtOff1->SetTitle("mtOff, L_{inst} = 0.6");
  mpf_nvert_mtOff1->Draw("COLZ");
  c_mpf_nvert->cd(6);
  mpf_nvert_mtOff2->SetTitle("mtOff, L_{inst} = 0.8");
  mpf_nvert_mtOff2->Draw("COLZ");
  c_mpf_nvert->cd(7);
  mpf_nvert_mtOff3->SetTitle("mtOff, L_{inst} = 1.2");
  mpf_nvert_mtOff3->Draw("COLZ");
  c_mpf_nvert->cd(8);
  mpf_nvert_mtOff_data->SetTitle("DATA RunBCD");
  mpf_nvert_mtOff_data->Draw("COLZ");


  TCanvas* c_rrel_nvert = new TCanvas("c_rrel_nvert","c_rrel_nvert",1000,600);
  c_rrel_nvert->Divide(4,2);


  //Draw
  c_rrel_nvert->cd(1);
  rrel_nvert_HIP1->SetTitle("mtOn, L_{inst} = 0.6");
  rrel_nvert_HIP1->Draw("COLZ");
  c_rrel_nvert->cd(2);
  rrel_nvert_HIP2->SetTitle("mtOn, L_{inst} = 0.8");
  rrel_nvert_HIP2->Draw("COLZ");
  c_rrel_nvert->cd(3);
  rrel_nvert_HIP3->SetTitle("mtOn, L_{inst} = 1.2");
  rrel_nvert_HIP3->Draw("COLZ");
  c_rrel_nvert->cd(4);
  rrel_nvert_HIP_data->SetTitle("DATA with mitigation");
  rrel_nvert_HIP_data->Draw("COLZ");
  c_rrel_nvert->cd(5);
  rrel_nvert_mtOff1->SetTitle("mtOff, L_{inst} = 0.6");
  rrel_nvert_mtOff1->Draw("COLZ");
  c_rrel_nvert->cd(6);
  rrel_nvert_mtOff2->SetTitle("mtOff, L_{inst} = 0.8");
  rrel_nvert_mtOff2->Draw("COLZ");
  c_rrel_nvert->cd(7);
  rrel_nvert_mtOff3->SetTitle("mtOff, L_{inst} = 1.2");
  rrel_nvert_mtOff3->Draw("COLZ");
  c_rrel_nvert->cd(8);
  rrel_nvert_mtOff_data->SetTitle("DATA RunBCD");
  rrel_nvert_mtOff_data->Draw("COLZ");



 




  TCanvas* c_Nj_eta = new TCanvas("c_Nj_eta","c_Nj_eta",1000,600);
  c_Nj_eta->Divide(4,2);

  //Draw
  c_Nj_eta->cd(1);
  njets_eta_HIP1->SetTitle("mtOn, L_{inst} = 0.6");
  njets_eta_HIP1->Draw("COLZ");
  c_Nj_eta->cd(2);
  njets_eta_HIP2->SetTitle("mtOn, L_{inst} = 0.8");
  njets_eta_HIP2->Draw("COLZ");
  c_Nj_eta->cd(3);
  njets_eta_HIP3->SetTitle("mtOn, L_{inst} = 1.2");
  njets_eta_HIP3->Draw("COLZ");
  c_Nj_eta->cd(4);
  njets_eta_HIP_data->SetTitle("DATA with mitigation");
  njets_eta_HIP_data->Draw("COLZ");
  c_Nj_eta->cd(5);
  njets_eta_mtOff1->SetTitle("mtOff, L_{inst} = 0.6");
  njets_eta_mtOff1->Draw("COLZ");
  c_Nj_eta->cd(6);
  njets_eta_mtOff2->SetTitle("mtOff, L_{inst} = 0.8");
  njets_eta_mtOff2->Draw("COLZ");
  c_Nj_eta->cd(7);
  njets_eta_mtOff3->SetTitle("mtOff, L_{inst} = 1.2");
  njets_eta_mtOff3->Draw("COLZ");
  c_Nj_eta->cd(8);
  njets_eta_mtOff_data->SetTitle("DATA RunBCD");
  njets_eta_mtOff_data->Draw("COLZ");


  TCanvas* c_Nj_pt = new TCanvas("c_Nj_pt","c_Nj_pt",1000,600);
  c_Nj_pt->Divide(4,2);

  //Draw
  c_Nj_pt->cd(1);
  njets_pt_HIP1->SetTitle("mtOn, L_{inst} = 0.6");
  njets_pt_HIP1->Draw("COLZ");
  c_Nj_pt->cd(2);
  njets_pt_HIP2->SetTitle("mtOn, L_{inst} = 0.8");
  njets_pt_HIP2->Draw("COLZ");
  c_Nj_pt->cd(3);
  njets_pt_HIP3->SetTitle("mtOn, L_{inst} = 1.2");
  njets_pt_HIP3->Draw("COLZ");
  c_Nj_pt->cd(4);
  njets_pt_HIP_data->SetTitle("DATA with mitigation");
  njets_pt_HIP_data->Draw("COLZ");
  c_Nj_pt->cd(5);
  njets_pt_mtOff1->SetTitle("mtOff, L_{inst} = 0.6");
  njets_pt_mtOff1->Draw("COLZ");
  c_Nj_pt->cd(6);
  njets_pt_mtOff2->SetTitle("mtOff, L_{inst} = 0.8");
  njets_pt_mtOff2->Draw("COLZ");
  c_Nj_pt->cd(7);
  njets_pt_mtOff3->SetTitle("mtOff, L_{inst} = 1.2");
  njets_pt_mtOff3->Draw("COLZ");
  c_Nj_pt->cd(8);
  njets_pt_mtOff_data->SetTitle("DATA RunBCD");
  njets_pt_mtOff_data->Draw("COLZ");


 //normalize each case (mtOn/Off) to the bin contents of L=0.6
  //mpf / rrel vs eta

  //clone
  TH2F* njets_eta_HIP1_norm = (TH2F*)njets_eta_HIP1->Clone("njets_eta_HIP1_norm");
  TH2F* njets_eta_HIP2_norm = (TH2F*)njets_eta_HIP2->Clone("njets_eta_HIP2_norm");
  TH2F* njets_eta_HIP3_norm = (TH2F*)njets_eta_HIP3->Clone("njets_eta_HIP3_norm");
  TH2F* njets_eta_mtOff1_norm = (TH2F*)njets_eta_mtOff1->Clone("njets_eta_mtOff1_norm");
  TH2F* njets_eta_mtOff2_norm = (TH2F*)njets_eta_mtOff2->Clone("njets_eta_mtOff2_norm");
  TH2F* njets_eta_mtOff3_norm = (TH2F*)njets_eta_mtOff3->Clone("njets_eta_mtOff3_norm");
  TH2F* njets_pt_HIP1_norm = (TH2F*)njets_pt_HIP1->Clone("njets_pt_HIP1_norm");
  TH2F* njets_pt_HIP2_norm = (TH2F*)njets_pt_HIP2->Clone("njets_pt_HIP2_norm");
  TH2F* njets_pt_HIP3_norm = (TH2F*)njets_pt_HIP3->Clone("njets_pt_HIP3_norm");
  TH2F* njets_pt_mtOff1_norm = (TH2F*)njets_pt_mtOff1->Clone("njets_pt_mtOff1_norm");
  TH2F* njets_pt_mtOff2_norm = (TH2F*)njets_pt_mtOff2->Clone("njets_pt_mtOff2_norm");
  TH2F* njets_pt_mtOff3_norm = (TH2F*)njets_pt_mtOff3->Clone("njets_pt_mtOff3_norm");
 
  njets_eta_HIP2_norm->Divide(njets_eta_HIP1);
  njets_eta_HIP3_norm->Divide(njets_eta_HIP1);
  njets_eta_HIP1_norm->Divide(njets_eta_HIP1);
  njets_eta_mtOff2_norm->Divide(njets_eta_mtOff1);
  njets_eta_mtOff3_norm->Divide(njets_eta_mtOff1);
  njets_eta_mtOff1_norm->Divide(njets_eta_mtOff1);
  njets_pt_HIP2_norm->Divide(njets_pt_HIP1);
  njets_pt_HIP3_norm->Divide(njets_pt_HIP1);
  njets_pt_HIP1_norm->Divide(njets_pt_HIP1);
  njets_pt_mtOff2_norm->Divide(njets_pt_mtOff1);
  njets_pt_mtOff3_norm->Divide(njets_pt_mtOff1);
  njets_pt_mtOff1_norm->Divide(njets_pt_mtOff1);


 TCanvas* c_njets_eta_norm = new TCanvas();
  c_njets_eta_norm->Divide(3,2);

  //Cosmetics
  njets_eta_HIP1_norm->GetZaxis()->SetRangeUser(0,2);
  njets_eta_HIP2_norm->GetZaxis()->SetRangeUser(0,2);
  njets_eta_HIP3_norm->GetZaxis()->SetRangeUser(0,2);
  njets_eta_mtOff1_norm->GetZaxis()->SetRangeUser(0,2);
  njets_eta_mtOff2_norm->GetZaxis()->SetRangeUser(0,2);
  njets_eta_mtOff3_norm->GetZaxis()->SetRangeUser(0,2);

  //Draw
  c_njets_eta_norm->cd(1);
  njets_eta_HIP1_norm->SetTitle("mtOn, L_{inst} = 0.6");
  njets_eta_HIP1_norm->Draw("COLZ");
  c_njets_eta_norm->cd(2);
  njets_eta_HIP2_norm->SetTitle("mtOn, L_{inst} = 0.8");
  njets_eta_HIP2_norm->Draw("COLZ");
  c_njets_eta_norm->cd(3);
  njets_eta_HIP3_norm->SetTitle("mtOn, L_{inst} = 1.2");
  njets_eta_HIP3_norm->Draw("COLZ");
  c_njets_eta_norm->cd(4);
  njets_eta_mtOff1_norm->SetTitle("mtOff, L_{inst} = 0.6");
  njets_eta_mtOff1_norm->Draw("COLZ");
  c_njets_eta_norm->cd(5);
  njets_eta_mtOff2_norm->SetTitle("mtOff, L_{inst} = 0.8");
  njets_eta_mtOff2_norm->Draw("COLZ");
  c_njets_eta_norm->cd(6);
  njets_eta_mtOff3_norm->SetTitle("mtOff, L_{inst} = 1.2");
  njets_eta_mtOff3_norm->Draw("COLZ");


  TCanvas* c_njets_pt_norm = new TCanvas();
  c_njets_pt_norm->Divide(3,2);

  //Cosmetics
  njets_pt_HIP1_norm->GetZaxis()->SetRangeUser(0,2);
  njets_pt_HIP2_norm->GetZaxis()->SetRangeUser(0,2);
  njets_pt_HIP3_norm->GetZaxis()->SetRangeUser(0,2);
  njets_pt_mtOff1_norm->GetZaxis()->SetRangeUser(0,2);
  njets_pt_mtOff2_norm->GetZaxis()->SetRangeUser(0,2);
  njets_pt_mtOff3_norm->GetZaxis()->SetRangeUser(0,2);

  //Draw
  c_njets_pt_norm->cd(1);
  njets_pt_HIP1_norm->SetTitle("mtOn, L_{inst} = 0.6");
  njets_pt_HIP1_norm->Draw("COLZ");
  c_njets_pt_norm->cd(2);
  njets_pt_HIP2_norm->SetTitle("mtOn, L_{inst} = 0.8");
  njets_pt_HIP2_norm->Draw("COLZ");
  c_njets_pt_norm->cd(3);
  njets_pt_HIP3_norm->SetTitle("mtOn, L_{inst} = 1.2");
  njets_pt_HIP3_norm->Draw("COLZ");
  c_njets_pt_norm->cd(4);
  njets_pt_mtOff1_norm->SetTitle("mtOff, L_{inst} = 0.6");
  njets_pt_mtOff1_norm->Draw("COLZ");
  c_njets_pt_norm->cd(5);
  njets_pt_mtOff2_norm->SetTitle("mtOff, L_{inst} = 0.8");
  njets_pt_mtOff2_norm->Draw("COLZ");
  c_njets_pt_norm->cd(6);
  njets_pt_mtOff3_norm->SetTitle("mtOff, L_{inst} = 1.2");
  njets_pt_mtOff3_norm->Draw("COLZ");



  //normalize each case (mtOn/Off) to the bin contents of L=0.6
  //mpf / rrel vs eta

  //clone
  TH2F* mpf_eta_HIP1_norm = (TH2F*)mpf_eta_HIP1->Clone("mpf_eta_HIP1_norm");
  TH2F* mpf_eta_HIP2_norm = (TH2F*)mpf_eta_HIP2->Clone("mpf_eta_HIP2_norm");
  TH2F* mpf_eta_HIP3_norm = (TH2F*)mpf_eta_HIP3->Clone("mpf_eta_HIP3_norm");
  TH2F* mpf_eta_mtOff1_norm = (TH2F*)mpf_eta_mtOff1->Clone("mpf_eta_mtOff1_norm");
  TH2F* mpf_eta_mtOff2_norm = (TH2F*)mpf_eta_mtOff2->Clone("mpf_eta_mtOff2_norm");
  TH2F* mpf_eta_mtOff3_norm = (TH2F*)mpf_eta_mtOff3->Clone("mpf_eta_mtOff3_norm");
  TH2F* rrel_eta_HIP1_norm = (TH2F*)rrel_eta_HIP1->Clone("rrel_eta_HIP1_norm");
  TH2F* rrel_eta_HIP2_norm = (TH2F*)rrel_eta_HIP2->Clone("rrel_eta_HIP2_norm");
  TH2F* rrel_eta_HIP3_norm = (TH2F*)rrel_eta_HIP3->Clone("rrel_eta_HIP3_norm");
  TH2F* rrel_eta_mtOff1_norm = (TH2F*)rrel_eta_mtOff1->Clone("rrel_eta_mtOff1_norm");
  TH2F* rrel_eta_mtOff2_norm = (TH2F*)rrel_eta_mtOff2->Clone("rrel_eta_mtOff2_norm");
  TH2F* rrel_eta_mtOff3_norm = (TH2F*)rrel_eta_mtOff3->Clone("rrel_eta_mtOff3_norm");
 
  mpf_eta_HIP2_norm->Divide(mpf_eta_HIP1);
  mpf_eta_HIP3_norm->Divide(mpf_eta_HIP1);
  mpf_eta_HIP1_norm->Divide(mpf_eta_HIP1);
  mpf_eta_mtOff2_norm->Divide(mpf_eta_mtOff1);
  mpf_eta_mtOff3_norm->Divide(mpf_eta_mtOff1);
  mpf_eta_mtOff1_norm->Divide(mpf_eta_mtOff1);
  rrel_eta_HIP2_norm->Divide(rrel_eta_HIP1);
  rrel_eta_HIP3_norm->Divide(rrel_eta_HIP1);
  rrel_eta_HIP1_norm->Divide(rrel_eta_HIP1);
  rrel_eta_mtOff2_norm->Divide(rrel_eta_mtOff1);
  rrel_eta_mtOff3_norm->Divide(rrel_eta_mtOff1);
  rrel_eta_mtOff1_norm->Divide(rrel_eta_mtOff1);


 TCanvas* c_mpf_norm = new TCanvas();
  c_mpf_norm->Divide(3,2);

  //Cosmetics
  mpf_eta_HIP1_norm->GetZaxis()->SetRangeUser(0,2);
  mpf_eta_HIP2_norm->GetZaxis()->SetRangeUser(0,2);
  mpf_eta_HIP3_norm->GetZaxis()->SetRangeUser(0,2);
  mpf_eta_mtOff1_norm->GetZaxis()->SetRangeUser(0,2);
  mpf_eta_mtOff2_norm->GetZaxis()->SetRangeUser(0,2);
  mpf_eta_mtOff3_norm->GetZaxis()->SetRangeUser(0,2);

  //Draw
  c_mpf_norm->cd(1);
  mpf_eta_HIP1_norm->SetTitle("mtOn, L_{inst} = 0.6");
  mpf_eta_HIP1_norm->Draw("COLZ");
  c_mpf_norm->cd(2);
  mpf_eta_HIP2_norm->SetTitle("mtOn, L_{inst} = 0.8");
  mpf_eta_HIP2_norm->Draw("COLZ");
  c_mpf_norm->cd(3);
  mpf_eta_HIP3_norm->SetTitle("mtOn, L_{inst} = 1.2");
  mpf_eta_HIP3_norm->Draw("COLZ");
  c_mpf_norm->cd(4);
  mpf_eta_mtOff1_norm->SetTitle("mtOff, L_{inst} = 0.6");
  mpf_eta_mtOff1_norm->Draw("COLZ");
  c_mpf_norm->cd(5);
  mpf_eta_mtOff2_norm->SetTitle("mtOff, L_{inst} = 0.8");
  mpf_eta_mtOff2_norm->Draw("COLZ");
  c_mpf_norm->cd(6);
  mpf_eta_mtOff3_norm->SetTitle("mtOff, L_{inst} = 1.2");
  mpf_eta_mtOff3_norm->Draw("COLZ");


  TCanvas* c_rrel_norm = new TCanvas();
  c_rrel_norm->Divide(3,2);

  //Cosmetics
  rrel_eta_HIP1_norm->GetZaxis()->SetRangeUser(0,2);
  rrel_eta_HIP2_norm->GetZaxis()->SetRangeUser(0,2);
  rrel_eta_HIP3_norm->GetZaxis()->SetRangeUser(0,2);
  rrel_eta_mtOff1_norm->GetZaxis()->SetRangeUser(0,2);
  rrel_eta_mtOff2_norm->GetZaxis()->SetRangeUser(0,2);
  rrel_eta_mtOff3_norm->GetZaxis()->SetRangeUser(0,2);

  //Draw
  c_rrel_norm->cd(1);
  rrel_eta_HIP1_norm->SetTitle("mtOn, L_{inst} = 0.6");
  rrel_eta_HIP1_norm->Draw("COLZ");
  c_rrel_norm->cd(2);
  rrel_eta_HIP2_norm->SetTitle("mtOn, L_{inst} = 0.8");
  rrel_eta_HIP2_norm->Draw("COLZ");
  c_rrel_norm->cd(3);
  rrel_eta_HIP3_norm->SetTitle("mtOn, L_{inst} = 1.2");
  rrel_eta_HIP3_norm->Draw("COLZ");
  c_rrel_norm->cd(4);
  rrel_eta_mtOff1_norm->SetTitle("mtOff, L_{inst} = 0.6");
  rrel_eta_mtOff1_norm->Draw("COLZ");
  c_rrel_norm->cd(5);
  rrel_eta_mtOff2_norm->SetTitle("mtOff, L_{inst} = 0.8");
  rrel_eta_mtOff2_norm->Draw("COLZ");
  c_rrel_norm->cd(6);
  rrel_eta_mtOff3_norm->SetTitle("mtOff, L_{inst} = 1.2");
  rrel_eta_mtOff3_norm->Draw("COLZ");



 //normalize each case (mtOn/Off) to the bin contents of L=0.6
  //mpf / rrel vs nvert

  //clone 
  TH2F* mpf_nvert_HIP1_norm = (TH2F*)mpf_nvert_HIP1->Clone("mpf_nvert_HIP1_norm");
  TH2F* mpf_nvert_HIP2_norm = (TH2F*)mpf_nvert_HIP2->Clone("mpf_nvert_HIP2_norm");
  TH2F* mpf_nvert_HIP3_norm = (TH2F*)mpf_nvert_HIP3->Clone("mpf_nvert_HIP3_norm");
  TH2F* mpf_nvert_mtOff1_norm = (TH2F*)mpf_nvert_mtOff1->Clone("mpf_nvert_mtOff1_norm");
  TH2F* mpf_nvert_mtOff2_norm = (TH2F*)mpf_nvert_mtOff2->Clone("mpf_nvert_mtOff2_norm");
  TH2F* mpf_nvert_mtOff3_norm = (TH2F*)mpf_nvert_mtOff3->Clone("mpf_nvert_mtOff3_norm");
  TH2F* rrel_nvert_HIP1_norm = (TH2F*)rrel_nvert_HIP1->Clone("rrel_nvert_HIP1_norm");
  TH2F* rrel_nvert_HIP2_norm = (TH2F*)rrel_nvert_HIP2->Clone("rrel_nvert_HIP2_norm");
  TH2F* rrel_nvert_HIP3_norm = (TH2F*)rrel_nvert_HIP3->Clone("rrel_nvert_HIP3_norm");
  TH2F* rrel_nvert_mtOff1_norm = (TH2F*)rrel_nvert_mtOff1->Clone("rrel_nvert_mtOff1_norm");
  TH2F* rrel_nvert_mtOff2_norm = (TH2F*)rrel_nvert_mtOff2->Clone("rrel_nvert_mtOff2_norm");
  TH2F* rrel_nvert_mtOff3_norm = (TH2F*)rrel_nvert_mtOff3->Clone("rrel_nvert_mtOff3_norm");
 
  mpf_nvert_HIP2_norm->Divide(mpf_nvert_HIP1);
  mpf_nvert_HIP3_norm->Divide(mpf_nvert_HIP1);
  mpf_nvert_HIP1_norm->Divide(mpf_nvert_HIP1);
  mpf_nvert_mtOff2_norm->Divide(mpf_nvert_mtOff1);
  mpf_nvert_mtOff3_norm->Divide(mpf_nvert_mtOff1);
  mpf_nvert_mtOff1_norm->Divide(mpf_nvert_mtOff1);
  rrel_nvert_HIP2_norm->Divide(rrel_nvert_HIP1);
  rrel_nvert_HIP3_norm->Divide(rrel_nvert_HIP1);
  rrel_nvert_HIP1_norm->Divide(rrel_nvert_HIP1);
  rrel_nvert_mtOff2_norm->Divide(rrel_nvert_mtOff1);
  rrel_nvert_mtOff3_norm->Divide(rrel_nvert_mtOff1);
  rrel_nvert_mtOff1_norm->Divide(rrel_nvert_mtOff1);


 TCanvas* c_mpf_nvert_norm = new TCanvas();
  c_mpf_nvert_norm->Divide(3,2);

  //Cosmetics
  mpf_nvert_HIP1_norm->GetZaxis()->SetRangeUser(0,2);
  mpf_nvert_HIP2_norm->GetZaxis()->SetRangeUser(0,2);
  mpf_nvert_HIP3_norm->GetZaxis()->SetRangeUser(0,2);
  mpf_nvert_mtOff1_norm->GetZaxis()->SetRangeUser(0,2);
  mpf_nvert_mtOff2_norm->GetZaxis()->SetRangeUser(0,2);
  mpf_nvert_mtOff3_norm->GetZaxis()->SetRangeUser(0,2);

  //Draw
  c_mpf_nvert_norm->cd(1);
  mpf_nvert_HIP1_norm->SetTitle("mtOn, L_{inst} = 0.6");
  mpf_nvert_HIP1_norm->Draw("COLZ");
  c_mpf_nvert_norm->cd(2);
  mpf_nvert_HIP2_norm->SetTitle("mtOn, L_{inst} = 0.8");
  mpf_nvert_HIP2_norm->Draw("COLZ");
  c_mpf_nvert_norm->cd(3);
  mpf_nvert_HIP3_norm->SetTitle("mtOn, L_{inst} = 1.2");
  mpf_nvert_HIP3_norm->Draw("COLZ");
  c_mpf_nvert_norm->cd(4);
  mpf_nvert_mtOff1_norm->SetTitle("mtOff, L_{inst} = 0.6");
  mpf_nvert_mtOff1_norm->Draw("COLZ");
  c_mpf_nvert_norm->cd(5);
  mpf_nvert_mtOff2_norm->SetTitle("mtOff, L_{inst} = 0.8");
  mpf_nvert_mtOff2_norm->Draw("COLZ");
  c_mpf_nvert_norm->cd(6);
  mpf_nvert_mtOff3_norm->SetTitle("mtOff, L_{inst} = 1.2");
  mpf_nvert_mtOff3_norm->Draw("COLZ");


  TCanvas* c_rrel_nvert_norm = new TCanvas();
  c_rrel_nvert_norm->Divide(3,2);

  //Cosmetics
  rrel_nvert_HIP1_norm->GetZaxis()->SetRangeUser(0,2);
  rrel_nvert_HIP2_norm->GetZaxis()->SetRangeUser(0,2);
  rrel_nvert_HIP3_norm->GetZaxis()->SetRangeUser(0,2);
  rrel_nvert_mtOff1_norm->GetZaxis()->SetRangeUser(0,2);
  rrel_nvert_mtOff2_norm->GetZaxis()->SetRangeUser(0,2);
  rrel_nvert_mtOff3_norm->GetZaxis()->SetRangeUser(0,2);

  //Draw
  c_rrel_nvert_norm->cd(1);
  rrel_nvert_HIP1_norm->SetTitle("mtOn, L_{inst} = 0.6");
  rrel_nvert_HIP1_norm->Draw("COLZ");
  c_rrel_nvert_norm->cd(2);
  rrel_nvert_HIP2_norm->SetTitle("mtOn, L_{inst} = 0.8");
  rrel_nvert_HIP2_norm->Draw("COLZ");
  c_rrel_nvert_norm->cd(3);
  rrel_nvert_HIP3_norm->SetTitle("mtOn, L_{inst} = 1.2");
  rrel_nvert_HIP3_norm->Draw("COLZ");
  c_rrel_nvert_norm->cd(4);
  rrel_nvert_mtOff1_norm->SetTitle("mtOff, L_{inst} = 0.6");
  rrel_nvert_mtOff1_norm->Draw("COLZ");
  c_rrel_nvert_norm->cd(5);
  rrel_nvert_mtOff2_norm->SetTitle("mtOff, L_{inst} = 0.8");
  rrel_nvert_mtOff2_norm->Draw("COLZ");
  c_rrel_nvert_norm->cd(6);
  rrel_nvert_mtOff3_norm->SetTitle("mtOff, L_{inst} = 1.2");
  rrel_nvert_mtOff3_norm->Draw("COLZ");
  


  /*
  c_eta->SaveAs(path+"plots/HIP_Mitigation_Comparison_Eta_Jets.eps");
  c_eta->SaveAs(path+"plots/HIP_Mitigation_Comparison_Eta_Jets.pdf");
  c_eta_HIP->SaveAs(path+"plots/HIP_Mitigation_On_Eta_Jets.eps");
  c_eta_HIP->SaveAs(path+"plots/HIP_Mitigation_On_Eta_Jets.pdf");
  c_eta_mtOff->SaveAs(path+"plots/HIP_Mitigation_Off_Eta_Jets.eps");
  c_eta_mtOff->SaveAs(path+"plots/HIP_Mitigation_Off_Eta_Jets.pdf");

  c_eta_norm->SaveAs(path+"plots/HIP_Mitigation_Comparison_Eta_Jets_norm.eps");
  c_eta_norm->SaveAs(path+"plots/HIP_Mitigation_Comparison_Eta_Jets_norm.pdf");
  c_eta_norm_HIP->SaveAs(path+"plots/HIP_Mitigation_On_Eta_Jets_norm.eps");
  c_eta_norm_HIP->SaveAs(path+"plots/HIP_Mitigation_On_Eta_Jets_norm.pdf");
  c_eta_norm_mtOff->SaveAs(path+"plots/HIP_Mitigation_Off_Eta_Jets_norm.eps");
  c_eta_norm_mtOff->SaveAs(path+"plots/HIP_Mitigation_Off_Eta_Jets_norm.pdf");

  c_pt->SaveAs(path+"plots/HIP_Mitigation_Comparison_Pt_Jets.eps");
  c_pt->SaveAs(path+"plots/HIP_Mitigation_Comparison_Pt_Jets.pdf");
  c_pt_HIP->SaveAs(path+"plots/HIP_Mitigation_On_Pt_Jets.eps");
  c_pt_HIP->SaveAs(path+"plots/HIP_Mitigation_On_Pt_Jets.pdf");
  c_pt_mtOff->SaveAs(path+"plots/HIP_Mitigation_Off_Pt_Jets.eps");
  c_pt_mtOff->SaveAs(path+"plots/HIP_Mitigation_Off_Pt_Jets.pdf");

  c_pt_norm->SaveAs(path+"plots/HIP_Mitigation_Comparison_Pt_Jets_norm.eps");
  c_pt_norm->SaveAs(path+"plots/HIP_Mitigation_Comparison_Pt_Jets_norm.pdf");
  c_eta_norm_HIP->SaveAs(path+"plots/HIP_Mitigation_On_Pt_Jets_norm.eps");
  c_pt_norm_HIP->SaveAs(path+"plots/HIP_Mitigation_On_Pt_Jets_norm.pdf");
  c_pt_norm_mtOff->SaveAs(path+"plots/HIP_Mitigation_Off_Pt_Jets_norm.eps");
  c_pt_norm_mtOff->SaveAs(path+"plots/HIP_Mitigation_Off_Pt_Jets_norm.pdf");
*/

  c_mpf->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_MPF_vs_Eta_Probe.eps");
  c_mpf->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_MPF_vs_Eta_Probe.pdf");
  c_rrel->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_R_rel_vs_Eta_Probe.eps");
  c_rrel->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_R_rel_vs_Eta_Probe.pdf");
  c_mpf_nvert->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_MPF_vs_Nvert_Probe.eps");
  c_mpf_nvert->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_MPF_vs_Nvert_Probe.pdf");
  c_rrel_nvert->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_R_rel_vs_Nvert_Probe.eps");
  c_rrel_nvert->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_R_rel_vs_Nvert_Probe.pdf");
  c_Nj_pt->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_NJets_vs_Pt_Probe.eps");
  c_Nj_pt->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_NJets_vs_Pt_Probe.pdf");
  c_Nj_eta->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_NJets_vs_Eta_Probe.eps");
  c_Nj_eta->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_NJets_vs_Eta_Probe.pdf");

  c_mpf_norm->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_MPF_vs_Eta_Probe_norm.eps");
  c_mpf_norm->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_MPF_vs_Eta_Probe_norm.pdf");
  c_rrel_norm->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_R_rel_vs_Eta_Probe_norm.eps");
  c_rrel_norm->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_R_rel_vs_Eta_Probe_norm.pdf");
  c_mpf_nvert_norm->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_MPF_vs_Nvert_Probe_norm.eps");
  c_mpf_nvert_norm->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_MPF_vs_Nvert_Probe_norm.pdf");
  c_rrel_nvert_norm->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_R_rel_vs_Nvert_Probe_norm.eps");
  c_rrel_nvert_norm->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_R_rel_vs_Nvert_Probe_norm.pdf");
  c_njets_pt_norm->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_NJets_vs_Pt_Probe_norm.eps");
  c_njets_pt_norm->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_NJets_vs_Pt_Probe_norm.pdf");
  c_njets_eta_norm->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_NJets_vs_Eta_Probe_norm.eps");
  c_njets_eta_norm->SaveAs(path+"plots/HIP_Mitigation_" + jettag +"_NJets_vs_Eta_Probe_norm.pdf");

  //save single plots as well
  TCanvas* c_tmp = new TCanvas();
  mpf_eta_HIP1->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p6_MPF_vs_Eta_Probe.pdf");
  mpf_eta_HIP2->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p8_MPF_vs_Eta_Probe.pdf");
  mpf_eta_HIP3->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_1p2_MPF_vs_Eta_Probe.pdf");
  mpf_eta_mtOff1->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p6_MPF_vs_Eta_Probe.pdf");
  mpf_eta_mtOff2->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p8_MPF_vs_Eta_Probe.pdf");
  mpf_eta_mtOff3->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_1p2_MPF_vs_Eta_Probe.pdf");


  rrel_eta_HIP1->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p6_R_rel_vs_Eta_Probe.pdf");
  rrel_eta_HIP2->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p8_R_rel_vs_Eta_Probe.pdf");
  rrel_eta_HIP3->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_1p2_R_rel_vs_Eta_Probe.pdf");
  rrel_eta_mtOff1->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p6_R_rel_vs_Eta_Probe.pdf");
  rrel_eta_mtOff2->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p8_R_rel_vs_Eta_Probe.pdf");
  rrel_eta_mtOff3->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_1p2_R_rel_vs_Eta_Probe.pdf");

  mpf_nvert_HIP1->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p6_MPF_vs_Nvert_Probe.pdf");
  mpf_nvert_HIP2->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p8_MPF_vs_Nvert_Probe.pdf");
  mpf_nvert_HIP3->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_1p2_MPF_vs_Nvert_Probe.pdf");
  mpf_nvert_mtOff1->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p6_MPF_vs_Nvert_Probe.pdf");
  mpf_nvert_mtOff2->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p8_MPF_vs_Nvert_Probe.pdf");
  mpf_nvert_mtOff3->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_1p2_MPF_vs_Nvert_Probe.pdf");


  rrel_nvert_HIP1->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p6_R_rel_vs_Nvert_Probe.pdf");
  rrel_nvert_HIP2->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p8_R_rel_vs_Nvert_Probe.pdf");
  rrel_nvert_HIP3->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_1p2_R_rel_vs_Nvert_Probe.pdf");
  rrel_nvert_mtOff1->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p6_R_rel_vs_Nvert_Probe.pdf");
  rrel_nvert_mtOff2->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p8_R_rel_vs_Nvert_Probe.pdf");
  rrel_nvert_mtOff3->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_1p2_R_rel_vs_Nvert_Probe.pdf");

  njets_eta_HIP1->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p6_NJets_vs_Eta_Probe.pdf");
  njets_eta_HIP2->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p8_NJets_vs_Eta_Probe.pdf");
  njets_eta_HIP3->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_1p2_NJets_vs_Eta_Probe.pdf");
  njets_eta_mtOff1->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p6_NJets_vs_Eta_Probe.pdf");
  njets_eta_mtOff2->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p8_NJets_vs_Eta_Probe.pdf");
  njets_eta_mtOff3->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_1p2_NJets_vs_Eta_Probe.pdf");

  njets_pt_HIP1->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p6_NJets_vs_Pt_Probe.pdf");
  njets_pt_HIP2->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p8_NJets_vs_Pt_Probe.pdf");
  njets_pt_HIP3->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_1p2_NJets_vs_Pt_Probe.pdf");
  njets_pt_mtOff1->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p6_NJets_vs_Pt_Probe.pdf");
  njets_pt_mtOff2->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p8_NJets_vs_Pt_Probe.pdf");
  njets_pt_mtOff3->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_1p2_NJets_vs_Pt_Probe.pdf");

  //normalized
  mpf_eta_HIP1_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p6_MPF_vs_Eta_Probe_norm.pdf");
  mpf_eta_HIP2_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p8_MPF_vs_Eta_Probe_norm.pdf");
  mpf_eta_HIP3_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_1p2_MPF_vs_Eta_Probe_norm.pdf");
  mpf_eta_mtOff1_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p6_MPF_vs_Eta_Probe_norm.pdf");
  mpf_eta_mtOff2_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p8_MPF_vs_Eta_Probe_norm.pdf");
  mpf_eta_mtOff3_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_1p2_MPF_vs_Eta_Probe_norm.pdf");


  rrel_eta_HIP1_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p6_R_rel_vs_Eta_Probe_norm.pdf");
  rrel_eta_HIP2_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p8_R_rel_vs_Eta_Probe_norm.pdf");
  rrel_eta_HIP3_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_1p2_R_rel_vs_Eta_Probe_norm.pdf");
  rrel_eta_mtOff1_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p6_R_rel_vs_Eta_Probe_norm.pdf");
  rrel_eta_mtOff2_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p8_R_rel_vs_Eta_Probe_norm.pdf");
  rrel_eta_mtOff3_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_1p2_R_rel_vs_Eta_Probe_norm.pdf");

  mpf_nvert_HIP1_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p6_MPF_vs_Nvert_Probe_norm.pdf");
  mpf_nvert_HIP2_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p8_MPF_vs_Nvert_Probe_norm.pdf");
  mpf_nvert_HIP3_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_1p2_MPF_vs_Nvert_Probe_norm.pdf");
  mpf_nvert_mtOff1_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p6_MPF_vs_Nvert_Probe_norm.pdf");
  mpf_nvert_mtOff2_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p8_MPF_vs_Nvert_Probe_norm.pdf");
  mpf_nvert_mtOff3_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_1p2_MPF_vs_Nvert_Probe_norm.pdf");


  rrel_nvert_HIP1_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p6_R_rel_vs_Nvert_Probe_norm.pdf");
  rrel_nvert_HIP2_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p8_R_rel_vs_Nvert_Probe_norm.pdf");
  rrel_nvert_HIP3_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_1p2_R_rel_vs_Nvert_Probe_norm.pdf");
  rrel_nvert_mtOff1_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p6_R_rel_vs_Nvert_Probe_norm.pdf");
  rrel_nvert_mtOff2_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p8_R_rel_vs_Nvert_Probe_norm.pdf");
  rrel_nvert_mtOff3_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_1p2_R_rel_vs_Nvert_Probe_norm.pdf");

  njets_eta_HIP1_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p6_NJets_vs_Eta_Probe_norm.pdf");
  njets_eta_HIP2_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p8_NJets_vs_Eta_Probe_norm.pdf");
  njets_eta_HIP3_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_1p2_NJets_vs_Eta_Probe_norm.pdf");
  njets_eta_mtOff1_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p6_NJets_vs_Eta_Probe_norm.pdf");
  njets_eta_mtOff2_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p8_NJets_vs_Eta_Probe_norm.pdf");
  njets_eta_mtOff3_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_1p2_NJets_vs_Eta_Probe_norm.pdf");

  njets_pt_HIP1_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p6_NJets_vs_Pt_Probe_norm.pdf");
  njets_pt_HIP2_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_0p8_NJets_vs_Pt_Probe_norm.pdf");
  njets_pt_HIP3_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_1p2_NJets_vs_Pt_Probe_norm.pdf");
  njets_pt_mtOff1_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p6_NJets_vs_Pt_Probe_norm.pdf");
  njets_pt_mtOff2_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_0p8_NJets_vs_Pt_Probe_norm.pdf");
  njets_pt_mtOff3_norm->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_1p2_NJets_vs_Pt_Probe_norm.pdf");

  mpf_eta_HIP_data->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_DATA_MPF_vs_Eta_Probe.pdf");
  mpf_eta_mtOff_data->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_DATA_MPF_vs_Eta_Probe.pdf");
  rrel_eta_HIP_data->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_DATA_R_rel_vs_Eta_Probe.pdf");
  rrel_eta_mtOff_data->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_DATA_R_rel_vs_Eta_Probe.pdf");
  njets_eta_HIP_data->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_DATA_NJets_vs_Eta_Probe.pdf");
  njets_eta_mtOff_data->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_DATA_NJets_vs_Eta_Probe.pdf");
  njets_pt_HIP_data->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_On_" + jettag +"_DATA_NJets_vs_Pt_Probe.pdf");
  njets_pt_mtOff_data->Draw("COLZ");
  c_tmp->SaveAs(path+"plots/HIP_Mitigation_Off_" + jettag +"_DATA_NJets_vs_Pt_Probe.pdf");

}
