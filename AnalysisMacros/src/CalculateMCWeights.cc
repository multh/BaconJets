#include "../include/CorrectionObject.h"
#include "../include/parameters.h"
#include "../include/useful_functions.h"

#include <TStyle.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TFile.h>
#include <TLine.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TProfile.h>



using namespace std;

void CorrectionObject::CalculateMCWeights(){
  cout << "--------------- Starting CalculateMCWeights() ---------------" << endl << endl;
  gStyle->SetOptStat(0);

  const int n_pt_bins = 1000;

  //Files
  TFile* f_mc   = new TFile(CorrectionObject::_basepath + CorrectionObject::_collection + "/DefaultWeights/uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_AK4CHS_Full.root","READ");
  TFile* f_data = new TFile(CorrectionObject::_basepath + CorrectionObject::_collection + "/DefaultWeights/uhh2.AnalysisModuleRunner.DATA.DATA_Run" + CorrectionObject::_runnr  + "_AK4CHS.root","READ");

 //pT_ave-histograms for MC & DATA

  TH2D* h_pt_ave_mc = new TH2D("pt_ave_mc","pt_ave mc;p_{T}^{ave};|#eta|", n_pt_bins,0,5000,n_eta-1, eta_bins);
  TH2D* h_pt_ave_data = new TH2D("pt_ave_data","pt_ave data;p_{T}^{ave};|#eta|", n_pt_bins,0,5000,n_eta-1,eta_bins);
  TH2D* h_pt_ave_mc_scaled = new TH2D("pt_ave_mc_scaled","pt_ave mc scaled;p_{T}^{ave};|#eta|", n_pt_bins,0,5000,n_eta-1,eta_bins);
  TH1D* h1_pt_ave_mc_scaled = new TH1D("h1_pt_ave_mc_scaled", "pt_ave mc scaled;p_{T}^{ave};entries", n_pt_bins, 0, 5000); //cross-check 1d
  TH1D* h1_pt_ave_data = new TH1D("h1_pt_ave_data", "pt_ave data;p_{T}^{ave};entries", n_pt_bins, 0, 5000);                //cross-check 1d

  //normalization histograms for MC & DATA
  double bins[16] = {51, 73, 95, 100, 126, 152, 163, 230, 250, 299, 319, 365, 433, 453, 566, 10000};
  double bins2[16] = {51, 73, 95, 100, 126, 152, 163, 230, 250, 299, 319, 365, 433, 453, 566, 1000};
  TH2D* h_pt_ave_binned_mc = new TH2D("pt_ave_binned_mc","pt_ave binned mc;p_{T}^{ave};|#eta|", 15, bins,n_eta-1,eta_bins);
  TH2D* h_pt_ave_binned_data = new TH2D("pt_ave_binned_data","pt_ave binned data;p_{T}^{ave};|#eta|", 15, bins,n_eta-1,eta_bins);
  TH2D* h_pt_ave_binned_mc_scaled = new TH2D("pt_ave_binned_mc_scaled","pt_ave binned mc scaled;p_{T}^{ave};|#eta|", 15, bins,n_eta-1,eta_bins);
  TH1D* h1_pt_ave_binned_mc_scaled = new TH1D("h1_pt_ave_binned_mc_scaled", "pt_ave mc scaled;p_{T}^{ave};entries", 15, bins2);   //cross-check 1d
  TH1D* h1_pt_ave_binned_data = new TH1D("h1_pt_ave_binned_data", "pt_ave data;p_{T}^{ave};entries", 15, bins2);                  //cross-check 1d 


  //Fill histograms
  TTreeReader myReader_data("AnalysisTree", f_data);
  TTreeReaderValue<Float_t> pt_ave_data(myReader_data, "pt_ave");
  TTreeReaderValue<Float_t> weight_data(myReader_data, "weight");
  TTreeReaderValue<Float_t> alpha_data(myReader_data, "alpha");
  TTreeReaderValue<Float_t> probejet_eta_data(myReader_data, "probejet_eta");
  while (myReader_data.Next()){
    if(*alpha_data >= 0.3) continue;
    h_pt_ave_data->Fill(*pt_ave_data, fabs(*probejet_eta_data), *weight_data);
    h_pt_ave_binned_data->Fill(*pt_ave_data, fabs(*probejet_eta_data), *weight_data);
    h1_pt_ave_data->Fill(*pt_ave_data, *weight_data);
    h1_pt_ave_binned_data->Fill(*pt_ave_data, *weight_data);
  }

  TTreeReader myReader_mc("AnalysisTree", f_mc);
  TTreeReaderValue<Float_t> pt_ave_mc(myReader_mc, "pt_ave");
  TTreeReaderValue<Float_t> weight_mc(myReader_mc, "weight");
  TTreeReaderValue<Float_t> alpha_mc(myReader_mc, "alpha");
  TTreeReaderValue<Float_t> probejet_eta_mc(myReader_mc, "probejet_eta");
  while (myReader_mc.Next()){
    if(*alpha_mc >= 0.3) continue;
    h_pt_ave_mc->Fill(*pt_ave_mc, fabs(*probejet_eta_mc), *weight_mc);
    h_pt_ave_binned_mc->Fill(*pt_ave_mc, fabs(*probejet_eta_mc), *weight_mc);
  }

  //Calculate scale factors
  TH2D* SF =  (TH2D*)h_pt_ave_data->Clone();
  TH2D* for_SF_mc =  (TH2D*)h_pt_ave_mc->Clone();
  
  for(int i=0; i<SF->GetNbinsX(); i++){
    for(int j=0; j<n_eta-1; j++){
      double content = 0;
      double data = SF->GetBinContent(i+1,j+1);
      double mc = for_SF_mc->GetBinContent(i+1,j+1);
      if(mc > 0) content = data/mc;
      else content = 0;
      SF->SetBinContent(i+1, j+1, content);
      //cout << "lower pt: " << 5+i*5 << ", eta-bin: " << j+1 << " DATA: " << h_pt_ave_data->GetBinContent(i+1,j+1) << ", MC: " << for_SF_mc->GetBinContent(i+1,j+1) << ", content: " << content << endl;
    }
  }


  TTreeReader myReader_mc2("AnalysisTree", f_mc);
  TTreeReaderValue<Float_t> pt_ave_mc2(myReader_mc2, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_mc2(myReader_mc2, "probejet_eta");
  TTreeReaderValue<Float_t> weight_mc2(myReader_mc2, "weight");
  TTreeReaderValue<Float_t> alpha_mc2(myReader_mc2, "alpha");

  int ind = 0;
  while (myReader_mc2.Next()){
    if(*alpha_mc2 >= 0.3) continue;
    double right_SF = 0;
    double pt = *pt_ave_mc2;
    double eta = fabs(*probejet_eta_mc2);
 

    int idx_x = 0;
    int idx_y = 0;
    while(pt > idx_x*5){ 
      idx_x++;
    }
    while(eta > eta_bins[idx_y]){
      idx_y++;
    }
    
    //Now i is the number of the bin, whose weight has to be taken for this event
    right_SF = SF->GetBinContent(idx_x,idx_y);

    double right_weight = *weight_mc2 * right_SF;

    h_pt_ave_mc_scaled->Fill(pt, eta, right_weight);
    h_pt_ave_binned_mc_scaled->Fill(pt, eta, right_weight);
    h1_pt_ave_mc_scaled->Fill(pt, right_weight);
    h1_pt_ave_binned_mc_scaled->Fill(pt, right_weight);

    ind++;
  } //looping over events for the 2nd time

  for(int i=0; i<n_pt_bins; i++){
    for(int j=0; j<n_eta-1; j++){
      //cout << "Data: " << h_pt_ave_data->GetBinContent(i+1, j+1) << ", MC: " << h_pt_ave_mc_scaled->GetBinContent(i+1, j+1) << endl;
      //cout << "Data: " << h_pt_ave_data->GetBinContent(i+1, j+1) << ", SF: " << SF->GetBinContent(i+1, j+1) << endl << endl;
      if(fabs(h_pt_ave_mc_scaled->GetBinContent(i+1, j+1) - h_pt_ave_data->GetBinContent(i+1, j+1)) > 0.05) cout << " MC: " << h_pt_ave_mc_scaled->GetBinContent(i+1, j+1) << ", DATA: " << h_pt_ave_data->GetBinContent(i+1, j+1) << endl;//throw runtime_error("In CalculateMCWeights.cc: Scaled MC bin-content does not match the bin-content in data.");
    }
  }

  TCanvas* c1 = new TCanvas();
  h1_pt_ave_mc_scaled->SetLineWidth(2);
  h1_pt_ave_mc_scaled->GetXaxis()->SetTitle("p_{T}^{ave} [GeV]");
  h1_pt_ave_mc_scaled->GetXaxis()->SetRangeUser(20,1000);
  h1_pt_ave_mc_scaled->Draw("HIST");
  h1_pt_ave_data->SetLineColor(2);
  h1_pt_ave_data->SetMarkerColor(2);
  h1_pt_ave_data->SetMarkerSize(0.5);
  h1_pt_ave_data->SetMarkerStyle(20);
  h1_pt_ave_data->Draw("P SAME");

  TLegend* l1 = new TLegend(0.52,0.7,0.9,0.9);    
  l1->AddEntry(h1_pt_ave_mc_scaled,"MC after reweighting","l");
  l1->AddEntry(h1_pt_ave_data,"DATA","p");
  l1->Draw();

  c1->SaveAs(CorrectionObject::_basepath + CorrectionObject::_collection + "/MC_scaled_PtEta_Fine.pdf");


  TCanvas* c2 = new TCanvas();
  h1_pt_ave_binned_mc_scaled->SetLineWidth(2);
  h1_pt_ave_binned_mc_scaled->GetXaxis()->SetTitle("p_{T}^{ave} [GeV]");
  h1_pt_ave_binned_mc_scaled->GetXaxis()->SetRangeUser(20,1000);
  h1_pt_ave_binned_mc_scaled->Draw("HIST");
  h1_pt_ave_binned_data->SetLineColor(2);
  h1_pt_ave_binned_data->SetMarkerColor(2);
  h1_pt_ave_binned_data->SetMarkerStyle(20);
  h1_pt_ave_binned_data->Draw("P SAME");

  TLegend* l2 = new TLegend(0.52,0.2,0.90,0.4);    
  l2->AddEntry(h1_pt_ave_binned_mc_scaled,"MC after reweighting","l");
  l2->AddEntry(h1_pt_ave_binned_data,"DATA","p");
  l2->Draw();

  c2->SaveAs(CorrectionObject::_basepath + CorrectionObject::_collection + "/MC_scaled_PtEta_Fine_binned.pdf");


  //TFile* out = new TFile(CorrectionObject::_basepath + CorrectionObject::_collection + "/MC_ReWeights_Run" + CorrectionObject::_runnr  + ".root","RECREATE");
  //SF->Write();
  //out->Close();

  //delete out;
  delete l2;
  delete c2;
  delete l1;
  delete c1;
  delete for_SF_mc;
  delete SF;
  delete h_pt_ave_binned_mc_scaled;
  delete h_pt_ave_binned_data;
  delete h_pt_ave_binned_mc;
  delete h_pt_ave_mc_scaled;
  delete h_pt_ave_data;
  delete h_pt_ave_mc;
  delete f_data;
  delete f_mc;

}
