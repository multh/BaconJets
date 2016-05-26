#include "header.h"

void alpha(TString path){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // get the residual root file for MPF and pt balance
  TFile* mpf = new TFile(path+"Histo_KFSR_MPF_L1.root","READ");
  TH1D* ampf = mpf->Get("kfsr_mpf");
  TFile* dijet = new TFile(path+"Histo_KFSR_DiJet_L1.root","READ"); 
  TH1D* adijet = dijet->Get("kfsr_dijet");



  // create histo for the normalization
  TH1D* ampf_norm=(TH1D*)ampf->Clone();
  TH1D* adijet_norm=(TH1D*)adijet->Clone();

  // create norm factor
  double norm_factor_MPF;
  double norm_factor_DiJet;

  for (int i=1; i<n_etabarr+1; i++){
    norm_factor_MPF+=ampf->GetBinContent(i);
    norm_factor_DiJet+=adijet->GetBinContent(i);
  }

  norm_factor_MPF = norm_factor_MPF/n_etabarr;
  norm_factor_DiJet = norm_factor_DiJet/n_etabarr;

  cout << "normalization factor mpf: " << norm_factor_MPF << endl;
  cout << "normalization factor pt:  " << norm_factor_DiJet << endl;

  // scale with norm factor
  ampf_norm->Scale(1/norm_factor_MPF);
  adijet_norm->Scale(1/norm_factor_DiJet);

  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,5.191,1);

  gROOT->LoadMacro("tdrstyle_mod14.C");

  ampf_norm->SetLineWidth(2);
  ampf_norm->SetLineColor(kRed+1);
  ampf_norm->Draw("E1");
  adijet_norm->SetLineWidth(2);
  adijet_norm->SetLineColor(kBlue+1);
  adijet_norm->Draw("E1 SAME");

  // Draw results
  TH1D *h = new TH1D("h",";|#eta|;k_{FSR}",41,0,5.191);
  h->SetMaximum(1.06);
  h->SetMinimum(0.98);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  lumi_13TeV = "2.11 fb^{-1}";
  TCanvas *c1 = tdrCanvas("c1",h,4,0,kSquare);
  //  TCanvas * c = new TCanvas("c", "c", 1000,1000);
  //gPad->SetGridy();

  ampf_norm->Draw("E1 SAME");
  adijet_norm->Draw("E1 SAME");
  line->SetLineStyle(2);
  line->Draw("SAME");

  TLegend *leg1 = tdrLeg(0.20,0.59,0.43,0.75);
  //leg1 -> SetHeader("AK4PFchs");
  leg1 -> AddEntry(ampf_norm, "MPF","L");
  leg1 -> AddEntry(adijet_norm, "Pt","L");
  leg1->Draw();

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  
  tex->DrawLatex(0.47,0.87,"Anti-k_{t} R = 0.4, PF+CHS");
  

  c1->SaveAs("KFSR_25ns_2p11fb.pdf");



}
