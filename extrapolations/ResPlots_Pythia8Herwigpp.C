#include "header.h"
#include "tdrstyle_mod14.C"

void ResPlots_Pythia8Herwigpp(TString pathPythia="/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_noReweight/pythia8/",TString pathHerwigpp="/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_noReweight/herwigpp/"){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // get the residual root file for MPF and pt balance with pythia8
  TFile* resmpfPythia = new TFile(pathPythia+"Histo_Res_MPF_L1.root","READ");
  TFile* resdijetPythia = new TFile(pathPythia+"Histo_Res_DiJet_L1.root","READ"); 
  TH1D* consmpfPythia = (TH1D*)resmpfPythia->Get("res_const_mpf");
  TH1D* consdijetPythia = (TH1D*)resdijetPythia->Get("res_const_dijet");

// get the residual root file for MPF and pt balance with herwigpp
  TFile* resmpfHerwigpp = new TFile(pathHerwigpp+"Histo_Res_MPF_L1.root","READ");
  TFile* resdijetHerwigpp = new TFile(pathHerwigpp+"Histo_Res_DiJet_L1.root","READ"); 
  TH1D* consmpfHerwigpp = (TH1D*)resmpfHerwigpp->Get("res_const_mpf");
  TH1D* consdijetHerwigpp = (TH1D*)resdijetHerwigpp->Get("res_const_dijet");


  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,5.191,1);
  consmpfPythia->SetLineWidth(3);
  consmpfPythia->SetLineColor(kRed);
  consdijetPythia->SetLineWidth(3);
  consdijetPythia->SetLineColor(kBlue);
  consmpfHerwigpp->SetLineWidth(3);
  consmpfHerwigpp->SetLineStyle(3);
  consmpfHerwigpp->SetLineColor(kRed);
  consdijetHerwigpp->SetLineStyle(3);
  consdijetHerwigpp->SetLineWidth(3);
  consdijetHerwigpp->SetLineColor(kBlue);
  // Draw results
  TH1D *h = new TH1D("h",";|#eta|;Relative correction",41,0,5.191);
  bool kSquare = true;
  TCanvas *c1 = tdrCanvas("c1",h,4,0,kSquare);
  consmpfPythia->GetYaxis()->SetRangeUser(0.95,1.15);
  consmpfPythia->GetYaxis()->SetTitle("Relative correction");
  consmpfPythia->GetXaxis()->SetTitleSize(0.05);
  consmpfPythia->GetXaxis()->SetTitle("|#eta|");
  consmpfPythia->GetYaxis()->SetTitleSize(0.05);
  consmpfPythia->Draw("E1");
  consdijetPythia->Draw("E1 SAME");
  consmpfHerwigpp->Draw("E1 SAME");
  consdijetHerwigpp->Draw("E1 SAME");
  line->SetLineStyle(2);
  line->Draw("SAME");

 

  TLegend *leg1 = tdrLeg(0.18,0.59,0.40,0.83);
  leg1 -> AddEntry(consmpfPythia, "MPF (pythia8)","L");
  leg1 -> AddEntry(consmpfHerwigpp, "MPF (Herwigpp)","L");
  leg1 -> AddEntry(consdijetPythia, "Pt (pythia8)","L");
  leg1 -> AddEntry(consdijetHerwigpp, "Pt (Herwigpp)","L");
  leg1->SetTextSize(0.04);
  leg1->Draw();

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  
  tex->DrawLatex(0.47,0.87,"Anti-k_{t} R = 0.4, PF+CHS");
  c1->SaveAs("plots/L2Res_pythia8_vs_herwigpp_25ns_2p11fb.pdf");
  





}
