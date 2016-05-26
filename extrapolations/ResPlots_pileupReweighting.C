#include "header.h"
#include "tdrstyle_mod14.C"

void ResPlots_pileupReweighting(TString pathNo="/nfs/dust/cms/user/karavdia/JEC_76X/RESULTS_noReweight/",TString pathNonimal="/nfs/dust/cms/user/karavdia/JEC_76X/RESULTS_69mb/", TString pathBelow="/nfs/dust/cms/user/karavdia/JEC_76X/RESULTS_58mb/", TString pathAbove="/nfs/dust/cms/user/karavdia/JEC_76X/RESULTS_80mb/"){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);




  // get the residual root file for MPF and pt balance without pile up correction
  TFile* resmpf = new TFile(pathNo+"Histo_Res_MPF_L1.root","READ");
  TH1D* consmpf = (TH1D*)resmpf->Get("res_const_mpf");
  TH1D* ptmpf = (TH1D*)resmpf->Get("res_logpt_mpf");
  TFile* resdijet = new TFile(pathNo+"Histo_Res_DiJet_L1.root","READ"); 
  TH1D* consdijet = (TH1D*)resdijet->Get("res_const_dijet");
  TH1D* ptdijet = (TH1D*)resdijet->Get("res_logpt_dijet");


  //get the residual root file for MPF and pt balance with nominal (69 mb) correction
  TFile* resmpf_nom = new TFile(pathNonimal+"Histo_Res_MPF_L1.root","READ");
  TH1D* consmpf_nom = (TH1D*)resmpf_nom->Get("res_const_mpf");
  TH1D* ptmpf_nom = (TH1D*)resmpf_nom->Get("res_logpt_mpf");
  TFile* resdijet_nom = new TFile(pathNonimal+"Histo_Res_DiJet_L1.root","READ"); 
  TH1D* consdijet_nom = (TH1D*)resdijet_nom->Get("res_const_dijet");
  TH1D* ptdijet_nom = (TH1D*)resdijet_nom->Get("res_logpt_dijet");

//get the residual root file for MPF and pt balance with below nominal (58 mb) correction
  TFile* resmpf_b = new TFile(pathBelow+"Histo_Res_MPF_L1.root","READ");
  TH1D* consmpf_b = (TH1D*)resmpf_b->Get("res_const_mpf");
  TH1D* ptmpf_b = (TH1D*)resmpf_b->Get("res_logpt_mpf");
  TFile* resdijet_b = new TFile(pathBelow+"Histo_Res_DiJet_L1.root","READ"); 
  TH1D* consdijet_b = (TH1D*)resdijet_b->Get("res_const_dijet");
  TH1D* ptdijet_b = (TH1D*)resdijet_b->Get("res_logpt_dijet");

//get the residual root file for MPF and pt balance with above nominal (80 mb) correction
  TFile* resmpf_a = new TFile(pathAbove+"Histo_Res_MPF_L1.root","READ");
  TH1D* consmpf_a = (TH1D*)resmpf_a->Get("res_const_mpf");
  TH1D* ptmpf_a = (TH1D*)resmpf_a->Get("res_logpt_mpf");
  TFile* resdijet_a = new TFile(pathAbove+"Histo_Res_DiJet_L1.root","READ"); 
  TH1D* consdijet_a = (TH1D*)resdijet_a->Get("res_const_dijet");
  TH1D* ptdijet_a = (TH1D*)resdijet_a->Get("res_logpt_dijet");


  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,5.191,1);
  consmpf->SetLineWidth(3);
  consmpf->SetLineColor(kRed);
  consdijet->SetLineWidth(3);
  consdijet->SetLineColor(kBlue);
  consmpf_nom->SetLineWidth(4);
  consmpf_nom->SetLineColor(kMagenta);
  consdijet_nom->SetLineWidth(4);
  consdijet_nom->SetLineColor(kCyan);
  consmpf_b->SetLineWidth(2);
  consmpf_b->SetLineColor(kOrange+3);
  consdijet_b->SetLineWidth(2);
  consdijet_b->SetLineColor(kGreen+3);
  consmpf_a->SetLineWidth(2);
  consmpf_a->SetLineStyle(2);
  consmpf_a->SetLineColor(kOrange+7);
  consdijet_a->SetLineWidth(2);
  consdijet_a->SetLineStyle(2);
  consdijet_a->SetLineColor(kGreen);

  // consmpf->Draw("E1");
  // consmpf->GetXaxis()->SetTitle("|#eta|");
  // consmpf->GetXaxis()->SetTitleSize(0.05);
  // consmpf->GetXaxis()->SetTitleOffset(0.99);
  // consmpf->GetYaxis()->SetTitle("Relative correction");
  // consmpf->GetYaxis()->SetTitleSize(0.05);
  // consmpf->GetYaxis()->SetTitleOffset(1.40);

  // Draw results
  TH1D *h = new TH1D("h",";|#eta|;Relative correction",41,0,5.191);
  lumi_13TeV = "2.11 fb^{-1}";
  bool kSquare = true;
  TCanvas *c1 = tdrCanvas("c1",h,4,0,kSquare);

  // h->SetMaximum(1.35);
  //  h->SetMinimum(0.95);
  h->SetMaximum(1.1);
  h->SetMinimum(0.9);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
 
 
  //  TCanvas * c = new TCanvas("c", "c", 1000,1000);
  //gPad->SetGridy();

  // consmpf_norm->Draw("E1 SAME");
  // ptmpf_norm->Draw("E1 SAME");
  consmpf->Draw("E1 SAME");
  consdijet->Draw("E1 SAME");
  consmpf_nom->Draw("E1 SAME");
  consdijet_nom->Draw("E1 SAME");
  consmpf_b->Draw("E1 SAME");
  consdijet_b->Draw("E1 SAME");
  consmpf_a->Draw("E1 SAME");
  consdijet_a->Draw("E1 SAME");
  line->SetLineStyle(2);
  line->Draw("SAME");

 

  TLegend *leg1 = tdrLeg(0.17,0.19,0.40,0.33);
  leg1 -> AddEntry(consmpf, "MPF (no reweight)","L");
  leg1 -> AddEntry(consmpf_nom, "MPF (69mb)","L");
  leg1 -> AddEntry(consmpf_b, "MPF (58mb)","L");
  leg1 -> AddEntry(consmpf_a, "MPF (80mb)","L");
  leg1->SetTextSize(0.04);
  leg1->Draw();

 TLegend *leg2 = tdrLeg(0.17,0.69,0.40,0.83);
  leg2 -> AddEntry(consdijet, "Pt (no reweight)","L");
  leg2 -> AddEntry(consdijet_nom, "Pt (69mb)","L");
  leg2 -> AddEntry(consdijet_b, "Pt (58mb)","L");
  leg2 -> AddEntry(consdijet_a, "Pt (80mb)","L");
  leg2->SetTextSize(0.04);
  leg2->Draw();

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  
  tex->DrawLatex(0.47,0.87,"Anti-k_{t} R = 0.4, PF+CHS");
  c1->SaveAs("plots/L2ResDiff_25ns_2p11fb.pdf");
  





}
