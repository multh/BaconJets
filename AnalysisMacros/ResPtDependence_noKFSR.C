#include "header.h"
#include "tdrstyle_mod14.C"

void ResPtDependence_noKFSR(TString path, double al_cut=0.2){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // get the residual root file for diff Pt
  TFile* res_nominal_file = new TFile(path+"Histo_Res_COMB_L1.root","READ");
  TFile* res_central_file = new TFile(path+"Histo_Res_COMB_L1_central.root","READ");
  TFile* res_down_file = new TFile(path+"Histo_Res_COMB_L1_down.root","READ");
  TFile* res_up_file = new TFile(path+"Histo_Res_COMB_L1_up.root","READ");
  TFile* res_doubleup_file = new TFile(path+"Histo_Res_COMB_L1_doubleup.root","READ");

  TH1D *res_nominal = (TH1D*)res_nominal_file->Get("res_logpt_comb");
  TH1D *res_central = (TH1D*)res_central_file->Get("res_logpt_comb");
  TH1D *res_down = (TH1D*)res_down_file->Get("res_logpt_comb");
  TH1D *res_up = (TH1D*)res_up_file->Get("res_logpt_comb");
  TH1D *res_doubleup = (TH1D*)res_doubleup_file->Get("res_logpt_comb");
  
   
  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,5.191,1);
  res_nominal->SetLineWidth(2);
  res_central->SetLineWidth(2);
  res_down->SetLineWidth(2);
  res_up->SetLineWidth(2);
  res_doubleup->SetLineWidth(2);
  res_nominal->SetLineColor(kRed+1);
  res_central->SetLineColor(kBlack);
  res_down->SetLineColor(kBlue+1);
  res_up->SetLineColor(kGreen+1);
  res_doubleup->SetLineColor(kOrange+7);

  res_nominal->GetXaxis()->SetTitle("|#eta|");
  res_nominal->GetXaxis()->SetTitleSize(0.05);
  res_nominal->GetXaxis()->SetTitleOffset(0.99);
  res_nominal->GetYaxis()->SetTitle("Relative correction");
  res_nominal->GetYaxis()->SetTitleSize(0.05);
  res_nominal->GetYaxis()->SetTitleOffset(1.3);

  // Draw results
  TH1D *h = new TH1D("h",";|#eta|;Relative correction",41,0,5.191);
  // h->SetMaximum(1.35);
  // h->SetMinimum(0.95);
  h->SetMaximum(1.2);
  h->SetMinimum(0.8);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  lumi_13TeV = "2.11 fb^{-1}";
  bool kSquare = true;

  TCanvas *c1 = tdrCanvas("c1",h,4,0,kSquare);
  res_nominal->GetYaxis()->SetRangeUser(0.8,1.3);
  res_nominal->Draw("E1");
  res_central->Draw("E1 SAME");
  res_down->Draw("E1 SAME");
  res_up->Draw("E1 SAME");
  res_doubleup->Draw("E1 SAME");
  line->SetLineStyle(2);
  line->Draw("SAME");

  TLegend *leg1 = tdrLeg(0.17,0.49,0.40,0.80);
  //  leg1 -> SetTitle("#bar{p_{T}}");
  leg1 -> AddEntry(res_down, "60 GeV","L");
  leg1 -> AddEntry(res_central, "120 GeV","L");
  leg1 -> AddEntry(res_up, "240 GeV","L");
  leg1 -> AddEntry(res_doubleup, "480 GeV","L");
  leg1 -> AddEntry(res_nominal, "Mean","L");
  leg1->Draw();

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  
  tex->DrawLatex(0.47,0.87,"Anti-k_{t} R = 0.4, PF+CHS");
  c1->SaveAs(path+"plots/RatioPtDepend_25ns_2p11fb_nokFSR.pdf");

  TFile* outputfilenorm = new TFile(path+"Histos_Res_COMB_L1_PtDepend.root","RECREATE");
  res_nominal->Write();
  res_central->Write();
  res_down->Write();
  res_up->Write();
  res_doubleup->Write();
  outputfilenorm->Write();
  outputfilenorm->Close();
}
