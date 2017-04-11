#include "header.h"
#include "tdrstyle_mod14.C"

void AllResPlots(TString path){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // get the residual root file for MPF and pt balance
  TFile* resmpf = new TFile(path+"Histo_Res_MPF_L1.root","READ");
  TH1D* consmpf = (TH1D*)resmpf->Get("res_const_mpf");
  TFile* resdijet = new TFile(path+"Histo_Res_DiJet_L1.root","READ"); 
  TH1D* consdijet = (TH1D*)resdijet->Get("res_const_dijet");
  consmpf->Print();
  consdijet->Print();

  // get the (R_{MC}/R_{DATA}) root file for MPF and pt balance
  TFile* rrmpf = new TFile(path+"Histo_ptave_MPF_L1.root","READ");
  TH1D* rrconstmpf = (TH1D*)rrmpf->Get("ptave_const_mpf");
  TFile* rrdijet = new TFile(path+"Histo_ptave_DiJet_L1.root","READ"); 
  TH1D* rrconstdijet = (TH1D*)rrdijet->Get("ptave_const_dijet");
  rrconstmpf->Print();
  rrconstdijet->Print();

  // get the kFSR root file for MPF and pt balance
  TFile* kfsrmpf = new TFile(path+"Histo_KFSR_MPF_L1.root","READ");
  TH1D* kfsrconstmpf = (TH1D*)kfsrmpf->Get("kfsr_mpf");
  TFile* kfsrdijet = new TFile(path+"Histo_KFSR_DiJet_L1.root","READ"); 
  TH1D* kfsrconstdijet = (TH1D*)kfsrdijet->Get("kfsr_dijet");
  kfsrconstmpf->Print();
  kfsrconstdijet->Print();

  // create histo for the normalization
  TH1D* consmpf_norm=(TH1D*)consmpf->Clone();
  TH1D* consdijet_norm=(TH1D*)consdijet->Clone();
  TH1D* rrconsmpf_norm=(TH1D*)rrconstmpf->Clone();
  TH1D* rrconsdijet_norm=(TH1D*)rrconstdijet->Clone();
  TH1D* kfsrconsmpf_norm=(TH1D*)kfsrconstmpf->Clone();
  TH1D* kfsrconsdijet_norm=(TH1D*)kfsrconstdijet->Clone();

  // create norm factor
  double norm_factor_MPF;
  double norm_factor_DiJet;
  double norm_factor_rrMPF;
  double norm_factor_rrDiJet;
  double norm_factor_kfsrMPF;
  double norm_factor_kfsrDiJet;

  for (int i=1; i<n_etabarr+1; i++){
    norm_factor_MPF+=consmpf->GetBinContent(i);
    norm_factor_DiJet+=consdijet->GetBinContent(i);
    norm_factor_rrMPF+=rrconstmpf->GetBinContent(i);
    norm_factor_rrDiJet+=rrconstdijet->GetBinContent(i);
    norm_factor_kfsrMPF+=kfsrconstmpf->GetBinContent(i);
    norm_factor_kfsrDiJet+=kfsrconstdijet->GetBinContent(i);
  }
  norm_factor_MPF = norm_factor_MPF/n_etabarr;
  norm_factor_DiJet = norm_factor_DiJet/n_etabarr;
  norm_factor_rrMPF = norm_factor_rrMPF/n_etabarr;
  norm_factor_rrDiJet = norm_factor_rrDiJet/n_etabarr;
  norm_factor_kfsrMPF = norm_factor_kfsrMPF/n_etabarr;
  norm_factor_kfsrDiJet = norm_factor_kfsrDiJet/n_etabarr;

  // scale with norm factor
  consmpf_norm->Scale(1/norm_factor_MPF);
  consdijet_norm->Scale(1/norm_factor_DiJet);
  rrconsmpf_norm->Scale(1/norm_factor_rrMPF);
  rrconsdijet_norm->Scale(1/norm_factor_rrDiJet);
  kfsrconsmpf_norm->Scale(1/norm_factor_kfsrMPF);
  kfsrconsdijet_norm->Scale(1/norm_factor_kfsrDiJet);

  TFile* outputfilenorm = new TFile(path+"Histos_Res_norm_L1.root","RECREATE");
  consmpf_norm->Write();
  consdijet_norm->Write();
  rrconsmpf_norm->Write();
  rrconsdijet_norm->Write();
  kfsrconsmpf_norm->Write();
  kfsrconsdijet_norm->Write();
  outputfilenorm->Write();
  outputfilenorm->Close();


  cout << "normalization factor (MPF): " << norm_factor_MPF << endl;
  cout << "normalization factor (ptb): " << norm_factor_DiJet << endl;

  cout << endl;
  for (int i=1; i<n_etabarr+1; i++){
    cout << consmpf_norm->GetBinContent(i) << endl;
  }


  
  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,5.191,1);
  //  gROOT->LoadMacro("tdrstyle_mod14.C");
  //  gROOT->ProcessLine(".L tdrstyle_mod14.C");
  consmpf_norm->SetLineWidth(2);
  consmpf_norm->SetLineColor(kRed+1);
  consdijet_norm->SetLineWidth(2);
  consdijet_norm->SetLineColor(kBlue+1);
  rrconsmpf_norm->SetLineWidth(2);
  rrconsmpf_norm->SetLineColor(kRed+1);
  rrconsdijet_norm->SetLineWidth(2);
  rrconsdijet_norm->SetLineColor(kBlue+1);
  kfsrconsmpf_norm->SetLineWidth(2);
  kfsrconsmpf_norm->SetLineColor(kRed+1);
  kfsrconsdijet_norm->SetLineWidth(2);
  kfsrconsdijet_norm->SetLineColor(kBlue+1);

  //  consdijet_norm->Draw("E1 SAME");
//  consmpf_norm->Draw("E1");
  consmpf_norm->GetXaxis()->SetTitle("|#eta|");
  consmpf_norm->GetXaxis()->SetTitleSize(0.05);
  consmpf_norm->GetXaxis()->SetTitleOffset(0.99);
  consmpf_norm->GetYaxis()->SetTitle("Relative correction");
  consmpf_norm->GetYaxis()->SetTitleSize(0.05);
  consmpf_norm->GetYaxis()->SetTitleOffset(1.40);

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
  consmpf_norm->GetYaxis()->SetRangeUser(0.95,1.15);
  consmpf_norm->Draw("E1");
  consdijet_norm->Draw("E1 SAME");
  line->SetLineStyle(2);
  line->Draw("SAME");

  TLegend *leg1 = tdrLeg(0.17,0.49,0.40,0.80);
  leg1 -> AddEntry(consmpf_norm, "MPF","L");
  leg1 -> AddEntry(consdijet_norm, "Pt","L");
  leg1->Draw();

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  
  tex->DrawLatex(0.47,0.87,"Anti-k_{t} R = 0.4, PF+CHS");
  

  c1->SaveAs("plots/L2Res_25ns_2p11fb.pdf");
  

  TCanvas *c2 = tdrCanvas("c2",h,4,0,kSquare);
  rrconsmpf_norm->GetYaxis()->SetTitle("(R^{MC}/R^{data})_{#alpha<0.2}");
  rrconsmpf_norm->GetYaxis()->SetTitleSize(0.05);
  rrconsmpf_norm->GetYaxis()->SetRangeUser(0.95,1.15);
  rrconsmpf_norm->GetXaxis()->SetTitle("|#eta|");
  rrconsmpf_norm->GetXaxis()->SetTitleSize(0.05);
  rrconsmpf_norm->Draw("E1");
  rrconsdijet_norm->Draw("E1 SAME");
  line->Draw("SAME");
  leg1->Draw();
  tex->DrawLatex(0.47,0.87,"Anti-k_{t} R = 0.4, PF+CHS");
  c2->SaveAs("plots/Ratio_25ns_2p11fb.pdf");

  TCanvas *c3 = tdrCanvas("c3",h,4,0,kSquare);
  kfsrconsmpf_norm->GetYaxis()->SetTitle("k_{FSR}");
  kfsrconsmpf_norm->GetYaxis()->SetTitleSize(0.05);
  kfsrconsmpf_norm->GetYaxis()->SetRangeUser(0.998,1.002);
  //kfsrconsmpf_norm->GetYaxis()->SetRangeUser(0.95,1.15);
  kfsrconsmpf_norm->GetXaxis()->SetTitle("|#eta|");
  kfsrconsmpf_norm->GetXaxis()->SetTitleSize(0.05);
  kfsrconsmpf_norm->Draw("E1");
  kfsrconsdijet_norm->Draw("E1 SAME");
  line->Draw("SAME");
  leg1->Draw();
  tex->DrawLatex(0.47,0.87,"Anti-k_{t} R = 0.4, PF+CHS");
  c3->SaveAs("plots/kFSR_25ns_2p11fb.pdf");
}
