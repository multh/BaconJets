#include "header.h"
#include "tdrstyle_mod14.C"

void AllResPlots_noKFSR(TString path){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // get the residual root file for MPF and pt balance
  TFile* resmpf = new TFile(path+"Histo_Res_MPF_L1.root","READ");
  TFile* resdijet = new TFile(path+"Histo_Res_DiJet_L1.root","READ"); 
  TFile* rescomb = new TFile(path+"Histo_Res_COMB_L1.root","READ"); 
  TH1D* consmpf = (TH1D*)resmpf->Get("res_logpt_mpf");
  TH1D* consdijet = (TH1D*)resdijet->Get("res_logpt_dijet");
  TH1D* conscomb = (TH1D*)rescomb->Get("res_logpt_comb");
  // consmpf->Print();
  // consdijet->Print();
  // conscomb->Print();

  // get the (R_{MC}/R_{DATA}) root file for MPF and pt balance
  TFile* rrmpf = new TFile(path+"Histo_ptave_MPF_L1.root","READ");
  rrmpf->Print();
  TFile* rrdijet = new TFile(path+"Histo_ptave_DiJet_L1.root","READ"); 
  rrdijet->Print();
  TFile* rrcomb = new TFile(path+"Histo_ptave_COMB_L1.root","READ"); 
  rrcomb->Print();
  TH1D* rrconstmpf = (TH1D*)rrmpf->Get("ptave_const_mpf");
  TH1D* rrconstdijet = (TH1D*)rrdijet->Get("ptave_const_dijet");
  TH1D* rrconstcomb = (TH1D*)rrcomb->Get("ptave_const_comb");
  TH1D* rrlogptmpf = (TH1D*)rrmpf->Get("ptave_logpt_mpf");
  TH1D* rrlogptdijet = (TH1D*)rrdijet->Get("ptave_logpt_dijet");
  TH1D* rrlogptcomb = (TH1D*)rrcomb->Get("ptave_logpt_comb");


  // create histo for the normalization
  TH1D* consmpf_norm=(TH1D*)consmpf->Clone();
  TH1D* consdijet_norm=(TH1D*)consdijet->Clone();
  TH1D* conscomb_norm=(TH1D*)conscomb->Clone();
  TH1D* rrconsmpf_norm=(TH1D*)rrconstmpf->Clone();
  TH1D* rrconsdijet_norm=(TH1D*)rrconstdijet->Clone();
  TH1D* rrconscomb_norm=(TH1D*)rrconstcomb->Clone();
  TH1D* rrlogptmpf_norm=(TH1D*)rrlogptmpf->Clone();
  TH1D* rrlogptdijet_norm=(TH1D*)rrlogptdijet->Clone();
  TH1D* rrlogptcomb_norm=(TH1D*)rrlogptcomb->Clone();

  // create norm factor
  double norm_factor_MPF;
  double norm_factor_DiJet;
  double norm_factor_COMB;
  double norm_factor_rrMPF;
  double norm_factor_rrDiJet;
  double norm_factor_rrCOMB;
  double norm_factor_rrMPFlog;
  double norm_factor_rrDiJetlog;
  double norm_factor_rrCOMBlog;


  for (int i=1; i<n_etabarr+1; i++){
    norm_factor_MPF+=consmpf->GetBinContent(i);
    norm_factor_DiJet+=consdijet->GetBinContent(i);
    norm_factor_COMB+=conscomb->GetBinContent(i);
    norm_factor_rrMPF+=rrconstmpf->GetBinContent(i);
    norm_factor_rrDiJet+=rrconstdijet->GetBinContent(i);
    norm_factor_rrCOMB+=rrconstcomb->GetBinContent(i);
    norm_factor_rrMPFlog+=rrlogptmpf->GetBinContent(i);
    norm_factor_rrDiJetlog+=rrlogptdijet->GetBinContent(i);
    norm_factor_rrCOMBlog+=rrlogptcomb->GetBinContent(i);

  }
  norm_factor_MPF = norm_factor_MPF/n_etabarr;
  norm_factor_DiJet = norm_factor_DiJet/n_etabarr;
  norm_factor_COMB = norm_factor_COMB/n_etabarr;
  norm_factor_rrMPF = norm_factor_rrMPF/n_etabarr;
  norm_factor_rrDiJet = norm_factor_rrDiJet/n_etabarr;
  norm_factor_rrCOMB = norm_factor_rrCOMB/n_etabarr;

  norm_factor_rrMPFlog = norm_factor_rrMPFlog/n_etabarr;
  norm_factor_rrDiJetlog = norm_factor_rrDiJetlog/n_etabarr;
  norm_factor_rrCOMBlog = norm_factor_rrCOMBlog/n_etabarr;

  // scale with norm factor
  consmpf_norm->Scale(1/norm_factor_MPF);
  consdijet_norm->Scale(1/norm_factor_DiJet);
  conscomb_norm->Scale(1/norm_factor_COMB);
  rrconsmpf_norm->Scale(1/norm_factor_rrMPF);
  rrconsdijet_norm->Scale(1/norm_factor_rrDiJet);
  rrconscomb_norm->Scale(1/norm_factor_rrCOMB);
  rrlogptmpf_norm->Scale(1/norm_factor_rrMPFlog);
  rrlogptdijet_norm->Scale(1/norm_factor_rrDiJetlog);
  rrlogptcomb_norm->Scale(1/norm_factor_rrCOMBlog);


  TFile* outputfilenorm = new TFile(path+"Histos_Res_norm_L1.root","RECREATE");
  consmpf_norm->Write();
  consdijet_norm->Write();
  conscomb_norm->Write();
  rrconsmpf_norm->Write();
  rrconsdijet_norm->Write();
  rrconscomb_norm->Write();
  rrlogptmpf_norm->Write();
  rrlogptdijet_norm->Write();
  rrlogptcomb_norm->Write();
  outputfilenorm->Write();
  outputfilenorm->Close();


  cout << "normalization factor (MPF): " << norm_factor_MPF << endl;
  cout << "normalization factor (ptb): " << norm_factor_DiJet << endl;
  cout << "normalization factor (comb): " << norm_factor_COMB << endl;

   
  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,5.191,1);
  //  gROOT->LoadMacro("tdrstyle_mod14.C");
  //  gROOT->ProcessLine(".L tdrstyle_mod14.C");
  consmpf_norm->SetLineWidth(2);
  consmpf_norm->SetLineColor(kRed+1);
  consdijet_norm->SetLineWidth(2);
  consdijet_norm->SetLineColor(kBlue+1);
  conscomb_norm->SetLineWidth(2);
  conscomb_norm->SetLineColor(kGreen+1);

  rrconsmpf_norm->SetLineWidth(2);
  rrconsmpf_norm->SetLineStyle(2);
  rrconsmpf_norm->SetLineColor(kRed+1);
  rrconsdijet_norm->SetLineWidth(2);
  rrconsdijet_norm->SetLineStyle(2);
  rrconsdijet_norm->SetLineColor(kBlue+1);
  rrconscomb_norm->SetLineWidth(2);
  rrconscomb_norm->SetLineStyle(2);
  rrconscomb_norm->SetLineColor(kGreen+1);
  rrlogptmpf_norm->SetLineWidth(2);
  rrlogptmpf_norm->SetLineColor(kRed+1);
  rrlogptdijet_norm->SetLineWidth(2);
  rrlogptdijet_norm->SetLineColor(kBlue+1);
  rrlogptcomb_norm->SetLineWidth(2);
  rrlogptcomb_norm->SetLineColor(kGreen+1);


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
  conscomb_norm->Draw("E1 SAME");
  line->SetLineStyle(2);
  line->Draw("SAME");

  TLegend *leg1 = tdrLeg(0.17,0.49,0.40,0.80);
  leg1 -> AddEntry(consmpf_norm, "MPF","L");
  leg1 -> AddEntry(consdijet_norm, "Pt","L");
  leg1 -> AddEntry(conscomb_norm, "COMB","L");
  leg1->Draw();

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  
  tex->DrawLatex(0.47,0.87,"Anti-k_{t} R = 0.4, PF+CHS");
  

  c1->SaveAs("plots/L2Res_25ns_2p11fb.pdf");
  

  TCanvas *c2 = tdrCanvas("c2",h,4,0,kSquare);
  rrconsmpf_norm->GetYaxis()->SetTitle("(R^{MC}/R^{data})_{#alpha<0.3}");
  rrconsmpf_norm->GetYaxis()->SetTitleSize(0.05);
  rrconsmpf_norm->GetYaxis()->SetRangeUser(0.95,1.15);
  rrconsmpf_norm->GetXaxis()->SetTitle("|#eta|");
  rrconsmpf_norm->GetXaxis()->SetTitleSize(0.05);
  rrconsmpf_norm->Draw("E1");
  rrconsdijet_norm->Draw("E1 SAME");
  rrconscomb_norm->Draw("E1 SAME");
  rrlogptmpf_norm->Draw("E1 SAME");
  rrlogptdijet_norm->Draw("E1 SAME");
  rrlogptcomb_norm->Draw("E1 SAME");
  line->Draw("SAME");
  TLegend *leg2 = tdrLeg(0.17,0.49,0.40,0.80);
  leg2 -> AddEntry(rrconsmpf_norm, "MPF FLAT","L");
  leg2 -> AddEntry(rrconsdijet_norm, "Pt FLAT","L");
  leg2 -> AddEntry(rrconscomb_norm, "COMB FLAT","L");
  leg2 -> AddEntry(rrlogptmpf_norm, "MPF LOGLIN","L");
  leg2 -> AddEntry(rrlogptdijet_norm, "Pt LOGLIN","L");
  leg2 -> AddEntry(rrlogptcomb_norm, "COMB LOGLIN","L");
  leg2->Draw();
  tex->DrawLatex(0.47,0.87,"Anti-k_{t} R = 0.4, PF+CHS");
  c2->SaveAs("plots/Ratio_25ns_2p11fb.pdf");
}
