#include "header.h"
#include "tdrstyle_mod14.C"

void ResPlots_addCOMB(TString path){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // get the residual root file for MPF and pt balance
  TFile* resmpf = new TFile(path+"Histo_Res_MPF_L1.root","READ");
  TH1D* consmpf = (TH1D*)resmpf->Get("res_const_mpf");
  TH1D* ptmpf = (TH1D*)resmpf->Get("res_logpt_mpf");
  TFile* resdijet = new TFile(path+"Histo_Res_DiJet_L1.root","READ"); 
  TH1D* consdijet = (TH1D*)resdijet->Get("res_const_dijet");
  TH1D* ptdijet = (TH1D*)resdijet->Get("res_logpt_dijet");

  TFile* rescomb = new TFile(path+"Histo_Res_COMB_L1.root","READ"); 
  TH1D* conscomb = (TH1D*)rescomb->Get("res_const_comb");
  TH1D* ptcomb = (TH1D*)rescomb->Get("res_logpt_comb");


  // create histo for the normalization
  TH1D* consmpf_norm=(TH1D*)consmpf->Clone();
  TH1D* consdijet_norm=(TH1D*)consdijet->Clone();
  TH1D* conscomb_norm=(TH1D*)conscomb->Clone();
  TH1D* ptmpf_norm=(TH1D*)ptmpf->Clone();
  TH1D* ptdijet_norm=(TH1D*)ptdijet->Clone();
  TH1D* ptcomb_norm=(TH1D*)ptcomb->Clone();

  // create norm factor
  double norm_factor_MPF;
  double norm_factor_DiJet;
  double norm_factor_comb;
  double norm_factor_ptMPF;
  double norm_factor_ptDiJet;
  double norm_factor_ptcomb;
  for (int i=1; i<n_etabarr+1; i++){
    norm_factor_MPF+=consmpf->GetBinContent(i);
    norm_factor_DiJet+=consdijet->GetBinContent(i);
    norm_factor_comb+=conscomb->GetBinContent(i);
    norm_factor_ptMPF+=ptmpf->GetBinContent(i);
    norm_factor_ptDiJet+=ptdijet->GetBinContent(i);
    norm_factor_ptcomb+=ptcomb->GetBinContent(i);
  }
  norm_factor_MPF = norm_factor_MPF/n_etabarr;
  norm_factor_DiJet = norm_factor_DiJet/n_etabarr;
  norm_factor_comb = norm_factor_comb/n_etabarr;
  norm_factor_ptMPF = norm_factor_ptMPF/n_etabarr;
  norm_factor_ptDiJet = norm_factor_ptDiJet/n_etabarr;
  norm_factor_ptcomb = norm_factor_ptcomb/n_etabarr;

  // scale with norm factor
  consmpf_norm->Scale(1/norm_factor_MPF);
  consdijet_norm->Scale(1/norm_factor_DiJet);
  conscomb_norm->Scale(1/norm_factor_comb);
  ptmpf_norm->Scale(1/norm_factor_ptMPF);
  ptdijet_norm->Scale(1/norm_factor_ptDiJet);
  ptcomb_norm->Scale(1/norm_factor_ptcomb);


  TFile* outputfilenorm = new TFile(path+"Histos_Res_norm_L1.root","RECREATE");
  consmpf_norm->Write();
  consdijet_norm->Write();
  conscomb_norm->Write();
  ptmpf_norm->Write();
  ptdijet_norm->Write();
  ptcomb_norm->Write();
  outputfilenorm->Write();
  outputfilenorm->Close();


  cout << "normalization factor (MPF): " << norm_factor_MPF << endl;
  cout << "normalization factor (ptb): " << norm_factor_DiJet << endl;
  cout << "normalization factor (comb): " << norm_factor_comb << endl;

  // cout << endl;
  // for (int i=1; i<n_etabarr+1; i++){
  //   cout << consmpf_norm->GetBinContent(i) << endl;
  // }


  // ofstream output, output_complete;
  // ofstream output_dijet, output_complete_dijet;
  // output.open("L2res_MPF_V3_norm.txt");
  // output_complete.open("L2res_MPF_complete_norm.txt");
  // output_dijet.open("L2res_DiJet_V3_norm.txt");
  // output_complete_dijet.open("L2res_DiJet_complete_norm.txt");
  // output  << "JetEta (abs)     " << " Correction " << "  statistical unc" << endl;
  // output_dijet  << "JetEta (abs)     " << " Correction " << "  statistical unc" << endl;
  // output_complete  << "JetEta (abs)     " << " Correction " << "  statistical unc" << endl;
  // output_complete_dijet  << "JetEta (abs)     " << " Correction " << "  statistical unc" << endl;

  // for (int j=1; j<n_eta; j++){
  //   output << std::setprecision(6)  << eta_range[j-1]<< "    " << eta_range[j] << "    " << consmpf_norm->GetBinContent(j) << "        " << consmpf_norm->GetBinError(j) << endl;
  //   output_dijet << std::setprecision(6)  << eta_range[j-1]<< "    " << eta_range[j] << "    " << consdijet_norm->GetBinContent(j) << "        " << consdijet_norm->GetBinError(j) << endl;
  //   output_complete << std::setprecision(6)  << eta_range[j-1]<< "    " << eta_range[j] << "    " << ptmpf_norm->GetBinContent(j) << "        " << ptmpf_norm->GetBinError(j) << endl;
  //   output_complete_dijet << std::setprecision(6)  << eta_range[j-1]<< "    " << eta_range[j] << "    " << ptdijet_norm->GetBinContent(j) << "        " << ptdijet_norm->GetBinError(j) << endl;
  // }

  
  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,5.191,1);
  //  gROOT->LoadMacro("tdrstyle_mod14.C");
  //  gROOT->ProcessLine(".L tdrstyle_mod14.C");
  consmpf_norm->SetLineWidth(2);
  consmpf_norm->SetLineColor(kRed+1);
  consmpf_norm->Draw("E1");
  ptmpf_norm->SetLineWidth(2);
  ptmpf_norm->SetLineStyle(2);
  ptmpf_norm->SetLineColor(kRed+1);
  ptmpf_norm->Draw("E1 SAME");
  consdijet_norm->SetLineWidth(2);
  consdijet_norm->SetLineColor(kBlue+1);
  consdijet_norm->Draw("E1 SAME");
  ptdijet_norm->SetLineWidth(2);
  ptdijet_norm->SetLineStyle(2);
  ptdijet_norm->SetLineColor(kBlue+1);
  ptdijet_norm->Draw("E1 SAME");
  conscomb_norm->SetLineWidth(2);
  conscomb_norm->SetLineColor(kGreen+4);
  conscomb_norm->Draw("E1 SAME");
  ptcomb_norm->SetLineWidth(2);
  ptcomb_norm->SetLineStyle(2);
  ptcomb_norm->SetLineColor(kGreen+4);
  ptcomb_norm->Draw("E1 SAME");


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
  h->SetMaximum(1.15);
  h->SetMinimum(0.95);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  lumi_13TeV = "2.11 fb^{-1}";
  bool kSquare = true;
  TCanvas *c1 = tdrCanvas("c1",h,4,0,kSquare);
  //  TCanvas * c = new TCanvas("c", "c", 1000,1000);
  //gPad->SetGridy();

  consmpf_norm->Draw("E1 SAME");
  ptmpf_norm->Draw("E1 SAME");
  consdijet_norm->Draw("E1 SAME");
  ptdijet_norm->Draw("E1 SAME");
  conscomb_norm->Draw("E1 SAME");
  ptcomb_norm->Draw("E1 SAME");
  line->SetLineStyle(2);
  line->Draw("SAME");

  TLegend *leg1 = tdrLeg(0.17,0.49,0.40,0.80);
  //leg1 -> SetHeader("AK4PFchs");
  leg1 -> AddEntry(consmpf_norm, "MPF FLAT","L");
  leg1 -> AddEntry(ptmpf_norm, "MPF LOGLIN","L");
  leg1 -> AddEntry(consdijet_norm, "Pt FLAT","L");
  leg1 -> AddEntry(ptdijet_norm, "Pt LOGLIN","L");
  leg1 -> AddEntry(conscomb_norm, "COMB FLAT","L");
  leg1 -> AddEntry(ptcomb_norm, "COMB LOGLIN","L");
  leg1->Draw();

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  
  tex->DrawLatex(0.47,0.87,"Anti-k_{t} R = 0.4, PF+CHS");
  

  c1->SaveAs("plots/L2Res_25ns_2p11fb_addCOMB.pdf");
  





}
