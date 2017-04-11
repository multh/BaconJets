#include "header.h"
#include "tdrstyle_mod14.C"

void ResPlots_OldNew(TString path){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // get the residual root file for MPF and pt balance
  TFile* resmpf = new TFile(path+"Histo_Res_MPF_L1.root","READ");
  TH1D* consmpf = (TH1D*)resmpf->Get("res_const_mpf");
  TH1D* ptmpf = (TH1D*)resmpf->Get("res_logpt_mpf");
  TFile* resdijet = new TFile(path+"Histo_Res_DiJet_L1.root","READ"); 
  TH1D* consdijet = (TH1D*)resdijet->Get("res_const_dijet");
  TH1D* ptdijet = (TH1D*)resdijet->Get("res_logpt_dijet");



// get OLD the residual root file for MPF and pt balance
  TFile* resmpf_old = new TFile("/nfs/dust/cms/user/kovalch/sFrame/JEC/V6/L2ResPlots/Histo_Res_MPF_25ns_data_V6_new_trigg60for80only2.root","READ");
  TH1D* consmpf_old = (TH1D*)resmpf_old->Get("res_const_mpf");
  TH1D* ptmpf_old = (TH1D*)resmpf_old->Get("res_logpt_mpf");
  TFile* resdijet_old = new TFile("/nfs/dust/cms/user/kovalch/sFrame/JEC/V6/L2ResPlots/Histo_Res_DiJet_25ns_data_V6_new_trigg60for80only2.root","READ"); 
  TH1D* consdijet_old = (TH1D*)resdijet_old->Get("res_const_dijet");
  TH1D* ptdijet_old = (TH1D*)resdijet_old->Get("res_logpt_dijet");


  // create histo for the normalization
  TH1D* consmpf_norm=(TH1D*)consmpf->Clone();
  TH1D* consdijet_norm=(TH1D*)consdijet->Clone();
  TH1D* ptmpf_norm=(TH1D*)ptmpf->Clone();
  TH1D* ptdijet_norm=(TH1D*)ptdijet->Clone();

  // create norm factor
  double norm_factor_MPF;
  double norm_factor_DiJet;
  double norm_factor_ptMPF;
  double norm_factor_ptDiJet;
  for (int i=1; i<n_etabarr+1; i++){
    norm_factor_MPF+=consmpf->GetBinContent(i);
    norm_factor_DiJet+=consdijet->GetBinContent(i);
    norm_factor_ptMPF+=ptmpf->GetBinContent(i);
    norm_factor_ptDiJet+=ptdijet->GetBinContent(i);
  }
  norm_factor_MPF = norm_factor_MPF/n_etabarr;
  norm_factor_DiJet = norm_factor_DiJet/n_etabarr;
  norm_factor_ptMPF = norm_factor_ptMPF/n_etabarr;
  norm_factor_ptDiJet = norm_factor_ptDiJet/n_etabarr;

  // scale with norm factor
  consmpf_norm->Scale(1/norm_factor_MPF);
  consdijet_norm->Scale(1/norm_factor_DiJet);
  ptmpf_norm->Scale(1/norm_factor_ptMPF);
  ptdijet_norm->Scale(1/norm_factor_ptDiJet);


  TFile* outputfilenorm = new TFile(path+"Histos_Res_norm_L1.root","RECREATE");
  consmpf_norm->Write();
  consdijet_norm->Write();
  ptmpf_norm->Write();
  ptdijet_norm->Write();
  outputfilenorm->Write();
  outputfilenorm->Close();


  cout << "normalization factor (MPF): " << norm_factor_MPF << endl;
  cout << "normalization factor (ptb): " << norm_factor_DiJet << endl;

  cout << endl;
  for (int i=1; i<n_etabarr+1; i++){
    cout << consmpf_norm->GetBinContent(i) << endl;
  }


  ofstream output, output_complete;
  ofstream output_dijet, output_complete_dijet;
  output.open("L2res_MPF_V3_norm.txt");
  output_complete.open("L2res_MPF_complete_norm.txt");
  output_dijet.open("L2res_DiJet_V3_norm.txt");
  output_complete_dijet.open("L2res_DiJet_complete_norm.txt");
  output  << "JetEta (abs)     " << " Correction " << "  statistical unc" << endl;
  output_dijet  << "JetEta (abs)     " << " Correction " << "  statistical unc" << endl;
  output_complete  << "JetEta (abs)     " << " Correction " << "  statistical unc" << endl;
  output_complete_dijet  << "JetEta (abs)     " << " Correction " << "  statistical unc" << endl;

  for (int j=1; j<n_eta; j++){
    output << std::setprecision(6)  << eta_range[j-1]<< "    " << eta_range[j] << "    " << consmpf_norm->GetBinContent(j) << "        " << consmpf_norm->GetBinError(j) << endl;
    output_dijet << std::setprecision(6)  << eta_range[j-1]<< "    " << eta_range[j] << "    " << consdijet_norm->GetBinContent(j) << "        " << consdijet_norm->GetBinError(j) << endl;
    output_complete << std::setprecision(6)  << eta_range[j-1]<< "    " << eta_range[j] << "    " << ptmpf_norm->GetBinContent(j) << "        " << ptmpf_norm->GetBinError(j) << endl;
    output_complete_dijet << std::setprecision(6)  << eta_range[j-1]<< "    " << eta_range[j] << "    " << ptdijet_norm->GetBinContent(j) << "        " << ptdijet_norm->GetBinError(j) << endl;
  }




  
  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,5.191,1);
  //  gROOT->LoadMacro("tdrstyle_mod14.C");
  //  gROOT->ProcessLine(".L tdrstyle_mod14.C");
  // consmpf_norm->SetLineWidth(2);
  // consmpf_norm->SetLineColor(kRed+1);
  // consmpf_norm->Draw("E1");
  // ptmpf_norm->SetLineWidth(2);
  // ptmpf_norm->SetLineStyle(2);
  // ptmpf_norm->SetLineColor(kRed+1);
  // ptmpf_norm->Draw("E1 SAME");
  consdijet_norm->SetLineWidth(2);
  consdijet_norm->SetLineColor(kGreen+1);
  consdijet_norm->Draw("E1");
  //  consdijet_norm->Draw("E1 SAME");
  ptdijet_norm->SetLineWidth(2);
  ptdijet_norm->SetLineStyle(2);
  ptdijet_norm->SetLineColor(kGreen+1);
  ptdijet_norm->Draw("E1 SAME");


  consmpf_norm->GetXaxis()->SetTitle("|#eta|");
  consmpf_norm->GetXaxis()->SetTitleSize(0.05);
  consmpf_norm->GetXaxis()->SetTitleOffset(0.99);
  consmpf_norm->GetYaxis()->SetTitle("Relative correction");
  consmpf_norm->GetYaxis()->SetTitleSize(0.05);
  consmpf_norm->GetYaxis()->SetTitleOffset(1.40);

  // Draw results
  TH1D *h = new TH1D("h",";|#eta|;Relative correction",41,0,5.191);
  // h->SetMaximum(1.35);
  //  h->SetMinimum(0.95);
  h->SetMaximum(1.4);
  h->SetMinimum(0.8);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  lumi_13TeV = "2.11 fb^{-1}";
  bool kSquare = true;
  TCanvas *c1 = tdrCanvas("c1",h,4,0,kSquare);
  //  TCanvas * c = new TCanvas("c", "c", 1000,1000);
  //gPad->SetGridy();

  // consmpf_norm->Draw("E1 SAME");
  // ptmpf_norm->Draw("E1 SAME");
  consdijet_norm->Draw("E1 SAME");
  ptdijet_norm->Draw("E1 SAME");
  line->SetLineStyle(2);
  line->Draw("SAME");

  consmpf_old->SetLineColor(kRed);
  consmpf_old->Draw("E1 SAME");
  ptmpf_old->SetLineColor(kRed);
  ptmpf_old->SetLineWidth(2);
  ptmpf_old->SetLineStyle(2);
  ptmpf_old->Draw("E1 SAME");
  consdijet_old->SetLineColor(kBlue+1);
  consdijet_old->SetLineWidth(2);
  consdijet_old->Draw("E1 SAME");
  ptdijet_old->SetLineColor(kBlue+1);
  ptdijet_old->SetLineWidth(2);
  ptdijet_old->SetLineStyle(2);
  ptdijet_old->Draw("E1 SAME");

  TLegend *leg1 = tdrLeg(0.17,0.49,0.40,0.80);
  //leg1 -> SetHeader("AK4PFchs");
  // leg1 -> AddEntry(consmpf_norm, "MPF FLAT","L");
  // leg1 -> AddEntry(ptmpf_norm, "MPF LOGLIN","L");
  leg1 -> AddEntry(consdijet_norm, "Pt FLAT (76X)","L");
  leg1 -> AddEntry(ptdijet_norm, "Pt LOGLIN (76X)","L");
  leg1 -> AddEntry(consdijet_old, "Pt FLAT (74X)","L");
  leg1 -> AddEntry(ptdijet_old, "Pt LOGLIN (74X)","L");
  leg1 -> AddEntry(consmpf_old, "MPF FLAT (74X)","L");
  leg1 -> AddEntry(ptmpf_old, "MPF LOGLIN (74X)","L");

  leg1->Draw();

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  
  tex->DrawLatex(0.47,0.87,"Anti-k_{t} R = 0.4, PF+CHS");
  

  c1->SaveAs("plots/L2Res_25ns_2p11fb.pdf");
  





}
