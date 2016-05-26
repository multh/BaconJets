#include "header.h"

void VarPtave(TString path){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // get the residual root file for MPF and pt balance
  TFile* resmpf = new TFile(path+"Histo_Res_MPF_L1_central.root","READ");
  TH1D* ptmpf = resmpf->Get("res_logpt_mpf");
  TFile* resmpf_half = new TFile(path+"Histo_Res_MPF_L1_down.root","READ"); 
  TH1D* ptmpf_half = resmpf_half->Get("res_logpt_mpf");
  TFile* resmpf_twice = new TFile(path+"Histo_Res_MPF_L1_up.root","READ"); 
  TH1D* ptmpf_twice = resmpf_twice->Get("res_logpt_mpf");
  TFile* resmpf_four = new TFile(path+"Histo_Res_MPF_L1_doubleup.root","READ"); 
  TH1D* ptmpf_four = resmpf_four->Get("res_logpt_mpf");


  // create histo for the normalization
  TH1D* ptmpf_norm=(TH1D*)ptmpf->Clone();
  TH1D* ptmpf_half_norm=(TH1D*)ptmpf_half->Clone();
  TH1D* ptmpf_twice_norm=(TH1D*)ptmpf_twice->Clone();
  TH1D* ptmpf_four_norm=(TH1D*)ptmpf_four->Clone();

  // create norm factor
  double norm_factor_ptMPF;
  double norm_factor_ptmpf_half;
  double norm_factor_ptmpf_twice;
  double norm_factor_ptmpf_four;
  for (int i=1; i<n_etabarr+1; i++){
    norm_factor_ptMPF+=ptmpf->GetBinContent(i);
    norm_factor_ptmpf_half+=ptmpf_half->GetBinContent(i);
    norm_factor_ptmpf_twice+=ptmpf_twice->GetBinContent(i);
    norm_factor_ptmpf_four+=ptmpf_four->GetBinContent(i);
  }
  norm_factor_ptMPF = norm_factor_ptMPF/n_etabarr;
  norm_factor_ptmpf_half = norm_factor_ptmpf_half/n_etabarr;
  norm_factor_ptmpf_twice = norm_factor_ptmpf_twice/n_etabarr;
  norm_factor_ptmpf_four = norm_factor_ptmpf_four/n_etabarr;

  // scale with norm factor
  ptmpf_norm->Scale(1/norm_factor_ptMPF);
  ptmpf_half_norm->Scale(1/norm_factor_ptmpf_half);
  ptmpf_twice_norm->Scale(1/norm_factor_ptmpf_twice);
  ptmpf_four_norm->Scale(1/norm_factor_ptmpf_four);

  TFile* outputfilenorm = new TFile(path+"Histos_Res_norm_L1.root","RECREATE");
  ptmpf_norm->Write();
  ptmpf_half_norm->Write();
  ptmpf_twice_norm->Write();
  ptmpf_four_norm->Write();
  outputfilenorm->Write();
  outputfilenorm->Close();


  //cout << "normalization factor (MPF): " << norm_factor_MPF << endl;
  //cout << "normalization factor (ptb): " << norm_factor_DiJet << endl;

  cout << endl;
  for (int i=1; i<n_etabarr+1; i++){
    cout << ptmpf_norm->GetBinContent(i) << endl;
  }


  TH1D* ptmpf_nominal=(TH1D*)ptmpf_norm->Clone();
  
  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,5.191,1);


  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  gROOT->LoadMacro("CMS_lumi.C");

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_13TeV  = "2.11 fb^{-1}"; // default is "19.7 fb^{-1}"
  //lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV 
  int iPos=11;

  
  // create a 'nice' plot..
  TCanvas * c = new TCanvas("c", "c", 1000,1000);
  gPad->SetTickx();
  gPad->SetTicky();
  ptmpf_norm->SetLineWidth(2);
  ptmpf_norm->SetLineColor(kRed+1);
  ptmpf_norm->SetMarkerColor(kRed+1);
  ptmpf_norm->SetMarkerStyle(21);
  ptmpf_norm->GetYaxis()->SetRangeUser(0.98,1.19);
  ptmpf_norm->GetXaxis()->SetTitle("|#eta|");
  ptmpf_norm->GetXaxis()->SetTitleSize(0.05);
  ptmpf_norm->GetXaxis()->SetTitleOffset(0.99);
  ptmpf_norm->GetYaxis()->SetTitle("Relative correction");
  ptmpf_norm->GetYaxis()->SetTitleSize(0.05);
  ptmpf_norm->GetYaxis()->SetTitleOffset(1.40);
  ptmpf_norm->Draw("E1");
  ptmpf_twice_norm->SetLineWidth(2);
  ptmpf_twice_norm->SetLineStyle(1);
  ptmpf_twice_norm->SetLineColor(kGreen+2);
  ptmpf_twice_norm->Draw("HIST SAME");
  ptmpf_half_norm->SetLineWidth(2);
  ptmpf_half_norm->SetLineStyle(1);
  ptmpf_half_norm->SetLineColor(kBlack);
  ptmpf_half_norm->Draw("HIST SAME");
  ptmpf_four_norm->SetLineWidth(2);
  ptmpf_four_norm->SetLineStyle(1);
  ptmpf_four_norm->SetLineColor(kBlue);
  line->SetLineStyle(kDashed);
  line->Draw("SAME");

  TLegend *leg1;
  leg1 = new TLegend(0.19,0.55,0.52,0.78,"","brNDC");//x+0.1
  leg1 -> SetBorderSize(0);
  leg1 -> SetTextSize(0.042);
  leg1 -> SetFillColor(10);
  leg1 -> SetLineColor(1);
  leg1 -> SetTextFont(62);
  leg1 -> SetHeader("p_{T}(jet)");
  leg1 -> AddEntry(ptmpf_half_norm, "0.5 #upoint #bar{p}_{T}","L");
  leg1 -> AddEntry(ptmpf_norm, "1 #upoint #bar{p}_{T}","L");
  leg1 -> AddEntry(ptmpf_twice_norm, "2 #upoint #bar{p}_{T}","L");
  leg1->Draw();
  TLegend *leg2;
  leg2 = new TLegend(0.51,0.85,0.91,0.94,"","brNDC");//x+0.1
  leg2 -> SetBorderSize(0);
  leg2 -> SetTextSize(0.042);
  leg2 -> SetFillColor(10);
  leg2 -> SetLineColor(1);
  leg2 -> SetTextFont(62);
  leg2 -> SetHeader("Anti-k_{T} R=0.4 PF+CHS");
  leg2 -> Draw();

  //  ptmpf_twice_norm->SetBinContent(16,1.);


  CMS_lumi( c, iPeriod, iPos );

  c->SaveAs("L2Res_25ns_2p11fb_ptave_systematics.pdf");
  
  gROOT->LoadMacro("tdrstyle_mod14.C");

  // Draw results
  TH1D *h = new TH1D("h",";|#eta|;Relative correction",41,0,5.191.);
  h->SetMaximum(1.19);
  h->SetMinimum(0.98);
  lumi_13TeV = "2.11 fb^{-1}";
  TCanvas *c1 = tdrCanvas("c1",h,4,0,kSquare);
  ptmpf_norm->Draw("SAMEP");
  ptmpf_half_norm->Draw("SAME HIST");
  ptmpf_twice_norm->Draw("SAME HIST");
  ptmpf_four_norm->Draw("SAME HIST");
  line->Draw("SAME");

  TLegend *leg1 = tdrLeg(0.17,0.36,0.50,0.77);
  leg1->SetHeader("p_{T}(jet)");
  leg1->AddEntry(ptmpf_half_norm, "60 GeV", "L");
  leg1->AddEntry(ptmpf_norm, "120 GeV", "P L");
  leg1->AddEntry(ptmpf_twice_norm, "240 GeV", "L");
  leg1->AddEntry(ptmpf_four_norm, "480 GeV", "L");
  leg1->Draw();

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  
  tex->DrawLatex(0.48,0.87,"Anti-k_{t} R = 0.4 PF+CHS");

  c1->SaveAs("L2Res_25ns_2p11fb_fixedPtave.pdf");


}
