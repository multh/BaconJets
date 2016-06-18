#include "header.h"
//#include "tdrstyle_mod14.C"
#include "tdrstyle_mod15.C"

void AllResPlots(TString path, TString txttag, TString jettag, TString variation, TString tag, double al_cut=0.2){

  int syst = 0;
  if(variation=="central") syst=0;
  if(variation=="down") syst=1;
  if(variation=="up") syst=2;
  if(variation=="doubleup") syst=3;
  if(variation=="nominal") syst=5;
  //  gROOT->LoadMacro("tdrstyle_mod15.C");
  //  setTDRStyle();
  // gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
    TString JetDescrib = "Anti-k_{t} R = 0.4, PF+CHS";
  //  TString JetDescrib = "Anti-k_{t} R = 0.4, PF+PUPPI";
  //  TString JetDescrib = "Anti-k_{t} R = 0.8, PF+CHS";
  //  TString JetDescrib = "Anti-k_{t} R = 0.8, PF+PUPPI";

  // get the residual root file for MPF and pt balance
  //  TFile* resmpf = new TFile(path+"Histo_Res_MPF_L1.root","READ");

  TFile* resmpf;
    if(syst==5){
      resmpf = new TFile(path+"Histo_Res_MPF_L1"+tag+".root","READ");
    }
    if(syst==0){
     resmpf = new TFile(path+"Histo_Res_MPF_L1_"+variation+".root","READ");
    }
    if(syst==1){
     resmpf = new TFile(path+"Histo_Res_MPF_L1_"+variation+".root","READ");
    }
    if(syst==2){
      resmpf = new TFile(path+"Histo_Res_MPF_L1_"+variation+".root","READ");
    }
    if(syst==3){
      resmpf = new TFile(path+"Histo_Res_MPF_L1_"+variation+".root","READ");
    }
  TH1D* consmpf = (TH1D*)resmpf->Get("res_const_mpf");
  TH1D* logptmpf = (TH1D*)resmpf->Get("res_logpt_mpf");

  //  TFile* resdijet = new TFile(path+"Histo_Res_DiJet_L1.root","READ"); 
  TFile* resdijet;
  if(syst==5){
    resdijet = new TFile(path+"Histo_Res_DiJet_L1"+tag+".root","READ");
  }
  if(syst==0){
    resdijet = new TFile(path+"Histo_Res_DiJet_L1_"+variation+".root","READ");
  }
  if(syst==1){
    resdijet = new TFile(path+"Histo_Res_DiJet_L1_"+variation+".root","READ");
  }
  if(syst==2){
    resdijet = new TFile(path+"Histo_Res_DiJet_L1_"+variation+".root","READ");
  }
  if(syst==3){
    resdijet = new TFile(path+"Histo_Res_DiJet_L1_"+variation+".root","READ");
  }
  TH1D* consdijet = (TH1D*)resdijet->Get("res_const_dijet");
  TH1D* logptdijet = (TH1D*)resdijet->Get("res_logpt_dijet");
  // consmpf->Print();
  // consdijet->Print();

  // get the (R_{MC}/R_{DATA}) root file for MPF and pt balance
  TFile* rrmpf;// = new TFile(path+"Histo_ptave_MPF_L1.root","READ");
 if(syst==5){
    rrmpf = new TFile(path+"Histo_ptave_MPF_L1"+tag+".root","READ");
  }
  if(syst==0){
    rrmpf = new TFile(path+"Histo_ptave_MPF_L1_"+variation+".root","READ");
  }
  if(syst==1){
    rrmpf = new TFile(path+"Histo_ptave_MPF_L1_"+variation+".root","READ");
  }
  if(syst==2){
    rrmpf = new TFile(path+"Histo_ptave_MPF_L1_"+variation+".root","READ");
  }
  if(syst==3){
    rrmpf = new TFile(path+"Histo_ptave_MPF_L1_"+variation+".root","READ");
  }
  TH1D* rrconstmpf = (TH1D*)rrmpf->Get("ptave_const_mpf");
  TH1D* rrlogptmpf = (TH1D*)rrmpf->Get("ptave_logpt_mpf");

  TFile* rrdijet;// = new TFile(path+"Histo_ptave_DiJet_L1.root","READ"); 
  if(syst==5){
    rrdijet = new TFile(path+"Histo_ptave_DiJet_L1"+tag+".root","READ");
  }
  if(syst==0){
    rrdijet = new TFile(path+"Histo_ptave_DiJet_L1_"+variation+".root","READ");
  }
  if(syst==1){
    rrdijet = new TFile(path+"Histo_ptave_DiJet_L1_"+variation+".root","READ");
  }
  if(syst==2){
    rrdijet = new TFile(path+"Histo_ptave_DiJet_L1_"+variation+".root","READ");
  }
  if(syst==3){
    rrdijet = new TFile(path+"Histo_ptave_DiJet_L1_"+variation+".root","READ");
  }
  TH1D* rrconstdijet = (TH1D*)rrdijet->Get("ptave_const_dijet");
  TH1D* rrlogptdijet = (TH1D*)rrdijet->Get("ptave_logpt_dijet");
  // rrconstmpf->Print();
  // rrconstdijet->Print();

  // get the kFSR root file for MPF and pt balance
  TFile* kfsrmpf = new TFile(path+"Histo_KFSR_MPF_L1.root","READ");
  TH1D* kfsrconstmpf = (TH1D*)kfsrmpf->Get("kfsr_mpf");
  TFile* kfsrdijet = new TFile(path+"Histo_KFSR_DiJet_L1.root","READ"); 
  TH1D* kfsrconstdijet = (TH1D*)kfsrdijet->Get("kfsr_dijet");


  // create histo for the normalization
  TH1D* consmpf_norm=(TH1D*)consmpf->Clone();
  TH1D* consdijet_norm=(TH1D*)consdijet->Clone();
  TH1D* rrconsmpf_norm=(TH1D*)rrconstmpf->Clone();
  TH1D* rrconsdijet_norm=(TH1D*)rrconstdijet->Clone();
  TH1D* logptmpf_norm=(TH1D*)logptmpf->Clone();
  TH1D* logptdijet_norm=(TH1D*)logptdijet->Clone();
  TH1D* rrlogptmpf_norm=(TH1D*)rrlogptmpf->Clone();
  TH1D* rrlogptdijet_norm=(TH1D*)rrlogptdijet->Clone();
  TH1D* kfsrconsmpf_norm=(TH1D*)kfsrconstmpf->Clone();
  TH1D* kfsrconsdijet_norm=(TH1D*)kfsrconstdijet->Clone();



  // create norm factor
  double norm_factor_MPF=0;
  double norm_factor_DiJet=0;
  double norm_factor_rrMPF=0;
  double norm_factor_rrDiJet=0;
  double norm_factor_MPF_logpt=0;
  double norm_factor_DiJet_logpt=0;
  double norm_factor_rrMPF_logpt=0;
  double norm_factor_rrDiJet_logpt=0;
  double norm_factor_kfsrMPF=0;
  double norm_factor_kfsrDiJet=0;

  for (int i=1; i<n_etabarr+1; i++){
    norm_factor_MPF+=consmpf->GetBinContent(i);
    norm_factor_DiJet+=consdijet->GetBinContent(i);
    norm_factor_rrMPF+=rrconstmpf->GetBinContent(i);
    norm_factor_rrDiJet+=rrconstdijet->GetBinContent(i);
    norm_factor_MPF_logpt+=logptmpf->GetBinContent(i);
    //   cout<<" "<<logptmpf->GetBinContent(i)<<" "<<norm_factor_MPF_logpt<<endl;
    norm_factor_DiJet_logpt+=logptdijet->GetBinContent(i);
    norm_factor_rrMPF_logpt+=rrlogptmpf->GetBinContent(i);
    norm_factor_rrDiJet_logpt+=rrlogptdijet->GetBinContent(i);
    norm_factor_kfsrMPF+=kfsrconstmpf->GetBinContent(i);
    norm_factor_kfsrDiJet+=kfsrconstdijet->GetBinContent(i);
  }
  norm_factor_MPF = norm_factor_MPF/n_etabarr;
  norm_factor_DiJet = norm_factor_DiJet/n_etabarr;
  norm_factor_rrMPF = norm_factor_rrMPF/n_etabarr;
  norm_factor_rrDiJet = norm_factor_rrDiJet/n_etabarr;
  norm_factor_MPF_logpt = norm_factor_MPF_logpt/n_etabarr;
  norm_factor_DiJet_logpt = norm_factor_DiJet_logpt/n_etabarr;
  norm_factor_rrMPF_logpt = norm_factor_rrMPF_logpt/n_etabarr;
  norm_factor_rrDiJet_logpt = norm_factor_rrDiJet_logpt/n_etabarr;
  norm_factor_kfsrMPF = norm_factor_kfsrMPF/n_etabarr;
  norm_factor_kfsrDiJet = norm_factor_kfsrDiJet/n_etabarr;

  // scale with norm factor
  consmpf_norm->Scale(1/norm_factor_MPF);
  consdijet_norm->Scale(1/norm_factor_DiJet);
  rrconsmpf_norm->Scale(1/norm_factor_rrMPF);
  rrconsdijet_norm->Scale(1/norm_factor_rrDiJet);
  logptmpf_norm->Scale(1/norm_factor_MPF_logpt);
  logptdijet_norm->Scale(1/norm_factor_DiJet_logpt);
  rrlogptmpf_norm->Scale(1/norm_factor_rrMPF_logpt);
  rrlogptdijet_norm->Scale(1/norm_factor_rrDiJet_logpt);
  kfsrconsmpf_norm->Scale(1/norm_factor_kfsrMPF);
  kfsrconsdijet_norm->Scale(1/norm_factor_kfsrDiJet);

  TFile* outputfilenorm = new TFile(path+"Histos_Res_norm_L1_"+jettag+"_"+txttag+"_"+variation+".root","RECREATE");
  consmpf_norm->Write();
  consdijet_norm->Write();
  rrconsmpf_norm->Write();
  rrconsdijet_norm->Write();
  logptmpf_norm->Write();
  logptdijet_norm->Write();
  rrlogptmpf_norm->Write();
  rrlogptdijet_norm->Write();
  kfsrconsmpf_norm->Write();
  kfsrconsdijet_norm->Write();
  outputfilenorm->Write();
  outputfilenorm->Close();


  cout << "normalization factor (MPF Flat): " << norm_factor_MPF << endl;
  cout << "normalization factor (ptb Flat): " << norm_factor_DiJet << endl;
  cout << "normalization factor (MPF loglin): " << norm_factor_MPF_logpt << endl;
  cout << "normalization factor (ptb loglin): " << norm_factor_DiJet_logpt << endl;

  // cout << endl;
  // for (int i=1; i<n_etabarr+1; i++){
  //   cout << consmpf_norm->GetBinContent(i) << endl;
  // }


  
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
  consmpf_norm->SetLineStyle(2);
  consdijet_norm->SetLineStyle(2);
  rrconsmpf_norm->SetLineStyle(2);
  rrconsdijet_norm->SetLineStyle(2);

  logptmpf_norm->SetLineWidth(2);
  logptmpf_norm->SetLineColor(kRed+1);
  logptdijet_norm->SetLineWidth(2);
  logptdijet_norm->SetLineColor(kBlue+1);
  rrlogptmpf_norm->SetLineWidth(2);
  rrlogptmpf_norm->SetLineColor(kRed+1);
  rrlogptdijet_norm->SetLineWidth(2);
  rrlogptdijet_norm->SetLineColor(kBlue+1);

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
  //  lumi_13TeV = "2.11 fb^{-1}";
  lumi_13TeV = "589.3 pb^{-1}";
  bool kSquare = true;

  TCanvas *c1 = tdrCanvas("c1",h,4,10,kSquare);
  //  consmpf_norm->GetYaxis()->SetRangeUser(0.91,1.15);
  consmpf_norm->GetYaxis()->SetRangeUser(0.81,1.15);
  consmpf_norm->Draw("E1 SAME");
  consdijet_norm->Draw("E1 SAME");
  logptmpf_norm->Draw("E1 SAME");
  logptdijet_norm->Draw("E1 SAME");
  line->SetLineStyle(2);
   line->Draw("SAME");

  // //  TLegend *leg1 = tdrLeg(0.17,0.49,0.40,0.80);
   TLegend *leg1 = tdrLeg(0.17,0.19,0.40,0.40);
  leg1 -> AddEntry(consmpf_norm, "MPF Flat","L");
  leg1 -> AddEntry(consdijet_norm, "Pt Flat","L");
  leg1 -> AddEntry(logptmpf_norm, "MPF Loglin","L");
  leg1 -> AddEntry(logptdijet_norm, "Pt Loglin","L");
  leg1->Draw();

   TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  
  tex->DrawLatex(0.47,0.87,JetDescrib);


  c1->SaveAs(path+"plots/L2Res_"+jettag+"_"+txttag+"_"+variation+".pdf");
  



  TCanvas *c2 = tdrCanvas("c2",h,4,10,kSquare);
  //  rrconsmpf_norm->GetYaxis()->SetTitle("(R^{MC}/R^{data})_{#alpha<0.3}");
  TString alVal;
  alVal.Form("%0.2f\n",al_cut);
  TString altitle = "{#alpha<"+alVal+"}";
  TString axistitle = "(R^{MC}/R^{data})_";
  axistitle +=altitle;
  rrconsmpf_norm->GetYaxis()->SetTitle(axistitle);
  rrconsmpf_norm->GetYaxis()->SetTitleSize(0.05);
  //  rrconsmpf_norm->GetYaxis()->SetRangeUser(0.91,1.15);
  rrconsmpf_norm->GetYaxis()->SetRangeUser(0.81,1.15);
  rrconsmpf_norm->GetXaxis()->SetTitle("|#eta|");
  rrconsmpf_norm->GetXaxis()->SetTitleSize(0.05);
  rrconsmpf_norm->Draw("E1 SAME");
  rrconsdijet_norm->Draw("E1 SAME");
  rrlogptmpf_norm->Draw("E1 SAME");
  rrlogptdijet_norm->Draw("E1 SAME");
  line->Draw("SAME");
  leg1->Draw();
  tex->DrawLatex(0.47,0.87,JetDescrib);
  c2->SaveAs(path+"plots/Ratio_"+jettag+"_"+txttag+"_"+variation+".pdf");

  TCanvas *c3 = tdrCanvas("c3",h,4,10,kSquare);
  kfsrconsmpf_norm->GetYaxis()->SetTitle("k_{FSR}");
  kfsrconsmpf_norm->GetYaxis()->SetTitleSize(0.05);
  //  kfsrconsmpf_norm->GetYaxis()->SetRangeUser(0.998,1.002);
  //  kfsrconsmpf_norm->GetYaxis()->SetRangeUser(0.91,1.15);
  kfsrconsmpf_norm->GetYaxis()->SetRangeUser(0.81,1.15);
  kfsrconsmpf_norm->GetXaxis()->SetTitle("|#eta|");
  kfsrconsmpf_norm->GetXaxis()->SetTitleSize(0.05);
 
  //
  TF1 *kfsr_fit_dijet = new TF1("kfsr_fit_dijet","[0]+([1]*TMath::CosH(x))/(1+[2]*TMath::CosH(x))",0,5.); 
  //  TF1 *kfsr_fit_dijet = new TF1("kfsr_fit_dijet","[0]+([1]*TMath::CosH(x))/(1+[2]*TMath::CosH(x))",0.,3.9); 
  //  kfsr_fit_dijet->SetParameters(0,0,0);
  //  kfsr_fit_dijet->SetParameters(0.,0,100.);//AK4CHS, phythia, Herwigg
  //  kfsr_fit_dijet->SetParameters(2.,-150.,145.); //AK4PUPPI, phythia
  //  kfsr_fit_dijet->SetParameters(1.,-1400.,100.);//AK4PUPPI, Herwigg; AK4CHS, Herwigg
  //  kfsr_fit_dijet->SetParameters(1.,-100.,100.); //AK8CHS
    kfsr_fit_dijet->SetParameters(1.,-200.,100.); //AK4CHS
  //  kfsr_fit_dijet->SetParLimits(0,1.,2.);
  //  kfsr_fit_dijet->SetParLimits(2,-100,1.);
  kfsr_fit_dijet->SetLineColor(kBlue+1);

  //  kfsrconsdijet_norm->Fit("kfsr_fit_dijet","SRMVEP");
  kfsrconsdijet_norm->Fit("kfsr_fit_dijet","SRMV");
  //Create a histogram to hold the confidence intervals
  TH1D *hint_dijet = (TH1D*)kfsrconsdijet_norm->Clone();
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_dijet);
   //Now the "hint" histogram has the fitted function values as the
   //bin contents and the confidence intervals as bin errors
   hint_dijet->SetStats(kFALSE);
   hint_dijet->SetFillColor(kBlue-10);
   
   TF1 *kfsr_fit_mpf = new TF1("kfsr_fit_mpf","[0]+([1]*TMath::CosH(x))/(1+[2]*TMath::CosH(x))",0,5.); 
   //   TF1 *kfsr_fit_mpf = new TF1("kfsr_fit_mpf","[0]+([1]*TMath::CosH(x))/(1+[2]*TMath::CosH(x))",0,3.1); 
   //   kfsr_fit_mpf->SetParameters(1,1,1);
   //   kfsr_fit_mpf->SetParameters(1.,0.,0.);
   kfsr_fit_mpf->SetParameters(0,0,100.);
   //   kfsr_fit_mpf->SetParameters(5,-100,100.);
   kfsr_fit_mpf->SetLineColor(kRed+1);
   kfsrconsmpf_norm->Fit("kfsr_fit_mpf","SR");
   //Create a histogram to hold the confidence intervals
  TH1D *hint_mpf = (TH1D*)kfsrconsmpf_norm->Clone();
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_mpf);
   //Now the "hint" histogram has the fitted function values as the
   //bin contents and the confidence intervals as bin errors
   hint_mpf->SetStats(kFALSE);
   hint_mpf->SetFillColor(kRed-10);

   hint_dijet->GetYaxis()->SetTitle("k_{FSR}");
   hint_dijet->GetYaxis()->SetTitleSize(0.05);
   //   hint_dijet->GetYaxis()->SetRangeUser(0.91,1.15);
   hint_dijet->GetYaxis()->SetRangeUser(0.81,1.15);
   hint_dijet->GetXaxis()->SetTitle("|#eta|");
   hint_dijet->GetXaxis()->SetTitleSize(0.05);

  hint_dijet->Draw("e3");
  hint_mpf->Draw("e2 SAME");
  kfsrconsmpf_norm->Draw("E1 SAME");
  kfsrconsdijet_norm->Draw("E1 SAME");


  line->Draw("SAME");

  //  TLegend *leg2 = tdrLeg(0.17,0.49,0.40,0.80);
  TLegend *leg2 = tdrLeg(0.17,0.19,0.40,0.30);
  leg2 -> AddEntry(kfsrconsmpf_norm, "MPF","L");
  leg2 -> AddEntry(kfsrconsdijet_norm, "Pt","L");
  leg2->Draw();

  tex->DrawLatex(0.47,0.87,JetDescrib);
  c3->SaveAs(path+"plots/kFSR_"+jettag+"_"+txttag+"_"+variation+".pdf");


  // use kFSR fit function instead of values in bins
  TH1D *h4 = new TH1D("h",";|#eta|;Relative correction",41,0,5.191);
  h4->SetMaximum(1.2);
  h4->SetMinimum(0.8);
  h4->GetXaxis()->SetTitleSize(0.05);
  h4->GetYaxis()->SetTitleSize(0.05);
  TCanvas *c4 = tdrCanvas("c4",h,4,10,kSquare);
  TH1D* rrconsmpf_norm_kFSR = (TH1D*)rrconsmpf_norm->Clone();
  TH1D* rrconsdijet_norm_kFSR = (TH1D*)rrconsdijet_norm->Clone();
  TH1D* rrlogptmpf_norm_kFSR = (TH1D*)rrlogptmpf_norm->Clone();
  TH1D* rrlogptdijet_norm_kFSR = (TH1D*)rrlogptdijet_norm->Clone();

  rrconsmpf_norm_kFSR->Multiply(hint_mpf);
  rrlogptmpf_norm_kFSR->Multiply(hint_mpf);
  rrconsdijet_norm_kFSR->Multiply(hint_dijet);
  rrlogptdijet_norm_kFSR->Multiply(hint_dijet);

  rrconsmpf_norm_kFSR->GetYaxis()->SetTitle("Relative correction");
  rrconsmpf_norm_kFSR->GetYaxis()->SetTitleSize(0.05);
  //  rrconsmpf_norm_kFSR->GetYaxis()->SetRangeUser(0.91,1.15);
  rrconsmpf_norm_kFSR->GetYaxis()->SetRangeUser(0.81,1.15);
  rrconsmpf_norm_kFSR->GetXaxis()->SetTitle("|#eta|");
  rrconsmpf_norm_kFSR->GetXaxis()->SetTitleSize(0.05);
  rrconsmpf_norm_kFSR->Draw("E1 SAME");
  rrconsdijet_norm_kFSR->Draw("E1 SAME");
  rrlogptmpf_norm_kFSR->Draw("E1 SAME");
  rrlogptdijet_norm_kFSR->Draw("E1 SAME");
  leg1->Draw();
  tex->DrawLatex(0.47,0.87,JetDescrib);
  line->Draw();
  // TLegend *leg1 = tdrLeg(0.17,0.49,0.40,0.80);
  // leg1 -> AddEntry(consmpf_norm, "MPF Flat","L");
  // leg1 -> AddEntry(consdijet_norm, "Pt Flat","L");
  // leg1 -> AddEntry(logptmpf_norm, "MPF Loglin","L");
  // leg1 -> AddEntry(logptdijet_norm, "Pt Loglin","L");
  // leg1->Draw();

  // TLatex *tex = new TLatex();
  // tex->SetNDC();
  // tex->SetTextSize(0.045);
  
  
  
  c4->SaveAs(path+"plots/L2Res_kFSRfit_"+jettag+"_"+txttag+"_"+variation+".pdf");

  TFile* mpf_out;// = path+"/Histo_Res_kFSRfit_MPF_L1_";  output_mpf += variation;  output_mpf += ".root";   TFile* mpf_out = new TFile(output_mpf,"RECREATE");
  if(syst==5){
    //    mpf_out = new TFile(path+"Histo_Res_kFSRfit_MPF_L1"+tag+".root","RECREATE");
    mpf_out = new TFile(path+"Histo_Res_kFSRfit_MPF_L1_"+jettag+"_"+txttag+"_"+variation+".root","RECREATE");
  }
  if(syst==0){
    mpf_out = new TFile(path+"Histo_Res_kFSRfit_MPF_L1_"+jettag+"_"+txttag+"_"+variation+".root","RECREATE");
  }
  if(syst==1){
    mpf_out = new TFile(path+"Histo_Res_kFSRfit_MPF_L1_"+jettag+"_"+txttag+"_"+variation+".root","RECREATE");
  }
  if(syst==2){
    mpf_out = new TFile(path+"Histo_Res_kFSRfit_MPF_L1_"+jettag+"_"+txttag+"_"+variation+".root","RECREATE");
  }
  if(syst==3){
    mpf_out = new TFile(path+"Histo_Res_kFSRfit_MPF_L1_"+jettag+"_"+txttag+"_"+variation+".root","RECREATE");
  }

  rrconsmpf_norm_kFSR->Write();
  rrlogptmpf_norm_kFSR->Write();
  mpf_out->Write();
  mpf_out->Close();

  // TString output_dijet = path+"/Histo_Res_kFSRfit_DiJet_L1_";
  // output_dijet += variation;
  // output_dijet += ".root";
  TFile* dijet_out;// = new TFile(output_dijet,"RECREATE");
  if(syst==5){
    //    dijet_out = new TFile(path+"Histo_Res_kFSRfit_DiJet_L1"+tag+".root","RECREATE");
   dijet_out = new TFile(path+"Histo_Res_kFSRfit_DiJet_L1_"+jettag+"_"+txttag+"_"+variation+".root","RECREATE");
  }
  if(syst==0){
    dijet_out = new TFile(path+"Histo_Res_kFSRfit_DiJet_L1_"+jettag+"_"+txttag+"_"+variation+".root","RECREATE");
  }
  if(syst==1){
   dijet_out = new TFile(path+"Histo_Res_kFSRfit_DiJet_L1_"+jettag+"_"+txttag+"_"+variation+".root","RECREATE");
  }
  if(syst==2){
    dijet_out = new TFile(path+"Histo_Res_kFSRfit_DiJet_L1_"+jettag+"_"+txttag+"_"+variation+".root","RECREATE");
  }
  if(syst==3){
    dijet_out = new TFile(path+"Histo_Res_kFSRfit_DiJet_L1_"+jettag+"_"+txttag+"_"+variation+".root","RECREATE");
  }
  rrconsdijet_norm_kFSR->Write();
  rrlogptdijet_norm_kFSR->Write();
  dijet_out->Write();
  dijet_out->Close();

  // TFile* resmpf = new TFile(path+"Histo_Res_MPF_L1.root","RECREATE");

  // TH1D* consmpf = (TH1D*)resmpf->Get("res_const_mpf");
  // TH1D* logptmpf = (TH1D*)resmpf->Get("res_logpt_mpf");
  // TFile* resdijet = new TFile(path+"Histo_Res_DiJet_L1.root","READ"); 
  // TH1D* consdijet = (TH1D*)resdijet->Get("res_const_dijet");
  // TH1D* logptdijet = (TH1D*)resdijet->Get("res_logpt_dijet");


}
