#include<TStyle.h>
#include<TLine.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include "tdrstyle_mod15.C"
#include <TFile.h>
#include <iostream>
#include <assert.h>
#include "../include/CorrectionObject.h"
#include <TLatex.h>
#include <TGraphAsymmErrors.h>
#include "../include/parameters.h"

using namespace std;

void CorrectionObject::L2ResOutput(){
  cout << "--------------- Starting L2ResOutput() ---------------" << endl << endl;

  TFile* f_Res_mpf = new TFile(CorrectionObject::_outpath+"Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","READ");
  TFile* f_Res_dijet = new TFile(CorrectionObject::_outpath+"Histo_Res_DiJet_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","READ");  

  TString JetDescrib;                                                                                                                            
  if (CorrectionObject::_collection=="AK4CHS") JetDescrib = "Anti-k_{t} R = 0.4, PF+CHS";
  if (CorrectionObject::_collection=="AK4Puppi") JetDescrib = "Anti-k_{t} R = 0.4, PF+PUPPI";
  if (CorrectionObject::_collection=="AK8CHS") JetDescrib = "Anti-k_{t} R = 0.8, PF+CHS";
  if (CorrectionObject::_collection=="AK8Puppi") JetDescrib = "Anti-k_{t} R = 0.8, PF+PUPPI"; 

  //plot results for "nominal" variation
 

  // get the (R_{MC}/R_{DATA}) hists for MPF and pt balance
  TH1D* pt_depend_const_mpf = (TH1D*)f_Res_mpf->Get("ptave_const_mpf");
  TH1D* pt_depend_logpt_mpf = (TH1D*)f_Res_mpf->Get("ptave_logpt_mpf");
  TH1D* pt_depend_const_dijet = (TH1D*)f_Res_dijet->Get("ptave_const_dijet");
  TH1D* pt_depend_logpt_dijet = (TH1D*)f_Res_dijet->Get("ptave_logpt_dijet");

  // get the kFSR hists for MPF and pt balance
  TH1D* kfsr_mpf = (TH1D*)f_Res_mpf->Get("kfsr_mpf");
  TH1D* kfsr_mpf_fit = (TH1D*)f_Res_mpf->Get("hist_kfsr_fit_mpf");
  TH1D* kfsr_dijet = (TH1D*)f_Res_dijet->Get("kfsr_dijet");
  TH1D* kfsr_dijet_fit = (TH1D*)f_Res_dijet->Get("hist_kfsr_fit_dijet");

  //get L2Res hists for MPF and pt balance
  TH1D* res_const_mpf_kfsrfit = (TH1D*)f_Res_mpf->Get("res_const_mpf");
  TH1D* res_const_dijet_kfsrfit = (TH1D*)f_Res_dijet->Get("res_const_dijet");
  TH1D* res_logpt_mpf_kfsrfit = (TH1D*)f_Res_mpf->Get("res_logpt_mpf");
  TH1D* res_logpt_dijet_kfsrfit = (TH1D*)f_Res_dijet->Get("res_logpt_dijet");


  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,5.191,1);
  res_const_mpf_kfsrfit->SetLineWidth(2);
  res_const_mpf_kfsrfit->SetLineColor(kRed+1);
  res_const_dijet_kfsrfit->SetLineWidth(2);
  res_const_dijet_kfsrfit->SetLineColor(kBlue+1);
  res_const_mpf_kfsrfit->SetLineStyle(2);
  res_const_dijet_kfsrfit->SetLineStyle(2);
  res_logpt_mpf_kfsrfit->SetLineWidth(2);
  res_logpt_mpf_kfsrfit->SetLineColor(kRed+1);
  res_logpt_dijet_kfsrfit->SetLineWidth(2);
  res_logpt_dijet_kfsrfit->SetLineColor(kBlue+1);

  pt_depend_const_mpf->SetLineWidth(2);
  pt_depend_const_mpf->SetLineColor(kRed+1);
  pt_depend_const_dijet->SetLineWidth(2);
  pt_depend_const_dijet->SetLineColor(kBlue+1);
  pt_depend_const_mpf->SetLineStyle(2);
  pt_depend_const_dijet->SetLineStyle(2);
  pt_depend_logpt_mpf->SetLineWidth(2);
  pt_depend_logpt_mpf->SetLineColor(kRed+1);
  pt_depend_logpt_dijet->SetLineWidth(2);
  pt_depend_logpt_dijet->SetLineColor(kBlue+1);

  kfsr_mpf->SetLineWidth(2);
  kfsr_mpf->SetLineColor(kRed+1);
  kfsr_dijet->SetLineWidth(2);
  kfsr_dijet->SetLineColor(kBlue+1);


  // Draw results
  TH1D *h = new TH1D("h",";|#eta|;Relative correction",41,0,5.191);
  h->SetMaximum(1.2);
  h->SetMinimum(0.8);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  lumi_13TeV = CorrectionObject::_lumitag;
  bool kSquare = true;

  TGraphAsymmErrors* gr_pt_depend_const_mpf = new TGraphAsymmErrors(pt_depend_const_mpf);
  TGraphAsymmErrors* gr_pt_depend_logpt_mpf = new TGraphAsymmErrors(pt_depend_logpt_mpf);
  TGraphAsymmErrors* gr_pt_depend_const_dijet = new TGraphAsymmErrors(pt_depend_const_dijet);
  TGraphAsymmErrors* gr_pt_depend_logpt_dijet = new TGraphAsymmErrors(pt_depend_logpt_dijet);
  

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);  
  tex->DrawLatex(0.45,0.87,JetDescrib);

  TLatex *tex_bars = new TLatex();
  tex_bars->SetNDC();
  tex_bars->SetTextSize(0.039);

  TLegend *leg1 = tdrLeg(0.17,0.19,0.40,0.40);
  leg1 -> AddEntry(pt_depend_const_mpf, "MPF Flat","L");
  leg1 -> AddEntry(pt_depend_const_dijet, "Pt Flat","L");
  leg1 -> AddEntry(pt_depend_logpt_mpf, "MPF Loglin","L");
  leg1 -> AddEntry(pt_depend_logpt_dijet, "Pt Loglin","L"); 

  TCanvas *c2 = tdrCanvas("c2",h,4,10,kSquare);

  TString alVal;
  alVal.Form("%0.2f\n",alpha_cut);
  TString altitle = "{#alpha<"+alVal+"}";
  TString axistitle = "(R^{MC}/R^{data})_";
  axistitle +=altitle;
  h->GetYaxis()->SetTitle(axistitle);
  h->GetYaxis()->SetRangeUser(0.81,1.15); 

  pt_depend_const_mpf->SetMarkerStyle(1);
  pt_depend_const_dijet->SetMarkerStyle(1);
  pt_depend_logpt_mpf->SetMarkerStyle(1);
  pt_depend_logpt_dijet->SetMarkerStyle(1);
  pt_depend_const_mpf->Draw("E1 SAME");
  pt_depend_const_dijet->Draw("E1 SAME");
  pt_depend_logpt_mpf->Draw("E1 SAME");
  pt_depend_logpt_dijet->Draw("E1 SAME");

  line->Draw("SAME");
  leg1->Draw();
  tex->DrawLatex(0.45,0.87,JetDescrib);
  c2->SaveAs(CorrectionObject::_outpath+"plots/Ratio_"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");


  TCanvas *c3 = tdrCanvas("c3",h,4,10,kSquare);
  h->GetYaxis()->SetTitle("k_{FSR}");
  h->GetYaxis()->SetRangeUser(0.81,1.15);
  kfsr_dijet_fit->SetMarkerStyle(1);
  kfsr_mpf_fit->SetMarkerStyle(1);
  kfsr_mpf->SetMarkerStyle(1);
  kfsr_dijet->SetMarkerStyle(1);
  kfsr_dijet_fit->Draw("E3 SAME");
  kfsr_mpf_fit->Draw("E3 SAME");
  kfsr_mpf->Draw("E1 SAME");
  kfsr_dijet->Draw("E1 SAME");
  line->Draw("SAME");

  TLegend *leg2 = tdrLeg(0.17,0.19,0.40,0.30);
  leg2 -> AddEntry(kfsr_mpf, "MPF","L");
  leg2 -> AddEntry(kfsr_dijet, "Pt","L");
  leg2->Draw();
  tex->DrawLatex(0.45,0.87,JetDescrib);
  c3->SaveAs(CorrectionObject::_outpath+"plots/kFSR_"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");


  TCanvas *c4 = tdrCanvas("L2res_kFSRfit",h,4,10,kSquare);
  h->GetYaxis()->SetTitle("Relative correction");
  h->GetYaxis()->SetRangeUser(0.8,1.2);
  res_const_mpf_kfsrfit->SetMarkerStyle(1);
  res_const_dijet_kfsrfit->SetMarkerStyle(1);
  res_logpt_mpf_kfsrfit->SetMarkerStyle(1);
  res_logpt_dijet_kfsrfit->SetMarkerStyle(1);
  res_const_mpf_kfsrfit->Draw("E1 SAME");
  res_const_dijet_kfsrfit->Draw("E1 SAME");
  res_logpt_mpf_kfsrfit->Draw("E1 SAME");
  res_logpt_dijet_kfsrfit->Draw("E1 SAME");
  leg1->Draw();
  tex->DrawLatex(0.45,0.87,JetDescrib);
  line->Draw();
  c4->SaveAs(CorrectionObject::_outpath+"plots/L2Res_kFSRfit_"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");


  //pt-dependence of L2Res (Money plot)

  TString var="";
  TH1D* res_logpt_mpf_kfsrfit_var[4];
  TH1D* res_logpt_dijet_kfsrfit_var[4];
  for(int i=0;i<4;i++){
    if(i==0) var="central"; 
    if(i==1) var="down";
    if(i==2) var="up";
    if(i==3) var="doubleup";

    TFile* f_Res_mpf_var;
    TFile* f_Res_dijet_var;
    f_Res_mpf_var = new TFile(CorrectionObject::_outpath+"Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_"+var+".root","READ");
    f_Res_dijet_var = new TFile(CorrectionObject::_outpath+"Histo_Res_DiJet_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_"+var+".root","READ"); 
    res_logpt_mpf_kfsrfit_var[i] = (TH1D*)f_Res_mpf_var->Get("res_logpt_mpf");
    res_logpt_dijet_kfsrfit_var[i] = (TH1D*)f_Res_dijet_var->Get("res_logpt_dijet");
  }

  TCanvas *c5 = tdrCanvas("L2res_logpt_MPF_kFSRfit_ptDepend",h,4,10,kSquare);
  res_logpt_mpf_kfsrfit->SetLineColor(kBlack);
  res_logpt_mpf_kfsrfit->Draw("E1 SAME");
  for(int i=0;i<4;i++){
    res_logpt_mpf_kfsrfit_var[i]->SetLineColor(kRed-3*i);
    res_logpt_mpf_kfsrfit_var[i]->SetMarkerColor(kRed-3*i);
    res_logpt_mpf_kfsrfit_var[i]->SetMarkerStyle(20+i);
    res_logpt_mpf_kfsrfit_var[i]->Draw("E1 SAME"); 
  }

  leg2 = tdrLeg(0.17,0.19,0.40,0.42); 
  leg2 -> AddEntry(res_logpt_mpf_kfsrfit_var[1] , "60 GeV","LP");  
  leg2 -> AddEntry(res_logpt_mpf_kfsrfit_var[0] , "120 GeV","LP");
  leg2 -> AddEntry(res_logpt_mpf_kfsrfit_var[2] , "240 GeV","LP"); 
  leg2 -> AddEntry(res_logpt_mpf_kfsrfit_var[3] , "480 GeV","LP");
  leg2 -> AddEntry(res_logpt_mpf_kfsrfit, "Nominal","LP");
  leg2->Draw();              
  tex->DrawLatex(0.45,0.87,JetDescrib);      
  c5->SaveAs(CorrectionObject::_outpath+"plots/L2Res_logpt_MPF_kFSRfit_"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");


  TCanvas *c6 = tdrCanvas("L2res_logpt_DiJet_kFSRfit_ptDepend",h,4,10,kSquare);
  res_logpt_dijet_kfsrfit->SetLineColor(kBlack);
  res_logpt_dijet_kfsrfit->Draw("E1 SAME");
  for(int i=0;i<4;i++){
    res_logpt_dijet_kfsrfit_var[i]->SetLineColor(kBlue+3*i);
    res_logpt_dijet_kfsrfit_var[i]->SetMarkerColor(kBlue+3*i);
    res_logpt_dijet_kfsrfit_var[i]->SetMarkerStyle(20+i);
    res_logpt_dijet_kfsrfit_var[i]->Draw("E1 SAME");  

  }
  leg2 = tdrLeg(0.17,0.19,0.40,0.42); 
  leg2 -> AddEntry(res_logpt_dijet_kfsrfit_var[1] , "60 GeV","LP");  
  leg2 -> AddEntry(res_logpt_dijet_kfsrfit_var[0] , "120 GeV","LP");
  leg2 -> AddEntry(res_logpt_dijet_kfsrfit_var[2] , "240 GeV","LP"); 
  leg2 -> AddEntry(res_logpt_dijet_kfsrfit_var[3] , "480 GeV","LP");
  leg2 -> AddEntry(res_logpt_dijet_kfsrfit, "Nominal","LP");
  leg2->Draw();
  tex->DrawLatex(0.45,0.87,JetDescrib);                    
  c6->SaveAs(CorrectionObject::_outpath+"plots/L2Res_logpt_DiJet_kFSRfit_"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");

  f_Res_mpf->Close();
  f_Res_dijet->Close();






















}
