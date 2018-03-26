#include<TStyle.h>
#include<TLine.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include "../include/tdrstyle_mod15.h"
#include <TFile.h>
#include <iostream>
#include <assert.h>
#include "../include/CorrectionObject.h"
#include <TLatex.h>
#include <TGraphAsymmErrors.h>
#include "../include/parameters.h"

using namespace std;

void CorrectionObject::L2ResOutput_CorrectExtrapolation(){
  cout << "--------------- Starting L2ResOutput_CorrectExtrapolation() ---------------" << endl << endl;
   cout<<"Before files"<<endl;
  TFile* f_Res_mpf = new TFile(CorrectionObject::_outpath+"Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","READ");
  TFile* f_Res_dijet = new TFile(CorrectionObject::_outpath+"Histo_Res_Rel_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","READ");  

  TFile* f_Res_mpf_old   = new TFile("/nfs/dust/cms/user/multh/JEC/2016ReReco/Residuals/Summer16_03Feb2017_V3/AK4CHS/MC_Reweighted_chsMET_NewSF/Run"+CorrectionObject::_runnr+"/Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","READ");
  TFile* f_Res_dijet_old = new TFile("/nfs/dust/cms/user/multh/JEC/2016ReReco/Residuals/Summer16_03Feb2017_V3/AK4CHS/MC_Reweighted_chsMET_NewSF/Run"+CorrectionObject::_runnr+"/Histo_Res_DiJet_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","READ");  
   cout<<"After files"<<endl;

  TString JetDescrib;                                                                                                                            
  if (CorrectionObject::_collection=="AK4CHS") JetDescrib = "Anti-k_{t} R = 0.4, PF+CHS";
  if (CorrectionObject::_collection=="AK4Puppi") JetDescrib = "Anti-k_{t} R = 0.4, PF+PUPPI";
  if (CorrectionObject::_collection=="AK8CHS") JetDescrib = "Anti-k_{t} R = 0.8, PF+CHS";
  if (CorrectionObject::_collection=="AK8Puppi") JetDescrib = "Anti-k_{t} R = 0.8, PF+PUPPI"; 


  //get L2Res hists for MPF and pt balance
  TH1D* res_const_mpf_extrapol = (TH1D*)f_Res_mpf->Get("res_const_mpf");
  TH1D* res_const_dijet_extrapol = (TH1D*)f_Res_dijet->Get("res_const_rel");
  TH1D* res_logpt_mpf_extrapol = (TH1D*)f_Res_mpf->Get("res_logpt_mpf");
  TH1D* res_logpt_dijet_extrapol = (TH1D*)f_Res_dijet->Get("res_logpt_rel");


  //get L2Res hists for MPF and pt balance with Old kFSR Fit values
  TH1D* res_const_mpf_kfsrfit = (TH1D*)f_Res_mpf_old->Get("res_const_mpf");
  TH1D* res_const_dijet_kfsrfit = (TH1D*)f_Res_dijet_old->Get("res_const_dijet");
  TH1D* res_logpt_mpf_kfsrfit = (TH1D*)f_Res_mpf_old->Get("res_logpt_mpf");
  TH1D* res_logpt_dijet_kfsrfit = (TH1D*)f_Res_dijet_old->Get("res_logpt_dijet");

  //get L2Res hists for MPF and pt balance with Old kFSR Histogram values
  TH1D* res_const_mpf_kfsrval = (TH1D*)f_Res_mpf_old->Get("res_const_mpf_val");
  TH1D* res_const_dijet_kfsrval = (TH1D*)f_Res_dijet_old->Get("res_const_dijet_val");
  TH1D* res_logpt_mpf_kfsrval = (TH1D*)f_Res_mpf_old->Get("res_logpt_mpf_val");
  TH1D* res_logpt_dijet_kfsrval = (TH1D*)f_Res_dijet_old->Get("res_logpt_dijet_val");

  //Draw Difference between new Method and kFSR Fit results
  TH1D* res_const_mpf_fit_diff = (TH1D*) res_const_mpf_kfsrfit->Clone("res_const_mpf_fit_diff");
  res_const_mpf_fit_diff->Add(res_const_mpf_kfsrfit,-1);
  TH1D* res_const_mpf_val_diff = (TH1D*) res_const_mpf_kfsrfit->Clone("res_const_mpf_val_diff");
  res_const_mpf_val_diff->Add(res_const_mpf_kfsrval,-1);
  TH1D* res_logpt_mpf_fit_diff = (TH1D*) res_logpt_mpf_kfsrfit->Clone("res_logpt_mpf_fit_diff");
  res_logpt_mpf_fit_diff->Add(res_logpt_mpf_kfsrfit,-1);
  TH1D* res_logpt_mpf_val_diff = (TH1D*) res_logpt_mpf_kfsrfit->Clone("res_logpt_mpf_val_diff");
  res_logpt_mpf_val_diff->Add(res_logpt_mpf_kfsrval,-1);

  //Draw Difference between new Method and kFSR HistVal results
  TH1D* res_const_dijet_fit_diff = (TH1D*) res_const_dijet_kfsrfit->Clone("res_const_dijet_fit_diff");
  res_const_dijet_fit_diff->Add(res_const_dijet_kfsrfit,-1);
  TH1D* res_const_dijet_val_diff = (TH1D*) res_const_dijet_kfsrfit->Clone("res_const_dijet_val_diff");
  res_const_dijet_val_diff->Add(res_const_dijet_kfsrval,-1);
  TH1D* res_logpt_dijet_fit_diff = (TH1D*) res_logpt_dijet_kfsrfit->Clone("res_logpt_dijet_fit_diff");
  res_logpt_dijet_fit_diff->Add(res_logpt_dijet_kfsrfit,-1);
  TH1D* res_logpt_dijet_val_diff = (TH1D*) res_logpt_dijet_kfsrfit->Clone("res_logpt_dijet_val_diff");
  res_logpt_dijet_val_diff->Add(res_logpt_dijet_kfsrval,-1);


 
  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0,1,5.191,1);
  TLine *line2 = new TLine(0,0,5.191,0);
  res_const_mpf_extrapol   ->SetLineWidth(2);
  res_const_mpf_extrapol   ->SetLineColor(kRed+1);
  res_const_dijet_extrapol ->SetLineWidth(2);
  res_const_dijet_extrapol ->SetLineColor(kBlue+1);
  res_const_mpf_extrapol   ->SetLineStyle(2);
  res_const_dijet_extrapol ->SetLineStyle(2);
  res_logpt_mpf_extrapol   ->SetLineWidth(2);
  res_logpt_mpf_extrapol   ->SetLineColor(kRed+1);
  res_logpt_dijet_extrapol ->SetLineWidth(2);
  res_logpt_dijet_extrapol ->SetLineColor(kBlue+1);

  //**************************Diff between new Method and kFSR with FIT Value **********************************
  res_const_mpf_fit_diff   ->SetLineWidth(2);
  res_const_mpf_fit_diff   ->SetLineColor(kRed+1);
  res_const_dijet_fit_diff ->SetLineWidth(2);
  res_const_dijet_fit_diff ->SetLineColor(kBlue+1);
  res_const_mpf_fit_diff   ->SetLineStyle(2);
  res_const_dijet_fit_diff ->SetLineStyle(2);
  res_logpt_mpf_fit_diff   ->SetLineWidth(2);
  res_logpt_mpf_fit_diff   ->SetLineColor(kRed+1);
  res_logpt_dijet_fit_diff ->SetLineWidth(2);
  res_logpt_dijet_fit_diff ->SetLineColor(kBlue+1);
  //*********************************************************************************
  
  //**************************Diff between new Method and kFSR with HIST Value **********************************
  res_const_mpf_val_diff   ->SetLineWidth(2);
  res_const_mpf_val_diff   ->SetLineColor(kRed+1);
  res_const_dijet_val_diff ->SetLineWidth(2);
  res_const_dijet_val_diff ->SetLineColor(kBlue+1);
  res_const_mpf_val_diff   ->SetLineStyle(2);
  res_const_dijet_val_diff ->SetLineStyle(2);
  res_logpt_mpf_val_diff   ->SetLineWidth(2);
  res_logpt_mpf_val_diff   ->SetLineColor(kRed+1);
  res_logpt_dijet_val_diff ->SetLineWidth(2);
  res_logpt_dijet_val_diff ->SetLineColor(kBlue+1);
  //*********************************************************************************

  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);  
  tex->DrawLatex(0.45,0.87,JetDescrib);

  TLatex *tex_bars = new TLatex();
  tex_bars->SetNDC();
  tex_bars->SetTextSize(0.039);
    
  TLatex *tex1 = new TLatex();
  tex1->SetNDC();
  tex1->SetTextSize(0.045); 
  tex1->DrawLatex(0.5,0.91,CorrectionObject::_lumitag+"(13TeV)"); 
  
  TCanvas *c1 = new TCanvas();
  c1 = new TCanvas("c1","L2res_kFSRextrapol",600,600);
  c1->DrawFrame(0,0.78,5.191,1.23,(";|#eta|;Relative corrections"));

  res_const_mpf_extrapol->SetMarkerStyle(1);
  res_const_dijet_extrapol->SetMarkerStyle(1);
  res_logpt_mpf_extrapol->SetMarkerStyle(1);
  res_logpt_dijet_extrapol->SetMarkerStyle(1);
  res_const_mpf_extrapol->Draw("E1 SAME");
  res_const_dijet_extrapol->Draw("E1 SAME");
  res_logpt_mpf_extrapol->Draw("E1 SAME");
  res_logpt_dijet_extrapol->Draw("E1 SAME");

  TLegend* leg1 = new TLegend(0.20,0.19,0.43,0.40,"","brNDC");
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetFillColor(10);
  leg1->SetLineColor(1);
  leg1->SetTextFont(42);
  
  leg1->AddEntry(res_const_mpf_extrapol, "MPF Flat","L");
  leg1->AddEntry(res_const_dijet_extrapol, "Pt Flat","L");
  leg1->AddEntry(res_logpt_mpf_extrapol, "MPF Loglin","L");
  leg1->AddEntry(res_logpt_dijet_extrapol, "Pt Loglin","L"); 
  

  leg1->Draw();
  tex ->DrawLatex(0.45,0.87,JetDescrib);
  line->Draw();
  tex1->DrawLatex(0.5,0.91,CorrectionObject::_lumitag+"(13TeV)"); 
  c1->SaveAs(CorrectionObject::_outpath+"plots/L2Res_kFSRextrapol_"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");



  TCanvas *c2 = new TCanvas();
  c2 = new TCanvas("c2","L2res_Diff_extrapol_fit",600,600);
  c2->DrawFrame(0,0.78,5.191,1.23,(";|#eta|;New - kFSR (Fit)"));
  res_const_mpf_fit_diff   ->SetMarkerStyle(1);
  res_const_dijet_fit_diff ->SetMarkerStyle(1);
  res_logpt_mpf_fit_diff   ->SetMarkerStyle(1);
  res_logpt_dijet_fit_diff ->SetMarkerStyle(1);

  res_const_mpf_fit_diff->Draw("E1 SAME");
  res_const_dijet_fit_diff->Draw("E1 SAME");
  res_logpt_mpf_fit_diff->Draw("E1 SAME");
  res_logpt_dijet_fit_diff->Draw("E1 SAME");

  leg1->Draw();
  tex ->DrawLatex(0.45,0.87,JetDescrib);
  line->Draw();
  tex1->DrawLatex(0.5,0.91,CorrectionObject::_lumitag+"(13TeV)"); 
  c2->SaveAs(CorrectionObject::_outpath+"plots/L2Res_diff_FitToExtrapol_"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");



  TCanvas *c3 = new TCanvas();
  c3 = new TCanvas("c3","L2res_Diff_extrapol_val",600,600);
  c3->DrawFrame(0,0.78,5.191,1.23,(";|#eta|;New - kFSR (Hist)"));
  res_const_mpf_val_diff   ->SetMarkerStyle(1);
  res_const_dijet_val_diff ->SetMarkerStyle(1);
  res_logpt_mpf_val_diff   ->SetMarkerStyle(1);
  res_logpt_dijet_val_diff ->SetMarkerStyle(1);

  res_const_mpf_val_diff->Draw("E1 SAME");
  res_const_dijet_val_diff->Draw("E1 SAME");
  res_logpt_mpf_val_diff->Draw("E1 SAME");
  res_logpt_dijet_val_diff->Draw("E1 SAME");

  leg1->Draw();
  tex ->DrawLatex(0.45,0.87,JetDescrib);
  line->Draw();
  tex1->DrawLatex(0.5,0.91,CorrectionObject::_lumitag+"(13TeV)"); 
  c3->SaveAs(CorrectionObject::_outpath+"plots/L2Res_diff_ValToExtrapol_"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");



  //pt-dependence of L2Res (Money plot)

  TString var="";
  TH1D* res_logpt_mpf_extrapol_pTdep[4];
  TH1D* res_logpt_dijet_extrapol_pTdep[4];
  
  TFile* f_Res_mpf_pTdep;
  TFile* f_Res_dijet_pTdep;

  for(int i=0;i<4;i++){
    if(i==0) var="central"; 
    if(i==1) var="down";
    if(i==2) var="up";
    if(i==3) var="doubleup";

    f_Res_mpf_pTdep = new TFile(CorrectionObject::_outpath+"Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_"+var+".root","READ");
    f_Res_dijet_pTdep = new TFile(CorrectionObject::_outpath+"Histo_Res_Rel_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_"+var+".root","READ");

    res_logpt_mpf_extrapol_pTdep[i]   = (TH1D*)f_Res_mpf_pTdep   ->Get("res_logpt_mpf");
    res_logpt_dijet_extrapol_pTdep[i] = (TH1D*)f_Res_dijet_pTdep ->Get("res_logpt_rel");
  }


  TCanvas *c4 = new TCanvas();
  c4 = new TCanvas("c4","L2res_logpt_MPF_extrapol_ptDepend",600,600);
  res_logpt_mpf_extrapol->SetLineColor(kBlack);
  res_logpt_mpf_extrapol->Draw("E1 SAME");
  for(int i=0;i<4;i++){
    res_logpt_mpf_extrapol_pTdep[i]->SetLineColor(kRed-3*i);
    res_logpt_mpf_extrapol_pTdep[i]->SetMarkerColor(kRed-3*i);
    res_logpt_mpf_extrapol_pTdep[i]->SetMarkerStyle(20+i);
    res_logpt_mpf_extrapol_pTdep[i]->Draw("E1 SAME"); 
  }
  
  TLegend* leg2 =  new TLegend(0.20,0.19,0.43,0.40,"","brNDC"); 
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.03);
  leg2->SetFillColor(10);
  leg2->SetLineColor(1);
  leg2->SetTextFont(42);
  leg2-> AddEntry(res_logpt_mpf_extrapol_pTdep[1] , "60 GeV","LP");  
  leg2-> AddEntry(res_logpt_mpf_extrapol_pTdep[0] , "120 GeV","LP");
  leg2-> AddEntry(res_logpt_mpf_extrapol_pTdep[2] , "240 GeV","LP"); 
  leg2-> AddEntry(res_logpt_mpf_extrapol_pTdep[3] , "480 GeV","LP");
  leg2-> AddEntry(res_logpt_mpf_extrapol, "Nominal","LP");
  leg2->Draw();              
  
  tex->DrawLatex(0.45,0.87,JetDescrib);  
  tex1->DrawLatex(0.5,0.91,CorrectionObject::_lumitag+"(13TeV)"); 
  c4->SaveAs(CorrectionObject::_outpath+"plots/L2Res_logpt_MPF_extrapol_"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");

 


  TCanvas *c5 = new TCanvas();
  c5 = new TCanvas("c5","L2res_logpt_Rel_extrapol_ptDepend",600,600);
  res_logpt_dijet_extrapol->SetLineColor(kBlack);
  res_logpt_dijet_extrapol->Draw("E1 SAME");
  for(int i=0;i<4;i++){
    res_logpt_dijet_extrapol_pTdep[i]->SetLineColor(kRed-3*i);
    res_logpt_dijet_extrapol_pTdep[i]->SetMarkerColor(kRed-3*i);
    res_logpt_dijet_extrapol_pTdep[i]->SetMarkerStyle(20+i);
    res_logpt_dijet_extrapol_pTdep[i]->Draw("E1 SAME"); 
  }

  TLegend* leg3 =  new TLegend(0.20,0.19,0.43,0.40,"","brNDC"); 
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.03);
  leg3->SetFillColor(10);
  leg3->SetLineColor(1);
  leg3->SetTextFont(42);

  leg3-> AddEntry(res_logpt_dijet_extrapol_pTdep[1] , "60 GeV","LP");  
  leg3-> AddEntry(res_logpt_dijet_extrapol_pTdep[0] , "120 GeV","LP");
  leg3-> AddEntry(res_logpt_dijet_extrapol_pTdep[2] , "240 GeV","LP"); 
  leg3-> AddEntry(res_logpt_dijet_extrapol_pTdep[3] , "480 GeV","LP");
  leg3-> AddEntry(res_logpt_dijet_extrapol, "Nominal","LP");
  leg3->Draw();              
 
  tex->DrawLatex(0.45,0.87,JetDescrib);  
  tex1->DrawLatex(0.5,0.91,CorrectionObject::_lumitag+"(13TeV)"); 
  c5->SaveAs(CorrectionObject::_outpath+"plots/L2Res_logpt_DiJet_extrapol_"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");
 



 

  cout << " ------------------- Finished L2Res(). Going to delete everything. ------------------------ " << endl;

  //delete everything
  delete c1;
  delete c2;
  delete c3;
  delete c4;
  delete c5;

  for(int i=0; i<4; i++){
    delete res_logpt_dijet_extrapol_pTdep[i];
    delete res_logpt_mpf_extrapol_pTdep[i];
  }

  delete f_Res_mpf_pTdep;
  delete f_Res_dijet_pTdep;

  delete tex_bars;
  delete tex;
  delete tex1; 
  
  delete line;
  delete line2;
 
  delete res_const_mpf_val_diff;
  delete res_const_dijet_val_diff;
  delete res_logpt_mpf_val_diff;
  delete res_logpt_dijet_val_diff;

  delete res_const_mpf_fit_diff;
  delete res_const_dijet_fit_diff;
  delete res_logpt_mpf_fit_diff;
  delete res_logpt_dijet_fit_diff;

  delete res_const_mpf_extrapol;
  delete res_const_dijet_extrapol;
  delete res_logpt_mpf_extrapol;
  delete res_logpt_dijet_extrapol;

  delete f_Res_dijet;
  delete f_Res_mpf;

  delete f_Res_mpf_old;
  delete f_Res_dijet_old;
}
