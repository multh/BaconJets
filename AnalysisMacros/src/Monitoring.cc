#include<TStyle.h>
#include<TLine.h>
#include <TH1.h>
#include <TH2.h>
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
#include <TProfile.h>
#include "../include/useful_functions.h"

using namespace std;

void CorrectionObject::Monitoring(){
  cout << "--------------- Starting Monitoring() ---------------" << endl << endl;

  TFile* f_monitoring_BCD = new TFile(CorrectionObject::_input_path+"/uhh2.AnalysisModuleRunner.DATA.DATA_RunBCD_AK4CHS.root");
  
  TFile* f_monitoring_EFearly = new TFile(CorrectionObject::_input_path+"/uhh2.AnalysisModuleRunner.DATA.DATA_RunEFearly_AK4CHS.root");

  TFile* f_monitoring_FlateG = new TFile(CorrectionObject::_input_path+"/uhh2.AnalysisModuleRunner.DATA.DATA_RunFlateG_AK4CHS.root","READ");
  TFile* f_monitoring_H = new TFile(CorrectionObject::_input_path+"/uhh2.AnalysisModuleRunner.DATA.DATA_RunH_AK4CHS.root","READ"); 
    
  TH2D *hist_A_BCD_asymmetry[n_eta-1][n_pt-1];
  TH2D *hist_A_EFearly_asymmetry[n_eta-1][n_pt-1];
  TH2D *hist_A_FlateG_asymmetry[n_eta-1][n_pt-1];
  TH2D *hist_A_H_asymmetry[n_eta-1][n_pt-1];

   
  TProfile *pr_A_BCD_asymmetry[n_eta-1][n_pt-1];
  TProfile *pr_A_EFearly_asymmetry[n_eta-1][n_pt-1];
  TProfile *pr_A_FlateG_asymmetry[n_eta-1][n_pt-1];
  TProfile *pr_A_H_asymmetry[n_eta-1][n_pt-1];

  TH2D *hist_B_BCD_asymmetry[n_eta-1][n_pt-1];
  TH2D *hist_B_EFearly_asymmetry[n_eta-1][n_pt-1];
  TH2D *hist_B_FlateG_asymmetry[n_eta-1][n_pt-1];
  TH2D *hist_B_H_asymmetry[n_eta-1][n_pt-1];
   
  TProfile *pr_B_BCD_asymmetry[n_eta-1][n_pt-1];
  TProfile *pr_B_EFearly_asymmetry[n_eta-1][n_pt-1];
  TProfile *pr_B_FlateG_asymmetry[n_eta-1][n_pt-1];
  TProfile *pr_B_H_asymmetry[n_eta-1][n_pt-1];


    for(int j=0; j<n_eta-1; j++){
      for(int i=0; i<n_pt-1; i++){

        hist_A_BCD_asymmetry[j][i]       = (TH2D*)f_monitoring_BCD    ->Get("Monitoring_Final/hist_data_A_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_pT_"+pt_range[i]+"_"+pt_range[i+1]);
        hist_A_EFearly_asymmetry[j][i]   = (TH2D*)f_monitoring_EFearly->Get("Monitoring_Final/hist_data_A_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_pT_"+pt_range[i]+"_"+pt_range[i+1]);
        hist_A_FlateG_asymmetry[j][i]    = (TH2D*)f_monitoring_FlateG ->Get("Monitoring_Final/hist_data_A_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_pT_"+pt_range[i]+"_"+pt_range[i+1]);
	hist_A_H_asymmetry[j][i]         = (TH2D*)f_monitoring_H      ->Get("Monitoring_Final/hist_data_A_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_pT_"+pt_range[i]+"_"+pt_range[i+1]);


	pr_A_BCD_asymmetry[j][i]    =(TProfile*)hist_A_BCD_asymmetry[j][i]    ->ProfileX(Form("profABCD%d%d",j,i));
	pr_A_EFearly_asymmetry[j][i]=(TProfile*)hist_A_EFearly_asymmetry[j][i]->ProfileX(Form("profAEFearly%d%d",j,i));
	pr_A_FlateG_asymmetry[j][i] =(TProfile*)hist_A_FlateG_asymmetry[j][i] ->ProfileX(Form("profAFlateG%d%d",j,i));
	pr_A_H_asymmetry[j][i]      =(TProfile*)hist_A_H_asymmetry[j][i]      ->ProfileX(Form("profAH%d%d",j,i));
	

        hist_B_BCD_asymmetry[j][i]       = (TH2D*)f_monitoring_BCD    ->Get("Monitoring_Final/hist_data_B_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_pT_"+pt_range[i]+"_"+pt_range[i+1]);
        hist_B_EFearly_asymmetry[j][i]   = (TH2D*)f_monitoring_EFearly->Get("Monitoring_Final/hist_data_B_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_pT_"+pt_range[i]+"_"+pt_range[i+1]);
        hist_B_FlateG_asymmetry[j][i]    = (TH2D*)f_monitoring_FlateG ->Get("Monitoring_Final/hist_data_B_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_pT_"+pt_range[i]+"_"+pt_range[i+1]);
	hist_B_H_asymmetry[j][i]         = (TH2D*)f_monitoring_H      ->Get("Monitoring_Final/hist_data_B_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_pT_"+pt_range[i]+"_"+pt_range[i+1]);

	pr_B_BCD_asymmetry[j][i]    =(TProfile*)hist_B_BCD_asymmetry[j][i]    ->ProfileX(Form("profBBCD%d%d",j,i));
	pr_B_EFearly_asymmetry[j][i]=(TProfile*)hist_B_EFearly_asymmetry[j][i]->ProfileX(Form("profBEFearly%d%d",j,i));
	pr_B_FlateG_asymmetry[j][i] =(TProfile*)hist_B_FlateG_asymmetry[j][i] ->ProfileX(Form("profBFlateG%d%d",j,i));
	pr_B_H_asymmetry[j][i]      =(TProfile*)hist_B_H_asymmetry[j][i]      ->ProfileX(Form("profBH%d%d",j,i));
      }
    }

   int n_lumi = pr_B_BCD_asymmetry[1][1]->GetNbinsX();
   double bin_width = pr_B_BCD_asymmetry[1][1]->GetXaxis()->GetBinWidth(1);


  double res_mpf_BCD[n_eta-1][n_pt-1][n_lumi-1];
  double res_mpf_EFearly[n_eta-1][n_pt-1][n_lumi-1];
  double res_mpf_FlateG[n_eta-1][n_pt-1][n_lumi-1];
  double res_mpf_H[n_eta-1][n_pt-1][n_lumi-1];

  double err_res_mpf_BCD[n_eta-1][n_pt-1][n_lumi-1];
  double err_res_mpf_EFearly[n_eta-1][n_pt-1][n_lumi-1];
  double err_res_mpf_FlateG[n_eta-1][n_pt-1][n_lumi-1];
  double err_res_mpf_H[n_eta-1][n_pt-1][n_lumi-1];

  TGraphErrors *mpf_res_BCD[n_eta-1][n_pt-1]; 
  TGraphErrors *mpf_res_EFearly[n_eta-1][n_pt-1]; 
  TGraphErrors *mpf_res_FlateG[n_eta-1][n_pt-1]; 
  TGraphErrors *mpf_res_H[n_eta-1][n_pt-1]; 

  double res_rel_BCD[n_eta-1][n_pt-1][n_lumi-1];
  double res_rel_EFearly[n_eta-1][n_pt-1][n_lumi-1];
  double res_rel_FlateG[n_eta-1][n_pt-1][n_lumi-1];
  double res_rel_H[n_eta-1][n_pt-1][n_lumi-1];

  double err_res_rel_BCD[n_eta-1][n_pt-1][n_lumi-1];
  double err_res_rel_EFearly[n_eta-1][n_pt-1][n_lumi-1];
  double err_res_rel_FlateG[n_eta-1][n_pt-1][n_lumi-1];
  double err_res_rel_H[n_eta-1][n_pt-1][n_lumi-1];

  TGraphErrors *rel_res_BCD[n_eta-1][n_pt-1]; 
  TGraphErrors *rel_res_EFearly[n_eta-1][n_pt-1]; 
  TGraphErrors *rel_res_FlateG[n_eta-1][n_pt-1]; 
  TGraphErrors *rel_res_H[n_eta-1][n_pt-1]; 



   double xbin_tgraph[n_lumi],zero[n_lumi];
   double xbin = 0;
   for(int i=0;i<n_lumi;i++){
     xbin = xbin + bin_width;
     xbin_tgraph[i] = xbin;
     zero[i] = 0;
   }

    for(int j=0; j<n_eta-1; j++){
      for(int i=0; i<n_pt-1; i++){
	for(int k =0; k<n_lumi; k++){
	  
	  res_rel_BCD[j][i][k]=0;
	  err_res_rel_BCD[j][i][k]=0;
	  res_rel_EFearly[j][i][k]=0;
	  err_res_rel_EFearly[j][i][k]=0;
	  res_rel_FlateG[j][i][k]=0;
	  err_res_rel_FlateG[j][i][k]=0;
	  res_rel_H[j][i][k]=0;
	  err_res_rel_H[j][i][k]=0;

	  res_mpf_BCD[j][i][k]=0;
	  err_res_mpf_BCD[j][i][k]=0;
	  res_mpf_EFearly[j][i][k]=0;
	  err_res_mpf_EFearly[j][i][k]=0;
	  res_mpf_FlateG[j][i][k]=0;
	  err_res_mpf_FlateG[j][i][k]=0;
	  res_mpf_H[j][i][k]=0;
	  err_res_mpf_H[j][i][k]=0;


	  if(hist_A_BCD_asymmetry[j][i]->GetEntries()>100){
	    res_rel_BCD[j][i][k] = (1+pr_A_BCD_asymmetry[j][i]->GetBinContent(k+1))/(1-pr_A_BCD_asymmetry[j][i]->GetBinContent(k+1));
	    err_res_rel_BCD[j][i][k] = 2/(pow((1-pr_A_BCD_asymmetry[j][i]->GetBinContent(k+1)),2)) * pr_A_BCD_asymmetry[j][i]->GetBinError(k+1);
	  }	 
	  if(hist_A_EFearly_asymmetry[j][i]->GetEntries()>100){
	    res_rel_EFearly[j][i][k] = (1+pr_A_EFearly_asymmetry[j][i]->GetBinContent(k+1))/(1-pr_A_EFearly_asymmetry[j][i]->GetBinContent(k+1));
	    err_res_rel_EFearly[j][i][k] = 2/(pow((1-pr_A_EFearly_asymmetry[j][i]->GetBinContent(k+1)),2)) * pr_A_EFearly_asymmetry[j][i]->GetBinError(k+1);
	  }
	  if(hist_A_FlateG_asymmetry[j][i]->GetEntries()>100){
	    res_rel_FlateG[j][i][k] = (1+pr_A_FlateG_asymmetry[j][i]->GetBinContent(k+1))/(1-pr_A_FlateG_asymmetry[j][i]->GetBinContent(k+1));
	    err_res_rel_FlateG[j][i][k] = 2/(pow((1-pr_A_FlateG_asymmetry[j][i]->GetBinContent(k+1)),2)) * pr_A_FlateG_asymmetry[j][i]->GetBinError(k+1);
	  }
	  if(hist_A_H_asymmetry[j][i]->GetEntries()>100){
	    res_rel_H[j][i][k] = (1+pr_A_H_asymmetry[j][i]->GetBinContent(k+1))/(1-pr_A_H_asymmetry[j][i]->GetBinContent(k+1));
	    err_res_rel_H[j][i][k] = 2/(pow((1-pr_A_H_asymmetry[j][i]->GetBinContent(k+1)),2)) * pr_A_H_asymmetry[j][i]->GetBinError(k+1);
	  }

	  rel_res_BCD[j][i]= new TGraphErrors(n_lumi, xbin_tgraph, res_rel_BCD[j][i], zero, err_res_rel_BCD[j][i]);
	  rel_res_BCD[j][i]=(TGraphErrors*)CleanEmptyPoints(rel_res_BCD[j][i]);
	  rel_res_EFearly[j][i]= new TGraphErrors(n_lumi, xbin_tgraph, res_rel_EFearly[j][i], zero, err_res_rel_EFearly[j][i]);
	  rel_res_EFearly[j][i]=(TGraphErrors*)CleanEmptyPoints(rel_res_EFearly[j][i]);
	  rel_res_FlateG[j][i]= new TGraphErrors(n_lumi, xbin_tgraph, res_rel_FlateG[j][i], zero, err_res_rel_FlateG[j][i]);
	  rel_res_FlateG[j][i]=(TGraphErrors*)CleanEmptyPoints(rel_res_FlateG[j][i]);
	  rel_res_H[j][i]= new TGraphErrors(n_lumi, xbin_tgraph, res_rel_H[j][i], zero, err_res_rel_H[j][i]);
	  rel_res_H[j][i]=(TGraphErrors*)CleanEmptyPoints(rel_res_H[j][i]);
	  
	  
	  if(hist_B_BCD_asymmetry[j][i]->GetEntries()>100){
	    res_mpf_BCD[j][i][k] = (1+pr_B_BCD_asymmetry[j][i]->GetBinContent(k+1))/(1-pr_B_BCD_asymmetry[j][i]->GetBinContent(k+1));
	    err_res_mpf_BCD[j][i][k] = 2/(pow((1-pr_B_BCD_asymmetry[j][i]->GetBinContent(k+1)),2)) * pr_B_BCD_asymmetry[j][i]->GetBinError(k+1);
	  }	 
	  if(hist_B_EFearly_asymmetry[j][i]->GetEntries()>100){
	    res_mpf_EFearly[j][i][k] = (1+pr_B_EFearly_asymmetry[j][i]->GetBinContent(k+1))/(1-pr_B_EFearly_asymmetry[j][i]->GetBinContent(k+1));
	    err_res_mpf_EFearly[j][i][k] = 2/(pow((1-pr_B_EFearly_asymmetry[j][i]->GetBinContent(k+1)),2)) * pr_B_EFearly_asymmetry[j][i]->GetBinError(k+1);
	  }
	  if(hist_B_FlateG_asymmetry[j][i]->GetEntries()>100){
	    res_mpf_FlateG[j][i][k] = (1+pr_B_FlateG_asymmetry[j][i]->GetBinContent(k+1))/(1-pr_B_FlateG_asymmetry[j][i]->GetBinContent(k+1));
	    err_res_mpf_FlateG[j][i][k] = 2/(pow((1-pr_B_FlateG_asymmetry[j][i]->GetBinContent(k+1)),2)) * pr_B_FlateG_asymmetry[j][i]->GetBinError(k+1);
	  }
	  if(hist_B_H_asymmetry[j][i]->GetEntries()>100){
	    res_mpf_H[j][i][k] = (1+pr_B_H_asymmetry[j][i]->GetBinContent(k+1))/(1-pr_B_H_asymmetry[j][i]->GetBinContent(k+1));
	    err_res_mpf_H[j][i][k] = 2/(pow((1-pr_B_H_asymmetry[j][i]->GetBinContent(k+1)),2)) * pr_B_H_asymmetry[j][i]->GetBinError(k+1);
	  }
	  
	  mpf_res_BCD[j][i]= new TGraphErrors(n_lumi, xbin_tgraph, res_mpf_BCD[j][i], zero, err_res_mpf_BCD[j][i]);
	  mpf_res_BCD[j][i]=(TGraphErrors*)CleanEmptyPoints(mpf_res_BCD[j][i]);
	  mpf_res_EFearly[j][i]= new TGraphErrors(n_lumi, xbin_tgraph, res_mpf_EFearly[j][i], zero, err_res_mpf_EFearly[j][i]);
	  mpf_res_EFearly[j][i]=(TGraphErrors*)CleanEmptyPoints(mpf_res_EFearly[j][i]);
	  mpf_res_FlateG[j][i]= new TGraphErrors(n_lumi, xbin_tgraph, res_mpf_FlateG[j][i], zero, err_res_mpf_FlateG[j][i]);
	  mpf_res_FlateG[j][i]=(TGraphErrors*)CleanEmptyPoints(mpf_res_FlateG[j][i]);
	  mpf_res_H[j][i]= new TGraphErrors(n_lumi, xbin_tgraph, res_mpf_H[j][i], zero, err_res_mpf_H[j][i]);
	  mpf_res_H[j][i]=(TGraphErrors*)CleanEmptyPoints(mpf_res_H[j][i]);
	}
  
  //dummy for tdrCanvas 
  TH1D *h = new TH1D("h",";dummy;",5000,0,40000);
  h->SetMaximum(-0.4);
  h->SetMinimum(0.4);
 
  TCanvas* c1 = new TCanvas();
  tdrCanvas(c1,"c1",h,4,10,true,CorrectionObject::_lumitag);

  TString alVal;
    alVal.Form("%0.2f\n",alpha_cut);
    TString axistitle_BCD = "RunBCD";
    TString axistitle_EFearly = "RunEFearly";
    TString axistitle_FlateG = "RunFlateG";
    TString axistitle_H = "RunH";
   
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 

    TLatex *tex1 = new TLatex();
    tex1->SetNDC();
    tex1->SetTextFont(42);
    tex1->SetTextSize(0.036); 
  
    h->GetXaxis()->SetTitle("Luminosity");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.80);
    h->GetXaxis()->SetLimits(0,35000);
    h->GetYaxis()->SetRangeUser(-0.4,0.4);
    h->GetYaxis()->SetTitle("<A>");

    pr_A_BCD_asymmetry[j][i]->SetLineColor(kBlue);
    pr_A_BCD_asymmetry[j][i]->Draw("SAME"); 

    pr_A_H_asymmetry[j][i]->SetLineColor(kRed);
    pr_A_H_asymmetry[j][i]->Draw("SAME");
     
    pr_A_EFearly_asymmetry[j][i]->SetLineColor(kGreen);
    pr_A_EFearly_asymmetry[j][i]->Draw("SAME");
    
    pr_A_FlateG_asymmetry[j][i]->SetLineColor(kViolet);
    pr_A_FlateG_asymmetry[j][i]->Draw("SAME");
 
    TLegend *leg_asymmetry;
    leg_asymmetry = new TLegend(0.35,0.72,0.51,0.90,"","brNDC");//x+0.1
    leg_asymmetry->SetBorderSize(0);
    leg_asymmetry->SetTextSize(0.036);
    leg_asymmetry->SetFillColor(10);
    leg_asymmetry->SetFillStyle(0);
    leg_asymmetry->SetLineColor(1);
    leg_asymmetry->SetTextFont(42);
    leg_asymmetry->SetHeader("A Response, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]+", #alpha<"+s_alpha_cut); //+", #alpha<"+s_alpha_cut
    leg_asymmetry->AddEntry(pr_A_BCD_asymmetry[j][i], axistitle_BCD,"lep"); 
    leg_asymmetry->AddEntry(pr_A_EFearly_asymmetry[j][i], axistitle_EFearly,"lep");
    leg_asymmetry->AddEntry(pr_A_FlateG_asymmetry[j][i], axistitle_FlateG,"lep");
    leg_asymmetry->AddEntry(pr_A_H_asymmetry[j][i], axistitle_H,"lep");
    leg_asymmetry->Draw();

    // tex->DrawLatex(0.53,0.95,CorrectionObject::_lumitag+"(13TeV)");
   tex1->DrawLatex(0.7,0.81, pt_range[i]+"<p_{T}<"+pt_range[i+1]);

    c1->SaveAs(CorrectionObject::_input_path+"Monitoring/A_symmetry_" + CorrectionObject::_generator_tag + "_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_pT_"+pt_range[i]+"_"+pt_range[i+1]+".pdf");
    delete c1;

 TCanvas* c2 = new TCanvas();
  tdrCanvas(c2,"c2",h,4,10,true,CorrectionObject::_lumitag);

    h->GetXaxis()->SetTitle("Luminosity");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.80);
    h->GetXaxis()->SetLimits(0,35000);
    h->GetYaxis()->SetRangeUser(-0.4,0.4);
    h->GetYaxis()->SetTitle("<B>");

   pr_B_BCD_asymmetry[j][i]->SetLineColor(kBlue);
    pr_B_BCD_asymmetry[j][i]->Draw("SAME"); 

    pr_B_H_asymmetry[j][i]->SetLineColor(kRed);
    pr_B_H_asymmetry[j][i]->Draw("SAME");
     
    pr_B_EFearly_asymmetry[j][i]->SetLineColor(kGreen);
    pr_B_EFearly_asymmetry[j][i]->Draw("SAME");
    
    pr_B_FlateG_asymmetry[j][i]->SetLineColor(kViolet);
    pr_B_FlateG_asymmetry[j][i]->Draw("SAME");

    TLegend *leg_bsymmetry;
    leg_bsymmetry = new TLegend(0.35,0.72,0.51,0.90,"","brNDC");//x+0.1
    leg_bsymmetry->SetBorderSize(0);
    leg_bsymmetry->SetTextSize(0.036);
    leg_bsymmetry->SetFillColor(10);
    leg_bsymmetry->SetFillStyle(0);
    leg_bsymmetry->SetLineColor(1);
    leg_bsymmetry->SetTextFont(42);
    leg_bsymmetry->SetHeader("B Response, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]); //+", #alpha<"+s_blpha_cut
    leg_bsymmetry->AddEntry(pr_B_BCD_asymmetry[j][i], axistitle_BCD,"lep"); 
    leg_bsymmetry->AddEntry(pr_B_EFearly_asymmetry[j][i], axistitle_EFearly,"lep");
    leg_bsymmetry->AddEntry(pr_B_FlateG_asymmetry[j][i], axistitle_FlateG,"lep");
    leg_bsymmetry->AddEntry(pr_B_H_asymmetry[j][i], axistitle_H,"lep");
    leg_bsymmetry->Draw();

    // tex->DrawLatex(0.53,0.95,CorrectionObject::_lumitag+"(13TeV)");
   tex1->DrawLatex(0.7,0.81, pt_range[i]+"<p_{T}<"+pt_range[i+1]);

    c2->SaveAs(CorrectionObject::_input_path+"/Monitoring/B_symmetry_" + CorrectionObject::_generator_tag + "_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_pT_"+pt_range[i]+"_"+pt_range[i+1]+".pdf");
    delete c2;



  TCanvas* c3 = new TCanvas();
  tdrCanvas(c3,"c3",h,4,10,true,CorrectionObject::_lumitag);

    h->GetXaxis()->SetTitle("Luminosity");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.80);
    h->GetXaxis()->SetLimits(0,35000);
    h->GetYaxis()->SetRangeUser(0.8,1.2);
    h->GetYaxis()->SetTitle("R_{pT-bal}");

    rel_res_BCD[j][i]->SetMarkerStyle(2);
    rel_res_BCD[j][i]->SetMarkerColor(kBlue);
    rel_res_BCD[j][i]->SetLineColor(kBlue);
    rel_res_BCD[j][i]->Draw("PE SAME"); 

    rel_res_H[j][i]->SetMarkerStyle(2);
    rel_res_H[j][i]->SetMarkerColor(kRed);
    rel_res_H[j][i]->SetLineColor(kRed);
    rel_res_H[j][i]->Draw("PE SAME");
     
    rel_res_EFearly[j][i]->SetMarkerStyle(2);
    rel_res_EFearly[j][i]->SetMarkerColor(kGreen);
    rel_res_EFearly[j][i]->SetLineColor(kGreen);
    rel_res_EFearly[j][i]->Draw("PE SAME");
    

    rel_res_FlateG[j][i]->SetMarkerStyle(2);
    rel_res_FlateG[j][i]->SetMarkerColor(kViolet);
    rel_res_FlateG[j][i]->SetLineColor(kViolet);
    rel_res_FlateG[j][i]->Draw("PE SAME");
 
    TLegend *leg_rel_res;
    leg_rel_res = new TLegend(0.35,0.72,0.51,0.90,"","brNDC");//x+0.1
    leg_rel_res->SetBorderSize(0);
    leg_rel_res->SetTextSize(0.036);
    leg_rel_res->SetFillColor(10);
    leg_rel_res->SetFillStyle(0);
    leg_rel_res->SetLineColor(1);
    leg_rel_res->SetTextFont(42);
    leg_rel_res->SetHeader("Rel Response, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]+", #alpha<"+s_alpha_cut); //+", #alpha<"+s_alpha_cut
    leg_rel_res->AddEntry(rel_res_BCD[j][i], axistitle_BCD,"lep"); 
    leg_rel_res->AddEntry(rel_res_EFearly[j][i], axistitle_EFearly,"lep");
    leg_rel_res->AddEntry(rel_res_FlateG[j][i], axistitle_FlateG,"lep");
    leg_rel_res->AddEntry(rel_res_H[j][i], axistitle_H,"lep");
    leg_rel_res->Draw();
   
    c3->SaveAs(CorrectionObject::_input_path+"Monitoring/Rel_res_" + CorrectionObject::_generator_tag + "_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_pT_"+pt_range[i]+"_"+pt_range[i+1]+".pdf");
    delete c3;

 TCanvas* c4 = new TCanvas();
  tdrCanvas(c4,"c4",h,4,10,true,CorrectionObject::_lumitag);

    h->GetXaxis()->SetTitle("Luminosity");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.80);
    h->GetXaxis()->SetLimits(0,35000);
    h->GetYaxis()->SetRangeUser(0.8,1.2);
    h->GetYaxis()->SetTitle("R_{MPF}");

    mpf_res_BCD[j][i]->SetMarkerStyle(2);
    mpf_res_BCD[j][i]->SetMarkerColor(kBlue);
    mpf_res_BCD[j][i]->SetLineColor(kBlue);
    mpf_res_BCD[j][i]->Draw("PE SAME"); 

    mpf_res_H[j][i]->SetMarkerStyle(2);
    mpf_res_H[j][i]->SetMarkerColor(kRed);
    mpf_res_H[j][i]->SetLineColor(kRed);
    mpf_res_H[j][i]->Draw("PE SAME");
     
    mpf_res_EFearly[j][i]->SetMarkerStyle(2);
    mpf_res_EFearly[j][i]->SetMarkerColor(kGreen);
    mpf_res_EFearly[j][i]->SetLineColor(kGreen);
    mpf_res_EFearly[j][i]->Draw("PE SAME");
    

    mpf_res_FlateG[j][i]->SetMarkerStyle(2);
    mpf_res_FlateG[j][i]->SetMarkerColor(kViolet);
    mpf_res_FlateG[j][i]->SetLineColor(kViolet);
    mpf_res_FlateG[j][i]->Draw("PE SAME");

    TLegend *leg_mpf_res;
    leg_mpf_res = new TLegend(0.35,0.72,0.51,0.90,"","brNDC");//x+0.1
    leg_mpf_res->SetBorderSize(0);
    leg_mpf_res->SetTextSize(0.036);
    leg_mpf_res->SetFillColor(10);
    leg_mpf_res->SetFillStyle(0);
    leg_mpf_res->SetLineColor(1);
    leg_mpf_res->SetTextFont(42);
    leg_mpf_res->SetHeader("MPF Response, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]); //+", #alpha<"+s_blpha_cut
    leg_mpf_res->AddEntry(mpf_res_BCD[j][i], axistitle_BCD,"lep"); 
    leg_mpf_res->AddEntry(mpf_res_EFearly[j][i], axistitle_EFearly,"lep");
    leg_mpf_res->AddEntry(mpf_res_FlateG[j][i], axistitle_FlateG,"lep");
    leg_mpf_res->AddEntry(mpf_res_H[j][i], axistitle_H,"lep");
    leg_mpf_res->Draw();

    // tex->DrawLatex(0.53,0.95,CorrectionObject::_lumitag+"(13TeV)");
   tex1->DrawLatex(0.7,0.81, pt_range[i]+"<p_{T}<"+pt_range[i+1]);

    c4->SaveAs(CorrectionObject::_input_path+"/Monitoring/MPF_res_" + CorrectionObject::_generator_tag + "_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_pT_"+pt_range[i]+"_"+pt_range[i+1]+".pdf");
    delete c4;


   }
   }
 }
 
