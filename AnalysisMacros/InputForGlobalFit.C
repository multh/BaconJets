// This script stores ratio between MC and DATA responses as function of pT in different eta bins and for different alpha values

#include "header.h"
#include "TGraph.h"
#include <iostream>
#include <cmath>
#include "TString.h"
#include "TH1D.h"

using namespace std;
void InputForGlobalFit(TString path, TFile* datafile, TFile* MCfile){
  const int min_number_events = 20;//required number of events in response histogram
  gStyle->SetOptFit(0);

  //calculate ratio in MC to DATA responses

  double mc_al_mpf[n_alpha_common][n_eta_common-1][n_pt-1]; 
  double err_mc_al_mpf[n_alpha_common][n_eta_common-1][n_pt-1]; 
  double mc_al_pTbal[n_alpha_common][n_eta_common-1][n_pt-1]; 
  double err_mc_al_pTbal[n_alpha_common][n_eta_common-1][n_pt-1]; 
  double data_al_mpf[n_alpha_common][n_eta_common-1][n_pt-1]; 
  double err_data_al_mpf[n_alpha_common][n_eta_common-1][n_pt-1]; 
  double data_al_pTbal[n_alpha_common][n_eta_common-1][n_pt-1]; 
  double err_data_al_pTbal[n_alpha_common][n_eta_common-1][n_pt-1]; 

  double ratio_al_mpf[n_alpha_common][n_eta_common-1][n_pt-1]; //ratio at pt,eta,alpha bins
  double err_ratio_al_mpf[n_alpha_common][n_eta_common-1][n_pt-1]; //error of ratio at pt,eta,alpha bins
  double ratio_al_pTbal[n_alpha_common][n_eta_common-1][n_pt-1]; //ratio at pt,eta,alpha bins
  double err_ratio_al_pTbal[n_alpha_common][n_eta_common-1][n_pt-1]; //error of ratio at pt,eta,alpha bins


  //Initialize with 0 values
  for(int i=0; i<n_alpha_common; i++){
    for(int j=0; j<n_eta_common-1; j++){
      for(int k=0; k<n_pt-1; k++){
	ratio_al_mpf[i][j][k] = 0;
	err_ratio_al_mpf[i][j][k] = 0;
	ratio_al_pTbal[i][j][k] = 0;
	err_ratio_al_pTbal[i][j][k] = 0;
	mc_al_mpf[i][j][k] = 0;
	err_mc_al_mpf[i][j][k] = 0; 
	mc_al_pTbal[i][j][k] = 0;  
	err_mc_al_pTbal[i][j][k] = 0;
	data_al_mpf[i][j][k] = 0; 
	err_data_al_mpf[i][j][k] = 0; 
	data_al_pTbal[i][j][k] = 0; 
	err_data_al_pTbal[i][j][k] = 0; 
      }
    }
  }

  for(int i=0; i<n_alpha_common; i++){
    for(int j=0; j<n_eta_common-1; j++){
      for(int k=0; k<n_pt-1; k++){
	cout<<"For alpha = "<<alpha_range_common[i]<<" eta = "<<eta_common_range[j]<<", "<<eta_common_range[j+1]<<" pT = "<<pt_range[k]<<", "<<pt_range[k+1]<<endl;
	pair<double,double> res_mc_mpf = Response(MCfile,alpha_bins_common[i],eta_common_bins[j],eta_common_bins[j+1],pt_bins[k],pt_bins[k+1],true);
	mc_al_mpf[i][j][k] = res_mc_mpf.first;
	err_mc_al_mpf[i][j][k] = res_mc_mpf.second;
	pair<double,double> res_mc_pTbal = Response(MCfile,alpha_bins_common[i],eta_common_bins[j],eta_common_bins[j+1],pt_bins[k],pt_bins[k+1],false);
	mc_al_pTbal[i][j][k] = res_mc_pTbal.first;  
	err_mc_al_pTbal[i][j][k] = res_mc_pTbal.second;
	pair<double,double> res_data_mpf = Response(datafile,alpha_bins_common[i],eta_common_bins[j],eta_common_bins[j+1],pt_bins[k],pt_bins[k+1],true);
	data_al_mpf[i][j][k] = res_data_mpf.first;
	err_data_al_mpf[i][j][k] = res_data_mpf.second;
	pair<double,double> res_data_pTbal = Response(datafile,alpha_bins_common[i],eta_common_bins[j],eta_common_bins[j+1],pt_bins[k],pt_bins[k+1],false);
	data_al_pTbal[i][j][k] = res_data_pTbal.first;  
	err_data_al_pTbal[i][j][k] = res_data_pTbal.second;

	pair<double,double> ratio_res_mpf = Rmc_to_Rdata(res_mc_mpf,res_data_mpf);
	ratio_al_mpf[i][j][k] = ratio_res_mpf.first;
	err_ratio_al_mpf[i][j][k] = ratio_res_mpf.second;
	pair<double,double> ratio_res_pTbal = Rmc_to_Rdata(res_mc_pTbal,res_data_pTbal);
	ratio_al_pTbal[i][j][k] = ratio_res_pTbal.first;
	err_ratio_al_pTbal[i][j][k] = ratio_res_pTbal.second;
      }
    }
  }


  //Store results in TGraphErrors
  double res_xbin_tgraph[n_pt-1];// = {(pt_bins[0]+pt_bins[1])/2
  double res_zero[n_pt-1];
  for(int i=0;i<n_pt-1;i++){
    res_xbin_tgraph[i]=(pt_bins[i]+pt_bins[i+1])/2;
    res_zero[i]=(pt_bins[i+1]-pt_bins[i])/2 ;
  }

  TGraphErrors *mpf_data[n_alpha_common-1][n_eta_common-1];
  TGraphErrors *rrel_data[n_alpha_common-1][n_eta_common-1];
  TGraphErrors *mpf_mc[n_alpha_common-1][n_eta_common-1];
  TGraphErrors *rrel_mc[n_alpha_common-1][n_eta_common-1];
  TGraphErrors *mpf_ratio[n_alpha_common-1][n_eta_common-1];
  TGraphErrors *rrel_ratio[n_alpha_common-1][n_eta_common-1];

  for(int i=0; i<n_alpha_common-1; i++){
    for(int j=0; j<n_eta_common-1; j++){
      mpf_data[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph,data_al_mpf[i][j],res_zero,err_data_al_mpf[i][j]);
      mpf_mc[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph,mc_al_mpf[i][j],res_zero,err_mc_al_mpf[i][j]);
      rrel_data[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph,data_al_pTbal[i][j],res_zero,err_data_al_pTbal[i][j]);
      rrel_mc[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph,mc_al_pTbal[i][j],res_zero,err_mc_al_pTbal[i][j]);
      mpf_ratio[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph,ratio_al_mpf[i][j],res_zero,err_ratio_al_mpf[i][j]);
      rrel_ratio[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph,ratio_al_pTbal[i][j],res_zero,err_ratio_al_pTbal[i][j]);

    }
  }

  // Cleaning for empty points (with low statistic)
  for(int i=0; i<n_alpha_common-1; i++){
    for(int j=0; j<n_eta_common-1; j++){
      mpf_data[i][j] = CleanEmptyPoints(mpf_data[i][j]);
      mpf_mc[i][j] = CleanEmptyPoints(mpf_mc[i][j]);
      rrel_data[i][j] = CleanEmptyPoints(rrel_data[i][j]);
      rrel_mc[i][j] = CleanEmptyPoints(rrel_mc[i][j]);
      mpf_ratio[i][j] = CleanEmptyPoints(mpf_ratio[i][j]);
      rrel_ratio[i][j] = CleanEmptyPoints(rrel_ratio[i][j]);

      mpf_data[i][j]->SetName("mpfchs_dijet_"+alpha_range_common[i]);
      mpf_mc[i][j]->SetName("mpfchs_dijet_"+alpha_range_common[i]);
      rrel_data[i][j]->SetName("ptchs_dijet_"+alpha_range_common[i]);
      rrel_mc[i][j]->SetName("ptchs_dijet_"+alpha_range_common[i]);
      mpf_ratio[i][j]->SetName("mpfchs_dijet_"+alpha_range_common[i]);
      rrel_ratio[i][j]->SetName("ptchs_dijet_"+alpha_range_common[i]);
      }
    }
 


  //Save results in root file
  TFile* outputfile = new TFile(path+"JECcombifile_Dijet.root","RECREATE");
  outputfile->Print();
 // TString eta_output[n_eta_common-1] = {"eta0000-0261", "eta0216-0522","eta0522-0783","eta0783-1044","eta1044-1305","eta1305-1653","eta1653-1930","eta1930-2172","eta2172-2322","eta2322-2500","eta2500-2650","eta2650-2853","eta2853-2964","eta2964-3139",
 // 				"eta3139-3489","eta3489-5191"};//TMP
  for(int i=0; i<n_alpha_common-1; i++){
    for(int j=0; j<n_eta_common-1; j++){
      if(i==0){
	outputfile->mkdir("ratio/"+eta_output[j]);
	outputfile->mkdir("data/"+eta_output[j]);
	outputfile->mkdir("mc/"+eta_output[j]);
      }
      outputfile->cd("ratio/"+eta_output[j]);
      mpf_ratio[i][j]->Write();
      rrel_ratio[i][j]->Write();
    
      outputfile->cd("data/"+eta_output[j]);
      mpf_data[i][j]->Write();
      rrel_data[i][j]->Write();
      
      outputfile->cd("mc/"+eta_output[j]);
      mpf_mc[i][j]->Write();
      rrel_mc[i][j]->Write();
    }
 }

  cout<<"Draw result for alpha = "<<alpha_range_common[3]<<" eta = "<<eta_common_range[n_eta_common-2]<<" "<<eta_common_range[n_eta_common-1]<<endl;
  mpf_ratio[3][n_eta_common-2]->Draw();

  outputfile->Write();
  outputfile->Close();

}
