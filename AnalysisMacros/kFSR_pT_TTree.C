// with this script kFSR extrapolations over alpha in bins of eta will be calculated for different pT bins
// works for MPF and pt-balance methods separetly

#include "header.h"
#include <iostream>
#include <cmath>
#include "TString.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
//#include "UsefulFunctions.C"
using namespace std;

void kFSR_pT_TTree(bool mpfMethod, TString path, TFile* datafile, TFile* MCfile, double al_cut=0.2,int nResponseBins=100){
  cout<<"alpha cut @"<<al_cut<<endl;
  TStyle* m_gStyle = new TStyle();

  m_gStyle->SetOptFit(0);
  // get ratio for MC to DATA responses
  double ratio_al_rel_r[n_pt-1][n_eta-1][n_alpha]; //ratio at pt,eta,alpha bins
  double err_ratio_al_rel_r[n_pt-1][n_eta-1][n_alpha]; //error of ratio at pt,eta,alpha bins
  double ratio_al_mpf_r[n_pt-1][n_eta-1][n_alpha]; //ratio at pt,eta,alpha bins
  double err_ratio_al_mpf_r[n_pt-1][n_eta-1][n_alpha]; //error of ratio at pt,eta,alpha bins
  TH1D *hdata_rel_r[n_pt-1][n_eta-1][n_alpha];// pT-balance responce for data
  TH1D *hdata_mpf_r[n_pt-1][n_eta-1][n_alpha];//MPF responce for data
  TH1D *hmc_rel_r[n_pt-1][n_eta-1][n_alpha];// pT-balance responce for MC
  TH1D *hmc_mpf_r[n_pt-1][n_eta-1][n_alpha];//MPF responce for MC
  int count = 0;
  TString name1 = "hist_data_rel_r_";
  TString name2 = "hist_data_mpf_r_";
  TString name3 = "hist_mc_rel_r_";
  TString name4 = "hist_mc_mpf_r_";
  for(int i=0; i<n_alpha; i++){
    for(int j=0; j<n_eta-1; j++){
      for(int k=0; k<n_pt-1; k++){
	ratio_al_rel_r[k][j][i] = 0;
	err_ratio_al_rel_r[k][j][i] = 0;
	ratio_al_mpf_r[k][j][i] = 0;
	err_ratio_al_mpf_r[k][j][i] = 0;
	TString name = name1; name+=count;
	hdata_rel_r[k][j][i] = new TH1D(name,"",nResponseBins, 0, 2.5);
	name = name2;name+=count;
	hdata_mpf_r[k][j][i] = new TH1D(name,"",nResponseBins, 0, 2.5);
	name = name3; name+=count;
	hmc_rel_r[k][j][i] = new TH1D(name,"",nResponseBins, 0, 2.5);
	name = name4; name+=count;
	hmc_mpf_r[k][j][i] = new TH1D(name,"",nResponseBins, 0, 2.5);
	count++;
      }
    }
  }


   // Create the tree reader and its data containers
   TTreeReader myReader_DATA("AnalysisTree", datafile);
   TTreeReaderValue<Float_t> pt_ave_data(myReader_DATA, "pt_ave");
   TTreeReaderValue<Float_t> probejet_eta_data(myReader_DATA, "probejet_eta");
   TTreeReaderValue<Float_t> alpha_data(myReader_DATA, "alpha");
   //   TTreeReaderValue<Float_t> alpha_data(myReader_DATA, "alpha_sum"); //TEST alpha_sum
   TTreeReaderValue<Float_t> rel_r_data(myReader_DATA, "rel_r");
   TTreeReaderValue<Float_t> mpf_r_data(myReader_DATA, "mpf_r");
   TTreeReaderValue<Float_t> weight_data(myReader_DATA, "weight");
   
   while (myReader_DATA.Next()) {
   //   for( int i=0;i<10;i++) {//TEST
   //     myReader_DATA.Next();//TEST
     //      cout<<"DATA point: alpha = "<<*alpha_data<<" eta_probe = "<<*probejet_eta_data<<" pT_ave = "<<*pt_ave_data<<endl;
     for(int k=0; k<n_pt-1; k++){
   	   if(*pt_ave_data<pt_bins[k] || *pt_ave_data>pt_bins[k+1]) continue;
	   for(int j=0; j<n_eta-1; j++){
	     if(fabs(*probejet_eta_data)>eta_bins[j+1] || fabs(*probejet_eta_data)<eta_bins[j]) continue;
	     //  if(*probejet_eta_data>eta_bins[j+1] || *probejet_eta_data<eta_bins[j]) continue;//TEST symmetry
	     for(int i=0; i<n_alpha; i++){
	       if(*alpha_data>alpha_bins[i]) continue;
	       else{
		 //	 std::cout<<"DATA: For alpha cut at "<<alpha_bins[i]<<" eta range = "<<eta_bins[j]<<", "<<eta_bins[j+1]
		 //		  <<" and pT range = "<<pt_bins[k]<<", "<<pt_bins[k+1]<<std::endl;
		 //	 std::cout<<"probe_eta = "<<*probejet_eta_data<<" pT_ave = "<<*pt_ave_data<<std::endl;
		 //		 cout<<"DATA response "<<*rel_r_data<<endl;
		 hdata_rel_r[k][j][i]->Fill(*rel_r_data,*weight_data);
		 hdata_mpf_r[k][j][i]->Fill(*mpf_r_data,*weight_data);
	       }
	   }
	 }
     }
   }

   TTreeReader myReader_MC("AnalysisTree", MCfile);
   TTreeReaderValue<Float_t> pt_ave_mc(myReader_MC, "pt_ave");
   TTreeReaderValue<Float_t> probejet_eta_mc(myReader_MC, "probejet_eta");
   TTreeReaderValue<Float_t> alpha_mc(myReader_MC, "alpha");
   //   TTreeReaderValue<Float_t> alpha_mc(myReader_MC, "alpha_sum");//TEST alpha_sum
   TTreeReaderValue<Float_t> rel_r_mc(myReader_MC, "rel_r");
   TTreeReaderValue<Float_t> mpf_r_mc(myReader_MC, "mpf_r");
   TTreeReaderValue<Float_t> weight_mc(myReader_MC, "weight");
   while (myReader_MC.Next()) {
     //     cout<<"*alpha_mc = "<<*alpha_mc<<endl;
     for(int k=0; k<n_pt-1; k++){
       if(*pt_ave_mc<pt_bins[k] || *pt_ave_mc>pt_bins[k+1]) continue;
       for(int j=0; j<n_eta-1; j++){
	 if(fabs(*probejet_eta_mc)>eta_bins[j+1] || fabs(*probejet_eta_mc)<eta_bins[j]) continue;
	 //	 if(*probejet_eta_mc>eta_bins[j+1] || *probejet_eta_mc<eta_bins[j]) continue; //TEST symmetry
   	 for(int i=0; i<n_alpha; i++){
   	   if(*alpha_mc>alpha_bins[i]) continue;
   	   else{
	     //	     cout<<"MC response "<<*rel_r_mc<<endl;
   	     hmc_rel_r[k][j][i]->Fill(*rel_r_mc,*weight_mc);
   	     hmc_mpf_r[k][j][i]->Fill(*mpf_r_mc,*weight_mc);
   	   }
   	 }
       }
     }
   }

   //TEST: normilise hists
   for(int i=0; i<n_alpha; i++){
     for(int j=0; j<n_eta-1; j++){
       for(int k=0; k<n_pt-1; k++){
	 double integ_rel_mc = hmc_rel_r[k][j][i]->Integral();
	 hmc_rel_r[k][j][i]->Scale(1./integ_rel_mc);
	 double integ_mpf_mc = hmc_mpf_r[k][j][i]->Integral();
	 hmc_mpf_r[k][j][i]->Scale(1./integ_mpf_mc);
	 double integ_rel_data = hdata_rel_r[k][j][i]->Integral();
	 hdata_rel_r[k][j][i]->Scale(1./integ_rel_data);
	 double integ_mpf_data = hdata_mpf_r[k][j][i]->Integral();
	 hdata_mpf_r[k][j][i]->Scale(1./integ_mpf_data);
       }
     }
   }

   for(int i=0; i<n_alpha; i++){
     for(int j=0; j<n_eta-1; j++){
      for(int k=0; k<n_pt-1; k++){
	if(fabs(alpha_bins[i]-al_cut)<1e-4){
	  if(k==0) cout<<eta_bins[j]<<" ";
	  //cout<<"& "<<hmc_rel_r[k][j][i]->GetEntries();//MC
	  cout<<"& "<<hdata_rel_r[k][j][i]->GetEntries();//DATA
	}
	pair<double,double> res_mc_rel_r,res_data_rel_r;
	pair<double,double> res_mc_mpf_r,res_data_mpf_r;
	res_mc_rel_r = GetValueAndError(hmc_rel_r[k][j][i]);
	res_data_rel_r = GetValueAndError(hdata_rel_r[k][j][i]);
	res_mc_mpf_r = GetValueAndError(hmc_mpf_r[k][j][i]);
	res_data_mpf_r = GetValueAndError(hdata_mpf_r[k][j][i]);

	pair<double,double> ratio_res_rel_r;
	if(res_mc_rel_r.first>0 && res_data_rel_r.first>0)
	  ratio_res_rel_r = Rmc_to_Rdata(res_mc_rel_r,res_data_rel_r);
	else 
	  ratio_res_rel_r.first = 0;
	pair<double,double> ratio_res_mpf_r;
	if(res_mc_mpf_r.first>0 && res_data_mpf_r.first>0)
	  ratio_res_mpf_r = Rmc_to_Rdata(res_mc_mpf_r,res_data_mpf_r);
	else 
	  ratio_res_mpf_r.first = 0;

	ratio_al_rel_r[k][j][i] = ratio_res_rel_r.first;
	err_ratio_al_rel_r[k][j][i] = ratio_res_rel_r.second;
	ratio_al_mpf_r[k][j][i] = ratio_res_mpf_r.first;
	err_ratio_al_mpf_r[k][j][i] = ratio_res_mpf_r.second;
      }
      if(fabs(alpha_bins[i]-al_cut)<1e-4) cout<<""<<endl;
     }
   }

   //Divide reponses by value at alpha=al_cut
   //find bin with alpha = al_cut
   int al_ref=0;
   for(int i=0; i<n_alpha; i++){
     if(fabs(alpha_bins[i]-al_cut)<1e-4) 
       al_ref=i;
   }
   // cout<<"alpha="<<al_cut<" for bin#"<<al_ref<<endl;

   // //Divide reponses by value at alpha=0.3
   // //find bin with alpha = 0.3
   // int al_ref=0;
   // for(int i=0; i<n_alpha; i++){
   //   if(fabs(alpha_bins[i]-0.3)<1e-4) al_ref=i;
   // }
   // cout<<"alpha=0.3 for bin#"<<al_ref<<endl;

  for(int j=0; j<n_eta-1; j++){
    for(int k=0; k<n_pt-1; k++){
      double norm_al02_rel_r = ratio_al_rel_r[k][j][al_ref];
      double err_norm_al02_rel_r = err_ratio_al_rel_r[k][j][al_ref];
      double norm_al02_mpf_r = ratio_al_mpf_r[k][j][al_ref];
      double err_norm_al02_mpf_r = err_ratio_al_mpf_r[k][j][al_ref];
      for(int i=0; i<n_alpha; i++){
  	if(norm_al02_rel_r>0){
  	  ratio_al_rel_r[k][j][i] =   ratio_al_rel_r[k][j][i]/norm_al02_rel_r;
  	  err_ratio_al_rel_r[k][j][i] = sqrt(abs(pow(err_ratio_al_rel_r[k][j][i],2)-pow(err_norm_al02_rel_r,2)));
  	}
  	if(norm_al02_mpf_r>0){
  	  ratio_al_mpf_r[k][j][i] =   ratio_al_mpf_r[k][j][i]/norm_al02_mpf_r;
  	  err_ratio_al_mpf_r[k][j][i] = sqrt(abs(pow(err_ratio_al_mpf_r[k][j][i],2)-pow(err_norm_al02_mpf_r,2)));
  	}
      }
    }
  }
  
  TLegend *leg1;
  leg1 = new TLegend(0.17,0.68,0.65,0.89,"","brNDC");//x+0.1
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetFillColor(10);
  leg1->SetLineColor(1);
  leg1->SetTextFont(42);
  leg1->SetNColumns(2);
  TGraphErrors *graph_rel_r[n_pt-1][n_eta-1];//set of points vs alpha
  TMultiGraph *pTgraph_rel_r[n_eta-1];//set of different pT bins in on eta bin
  TGraphErrors *graph_mpf_r[n_pt-1][n_eta-1];//set of points vs alpha
  TMultiGraph *pTgraph_mpf_r[n_eta-1];//set of different pT bins in on eta bin
  double xbin_tgraph[n_alpha],zero[n_alpha];
  for(int i=0;i<n_alpha;i++){
    xbin_tgraph[i] = alpha_bins[i];
    zero[i] = 0;
  }
  for(int j=0; j<n_eta-1; j++){
    pTgraph_rel_r[j] = new TMultiGraph();
    pTgraph_mpf_r[j] = new TMultiGraph();
    for(int k=0; k<n_pt-1; k++){
      graph_rel_r[k][j] = new TGraphErrors(n_alpha,xbin_tgraph,ratio_al_rel_r[k][j],zero,err_ratio_al_rel_r[k][j]);
      graph_rel_r[k][j] = (TGraphErrors*)CleanEmptyPoints(graph_rel_r[k][j]);
      graph_rel_r[k][j]->SetMarkerSize(1.5);
      graph_rel_r[k][j]->SetMarkerStyle(20);
      if(k<9){
      graph_rel_r[k][j]->SetMarkerColor(k+1);
      graph_rel_r[k][j]->SetLineColor(k+1);
      }
      else{
	graph_rel_r[k][j]->SetMarkerColor(k+19);
	graph_rel_r[k][j]->SetLineColor(k+19);
      }
      if(graph_rel_r[k][j]->GetN()>0) 
	pTgraph_rel_r[j]->Add(graph_rel_r[k][j]);
      TString pTbin_label = "";
      pTbin_label+=pt_bins[k];
      pTbin_label+=" < p_{T} < ";
      pTbin_label+=pt_bins[k+1];
      //      cout<<j<<" "<<k<<" "<<pTbin_label<<endl;
      if(j==0) leg1->AddEntry(graph_rel_r[k][j],pTbin_label,"epl");

      graph_mpf_r[k][j] = new TGraphErrors(n_alpha,xbin_tgraph,ratio_al_mpf_r[k][j],zero,err_ratio_al_mpf_r[k][j]);
      graph_mpf_r[k][j] = (TGraphErrors*)CleanEmptyPoints(graph_mpf_r[k][j]);
      graph_mpf_r[k][j]->SetMarkerSize(1.3);
      graph_mpf_r[k][j]->SetMarkerStyle(20);
      if(k<9){
      graph_mpf_r[k][j]->SetMarkerColor(k+1);
      graph_mpf_r[k][j]->SetLineColor(k+1);
      }
      else{
	graph_mpf_r[k][j]->SetMarkerColor(k+19);
	graph_mpf_r[k][j]->SetLineColor(k+19);
      }
      if(graph_mpf_r[k][j]->GetN()>0) 
	pTgraph_mpf_r[j]->Add(graph_mpf_r[k][j]);

    }
  }

// // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//   // // Plots

//   // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(alpha_bins[0]-0.01,1,alpha_bins[n_alpha-1]+0.01,1);
  // create output .dat file, including the kFSR extrapolation (alpha->0)
  FILE *fp_rel_r; FILE *fp_mpf_r; 
  TH1D* kFSR_MPF; TH1D* kFSR_DiJet;
  fp_mpf_r = fopen(path+"output/KFSR_MPF_extrapolation.dat","w");
  fp_rel_r = fopen(path+"output/KFSR_DiJet_extrapolation.dat","w");
  kFSR_MPF = new TH1D("kfsr_mpf","kfsr_mpf", n_eta-1,eta_bins);
  kFSR_DiJet = new TH1D("kfsr_dijet","kfsr_dijet", n_eta-1,eta_bins);
  TH1D* plotkfsr = new TH1D("kfsr","kfsr", n_eta-1,eta_bins);


  // Start with pT-balance
  //create plots
  TCanvas* a[n_eta-1];
  TString plotname[n_eta-1];
  TF1 *pol1[n_eta-1];
  for (int j=0; j<n_eta-1; j++){
    plotname[j]="dijet_kfsr_diffPt_eta_"+eta_range[j]+"_"+eta_range[j+1];
    a[j] = new TCanvas(plotname[j], plotname[j], 800,700);
    m_gStyle->SetOptTitle(0);
    pTgraph_rel_r[j]->Draw("AP");
    pTgraph_rel_r[j]->GetYaxis()->SetRangeUser(0.92,1.08);
    //    pTgraph_rel_r[j]->GetYaxis()->SetRangeUser(0.9,1.1);
    pTgraph_rel_r[j]->GetXaxis()->SetRangeUser(0,alpha_bins[n_alpha-1]+0.01);
    //    pTgraph_rel_r[j]->GetYaxis()->SetTitle("(R_{MC}/R_{DATA})/(R_{MC}/R_{DATA})_{#alpha<0.2}");
    pTgraph_rel_r[j]->GetYaxis()->SetTitle("(R_{MC}/R_{DATA})");
    pTgraph_rel_r[j]->GetXaxis()->SetTitle("cut on #alpha");
    //    pol1[j] = new TF1("pol1","pol1",0.19,0.36);  //TEST AK4
    pol1[j] = new TF1("pol1","pol1",0.09,0.36);  //TEST AK4
    //    pol1[j] = new TF1("pol1","pol1",0.19,0.46);  //TEST AK8CHS
    //    pol1[j] = new TF1("pol1","pol1",0.19,0.41); //TEST AK8PUPPI

    pol1[j]->SetParameters(0,0);
    pTgraph_rel_r[j]->Fit(pol1[j],"R");
    line->SetLineStyle(2);
    line->Draw("SAME");

    // //Divide by value@alpha cut
    // double y_al_cut = pol1[j]->GetParameter(0)+pol1[j]->GetParameter(1)*al_cut;
    //    cout<<"y_al_cut = "<<y_al_cut<<endl;
    // fill the output.dat file
    if (fp_rel_r!=NULL) {
      Float_t value = pol1[j]->GetParameter(0);
      //      Float_t value = pol1[j]->GetParameter(0)/y_al_cut;
      Float_t uncert = pol1[j]->GetParError(0);
      fprintf(fp_rel_r, "%f %f\n",value,uncert);
    }
    plotkfsr->SetBinContent(j+1,pol1[j]->GetParameter(0));
    plotkfsr->SetBinError(j+1,pol1[j]->GetParError(0));
    kFSR_DiJet->SetBinContent(j+1,pol1[j]->GetParameter(0));
    kFSR_DiJet->SetBinError(j+1,pol1[j]->GetParError(0));
    // plotkfsr->SetBinContent(j+1,pol1[j]->GetParameter(0)/y_al_cut);
    // plotkfsr->SetBinError(j+1,pol1[j]->GetParError(0));
    // kFSR_DiJet->SetBinContent(j+1,pol1[j]->GetParameter(0)/y_al_cut);
    // kFSR_DiJet->SetBinError(j+1,pol1[j]->GetParError(0));
    leg1->SetHeader("p_{T} balance, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
    leg1->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    tex->DrawLatex(0.64,0.91,"2.11fb^{-1} (13TeV)");


    TString chi2_loglin = "#chi^{2}/n.d.f = ";
    chi2_loglin += trunc(pol1[j]->GetChisquare());
    chi2_loglin +="/";
    chi2_loglin +=trunc(pol1[j]->GetNDF());
    TLatex *tex2 = new TLatex();
    tex2->SetNDC();
    tex2->SetTextSize(0.035); 
    tex2->DrawLatex(0.64,0.35,chi2_loglin);
    a[j]->Print(path+"plots/kFSR_Pt_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
  }
  fclose(fp_rel_r);
  // create output file including the kFSR plot
  TFile* outputfile_rel_r = new TFile(path+"Histo_KFSR_DiJet_L1.root","RECREATE");
  kFSR_DiJet->Write();
  outputfile_rel_r->Write();
  outputfile_rel_r->Close();


  // And now MPF results
 //create plots
  TCanvas* b[n_eta-1];
  for (int j=0; j<n_eta-1; j++){
    plotname[j]="mpf_kfsr_diffPt_eta_"+eta_range[j]+"_"+eta_range[j+1];
    b[j] = new TCanvas(plotname[j], plotname[j], 800,700);
    m_gStyle->SetOptTitle(0);
    pTgraph_mpf_r[j]->Draw("AP");
    pTgraph_mpf_r[j]->GetYaxis()->SetRangeUser(0.92,1.08);
    //    pTgraph_mpf_r[j]->GetYaxis()->SetRangeUser(0.9,1.1);
    pTgraph_mpf_r[j]->GetXaxis()->SetRangeUser(0,alpha_bins[n_alpha-1]+0.01);
    //    pTgraph_mpf_r[j]->GetYaxis()->SetTitle("(R_{MC}/R_{DATA})/(R_{MC}/R_{DATA})_{#alpha<0.2}");
    pTgraph_mpf_r[j]->GetYaxis()->SetTitle("(R_{MC}/R_{DATA})");
    pTgraph_mpf_r[j]->GetXaxis()->SetTitle("cut on #alpha");
    //    pol1[j] = new TF1("pol1","pol1",0.19,0.36);  //TEST AK4
    pol1[j] = new TF1("pol1","pol1",0.09,0.36);  //TEST AK4
    //    pol1[j] = new TF1("pol1","pol1",0.19,0.46);  //TEST AK8CHS
    //    pol1[j] = new TF1("pol1","pol1",0.19,0.41); //TEST AK8PUPPI
    pol1[j]->SetParameters(0,0);
    pTgraph_mpf_r[j]->Fit(pol1[j],"R");
    line->SetLineStyle(2);
    line->Draw("SAME");
    // //Divide by value@alpha cut
    // double y_al_cut = pol1[j]->GetParameter(0)+pol1[j]->GetParameter(1)*al_cut;
    //    cout<<"y_al_cut = "<<y_al_cut<<endl;
    // fill the output.dat file
    if (fp_mpf_r!=NULL) {
      Float_t value = pol1[j]->GetParameter(0);
      //   Float_t value = pol1[j]->GetParameter(0)/y_al_cut;
      Float_t uncert = pol1[j]->GetParError(0);
      fprintf(fp_mpf_r, "%f %f\n",value,uncert);
    }
    plotkfsr->SetBinContent(j+1,pol1[j]->GetParameter(0));
    //plotkfsr->SetBinContent(j+1,pol1[j]->GetParameter(0)/y_al_cut);
    plotkfsr->SetBinError(j+1,pol1[j]->GetParError(0));
    kFSR_MPF->SetBinContent(j+1,pol1[j]->GetParameter(0));
    // kFSR_MPF->SetBinContent(j+1,pol1[j]->GetParameter(0)/y_al_cut);
    kFSR_MPF->SetBinError(j+1,pol1[j]->GetParError(0));
    leg1->SetHeader("MPF, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
    leg1->Draw();
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    tex->DrawLatex(0.64,0.91,"2.11fb^{-1} (13TeV)");
    TString chi2_loglin = "#chi^{2}/n.d.f = ";
    chi2_loglin += trunc(pol1[j]->GetChisquare());
    chi2_loglin +="/";
    chi2_loglin +=trunc(pol1[j]->GetNDF());
    TLatex *tex2 = new TLatex();
    tex2->SetNDC();
    tex2->SetTextSize(0.035); 
    tex2->DrawLatex(0.64,0.35,chi2_loglin);

    //save the plots
      b[j]->Print(path+"plots/kFSR_MPF_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
  }
  fclose(fp_mpf_r);


  // create output file including the kFSR plot
    TFile* outputfile_mpf_r = new TFile(path+"Histo_KFSR_MPF_L1.root","RECREATE");
    kFSR_MPF->Write();
    outputfile_mpf_r->Write();
    outputfile_mpf_r->Close();

}
