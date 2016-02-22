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
#include "UsefulFunctions.C"
using namespace std;

void kFSR_pT(bool mpfMethod, TString path, TFile* datafile, TFile* MCfile){
  TStyle* m_gStyle = new TStyle();

  m_gStyle->SetOptFit(0);
  //datafile->Print();
  //MCfile->Print();

  // get ratio for MC to DATA responses
  double ratio_al[n_pt-1][n_eta-1][n_alpha]; //ratio at pt,eta,alpha bins
  double err_ratio_al[n_pt-1][n_eta-1][n_alpha]; //error of ratio at pt,eta,alpha bins
  for(int i=0; i<n_alpha; i++){
    for(int j=0; j<n_eta-1; j++){
      for(int k=0; k<n_pt-1; k++){
	ratio_al[k][j][i] = 0;
	err_ratio_al[k][j][i] = 0;
      }
    }
  }

  for(int i=0; i<n_alpha; i++){
    for(int j=0; j<n_eta-1; j++){
      for(int k=0; k<n_pt-1; k++){
	std::cout<<"For alpha cut at "<<alpha_bins[i]<<" eta range = "<<eta_bins[j]<<", "<<eta_bins[j+1]<<" and pT range = "<<pt_bins[k]<<", "<<pt_bins[k+1]<<std::endl;
	pair<double,double> res_mc =  Response(MCfile,alpha_bins[i],eta_bins[j],eta_bins[j+1],pt_bins[k],pt_bins[k+1],mpfMethod);
	pair<double,double> res_data =  Response(datafile,alpha_bins[i],eta_bins[j],eta_bins[j+1],pt_bins[k],pt_bins[k+1],mpfMethod);
	pair<double,double> ratio_res = Rmc_to_Rdata(res_mc,res_data);
	ratio_al[k][j][i] = ratio_res.first;
	err_ratio_al[k][j][i] = ratio_res.second;
      }
    }
  }

  //Divide reponses by value at alpha=0.3
  //find bin with alpha = 0.3
  int al_ref=0;
  for(int i=0; i<n_alpha; i++){
    if(fabs(alpha_bins[i]-0.3)<1e-4) al_ref=i;
  }
  cout<<"alpha=0.3 for bin#"<<al_ref<<endl;

  for(int j=0; j<n_eta-1; j++){
    for(int k=0; k<n_pt-1; k++){
      double norm_al02 = ratio_al[k][j][al_ref];
      double err_norm_al02 = err_ratio_al[k][j][al_ref];
      for(int i=0; i<n_alpha; i++){
  	if(norm_al02>0){
	  ratio_al[k][j][i] =   ratio_al[k][j][i]/norm_al02;
	  //	  cout<<"ratio_al["<<k<<"]["<<j<<"]["<<i<<"] = "<<ratio_al[k][j][i]<<endl;
	  err_ratio_al[k][j][i] = sqrt(abs(pow(err_ratio_al[k][j][i],2)-pow(err_norm_al02,2)));
	}
      }
    }
  }
  
  TLegend *leg1;
  leg1 = new TLegend(0.53,0.54,0.65,0.88,"","brNDC");//x+0.1
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.035);
  leg1->SetFillColor(10);
  leg1->SetLineColor(1);
  leg1->SetTextFont(42);

  TGraphErrors *graph[n_pt-1][n_eta-1];//set of points vs alpha
  TMultiGraph *pTgraph[n_eta-1];//set of different pT bins in on eta bin
  double xbin_tgraph[n_alpha],zero[n_alpha];
  for(int i=0;i<n_alpha;i++){
    xbin_tgraph[i] = alpha_bins[i];
    zero[i] = 0;
  }
  for(int j=0; j<n_eta-1; j++){
    pTgraph[j] = new TMultiGraph();
    for(int k=0; k<n_pt-1; k++){
      graph[k][j] = new TGraphErrors(n_alpha,xbin_tgraph,ratio_al[k][j],zero,err_ratio_al[k][j]);
      graph[k][j] = (TGraphErrors*)CleanEmptyPoints(graph[k][j]);
      //      graph[k][j]->SetMarkerStyle(20);
      // graph[k][j]->SetMarkerColor(kBlue-3+k);
      // graph[k][j]->SetLineColor(kBlue-3+k);
      graph[k][j]->SetMarkerSize(1.3);
      graph[k][j]->SetMarkerStyle(20+k);
      if(k<4){ //skip yellow color
      graph[k][j]->SetMarkerColor(k+1);
      graph[k][j]->SetLineColor(k+1);
      }
      else{
	graph[k][j]->SetMarkerColor(k+2);
	graph[k][j]->SetLineColor(k+2);
      }
      if(graph[k][j]->GetN()>0) 
	pTgraph[j]->Add(graph[k][j]);
      TString pTbin_label = "";
      pTbin_label+=pt_bins[k];
      pTbin_label+=" < p_{T} < ";
      pTbin_label+=pt_bins[k+1];
      //      cout<<j<<" "<<k<<" "<<pTbin_label<<endl;
      if(j==0) leg1->AddEntry(graph[k][j],pTbin_label,"epl");
    }
  }

// //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // // Plots

  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(alpha_bins[0]-0.01,1,alpha_bins[n_alpha-1]+0.01,1);

 
  // create output .dat file, including the kFSR extrapolation (alpha->0)
  FILE *fp; TH1D* kFSR_MPF; TH1D* kFSR_DiJet;
  if(mpfMethod){
    fp = fopen("output/KFSR_MPF_extrapolation.dat","w");
    kFSR_MPF = new TH1D("kfsr_mpf","kfsr_mpf", n_eta-1,eta_bins);
  }
  else{
    fp = fopen("output/KFSR_DiJet_extrapolation.dat","w");
    kFSR_DiJet = new TH1D("kfsr_dijet","kfsr_dijet", n_eta-1,eta_bins);
  }

  TH1D* plotkfsr = new TH1D("kfsr","kfsr", n_eta-1,eta_bins);


  


  //create plots
  TCanvas* a[n_eta-1];
  TString plotname[n_eta-1];
  TF1 *pol1[n_eta-1];
  for (int j=0; j<n_eta-1; j++){
  //  for (int j=0; j<1; j++){//TEST
    if(mpfMethod){
      plotname[j]="mpf_kfsr_diffPt_eta_"+eta_range[j]+"_"+eta_range[j+1];
    }
    else{
      plotname[j]="dijet_kfsr_diffPt_eta_"+eta_range[j]+"_"+eta_range[j+1];
    }
    a[j] = new TCanvas(plotname[j], plotname[j], 800,700);
    m_gStyle->SetOptTitle(0);
    pTgraph[j]->Draw("AP");
    //    pTgraph[j]->GetYaxis()->SetRangeUser(0.8,1.4);
    //    pTgraph[j]->GetYaxis()->SetRangeUser(0.9,1.1);
    pTgraph[j]->GetYaxis()->SetRangeUser(0.92,1.08);
    pTgraph[j]->GetYaxis()->SetTitle("(R_{MC}/R_{DATA})/(R_{MC}/R_{DATA})_{#alpha<0.2}");
    pTgraph[j]->GetXaxis()->SetTitle("cut on #alpha");

    pol1[j] = new TF1("pol1","pol1",0.14,0.46);  //TEST
    pol1[j]->SetParameters(0,0);
    pTgraph[j]->Fit(pol1[j],"R");
    line->SetLineStyle(2);
    line->Draw("SAME");
    // fill the output.dat file
    if (fp!=NULL) {
      Float_t value = pol1[j]->GetParameter(0);
      Float_t uncert = pol1[j]->GetParError(0);
      fprintf(fp, "%f %f\n",value,uncert);
    }
    plotkfsr->SetBinContent(j+1,pol1[j]->GetParameter(0));
    plotkfsr->SetBinError(j+1,pol1[j]->GetParError(0));
    if(mpfMethod){
      kFSR_MPF->SetBinContent(j+1,pol1[j]->GetParameter(0));
      kFSR_MPF->SetBinError(j+1,pol1[j]->GetParError(0));
    }
    else{
      kFSR_DiJet->SetBinContent(j+1,pol1[j]->GetParameter(0));
      kFSR_DiJet->SetBinError(j+1,pol1[j]->GetParError(0));
    }
   
    if(mpfMethod){
      leg1->SetHeader("MPF, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
    }
    else{
      leg1->SetHeader("p_{T} balance, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
    }
    // leg1->AddEntry(graph1[j], "R(MC)/R(DATA)","P");
    // leg1->AddEntry(pol1, "linear fit","L");
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
    if(mpfMethod){
      a[j]->Print("plots/kFSR_MPF_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
    }
    else{
      a[j]->Print("plots/kFSR_Pt_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
    }

  }
  fclose(fp);


  // create output file including the kFSR plot
  if(mpfMethod){
    TFile* outputfile = new TFile(path+"Histo_KFSR_MPF_L1.root","RECREATE");
    kFSR_MPF->Write();
    outputfile->Write();
    outputfile->Close();
  }
  else{
    TFile* outputfile = new TFile(path+"Histo_KFSR_DiJet_L1.root","RECREATE");
    kFSR_DiJet->Write();
    outputfile->Write();
    outputfile->Close();
  }


}
