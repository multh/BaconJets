// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// with this script pT extrapolations in bins of eta will be calculated
//
// kFSR is aproximated by function and kFSR correction is done by the 
// function values at each eta bin
//
// number of root files with historgrams and txt files for JEC db
// are produced
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


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
#include "TMath.h"
#include "TPaveStats.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"
#include "TVirtualFitter.h"
using namespace std;

void PTextrapolation_TTree_kFSRfit(bool mpfMethod, TString path, TFile* datafile, TFile* MCfile, TString txttag, TString jettag, TString variation, TString tag, double al_cut=0.2, int nResponseBins=100){
  gStyle->SetOptFit(000);


  // Tag for time dependence plots
  //TString tag = "";

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




  // systematic uncertainty: "down" means that 0.5*ptave will be used in the loglinear fit. "up" means that 2*ptave will be used in the loglinear fit
  int syst = 0;
  if(variation=="central") syst=0;
  if(variation=="down") syst=1;
  if(variation=="up") syst=2;
  if(variation=="doubleup") syst=3;
  if(variation=="nominal") syst=5;


  // get the histos for pt average
  TH1D* ptave_data[n_eta-1];
  int countPt = 0;
  TString namehist = "ptave_";
  for(int i=0; i<n_eta-1; i++){
    TString selection = "alpha<";
    //    TString selection = "alpha_sum<";
    selection+=al_cut;
    selection+=" && fabs(probejet_eta)<";
    selection+=eta_range[i+1];
    selection+=" && fabs(probejet_eta)>=";
    selection+=eta_range[i];
    TString var1 = "pt_ave";   
    ptave_data[i] = (TH1D*)GetHist(datafile, selection, var1, 300,0,3000)->Clone();
    TString namecur = namehist;namecur+=countPt;
    ptave_data[i]->SetName(namecur);
    countPt++;
  }




  // get ratio for MC to DATA responses
  double ratio_al_rel_r[n_pt-1][n_eta-1]; //ratio at pt,eta,alpha bins
  double err_ratio_al_rel_r[n_pt-1][n_eta-1]; //error of ratio at pt,eta,alpha bins
  double ratio_al_mpf_r[n_pt-1][n_eta-1]; //ratio at pt,eta,alpha bins
  double err_ratio_al_mpf_r[n_pt-1][n_eta-1]; //error of ratio at pt,eta,alpha bins
  TH1D *hdata_rel_r[n_pt-1][n_eta-1];// pT-balance responce for data
  TH1D *hdata_mpf_r[n_pt-1][n_eta-1];//MPF responce for data
  TH1D *hmc_rel_r[n_pt-1][n_eta-1];// pT-balance responce for MC
  TH1D *hmc_mpf_r[n_pt-1][n_eta-1];//MPF responce for MC
  int count = 0;
  TString name1 = "hist_data_rel_r_";
  TString name2 = "hist_data_mpf_r_";
  TString name3 = "hist_mc_rel_r_";
  TString name4 = "hist_mc_mpf_r_";
    for(int j=0; j<n_eta-1; j++){
      for(int k=0; k<n_pt-1; k++){
	ratio_al_rel_r[k][j] = 0;
	err_ratio_al_rel_r[k][j] = 0;
	ratio_al_mpf_r[k][j] = 0;
	err_ratio_al_mpf_r[k][j] = 0;
	TString name = name1; name+=count;
	hdata_rel_r[k][j] = new TH1D(name,"",nResponseBins, 0, 2.5);
	name = name2;name+=count;
	hdata_mpf_r[k][j] = new TH1D(name,"",nResponseBins, 0, 2.5);
	name = name3; name+=count;
	hmc_rel_r[k][j] = new TH1D(name,"",nResponseBins, 0, 2.5);
	name = name4; name+=count;
	hmc_mpf_r[k][j] = new TH1D(name,"",nResponseBins, 0, 2.5);
	count++;
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
     if(*alpha_data>al_cut) continue;
     for(int k=0; k<n_pt-1; k++){
   	   if(*pt_ave_data<pt_bins[k] || *pt_ave_data>pt_bins[k+1]) continue;
	   for(int j=0; j<n_eta-1; j++){
	     if(fabs(*probejet_eta_data)>eta_bins[j+1] || fabs(*probejet_eta_data)<eta_bins[j]) continue;
	     //  if(*probejet_eta_data>eta_bins[j+1] || *probejet_eta_data<eta_bins[j]) continue;//TEST 
	     else{
	       hdata_rel_r[k][j]->Fill(*rel_r_data,*weight_data);
	       hdata_mpf_r[k][j]->Fill(*mpf_r_data,*weight_data);
	     }
	   }
     }
   }

   TTreeReader myReader_MC("AnalysisTree", MCfile);
   TTreeReaderValue<Float_t> pt_ave_mc(myReader_MC, "pt_ave");
   TTreeReaderValue<Float_t> probejet_eta_mc(myReader_MC, "probejet_eta");
   TTreeReaderValue<Float_t> alpha_mc(myReader_MC, "alpha");
   //   TTreeReaderValue<Float_t> alpha_mc(myReader_MC, "alpha_sum"); //TEST alpha_sum
   TTreeReaderValue<Float_t> rel_r_mc(myReader_MC, "rel_r");
   TTreeReaderValue<Float_t> mpf_r_mc(myReader_MC, "mpf_r");
   TTreeReaderValue<Float_t> weight_mc(myReader_MC, "weight");
   while (myReader_MC.Next()) {
     if(*alpha_mc>al_cut) continue;
     for(int k=0; k<n_pt-1; k++){
       if(*pt_ave_mc<pt_bins[k] || *pt_ave_mc>pt_bins[k+1]) continue;
       for(int j=0; j<n_eta-1; j++){
	 if(fabs(*probejet_eta_mc)>eta_bins[j+1] || fabs(*probejet_eta_mc)<eta_bins[j]) continue;
	 //	 if(*probejet_eta_mc>eta_bins[j+1] || *probejet_eta_mc<eta_bins[j]) continue; //TEST
   	   else{
   	     hmc_rel_r[k][j]->Fill(*rel_r_mc,*weight_mc);
   	     hmc_mpf_r[k][j]->Fill(*mpf_r_mc,*weight_mc);
   	   }
       }
     }
   }


     for(int j=0; j<n_eta-1; j++){
      for(int k=0; k<n_pt-1; k++){
	std::cout<<"Ratio: For eta range = "<<eta_bins[j]<<", "<<eta_bins[j+1]<<" and pT range = "<<pt_bins[k]<<", "<<pt_bins[k+1]<<std::endl;
	pair<double,double> res_mc_rel_r,res_data_rel_r;
	pair<double,double> res_mc_mpf_r,res_data_mpf_r;
	res_mc_rel_r = GetValueAndError(hmc_rel_r[k][j]);
	res_data_rel_r = GetValueAndError(hdata_rel_r[k][j]);
	res_mc_mpf_r = GetValueAndError(hmc_mpf_r[k][j]);
	res_data_mpf_r = GetValueAndError(hdata_mpf_r[k][j]);


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
	ratio_al_rel_r[k][j] = ratio_res_rel_r.first;
	err_ratio_al_rel_r[k][j] = ratio_res_rel_r.second;
	ratio_al_mpf_r[k][j] = ratio_res_mpf_r.first;
	err_ratio_al_mpf_r[k][j] = ratio_res_mpf_r.second;
      }
    }
 




  // get ratio for MC to DATA responses
  double ratio_mpf[n_eta-1][n_pt-1]; //ratio at pt,eta bins for alpha = 0.3
  double err_ratio_mpf[n_eta-1][n_pt-1]; //error of ratio at pt,eta bins for alpha = 0.3
  for(int j=0; j<n_eta-1; j++){
    for(int k=0; k<n_pt-1; k++){
      if(mpfMethod){
	ratio_mpf[j][k] = ratio_al_mpf_r[k][j];
	err_ratio_mpf[j][k] = err_ratio_al_mpf_r[k][j];
      }
      else{
	ratio_mpf[j][k] = ratio_al_rel_r[k][j];
	err_ratio_mpf[j][k] = err_ratio_al_rel_r[k][j];
      }
    }
  }

  // create and fill tgrapherrors
  double xbin_tgraph[n_pt-1];
  double zero[n_pt-1];
  for(int i=0;i<n_pt-1;i++){
    xbin_tgraph[i]=(pt_bins[i]+pt_bins[i+1])/2;
    zero[i]=(pt_bins[i+1]-pt_bins[i])/2 ;
  }
  TGraphErrors *graph1_mpf[n_eta-1];
  TGraphErrors *graph1_dijet[n_eta-1];
 
 
  for(int j=0; j<n_eta-1; j++){
    graph1_mpf[j] = new TGraphErrors(n_pt-1, xbin_tgraph, ratio_mpf[j], zero, err_ratio_mpf[j]);
    graph1_mpf[j] = (TGraphErrors*)CleanEmptyPoints(graph1_mpf[j]);
  }

  

  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,pt_bins[n_pt-1]+10,1);

  // create a function for the loglinear fit
  TF1 * f1[n_eta-1];
  // create a function for the constant fit
  TF1 * f2[n_eta-1];

  // create output .dat file, including the kFSR extrapolation (alpha->0)
  FILE *fp; FILE *fp2; FILE *l2resfile;
  if(mpfMethod){
    fp = fopen(path+"output/pT_MPF_extrapolations.dat","w");
    fp2 = fopen(path+"output/pT_MPF_constantExtrapolation.dat","w");
    l2resfile = fopen(path+"output/L2Res_MPF.dat","w");
  }
  else{
    fp = fopen(path+"output/pT_DiJet_extrapolations.dat","w");
    fp2 = fopen(path+"output/pT_DiJet_constantExtrapolation.dat","w");
    l2resfile = fopen(path+"output/L2Res_DiJet.dat","w");
  }


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // Plots

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  TCanvas* asd[n_eta-1];
  TString plotname[n_eta-1];
  double Vcov[3][n_eta-1];//covarance matrix for log lin fit results

  for (int j=0; j<n_eta-1; j++){
    if(mpfMethod){
      plotname[j]="mpf_ptextra_eta_"+eta_range[j]+"_"+eta_range[j+1];
    }
    else{
      plotname[j]="dijet_ptextra_eta_"+eta_range[j]+"_"+eta_range[j+1];
    }
    asd[j] = new TCanvas(plotname[j], plotname[j], 800,700);
    gStyle->SetOptTitle(0);
    gPad->SetLogx();
    graph1_mpf[j]->SetMarkerColor(kBlue);
    graph1_mpf[j]->SetMarkerStyle(20);
    graph1_mpf[j]->SetLineColor(kBlue);
   
 
    f1[j] = new TF1(plotname[j]+"f1","[0]+[1]*TMath::Log(x)", 50 , pt_bins[n_pt-1]+10);
    //f1[j] = new TF1(plotname[j]+"f1","[0]+[1]*TMath::Log(x)", 70 , pt_bins[n_pt-1]+10);
    //    f1[j] = new TF1(plotname[j]+"f1","[0]+[1]*TMath::Log(x)", 70 , 570);
    f1[j]->SetParameters(1,0);
    f2[j] = new TF1(plotname[j]+"f2","pol0", 50 , pt_bins[n_pt-1]+10);
    //    f2[j]->SetParameters(1);
    //    f2[j] = new TF1(plotname[j]+"f2","pol0", 70 , pt_bins[n_pt-1]+10);
    f2[j]->SetLineColor(kBlue);
    f2[j]->SetLineStyle(3);

    TFitResultPtr fitloglin =  graph1_mpf[j]->Fit(plotname[j]+"f1","SRM");
    TMatrixDSym cov = fitloglin->GetCovarianceMatrix();
    Vcov[0][j] = cov(0,0);     Vcov[1][j] = cov(1,1);     Vcov[2][j] = cov(0,1);
    graph1_mpf[j]->Fit(plotname[j]+"f2","SRM + SAME");
    graph1_mpf[j]->Draw("AP");
    graph1_mpf[j]->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph1_mpf[j]->GetXaxis()->SetTitleSize(0.05);
    graph1_mpf[j]->GetXaxis()->SetTitleOffset(0.80);
    graph1_mpf[j]->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    //    graph1_mpf[j]->GetYaxis()->SetRangeUser(0.9,1.15);
    //   graph1_mpf[j]->GetYaxis()->SetRangeUser(0.9,1.25);
   graph1_mpf[j]->GetYaxis()->SetRangeUser(0.85,1.25);
 
    line->SetLineStyle(2);
    line->Draw("SAME");

    if (fp2!=NULL) {
      // getting the p0 parameter from the constant fit
      Float_t value = f2[j]->GetParameter(0);
      fprintf(fp2, "%f\n",value);
    }
    if (l2resfile!=NULL) {
      fprintf(l2resfile, "%f\n", eta_bins[j]);
    }


    asd[j]->Modified(); asd[j]->Update();
    TPaveStats *st = ((TPaveStats*)(graph1_mpf[j]->GetListOfFunctions()->FindObject("stats")));
    if (st) {
      st->SetTextColor(kRed);
      st->SetX1NDC(0.69); st->SetX2NDC(0.99);
      st->SetY1NDC(0.65); st->SetY2NDC(0.8);
    }
    st = ((TPaveStats*)(graph1_mpf[j]->GetListOfFunctions()->FindObject("stats")));
    if (st) {
      st->SetTextColor(graph1_mpf[j]->GetLineColor());
      st->SetX1NDC(0.69); st->SetX2NDC(0.99);
      st->SetY1NDC(0.85); st->SetY2NDC(0.95);
    }
    asd[j]->Modified(); asd[j]->Update();

    TLegend *leg1;
    leg1 = new TLegend(0.12,0.68,0.35,0.88,"","brNDC");//x+0.1
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.038);
    leg1->SetFillColor(10);
    leg1->SetLineColor(1);
    leg1->SetTextFont(42);
    if(mpfMethod){
      leg1->SetHeader("MPF #bar{p}_{T} extrapolation, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
    }
    else{
      leg1->SetHeader("p_{T} balance #bar{p}_{T} extrapolation, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
    }

    TString alVal;
    alVal.Form("%0.2f\n",al_cut);
    TString altitle = "{#alpha<"+alVal+"}";
    TString axistitle = "(R^{MC}/R^{data})_";
    axistitle +=altitle;
    leg1->AddEntry(graph1_mpf[j], axistitle,"P");
    leg1->AddEntry(f1[j], "loglinear fit","L");
    leg1->AddEntry(f2[j], "constant fit","L");
    leg1->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    //    tex->DrawLatex(0.64,0.91,"2.11fb^{-1} (13TeV)");
    tex->DrawLatex(0.64,0.91,"0.804fb^{-1} (13TeV)");
    //    tex->DrawLatex(0.64,0.91,"0.804fb^{-1} (13TeV)");

    TString chi2_loglin = "loglinear fit #chi^{2}/n.d.f = ";
    chi2_loglin += trunc(f1[j]->GetChisquare());
    chi2_loglin +="/";
    chi2_loglin +=trunc(f1[j]->GetNDF());
    TString chi2_const = "constant fit #chi^{2}/n.d.f = ";
    chi2_const+=trunc(f2[j]->GetChisquare());
    chi2_const+="/";
    chi2_const+=trunc(f2[j]->GetNDF());

    TLatex *tex2 = new TLatex();
    tex2->SetNDC();
    tex2->SetTextSize(0.035); 
    tex2->DrawLatex(0.51,0.73,chi2_loglin);
    tex2->DrawLatex(0.51,0.69,chi2_const);

    //save plots
    if(mpfMethod){
      asd[j]->Print(path+"plots/pTextrapolation_MPF_pT_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
    }
    else{
      asd[j]->Print(path+"plots/pTextrapolation_Pt_pT_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
    }
  }
  fclose(fp);
  fclose(fp2);


  // //==============================================================================================================================
  // // some couts for cross check
  // for(int i=0; i<n_eta-1; i++){
  //   cout << "mean value, eta " << eta_range[i] << " to " << eta_range[i+1] << ": " << ptave_data[i]->GetMean() << endl;
  //   //cout << "max value, eta " << eta_range[i] << " to " << eta_range[i+1] << ": " << ptave_data[i]->FindLastBinAbove() << endl;
  // }
  // for(int i=0; i<n_eta-1; i++){
  //   cout << "loglin fit value " << eta_range[i] << " to " << eta_range[i+1] << ": " << f1[i]->GetParameter(1) << endl;
  // }
  // for(int i=0; i<n_eta-1; i++){
  //   cout << "max value of pt " << eta_range[i] << " to " << eta_range[i+1] << ": " << ptave_data[i]->FindLastBinAbove(0.) << endl;
  // }

  // for(int i=0; i<n_eta-1; i++){
  //   for(int j=0; j<graph1_mpf[i]->GetN(); j++) {
  //     cout << "uncert, eta " << eta_range[i] << " to " << eta_range[i+1] << ", pt range " << pt_range[j] << " to " << pt_range[j+1] << ": "  << graph1_mpf[i]->GetEY()[j] << endl;
  //   }
  // }


  // //==============================================================================================================================
  // get the kFSR plots and calculate residuals
  if(mpfMethod){
    TFile* kfsr_mpf = new TFile(path+"Histo_KFSR_MPF_L1.root","READ");
    TH1D* hist_kfsr_mpf = (TH1D*)kfsr_mpf->Get("kfsr_mpf");
    //fit kFSR
    TF1 *kfsr_fit_mpf = new TF1("kfsr_fit_mpf","[0]+([1]*TMath::CosH(x))/(1+[2]*TMath::CosH(x))",0,5.);                                          
    kfsr_fit_mpf->SetParameters(0,0,100.);
    kfsr_fit_mpf->SetLineColor(kRed+1);
    hist_kfsr_mpf->Fit("kfsr_fit_mpf","SR");
    //Create a histogram to hold the confidence intervals                                                                                        
    TH1D *hist_kfsr_fit_mpf = (TH1D*)hist_kfsr_mpf->Clone();
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hist_kfsr_fit_mpf);
    //Now the "hist_kfsr_fit" histogram has the fitted function values as the
    //bin contents and the confidence intervals as bin errors
    hist_kfsr_fit_mpf->SetStats(kFALSE);
    hist_kfsr_fit_mpf->SetFillColor(kRed-10);     
    hist_kfsr_fit_mpf->SetName("hist_kfsr_fit_mpf");
    hist_kfsr_fit_mpf->SetTitle("kfsr fit for mpf");
    double flat[n_eta-1];
    double loglin[n_eta-1];
    for (int j=0; j<n_eta-1; j++){
      flat[j] = hist_kfsr_fit_mpf->GetBinContent(j+1)*f2[j]->GetParameter(0);
      loglin[j] = hist_kfsr_fit_mpf->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()) );
    }
    double flat_norm=0, loglin_norm=0;
    for (int j=0; j<n_etabarr; j++){
      flat_norm += flat[j];
      loglin_norm += loglin[j];
    }
    flat_norm = flat_norm/n_etabarr;
    loglin_norm = loglin_norm/n_etabarr;
    for (int j=0; j<n_eta-1; j++){
      flat[j] = flat[j]/flat_norm;
      loglin[j] = loglin[j]/loglin_norm;
      //cout << flat[j] << endl;
    }

    cout << "normalization factor (flat): " << flat_norm << endl;
    cout << "normalization factor (logl): " << loglin_norm << endl;

    

    TH1D* Residual_logpt_MPF = new TH1D("res_logpt_mpf","res_logpt_mpf", n_eta-1,eta_bins);
    TH1D* Residual_const_MPF = new TH1D("res_const_mpf","res_const_mpf", n_eta-1,eta_bins);
    TH1D* ptave_const_MPF = new TH1D("ptave_const_mpf","ptave_const_mpf", n_eta-1,eta_bins);
    TH1D* ptave_logpt_MPF = new TH1D("ptave_logpt_mpf","ptave_logpt_mpf", n_eta-1,eta_bins);
    // for (int j=n_eta-1; j>0; --j){
    //   if(syst==5){
    // 	f1[j]->Print();
    // 	ptave_logpt_MPF->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()));
    // 	ptave_logpt_MPF->SetBinError(j+1,Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(ptave_data[j]->GetMean()),2)+2*Vcov[2][j]*TMath::Log(ptave_data[j]->GetMean()));
    //   }
    //   ptave_const_MPF->SetBinContent(j+1,f2[j]->GetParameter(0));
    //   ptave_const_MPF->SetBinError(j+1,f2[j]->GetParError(0));                                                                                 
     
    //   //   ptave_logpt_MPF->SetBinError(j+1,f2[j]->GetParError(0));   
    // }
    ofstream output, output_loglin, uncerts, uncerts_loglin;
    output.open(path+"output/Spring16_25ns_MPF_FLAT_L2Residual_"+txttag+"_"+jettag+".txt");
    output_loglin.open(path+"output/Spring16_25ns_MPF_LOGLIN_L2Residual_"+txttag+"_"+jettag+".txt");
    uncerts.open(path+"output/Spring16_25ns_MPF_FLAT_L2Residual_"+txttag+"_"+jettag+".txt.STAT");
    uncerts_loglin.open(path+"output/Spring16_25ns_MPF_LOGLIN_L2Residual_"+txttag+"_"+jettag+".txt.STAT");
    output  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(208))/208))) Correction L2Relative}" << endl;
    output_loglin  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(208))/208))) Correction L2Relative}" << endl;
    uncerts << "{ 1 JetEta 1 JetPt [0] kFSR_err f0_err}" << endl;
    uncerts_loglin << "{ 1 JetEta 1 JetPt [0] kFSR_err cov(0,0) cov(1,1) cov(0,1) }" << endl;
    for (int j=n_eta-1; j>0; --j){
      output << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   55   " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "   " << 1/flat_norm << " " << hist_kfsr_fit_mpf->GetBinContent(j) << "   " << f2[j-1]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;
      output_loglin << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   55   " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "   " << 1/loglin_norm << " " << hist_kfsr_fit_mpf->GetBinContent(j) << "   " << f1[j-1]->GetParameter(0) << " " << f1[j-1]->GetParameter(1) << "   1 0.0000 0.0" << endl;
      // uncerts << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "       " << sqrt( pow(f2[j-1]->GetParameter(0)*hist_kfsr_fit_mpf->GetBinError(j),2)+pow( hist_kfsr_fit_mpf->GetBinContent(j)*f2[j-1]->GetParError(0),2)  ) / flat_norm  << endl;
      // uncerts_loglin << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "       " << sqrt(  pow(f1[j-1]->GetParameter(0)+f1[j-1]->GetParameter(1)*TMath::Log(ptave_data[j-1]->GetMean()),2)*pow(hist_kfsr_fit_mpf->GetBinError(j),2)  + pow(hist_kfsr_fit_mpf->GetBinContent(j),2)*(Vcov[0][j-1]+Vcov[1][j-1]*pow(TMath::Log(ptave_data[j-1]->GetMean()),2)+2*Vcov[2][j-1]*TMath::Log(ptave_data[j-1]->GetMean()))) / loglin_norm << endl;

      uncerts << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " 
	      << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "  " << hist_kfsr_fit_mpf->GetBinError(j-1) / flat_norm 
	      <<" "<< f2[j-1]->GetParError(0) << endl;

      uncerts_loglin << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" 
		     << eta_range[j-1] << "  3    55    " << ptave_data[j-1]->FindLastBinAbove(0.)*10 
		     << " " << hist_kfsr_fit_mpf->GetBinError(j-1) / loglin_norm  
		     << " " <<Vcov[0][j-1] <<" "<<Vcov[1][j-1]<<" "<<Vcov[2][j-1]<<endl; 

    }
  
    for (int j=0; j<n_eta-1; j++){
      output << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   55   " 
	     << ptave_data[j]->FindLastBinAbove(0.)*10 << "   " << 1/flat_norm << " " << hist_kfsr_fit_mpf->GetBinContent(j+1) 
	     << "   " << f2[j]->GetParameter(0) << " 0   1 0.0000 0.0"  << endl;
      output_loglin << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   55   " 
		    << ptave_data[j]->FindLastBinAbove(0.)*10 << "   " << 1/loglin_norm << " " 
		    << hist_kfsr_fit_mpf->GetBinContent(j+1) << "   " << f1[j]->GetParameter(0) 
		    << " " << f1[j]->GetParameter(1) << "   1 0.0000 0.0" << endl;
      // uncerts << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " << ptave_data[j]->FindLastBinAbove(0.)*10 << "       " << sqrt( pow(f2[j]->GetParameter(0)*hist_kfsr_fit_mpf->GetBinError(j+1),2)+pow( hist_kfsr_fit_mpf->GetBinContent(j+1)*f2[j]->GetParError(0),2)  ) / flat_norm  << endl;
      // uncerts_loglin << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " << ptave_data[j]->FindLastBinAbove(0.)*10 << "       " << sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()),2)*pow(hist_kfsr_fit_mpf->GetBinError(j+1),2)   + pow(hist_kfsr_fit_mpf->GetBinContent(j+1),2)*(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(ptave_data[j]->GetMean()),2)+2*Vcov[2][j]*TMath::Log(ptave_data[j]->GetMean()))) / loglin_norm << endl;

      uncerts << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " 
	      << ptave_data[j]->FindLastBinAbove(0.)*10 
	      <<" "<< hist_kfsr_fit_mpf->GetBinError(j)/flat_norm<< " " << f2[j]->GetParameter(0) <<endl;

      uncerts_loglin << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " 
		  << ptave_data[j]->FindLastBinAbove(0.)*10 << " " 
		  << hist_kfsr_fit_mpf->GetBinError(j)/loglin_norm<< " "<< Vcov[0][j]<<" "<< Vcov[1][j]<<" "<< Vcov[2][j] << endl; 

   double pave_value_corr;
      if(syst==5){
	pave_value_corr = ptave_data[j]->GetMean();/*
    	Residual_logpt_MPF->SetBinContent(j+1,hist_kfsr_fit_mpf->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()) ) );
    	Residual_logpt_MPF->SetBinError(j+1, sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()),2)*pow(hist_kfsr_fit_mpf->GetBinError(j+1),2)   + pow(hist_kfsr_fit_mpf->GetBinContent(j+1),2)*(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(ptave_data[j]->GetMean()),2)+2*Vcov[2][j]*TMath::Log(ptave_data[j]->GetMean()))));
	ptave_logpt_MPF->SetBinError(j+1,sqrt(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(ptave_data[j]->GetMean()),2)+2*Vcov[2][j]*TMath::Log(ptave_data[j]->GetMean())));
	ptave_logpt_MPF->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean())); */
      }
      if(syst==0){
	pave_value_corr = 120;
	/*	Residual_logpt_MPF->SetBinContent(j+1,hist_kfsr_fit_mpf->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(120) ) );
    	Residual_logpt_MPF->SetBinError(j+1, sqrt(pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(120),2)*pow(hist_kfsr_fit_mpf->GetBinError(j+1),2)   + pow(hist_kfsr_fit_mpf->GetBinContent(j+1),2)*(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(120),2)+2*Vcov[2][j]*TMath::Log(120))));
	ptave_logpt_MPF->SetBinError(j+1,sqrt(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(120.),2)+2*Vcov[2][j]*TMath::Log(120.)));
	ptave_logpt_MPF->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(120.));
	*/
      }
      if(syst==1){
	pave_value_corr = 60;
      }
      if(syst==2){
	pave_value_corr = 240;
	/*
    	Residual_logpt_MPF->SetBinContent(j+1,hist_kfsr_fit_mpf->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(240) ) );
    	Residual_logpt_MPF->SetBinError(j+1, sqrt(pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(240),2)*pow(hist_kfsr_fit_mpf->GetBinError(j+1),2)   + pow(hist_kfsr_fit_mpf->GetBinContent(j+1),2)*(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(240),2)+2*Vcov[2][j]*TMath::Log(240))));
	ptave_logpt_MPF->SetBinError(j+1,sqrt(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(120.),2)+2*Vcov[2][j]*TMath::Log(240.)));
	ptave_logpt_MPF->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(240.));      
	*/
}
      if(syst==3){
	pave_value_corr = 480;
	/*   
    	Residual_logpt_MPF->SetBinContent(j+1,hist_kfsr_fit_mpf->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(480) ) );
    	Residual_logpt_MPF->SetBinError(j+1, sqrt(pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(480),2)*pow(hist_kfsr_fit_mpf->GetBinError(j+1),2)   + pow(hist_kfsr_fit_mpf->GetBinContent(j+1),2)*(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(480),2)+2*Vcov[2][j]*TMath::Log(480))));
	ptave_logpt_MPF->SetBinError(j+1,sqrt(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(120.),2)+2*Vcov[2][j]*TMath::Log(120.)));
	ptave_logpt_MPF->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(120.));
	*/
}

    	Residual_logpt_MPF->SetBinContent(j+1,hist_kfsr_fit_mpf->GetBinContent(j+1)*(f1[j]->GetParameter(0)
										     +f1[j]->GetParameter(1)*TMath::Log(pave_value_corr)));
    	Residual_logpt_MPF->SetBinError(j+1, sqrt(pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(pave_value_corr),2)*pow(hist_kfsr_fit_mpf->GetBinError(j+1),2)   + pow(hist_kfsr_fit_mpf->GetBinContent(j+1),2)*(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(pave_value_corr),2)+2*Vcov[2][j]*TMath::Log(pave_value_corr))));
	ptave_logpt_MPF->SetBinError(j+1,sqrt(Vcov[0][j]
					      +Vcov[1][j]*pow(TMath::Log(pave_value_corr),2)+2*Vcov[2][j]*TMath::Log(pave_value_corr)));
	ptave_logpt_MPF->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(pave_value_corr));


      Residual_const_MPF->SetBinContent(j+1,hist_kfsr_fit_mpf->GetBinContent(j+1)*f2[j]->GetParameter(0));
      Residual_const_MPF->SetBinError(j+1, sqrt( pow(f2[j]->GetParameter(0)*hist_kfsr_fit_mpf->GetBinError(j+1),2)+pow( hist_kfsr_fit_mpf->GetBinContent(j+1)*f2[j]->GetParError(0),2)  )  );
      ptave_const_MPF->SetBinContent(j+1,f2[j]->GetParameter(0));
      ptave_const_MPF->SetBinError(j+1,f2[j]->GetParError(0));
     

    }




    TFile* outputfile;
    if(syst==5){
      outputfile = new TFile(path+"Histo_Res_MPF_L1_"+txttag+"_"+jettag+tag+".root","RECREATE");
    }
    if(syst<4){
     outputfile = new TFile(path+"Histo_Res_MPF_L1_"+txttag+"_"+jettag+"_"+variation+".root","RECREATE");
    }
    /*
    if(syst==1){
     outputfile = new TFile(path+"Histo_Res_MPF_L1_"+variation+".root","RECREATE");
    }
    if(syst==2){
      outputfile = new TFile(path+"Histo_Res_MPF_L1_"+variation+".root","RECREATE");
    }
    if(syst==3){
      outputfile = new TFile(path+"Histo_Res_MPF_L1_"+variation+".root","RECREATE");
    }
    */
    hist_kfsr_mpf->Write();
    hist_kfsr_fit_mpf->Write();                                                                                                                
    kfsr_fit_mpf->Write();
    ptave_const_MPF->Write();                                                                                                                    
    ptave_logpt_MPF->Write();
    /*    for (int j=0; j<n_eta-1; j++){                                                                                                               
      f1[j]->Write();                                                                                                                            
      f2[j]->Write();                                                                                                                            
      } */      
    Residual_logpt_MPF->Write();
    Residual_const_MPF->Write();
    outputfile->Write();
    outputfile->Close();
    // TFile* outputfile2 = new TFile(path+"Histo_ptave_MPF_L1.root","RECREATE");

    TFile* outputfile2;
    if(syst==5){
      outputfile2 = new TFile(path+"Histo_ptave_MPF_L1_"+txttag+"_"+jettag+tag+".root","RECREATE");
    }
    if(syst<4){
     outputfile2 = new TFile(path+"Histo_ptave_MPF_L1_"+txttag+"_"+jettag+"_"+variation+".root","RECREATE");
    }
    /*
    if(syst==1){
      outputfile2 = new TFile(path+"Histo_ptave_MPF_L1_"+variation+".root","RECREATE");
    }
    if(syst==2){
      outputfile2 = new TFile(path+"Histo_ptave_MPF_L1_"+variation+".root","RECREATE");
    }
    if(syst==3){
      outputfile2 = new TFile(path+"Histo_ptave_MPF_L1_"+variation+".root","RECREATE");
    }
    */
    ptave_const_MPF->Write();
    ptave_logpt_MPF->Write();
  for (int j=0; j<n_eta-1; j++){
    f1[j]->Write();
    f2[j]->Write();
  }
    outputfile2->Write();
    outputfile2->Close();
  }



  //==============================================================================================================================
  // DIJET balance
  else{
    TFile* kfsr_dijet = new TFile(path+"Histo_KFSR_DiJet_L1.root","READ");
    TH1D* hist_kfsr_dijet = (TH1D*)kfsr_dijet->Get("kfsr_dijet");
    //suda fit kFSR
    TF1 *kfsr_fit_dijet = new TF1("kfsr_fit_dijet","[0]+([1]*TMath::CosH(x))/(1+[2]*TMath::CosH(x))",0,5.); 
    kfsr_fit_dijet->SetParameters(0,0,200.); //Almost everything, but not AK4PUPPI herwigpp
    //    kfsr_fit_dijet->SetParameters(1.,1.,100.); //AK4CHS, AK4PUPPI herwigpp
    //    kfsr_fit_dijet->SetParameters(1.,-200.,100.); //AK4CHS, AK8PUPPI pythia8 
    //    kfsr_fit_dijet->SetParameters(0,0,100.); //AK4PUPPI

    //    kfsr_fit_dijet->SetParameters(10.,-1000.,200.); //AK8CHS herwigpp, AK8PUPPI herwigpp 
    kfsr_fit_dijet->SetLineColor(kBlue+1);
    hist_kfsr_dijet->Fit("kfsr_fit_dijet","SR");
    //Create a histogram to hold the confidence intervals                                                                                        
    TH1D *hist_kfsr_fit_dijet = (TH1D*)hist_kfsr_dijet->Clone();
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hist_kfsr_fit_dijet);
    //Now the "hist_kfsr_fit" histogram has the fitted function values as the
    //bin contents and the confidence intervals as bin errors
    hist_kfsr_fit_dijet->SetStats(kFALSE);
    hist_kfsr_fit_dijet->SetFillColor(kBlue-10);     
    hist_kfsr_fit_dijet->SetName("hist_kfsr_fit_dijet");
    hist_kfsr_fit_dijet->SetTitle("kfsr fit for dijet");
    double flat[n_eta-1];
    double loglin[n_eta-1];
    for (int j=0; j<n_eta-1; j++){
      flat[j] = hist_kfsr_fit_dijet->GetBinContent(j+1)*f2[j]->GetParameter(0);
      loglin[j] = hist_kfsr_fit_dijet->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()) );
    }
    double flat_norm = 0, loglin_norm = 0;
    for (int j=0; j<n_etabarr; j++){
      flat_norm += flat[j];
      loglin_norm += loglin[j];
    }
    flat_norm = flat_norm/n_etabarr;
    loglin_norm = loglin_norm/n_etabarr;
    for (int j=0; j<n_eta-1; j++){
      flat[j] = flat[j]/flat_norm;
      loglin[j] = loglin[j]/loglin_norm;
      //cout << flat[j] << endl;
    }

    cout << "normalization factor (flat): " << flat_norm << endl;
    cout << "normalization factor (logl): " << loglin_norm << endl;


    TH1D* Residual_logpt_DiJet = new TH1D("res_logpt_dijet","res_logpt_dijet", n_eta-1,eta_bins);
    TH1D* Residual_const_DiJet = new TH1D("res_const_dijet","res_const_dijet", n_eta-1,eta_bins);
    TH1D* ptave_const_DiJet = new TH1D("ptave_const_dijet","ptave_const_dijet", n_eta-1,eta_bins);
    TH1D* ptave_logpt_DiJet = new TH1D("ptave_logpt_dijet","ptave_logpt_dijet", n_eta-1,eta_bins);
    // for (int j=n_eta-1; j>0; --j){
    //   if(syst==5){
    //     ptave_logpt_DiJet->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()));
    //     ptave_logpt_DiJet->SetBinError(j+1,Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(ptave_data[j]->GetMean()),2)
    // 				       +2*Vcov[2][j]*TMath::Log(ptave_data[j]->GetMean()));                                                                                                                                    
    //   }                                                                                                                                          
    //   ptave_const_DiJet->SetBinContent(j+1,f2[j]->GetParameter(0));                                                                              
    //   ptave_const_DiJet->SetBinError(j+1,f2[j]->GetParError(0));                                                                                 
      
      //   ptave_logpt_MPF->SetBinError(j+1,f2[j]->GetParError(0));                                                                               
    // }                 
    ofstream output, output_loglin, uncerts, uncerts_loglin;
    output.open(path+"output/Spring16_25ns_pT_FLAT_L2Residual_"+txttag+"_"+jettag+".txt");
    output_loglin.open(path+"output/Spring16_25ns_pT_LOGLIN_L2Residual_"+txttag+"_"+jettag+".txt");
    uncerts.open(path+"output/Spring16_25ns_pT_FLAT_L2Residual_"+txttag+"_"+jettag+".txt.STAT");
    uncerts_loglin.open(path+"output/Spring16_25ns_pT_LOGLIN_L2Residual_"+txttag+"_"+jettag+".txt.STAT");
    output  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(208))/208))) Correction L2Relative}" << endl;
    output_loglin  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(208))/208))) Correction L2Relative}" << endl;
    uncerts << "{ 1 JetEta 1 JetPt [0] kFSR_err f0_err}" << endl;
    uncerts_loglin << "{ 1 JetEta 1 JetPt [0] kFSR_err cov(0,0) cov(1,1) cov(0,1) }" << endl;
    // uncerts << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}" << endl;
    // uncerts_loglin << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}" << endl;

    for (int j=n_eta-1; j>0; --j){
      output << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   55   " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "   " << 1/flat_norm << " " << hist_kfsr_fit_dijet->GetBinContent(j) << "   " << f2[j-1]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;
      output_loglin << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   55   " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "   " << 1/loglin_norm << " " << hist_kfsr_fit_dijet->GetBinContent(j) << "   " << f1[j-1]->GetParameter(0) << " " << f1[j-1]->GetParameter(1) <<  "   1 0.0000 0.0" << endl;
      // uncerts << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "       " << sqrt( pow(f2[j-1]->GetParameter(0)*hist_kfsr_fit_dijet->GetBinError(j),2)+pow( hist_kfsr_fit_dijet->GetBinContent(j)*f2[j-1]->GetParError(0),2)  ) / flat_norm  << endl;
      // uncerts_loglin << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "       " << sqrt(  pow(f1[j-1]->GetParameter(0)+f1[j-1]->GetParameter(1)*TMath::Log(ptave_data[j-1]->GetMean()),2)*pow(hist_kfsr_fit_dijet->GetBinError(j),2)  + pow(hist_kfsr_fit_dijet->GetBinContent(j),2)*(Vcov[0][j-1]+Vcov[1][j-1]*pow(TMath::Log(ptave_data[j-1]->GetMean()),2)+2*Vcov[2][j-1]*TMath::Log(ptave_data[j-1]->GetMean()))) / loglin_norm<< endl;
 
      uncerts << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    "
	      << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "  " << hist_kfsr_fit_dijet->GetBinError(j-1) / flat_norm
	      <<" "<< f2[j-1]->GetParError(0) << endl; 
      uncerts_loglin << fixed << std::setprecision(6) << "-" << eta_range[j] << " -"
		     << eta_range[j-1] << "  3    55    " << ptave_data[j-1]->FindLastBinAbove(0.)*10
		     << " " << hist_kfsr_fit_dijet->GetBinError(j-1) / loglin_norm<< " " <<Vcov[0][j-1] <<" "<<Vcov[1][j-1]<<" "<<Vcov[2][j-1]<<endl;   
   }


    for (int j=0; j<n_eta-1; j++){
      
      output << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   55   " << ptave_data[j]->FindLastBinAbove(0.)*10 << "   " << 1/flat_norm << " " << hist_kfsr_fit_dijet->GetBinContent(j+1) << "   " << f2[j]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;
      output_loglin << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   55   " << ptave_data[j]->FindLastBinAbove(0.)*10 << "   " << 1/loglin_norm << " " << hist_kfsr_fit_dijet->GetBinContent(j+1) << "   " << f1[j]->GetParameter(0) << " " << f1[j]->GetParameter(1) <<  "   1 0.0000 0.0" << endl;
      // uncerts << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " << ptave_data[j]->FindLastBinAbove(0.)*10 << "       " << sqrt( pow(f2[j]->GetParameter(0)*hist_kfsr_fit_dijet->GetBinError(j+1),2)+pow( hist_kfsr_fit_dijet->GetBinContent(j+1)*f2[j]->GetParError(0),2)  ) / flat_norm  << endl;
      // uncerts_loglin << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " << ptave_data[j]->FindLastBinAbove(0.)*10 << "       " << sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()),2)*pow(hist_kfsr_fit_dijet->GetBinError(j+1),2)   + pow(hist_kfsr_fit_dijet->GetBinContent(j+1),2)*(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(ptave_data[j]->GetMean()),2)+2*Vcov[2][j]*TMath::Log(ptave_data[j]->GetMean()))) / loglin_norm<< endl;

      uncerts << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    "
	      << ptave_data[j]->FindLastBinAbove(0.)*10<<" "<< hist_kfsr_fit_dijet->GetBinError(j)/flat_norm<< " " 
	      << f2[j]->GetParameter(0) <<endl;                                         
      uncerts_loglin << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] 
		     << "  3    55    "<< ptave_data[j]->FindLastBinAbove(0.)*10 << " "
		     << hist_kfsr_fit_dijet->GetBinError(j)/loglin_norm<< " "<< Vcov[0][j]<<" "<< Vcov[1][j]<<" "<< Vcov[2][j] << endl;     

      double pave_corr;
      if(syst==5){
	pave_corr = ptave_data[j]->GetMean();
	/*
    	Residual_logpt_DiJet->SetBinContent(j+1,hist_kfsr_fit_dijet->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()) ) );
    	Residual_logpt_DiJet->SetBinError(j+1, sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()),2)*pow(hist_kfsr_fit_dijet->GetBinError(j+1),2)   + pow(hist_kfsr_fit_dijet->GetBinContent(j+1),2)*(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(ptave_data[j]->GetMean()),2)+2*Vcov[2][j]*TMath::Log(ptave_data[j]->GetMean()))));

    	ptave_logpt_DiJet->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()));
    	ptave_logpt_DiJet->SetBinError(j+1,sqrt((Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(ptave_data[j]->GetMean()),2)+2*Vcov[2][j]*TMath::Log(ptave_data[j]->GetMean()))));
	*/
      }
      if(syst==0){
	pave_corr = 120.;
	/* 
   	Residual_logpt_DiJet->SetBinContent(j+1,hist_kfsr_fit_dijet->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(120) ) );
    	Residual_logpt_DiJet->SetBinError(j+1, sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(120),2)*pow(hist_kfsr_fit_dijet->GetBinError(j+1),2)   + pow(hist_kfsr_fit_dijet->GetBinContent(j+1),2)*(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(120),2)+2*Vcov[2][j]*TMath::Log(120))));
    	ptave_logpt_DiJet->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(120));
    	ptave_logpt_DiJet->SetBinError(j+1,sqrt((Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(120),2)+2*Vcov[2][j]*TMath::Log(120))));
	*/
      }
      if(syst==1){
	pave_corr = 60.;
	/*
    	Residual_logpt_DiJet->SetBinContent(j+1,hist_kfsr_fit_dijet->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(60) ) );
    	Residual_logpt_DiJet->SetBinError(j+1, sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(60),2)*pow(hist_kfsr_fit_dijet->GetBinError(j+1),2)   + pow(hist_kfsr_fit_dijet->GetBinContent(j+1),2)*(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(60),2)+2*Vcov[2][j]*TMath::Log(60))));
    	ptave_logpt_DiJet->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(60));
    	ptave_logpt_DiJet->SetBinError(j+1,sqrt((Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(60),2)+2*Vcov[2][j]*TMath::Log(60))));
	*/      
}
      if(syst==2){
	pave_corr = 240.;
	/*
    	Residual_logpt_DiJet->SetBinContent(j+1,hist_kfsr_fit_dijet->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(240) ) );
    	Residual_logpt_DiJet->SetBinError(j+1, sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(240),2)*pow(hist_kfsr_fit_dijet->GetBinError(j+1),2)   + pow(hist_kfsr_fit_dijet->GetBinContent(j+1),2)*(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(240),2)+2*Vcov[2][j]*TMath::Log(240))));
    	ptave_logpt_DiJet->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(240));
    	ptave_logpt_DiJet->SetBinError(j+1,sqrt((Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(240),2)+2*Vcov[2][j]*TMath::Log(240))));
	*/
      }
      if(syst==3){
	pave_corr = 480.;
	/*
    	Residual_logpt_DiJet->SetBinContent(j+1,hist_kfsr_fit_dijet->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(480) ) );
    	Residual_logpt_DiJet->SetBinError(j+1, sqrt( pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(480),2)*pow(hist_kfsr_fit_dijet->GetBinError(j+1),2)   + pow(hist_kfsr_fit_dijet->GetBinContent(j+1),2)*(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(480),2)+2*Vcov[2][j]*TMath::Log(480))));
    	ptave_logpt_DiJet->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(480));
    	ptave_logpt_DiJet->SetBinError(j+1,sqrt((Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(480),2)+2*Vcov[2][j]*TMath::Log(480))));
	*/      
}


    	Residual_logpt_DiJet->SetBinContent(j+1,hist_kfsr_fit_dijet->GetBinContent(j+1)*(f1[j]->GetParameter(0)
										     +f1[j]->GetParameter(1)*TMath::Log(pave_corr)));
    	Residual_logpt_DiJet->SetBinError(j+1, sqrt(pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(pave_corr),2)*pow(hist_kfsr_fit_dijet->GetBinError(j+1),2)+ pow(hist_kfsr_fit_dijet->GetBinContent(j+1),2)*(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(pave_corr),2)+2*Vcov[2][j]*TMath::Log(pave_corr))));

    	ptave_logpt_DiJet->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(pave_corr));
    	ptave_logpt_DiJet->SetBinError(j+1,sqrt((Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(ptave_data[j]->GetMean()),2)+2*Vcov[2][j]*TMath::Log(pave_corr))));

      Residual_const_DiJet->SetBinContent(j+1,hist_kfsr_fit_dijet->GetBinContent(j+1)*f2[j]->GetParameter(0));
      Residual_const_DiJet->SetBinError(j+1, sqrt( pow(f2[j]->GetParameter(0)*hist_kfsr_fit_dijet->GetBinError(j+1),2)+pow(hist_kfsr_fit_dijet->GetBinContent(j+1)*f2[j]->GetParError(0),2)  )  );
      ptave_const_DiJet->SetBinContent(j+1,f2[j]->GetParameter(0));
      ptave_const_DiJet->SetBinError(j+1,f2[j]->GetParError(0));
    }





    TFile* outputfile;
    if(syst==5){
      outputfile = new TFile(path+"Histo_Res_DiJet_L1_"+txttag+"_"+jettag+tag+".root","RECREATE");
    }
    if(syst<4){
     outputfile = new TFile(path+"Histo_Res_DiJet_L1_"+txttag+"_"+jettag+"_"+variation+".root","RECREATE");
    }
    /*
    if(syst==1){
      outputfile = new TFile(path+"Histo_Res_DiJet_L1_"+txttag+"_"+jettag+"_"+variation+".root","RECREATE");
    }
    if(syst==2){
      outputfile = new TFile(path+"Histo_Res_DiJet_L1_"+variation+".root","RECREATE");
    }
    if(syst==3){
      outputfile = new TFile(path+"Histo_Res_DiJet_L1_"+variation+".root","RECREATE");
    }
    */
    hist_kfsr_dijet->Write();
    hist_kfsr_fit_dijet->Write();
    kfsr_fit_dijet->Write();
    ptave_const_DiJet->Write();
    ptave_logpt_DiJet->Write();
    /*  for (int j=0; j<n_eta-1; j++){
      f1[j]->Write();
      f2[j]->Write();
      } */              
    Residual_logpt_DiJet->Write();
    Residual_const_DiJet->Write();
    outputfile->Write();
    outputfile->Close();
    //    TFile* outputfile2 = new TFile(path+"Histo_ptave_DiJet_L1.root","RECREATE");
    TFile* outputfile2;
    if(syst==5){
      outputfile2 = new TFile(path+"Histo_ptave_DiJet_L1_"+txttag+"_"+jettag+tag+".root","RECREATE");
    }
    if(syst<4){
     outputfile2 = new TFile(path+"Histo_ptave_DiJet_L1_"+txttag+"_"+jettag+"_"+variation+".root","RECREATE");
    }
    /*
    if(syst==1){
      outputfile2 = new TFile(path+"Histo_ptave_DiJet_L1_"+variation+".root","RECREATE");
    }
    if(syst==2){
      outputfile2 = new TFile(path+"Histo_ptave_DiJet_L1_"+variation+".root","RECREATE");
    }
    if(syst==3){
      outputfile2 = new TFile(path+"Histo_ptave_DiJet_L1_"+variation+".root","RECREATE");
    }
    */

    ptave_const_DiJet->Write();
    ptave_logpt_DiJet->Write();
    for (int j=0; j<n_eta-1; j++){
      f1[j]->Write();
      f2[j]->Write();
    }             
    outputfile2->Write();
    outputfile2->Close();
  }






}
