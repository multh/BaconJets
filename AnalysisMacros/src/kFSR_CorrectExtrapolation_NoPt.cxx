#include "../include/CorrectionObject.h"
#include "../include/parameters.h"
#include "../include/useful_functions.h"

#include <TStyle.h>
#include <TH1.h>
#include <TH1D.h>
#include <TString.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TFile.h>
#include <TLine.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TGraphErrors.h>
#include <TMatrixD.h>
#include <assert.h> 
#include <TMinuit.h>
#include <TVirtualFitter.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>
#include <TMath.h>
#include <TPaveStats.h>


using namespace std;

struct fit_data{
  // x, y values:
  std::vector<double> x_val, y_val, y_val_err;
  // variance-covariance matrix for y values:
  TMatrixD y_cov;
  // inverted cov matrix; calculated by chi2_linear "on demand".
  TMatrixD y_cov_inv;
  
  void reset(){
    x_val.clear();
    y_val.clear();
    y_val_err.clear();
    y_cov.ResizeTo(0,0);
    y_cov_inv.ResizeTo(0,0);
  }
  
  void CheckPoints(){
    std::vector<int> RemovedPoints;
    TMatrixD y_cov_new;
    int j = 0;

    for(int i = 0; i < y_val.size(); i++) {
      if( y_val.at(i) == 0 || y_val_err.at(i) == 0) {
	x_val.erase(x_val.begin()+i);
	y_val.erase(y_val.begin()+i);
	y_val_err.erase(y_val_err.begin()+i);
	RemovedPoints.push_back(j);
	i = i-1;
      }
      j++;
    }    
 
    y_cov_new.ResizeTo(x_val.size(),x_val.size());
    for( int i=0; i < x_val.size(); i++) {
      for(int k= 0; k < x_val.size(); k++) {
	y_cov_new(i,k) = y_cov(i+RemovedPoints.size(),k+RemovedPoints.size());
      }
    }
    y_cov.ResizeTo(0,0);
    y_cov.ResizeTo(x_val.size(),x_val.size());
    y_cov = y_cov_new;  
  }
};

fit_data data;


void chi2_linear(Int_t& npar, Double_t* grad, Double_t& fval, Double_t* p, Int_t status){
    if(data.y_cov_inv.GetNcols()==0){
        double dummy;
        int ncols = data.y_cov.GetNcols();
        data.y_cov_inv.ResizeTo(ncols, ncols);
	
        data.y_cov_inv = data.y_cov.Invert(&dummy);
	
//	data.y_cov_inv.Print();
	
    }
    const size_t ndata = data.x_val.size(); // number of data points in x,y graph to fit to

    std::vector<double> delta_y(ndata);

    for(size_t i=0; i<ndata; ++i){

        delta_y[i] = data.x_val[i]*p[0] + p[1] - data.y_val[i];
    }
    // now calculate the chi2, i.e.
    //  dy^T * C^{-1} * dy
    // where C is the variance--covariance matrix and dy = (y_data - y_pred)
    // This could probably be implemented in ROOT, but it's so simple, we just do it here:
    fval = 0.0;
    for(size_t i=0; i<ndata; ++i){
        for(size_t j=0; j<ndata; ++j){
	  fval += delta_y[i] * delta_y[j] * data.y_cov_inv(i,j);
        }
    }
}



void make_lin_fit(double & slope, double & d_slope, double & offset, double & d_offset){
  TMinuit *min = new TMinuit();
  min->SetPrintLevel(-1);
  //  min->SetPrintLevel(0);
  int err = min->DefineParameter(0, "slope", slope, d_slope, -2*fabs(slope),2*fabs(slope));
  assert(err==0);
  err = min->DefineParameter(1, "offset", offset, d_offset, -2*fabs(offset),2*fabs(offset));
  assert(err==0);
  
  min->SetFCN(chi2_linear);
  
  double arglist[10];
  int ierrflag =0;
  
  arglist[0] =1;
  min->mnexcm("SET ERR", arglist ,1,ierrflag);
  
  arglist[0]=100;
  arglist[1]=1.;
  
  min->mnexcm("MIGRAD", arglist ,2,ierrflag);
  
  min->mnmigr();
  
  TMatrixD matrix0(2,2);
  min->mnemat(matrix0.GetMatrixArray(),2);
  //matrix0.Print();
  
  min->GetParameter(0, slope, d_slope);
  min->GetParameter(1, offset, d_offset);
}


TGraphErrors* Fit_User(TH1D** hist){
  
  vector<double> alpha;
  //  alpha.push_back(0.1); 
  alpha.push_back(0.15); 
  alpha.push_back(0.2); 
  alpha.push_back(0.25); 
  alpha.push_back(0.3);
  alpha.push_back(0.35); 
  alpha.push_back(0.4);
  //alpha.push_back(0.45);
  
  
  vector<double> mean_value;
  vector<double> mean_error;
  vector<double> std_dev;

  vector<double> mean_graph;
  vector<double> error_graph;
  
  for(int i=0; i<alpha.size();i++){
    if(hist[i+2]->GetEntries()>100){  
      mean_value.push_back(hist[i+2]->GetMean());
      mean_error.push_back((hist[i+2]->GetMeanError()));
      std_dev.push_back((hist[i+2]->GetStdDev()));
    }
    else{
      mean_value.push_back(0);
      mean_error.push_back(0);
      std_dev.push_back(0);
    }
  }
  
  for(int i=0; i<n_alpha;i++){
    if(hist[i]->GetEntries()>100){  
      mean_graph.push_back(hist[i]->GetMean());
      error_graph.push_back((hist[i]->GetMeanError()));
    }
    else{
      mean_graph.push_back(0);
      error_graph.push_back(0);
    }
  }


  TMatrixD Cov_alpha;
  Cov_alpha.ResizeTo(alpha.size(), alpha.size());
  
  // fill covariance matrix for data and mc
  for(int ialpha=0; ialpha < alpha.size(); ++ialpha){
    for (Int_t jalpha =0; jalpha < alpha.size(); jalpha++){
      if( ialpha <= jalpha ) {
	Cov_alpha(ialpha, jalpha) = pow(std_dev.at(ialpha),2)*(pow(mean_error.at(jalpha),2)/pow(std_dev.at(jalpha),2)); 
      }
      else {
	Cov_alpha(ialpha, jalpha) = pow(std_dev.at(jalpha),2)*(pow(mean_error.at(ialpha),2)/pow(std_dev.at(ialpha),2));
      }
    }
  }
  
  
  // fit linear extrapolation function
  TF1 *lin_extrapol = new TF1("lin_extrapol","[0]+[1]*x",0,alpha.back()+0.05); 
  
  
  data.reset();
  data.x_val = alpha;
  data.y_val = mean_value;
  data.y_val_err = mean_error;
  data.y_cov.ResizeTo(alpha.size(),alpha.size());
  data.y_cov = Cov_alpha;
  
  data.CheckPoints();

  //    data.y_cov.Print();
  
  double slope = (mean_value[alpha.size()-1] - mean_value[alpha.size()-3])/(alpha.at(alpha.size()-1) - alpha.at(alpha.size()-3));
  double d_slope = abs(slope*0.05);
  double offset = mean_value[alpha.size()-1] - (slope*alpha.at(alpha.size()-3));
  double d_offset = abs(offset*0.05);
  
  make_lin_fit(slope, d_slope, offset, d_offset);
  
  //Set Parameters
  if(data.y_cov.GetNcols()>3){
    lin_extrapol->SetParameter(0, offset);
    lin_extrapol->SetParError(0, d_offset);
    lin_extrapol->SetParameter(1, slope);
    lin_extrapol->SetParError(1, d_slope);
  }
  else{
    lin_extrapol->SetParameter(0, 0);
    lin_extrapol->SetParError(0, 0);
    lin_extrapol->SetParameter(1, 0);
    lin_extrapol->SetParError(1, 0);
  }  
  
  //Setup Graph
  double xbin_tgraph[n_alpha],zero[n_alpha], Mean[n_alpha], Error[n_alpha];
  for(int i=0;i<n_alpha;i++){
    xbin_tgraph[i] = alpha_bins[i];
    zero[i] = 0;
    Mean[i] = mean_graph[i];
    Error[i]= error_graph[i];
  }
  
  TGraphErrors* graph = new TGraphErrors(n_alpha,xbin_tgraph,Mean,zero,Error);
  graph = (TGraphErrors*)CleanEmptyPoints(graph);
  graph->GetListOfFunctions()->Add(lin_extrapol);
  
  return graph;
}


  

void CorrectionObject::kFSR_CorrectExtrapolation(){
  cout << "--------------- Starting kFSR() ---------------" << endl << endl;
  //  TStyle* m_gStyle = new TStyle();
  // m_gStyle->SetOptFit(0);
  TH1::SetDefaultSumw2(kTRUE);

  int countPt = 0;
  TH1D* ptave_data[n_eta-1];
  TString namehist = "ptave_";
  for(int i=0; i<n_eta-1; i++){
    TString selection = "alpha<";
    selection+=alpha_cut;
    selection+=" && fabs(probejet_eta)<";
    selection+=eta_range[i+1];
    selection+=" && fabs(probejet_eta)>=";
    selection+=eta_range[i];
    TString var1 = "pt_ave";
    ptave_data[i] = (TH1D*)GetHist(CorrectionObject::_DATAFile, selection, var1, 300,0,3000)->Clone();
    TString namecur = namehist;
    namecur += countPt;
    ptave_data[i]->SetName(namecur);
    countPt++;
  }


  TH1D *hdata_asymmetry[n_eta-1][n_alpha];
  TH1D *hdata_B[n_eta-1][n_alpha];
  TH1D *hmc_asymmetry[n_eta-1][n_alpha];
  TH1D *hmc_B[n_eta-1][n_alpha];

  int count = 0;

  TString name1 = "hist_data_asymmetry_";
  TString name2 = "hist_data_B_";
  TString name3 = "hist_mc_asymmetry_";
  TString name4 = "hist_mc_B_";

  for(int i=0; i<n_alpha; i++){
    for(int j=0; j<n_eta-1; j++){
      TString name = name1; name+=count;
      hdata_asymmetry[j][i] = new TH1D(name,"A in DATA; A ; Events ",nResponseBins, -1.2, 1.2);
      name = name2;name+=count;
      hdata_B[j][i]         = new TH1D(name,"B in DATA; B ; Events",nResponseBins, -1.2, 1.2);
      name = name3; name+=count;
      hmc_asymmetry[j][i]   = new TH1D(name,"A in MC; A ; Events",nResponseBins, -1.2, 1.2);
      name = name4; name+=count;
      hmc_B[j][i]           = new TH1D(name,"B in MC;B;Events", nResponseBins, -1.2, 1.2);
      
      count++;
    }
  }
  
  TH1D *hdata_ptave[n_eta-1][n_pt-1];//pt-ave in each bin of pT_ave in bins of eta

  for(int k=0; k<n_eta-1; k++){
    for(int j=0; j<n_pt-1; j++){
      TString name = "hist_data_pt_ave_eta_"+eta_range2[k]+"_"+eta_range2[k+1]+"pt_"+pt_range[j]+"_"+pt_range[j+1];
      hdata_ptave[k][j] = new TH1D(name,"",3000,0,3000); //only used for GetMean and GetStdDev in the pT extrapolations
    }
  }


  cout << "Set up a total of " << count << " histograms." << endl;
    
  // Get relevant quantities from DATA, loop over events
  TTreeReader myReader_DATA("AnalysisTree", CorrectionObject::_DATAFile);
  TTreeReaderValue<Float_t> pt_ave_data(myReader_DATA, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_data(myReader_DATA, "probejet_eta");
  TTreeReaderValue<Float_t> alpha_data(myReader_DATA, "alpha");
  TTreeReaderValue<Float_t> asymmetry_data(myReader_DATA, "asymmetry");
  TTreeReaderValue<Float_t> B_data(myReader_DATA, "B");   
  TTreeReaderValue<Float_t> weight_data(myReader_DATA, "weight");
  int idx = 0;
  
  cout << "starting to loop over DATA events." << endl;
  
  while (myReader_DATA.Next()) {
    for(int j=0; j<n_eta-1; j++){
      if(fabs(*probejet_eta_data)>eta_bins[j+1] || fabs(*probejet_eta_data)<eta_bins[j]) continue;
      for(int i=0; i<n_alpha; i++){
	if(*alpha_data>alpha_bins[i]) continue;
	else{
	  hdata_asymmetry[j][i]->Fill(*asymmetry_data,*weight_data);
	  hdata_B[j][i]->Fill(*B_data,*weight_data);
	  idx++;
	}
      }
      for(int k=0; k<n_pt-1; k++){
	if(*pt_ave_data>pt_bins[k+1] || *pt_ave_data<pt_bins[k]) continue;
	hdata_ptave[j][k]->Fill(*pt_ave_data,*weight_data);
      }
    }
  }
  cout << "Finished running over DATA events. Read in a total of " << idx << " events." << endl;
  cout << "starting to loop over MC events." << endl;

  // Get relevant quantities from MC, loop over events
  TTreeReader myReader_MC("AnalysisTree", CorrectionObject::_MCFile);
  TTreeReaderValue<Float_t> pt_ave_mc(myReader_MC, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_mc(myReader_MC, "probejet_eta");
  TTreeReaderValue<Float_t> alpha_mc(myReader_MC, "alpha");
  TTreeReaderValue<Float_t> asymmetry_mc(myReader_MC, "asymmetry");
  TTreeReaderValue<Float_t> B_mc(myReader_MC, "B");
  TTreeReaderValue<Float_t> weight_mc(myReader_MC, "weight");
  idx = 0;
  
  while (myReader_MC.Next()) {
    for(int j=0; j<n_eta-1; j++){
      if(fabs(*probejet_eta_mc)>eta_bins[j+1] || fabs(*probejet_eta_mc)<eta_bins[j]) continue;
      for(int i=0; i<n_alpha; i++){
	if(*alpha_mc>alpha_bins[i]) continue;
	else{
	  hmc_asymmetry[j][i]->Fill(*asymmetry_mc,*weight_mc);
	  hmc_B[j][i]->Fill(*B_mc,*weight_mc);
	  idx++;
	}
      }
    }
  }
  
  cout << "Finished running over MC events. Read in a total of " << idx << " events." << endl;

 /*
 
  for(int j=0; j<n_eta-1; j++){
    for(int k=0; k<2; k++){
      for(int i=0; i<n_alpha; i++){
	//print TH2D histos
	TCanvas* c_dummy1 = new TCanvas();
	hdata_asymmetry[j][k][i]->GetYaxis()->SetTitleSize(0.048);
	hdata_asymmetry[j][k][i]->GetYaxis()->SetTitleOffset(0.6);
	hdata_asymmetry[j][k][i]->GetXaxis()->SetTitleSize(0.048);
	hdata_asymmetry[j][k][i]->Draw();
	c_dummy1->SaveAs(CorrectionObject::_outpath+"plots/control/TH1_A_DATA_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_pT_"+pt_range[k]+"_"+pt_range[k+1]+"_"+alpha_range[i]+".pdf");
	delete c_dummy1;
	TCanvas* c_dummy2 = new TCanvas();
	hmc_asymmetry[j][k][i]->GetYaxis()->SetTitleSize(0.048);
	hmc_asymmetry[j][k][i]->GetYaxis()->SetTitleOffset(0.6);
	hmc_asymmetry[j][k][i]->GetXaxis()->SetTitleSize(0.05);
	hmc_asymmetry[j][k][i]->Draw();
	c_dummy2->SaveAs(CorrectionObject::_outpath+"plots/control/TH1_A_MC_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_pT_"+pt_range[k]+"_"+pt_range[k+1]+"_"+alpha_range[i]+".pdf");
	delete c_dummy2;
	TCanvas* c_dummy3 = new TCanvas();
	hdata_B[j][k][i]->GetYaxis()->SetTitleSize(0.048);
	hdata_B[j][k][i]->GetYaxis()->SetTitleOffset(0.6);
	hdata_B[j][k][i]->GetXaxis()->SetTitleSize(0.048);
	hdata_B[j][k][i]->Draw();
	c_dummy3->SaveAs(CorrectionObject::_outpath+"plots/control/TH1_B_DATA_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_pT_"+pt_range[k]+"_"+pt_range[k+1]+"_"+alpha_range[i]+".pdf");
	delete c_dummy3;
	TCanvas* c_dummy4 = new TCanvas();
	hmc_B[j][k][i]->GetYaxis()->SetTitleSize(0.048);
	hmc_B[j][k][i]->GetYaxis()->SetTitleOffset(0.6);
	hmc_B[j][k][i]->GetXaxis()->SetTitleSize(0.048);
	hmc_B[j][k][i]->Draw();
	c_dummy4->SaveAs(CorrectionObject::_outpath+"plots/control/TH1_B_MC_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_pT_"+pt_range[k]+"_"+pt_range[k+1]+"_"+alpha_range[i]+".pdf");
	delete c_dummy4;
      }
    }
  }
 */
  
  
  TGraphErrors *graph_rel_mc[n_eta-1];
  TGraphErrors *graph_rel_data[n_eta-1];
  
  TGraphErrors *graph_mpf_mc[n_eta-1];
  TGraphErrors *graph_mpf_data[n_eta-1];
  
  TF1 *Rel_data_Extrapol[n_eta-1];
  TF1 *Rel_mc_Extrapol[n_eta-1];
  TF1 *MPF_data_Extrapol[n_eta-1];
  TF1 *MPF_mc_Extrapol[n_eta-1];


  
  double xbin_tgraph[n_eta-1],zero[n_eta-1];
  double MPF_data_mean_ex[n_eta-1], MPF_data_error_ex[n_eta-1];
  double MPF_mc_mean_ex[n_eta-1], MPF_mc_error_ex[n_eta-1];
  double Rel_data_mean_ex[n_eta-1], Rel_data_error_ex[n_eta-1];
  double Rel_mc_mean_ex[n_eta-1], Rel_mc_error_ex[n_eta-1];

  /*
  TGraphErrors *graph_mpf_mc_ex[n_eta-1];
  TGraphErrors *graph_mpf_data_ex[n_eta-1];
  TGraphErrors *graph_rel_mc_ex[n_eta-1];
  TGraphErrors *graph_rel_data_ex[n_eta-1];
  */
  double Ratio_Res_Rel[n_eta-1], Error_Res_Rel[n_eta-1];
  double Ratio_Res_MPF[n_eta-1], Error_Res_MPF[n_eta-1];

  /*
  TGraphErrors *graph_ratio_res_rel[n_eta-1];
  TGraphErrors *graph_ratio_res_mpf[n_eta-1];
  
  TString plotname_rel[n_eta-1];
  double Vcov_rel[3][n_eta-1];

  TString plotname_mpf[n_eta-1];
  double Vcov_mpf[3][n_eta-1];

  // create a function for the loglinear fit
  TF1 * f1_rel[n_eta-1];
  // create a function for the constant fit
  TF1 * f2_rel[n_eta-1];

  // create a function for the loglinear fit
  TF1 * f1_mpf[n_eta-1];
  // create a function for the constant fit
  TF1 * f2_mpf[n_eta-1];
  
  TH1D* h_chi2_loglin_mpf = new TH1D("h_chi2_loglin_mpf", "Chi2/ndf for each eta bin;|#eta|;#chi^{2}/ndf", n_eta-1, eta_bins);
  TH1D* h_chi2_const_mpf =  new TH1D("h_chi2_const_mpf", "Chi2/ndf for each eta bin;|#eta|;#chi^{2}/ndf", n_eta-1, eta_bins);

  TH1D* h_chi2_loglin_rel = new TH1D("h_chi2_loglin_rel", "Chi2/ndf for each eta bin;|#eta|;#chi^{2}/ndf", n_eta-1, eta_bins);
  TH1D* h_chi2_const_rel =  new TH1D("h_chi2_const_rel", "Chi2/ndf for each eta bin;|#eta|;#chi^{2}/ndf", n_eta-1, eta_bins);
  */
  TH1D* h_kFSR_mpf = new TH1D("h_kFSR_mpf", "kFSR MPF;|#eta|; kFSR_{#alpha->0}", n_eta-1, eta_bins);
  TH1D* h_kFSR_rel =  new TH1D("h_kFSR_rel", "kFSR DiJet;|#eta|; kFSR_{#alpha->0}", n_eta-1, eta_bins);

  for(int j=0; j<n_eta-1; j++){
      
      graph_rel_mc[j]=(TGraphErrors*)Fit_User(hmc_asymmetry[j]);
      graph_mpf_mc[j]=(TGraphErrors*)Fit_User(hmc_B[j]);
      graph_mpf_data[j]=(TGraphErrors*)Fit_User(hdata_B[j]);
      graph_rel_data[j]=(TGraphErrors*)Fit_User(hdata_asymmetry[j]);
      
      Rel_data_Extrapol[j] = graph_rel_data[j]->GetFunction("lin_extrapol");
      Rel_mc_Extrapol[j]   = graph_rel_mc[j]->GetFunction("lin_extrapol");

      MPF_data_Extrapol[j] = graph_mpf_data[j]->GetFunction("lin_extrapol");
      MPF_mc_Extrapol[j]   = graph_mpf_mc[j]->GetFunction("lin_extrapol");


      TLegend *leg1;
      leg1 = new TLegend(0.15,0.68,0.35,0.88,"","brNDC");//x+0.1
      leg1->SetBorderSize(0);
      leg1->SetTextSize(0.03);
      leg1->SetFillColor(10);
      leg1->SetLineColor(1);
      leg1->SetTextFont(42);


      // draw extrapolations  mc
      TCanvas *c1 = new TCanvas("c1","",600,600);
      c1->DrawFrame(0,-0.29,0.5,0.19,(";Threshold #alpha_{max};#mu_{A}"));
      graph_rel_mc[j]->GetYaxis()->SetRangeUser(0.8, 1.2);
      graph_rel_mc[j]->SetMarkerStyle(20);
      graph_rel_mc[j]->SetMarkerColor(kRed+1);
      graph_rel_mc[j]->SetLineColor(kRed+1);
      
      graph_rel_mc[j]->Draw("P");
      
      TF1* Temp = new TF1();
      graph_rel_mc[j]->GetFunction("lin_extrapol")->SetLineColor(kRed+1);
      graph_rel_mc[j]->GetFunction("lin_extrapol")->SetLineStyle(2);
      
      Temp=(TF1*) graph_rel_mc[j]->GetFunction("lin_extrapol")->Clone();
      Temp->SetRange(0.1,1);
      Temp->SetLineStyle(1);
      Temp->Draw("same");

      leg1->SetHeader("p_{T}-balance, #alpha extrapolation, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
      leg1->AddEntry(graph_rel_mc[j], "#mu_{A, MC}","P");
      leg1->AddEntry(Temp, "linear fit","L");
      leg1->Draw();

      TLatex *tex1 = new TLatex();
      tex1->SetNDC();
      tex1->SetTextSize(0.045); 
      tex1->DrawLatex(0.5,0.91,CorrectionObject::_lumitag+"(13TeV)"); 

      c1->Print(CorrectionObject::_outpath+"plots/kFSR_Rel_MC_Extrapol_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
      delete c1;
      delete tex1;
      leg1->Clear();


      // draw extrapolations mc
      TCanvas *c2 = new TCanvas("c2","",600,600);
      c2->DrawFrame(0,-0.29,0.5,0.19,(";Threshold #alpha_{max};#mu_{B}"));
      graph_mpf_mc[j]->GetYaxis()->SetRangeUser(0.8, 1.2);
      graph_mpf_mc[j]->SetMarkerStyle(20);
      graph_mpf_mc[j]->SetMarkerColor(kRed+1);
      graph_mpf_mc[j]->SetLineColor(kRed+1);
      
      graph_mpf_mc[j]->Draw("P");

      graph_mpf_mc[j]->GetFunction("lin_extrapol")->SetLineColor(kRed+1);
      graph_mpf_mc[j]->GetFunction("lin_extrapol")->SetLineStyle(2);
      
      Temp=(TF1*) graph_mpf_mc[j]->GetFunction("lin_extrapol")->Clone();
      Temp->SetRange(0.15,1);
      Temp->SetLineStyle(1);
      Temp->Draw("same");

      leg1->SetHeader("MPF, #alpha extrapolation, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
      leg1->AddEntry(graph_mpf_mc[j], "#mu_{B, MC}","P");
      leg1->AddEntry(Temp, "linear fit","L");
      leg1->Draw();

      TLatex *tex2 = new TLatex();
      tex2->SetNDC();
      tex2->SetTextSize(0.045); 
      tex2->DrawLatex(0.5,0.91,CorrectionObject::_lumitag+"(13TeV)"); 
      
      c2->Print(CorrectionObject::_outpath+"plots/kFSR_MPF_MC_Extrapol_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
      delete c2;
      delete tex2;
      delete Temp;
      leg1->Clear();
      


      // draw extrapolations data 
      TCanvas *c3 = new TCanvas("c3","",600,600);
      c3->DrawFrame(0,-0.29,0.5,0.19,(";Threshold #alpha_{max};#mu_{A}"));
      graph_rel_data[j]->GetYaxis()->SetRangeUser(0.8, 1.2);
      graph_rel_data[j]->SetMarkerStyle(20);
      graph_rel_data[j]->SetMarkerColor(kRed+1);
      graph_rel_data[j]->SetLineColor(kRed+1);
      
      graph_rel_data[j]->Draw("P");

      graph_rel_data[j]->GetFunction("lin_extrapol")->SetLineColor(kRed+1);
      graph_rel_data[j]->GetFunction("lin_extrapol")->SetLineStyle(2);

      Temp=(TF1*) graph_rel_data[j]->GetFunction("lin_extrapol")->Clone();
      Temp->SetRange(0.15,1);
      Temp->SetLineStyle(1);
      Temp->Draw("same");

      leg1->SetHeader("p_{T} bal, #alpha extrapolation, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
      leg1->AddEntry(graph_rel_data[j], "#mu_{A, Data}","P");
      leg1->AddEntry(Temp, "linear fit","L");
      leg1->Draw();

      TLatex *tex3 = new TLatex();
      tex3->SetNDC();
      tex3->SetTextSize(0.045); 
      tex3->DrawLatex(0.5,0.91,CorrectionObject::_lumitag+"(13TeV)"); 

      c3->Print(CorrectionObject::_outpath+"plots/kFSR_Rel_DATA_Extrapol_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
      delete c3;
      delete tex3;
      delete Temp;
      leg1->Clear();
      
      
      // draw extrapolations data + data
      TCanvas *c4 = new TCanvas("c4","",600,600);
      c4->DrawFrame(0,-0.29,0.5,0.19,(";Threshold #alpha_{max};#mu_{B}"));
      graph_mpf_data[j]->GetYaxis()->SetRangeUser(0.8, 1.2);
      graph_mpf_data[j]->SetMarkerStyle(20);
      graph_mpf_data[j]->SetMarkerColor(kRed+1);
      graph_mpf_data[j]->SetLineColor(kRed+1);
      
      graph_mpf_data[j]->Draw("P");
      
      graph_mpf_data[j]->GetFunction("lin_extrapol")->SetLineColor(kRed+1);
      graph_mpf_data[j]->GetFunction("lin_extrapol")->SetLineStyle(2);
      
      Temp=(TF1*) graph_mpf_data[j]->GetFunction("lin_extrapol")->Clone();
      Temp->SetRange(0.15,1);
      Temp->SetLineStyle(1);
      Temp->Draw("same");

      leg1->SetHeader("MPF, #alpha extrapolation, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
      leg1->AddEntry(graph_mpf_data[j], "#mu_{B, Data}","P");
      leg1->AddEntry(Temp, "linear fit","L");
      leg1->Draw();

      TLatex *tex4 = new TLatex();
      tex4->SetNDC();
      tex4->SetTextSize(0.045); 
      tex4->DrawLatex(0.5,0.91,CorrectionObject::_lumitag+"(13TeV)"); 

      c4->Print(CorrectionObject::_outpath+"plots/kFSR_MPF_DATA_Extrapol_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
      delete c4;
      delete tex4;
      delete Temp;
      leg1->Clear();

      Rel_data_mean_ex[j] = (1+Rel_data_Extrapol[j]->GetParameter(0))/(1-Rel_data_Extrapol[j]->GetParameter(0)); 
      if(Rel_data_Extrapol[j]->GetParameter(0)==0) Rel_data_mean_ex[j] = 0;
      Rel_data_error_ex[j]= 2/(pow((Rel_data_Extrapol[j]->GetParameter(0)-1),2)) * Rel_data_Extrapol[j]->GetParError(0); 

      Rel_mc_mean_ex[j] = (1+Rel_mc_Extrapol[j]->GetParameter(0))/(1-Rel_mc_Extrapol[j]->GetParameter(0)); 
      if(Rel_mc_Extrapol[j]->GetParameter(0)==0) Rel_mc_mean_ex[j] = 0;
      Rel_mc_error_ex[j]= 2/(pow((Rel_mc_Extrapol[j]->GetParameter(0)-1),2)) * Rel_mc_Extrapol[j]->GetParError(0); 

 
      MPF_data_mean_ex[j] = (1+MPF_data_Extrapol[j]->GetParameter(0))/(1-MPF_data_Extrapol[j]->GetParameter(0)); 
      if(MPF_data_Extrapol[j]->GetParameter(0)==0) MPF_data_mean_ex[j] = 0;
      MPF_data_error_ex[j]= 2/(pow((MPF_data_Extrapol[j]->GetParameter(0)-1),2)) * MPF_data_Extrapol[j]->GetParError(0); 

      MPF_mc_mean_ex[j] = (1+MPF_mc_Extrapol[j]->GetParameter(0))/(1-MPF_mc_Extrapol[j]->GetParameter(0)); 
      if(MPF_mc_Extrapol[j]->GetParameter(0)==0) MPF_mc_mean_ex[j] = 0;
      MPF_mc_error_ex[j]= 2/(pow((MPF_mc_Extrapol[j]->GetParameter(0)-1),2)) * MPF_mc_Extrapol[j]->GetParError(0);

      if(Rel_data_mean_ex[j] > 0){
	Ratio_Res_Rel[j] = Rel_mc_mean_ex[j]/Rel_data_mean_ex[j];
	Error_Res_Rel[j] = sqrt(pow(1/Rel_data_mean_ex[j]*Rel_mc_error_ex[j],2)+pow(Rel_mc_mean_ex[j]/(Rel_data_mean_ex[j]*Rel_data_mean_ex[j])*Rel_data_error_ex[j],2));
	}
      else{
	Ratio_Res_Rel[j] = 0;
	Error_Res_Rel[j] = 0;
      }

      if(MPF_data_mean_ex[j] > 0){
	Ratio_Res_MPF[j] = MPF_mc_mean_ex[j]/MPF_data_mean_ex[j];
	Error_Res_MPF[j] = sqrt(pow(MPF_mc_error_ex[j]/MPF_data_mean_ex[j],2)+pow(MPF_mc_mean_ex[j]/(MPF_data_mean_ex[j]*MPF_data_mean_ex[j])*MPF_data_error_ex[j],2));
      }
      else{
	Ratio_Res_MPF[j] = 0;
	Error_Res_MPF[j] = 0;
      }

      double Ratio_Res_Rel_kFSR = Ratio_Res_Rel[j];
      double Ratio_Res_MPF_kFSR = Ratio_Res_MPF[j];
      double Error_Res_Rel_kFSR = Error_Res_Rel[j];
      double Error_Res_MPF_kFSR = Error_Res_MPF[j];

      h_kFSR_rel->SetBinContent(j+1,Ratio_Res_Rel_kFSR); //to make sure to fill the right eta-bin...
      h_kFSR_rel->SetBinError(j+1, Error_Res_Rel_kFSR);
      h_kFSR_mpf->SetBinContent(j+1,Ratio_Res_MPF_kFSR);   //to make sure to fill the right eta-bin...
      h_kFSR_mpf->SetBinError(j+1, Error_Res_MPF_kFSR);


  }


  TCanvas* c_kFSR_mpf = new TCanvas();
  h_kFSR_mpf->GetYaxis()->SetRangeUser(0.85, 1.15);
  h_kFSR_mpf->Draw("P");
  c_kFSR_mpf ->SaveAs(CorrectionObject::_outpath+"plots/kFSR_MPF_Comparison_"+CorrectionObject::_generator_tag+".pdf");
  delete c_kFSR_mpf;
  delete h_kFSR_mpf;


  TCanvas* c_kFSR_rel = new TCanvas();
  h_kFSR_rel->GetYaxis()->SetRangeUser(0.85, 1.15);
  h_kFSR_rel->Draw("P");
  c_kFSR_rel ->SaveAs(CorrectionObject::_outpath+"plots/kFSR_REL_Comparison_"+CorrectionObject::_generator_tag+".pdf");
  delete c_kFSR_rel;
  delete h_kFSR_rel;
	



      /*
    graph_rel_data_ex[j] = new TGraphErrors(n_pt,xbin_tgraph[j],Rel_data_mean_ex[j],zero[j],Rel_data_error_ex[j]);
    graph_rel_data_ex[j] = (TGraphErrors*)CleanEmptyPoints(graph_rel_data_ex[j]);

    graph_rel_mc_ex[j] = new TGraphErrors(n_pt,xbin_tgraph[j],Rel_mc_mean_ex[j],zero[j],Rel_mc_error_ex[j]);
    graph_rel_mc_ex[j] = (TGraphErrors*)CleanEmptyPoints(graph_rel_mc_ex[j]);

    graph_mpf_data_ex[j] = new TGraphErrors(n_pt,xbin_tgraph[j],MPF_data_mean_ex[j],zero[j],MPF_data_error_ex[j]);
    graph_mpf_data_ex[j] = (TGraphErrors*)CleanEmptyPoints(graph_mpf_data_ex[j]);

    graph_mpf_mc_ex[j] = new TGraphErrors(n_pt,xbin_tgraph[j],MPF_mc_mean_ex[j],zero[j],MPF_mc_error_ex[j]);
    graph_mpf_mc_ex[j] = (TGraphErrors*)CleanEmptyPoints(graph_mpf_mc_ex[j]);




    graph_ratio_res_rel[j] = new TGraphErrors(n_pt, xbin_tgraph[j], Ratio_Res_Rel[j],zero[j],Error_Res_Rel[j]);
    graph_ratio_res_rel[j] = (TGraphErrors*)CleanEmptyPoints(graph_ratio_res_rel[j]);

    graph_ratio_res_mpf[j] = new TGraphErrors(n_pt, xbin_tgraph[j], Ratio_Res_MPF[j],zero[j],Error_Res_MPF[j]);
    graph_ratio_res_mpf[j] = (TGraphErrors*)CleanEmptyPoints(graph_ratio_res_mpf[j]);



    plotname_rel[j]="rel_ptextra_"+CorrectionObject::_generator_tag+"_eta_"+eta_range[j]+"_"+eta_range[j+1];
    //Log Linear/Constant fit for pT extrapolation
    f1_rel[j] = new TF1(plotname_rel[j]+"f1","[0]+[1]*TMath::Log(x)", 60 , 300);
    f1_rel[j]->SetParameters(1,0);

    f2_rel[j] = new TF1(plotname_rel[j]+"f2","pol0", 60 , 300);
    f2_rel[j]->SetLineColor(kBlue);
    f2_rel[j]->SetLineStyle(3);
   
    //Do the fit!
    TFitResultPtr fitloglin_rel;
    if(j==17){
      double slope_rel = f1_rel[j-1]->GetParameter(1.);
      f1_rel[j]->FixParameter(1,slope_rel); //fixing slope for last eta bin (Vidyo Meeting 07.07.2017)
    }
    fitloglin_rel =  graph_ratio_res_rel[j]->Fit(plotname_rel[j]+"f1","SMR","",51,1200);
    TMatrixDSym cov_rel = fitloglin_rel->GetCovarianceMatrix();
    Vcov_rel[0][j] = cov_rel(0,0);
    Vcov_rel[1][j] = cov_rel(1,1);     
    Vcov_rel[2][j] = cov_rel(0,1);
    
    graph_ratio_res_rel[j]->Fit(plotname_rel[j]+"f2","SMR + SAME","",51,1200); 
    

    plotname_mpf[j]="mpf_ptextra_"+CorrectionObject::_generator_tag+"_eta_"+eta_range[j]+"_"+eta_range[j+1];
    //Log Linear/Constant fit for pT extrapolation
    f1_mpf[j] = new TF1(plotname_mpf[j]+"f1","[0]+[1]*TMath::Log(x)", 51 , 1200);
    f1_mpf[j]->SetParameters(1,0);

    f2_mpf[j] = new TF1(plotname_mpf[j]+"f2","pol0", 510 , 1200);
    f2_mpf[j]->SetLineColor(kBlue);
    f2_mpf[j]->SetLineStyle(3);
   
    //Do the fit!
    TFitResultPtr fitloglin_mpf;
    if(j==17){
      double slope_mpf = f1_mpf[j-1]->GetParameter(1.);
      f1_mpf[j]->FixParameter(1,slope_mpf); //fixing slope for last eta bin (Vidyo Meeting 07.07.2017)
    }
    fitloglin_mpf =  graph_ratio_res_mpf[j]->Fit(plotname_mpf[j]+"f1","SMR","", 51, 1200);
    TMatrixDSym cov_mpf = fitloglin_mpf->GetCovarianceMatrix();
    Vcov_mpf[0][j] = cov_mpf(0,0);
    Vcov_mpf[1][j] = cov_mpf(1,1);     
    Vcov_mpf[2][j] = cov_mpf(0,1);
    
    graph_ratio_res_mpf[j]->Fit(plotname_mpf[j]+"f2","SMR + SAME","", 51, 1200); 


    double chi2ndf_loglin_mpf = f1_mpf[j]->GetChisquare() / f1_mpf[j]->GetNDF();
    double chi2ndf_const_mpf  = f2_mpf[j]->GetChisquare() / f2_mpf[j]->GetNDF();
    h_chi2_loglin_mpf->SetBinContent(j+1,chi2ndf_loglin_mpf); //to make sure to fill the right eta-bin...
    h_chi2_const_mpf ->SetBinContent(j+1,chi2ndf_const_mpf);   //to make sure to fill the right eta-bin...

    double chi2ndf_loglin_rel = f1_rel[j]->GetChisquare() / f1_rel[j]->GetNDF();
    double chi2ndf_const_rel  = f2_rel[j]->GetChisquare() / f2_rel[j]->GetNDF();
    h_chi2_loglin_rel->SetBinContent(j+1,chi2ndf_loglin_rel); //to make sure to fill the right eta-bin...
    h_chi2_const_rel ->SetBinContent(j+1,chi2ndf_const_rel);   //to make sure to fill the right eta-bin...




    TLegend *leg2;
    leg2 = new TLegend(0.15,0.68,0.35,0.88,"","brNDC");//x+0.1
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.038);
    leg2->SetFillColor(10);
    leg2->SetLineColor(1);
    leg2->SetTextFont(42);

    // draw extrapolations  data
    TCanvas *c5 = new TCanvas("c5","",850,700);
    c5->DrawFrame(30,0.5,2100,1.5,(";#bar{p}_{T};Response"));
    gPad->SetLogx();
    graph_rel_data_ex[j]->SetMarkerStyle(20);
    graph_rel_data_ex[j]->SetMarkerColor(kBlack);
    graph_rel_data_ex[j]->SetLineColor(kBlack);
    
    graph_rel_data_ex[j]->Draw("P");
    
    graph_rel_mc_ex[j]->SetMarkerStyle(20);
    graph_rel_mc_ex[j]->SetMarkerColor(kRed);
    graph_rel_mc_ex[j]->SetLineColor(kRed);
    
    graph_rel_mc_ex[j]->Draw("P SAME");

    leg2->SetHeader("p_{T}-balance, Response, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
    leg2->AddEntry(graph_rel_data_ex[j], "R_{Data, #alpha -> 0}","P");
    leg2->AddEntry(graph_rel_mc_ex[j], "R_{MC, #alpha -> 0}","P");
    leg2->Draw();

    TLatex *tex5 = new TLatex();
    tex5->SetNDC();
    tex5->SetTextSize(0.045); 
    tex5->DrawLatex(0.5,0.91,CorrectionObject::_lumitag+"(13TeV)"); 

    c5->Print(CorrectionObject::_outpath+"plots/Response_Rel_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
    delete c5;
    delete tex5;
    leg2->Clear();
    

    // draw extrapolations  data
    TCanvas *c6 = new TCanvas("c6","",850,700);
    c6->DrawFrame(30,0.5,2100,1.5,(";#bar{p}_{T};Response"));
    gPad->SetLogx();
    graph_mpf_data_ex[j]->SetMarkerStyle(20);
    graph_mpf_data_ex[j]->SetMarkerColor(kBlack);
    graph_mpf_data_ex[j]->SetLineColor(kBlack);
    
    graph_mpf_data_ex[j]->Draw("P");
    
    graph_mpf_mc_ex[j]->SetMarkerStyle(20);
    graph_mpf_mc_ex[j]->SetMarkerColor(kRed);
    graph_mpf_mc_ex[j]->SetLineColor(kRed);
    
    graph_mpf_mc_ex[j]->Draw("P SAME");

    leg2->SetHeader("MPF, Response, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
    leg2->AddEntry(graph_mpf_data_ex[j], "R_{Data, #alpha -> 0}","P");
    leg2->AddEntry(graph_mpf_mc_ex[j], "R_{MC, #alpha -> 0}","P");
    leg2->Draw();

    TLatex *tex6 = new TLatex();
    tex6->SetNDC();
    tex6->SetTextSize(0.045); 
    tex6->DrawLatex(0.5,0.91,CorrectionObject::_lumitag+"(13TeV)"); 

    c6->Print(CorrectionObject::_outpath+"plots/Response_MPF_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
    delete c6;
    delete tex6;
    leg2->Clear();


    // draw extrapolations  data
    TCanvas *c7 = new TCanvas("c7","",850,700);
    c7->DrawFrame(30,0.83,2100,1.27,(";#bar{p}_{T} [GeV];(R^{MC}/R^{Data})_{#alpha -> 0}"));
    gPad->SetLogx();
    graph_ratio_res_rel[j]->SetMarkerStyle(20);
    graph_ratio_res_rel[j]->SetMarkerColor(kBlue);
    graph_ratio_res_rel[j]->SetLineColor(kBlue);
    
    graph_ratio_res_rel[j]->Draw("P");

    leg2->SetHeader("p_{T}-balance, p_{T} extrapolation, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
    leg2->AddEntry(graph_ratio_res_rel[j], "(R^{MC}/R^{Data})_{#alpha -> 0}","P");
    leg2->AddEntry(f1_rel[j], "loglinear fit","L");
    leg2->AddEntry(f2_rel[j], "constant fit","L");
    leg2->Draw();

    TLatex *tex7 = new TLatex();
    tex7->SetNDC();
    tex7->SetTextSize(0.045); 
    tex7->DrawLatex(0.5,0.91,CorrectionObject::_lumitag+"(13TeV)"); 

    TLatex *tex_chi = new TLatex();
    TString chi2_loglin = "loglinear fit #chi^{2}/n.d.f = ";
    chi2_loglin += trunc(f1_rel[j]->GetChisquare());
    chi2_loglin +="/";
    chi2_loglin +=trunc(f1_rel[j]->GetNDF());
    TString chi2_const = "constant fit #chi^{2}/n.d.f = ";
    chi2_const+=trunc(f2_rel[j]->GetChisquare());
    chi2_const+="/";
    chi2_const+=trunc(f2_rel[j]->GetNDF());
    
    tex_chi->SetNDC();
    tex_chi->SetTextSize(0.035); 
    tex_chi->DrawLatex(0.51,0.73,chi2_loglin);
    tex_chi->DrawLatex(0.51,0.69,chi2_const);
    
    c7->Print(CorrectionObject::_outpath+"plots/pTextrapolation_Rel_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
    delete c7;
    delete tex7;
    delete tex_chi;
    leg2->Clear();


    // draw extrapolations  data
    TCanvas *c8 = new TCanvas("c8","",850,700);
    c8->DrawFrame(30,0.83,2100,1.27,(";#bar{p}_{T} [GeV];(R^{Data}/R^{MC})_{#alpha -> 0}"));
    gPad->SetLogx();
    graph_ratio_res_mpf[j]->SetMarkerStyle(20);
    graph_ratio_res_mpf[j]->SetMarkerColor(kBlue);
    graph_ratio_res_mpf[j]->SetLineColor(kBlue);
    
    graph_ratio_res_mpf[j]->Draw("P");

    leg2->SetHeader("MPF, p_{T} extrapolation, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
    leg2->AddEntry(graph_ratio_res_mpf[j], "(R^{Data}/R^{MC})_{#alpha -> 0}","P");
    leg2->Draw();

    TLatex *tex8 = new TLatex();
    tex8->SetNDC();
    tex8->SetTextSize(0.045); 
    tex8->DrawLatex(0.5,0.91,CorrectionObject::_lumitag+"(13TeV)"); 

    TLatex *tex_chi2 = new TLatex();
    TString chi2_loglin1 = "loglinear fit #chi^{2}/n.d.f = ";
    chi2_loglin1 += trunc(f1_mpf[j]->GetChisquare());
    chi2_loglin1 +="/";
    chi2_loglin1 +=trunc(f1_mpf[j]->GetNDF());
    TString chi2_const1 = "constant fit #chi^{2}/n.d.f = ";
    chi2_const1+=trunc(f2_mpf[j]->GetChisquare());
    chi2_const1+="/";
    chi2_const1+=trunc(f2_mpf[j]->GetNDF());
    
    tex_chi2->SetNDC();
    tex_chi2->SetTextSize(0.035); 
    tex_chi2->DrawLatex(0.51,0.73,chi2_loglin);
    tex_chi2->DrawLatex(0.51,0.69,chi2_const);



    c8->Print(CorrectionObject::_outpath+"plots/pTextrapolation_MPF_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
    delete c8;
    delete tex8;
    delete tex_chi2;
    leg2->Clear();

  }
  

  h_chi2_loglin_mpf->GetYaxis()->SetRangeUser(0,20);
  h_chi2_const_mpf ->GetYaxis()->SetRangeUser(0,20);

  h_chi2_loglin_rel->GetYaxis()->SetRangeUser(0,20);
  h_chi2_const_rel ->GetYaxis()->SetRangeUser(0,20);


  TCanvas* c_chi2_loglin_mpf = new TCanvas();
  h_chi2_loglin_mpf->Draw("HIST");
  c_chi2_loglin_mpf ->SaveAs(CorrectionObject::_outpath+"plots/pTextrapolation_MPF_chi2ndf_loglin_"+CorrectionObject::_generator_tag+".pdf");
  delete c_chi2_loglin_mpf;
  delete h_chi2_loglin_mpf;

  TCanvas* c_chi2_const_mpf = new TCanvas();
  h_chi2_const_mpf->Draw("HIST");
  c_chi2_const_mpf ->SaveAs(CorrectionObject::_outpath+"plots/pTextrapolation_MPF_chi2ndf_const_"+CorrectionObject::_generator_tag+".pdf");
  delete c_chi2_const_mpf;
  delete h_chi2_const_mpf;

  TCanvas* c_chi2_loglin_rel = new TCanvas();
  h_chi2_loglin_rel->Draw("HIST");
  c_chi2_loglin_rel ->SaveAs(CorrectionObject::_outpath+"plots/pTextrapolation_Rel_chi2ndf_loglin_"+CorrectionObject::_generator_tag+".pdf");
  delete c_chi2_loglin_rel;
  delete h_chi2_loglin_rel;

  TCanvas* c_chi2_const_rel = new TCanvas();
  h_chi2_const_rel->Draw("HIST");
  c_chi2_const_rel ->SaveAs(CorrectionObject::_outpath+"plots/pTextrapolation_Rel_chi2ndf_const_"+CorrectionObject::_generator_tag+".pdf");
  delete c_chi2_const_rel;
  delete h_chi2_const_rel;




  double flat_mpf[n_eta-1];
  double loglin_mpf[n_eta-1];

  double flat_norm_mpf;
  double loglin_norm_mpf;
  
  double flat_rel[n_eta-1];
  double loglin_rel[n_eta-1];

  double flat_norm_rel;
  double loglin_norm_rel;
  
  flat_norm_mpf=0;
  loglin_norm_mpf=0;
  for (int j=0; j<n_eta-1; j++){
    flat_mpf[j] = f2_mpf[j]->GetParameter(0);
    loglin_mpf[j] = f1_mpf[j]->GetParameter(0)+f1_mpf[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()) ;
  } 
  for (int j=0; j<n_etabarr; j++){
    flat_norm_mpf += flat_mpf[j];
    loglin_norm_mpf += loglin_mpf[j];
  }
  flat_norm_mpf = flat_norm_mpf/n_etabarr;
  loglin_norm_mpf = loglin_norm_mpf/n_etabarr;
  for (int j=0; j<n_eta-1; j++){
      flat_mpf[j] = flat_mpf[j]/flat_norm_mpf;
      loglin_mpf[j] = loglin_mpf[j]/loglin_norm_mpf;
  }
  
  flat_norm_rel=0;
  loglin_norm_rel=0;
  for (int j=0; j<n_eta-1; j++){
    flat_rel[j] = f2_rel[j]->GetParameter(0);
    loglin_rel[j] = f1_rel[j]->GetParameter(0)+f1_rel[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()) ;
  } 
  for (int j=0; j<n_etabarr; j++){
    flat_norm_rel += flat_rel[j];
    loglin_norm_rel += loglin_rel[j];
  }
  flat_norm_rel = flat_norm_rel/n_etabarr;
  loglin_norm_rel = loglin_norm_rel/n_etabarr;
  for (int j=0; j<n_eta-1; j++){
    flat_rel[j] = flat_rel[j]/flat_norm_rel;
    loglin_rel[j] = loglin_rel[j]/loglin_norm_rel;
  }


  TH1D* Residual_logpt_Rel = new TH1D("res_logpt_rel","res_logpt_rel", n_eta-1,eta_bins);
  TH1D* Residual_const_Rel = new TH1D("res_const_rel","res_const_rel", n_eta-1,eta_bins);
  
  TH1D* Residual_logpt_MPF = new TH1D("res_logpt_mpf","res_logpt_mpf", n_eta-1,eta_bins);
  TH1D* Residual_const_MPF = new TH1D("res_const_mpf","res_const_mpf", n_eta-1,eta_bins);
  
  
  for(int k=0; k<5; k++){
    for (int j=0; j<n_eta-1; j++){
      
      double pave_value_corr = 0;
      int reject = 1;
      if(k==0) pave_value_corr = ptave_data[j]->GetMean();
      
      if(k==1) pave_value_corr = 120;
      if(k==2) pave_value_corr = 60;
      if(k==3) pave_value_corr = 240;
      if(k==4) pave_value_corr = 480;
      
  */
      //**************************************  Test reject high pT bins without data ***************************
      /*
      double ptave_value = 0;
      ptave_value = ptave_data[j]->FindLastBinAbove(100.);
      if(k==1 && ptave_value > 120) reject = 1;
      else{ if(k==1) reject = 1000;}
      if(k==2 && ptave_value > 60) reject = 1;
      else{ if(k==2)reject = 1000;}
      if(k==3 && ptave_value > 240) reject = 1;
      else{ if(k==3) reject = 1000;}
      if(k==4 && ptave_value > 480) reject = 1;
      else{if(k==4) reject = 1000;}
      */
      //************************************************************************************************************ 
      /*
      Residual_logpt_Rel->SetBinContent(j+1,(f1_rel[j]->GetParameter(0)+f1_rel[j]->GetParameter(1)*TMath::Log(pave_value_corr))/loglin_norm_rel);
      Residual_logpt_Rel->SetBinError(j+1, sqrt(Vcov_rel[0][j]+Vcov_rel[1][j]*pow(TMath::Log(pave_value_corr),2)+2*Vcov_rel[2][j]*TMath::Log(pave_value_corr))/loglin_norm_rel);
      
      Residual_const_Rel->SetBinContent(j+1,f2_rel[j]->GetParameter(0)/flat_norm_rel);
      Residual_const_Rel->SetBinError(j+1, f2_rel[j]->GetParError(0)/flat_norm_rel);


      Residual_logpt_MPF->SetBinContent(j+1,(f1_mpf[j]->GetParameter(0)+f1_mpf[j]->GetParameter(1)*TMath::Log(pave_value_corr))/loglin_norm_mpf);
      Residual_logpt_MPF->SetBinError(j+1, sqrt(Vcov_mpf[0][j]+Vcov_mpf[1][j]*pow(TMath::Log(pave_value_corr),2)+2*Vcov_mpf[2][j]*TMath::Log(pave_value_corr))/loglin_norm_mpf);
      
      Residual_const_MPF->SetBinContent(j+1,f2_mpf[j]->GetParameter(0)/flat_norm_mpf);
      Residual_const_MPF->SetBinError(j+1, f2_mpf[j]->GetParError(0)/flat_norm_mpf);


    }
    
    //File containing the results with kFSR correction applied
    TFile* outputfile;
    TString variation;
    TFile* outputfile2;

    if(k==0){
      outputfile = new TFile(CorrectionObject::_outpath+"Histo_Res_Rel_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","RECREATE");
    }
    if(k>0){
      if(k==1) variation = "central";
      if(k==2) variation = "down";
      if(k==3) variation = "up";
      if(k==4) variation = "doubleup";
      outputfile = new TFile(CorrectionObject::_outpath+"Histo_Res_Rel_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_"+variation+".root","RECREATE");
    }
    
    Residual_logpt_Rel->Write();
    Residual_const_Rel->Write();
    outputfile->Write();
    outputfile->Close();
    
    delete outputfile;

    if(k==0){
      outputfile2 = new TFile(CorrectionObject::_outpath+"Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","RECREATE");
    }
    if(k>0){
      if(k==1) variation = "central";
      if(k==2) variation = "down";
      if(k==3) variation = "up";
      if(k==4) variation = "doubleup";
      outputfile2 = new TFile(CorrectionObject::_outpath+"Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_"+variation+".root","RECREATE");
    }
    
    Residual_logpt_MPF->Write();
    Residual_const_MPF->Write();
    outputfile2->Write();
    outputfile2->Close();
    
    delete outputfile2;


  }
      */
  cout << "+++++++++++++++++ Finished kFSR() +++++++++++++++++++" << endl;
  
  //delete everything
  
  for(int i=0; i<n_alpha; i++){
    for(int j=0; j<n_eta-1; j++){
      delete hdata_asymmetry[j][i];
      delete hdata_B[j][i];
      delete hmc_asymmetry[j][i];
      delete hmc_B[j][i];
    }
  }
  /*
  for (int j=0; j<n_eta-1; j++) {
    //    delete graph_ratio_res_mpf[j];
    //    delete graph_ratio_res_rel[j];

    delete f1_mpf[j];
    delete f2_mpf[j];
    delete f1_rel[j];
    delete f2_rel[j];
  }
  */
  cout << "++++++++++++ Deleted everything in kFSR(), exiting ++++++++++++++" << endl;
  
}
