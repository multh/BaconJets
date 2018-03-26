// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// with this script pT extrapolations in bins of eta will be calculated
//
// kFSR is aproximated by function and kFSR correction is done by the 
// function values at each eta bin
//
// number of root files with historgrams and txt files for JEC db
// are produced
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "../include/parameters.h"
#include "../include/useful_functions.h"
#include "../include/CorrectionObject.h"
#include "../include/tdrstyle_mod15.h"

#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TMatrixDSym.h>
#include <TPaveStats.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TVirtualFitter.h>
#include <TMath.h>
#include <TFile.h>
#include <TProfile.h>


using namespace std;

void CorrectionObject::Pt_Extrapolation_Alternative_CorrectFormulae(bool mpfMethod){
  cout << "--------------- Starting Pt_Extrapolation() ---------------" << endl << endl;

  TStyle* m_gStyle = new TStyle();
  m_gStyle->SetOptFit(000);

  string output_header = "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*0.021*(-1.+1./(1.+exp(-(TMath::Log(x)-5.030)/0.395))))) Correction L2Relative}";

  // fill the histos for pt average in bins of eta
  TH1D* ptave_data[n_eta-1];
  TH1D* ptave_mc[n_eta-1];

  TH1D* ptave_data_eta = new TH1D("ptave_data_eta", " ;pt_ave;#eta", n_eta-1, eta_bins);
  TH1D* ptave_mc_eta   = new TH1D("ptave_mc_eta", " ;pt_ave; #eta", n_eta-1, eta_bins);

  int countPt = 0;
  TString namehist = "ptave_";
  for(int i=0; i<n_eta-1; i++){
    TString selection = "alpha<";
    selection+=alpha_cut;
    selection+=" && fabs(probejet_eta)<";
    selection+=eta_range[i+1];
    selection+=" && fabs(probejet_eta)>=";
    selection+=eta_range[i];
    TString var1 = "pt_ave";
    ptave_mc[i] = (TH1D*)GetHist(CorrectionObject::_MCFile, selection, var1, 300,0,3000)->Clone();
    ptave_data[i] = (TH1D*)GetHist(CorrectionObject::_DATAFile, selection, var1, 300,0,3000)->Clone();
    
    ptave_data_eta->SetBinContent(i+1,ptave_data[i]->GetMean());
    ptave_mc_eta->SetBinContent(i+1, ptave_mc[i]->GetMean());

    TString namecur = namehist;
    namecur += countPt;
    ptave_data[i]->SetName(namecur);
    countPt++;
  }

  TLine *line_ptave = new TLine(0.,1,eta_bins[n_eta-1],1);

  TCanvas* c_ptave_data = new TCanvas("ptave_data","ptave_data", 700, 500);
  ptave_data_eta->Draw();
  line_ptave->Draw("SAME");
  c_ptave_data -> Print(CorrectionObject::_outpath+"plots/pTave_data.pdf");
  delete c_ptave_data;

  TCanvas* c_ptave_mc = new TCanvas("ptave_mc","ptave_mc", 700, 500);
  ptave_mc_eta->Draw();
  line_ptave->Draw("SAME");
  c_ptave_mc -> Print(CorrectionObject::_outpath+"plots/pTave_mc.pdf");
  delete c_ptave_mc;
  


  //Set up histos for ratios of responses
  double ratio_al_rel_r[n_pt-1][n_eta-1]; //ratio at pt,eta,alpha bins
  double err_ratio_al_rel_r[n_pt-1][n_eta-1]; //error of ratio at pt,eta,alpha bins
  double ratio_al_mpf_r[n_pt-1][n_eta-1]; //ratio at pt,eta,alpha bins
  double err_ratio_al_mpf_r[n_pt-1][n_eta-1]; //error of ratio at pt,eta,alpha bins

  TProfile *pr_data_asymmetry[n_eta-1];// pT-balance response for data  
  TProfile *pr_data_B[n_eta-1];//MPF response for data
  TProfile *pr_mc_asymmetry[n_eta-1];// pT-balanse responce for MC  
  TProfile *pr_mc_B[n_eta-1];//MPF response for MC

  TH2D *hdata_asymmetry[n_eta-1];
  TH2D *hdata_B[n_eta-1];
  TH2D *hmc_asymmetry[n_eta-1];
  TH2D *hmc_B[n_eta-1];

  int n_entries_mc[n_eta-1][n_pt-1];
  int n_entries_data[n_eta-1][n_pt-1];

  TH1D *hdata_ptave[n_pt-1][n_eta-1];//pt-ave in each bin of pT_ave in bins of eta
  
  int count = 0;

  TString name1 = "hist_data_asymmetry_";
  TString name2 = "hist_data_B_";
  TString name3 = "hist_mc_asymmetry_";
  TString name4 = "hist_mc_B_";
  TString name5 = "hist_data_pt_ave";

  for(int j=0; j<n_eta-1; j++){
    TString eta_name = "eta_"+eta_range2[j]+"_"+eta_range2[j+1];
    
    TString name = name1 + eta_name;
    hdata_asymmetry[j] = new TH2D(name,"",n_pt-1,pt_bins,nResponseBins, -1.2, 1.2);
    name = name2 + eta_name;
    hdata_B[j] = new TH2D(name,"",n_pt-1,pt_bins,nResponseBins, -1.2, 1.2);
    name = name3 + eta_name;
    hmc_asymmetry[j] = new TH2D(name,"",n_pt-1,pt_bins,nResponseBins, -1.2, 1.2);
    name = name4 + eta_name;
    hmc_B[j] = new TH2D(name,"",n_pt-1,pt_bins,nResponseBins, -1.2, 1.2);
     
    for(int k=0; k<n_pt-1; k++){
      TString pt_name = "pt_"+pt_range[k]+"_"+pt_range[k+1];
      name = name5 + eta_name + "_" + pt_name; 
      hdata_ptave[k][j] = new TH1D(name,"",3000,0,3000); //only used for GetMean and GetStdDev in the pT extrapolations
      count++;
      ratio_al_rel_r[k][j] = 0;
      err_ratio_al_rel_r[k][j] = 0;
      ratio_al_mpf_r[k][j] = 0;
      err_ratio_al_mpf_r[k][j] = 0;
      n_entries_mc[j][k] = 0;
      n_entries_data[j][k] = 0;
    }
  }  




  //Get relevant information from DATA, loop over DATA events
  TTreeReader myReader_DATA("AnalysisTree", CorrectionObject::_DATAFile);
  TTreeReaderValue<Float_t> pt_ave_data(myReader_DATA, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_data(myReader_DATA, "probejet_eta");
  TTreeReaderValue<Float_t> alpha_data(myReader_DATA, "alpha");
  TTreeReaderValue<Float_t> rel_r_data(myReader_DATA, "rel_r");
  TTreeReaderValue<Float_t> mpf_r_data(myReader_DATA, "mpf_r");
  TTreeReaderValue<Float_t> asymmetry_data(myReader_DATA, "asymmetry");
  TTreeReaderValue<Float_t> B_data(myReader_DATA, "B");   
  TTreeReaderValue<Float_t> weight_data(myReader_DATA, "weight");
   
  while (myReader_DATA.Next()) {
    if(*alpha_data>alpha_cut) continue;
    for(int j=0; j<n_eta-1; j++){
      if(fabs(*probejet_eta_data)>eta_bins[j+1] || fabs(*probejet_eta_data)<eta_bins[j]) continue;
      else{
	hdata_asymmetry[j]->Fill(*pt_ave_data,*asymmetry_data,*weight_data);
	hdata_B[j]->Fill(*pt_ave_data,*B_data,*weight_data);

	for(int k=0; k<n_pt-1; k++){
	  if(*pt_ave_data<pt_bins[k] || *pt_ave_data>pt_bins[k+1]) continue;
	  hdata_ptave[k][j]->Fill(*pt_ave_data,*weight_data);
	  n_entries_data[j][k]++;
	}
      }
    }
  }

  //Get relevant information from MC, loop over MC events 
  TTreeReader myReader_MC("AnalysisTree", CorrectionObject::_MCFile);
  TTreeReaderValue<Float_t> pt_ave_mc(myReader_MC, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_mc(myReader_MC, "probejet_eta");
  TTreeReaderValue<Float_t> alpha_mc(myReader_MC, "alpha");
  TTreeReaderValue<Float_t> rel_r_mc(myReader_MC, "rel_r");
  TTreeReaderValue<Float_t> mpf_r_mc(myReader_MC, "mpf_r");
  TTreeReaderValue<Float_t> asymmetry_mc(myReader_MC, "asymmetry");
  TTreeReaderValue<Float_t> B_mc(myReader_MC, "B");
  TTreeReaderValue<Float_t> weight_mc(myReader_MC, "weight");

  while (myReader_MC.Next()) {
    if(*alpha_mc>alpha_cut) continue;
    for(int j=0; j<n_eta-1; j++){
      if(fabs(*probejet_eta_mc)>eta_bins[j+1] || fabs(*probejet_eta_mc)<eta_bins[j]) continue;
      else{
	hmc_asymmetry[j]->Fill(*pt_ave_mc,*asymmetry_mc,*weight_mc);
	hmc_B[j]->Fill(*pt_ave_mc,*B_mc,*weight_mc);
	for(int k=0; k<n_pt-1; k++){
	  if(*pt_ave_mc<pt_bins[k] || *pt_ave_mc>pt_bins[k+1]) continue;
	  n_entries_mc[j][k]++;
	}
      }
    }
  }


  //Check number of entries in each bin
  bool enough_entries[n_eta-1][n_pt-1];
  for(int j=0; j<n_eta-1; j++){
    for(int k=0; k<n_pt-1; k++){
      enough_entries[j][k] = false;
      if(n_entries_mc[j][k] > 100 && n_entries_data[j][k] > 100) enough_entries[j][k] = true;
    }
  }

  //build profiles out of asymmetry and B 2d-histos to get <A> and <B> as a function of pT in bins of eta,alpha
  for(int j=0; j<n_eta-1; j++){
    //build profiles
    pr_data_asymmetry[j] = (TProfile*)hdata_asymmetry[j]->ProfileX();
    pr_data_B[j]         = (TProfile*)hdata_B[j]->ProfileX();
    pr_mc_asymmetry[j]   = (TProfile*)hmc_asymmetry[j]->ProfileX();
    pr_mc_B[j]           = (TProfile*)hmc_B[j]->ProfileX();
  }
   
  //R_MC and R_DATA overlaid in the same plot as a function of pT, in bins of |eta|
  double val_rel_mc[n_eta-1][n_pt-1]; //value at pt,eta
  double rel_mc_err[n_eta-1][n_pt-1]; //error of ratio at pt,eta
  double val_mpf_mc[n_eta-1][n_pt-1]; //ratio at pt,eta
  double mpf_mc_err[n_eta-1][n_pt-1]; //error of ratio at pt,eta
  double val_rel_data[n_eta-1][n_pt-1]; //value at pt,eta
  double rel_data_err[n_eta-1][n_pt-1]; //error of ratio at pt,eta
  double val_mpf_data[n_eta-1][n_pt-1]; //ratio at pt,eta
  double mpf_data_err[n_eta-1][n_pt-1]; //error of ratio at pt,eta


  //calculate response from <A> and <B> in bins of pt,eta,alpha
  //gaussian error propagation from errors on <A> and <B>
  for(int k=0; k<n_pt-1; k++){
    for(int j=0; j<n_eta-1; j++){

      //responses for data, MC separately. Only for bins with >= 100 entries
      double mpf_mc = (1+pr_mc_B[j]->GetBinContent(k+1))/(1-pr_mc_B[j]->GetBinContent(k+1));
      if(!enough_entries[j][k] || pr_mc_B[j]->GetBinContent(k+1)==0) mpf_mc = 0;
      val_mpf_mc[j][k] = mpf_mc;

      double mpf_data = (1+pr_data_B[j]->GetBinContent(k+1))/(1-pr_data_B[j]->GetBinContent(k+1));
      if(!enough_entries[j][k] || pr_data_B[j]->GetBinContent(k+1)==0) mpf_data = 0;
      val_mpf_data[j][k] = mpf_data;

      double rel_mc = (1+pr_mc_asymmetry[j]->GetBinContent(k+1))/(1-pr_mc_asymmetry[j]->GetBinContent(k+1));
      if(!enough_entries[j][k] || pr_mc_asymmetry[j]->GetBinContent(k+1)==0) rel_mc = 0;
      val_rel_mc[j][k] = rel_mc;

      double rel_data = (1+pr_data_asymmetry[j]->GetBinContent(k+1))/(1-pr_data_asymmetry[j]->GetBinContent(k+1));
      if(!enough_entries[j][k] || pr_data_asymmetry[j]->GetBinContent(k+1)==0) rel_data = 0;
      val_rel_data[j][k] = rel_data;

      double err_mpf_mc = 2/(pow((1-pr_mc_B[j]->GetBinContent(k+1)),2)) * pr_mc_B[j]->GetBinError(k+1);
      if(!enough_entries[j][k] || pr_mc_B[j]->GetBinContent(k+1)==0) err_mpf_mc = 0;
      mpf_mc_err[j][k] = err_mpf_mc;

      double err_mpf_data = 2/(pow((1-pr_data_B[j]->GetBinContent(k+1)),2)) * pr_data_B[j]->GetBinError(k+1);
      if(!enough_entries[j][k] || pr_data_B[j]->GetBinContent(k+1)==0) err_mpf_data = 0;
      mpf_data_err[j][k] = err_mpf_data;

      double err_rel_mc = 2/(pow((1-pr_mc_asymmetry[j]->GetBinContent(k+1)),2)) * pr_mc_asymmetry[j]->GetBinError(k+1);
      if(!enough_entries[j][k] || pr_mc_asymmetry[j]->GetBinContent(k+1)==0) err_rel_mc = 0;
      rel_mc_err[j][k] = err_rel_mc;

      double err_rel_data = 2/(pow((1-pr_data_asymmetry[j]->GetBinContent(k+1)),2)) * pr_data_asymmetry[j]->GetBinError(k+1);
      if(!enough_entries[j][k] || pr_data_asymmetry[j]->GetBinContent(k+1)==0) err_rel_data = 0;
      rel_data_err[j][k] = err_rel_data;

      //ratio of responses, again gaussian error propagation
      if(rel_data > 0) ratio_al_rel_r[k][j] = rel_mc/rel_data;
      else ratio_al_rel_r[k][j] = 0;
      err_ratio_al_rel_r[k][j] = sqrt(pow(1/rel_data*err_rel_mc,2) + pow(rel_mc/(rel_data*rel_data)*err_rel_data,2));
      if(mpf_data > 0) ratio_al_mpf_r[k][j] = mpf_mc/mpf_data;
      else ratio_al_mpf_r[k][j] = 0;
      err_ratio_al_mpf_r[k][j] = sqrt(pow(1/mpf_data*err_mpf_mc,2) + pow(mpf_mc/(mpf_data*mpf_data)*err_mpf_data,2));

    }
  }

  //dummy for tdrCanvas
  TH1D *h = new TH1D("h",";dummy;",41,0,5.191);
  h->SetMaximum(1.2);
  h->SetMinimum(0.8);

  TCanvas* c_0 = new TCanvas();
  tdrCanvas(c_0,"c_0",h,4,10,true,CorrectionObject::_lumitag);


  for(int i=0; i<n_eta-1; i++){
    //Create and fill TGraphErrors
    double xbin_tgraph[n_pt-1];
    double zero[n_pt-1];
    for(int i=0;i<n_pt-1;i++){
      xbin_tgraph[i]=(pt_bins[i]+pt_bins[i+1])/2;
      zero[i]=(pt_bins[i+1]-pt_bins[i])/2 ;
    }
 

    TGraphErrors *graph_mpf_mc   = new TGraphErrors(n_pt-1, xbin_tgraph, val_mpf_mc[i], zero, mpf_mc_err[i]);
    TGraphErrors *graph_mpf_data = new TGraphErrors(n_pt-1, xbin_tgraph, val_mpf_data[i], zero, mpf_data_err[i]);
    TGraphErrors *graph_rel_mc   = new TGraphErrors(n_pt-1, xbin_tgraph, val_rel_mc[i], zero, rel_mc_err[i]);
    TGraphErrors *graph_rel_data = new TGraphErrors(n_pt-1, xbin_tgraph, val_rel_data[i], zero, rel_data_err[i]);
    graph_mpf_mc   = (TGraphErrors*)CleanEmptyPoints(graph_mpf_mc);
    graph_mpf_data = (TGraphErrors*)CleanEmptyPoints(graph_mpf_data);
    graph_rel_mc   = (TGraphErrors*)CleanEmptyPoints(graph_rel_mc);
    graph_rel_data = (TGraphErrors*)CleanEmptyPoints(graph_rel_data);


    graph_mpf_mc->SetTitle("");
    graph_mpf_mc->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph_mpf_mc->GetXaxis()->SetTitleSize(0.05);
    graph_mpf_mc->GetXaxis()->SetTitleOffset(0.80);
    graph_mpf_mc->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    graph_mpf_mc->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_mpf_mc->SetMarkerColor(kRed);
    graph_mpf_mc->SetMarkerStyle(20);
    graph_mpf_mc->SetLineColor(kRed);

    graph_mpf_data->SetTitle("");
    graph_mpf_data->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph_mpf_data->GetXaxis()->SetTitleSize(0.05);
    graph_mpf_data->GetXaxis()->SetTitleOffset(0.80);
    graph_mpf_data->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    graph_mpf_data->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_mpf_data->SetMarkerColor(kBlack);
    graph_mpf_data->SetMarkerStyle(20);
    graph_mpf_data->SetLineColor(kBlack);

    graph_rel_mc->SetTitle("");
    graph_rel_mc->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph_rel_mc->GetXaxis()->SetTitleSize(0.05);
    graph_rel_mc->GetXaxis()->SetTitleOffset(0.80);
    graph_rel_mc->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    graph_rel_mc->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_rel_mc->SetMarkerColor(kRed);
    graph_rel_mc->SetMarkerStyle(20);
    graph_rel_mc->SetLineColor(kRed);

    graph_rel_data->SetTitle("");
    graph_rel_data->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph_rel_data->GetXaxis()->SetTitleSize(0.05);
    graph_rel_data->GetXaxis()->SetTitleOffset(0.80);
    graph_rel_data->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    graph_rel_data->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_rel_data->SetMarkerColor(kBlack);
    graph_rel_data->SetMarkerStyle(20);
    graph_rel_data->SetLineColor(kBlack);



    TString alVal;
    alVal.Form("%0.2f\n",alpha_cut);
    TString altitle = "{#alpha<"+alVal+"}";
    TString axistitle_mc = "R^{MC}_";
    TString axistitle_data = "R^{Data}_";
    axistitle_mc   += altitle;
    axistitle_data += altitle;

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 



    TCanvas* c_mpf = new TCanvas();
    tdrCanvas(c_mpf,"c_mpf",h,4,10,true,CorrectionObject::_lumitag);
    h->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.80);
    h->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    h->GetYaxis()->SetRangeUser(0.50,1.50);
    graph_mpf_mc->Draw("P SAME");
    graph_mpf_data->Draw("P SAME");
    gPad->SetLogx();

    TLegend *leg_mpf;
    leg_mpf = new TLegend(0.35,0.72,0.51,0.92,"","brNDC");//x+0.1
    leg_mpf->SetBorderSize(0);
    leg_mpf->SetTextSize(0.038);
    leg_mpf->SetFillColor(10);
    leg_mpf->SetFillStyle(0);
    leg_mpf->SetLineColor(1);
    leg_mpf->SetTextFont(42);
    leg_mpf->SetHeader("MPF response, "+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]+", #alpha<"+s_alpha_cut);
    leg_mpf->AddEntry(graph_mpf_mc, axistitle_mc,"P");
    leg_mpf->AddEntry(graph_mpf_data, axistitle_data,"P");
    leg_mpf->Draw();

    //tex->DrawLatex(0.53,0.91,CorrectionObject::_lumitag+"(13TeV)");

    c_mpf->SaveAs(CorrectionObject::_outpath+"plots/control/MPF_Response_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + ".pdf");


    TCanvas* c_rel = new TCanvas();
    tdrCanvas(c_rel,"c_rel",h,4,10,true,CorrectionObject::_lumitag);
    h->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.80);
    h->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    h->GetYaxis()->SetRangeUser(0.50,1.50);
    graph_rel_mc->Draw("P SAME");
    graph_rel_data->Draw("P SAME");
    gPad->SetLogx();

    TLegend *leg_rel;
    leg_rel = new TLegend(0.35,0.72,0.51,0.92,"","brNDC");//x+0.1
    leg_rel->SetBorderSize(0);
    leg_rel->SetTextSize(0.038);
    leg_rel->SetFillColor(10);
    leg_rel->SetFillStyle(0);
    leg_rel->SetLineColor(1);
    leg_rel->SetTextFont(42);
    leg_rel->SetHeader("p_{T}-balance response, "+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]); 
    leg_rel->AddEntry(graph_rel_mc, axistitle_mc,"P");
    leg_rel->AddEntry(graph_rel_data, axistitle_data,"P");
    leg_rel->Draw();

    //tex->DrawLatex(0.53,0.91,CorrectionObject::_lumitag+"(13TeV)");

    c_rel->SaveAs(CorrectionObject::_outpath+"plots/control/Rel_Response_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + ".pdf");


    //delete leg_rel;
    delete c_rel;
    //delete leg_mpf;
    delete c_mpf;
    delete tex;
    delete graph_rel_data;
    delete graph_rel_mc;
    delete graph_mpf_data;
    delete graph_mpf_mc;
  }


  // get ratio for MC to DATA responses
  //BAD NOMENCLATURE
  double ratio_mpf[n_eta-1][n_pt-1];     //ratio at pt,eta bins for alpha = 0.3
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


  //Create and fill TGraphErrors
  // double xbin_tgraph[n_pt-1];
  // double zero[n_pt-1];
  double xbin_tgraph[n_eta-1][n_pt-1];
  double zero[n_eta-1][n_pt-1];
  for(int i=0;i<n_pt-1;i++){
    for(int j=0; j<n_eta-1; j++){
   
      xbin_tgraph[j][i] = hdata_ptave[i][j]->GetMean();
      zero[j][i] = hdata_ptave[i][j]->GetStdDev();
    }
  }
  TGraphErrors *graph1_mpf[n_eta-1];
 
  bool graph_filled[n_eta-1];
  for(int j=0; j<n_eta-1; j++) graph_filled[j] = false;
  for(int j=0; j<n_eta-1; j++){
    graph1_mpf[j] = new TGraphErrors(n_pt-1, xbin_tgraph[j], ratio_mpf[j], zero[j], err_ratio_mpf[j]);
    graph1_mpf[j] = (TGraphErrors*)CleanEmptyPoints(graph1_mpf[j]); 
    if(graph1_mpf[j]->GetN() > 0) graph_filled[j] = true;
    // cout << "Eta bin no: " << j << ", graph filled?: " << graph_filled[j] << endl;
  }

  for(int k=0; k<n_pt-1; k++){
    double x = 0;
    double val = 0;
    graph1_mpf[13]->GetPoint(k,x,val);
 }

  //Mikko's request: delete super-high data point in 2.8-2.9 bin (j==13) above ~400 GeV pT (point no 6 & 7 [c++])
// <<<<<<< HEAD

   if(graph1_mpf[13]->GetN() == 9){
     graph1_mpf[13]->RemovePoint(8);
     graph1_mpf[13]->RemovePoint(7);
    }
    else if(graph1_mpf[13]->GetN() == 8){
     graph1_mpf[13]->RemovePoint(7);
    }
    if(graph1_mpf[14]->GetN() == 9){
       graph1_mpf[14]->RemovePoint(8);
    }

    if(graph1_mpf[15]->GetN() == 8){
      graph1_mpf[15]->RemovePoint(7);
    }
    
    if(graph1_mpf[16]->GetN() == 6){
      graph1_mpf[16]->RemovePoint(5);
    }
    

  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,pt_bins[n_pt-1]+10,1);

  // create a function for the loglinear fit
  TF1 * f1[n_eta-1];
  // create a function for the constant fit
  TF1 * f2[n_eta-1];



  // create output .dat file, including the kFSR extrapolation (alpha->0)
  FILE *fp, *fp2, *l2resfile;
  if(mpfMethod){
    fp = fopen(CorrectionObject::_outpath+"output/pT_MPF_"+CorrectionObject::_generator_tag+"_extrapolations.dat","w");
    fp2 = fopen(CorrectionObject::_outpath+"output/pT_MPF_"+CorrectionObject::_generator_tag+"_constantExtrapolation.dat","w");
    l2resfile = fopen(CorrectionObject::_outpath+"output/L2Res_MPF_"+CorrectionObject::_generator_tag+".dat","w");
  }
  else{
    fp = fopen(CorrectionObject::_outpath+"output/pT_DiJet_"+CorrectionObject::_generator_tag+"_extrapolations.dat","w");
    fp2 = fopen(CorrectionObject::_outpath+"output/pT_DiJet_"+CorrectionObject::_generator_tag+"_constantExtrapolation.dat","w");
    l2resfile = fopen(CorrectionObject::_outpath+"output/L2Res_DiJet_"+CorrectionObject::_generator_tag+".dat","w");
  }



  /* +++++++++++++++++++++++ Plots +++++++++++++++++++ */
  if(mpfMethod){
    TFile* pT_extrapolation_mpf_out = new TFile(CorrectionObject::_outpath+"plots/pT_extrapolation_mpf.root","RECREATE");}
  else{   TFile* pT_extrapolation_pt_out = new TFile(CorrectionObject::_outpath+"plots/pT_extrapolation_pt.root","RECREATE");}

  TCanvas* asd[n_eta-1];
  TString plotname[n_eta-1];
  double Vcov[3][n_eta-1];//covariance matrix for log lin fit results
  TH1D* h_chi2_loglin = new TH1D("h_chi2_loglin", "Chi2/ndf for each eta bin;|#eta|;#chi^{2}/ndf", n_eta-1, eta_bins);
  TH1D* h_chi2_const =  new TH1D("h_chi2_const", "Chi2/ndf for each eta bin;|#eta|;#chi^{2}/ndf", n_eta-1, eta_bins);
  
  for (int j=0; j<n_eta-1; j++){
    if(mpfMethod){
      plotname[j]="mpf_ptextra_"+CorrectionObject::_generator_tag+"_eta_"+eta_range[j]+"_"+eta_range[j+1];
    }
    else{
      plotname[j]="dijet_ptextra_"+CorrectionObject::_generator_tag+"_eta_"+eta_range[j]+"_"+eta_range[j+1];
    }
    asd[j] = new TCanvas(plotname[j], plotname[j], 800,700);
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(0);
    gPad->SetLogx();
    graph1_mpf[j]->SetMarkerColor(kBlue);
    graph1_mpf[j]->SetMarkerStyle(20);
    graph1_mpf[j]->SetLineColor(kBlue);


    //Log Linear/Constant fit for pT extrapolation
    f1[j] = new TF1(plotname[j]+"f1","[0]+[1]*TMath::Log(x)", 50 , pt_bins[n_pt-1]+10);
    f1[j]->SetParameters(1,0);

    f2[j] = new TF1(plotname[j]+"f2","pol0", 50 , pt_bins[n_pt-1]+10);
    f2[j]->SetLineColor(kBlue);
    f2[j]->SetLineStyle(3);
   
    //Do the fit!
    if(graph_filled[j]){
    TFitResultPtr fitloglin;
    
     if(j==17){
      cout<<"Fixed slope parameter: "<< f1[j-1]->GetParameter(1.) <<endl;
      double slope = f1[j-1]->GetParameter(1.);
      f1[j]->FixParameter(1,slope); //fixing slope for last eta bin (Vidyo Meeting 07.07.2017)
      }
    
      fitloglin =  graph1_mpf[j]->Fit(plotname[j]+"f1","SM","",55,1200);
      TMatrixDSym cov = fitloglin->GetCovarianceMatrix();
      Vcov[0][j] = cov(0,0);
      Vcov[1][j] = cov(1,1);     
      Vcov[2][j] = cov(0,1);
    }
    else{
      f1[j]->SetParameters(0.0001,0.0001);
      f1[j]->SetParError(0,0.0001);
      f1[j]->SetParError(1,0.0001);
      Vcov[0][j] = 0.0001;
      Vcov[1][j] = 0.0001;     
      Vcov[2][j] = 0.0001;
    }
    
    if(graph_filled[j]){
      graph1_mpf[j]->Fit(plotname[j]+"f2","SM + SAME","",55,1200); }
    else{
      f2[j]->SetParameter(0,0.0001);
      f2[j]->SetParError(0,0.0001);
    }
   
    graph1_mpf[j]->Draw("AP");
    if(graph_filled[j]){
      graph1_mpf[j]->SetTitle("");
      graph1_mpf[j]->GetYaxis()->SetTitle("(R_{MC}/R_{DATA})_{#alpha<0.3}");
      graph1_mpf[j]->GetYaxis()->SetTitleSize(0.045);
      graph1_mpf[j]->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
      graph1_mpf[j]->GetXaxis()->SetTitleSize(0.045);
      graph1_mpf[j]->GetXaxis()->SetTitleOffset(0.82);
      graph1_mpf[j]->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
      graph1_mpf[j]->GetYaxis()->SetRangeUser(0.85,1.25);
    }
 
    line->SetLineStyle(2);
    line->Draw("SAME");
 
    //Store the p0 parameter from const fit
    if (fp2!=NULL) {
      // getting the p0 parameter from the constant fit
      Float_t value = f2[j]->GetParameter(0);
      fprintf(fp2, "%f\n",value);
    }
    if (l2resfile!=NULL) {
      fprintf(l2resfile, "%f\n", eta_bins[j]);
    }

    //Cosmetics
    asd[j]->Modified(); 
    asd[j]->Update();
   
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
    asd[j]->Modified(); 
    asd[j]->Update();


    //Set up legend, prepare canvas layout
    TLegend *leg1;
    leg1 = new TLegend(0.2,0.68,0.35,0.88,"","brNDC");//x+0.1
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.039);
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
    alVal.Form("%0.2f\n",alpha_cut);
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
    tex->DrawLatex(0.5,0.96,CorrectionObject::_lumitag+"(13TeV)"); 

    TLatex *tex2 = new TLatex();
    if(graph_filled[j]){    
      TString chi2_loglin = "loglinear fit #chi^{2}/n.d.f = ";
      chi2_loglin += trunc(f1[j]->GetChisquare());
      chi2_loglin +="/";
      chi2_loglin +=trunc(f1[j]->GetNDF());
      TString chi2_const = "constant fit #chi^{2}/n.d.f = ";
      chi2_const+=trunc(f2[j]->GetChisquare());
      chi2_const+="/";
      chi2_const+=trunc(f2[j]->GetNDF());

      tex2->SetNDC();
      tex2->SetTextSize(0.035); 
      tex2->DrawLatex(0.51,0.73,chi2_loglin);
      tex2->DrawLatex(0.51,0.69,chi2_const);

      double chi2ndf_loglin = f1[j]->GetChisquare() / f1[j]->GetNDF();
      double chi2ndf_const =  f2[j]->GetChisquare() / f2[j]->GetNDF();
      h_chi2_loglin->SetBinContent(j+1,chi2ndf_loglin); //to make sure to fill the right eta-bin...
      h_chi2_const->SetBinContent(j+1,chi2ndf_const);   //to make sure to fill the right eta-bin...
    }

    h_chi2_loglin->GetYaxis()->SetRangeUser(0,20);
    h_chi2_const->GetYaxis()->SetRangeUser(0,20);
 

   //Store plots
   if(mpfMethod){
      graph1_mpf[j]->Write("pTextrapolation_MPF_"+CorrectionObject::_generator_tag+"_pT_"+eta_range2[j]+"_"+eta_range2[j+1]);
       asd[j]->Print(CorrectionObject::_outpath+"plots/pTextrapolation_MPF_"+CorrectionObject::_generator_tag+"_pT_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
    }
    else{
      graph1_mpf[j]->Write("pTextrapolation_Pt_"+CorrectionObject::_generator_tag+"_pT_"+eta_range2[j]+"_"+eta_range2[j+1]);
       asd[j]->Print(CorrectionObject::_outpath+"plots/pTextrapolation_Pt_"+CorrectionObject::_generator_tag+"_pT_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
    }
    delete tex2;
    delete tex;
    delete leg1;
 
  }
  
  TCanvas* c_chi2_loglin = new TCanvas();
  h_chi2_loglin->Draw("HIST");
  if(mpfMethod)c_chi2_loglin->SaveAs(CorrectionObject::_outpath+"plots/pTextrapolation_MPF_chi2ndf_loglin_"+CorrectionObject::_generator_tag+".pdf");
  else c_chi2_loglin->SaveAs(CorrectionObject::_outpath+"plots/pTextrapolation_Pt_chi2ndf_loglin_"+CorrectionObject::_generator_tag+".pdf");
  delete c_chi2_loglin;
  delete h_chi2_loglin;

  TCanvas* c_chi2_const = new TCanvas();
  h_chi2_const->Draw("HIST");
  if(mpfMethod) c_chi2_const->SaveAs(CorrectionObject::_outpath+"plots/pTextrapolation_MPF_chi2ndf_const_"+CorrectionObject::_generator_tag+".pdf");
  else c_chi2_const->SaveAs(CorrectionObject::_outpath+"plots/pTextrapolation_Pt_chi2ndf_const_"+CorrectionObject::_generator_tag+".pdf");
  delete c_chi2_const;
  delete h_chi2_const;

  fclose(fp);
  fclose(fp2);


  /* ++++++++++++++++++++++++++ Calculate L2Residuals MPF ++++++++++++++++++++++++++++++ */

 
  // get the kFSR file
  TCanvas* c_kfsr_fit = new TCanvas("c_kfsr_fit", "c_kfsr_fit",1);
  c_kfsr_fit->Update();
  if(mpfMethod){
    TFile* kfsr_mpf;
    TH1D* hist_kfsr_fit_mpf;
    TH1D* hist_kfsr_mpf;
    
    if(CorrectionObject::_runnr != "BCDEFGH"){
        kfsr_mpf = new TFile(CorrectionObject::_input_path+"RunBCDEFGH/Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_AK4PFchs.root","READ");
	hist_kfsr_fit_mpf = (TH1D*)kfsr_mpf->Get("hist_kfsr_fit_mpf");
	hist_kfsr_mpf = (TH1D*)kfsr_mpf->Get("kfsr_mpf");
    }
    
     else{
    kfsr_mpf = new TFile(CorrectionObject::_outpath+"Histo_KFSR_MPF_"+CorrectionObject::_generator_tag+"_L1.root","READ");
    hist_kfsr_mpf = (TH1D*)kfsr_mpf->Get("kfsr_mpf");

    //fit the kFSR values
    TF1 *kfsr_fit_mpf = new TF1("kfsr_fit_mpf","[0]+([1]*TMath::CosH(x))/(1+[2]*TMath::CosH(x))",0,5.);   //Range: 0,5. by default

    bool fit_fullrange = false;
    bool fit_285       = false;
 
    //VERY fragile fit, carefully set initial values
    if(CorrectionObject::_generator == "pythia"){

      if(CorrectionObject::_collection == "AK4CHS"){
	if(CorrectionObject::_runnr == "BCDEFGH"){
	  if(!CorrectionObject::_closuretest) kfsr_fit_mpf->SetParameters(0.9,3.9,110); //RES
	  else kfsr_fit_mpf->SetParameters(1,2,50); //CLOSURETEST
	  fit_fullrange = true;
	}
      }
    }    
    else throw runtime_error("PTextrapolation, MPF kFSR-fit: Invalid generator specified.");

    //Finally perform the fit!
    //Be carefull with the fit range!
    kfsr_fit_mpf->SetLineColor(kRed+1);
    std::cout<<"!!! kFSR MPF fit !!!"<<std::endl;

    if(!fit_fullrange && !fit_285) hist_kfsr_mpf->Fit("kfsr_fit_mpf","SR","",0,3.14);
    //   if(!fit_fullrange) hist_kfsr_mpf->Fit("kfsr_fit_mpf","SR","",0,2.4);
    else if(!fit_fullrange && fit_285) hist_kfsr_mpf->Fit("kfsr_fit_mpf","SR","",0,2.85);
    else hist_kfsr_mpf->Fit("kfsr_fit_mpf","SR","",0,5.19);

    //Create a histogram to hold the confidence intervals                                                                                        
    hist_kfsr_fit_mpf = (TH1D*)hist_kfsr_mpf->Clone();
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hist_kfsr_fit_mpf);
    //Now the "hist_kfsr_fit" histogram has the fitted function values as the
    //bin contents and the confidence intervals as bin errors

    hist_kfsr_fit_mpf->SetStats(kFALSE);
    hist_kfsr_fit_mpf->SetFillColor(kRed-10);     
    hist_kfsr_fit_mpf->SetName("hist_kfsr_fit_mpf");
    hist_kfsr_fit_mpf->SetTitle("kfsr fit for mpf");
     }
    //  }



    double flat[n_eta-1];
    double loglin[n_eta-1];

    double flat_var[n_eta-1];
    double loglin_var[n_eta-1];

    double flat_norm;
    double loglin_norm;

    double flat_norm_var;
    double loglin_norm_var;


    for (int j=0; j<n_eta-1; j++){
      flat[j] = hist_kfsr_fit_mpf->GetBinContent(j+1)*f2[j]->GetParameter(0);
      loglin[j] = hist_kfsr_fit_mpf->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()) );
    }
    flat_norm=0;
    loglin_norm=0;
    for (int j=0; j<n_etabarr; j++){
      flat_norm += flat[j];
      loglin_norm += loglin[j];
    }
    flat_norm = flat_norm/n_etabarr;
    loglin_norm = loglin_norm/n_etabarr;
    for (int j=0; j<n_eta-1; j++){
      flat[j] = flat[j]/flat_norm;
      loglin[j] = loglin[j]/loglin_norm;
    }

    for (int j=0; j<n_eta-1; j++){
      flat_var[j] = hist_kfsr_mpf->GetBinContent(j+1)*f2[j]->GetParameter(0);
      loglin_var[j] = hist_kfsr_mpf->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()) );
    }

    flat_norm_var=0;
    loglin_norm_var=0;

    for (int j=0; j<n_etabarr; j++){
      flat_norm_var += flat_var[j];
      loglin_norm_var += loglin_var[j];
    }
    flat_norm_var = flat_norm_var/n_etabarr;
    loglin_norm_var = loglin_norm_var/n_etabarr;
    for (int j=0; j<n_eta-1; j++){
      flat_var[j] = flat_var[j]/flat_norm_var;
      loglin_var[j] = loglin_var[j]/loglin_norm_var;
    }
 

    //Histograms to hold the output
    TH1D* Residual_logpt_MPF = new TH1D("res_logpt_mpf","res_logpt_mpf", n_eta-1,eta_bins);
    TH1D* Residual_const_MPF = new TH1D("res_const_mpf","res_const_mpf", n_eta-1,eta_bins);
    TH1D* ptave_const_MPF    = new TH1D("ptave_const_mpf","ptave_const_mpf", n_eta-1,eta_bins);
    TH1D* ptave_logpt_MPF    = new TH1D("ptave_logpt_mpf","ptave_logpt_mpf", n_eta-1,eta_bins);

    TH1D* Residual_logpt_MPF_val = new TH1D("res_logpt_mpf_val","res_logpt_mpf_val", n_eta-1,eta_bins);
    TH1D* Residual_const_MPF_val = new TH1D("res_const_mpf_val","res_const_mpf_val", n_eta-1,eta_bins);
    
    string kFSR_method[2] = {"kFSR_Fit", "kFSR_Hist"};    

    for(int i=0; i<2;i++){ 

      ofstream output, output_loglin, uncerts, uncerts_loglin, output_hybrid, uncerts_hybrid, output_hybrid_barrel, uncerts_hybrid_barrel;
      output.open(CorrectionObject::_outpath+"output/"+kFSR_method[i]+"/Summer16_07Aug2017_MPF_FLAT_L2Residual_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt");
      output_loglin.open(CorrectionObject::_outpath+"output/"+kFSR_method[i]+"/Summer16_07Aug2017_MPF_LOGLIN_L2Residual_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt");
      uncerts.open(CorrectionObject::_outpath+"output/"+kFSR_method[i]+"/Summer16_07Aug2017_MPF_FLAT_L2Residual_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt.STAT");
      uncerts_loglin.open(CorrectionObject::_outpath+"output/"+kFSR_method[i]+"/Summer16_07Aug2017_MPF_LOGLIN_L2Residual_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt.STAT");
      
      output_hybrid.open(CorrectionObject::_outpath+"output/"+kFSR_method[i]+"/Summer16_07Aug2017_MPF_Hybrid_L2Residual_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt");
      uncerts_hybrid.open(CorrectionObject::_outpath+"output/"+kFSR_method[i]+"/Summer16_07Aug2017_MPF_Hybrid_L2Residual_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt.STAT");
      
      output_hybrid_barrel.open(CorrectionObject::_outpath+"output/"+kFSR_method[i]+"/Summer16_07Aug2017_MPF_Hybrid_Barrel_L2Residual_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt");
      uncerts_hybrid_barrel.open(CorrectionObject::_outpath+"output/"+kFSR_method[i]+"/Summer16_07Aug2017_MPF_Hybrid_Barrel_L2Residual_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt.STAT");
      
      //   output  <<"{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*0.021*(-1.+1./(1.+TMath::Exp(-(TMath::Log(x)-5.030)/0.395))))) Correction L2Relative}" << endl;  //Old Header (Stand: 23.Jan2018)


      output  <<"{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*0.021*(-1.+1./(1.+exp(-(TMath::Log(x)-5.030)/0.395))))) Correction L2Relative}" << endl;
      
      output_loglin  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*0.021*(-1.+1./(1.+exp(-(TMath::Log(x)-5.030)/0.395))))) Correction L2Relative}" << endl;
      uncerts << "{ 1 JetEta 1 JetPt [0] kFSR_err f0_err}" << endl;
      //uncerts_loglin << "{ 1 JetEta 1 JetPt [0] kFSR_err cov(0,0) cov(1,1) cov(0,1) }" << endl;
      uncerts_loglin << "{ 1 JetEta 1 JetPt sqrt(fabs([0]*[0]+[1]+[2]*log(x)*log(x)+2*[3]*log(x))) Correction L2Relative}" << endl;
      
      output_hybrid  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*0.021*(-1.+1./(1.+exp(-(TMath::Log(x)-5.030)/0.395))))) Correction L2Relative}" << endl;

      uncerts_hybrid << "{ 1 JetEta 1 JetPt Hybrid Method: Use const fit for 2.65 < |eta| < 3.1, else Log-lin}" << endl;
      
      output_hybrid_barrel  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*0.021*(-1.+1./(1.+exp(-(TMath::Log(x)-5.030)/0.395))))) Correction L2Relative}" << endl;

      uncerts_hybrid_barrel << "{ 1 JetEta 1 JetPt Hybrid Method: Use const fit for 2.65 < |eta| < 3.1, else Log-lin}" << endl;
      

      double norm_flat=0;
      double norm_loglin=0;

      TH1D*  kfsr_FitVal_mpf;
      //!!!Check if fit or hist values are used!!!
      if(kFSR_method[i]=="kFSR_Fit"){
	cout<<"!kFSR: Fit Values used for txt-files and plots!"<<endl;
	kfsr_FitVal_mpf = (TH1D*) hist_kfsr_fit_mpf->Clone();
	norm_flat = flat_norm;
	norm_loglin = loglin_norm;
      }
      else{
	cout<<"!kFSR: Hist Values used for txt-files and plots!"<<endl;
	kfsr_FitVal_mpf = (TH1D*) hist_kfsr_mpf->Clone();
	norm_flat = flat_norm_var;
	norm_loglin = loglin_norm_var;
      }
      
      
      for (int j=n_eta-1; j>0; --j){
	
	output << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   30   " << ptave_data[j-1]->FindLastBinAbove(100.)*10 << "   " << 1/norm_flat << " " << kfsr_FitVal_mpf->GetBinContent(j) << "   " << f2[j-1]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;
	
	output_loglin << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   30   " << ptave_data[j-1]->FindLastBinAbove(100.)*10 << "   " << 1/norm_loglin << " " << kfsr_FitVal_mpf->GetBinContent(j) << "   " << f1[j-1]->GetParameter(0) << " " << f1[j-1]->GetParameter(1) << "   1 0.0000 0.0" << endl;
	
	uncerts << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " 
		<< ptave_data[j-1]->FindLastBinAbove(100.)*10 << "  " << kfsr_FitVal_mpf->GetBinError(j-1) / norm_flat 
		<<" "<< f2[j-1]->GetParError(0) << endl;
	
	uncerts_loglin << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" 
		       << eta_range[j-1] << "  6    30    " << ptave_data[j-1]->FindLastBinAbove(100.)*10 
		       << " " << kfsr_FitVal_mpf->GetBinError(j-1) / norm_loglin  
		       << " " <<Vcov[0][j-1] <<" "<<Vcov[1][j-1]<<" "<<Vcov[2][j-1]<<endl; 
	
	if(j>12 && j<16){
	  output_hybrid << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   30   " << ptave_data[j-1]->FindLastBinAbove(100.)*10 << "   " << 1/norm_flat << " " << kfsr_FitVal_mpf->GetBinContent(j) << "   " << f2[j-1]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;
	  
	  uncerts_hybrid << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    30    " 
			 << ptave_data[j-1]->FindLastBinAbove(100.)*10 << "  " << kfsr_FitVal_mpf->GetBinError(j-1) / norm_flat 
			 <<" "<< f2[j-1]->GetParError(0) << endl;
	}
	else{
	  output_hybrid <<  fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   30   " << ptave_data[j-1]->FindLastBinAbove(100.)*10 << "   " << 1/norm_loglin << " " << kfsr_FitVal_mpf->GetBinContent(j) << "   " << f1[j-1]->GetParameter(0) << " " << f1[j-1]->GetParameter(1) << "   1 0.0000 0.0" << endl;
	  
	  uncerts_hybrid << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" 
			 << eta_range[j-1] << "  6    30    " << ptave_data[j-1]->FindLastBinAbove(100.)*10 
			 << " " << kfsr_FitVal_mpf->GetBinError(j-1) / norm_loglin  
			 << " " <<Vcov[0][j-1] <<" "<<Vcov[1][j-1]<<" "<<Vcov[2][j-1]<<endl; 
	}
	
	
	if(j<6){
	  output_hybrid_barrel << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   30   " << ptave_data[j-1]->FindLastBinAbove(100.)*10 << "   " << 1/norm_flat << " " << kfsr_FitVal_mpf->GetBinContent(j) << "   " << f2[j-1]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;
	  
	  uncerts_hybrid_barrel << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    30    " 
				<< ptave_data[j-1]->FindLastBinAbove(100.)*10 << "  " << kfsr_FitVal_mpf->GetBinError(j-1) / norm_flat 
				<<" "<< f2[j-1]->GetParError(0) << endl;
	}
	else{
	  output_hybrid_barrel <<  fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   30   " << ptave_data[j-1]->FindLastBinAbove(100.)*10 << "   " << 1/norm_flat << " " << kfsr_FitVal_mpf->GetBinContent(j) << "   " << f1[j-1]->GetParameter(0) << " " << f1[j-1]->GetParameter(1) << "   1 0.0000 0.0" << endl;
	  
	  uncerts_hybrid_barrel << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" 
				<< eta_range[j-1] << "  6    30    " << ptave_data[j-1]->FindLastBinAbove(100.)*10 
				<< " " << kfsr_FitVal_mpf->GetBinError(j-1) / norm_flat  
				<< " " <<Vcov[0][j-1] <<" "<<Vcov[1][j-1]<<" "<<Vcov[2][j-1]<<endl; 
	}
      }
      
      
      for (int j=0; j<n_eta-1; j++){
	
	
	output << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   30   " 
	       << ptave_data[j]->FindLastBinAbove(100.)*10 << "   " << 1/norm_flat << " " << kfsr_FitVal_mpf->GetBinContent(j+1) 
	       << "   " << f2[j]->GetParameter(0) << " 0   1 0.0000 0.0"  << endl;
	
	output_loglin << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   30   " 
		      << ptave_data[j]->FindLastBinAbove(100.)*10 << "   " << 1/norm_loglin << " " 
		      << kfsr_FitVal_mpf->GetBinContent(j+1) << "   " << f1[j]->GetParameter(0) 
		      << " " << f1[j]->GetParameter(1) << "   1 0.0000 0.0" << endl;
	
	uncerts << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    30    " 
		<< ptave_data[j]->FindLastBinAbove(100.)*10 
		<<" "<< kfsr_FitVal_mpf->GetBinError(j)/norm_flat<< " " << f2[j]->GetParameter(0) <<endl;
	
	uncerts_loglin << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  6    30    " 
		       << ptave_data[j]->FindLastBinAbove(100.)*10 << " " 
		       << kfsr_FitVal_mpf->GetBinError(j)/norm_loglin<< " "<< Vcov[0][j]<<" "<< Vcov[1][j]<<" "<< Vcov[2][j] << endl; 
	
	if(j>11 && j<15){
	  output_hybrid << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   30   " 
			<< ptave_data[j]->FindLastBinAbove(100.)*10 << "   " << 1/norm_flat << " " << kfsr_FitVal_mpf->GetBinContent(j+1) 
			<< "   " << f2[j]->GetParameter(0) << " 0   1 0.0000 0.0"  << endl;
	  
	  uncerts_hybrid << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    30    " 
			 << ptave_data[j]->FindLastBinAbove(100.)*10 
			 <<" "<< kfsr_FitVal_mpf->GetBinError(j)/norm_flat<< " " << f2[j]->GetParameter(0) <<endl;
	}
	else{
	  output_hybrid << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   30   " 
			<< ptave_data[j]->FindLastBinAbove(100.)*10 << "   " << 1/norm_loglin << " " 
			<< kfsr_FitVal_mpf->GetBinContent(j+1) << "   " << f1[j]->GetParameter(0) 
			<< " " << f1[j]->GetParameter(1) << "   1 0.0000 0.0"<< endl; 
	  
	  uncerts_hybrid << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  6    30    " 
			 << ptave_data[j]->FindLastBinAbove(100.)*10 << " " 
			 << kfsr_FitVal_mpf->GetBinError(j)/norm_loglin<< " "<< Vcov[0][j]<<" "<< Vcov[1][j]<<" "<< Vcov[2][j] << endl; 
	}
	
	if(j<5){
	  output_hybrid_barrel << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   30   " 
			       << ptave_data[j]->FindLastBinAbove(100.)*10 << "   " << 1/norm_flat << " " << kfsr_FitVal_mpf->GetBinContent(j+1) 
			       << "   " << f2[j]->GetParameter(0) << " 0   1 0.0000 0.0"  << endl;
	  
	  uncerts_hybrid_barrel << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    30    " 
				<< ptave_data[j]->FindLastBinAbove(100.)*10 
				<<" "<< kfsr_FitVal_mpf->GetBinError(j)/norm_flat<< " " << f2[j]->GetParameter(0) <<endl;
	}
	else{
	  output_hybrid_barrel << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   30   " 
			       << ptave_data[j]->FindLastBinAbove(100.)*10 << "   " << 1/norm_flat << " " 
			       << kfsr_FitVal_mpf->GetBinContent(j+1) << "   " << f1[j]->GetParameter(0) 
			       << " " << f1[j]->GetParameter(1) << "   1 0.0000 0.0"<< endl; 
	  
	  uncerts_hybrid_barrel << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  6    30    " 
				<< ptave_data[j]->FindLastBinAbove(100.)*10 << " " 
				<< kfsr_FitVal_mpf->GetBinError(j)/norm_flat<< " "<< Vcov[0][j]<<" "<< Vcov[1][j]<<" "<< Vcov[2][j] << endl; 
	}
	
      }
    }

    //pT-dependence plots, setting bin contents and errors
    for(int k=0; k<5; k++){
      for (int j=0; j<n_eta-1; j++){

	double pave_value_corr = 0;
	int reject = 1;
	if(graph_filled[j]) {if(k==0) pave_value_corr = ptave_data[j]->GetMean();}
	else{
	  if(k==0) pave_value_corr = 0.0001;
	}
	
	if(k==1) pave_value_corr = 120;
	if(k==2) pave_value_corr = 60;
	if(k==3) pave_value_corr = 240;
	if(k==4) pave_value_corr = 480;
	
	//**************************************  Test reject high pT bins without data ***************************
	 
	double ptave_value = 0;
	if(graph_filled[j]) { ptave_value = ptave_data[j]->FindLastBinAbove(100.)*10;
	if(k==1 && ptave_value > 120) reject = 1;
	else{ if(k==1) reject = 1000;}
	if(k==2 && ptave_value > 60) reject = 1;
	else{ if(k==2)reject = 1000;}
	if(k==3 && ptave_value > 240) reject = 1;
	else{ if(k==3) reject = 1000;}
	if(k==4 && ptave_value > 480) reject = 1;
	else{if(k==4) reject = 1000;}
	}
	else{
	reject = 1000;
	}
	
	//************************************************************************************************************ 
  
                                
	//Residual_logpt_MPF->SetBinContent(j+1,hist_kfsr_fit_mpf->GetBinContent(j+1)*(f1[j]->GetParameter(0) +f1[j]->GetParameter(1)*TMath::Log(pave_value_corr)))
	//in above formula, no normalization has been applied!
	Residual_logpt_MPF->SetBinContent(j+1,reject*hist_kfsr_fit_mpf->GetBinContent(j+1)*(f1[j]->GetParameter(0)
										     +f1[j]->GetParameter(1)*TMath::Log(pave_value_corr))/loglin_norm);
	Residual_logpt_MPF->SetBinError(j+1, sqrt(pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(pave_value_corr),2)*pow(hist_kfsr_fit_mpf->GetBinError(j+1),2)   + pow(hist_kfsr_fit_mpf->GetBinContent(j+1),2)*(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(pave_value_corr),2)+2*Vcov[2][j]*TMath::Log(pave_value_corr)))/loglin_norm);
	ptave_logpt_MPF->SetBinError(j+1,sqrt(Vcov[0][j]
					      +Vcov[1][j]*pow(TMath::Log(pave_value_corr),2)+2*Vcov[2][j]*TMath::Log(pave_value_corr)));
	ptave_logpt_MPF->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(pave_value_corr));

                                                
	//Residual_const_MPF->SetBinContent(j+1,hist_kfsr_fit_mpf->GetBinContent(j+1)*f2[j]->GetParameter(0));
	//in above formula, no normalization has been applied!
	Residual_const_MPF->SetBinContent(j+1,reject*hist_kfsr_fit_mpf->GetBinContent(j+1)*f2[j]->GetParameter(0)/flat_norm);
	Residual_const_MPF->SetBinError(j+1, sqrt( pow(f2[j]->GetParameter(0)*hist_kfsr_fit_mpf->GetBinError(j+1),2)+pow( hist_kfsr_fit_mpf->GetBinContent(j+1)*f2[j]->GetParError(0),2)  ) /flat_norm );
	ptave_const_MPF->SetBinContent(j+1,f2[j]->GetParameter(0));
	ptave_const_MPF->SetBinError(j+1,f2[j]->GetParError(0));



	//******************************************************  USE kFSR Values/ NO FIT VALUES ************************************************************************************************
	Residual_logpt_MPF_val->SetBinContent(j+1,reject*hist_kfsr_mpf->GetBinContent(j+1)*(f1[j]->GetParameter(0)
										     +f1[j]->GetParameter(1)*TMath::Log(pave_value_corr))/loglin_norm_var);
	Residual_logpt_MPF_val->SetBinError(j+1, sqrt(pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(pave_value_corr),2)*pow(hist_kfsr_mpf->GetBinError(j+1),2)   + pow(hist_kfsr_mpf->GetBinContent(j+1),2)*(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(pave_value_corr),2)+2*Vcov[2][j]*TMath::Log(pave_value_corr)))/loglin_norm_var);
                                                
	//Residual_const_MPF->SetBinContent(j+1,hist_kfsr_mpf->GetBinContent(j+1)*f2[j]->GetParameter(0));
	//in above formula, no normalization has been applied!
	Residual_const_MPF_val->SetBinContent(j+1,reject*hist_kfsr_mpf->GetBinContent(j+1)*f2[j]->GetParameter(0)/flat_norm_var);
	Residual_const_MPF_val->SetBinError(j+1, sqrt( pow(f2[j]->GetParameter(0)*hist_kfsr_mpf->GetBinError(j+1),2)+pow( hist_kfsr_mpf->GetBinContent(j+1)*f2[j]->GetParError(0),2)  ) /flat_norm_var);

	//**********************************************************************************************************************************************************************************************
      }
    
      //File containing the results with kFSR correction applied
      TFile* outputfile;
      TString variation2;
      if(k==0){
	outputfile = new TFile(CorrectionObject::_outpath+"Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","RECREATE");
      }
      if(k>0){
	if(k==1) variation2 = "central";
	if(k==2) variation2 = "down";
	if(k==3) variation2 = "up";
	if(k==4) variation2 = "doubleup";
	outputfile = new TFile(CorrectionObject::_outpath+"Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_"+variation2+".root","RECREATE");
      }

      hist_kfsr_mpf->Write();
      hist_kfsr_fit_mpf->Write();                                                                                                                
      //kfsr_fit_mpf->Write();
      ptave_const_MPF->Write();                                                                                                                    
      ptave_logpt_MPF->Write();
      
      Residual_logpt_MPF->Write();
      Residual_const_MPF->Write();
      Residual_logpt_MPF_val->Write();
      Residual_const_MPF_val->Write();
      outputfile->Write();
      outputfile->Close();

      //File containing the pure ratios, no kFSR correction applied
      TFile* outputfile2;
      if(k==0){
	outputfile2 = new TFile(CorrectionObject::_outpath+"Histo_ptave_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","RECREATE");
      }
      if(k>0){
	if(k==1) variation2 = "central";
	if(k==2) variation2 = "down";
	if(k==3) variation2 = "up";
	if(k==4) variation2 = "doubleup";
	outputfile2 = new TFile(CorrectionObject::_outpath+"Histo_ptave_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_"+variation2+".root","RECREATE");
      }
    
      ptave_const_MPF->Write();
      ptave_logpt_MPF->Write();
      for (int j=0; j<n_eta-1; j++){
	f1[j]->Write();
	f2[j]->Write();
      }
      outputfile2->Write();
      outputfile2->Close();

      //delete stuff
      delete outputfile2;
      delete outputfile;
    }

    //delete stuff
    delete Residual_logpt_MPF;
    delete Residual_const_MPF;
    delete ptave_const_MPF;
    delete ptave_logpt_MPF;
    delete hist_kfsr_fit_mpf;
    delete hist_kfsr_mpf;
    // delete kfsr_fit_mpf;
    delete kfsr_mpf;
  } //if(mpfMethod) ends here


  /* ++++++++++++++++++++++++++ Calculate L2Residuals PT balance ++++++++++++++++++++++++++++++ */

  else{
    cout<<"HELLO EVERYBODY!"<<endl;

    TFile* kfsr_dijet;
    TH1D* hist_kfsr_fit_dijet;
    TH1D* hist_kfsr_dijet;

   
    if(CorrectionObject::_runnr != "BCDEFGH"){
        kfsr_dijet = new TFile(CorrectionObject::_input_path+"RunBCDEFGH/Histo_Res_DiJet_L1_"+CorrectionObject::_generator_tag+"_AK4PFchs.root","READ");
	hist_kfsr_fit_dijet = (TH1D*)kfsr_dijet->Get("hist_kfsr_fit_dijet");
	hist_kfsr_dijet = (TH1D*)kfsr_dijet->Get("kfsr_dijet");
    }
    
    else{
    kfsr_dijet = new TFile(CorrectionObject::_outpath+"Histo_KFSR_DiJet_"+CorrectionObject::_generator_tag+"_L1.root","READ");
    hist_kfsr_dijet = (TH1D*)kfsr_dijet->Get("kfsr_dijet");
    
    //kFSR fit function
    TF1 *kfsr_fit_dijet = new TF1("kfsr_fit_dijet","[0]+([1]*TMath::CosH(x))/(1+[2]*TMath::CosH(x))",0,5.); //Range: 0,5. by default
    
    //Be carefull with fit range!
    bool fit_fullrange = false;
    bool fit_285 = false;
  
    //Very fragile fit! tune the initial values for the fit
    ///Set Parameter for kFSR Fit (pT-Balance)
    if(CorrectionObject::_generator == "pythia"){
       if(CorrectionObject::_collection == "AK4CHS"){ 
	if(CorrectionObject::_runnr == "BCDEFGH"){
	  if(!CorrectionObject::_closuretest) kfsr_fit_dijet->SetParameters(-10,3600,300); //RES
	  else kfsr_fit_dijet->SetParameters(1,500,150); //CLOSURETEST
	  //fit_fullrange = true;
	  fit_285 = true;
	}
       }
    }
    
    else throw runtime_error("PTextrapolation: Invalid generator specified.");

    //Finally perform the fit 
    kfsr_fit_dijet->SetLineColor(kBlue+1);
    std::cout<<"------ !!! kFSR pT-balance fit !!! ------"<<std::endl;
    if(fit_fullrange) hist_kfsr_dijet->Fit("kfsr_fit_dijet","SR","",0,5.);
    else if(!fit_285) hist_kfsr_dijet->Fit("kfsr_fit_dijet","SR","",0,2.65);
    else              hist_kfsr_dijet->Fit("kfsr_fit_dijet","SR","",0,2.85);
 
    //Create a histogram to hold the confidence intervals  
    hist_kfsr_fit_dijet = (TH1D*)hist_kfsr_dijet->Clone();
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hist_kfsr_fit_dijet);
    //Now the "hist_kfsr_fit" histogram has the fitted function values as the
    //bin contents and the confidence intervals as bin errors

    hist_kfsr_fit_dijet->SetStats(kFALSE);
    hist_kfsr_fit_dijet->SetFillColor(kBlue-10);     
    hist_kfsr_fit_dijet->SetName("hist_kfsr_fit_dijet");
    hist_kfsr_fit_dijet->SetTitle("kfsr fit for dijet");
    }
  

    double flat[n_eta-1];
    double loglin[n_eta-1];

    double flat_var[n_eta-1];
    double loglin_var[n_eta-1];

    double flat_norm;
    double loglin_norm;

    double flat_norm_var;
    double loglin_norm_var;


    for (int j=0; j<n_eta-1; j++){
      flat[j] = hist_kfsr_fit_dijet->GetBinContent(j+1)*f2[j]->GetParameter(0);
      loglin[j] = hist_kfsr_fit_dijet->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()) );
    }

    flat_norm = 0, loglin_norm = 0;
    for (int j=0; j<n_etabarr; j++){
      flat_norm += flat[j];
      loglin_norm += loglin[j];
    }
    flat_norm = flat_norm/n_etabarr;
    loglin_norm = loglin_norm/n_etabarr;
    for (int j=0; j<n_eta-1; j++){
      flat[j] = flat[j]/flat_norm;
      loglin[j] = loglin[j]/loglin_norm;
    }

    for (int j=0; j<n_eta-1; j++){
      flat_var[j] = hist_kfsr_dijet->GetBinContent(j+1)*f2[j]->GetParameter(0);
      loglin_var[j] = hist_kfsr_dijet->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()) );
    }

    flat_norm_var = 0, loglin_norm_var = 0;
    for (int j=0; j<n_etabarr; j++){
      flat_norm_var += flat_var[j];
      loglin_norm_var += loglin_var[j];
    }
    flat_norm_var = flat_norm_var/n_etabarr;
    loglin_norm_var = loglin_norm_var/n_etabarr;
    for (int j=0; j<n_eta-1; j++){
      flat_var[j] = flat_var[j]/flat_norm_var;
      loglin_var[j] = loglin_var[j]/loglin_norm_var;
    }
    //   }
 
    //Histograms holding the output
    TH1D* Residual_logpt_DiJet = new TH1D("res_logpt_dijet","res_logpt_dijet", n_eta-1,eta_bins);
    TH1D* Residual_const_DiJet = new TH1D("res_const_dijet","res_const_dijet", n_eta-1,eta_bins);
    TH1D* ptave_const_DiJet    = new TH1D("ptave_const_dijet","ptave_const_dijet", n_eta-1,eta_bins);
    TH1D* ptave_logpt_DiJet    = new TH1D("ptave_logpt_dijet","ptave_logpt_dijet", n_eta-1,eta_bins);

    TH1D* Residual_logpt_DiJet_val = new TH1D("res_logpt_dijet_val","res_logpt_dijet_val", n_eta-1,eta_bins);
    TH1D* Residual_const_DiJet_val = new TH1D("res_const_dijet_val","res_const_dijet_val", n_eta-1,eta_bins);

    string kFSR_method[2] = {"kFSR_Fit", "kFSR_Hist"};    

    for(int i=0; i<2;i++){ 

    ofstream output, output_loglin, uncerts, uncerts_loglin, output_hybrid, uncerts_hybrid, output_hybrid_barrel, uncerts_hybrid_barrel;
    output.open(CorrectionObject::_outpath+"output/"+kFSR_method[i]+"/Summer16_07Aug2017_pT_FLAT_L2Residual_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt");
    output_loglin.open(CorrectionObject::_outpath+"output/"+kFSR_method[i]+"/Summer16_07Aug2017_pT_LOGLIN_L2Residual_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt");
    uncerts.open(CorrectionObject::_outpath+"output/"+kFSR_method[i]+"/Summer16_07Aug2017_pT_FLAT_L2Residual_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt.STAT");
    uncerts_loglin.open(CorrectionObject::_outpath+"output/"+kFSR_method[i]+"/Summer16_07Aug2017_pT_LOGLIN_L2Residual_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt.STAT");

    output_hybrid.open(CorrectionObject::_outpath+"output/"+kFSR_method[i]+"/Summer16_07Aug2017_pT_Hybrid_L2Residual_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt");
    uncerts_hybrid.open(CorrectionObject::_outpath+"output/"+kFSR_method[i]+"/Summer16_07Aug2017_pT_Hybrid_L2Residual_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt.STAT");

    output_hybrid_barrel.open(CorrectionObject::_outpath+"output/"+kFSR_method[i]+"/Summer16_07Aug2017_pT_Hybrid_Barrel_L2Residual_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt");
    uncerts_hybrid_barrel.open(CorrectionObject::_outpath+"output/"+kFSR_method[i]+"/Summer16_07Aug2017_pT_Hybrid_Barrel_L2Residual_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt.STAT");

  
    output  << "{ 1 JetEta 1 JetPt[2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*0.021*(-1.+1./(1.+exp(-(TMath::Log(x)-5.030)/0.395))))) Correction L2Relative}" << endl;

    output_loglin  << "{ 1 JetEta 1 JetPt[2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*0.021*(-1.+1./(1.+exp(-(TMath::Log(x)-5.030)/0.395))))) Correction L2Relative}" << endl;
    uncerts << "{ 1 JetEta 1 JetPt [0] kFSR_err f0_err}" << endl;

    uncerts_loglin << "{ 1 JetEta 1 JetPt sqrt(fabs([0]*[0]+[1]+[2]*log(x)*log(x)+2*[3]*log(x))) Correction L2Relative}" << endl;

    output_hybrid  << "{ 1 JetEta 1 JetPt[2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*0.021*(-1.+1./(1.+exp(-(TMath::Log(x)-5.030)/0.395))))) Correction L2Relative}" << endl;

    uncerts_hybrid << "{ 1 JetEta 1 JetPt Hybrid Method: Use const fit for 2.65 < |eta| < 3.1, else Log-lin}" << endl;

    output_hybrid_barrel  << "{ 1 JetEta 1 JetPt[2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*0.021*(-1.+1./(1.+exp(-(TMath::Log(x)-5.030)/0.395))))) Correction L2Relative}" << endl;

    uncerts_hybrid_barrel << "{ 1 JetEta 1 JetPt Hybrid Method: Use const fit for 2.65 < |eta| < 3.1, else Log-lin}" << endl;

 
    double norm_flat=0;
    double norm_loglin=0;
    
    TH1D* kfsr_FitVal_dijet;

    if(kFSR_method[i]=="kFSR_Fit"){
      cout<<"!kFSR: Fit Values used for txt-files and plots!"<<endl;
      kfsr_FitVal_dijet = (TH1D*) hist_kfsr_fit_dijet->Clone();
      norm_flat = flat_norm;
      norm_loglin = loglin_norm;
    }
    else{
      cout<<"!kFSR: Hist Values used for txt-files and plots!"<<endl;
      kfsr_FitVal_dijet = (TH1D*) hist_kfsr_dijet->Clone();
      norm_flat = flat_norm_var;
      norm_loglin = loglin_norm_var;
    }



    for (int j=n_eta-1; j>0; --j){
      output << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   30   " << ptave_data[j-1]->FindLastBinAbove(100.)*10 << "   " << 1/norm_flat << " " << kfsr_FitVal_dijet->GetBinContent(j) << "   " << f2[j-1]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;
      
      output_loglin << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   30   " << ptave_data[j-1]->FindLastBinAbove(100.)*10 << "   " << 1/norm_loglin << " " << kfsr_FitVal_dijet->GetBinContent(j) << "   " << f1[j-1]->GetParameter(0) << " " << f1[j-1]->GetParameter(1) << "   1 0.0000 0.0" << endl;

      uncerts << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    30    " 
	      << ptave_data[j-1]->FindLastBinAbove(100.)*10 << "  " << kfsr_FitVal_dijet->GetBinError(j-1) / norm_flat 
	      <<" "<< f2[j-1]->GetParError(0) << endl;

      uncerts_loglin << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" 
		     << eta_range[j-1] << "  6    30    " << ptave_data[j-1]->FindLastBinAbove(100.)*10 
		     << " " << kfsr_FitVal_dijet->GetBinError(j-1) / norm_loglin  
		     << " " <<Vcov[0][j-1] <<" "<<Vcov[1][j-1]<<" "<<Vcov[2][j-1]<<endl; 
      
      
      if(j>12 && j<16){
	output_hybrid << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   30   " << ptave_data[j-1]->FindLastBinAbove(100.)*10 << "   " << 1/norm_flat << " " << kfsr_FitVal_dijet->GetBinContent(j) << "   " << f2[j-1]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;

	uncerts_hybrid << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    30    " 
	      << ptave_data[j-1]->FindLastBinAbove(100.)*10 << "  " << kfsr_FitVal_dijet->GetBinError(j-1) / norm_flat 
	      <<" "<< f2[j-1]->GetParError(0) << endl;
      }
      else{
	output_hybrid <<  fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   30   " << ptave_data[j-1]->FindLastBinAbove(100.)*10 << "   " << 1/norm_loglin << " " << kfsr_FitVal_dijet->GetBinContent(j) << "   " << f1[j-1]->GetParameter(0) << " " << f1[j-1]->GetParameter(1) << "   1 0.0000 0.0" << endl;

        uncerts_hybrid << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" 
		     << eta_range[j-1] << "  6    30    " << ptave_data[j-1]->FindLastBinAbove(100.)*10 
		     << " " << kfsr_FitVal_dijet->GetBinError(j-1) / norm_loglin  
		     << " " <<Vcov[0][j-1] <<" "<<Vcov[1][j-1]<<" "<<Vcov[2][j-1]<<endl; 
      }
      

      if(j<6){
	output_hybrid_barrel << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   30   " << ptave_data[j-1]->FindLastBinAbove(100.)*10 << "   " << 1/norm_flat << " " << kfsr_FitVal_dijet->GetBinContent(j) << "   " << f2[j-1]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;

	uncerts_hybrid_barrel << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    30    " 
	      << ptave_data[j-1]->FindLastBinAbove(100.)*10 << "  " << kfsr_FitVal_dijet->GetBinError(j-1) / norm_flat 
	      <<" "<< f2[j-1]->GetParError(0) << endl;
      }
      else{
	output_hybrid_barrel <<  fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   30   " << ptave_data[j-1]->FindLastBinAbove(100.)*10 << "   " << 1/norm_flat << " " << kfsr_FitVal_dijet->GetBinContent(j) << "   " << f1[j-1]->GetParameter(0) << " " << f1[j-1]->GetParameter(1) << "   1 0.0000 0.0" << endl;

        uncerts_hybrid_barrel << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" 
		     << eta_range[j-1] << "  6    30    " << ptave_data[j-1]->FindLastBinAbove(100.)*10 
		     << " " << kfsr_FitVal_dijet->GetBinError(j-1) / norm_flat  
		     << " " <<Vcov[0][j-1] <<" "<<Vcov[1][j-1]<<" "<<Vcov[2][j-1]<<endl; 
      }
    }


    for (int j=0; j<n_eta-1; j++){
 
      output << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   55   " 
	     << ptave_data[j]->FindLastBinAbove(100.)*10 << "   " << 1/norm_flat << " " << kfsr_FitVal_dijet->GetBinContent(j+1) 
	     << "   " << f2[j]->GetParameter(0) << " 0   1 0.0000 0.0"  << endl;

      output_loglin << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   30   " 
		    << ptave_data[j]->FindLastBinAbove(100.)*10 << "   " << 1/norm_loglin << " " 
		    << kfsr_FitVal_dijet->GetBinContent(j+1) << "   " << f1[j]->GetParameter(0) 
		    << " " << f1[j]->GetParameter(1) << "   1 0.0000 0.0" << endl;

      uncerts << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " 
	      << ptave_data[j]->FindLastBinAbove(100.)*10 
	      <<" "<< kfsr_FitVal_dijet->GetBinError(j)/norm_flat << " " << f2[j]->GetParameter(0) <<endl;

      uncerts_loglin << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  6    30    " 
		     << ptave_data[j]->FindLastBinAbove(100.)*10 << " " 
		     << kfsr_FitVal_dijet->GetBinError(j)/norm_loglin << " "<< Vcov[0][j]<<" "<< Vcov[1][j]<<" "<< Vcov[2][j] << endl; 

      
      if(j>11 && j<15){
	output_hybrid << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   30   " 
	     << ptave_data[j]->FindLastBinAbove(100.)*10 << "   " << 1/norm_flat << " " << kfsr_FitVal_dijet->GetBinContent(j+1) 
	     << "   " << f2[j]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;

	uncerts_hybrid << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    30    " 
	      << ptave_data[j]->FindLastBinAbove(100.)*10 
	      <<" "<< kfsr_FitVal_dijet->GetBinError(j)/norm_flat << " " << f2[j]->GetParameter(0) <<endl;
      }
      else{
	output_hybrid << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   30   " 
		    << ptave_data[j]->FindLastBinAbove(100.)*10 << "   " << 1/norm_loglin << " " 
		    << kfsr_FitVal_dijet->GetBinContent(j+1) << "   " << f1[j]->GetParameter(0) 
		    << " " << f1[j]->GetParameter(1) << "   1 0.0000 0.0" << endl; 

        uncerts_hybrid << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  6    30    " 
		     << ptave_data[j]->FindLastBinAbove(100.)*10 << " " 
		     << kfsr_FitVal_dijet->GetBinError(j)/norm_loglin << " "<< Vcov[0][j]<<" "<< Vcov[1][j]<<" "<< Vcov[2][j] << endl; 
      }
      

      if(j<5){
	output_hybrid_barrel << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   30   " 
	     << ptave_data[j]->FindLastBinAbove(100.)*10 << "   " << 1/norm_flat << " " << kfsr_FitVal_dijet->GetBinContent(j+1) 
	     << "   " << f2[j]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;

	uncerts_hybrid_barrel << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    30    " 
	      << ptave_data[j]->FindLastBinAbove(100.)*10 
	      <<" "<< kfsr_FitVal_dijet->GetBinError(j)/norm_flat << " " << f2[j]->GetParameter(0) <<endl;
      }
      else{
	output_hybrid_barrel << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   30   " 
		    << ptave_data[j]->FindLastBinAbove(100.)*10 << "   " << 1/norm_flat << " " 
		    << kfsr_FitVal_dijet->GetBinContent(j+1) << "   " << f1[j]->GetParameter(0) 
		    << " " << f1[j]->GetParameter(1) << "   1 0.0000 0.0" << endl; 

        uncerts_hybrid_barrel << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  6    30    " 
		     << ptave_data[j]->FindLastBinAbove(100.)*10 << " " 
		     << kfsr_FitVal_dijet->GetBinError(j)/norm_flat << " "<< Vcov[0][j]<<" "<< Vcov[1][j]<<" "<< Vcov[2][j] << endl; 
      }
    }
    }


    for(int k=0; k<5; k++){  
      for (int j=0; j<n_eta-1; j++){

	double pave_corr;
	int reject = 1;
	if(graph_filled[j]) {if(k==0) pave_corr = ptave_data[j]->GetMean();}
	else{
	  if(k==0) pave_corr = 0.0001;
	}
	
	if(k==1) pave_corr = 120;
	if(k==2) pave_corr = 60;
	if(k==3) pave_corr = 240;
	if(k==4) pave_corr = 480;
	
	
	//**************************************  reject high pT bins without data ***************************
	double ptave_value = 0;
	if(graph_filled[j]) { ptave_value = ptave_data[j]->FindLastBinAbove(100.)*10;
	if(k==1 && ptave_value > 120) reject = 1;
	else{ if(k==1) reject = 1000;}
	if(k==2 && ptave_value > 60) reject = 1;
	else{ if(k==2)reject = 1000;}
	if(k==3 && ptave_value > 240) reject = 1;
	else{ if(k==3) reject = 1000;}
	if(k==4 && ptave_value > 480) reject = 1;
	else{if(k==4) reject = 1000;}
	}
	else{
	reject = 1000;
	}
	
	//************************************************************************************************************ 

	Residual_logpt_DiJet->SetBinContent(j+1,reject*hist_kfsr_fit_dijet->GetBinContent(j+1)*(f1[j]->GetParameter(0)
											 +f1[j]->GetParameter(1)*TMath::Log(pave_corr))/loglin_norm);
	Residual_logpt_DiJet->SetBinError(j+1, sqrt(pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(pave_corr),2)*pow(hist_kfsr_fit_dijet->GetBinError(j+1),2)+ pow(hist_kfsr_fit_dijet->GetBinContent(j+1),2)*(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(pave_corr),2)+2*Vcov[2][j]*TMath::Log(pave_corr)))/loglin_norm);

	ptave_logpt_DiJet->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(pave_corr));
	ptave_logpt_DiJet->SetBinError(j+1,sqrt((Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(ptave_data[j]->GetMean()),2)+2*Vcov[2][j]*TMath::Log(pave_corr))));

	Residual_const_DiJet->SetBinContent(j+1,reject*hist_kfsr_fit_dijet->GetBinContent(j+1)*f2[j]->GetParameter(0)/flat_norm);
	Residual_const_DiJet->SetBinError(j+1, sqrt( pow(f2[j]->GetParameter(0)*hist_kfsr_fit_dijet->GetBinError(j+1),2)+pow(hist_kfsr_fit_dijet->GetBinContent(j+1)*f2[j]->GetParError(0),2)  )/flat_norm  );
	ptave_const_DiJet->SetBinContent(j+1,f2[j]->GetParameter(0));
	ptave_const_DiJet->SetBinError(j+1,f2[j]->GetParError(0));


	//*********************************************************************** TEST: kFSR VALUES/ NO FIT USED
	Residual_logpt_DiJet_val->SetBinContent(j+1,reject*hist_kfsr_dijet->GetBinContent(j+1)*(f1[j]->GetParameter(0)
											 +f1[j]->GetParameter(1)*TMath::Log(pave_corr))/loglin_norm_var);
	Residual_logpt_DiJet_val->SetBinError(j+1, sqrt(pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(pave_corr),2)*pow(hist_kfsr_dijet->GetBinError(j+1),2)+ pow(hist_kfsr_dijet->GetBinContent(j+1),2)*(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(pave_corr),2)+2*Vcov[2][j]*TMath::Log(pave_corr)))/loglin_norm_var);

	Residual_const_DiJet_val->SetBinContent(j+1,reject*hist_kfsr_dijet->GetBinContent(j+1)*f2[j]->GetParameter(0)/flat_norm_var);
	Residual_const_DiJet_val->SetBinError(j+1, sqrt( pow(f2[j]->GetParameter(0)*hist_kfsr_dijet->GetBinError(j+1),2)+pow(hist_kfsr_dijet->GetBinContent(j+1)*f2[j]->GetParError(0),2)  )/flat_norm_var);
	//*************************************************************************************************************************************************
      }


      //File containing the corrections, with kFSR correction applied
      TFile* outputfile;
      TString variation2;
      if(k==0){
	outputfile = new TFile(CorrectionObject::_outpath+"Histo_Res_DiJet_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","RECREATE");
      }
      if(k>0){
	if(k==1) variation2 = "central";
	if(k==2) variation2 = "down";
	if(k==3) variation2 = "up";
	if(k==4) variation2 = "doubleup";
	outputfile = new TFile(CorrectionObject::_outpath+"Histo_Res_DiJet_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_"+variation2+".root","RECREATE");
      }
    
      hist_kfsr_dijet->Write();
      hist_kfsr_fit_dijet->Write();
      //kfsr_fit_dijet->Write();
      ptave_const_DiJet->Write();
      ptave_logpt_DiJet->Write();
                  
      Residual_logpt_DiJet->Write();
      Residual_const_DiJet->Write();
      Residual_logpt_DiJet_val->Write();
      Residual_const_DiJet_val->Write();
      outputfile->Write();
      outputfile->Close();


      //Contains the pT-Extrapolations
      TFile* outputfile2;
      if(k==0){
	outputfile2 = new TFile(CorrectionObject::_outpath+"Histo_ptave_DiJet_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","RECREATE");
      }
      if(k>0){
	if(k==1) variation2 = "central";
	if(k==2) variation2 = "down";
	if(k==3) variation2 = "up";
	if(k==4) variation2 = "doubleup";
	outputfile2 = new TFile(CorrectionObject::_outpath+"Histo_ptave_DiJet_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_"+variation2+".root","RECREATE");
      }
    

      ptave_const_DiJet->Write();
      ptave_logpt_DiJet->Write();
      for (int j=0; j<n_eta-1; j++){
	f1[j]->Write();
	f2[j]->Write();
      }             
      outputfile2->Write();
      outputfile2->Close();

      //delete stuff
      delete outputfile2;
      delete outputfile;
    }



    //delete stuff
    delete ptave_logpt_DiJet;
    delete ptave_const_DiJet;
    delete Residual_const_DiJet;
    delete Residual_logpt_DiJet;

    // delete hist_kfsr_fit_dijet;
    // delete hist_kfsr_dijet;
    // delete kfsr_fit_dijet;
    // delete kfsr_dijet;
  } //end of (if(mpfMethod...else))




//delete everything
  delete c_kfsr_fit;

  for (int j=0; j<n_eta-1; j++) {
    delete asd[j];
    delete f1[j];
    delete f2[j];
  }

  delete line;

  for(int j=0; j<n_eta-1; j++){
    delete graph1_mpf[j];
  }
 
  for(int j=0; j<n_eta-1; j++){
    delete pr_data_asymmetry[j];
    delete pr_data_B[j];
    delete pr_mc_asymmetry[j];
    delete pr_mc_B[j];
    delete hdata_asymmetry[j];
    delete hdata_B[j];
    delete hmc_asymmetry[j];
    delete hmc_B[j];
  }
     
  for(int i=0; i<n_eta-1; i++) delete ptave_data[i];
  delete m_gStyle;

}