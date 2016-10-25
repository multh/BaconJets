#include "header.h"
//#include "tdrstyle_mod14.C"
//#include "tdrstyle_mod15.C"

void L2ResOutput(bool divide_by_lumi, bool use_sys_hists , int current_lumibin, TString path, TString path_up, TString path_down, TString txttag, TString jettag, TString tag, TString lumitag, double al_cut=0.2){

  if(divide_by_lumi){
    cout << "Plots are made in bins of instantaneous lumi, this time for " << lumi_bins[current_lumibin-1] << " < lumi < " << lumi_bins[current_lumibin] << endl;
  }

  if(divide_by_lumi) {use_sys_hists = false; cout << "Warning, for lumi-binned plots no JER uncertainties are drawn -- too much effort." << endl;}

  TFile* f_Res_mpf, *f_Res_mpf_up, *f_Res_mpf_down;
  if(divide_by_lumi) f_Res_mpf = new TFile(path+"Histo_Res_MPF_L1_"+txttag+"_"+jettag+tag+"_Lumi_"+lumi_range[current_lumibin-1]+"_"+lumi_range[current_lumibin]+".root","READ");   
  else f_Res_mpf = new TFile(path+"Histo_Res_MPF_L1_"+txttag+"_"+jettag+tag+".root","READ");   
  if(use_sys_hists){
    f_Res_mpf_up = new TFile(path_up+"Histo_Res_MPF_L1_"+txttag+"_"+jettag+tag+".root","READ");
    f_Res_mpf_down = new TFile(path_down+"Histo_Res_MPF_L1_"+txttag+"_"+jettag+tag+".root","READ");
  }
  

  TFile* f_Res_dijet, *f_Res_dijet_up, *f_Res_dijet_down;
  if(divide_by_lumi) f_Res_dijet = new TFile(path+"Histo_Res_DiJet_L1_"+txttag+"_"+jettag+tag+"_Lumi_"+lumi_range[current_lumibin-1]+"_"+lumi_range[current_lumibin]+".root","READ");
  else f_Res_dijet = new TFile(path+"Histo_Res_DiJet_L1_"+txttag+"_"+jettag+tag+".root","READ");  
  if(use_sys_hists){
    f_Res_dijet_up = new TFile(path_up+"Histo_Res_DiJet_L1_"+txttag+"_"+jettag+tag+".root","READ");
    f_Res_dijet_down = new TFile(path_down+"Histo_Res_DiJet_L1_"+txttag+"_"+jettag+tag+".root","READ");
  }


  TString JetDescrib;                                                                                                                            
  if (jettag=="AK4PFchs") JetDescrib = "Anti-k_{t} R = 0.4, PF+CHS";
  if (jettag=="AK4PFpuppi") JetDescrib = "Anti-k_{t} R = 0.4, PF+PUPPI";
  if (jettag=="AK8PFchs") JetDescrib = "Anti-k_{t} R = 0.8, PF+CHS";
  if (jettag=="AK8PFpuppi") JetDescrib = "Anti-k_{t} R = 0.8, PF+PUPPI"; 

  //plot results for "nominal" variation
 

  // get the (R_{MC}/R_{DATA}) hists for MPF and pt balance
  TH1D* pt_depend_const_mpf = (TH1D*)f_Res_mpf->Get("ptave_const_mpf");
  TH1D* pt_depend_logpt_mpf = (TH1D*)f_Res_mpf->Get("ptave_logpt_mpf");
  TH1D* pt_depend_const_mpf_up, *pt_depend_const_mpf_down, *pt_depend_logpt_mpf_up, *pt_depend_logpt_mpf_down;
  if(use_sys_hists){
    pt_depend_const_mpf_up = (TH1D*)f_Res_mpf_up->Get("ptave_const_mpf");
    pt_depend_const_mpf_down = (TH1D*)f_Res_mpf_down->Get("ptave_const_mpf");
    pt_depend_logpt_mpf_up = (TH1D*)f_Res_mpf_up->Get("ptave_logpt_mpf");
    pt_depend_logpt_mpf_down = (TH1D*)f_Res_mpf_down->Get("ptave_logpt_mpf");
  }
 
  TH1D* pt_depend_const_dijet = (TH1D*)f_Res_dijet->Get("ptave_const_dijet");
  TH1D* pt_depend_logpt_dijet = (TH1D*)f_Res_dijet->Get("ptave_logpt_dijet");
  TH1D* pt_depend_const_dijet_up, *pt_depend_const_dijet_down, *pt_depend_logpt_dijet_up, *pt_depend_logpt_dijet_down;
  if(use_sys_hists){
    pt_depend_const_dijet_up = (TH1D*)f_Res_dijet_up->Get("ptave_const_dijet");
    pt_depend_const_dijet_down = (TH1D*)f_Res_dijet_down->Get("ptave_const_dijet");
    pt_depend_logpt_dijet_up = (TH1D*)f_Res_dijet_up->Get("ptave_logpt_dijet");
    pt_depend_logpt_dijet_down = (TH1D*)f_Res_dijet_down->Get("ptave_logpt_dijet");
 }

  // get the kFSR hists for MPF and pt balance
  TH1D* kfsr_mpf = (TH1D*)f_Res_mpf->Get("kfsr_mpf");
  TH1D* kfsr_mpf_fit = (TH1D*)f_Res_mpf->Get("hist_kfsr_fit_mpf");
  TH1D* kfsr_dijet = (TH1D*)f_Res_dijet->Get("kfsr_dijet");
  TH1D* kfsr_dijet_fit = (TH1D*)f_Res_dijet->Get("hist_kfsr_fit_dijet");
  TH1D* kfsr_mpf_fit_up, *kfsr_mpf_fit_down, *kfsr_mpf_up, *kfsr_mpf_down, *kfsr_dijet_fit_up, *kfsr_dijet_fit_down, *kfsr_dijet_up, *kfsr_dijet_down;
  if(use_sys_hists){
    kfsr_mpf_fit_up = (TH1D*)f_Res_mpf_up->Get("hist_kfsr_fit_mpf");
    kfsr_mpf_fit_down = (TH1D*)f_Res_mpf_down->Get("hist_kfsr_fit_mpf");
    kfsr_mpf_up = (TH1D*)f_Res_mpf_up->Get("hist_kfsr_mpf");
    kfsr_mpf_down = (TH1D*)f_Res_mpf_down->Get("hist_kfsr_mpf");
    kfsr_dijet_fit_up = (TH1D*)f_Res_dijet_up->Get("hist_kfsr_fit_dijet");
    kfsr_dijet_fit_down = (TH1D*)f_Res_dijet_down->Get("hist_kfsr_fit_dijet");
    kfsr_dijet_up = (TH1D*)f_Res_dijet_up->Get("hist_kfsr_dijet");
    kfsr_dijet_down = (TH1D*)f_Res_dijet_down->Get("hist_kfsr_dijet");
  }

  //get L2Res hists for MPF and pt balance
  TH1D* res_const_mpf_kfsrfit = (TH1D*)f_Res_mpf->Get("res_const_mpf");
  TH1D* res_const_dijet_kfsrfit = (TH1D*)f_Res_dijet->Get("res_const_dijet");
  TH1D* res_logpt_mpf_kfsrfit = (TH1D*)f_Res_mpf->Get("res_logpt_mpf");
  TH1D* res_logpt_dijet_kfsrfit = (TH1D*)f_Res_dijet->Get("res_logpt_dijet");
  TH1D* res_const_mpf_kfsrfit_up, *res_const_mpf_kfsrfit_down, *res_logpt_mpf_kfsrfit_up, *res_logpt_mpf_kfsrfit_down, *res_const_dijet_kfsrfit_up, *res_const_dijet_kfsrfit_down, *res_logpt_dijet_kfsrfit_up, *res_logpt_dijet_kfsrfit_down;
  if(use_sys_hists){
    res_const_mpf_kfsrfit_up = (TH1D*)f_Res_mpf_up->Get("res_const_mpf");
    res_const_mpf_kfsrfit_down = (TH1D*)f_Res_mpf_down->Get("res_const_mpf");
    res_logpt_mpf_kfsrfit_up = (TH1D*)f_Res_mpf_up->Get("res_logpt_mpf");
    res_logpt_mpf_kfsrfit_down = (TH1D*)f_Res_mpf_down->Get("res_logpt_mpf");
    res_const_dijet_kfsrfit_up = (TH1D*)f_Res_dijet_up->Get("res_const_dijet");
    res_const_dijet_kfsrfit_down = (TH1D*)f_Res_dijet_down->Get("res_const_dijet");
    res_logpt_dijet_kfsrfit_up = (TH1D*)f_Res_dijet_up->Get("res_logpt_dijet");
    res_logpt_dijet_kfsrfit_down = (TH1D*)f_Res_dijet_down->Get("res_logpt_dijet");
  }
 
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


//   // Draw results
  TH1D *h = new TH1D("h",";|#eta|;Relative correction",41,0,5.191);
  h->SetMaximum(1.2);
  h->SetMinimum(0.8);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  //  lumi_13TeV = "Run2015  2.1 fb^{-1}";
  //lumi_13TeV = "Run2016  12.9 fb^{-1}";
  lumi_13TeV = lumitag;
  //lumi_13TeV = "589.3 pb^{-1}";
  bool kSquare = true;

  TGraphAsymmErrors* gr_pt_depend_const_mpf = new TGraphAsymmErrors(pt_depend_const_mpf);
  TGraphAsymmErrors* gr_pt_depend_logpt_mpf = new TGraphAsymmErrors(pt_depend_logpt_mpf);
  TGraphAsymmErrors* gr_pt_depend_const_dijet = new TGraphAsymmErrors(pt_depend_const_dijet);
  TGraphAsymmErrors* gr_pt_depend_logpt_dijet = new TGraphAsymmErrors(pt_depend_logpt_dijet);
  if(use_sys_hists){
    const int n_bins_x = pt_depend_const_mpf->GetNbinsX();
    double mpf_const_up[n_bins_x], mpf_const_down[n_bins_x], mpf_logpt_up[n_bins_x], mpf_logpt_down[n_bins_x], dijet_const_up[n_bins_x], dijet_const_down[n_bins_x], dijet_logpt_up[n_bins_x], dijet_logpt_down[n_bins_x];
    for(int i=1; i<pt_depend_const_mpf->GetNbinsX()+1; i++){
      //calculate systematic errors
      double mpf_const_sys_up=0, mpf_const_sys_down=0, mpf_logpt_sys_up=0, mpf_logpt_sys_down=0, dijet_const_sys_up=0, dijet_const_sys_down=0, dijet_logpt_sys_up=0, dijet_logpt_sys_down=0;

      //start with mpf_const
      if(pt_depend_const_mpf->GetBinContent(i) > pt_depend_const_mpf_up->GetBinContent(i) && pt_depend_const_mpf->GetBinContent(i) < pt_depend_const_mpf_down->GetBinContent(i)) {
	mpf_const_sys_up = pt_depend_const_mpf_down->GetBinContent(i) - pt_depend_const_mpf->GetBinContent(i);
	mpf_const_sys_down = pt_depend_const_mpf->GetBinContent(i) - pt_depend_const_mpf_up->GetBinContent(i); 
      } 

      else if(pt_depend_const_mpf->GetBinContent(i) < pt_depend_const_mpf_up->GetBinContent(i) && pt_depend_const_mpf->GetBinContent(i) > pt_depend_const_mpf_down->GetBinContent(i)){
	mpf_const_sys_up = pt_depend_const_mpf_up->GetBinContent(i) - pt_depend_const_mpf->GetBinContent(i);
	mpf_const_sys_down = pt_depend_const_mpf->GetBinContent(i) - pt_depend_const_mpf_down->GetBinContent(i); 
      }

      else if(pt_depend_const_mpf->GetBinContent(i) > pt_depend_const_mpf_up->GetBinContent(i) && pt_depend_const_mpf->GetBinContent(i) > pt_depend_const_mpf_down->GetBinContent(i)){
	if(pt_depend_const_mpf_up->GetBinContent(i) < pt_depend_const_mpf_down->GetBinContent(i)) mpf_const_sys_down = pt_depend_const_mpf->GetBinContent(i) - pt_depend_const_mpf_up->GetBinContent(i);
	else mpf_const_sys_down = pt_depend_const_mpf->GetBinContent(i) - pt_depend_const_mpf_down->GetBinContent(i); 
      }

      else if(pt_depend_const_mpf->GetBinContent(i) < pt_depend_const_mpf_up->GetBinContent(i) && pt_depend_const_mpf->GetBinContent(i) < pt_depend_const_mpf_down->GetBinContent(i)){
	if(pt_depend_const_mpf_up->GetBinContent(i) < pt_depend_const_mpf_down->GetBinContent(i)) mpf_const_sys_up = pt_depend_const_mpf_down->GetBinContent(i) - pt_depend_const_mpf->GetBinContent(i);
	else mpf_const_sys_up = pt_depend_const_mpf_up->GetBinContent(i) - pt_depend_const_mpf->GetBinContent(i); 
      }

      else {throw runtime_error("what did i do?! - mpf const");}
 
      //again for mpf_logpt
      if(pt_depend_logpt_mpf->GetBinContent(i) > pt_depend_logpt_mpf_up->GetBinContent(i) && pt_depend_logpt_mpf->GetBinContent(i) < pt_depend_logpt_mpf_down->GetBinContent(i)) {
	mpf_logpt_sys_up = pt_depend_logpt_mpf_down->GetBinContent(i) - pt_depend_logpt_mpf->GetBinContent(i);
	mpf_logpt_sys_down = pt_depend_logpt_mpf->GetBinContent(i) - pt_depend_logpt_mpf_up->GetBinContent(i);
      } 
      else if(pt_depend_logpt_mpf->GetBinContent(i) < pt_depend_logpt_mpf_up->GetBinContent(i) && pt_depend_logpt_mpf->GetBinContent(i) > pt_depend_logpt_mpf_down->GetBinContent(i)){
	mpf_logpt_sys_up = pt_depend_logpt_mpf_up->GetBinContent(i) - pt_depend_logpt_mpf->GetBinContent(i);
	mpf_logpt_sys_down = pt_depend_logpt_mpf->GetBinContent(i) - pt_depend_logpt_mpf_down->GetBinContent(i);
      }
      else if(pt_depend_logpt_mpf->GetBinContent(i) > pt_depend_logpt_mpf_up->GetBinContent(i) && pt_depend_logpt_mpf->GetBinContent(i) > pt_depend_logpt_mpf_down->GetBinContent(i)){
	if(pt_depend_logpt_mpf_up->GetBinContent(i) < pt_depend_logpt_mpf_down->GetBinContent(i)) mpf_logpt_sys_down = pt_depend_logpt_mpf->GetBinContent(i) - pt_depend_logpt_mpf_up->GetBinContent(i);
	else mpf_logpt_sys_down = pt_depend_logpt_mpf->GetBinContent(i) - pt_depend_logpt_mpf_down->GetBinContent(i);
      }
      else if(pt_depend_logpt_mpf->GetBinContent(i) < pt_depend_logpt_mpf_up->GetBinContent(i) && pt_depend_logpt_mpf->GetBinContent(i) < pt_depend_logpt_mpf_down->GetBinContent(i)){
	if(pt_depend_logpt_mpf_up->GetBinContent(i) < pt_depend_logpt_mpf_down->GetBinContent(i)) mpf_logpt_sys_up = pt_depend_logpt_mpf_down->GetBinContent(i) - pt_depend_logpt_mpf->GetBinContent(i);
	else mpf_logpt_sys_up = pt_depend_logpt_mpf_up->GetBinContent(i) - pt_depend_logpt_mpf->GetBinContent(i);
      }
      else {throw runtime_error("what did i do - mpf logpt?!");}
 
      //dijet_const
      if(pt_depend_const_dijet->GetBinContent(i) > pt_depend_const_dijet_up->GetBinContent(i) && pt_depend_const_dijet->GetBinContent(i) < pt_depend_const_dijet_down->GetBinContent(i)) {
	dijet_const_sys_up = pt_depend_const_dijet_down->GetBinContent(i) - pt_depend_const_dijet->GetBinContent(i);
	dijet_const_sys_down = pt_depend_const_dijet->GetBinContent(i) - pt_depend_const_dijet_up->GetBinContent(i);
      } 
      else if(pt_depend_const_dijet->GetBinContent(i) < pt_depend_const_dijet_up->GetBinContent(i) && pt_depend_const_dijet->GetBinContent(i) > pt_depend_const_dijet_down->GetBinContent(i)){
	dijet_const_sys_up = pt_depend_const_dijet_up->GetBinContent(i) - pt_depend_const_dijet->GetBinContent(i);
	dijet_const_sys_down = pt_depend_const_dijet->GetBinContent(i) - pt_depend_const_dijet_down->GetBinContent(i);
      }
      else if(pt_depend_const_dijet->GetBinContent(i) > pt_depend_const_dijet_up->GetBinContent(i) && pt_depend_const_dijet->GetBinContent(i) > pt_depend_const_dijet_down->GetBinContent(i)){
	if(pt_depend_const_dijet_up->GetBinContent(i) < pt_depend_const_dijet_down->GetBinContent(i)) dijet_const_sys_down = pt_depend_const_dijet->GetBinContent(i) - pt_depend_const_dijet_up->GetBinContent(i);
	else dijet_const_sys_down = pt_depend_const_dijet->GetBinContent(i) - pt_depend_const_dijet_down->GetBinContent(i);
      }
      else if(pt_depend_const_dijet->GetBinContent(i) < pt_depend_const_dijet_up->GetBinContent(i) && pt_depend_const_dijet->GetBinContent(i) < pt_depend_const_dijet_down->GetBinContent(i)){
	if(pt_depend_const_dijet_up->GetBinContent(i) < pt_depend_const_dijet_down->GetBinContent(i)) dijet_const_sys_up = pt_depend_const_dijet_down->GetBinContent(i) - pt_depend_const_dijet->GetBinContent(i);
	else dijet_const_sys_up = pt_depend_const_dijet_up->GetBinContent(i) - pt_depend_const_dijet->GetBinContent(i);
      }
      else {throw runtime_error("what did i do?! - dijet const");}

      //again for dijet_logpt
      if(pt_depend_logpt_dijet->GetBinContent(i) > pt_depend_logpt_dijet_up->GetBinContent(i) && pt_depend_logpt_dijet->GetBinContent(i) < pt_depend_logpt_dijet_down->GetBinContent(i)) {
	dijet_logpt_sys_up = pt_depend_logpt_dijet_down->GetBinContent(i) - pt_depend_logpt_dijet->GetBinContent(i);
	dijet_logpt_sys_down = pt_depend_logpt_dijet->GetBinContent(i) - pt_depend_logpt_dijet_up->GetBinContent(i);
      } 
      else if(pt_depend_logpt_dijet->GetBinContent(i) < pt_depend_logpt_dijet_up->GetBinContent(i) && pt_depend_logpt_dijet->GetBinContent(i) > pt_depend_logpt_dijet_down->GetBinContent(i)){
	dijet_logpt_sys_up = pt_depend_logpt_dijet_up->GetBinContent(i) - pt_depend_logpt_dijet->GetBinContent(i);
	dijet_logpt_sys_down = pt_depend_logpt_dijet->GetBinContent(i) - pt_depend_logpt_dijet_down->GetBinContent(i);
      }
      else if(pt_depend_logpt_dijet->GetBinContent(i) > pt_depend_logpt_dijet_up->GetBinContent(i) && pt_depend_logpt_dijet->GetBinContent(i) > pt_depend_logpt_dijet_down->GetBinContent(i)){
	if(pt_depend_logpt_dijet_up->GetBinContent(i) < pt_depend_logpt_dijet_down->GetBinContent(i)) dijet_logpt_sys_down = pt_depend_logpt_dijet->GetBinContent(i) - pt_depend_logpt_dijet_up->GetBinContent(i);
	else dijet_logpt_sys_down = pt_depend_logpt_dijet->GetBinContent(i) - pt_depend_logpt_dijet_down->GetBinContent(i);
      }
      else if(pt_depend_logpt_dijet->GetBinContent(i) < pt_depend_logpt_dijet_up->GetBinContent(i) && pt_depend_logpt_dijet->GetBinContent(i) < pt_depend_logpt_dijet_down->GetBinContent(i)){
	if(pt_depend_logpt_dijet_up->GetBinContent(i) < pt_depend_logpt_dijet_down->GetBinContent(i)) dijet_logpt_sys_up = pt_depend_logpt_dijet_down->GetBinContent(i) - pt_depend_logpt_dijet->GetBinContent(i);
	else dijet_logpt_sys_up = pt_depend_logpt_dijet_up->GetBinContent(i) - pt_depend_logpt_dijet->GetBinContent(i);
      }
      else {throw runtime_error("what did i do - dijet logpt?!");}
 

      //Set arrays
      mpf_const_up[i-1] = mpf_const_sys_up;
      mpf_const_down[i-1] = mpf_const_sys_down;
      mpf_logpt_up[i-1] = mpf_logpt_sys_up;
      mpf_logpt_down[i-1] = mpf_logpt_sys_down;
      dijet_const_up[i-1] = dijet_const_sys_up;
      dijet_const_down[i-1] = dijet_const_sys_down;
      dijet_logpt_up[i-1] = dijet_logpt_sys_up;
      dijet_logpt_down[i-1] = dijet_logpt_sys_down;
      cout << "bin " << i << "---" << "mpf_const_sys_up: " << mpf_const_sys_up << ", mpf_const_sys_down: " << mpf_const_sys_down << ", mpf_logpt_sys_up:" << mpf_logpt_sys_up << ", mpf_logpt_sys_down:" << mpf_logpt_sys_down << ", dijet_const_sys_up:" <<dijet_const_sys_up << ", dijet_const_sys_down:" <<  dijet_const_sys_down << ", dijet_logpt_sys_up:" << dijet_logpt_sys_up << ", dijet_logpt_sys_down:" << dijet_logpt_sys_down << endl;


      cout << "for logpt dijet: nominal: " << pt_depend_logpt_dijet->GetBinContent(i) << ", up = " << pt_depend_logpt_dijet_up->GetBinContent(i) << ", down: " << pt_depend_logpt_dijet_down->GetBinContent(i) << endl << endl;
    }
 

    for(int i=0; i<n_bins_x; i++){
      //Calculate total error
      gr_pt_depend_const_mpf->SetPointEYhigh(i,sqrt(pow(mpf_const_up[i],2) + pow(gr_pt_depend_const_mpf->GetErrorYhigh(i),2)));
      gr_pt_depend_const_mpf->SetPointEYlow(i,sqrt(pow(mpf_const_down[i],2) + pow(gr_pt_depend_const_mpf->GetErrorYlow(i),2)));
      gr_pt_depend_const_dijet->SetPointEYhigh(i,sqrt(pow(dijet_const_up[i],2) + pow(gr_pt_depend_const_dijet->GetErrorYhigh(i),2)));
      gr_pt_depend_const_dijet->SetPointEYlow(i,sqrt(pow(dijet_const_down[i],2) + pow(gr_pt_depend_const_dijet->GetErrorYlow(i),2)));
      gr_pt_depend_logpt_mpf->SetPointEYhigh(i,sqrt(pow(mpf_logpt_up[i],2) + pow(gr_pt_depend_logpt_mpf->GetErrorYhigh(i),2)));
      gr_pt_depend_logpt_mpf->SetPointEYlow(i,sqrt(pow(mpf_logpt_down[i],2) + pow(gr_pt_depend_logpt_mpf->GetErrorYlow(i),2)));
      gr_pt_depend_logpt_dijet->SetPointEYhigh(i,sqrt(pow(dijet_logpt_up[i],2) + pow(gr_pt_depend_logpt_dijet->GetErrorYhigh(i),2)));
      gr_pt_depend_logpt_dijet->SetPointEYlow(i,sqrt(pow(dijet_logpt_down[i],2) + pow(gr_pt_depend_logpt_dijet->GetErrorYlow(i),2)));

      cout << "actual errors (down only)" << endl << "bin " << i << "gr_pt_depend_logpt_dijet" << gr_pt_depend_logpt_dijet->GetErrorYlow(i) << endl;
    }
  }//if(use_sys_hists)


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
  alVal.Form("%0.2f\n",al_cut);
  TString altitle = "{#alpha<"+alVal+"}";
  TString axistitle = "(R^{MC}/R^{data})_";
  axistitle +=altitle;
  h->GetYaxis()->SetTitle(axistitle);
  h->GetYaxis()->SetRangeUser(0.81,1.15);

  if(use_sys_hists){
    //pt_depend_const_mpf->Draw("E1 SAME");
    //pt_depend_const_dijet->Draw("E1 SAME");
    pt_depend_logpt_mpf->Draw("E1 SAME");
    pt_depend_logpt_dijet->Draw("E1 SAME");

    //gr_pt_depend_const_mpf->Draw("P SAME");
    gr_pt_depend_logpt_dijet->Draw("P SAME");
    gr_pt_depend_logpt_mpf->Draw("P SAME");
    //gr_pt_depend_const_dijet->Draw("P SAME");
  }
  else{
    pt_depend_const_mpf->SetMarkerStyle(1);
    pt_depend_const_dijet->SetMarkerStyle(1);
    pt_depend_logpt_mpf->SetMarkerStyle(1);
    pt_depend_logpt_dijet->SetMarkerStyle(1);
    pt_depend_const_mpf->Draw("E1 SAME");
    pt_depend_const_dijet->Draw("E1 SAME");
    pt_depend_logpt_mpf->Draw("E1 SAME");
    pt_depend_logpt_dijet->Draw("E1 SAME");
  }

  line->Draw("SAME");
  leg1->Draw();
  tex->DrawLatex(0.45,0.87,JetDescrib);
  if(use_sys_hists){
    tex_bars->DrawLatex(0.45,0.37,"Inner: stat");
    tex_bars->DrawLatex(0.45,0.32,"Outer: stat #oplus JER");
  }
  if(divide_by_lumi)c2->SaveAs(path+"plots/Ratio_"+jettag+"_"+txttag+"_Lumi_"+lumi_range[current_lumibin-1]+"_"+lumi_range[current_lumibin]+".pdf"); 
  else c2->SaveAs(path+"plots/Ratio_"+jettag+"_"+txttag+".pdf");


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

  //TLegend *leg2 = tdrLeg(0.17,0.49,0.40,0.80);
  TLegend *leg2 = tdrLeg(0.17,0.19,0.40,0.30);
  leg2 -> AddEntry(kfsr_mpf, "MPF","L");
  leg2 -> AddEntry(kfsr_dijet, "Pt","L");
  leg2->Draw();
  tex->DrawLatex(0.45,0.87,JetDescrib);
  if(divide_by_lumi) c3->SaveAs(path+"plots/kFSR_"+jettag+"_"+txttag+"_Lumi_"+lumi_range[current_lumibin-1]+"_"+lumi_range[current_lumibin]+".pdf");
  else c3->SaveAs(path+"plots/kFSR_"+jettag+"_"+txttag+".pdf");


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
  if(divide_by_lumi) c4->SaveAs(path+"plots/L2Res_kFSRfit_"+jettag+"_"+txttag+"_Lumi_"+lumi_range[current_lumibin-1]+"_"+lumi_range[current_lumibin]+".pdf");
  else c4->SaveAs(path+"plots/L2Res_kFSRfit_"+jettag+"_"+txttag+".pdf");

  //pt-dependence of L2Res (Money plot)

  TString var="";
  TH1D* res_logpt_mpf_kfsrfit_var[4];
  TH1D* res_logpt_dijet_kfsrfit_var[4];
  for(int i=0;i<4;i++){
    if(i==0) var="central"; 
    if(i==1) var="down";
    if(i==2) var="up";
    if(i==3) var="doubleup";
    //    cout<<var<<endl; 
    TFile* f_Res_mpf_var;
    TFile* f_Res_dijet_var;
    if(divide_by_lumi) f_Res_mpf_var = new TFile(path+"Histo_Res_MPF_L1_"+txttag+"_"+jettag+"_"+var+"_Lumi_"+lumi_range[current_lumibin-1]+"_"+lumi_range[current_lumibin]+".root","READ"); 
    else f_Res_mpf_var = new TFile(path+"Histo_Res_MPF_L1_"+txttag+"_"+jettag+"_"+var+".root","READ");   
    if(divide_by_lumi) f_Res_dijet_var = new TFile(path+"Histo_Res_DiJet_L1_"+txttag+"_"+jettag+"_"+var+"_Lumi_"+lumi_range[current_lumibin-1]+"_"+lumi_range[current_lumibin]+".root","READ");
    else f_Res_dijet_var = new TFile(path+"Histo_Res_DiJet_L1_"+txttag+"_"+jettag+"_"+var+".root","READ"); 
 
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
  if(divide_by_lumi) c5->SaveAs(path+"plots/L2Res_logpt_MPF_kFSRfit_"+jettag+"_"+txttag+"_Lumi_"+lumi_range[current_lumibin-1]+"_"+lumi_range[current_lumibin]+".pdf");
  else c5->SaveAs(path+"plots/L2Res_logpt_MPF_kFSRfit_"+jettag+"_"+txttag+".pdf");


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
  if(divide_by_lumi) c6->SaveAs(path+"plots/L2Res_logpt_DiJet_kFSRfit_"+jettag+"_"+txttag+"_Lumi_"+lumi_range[current_lumibin-1]+"_"+lumi_range[current_lumibin]+".pdf");
  else c6->SaveAs(path+"plots/L2Res_logpt_DiJet_kFSRfit_"+jettag+"_"+txttag+".pdf");




  f_Res_mpf->Close();
  //f_Res_mpf_up->Close(); 
  //f_Res_mpf_down->Close();
  f_Res_dijet->Close();
  //f_Res_dijet_up->Close();
  //f_Res_dijet_down->Close();


}


