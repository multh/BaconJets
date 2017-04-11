#include "../include/parameters.h"
#include "../include/useful_functions.h"
#include "../include/CorrectionObject.h"
#include "../include/tdrstyle_mod15.h"

#include <TStyle.h>
#include <TF1.h>
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


using namespace std;

void CorrectionObject::FinalControlPlots_CorrectFormulae(){
  cout << "--------------- Starting FinalControlPlots_CorrectFormulae() ---------------" << endl << endl;
  gStyle->SetOptStat(0);

  //Table with number of events in each pT- and eta-bin
  
  //Set up histos for ratios of responses
  TH1D *hdata_asymmetry[n_pt-1][n_eta-1]; // A for data
  TH1D *hdata_B[n_pt-1][n_eta-1];         // B for data
  TH1D *hdata_METoverJetsPt[n_pt-1][n_eta-1];         // MET/sum_jets_pt for data
  TH1D *hmc_asymmetry[n_pt-1][n_eta-1];   // A for MC
  TH1D *hmc_B[n_pt-1][n_eta-1];           // B for MC
  TH1D *hmc_METoverJetsPt[n_pt-1][n_eta-1];         // MET/sum_jets_pt for MC
  TH1D* hmc_pt_ave[n_eta-1];              // pt_ave for MC
  TH1D* hdata_pt_ave[n_eta-1];            // pt_ave for data


  int count = 0;
  TString name1 = "hist_data_A_";
  TString name2 = "hist_data_B_";
  TString name3 = "hist_mc_A_";
  TString name4 = "hist_mc_B_";
  TString name5 = "hist_data_METoverJetsPt_";
  TString name6 = "hist_mc_METoverJetsPt_";
  for(int j=0; j<n_eta-1; j++){
      TString eta_name = "eta_"+eta_range2[j]+"_"+eta_range2[j+1];
    for(int k=0; k<n_pt-1; k++){
      TString pt_name = "pt_"+pt_range[k]+"_"+pt_range[k+1];

      
      TString name = name1 + eta_name + "_" + pt_name; 
      hdata_asymmetry[k][j] = new TH1D(name,"",nResponseBins, -1.2, 1.2);
      name = name2 + eta_name + "_" + pt_name;
      hdata_B[k][j] = new TH1D(name,"",nResponseBins, -1.2, 1.2);
      name = name5 + eta_name + "_" + pt_name;
      hdata_METoverJetsPt[k][j] = new TH1D(name,"",60,0,1.2);
      // hdata_METoverJetsPt[k][j]->Print();
      name = name3 + eta_name + "_" + pt_name;
      hmc_asymmetry[k][j] = new TH1D(name,"",nResponseBins, -1.2, 1.2);
      name = name4 + eta_name + "_" + pt_name;
      hmc_B[k][j] = new TH1D(name,"",nResponseBins, -1.2, 1.2);
      name = name6 + eta_name + "_" + pt_name;
      hmc_METoverJetsPt[k][j] = new TH1D(name,"",50,0,1.2);
      //      hmc_METoverJetsPt[k][j]->Print();
      /*
      TString name = name1 + eta_name + "_" + pt_name; 
      hdata_asymmetry[k][j] = new TH1D(name,"",nResponseBins, 0, 2.5);
      name = name2 + eta_name + "_" + pt_name;
      hdata_B[k][j] = new TH1D(name,"",nResponseBins, 0, 2.5);
      name = name3 + eta_name + "_" + pt_name;
      hmc_asymmetry[k][j] = new TH1D(name,"",nResponseBins, 0, 2.5);
      name = name4 + eta_name + "_" + pt_name;
      hmc_B[k][j] = new TH1D(name,"",nResponseBins, 0, 2.5);
      
      TString name = name1 + eta_name + "_" + pt_name; 
      hdata_asymmetry[k][j] = new TH1D(name,"",nResponseBins, -2.5, 0);
      name = name2 + eta_name + "_" + pt_name;
      hdata_B[k][j] = new TH1D(name,"",nResponseBins, -2.5, 0);
      name = name3 + eta_name + "_" + pt_name;
      hmc_asymmetry[k][j] = new TH1D(name,"",nResponseBins, -2.5, 0);
      name = name4 + eta_name + "_" + pt_name;
      hmc_B[k][j] = new TH1D(name,"",nResponseBins, -2.5, 0);
      */
      count++;
    }

    //define pt_ave[eta] histos
    TString name_pt_ave = "hist_pt_ave_";
    hmc_pt_ave[j] = new TH1D(name_pt_ave+"MC_"+eta_name,"",1000,0,5000);
    hdata_pt_ave[j] = new TH1D(name_pt_ave+"data_"+eta_name,"",1000,0,5000);

  }


  //Get relevant information from DATA, loop over DATA events
  TTreeReader myReader_DATA("AnalysisTree", CorrectionObject::_DATAFile);
  TTreeReaderValue<Float_t> pt_ave_data(myReader_DATA, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_data(myReader_DATA, "probejet_eta");
  TTreeReaderValue<Float_t> probejet_pt_data(myReader_DATA, "probejet_pt");
  TTreeReaderValue<Float_t> barreljet_pt_data(myReader_DATA, "barreljet_pt");
  TTreeReaderValue<Float_t> alpha_data(myReader_DATA, "alpha");
  TTreeReaderValue<Float_t> asymmetry_data(myReader_DATA, "asymmetry");
  TTreeReaderValue<Float_t> B_data(myReader_DATA, "B");
  TTreeReaderValue<Float_t> weight_data(myReader_DATA, "weight");
  TTreeReaderValue<Float_t> MET_data(myReader_DATA, "MET");
  TTreeReaderValue<Float_t> sum_jets_pt_data(myReader_DATA, "sum_jets_pt");
   
  while (myReader_DATA.Next()) {
    if(*alpha_data>alpha_cut) continue;
    //fill histos in bins of eta
    for(int i=0; i<n_eta-1; i++){
      if(fabs(*probejet_eta_data)<eta_bins[i+1] && fabs(*probejet_eta_data)>=eta_bins[i]){
	hdata_pt_ave[i]->Fill(*pt_ave_data,*weight_data);
      }
    }

    //fill histos in bins of pt and eta
    for(int k=0; k<n_pt-1; k++){
      if(*pt_ave_data<pt_bins[k] || *pt_ave_data>pt_bins[k+1]) continue;
      for(int j=0; j<n_eta-1; j++){
	if(fabs(*probejet_eta_data)>eta_bins[j+1] || fabs(*probejet_eta_data)<eta_bins[j]) continue;
	else{
	  hdata_asymmetry[k][j]->Fill(*asymmetry_data,*weight_data);
	  hdata_B[k][j]->Fill(*B_data,*weight_data);
	  hdata_METoverJetsPt[k][j]->Fill((*MET_data)/(*sum_jets_pt_data+*probejet_pt_data+*barreljet_pt_data),*weight_data);
	}
      }
    }
  }


  //Get relevant information from MC, loop over MC events 
  TTreeReader myReader_MC("AnalysisTree", CorrectionObject::_MCFile);
  TTreeReaderValue<Float_t> pt_ave_mc(myReader_MC, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_mc(myReader_MC, "probejet_eta");
  TTreeReaderValue<Float_t> probejet_pt_mc(myReader_MC, "probejet_pt");
  TTreeReaderValue<Float_t> barreljet_pt_mc(myReader_MC, "barreljet_pt");
  TTreeReaderValue<Float_t> alpha_mc(myReader_MC, "alpha");
  TTreeReaderValue<Float_t> asymmetry_mc(myReader_MC, "asymmetry");
  TTreeReaderValue<Float_t> B_mc(myReader_MC, "B");
  TTreeReaderValue<Float_t> weight_mc(myReader_MC, "weight");
  TTreeReaderValue<Float_t> MET_mc(myReader_MC, "MET");
  TTreeReaderValue<Float_t> sum_jets_pt_mc(myReader_MC, "sum_jets_pt");
  while (myReader_MC.Next()) {
    if(*alpha_mc>alpha_cut) continue;
    //fill histos in bins of eta
    for(int i=0; i<n_eta-1; i++){
      if(fabs(*probejet_eta_mc)<eta_bins[i+1] && fabs(*probejet_eta_mc)>=eta_bins[i]){
	hmc_pt_ave[i]->Fill(*pt_ave_mc,*weight_mc);
      }
    }

    //fill histos in bins of pt and eta
    for(int k=0; k<n_pt-1; k++){
      if(*pt_ave_mc<pt_bins[k] || *pt_ave_mc>pt_bins[k+1]) continue;
      for(int j=0; j<n_eta-1; j++){
	if(fabs(*probejet_eta_mc)>eta_bins[j+1] || fabs(*probejet_eta_mc)<eta_bins[j]) continue;
	else{
	  hmc_asymmetry[k][j]->Fill(*asymmetry_mc,*weight_mc);
	  hmc_B[k][j]->Fill(*B_mc,*weight_mc);
	  hmc_METoverJetsPt[k][j]->Fill((*MET_mc)/(*sum_jets_pt_mc+*probejet_pt_mc+*barreljet_pt_mc),*weight_mc);
	  //	  hmc_METoverJetsPt[k][j]->Print();
	}
      }
    }
  }
  




 

  ofstream output;
  output.open(CorrectionObject::_outpath+"plots/control/Number_Events_Pt_Eta_bins_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt");
    

  output << "Number of events in each bin for MC" << endl;
  output << "|Eta|:          ";
  double n_tot_MC = 0;
  double n_tot_DATA = 0;
  for(int i=0; i<n_eta; i++) {
    if(i != n_eta-1) output << eta_range[i] << " -- ";
    else output << eta_range[i] << endl;
  }
  for(int i=0; i<n_pt-1; i++){
    if(i==0) output << "pT = ["  << fixed << setprecision(0) << pt_bins[i] << "," << setprecision(0) << pt_bins[i+1] << "]  :    ";
    else if(i==1) output << "pT = ["  << fixed << setprecision(0) << pt_bins[i] << "," << setprecision(0) << pt_bins[i+1] << "] :    ";
    else output << "pT = ["  << fixed << setprecision(0) << pt_bins[i] << "," << setprecision(0) << pt_bins[i+1] << "]:    ";

    for(int j=0; j<n_eta-1; j++){
      if(j!=n_eta-2){
	if(hmc_B[i][j]->GetEntries()/1000 < 0.01)     output << hmc_B[i][j]->GetEntries() << "      - "; //<1000
	else if(hmc_B[i][j]->GetEntries()/1000 < 0.1) output << hmc_B[i][j]->GetEntries() << "     - "; //<1000
	else if(hmc_B[i][j]->GetEntries()/1000 < 1)   output << hmc_B[i][j]->GetEntries() << "    - "; //<1000
	else if(hmc_B[i][j]->GetEntries()/1000 <10)   output << hmc_B[i][j]->GetEntries() << "   - "; //<10000
	else if(hmc_B[i][j]->GetEntries()/1000 <100)  output << hmc_B[i][j]->GetEntries() << "  - ";
	else                                              output << hmc_B[i][j]->GetEntries() << " - ";
      }
      else output << hmc_B[i][j]->GetEntries() << endl;
      n_tot_MC+= hmc_B[i][j]->GetEntries();
    }

  }
  output << endl << endl << "Total number of events in MC: " << n_tot_MC << endl;

  output << endl << endl << endl << endl << "Number of events in each bin for DATA" << endl;
  output << "|Eta|:          ";
  for(int i=0; i<n_eta; i++) {
    if(i != n_eta-1) output << eta_range[i] << " -- ";
    else output << eta_range[i] << endl;
  }
  for(int i=0; i<n_pt-1; i++){
    if(i==0) output << "pT = ["  << fixed << setprecision(0) << pt_bins[i] << "," << setprecision(0) << pt_bins[i+1] << "]  :    ";
    else if(i==1) output << "pT = ["  << fixed << setprecision(0) << pt_bins[i] << "," << setprecision(0) << pt_bins[i+1] << "] :    ";
    else output << "pT = ["  << fixed << setprecision(0) << pt_bins[i] << "," << setprecision(0) << pt_bins[i+1] << "]:    ";

    for(int j=0; j<n_eta-1; j++){
      if(j!=n_eta-2){
	if(hdata_B[i][j]->GetEntries()/1000 < 0.01)     output << hdata_B[i][j]->GetEntries() << "      - "; //<1000
	else if(hdata_B[i][j]->GetEntries()/1000 < 0.1) output << hdata_B[i][j]->GetEntries() << "     - "; //<1000
	else if(hdata_B[i][j]->GetEntries()/1000 < 1)   output << hdata_B[i][j]->GetEntries() << "    - "; //<1000
	else if(hdata_B[i][j]->GetEntries()/1000 <10)   output << hdata_B[i][j]->GetEntries() << "   - "; //<10000
	else if(hdata_B[i][j]->GetEntries()/1000 <100)  output << hdata_B[i][j]->GetEntries() << "  - ";
	else                                                output << hdata_B[i][j]->GetEntries() << " - ";
      }
      else output << hdata_B[i][j]->GetEntries() << endl;
      n_tot_DATA += hdata_B[i][j]->GetEntries();
    }

  }
  output << endl << endl << "Total number of events in DATA: " << n_tot_DATA << endl;



  // Dump 1-d distributions of A and B in bins of pT, eta

  TFile* test_out_mc_B = new TFile(CorrectionObject::_outpath+"plots/control/B_1d_mc.root","RECREATE");
  for(int j=0; j<n_eta-1; j++){
    for(int k=0; k<n_pt-1; k++){
      hmc_B[k][j]->Write();
      hmc_METoverJetsPt[k][j]->Write();
    }
  }
  test_out_mc_B->Close();
  delete test_out_mc_B;

  TFile* test_out_data_B = new TFile(CorrectionObject::_outpath+"plots/control/B_1d_data.root","RECREATE");
  for(int j=0; j<n_eta-1; j++){
    for(int k=0; k<n_pt-1; k++){
      hdata_B[k][j]->Write();
      hdata_METoverJetsPt[k][j]->Write();
    }
  }
  test_out_data_B->Close();
  delete test_out_data_B;

  TFile* test_out_mc_A = new TFile(CorrectionObject::_outpath+"plots/control/A_1d_mc.root","RECREATE");
  for(int j=0; j<n_eta-1; j++){
    for(int k=0; k<n_pt-1; k++){
      hmc_asymmetry[k][j]->Write();
    }
  }
  test_out_mc_A->Close();
  delete test_out_mc_A;

  TFile* test_out_data_A = new TFile(CorrectionObject::_outpath+"plots/control/A_1d_data.root","RECREATE");
  for(int j=0; j<n_eta-1; j++){
    for(int k=0; k<n_pt-1; k++){
      hdata_asymmetry[k][j]->Write();
    }
  }
  test_out_data_A->Close();
  delete test_out_data_A;




  //R_MC and R_DATA overlaid in the same plot as a function of pT, in bins of |eta|
  double val_rel_mc[n_eta-1][n_pt-1]; //value at pt,eta
  double err_rel_mc[n_eta-1][n_pt-1]; //error of ratio at pt,eta
  double val_mpf_mc[n_eta-1][n_pt-1]; //ratio at pt,eta
  double err_mpf_mc[n_eta-1][n_pt-1]; //error of ratio at pt,eta
  double val_rel_data[n_eta-1][n_pt-1]; //value at pt,eta
  double err_rel_data[n_eta-1][n_pt-1]; //error of ratio at pt,eta
  double val_mpf_data[n_eta-1][n_pt-1]; //ratio at pt,eta
  double err_mpf_data[n_eta-1][n_pt-1]; //error of ratio at pt,eta



  for(int i=0; i<n_eta-1; i++){
    for(int j=0; j<n_pt-1; j++){

      //get <A> and error on <A>
      pair <double,double> A_mc = GetValueAndError(hmc_asymmetry[j][i]);
      pair <double,double> A_data = GetValueAndError(hdata_asymmetry[j][i]);
      pair <double,double> B_mc = GetValueAndError(hmc_B[j][i]);
      pair <double,double> B_data = GetValueAndError(hdata_B[j][i]);

      //build MPF and pt_bal and their errors


      pair<double,double> res_mc_rel_r,res_data_rel_r;
      pair<double,double> res_mc_mpf_r,res_data_mpf_r;
      res_mc_mpf_r.first = (1+B_mc.first)/(1-B_mc.first);
      res_mc_mpf_r.second = 2/(pow((1-B_mc.first),2)) * B_mc.second;
      res_data_mpf_r.first = (1+B_data.first)/(1-B_data.first);
      res_data_mpf_r.second = 2/(pow((1-B_data.first),2)) * B_data.second;
      res_mc_rel_r.first = (1+A_mc.first)/(1-A_mc.first);
      res_mc_rel_r.second = 2/(pow((1-A_mc.first),2)) * A_mc.second;
      res_data_rel_r.first = (1+A_data.first)/(1-A_data.first);
      res_data_rel_r.second = 2/(pow((1-A_data.first),2)) * A_data.second;


      val_rel_mc[i][j] = res_mc_rel_r.first;
      err_rel_mc[i][j] = res_mc_rel_r.second;
      val_mpf_mc[i][j] = res_mc_mpf_r.first;
      err_mpf_mc[i][j] = res_mc_mpf_r.second;
      val_rel_data[i][j] = res_data_rel_r.first;
      err_rel_data[i][j] = res_data_rel_r.second;
      val_mpf_data[i][j] = res_data_mpf_r.first;
      err_mpf_data[i][j] = res_data_mpf_r.second;
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


 

    TGraphErrors *graph_mpf_mc   = new TGraphErrors(n_pt-1, xbin_tgraph, val_mpf_mc[i], zero, err_mpf_mc[i]);
    TGraphErrors *graph_mpf_data = new TGraphErrors(n_pt-1, xbin_tgraph, val_mpf_data[i], zero, err_mpf_data[i]);
    TGraphErrors *graph_rel_mc   = new TGraphErrors(n_pt-1, xbin_tgraph, val_rel_mc[i], zero, err_rel_mc[i]);
    TGraphErrors *graph_rel_data = new TGraphErrors(n_pt-1, xbin_tgraph, val_rel_data[i], zero, err_rel_data[i]);
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
    //    TString axistitle_mc = "R^{PromtReco}_";
    //    TString axistitle_mc = "R^{23SepReReco}_";
    //    TString axistitle_mc = "R^{03FebReMINIAOD}_";
    TString axistitle_data = "R^{DATA}_";
    //    TString axistitle_data = "R^{preLegacy}_";
    // TString axistitle_data = "R^{23SepReReco}_";
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





  //Plot 1d response distributions in a particular eta-bin for different pt-bins onto a single canvas

  //Get histo files
  TFile* f_mpf_mc = new TFile(CorrectionObject::_outpath+"plots/control/B_1d_mc.root","READ");
  TFile* f_mpf_data = new TFile(CorrectionObject::_outpath+"plots/control/B_1d_data.root","READ");
  TFile* f_rel_mc = new TFile(CorrectionObject::_outpath+"plots/control/A_1d_mc.root","READ");
  TFile* f_rel_data = new TFile(CorrectionObject::_outpath+"plots/control/A_1d_data.root","READ");
  for(int i=0; i<n_eta-1; i++){
    TString eta_name = "eta_"+eta_range2[i]+"_"+eta_range2[i+1];
    
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    TString text = eta_range[i] + " < |#eta| < " + eta_range[i+1];

    TLatex *tex_lumi = new TLatex();
    tex_lumi->SetNDC();
    tex_lumi->SetTextSize(0.045); 
    

    TCanvas* c1 = new TCanvas();
    tdrCanvas(c1,"c1",h,4,10,kSquare,"MC");
    TLegend leg1 = tdrLeg(0.62,0.46,0.85,0.81);
    TH1D* htemp_mpf_mc;
    for(int j=0; j<n_pt-1; j++){
    //    for(int j=0; j<5; j++){ //TEST
    //    for(int j=5; j<n_pt-1; j++){//TEST
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_mpf_mc = "hist_mc_B_"+eta_name+"_"+pt_name;
      htemp_mpf_mc = (TH1D*)f_mpf_mc->Get(name_mpf_mc);
      int n_ev = htemp_mpf_mc->GetEntries();
      if(htemp_mpf_mc->Integral() > 0)htemp_mpf_mc->Scale(1/htemp_mpf_mc->Integral());
      h->GetXaxis()->SetTitle("B");
      h->GetYaxis()->SetTitle("Normalized entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(-1.2,1.2);
      h->SetMaximum(0.1);
      if(j<9) htemp_mpf_mc->SetLineColor(j+1);
      else    htemp_mpf_mc->SetLineColor(j+31);
      htemp_mpf_mc->SetLineWidth(3);
      if(n_ev>100) htemp_mpf_mc->Draw("HIST SAME");
      leg1.AddEntry(htemp_mpf_mc, legname, "l");
    }

    leg1.Draw();
    tex->DrawLatex(0.52,0.85,"MC, " + text);
    //tex_lumi->DrawLatex(0.60,0.91,"MC");
    c1->SaveAs(CorrectionObject::_outpath+"plots/control/B_NormDistribution_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + ".pdf");


    TCanvas* c2 = new TCanvas();
    tdrCanvas(c2,"c2",h,4,10,kSquare,CorrectionObject::_lumitag);
    TLegend leg2 = tdrLeg(0.62,0.46,0.85,0.81);
    TH1D* htemp_mpf_data;
    for(int j=0; j<n_pt-1; j++){
    //    for(int j=0; j<5; j++){ //TEST
    //    for(int j=5; j<n_pt-1; j++){//TEST
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_mpf_data = "hist_data_B_"+eta_name+"_"+pt_name;
      htemp_mpf_data = (TH1D*)f_mpf_data->Get(name_mpf_data);
      int n_ev = htemp_mpf_data->GetEntries();
      if(htemp_mpf_data->Integral() > 0)htemp_mpf_data->Scale(1/htemp_mpf_data->Integral());
      h->GetXaxis()->SetTitle("B");
      h->GetYaxis()->SetTitle("Normalized entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(-1.2,1.2);
      h->SetMaximum(0.1);
      if(j<9) htemp_mpf_data->SetLineColor(j+1);
      else    htemp_mpf_data->SetLineColor(j+31);
      htemp_mpf_data->SetLineWidth(3);
      if(n_ev>100) htemp_mpf_data->Draw("HIST SAME");
      leg2.AddEntry(htemp_mpf_data, legname ,"l");
    }

    leg2.Draw();
    tex->DrawLatex(0.52,0.85,"Data, " + text);
    //tex_lumi->DrawLatex(0.50,0.91,CorrectionObject::_lumitag+"(13TeV)");
    c2->SaveAs(CorrectionObject::_outpath+"plots/control/B_NormDistribution_DATA_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + ".pdf");

    TCanvas* c3 = new TCanvas();
    tdrCanvas(c3,"c3",h,4,10,kSquare,"MC");
    TLegend leg3 = tdrLeg(0.62,0.46,0.85,0.81);
    TH1D* htemp_rel_mc;
    for(int j=0; j<n_pt-1; j++){
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_rel_mc = "hist_mc_A_"+eta_name+"_"+pt_name;
      htemp_rel_mc = (TH1D*)f_rel_mc->Get(name_rel_mc);
      int n_ev =  htemp_rel_mc->GetEntries();
      if(htemp_rel_mc->Integral() > 0)htemp_rel_mc->Scale(1/htemp_rel_mc->Integral());
      h->GetXaxis()->SetTitle("A");
      h->GetYaxis()->SetTitle("Normalized entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(-1.2,1.2);
      h->SetMaximum(0.1);
      if(j<9) htemp_rel_mc->SetLineColor(j+1);
      else    htemp_rel_mc->SetLineColor(j+31);
      htemp_rel_mc->SetLineWidth(3);
      if(n_ev>100) htemp_rel_mc->Draw("HIST SAME");
      leg3.AddEntry(htemp_rel_mc, legname);
    }

    leg3.Draw();
    tex->DrawLatex(0.52,0.85,"MC, " + text);
    //tex_lumi->DrawLatex(0.6,0.91,"MC");
    c3->SaveAs(CorrectionObject::_outpath+"plots/control/A_NormDistribution_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + ".pdf");

    TCanvas* c4 = new TCanvas();
    tdrCanvas(c4,"c4",h,4,10,kSquare,CorrectionObject::_lumitag);
    TLegend leg4 = tdrLeg(0.62,0.46,0.85,0.81);
    TH1D* htemp_rel_data;
    for(int j=0; j<n_pt-1; j++){
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_rel_data = "hist_data_A_"+eta_name+"_"+pt_name;
      htemp_rel_data = (TH1D*)f_rel_data->Get(name_rel_data);
      int n_ev = htemp_rel_data->GetEntries();
      if(htemp_rel_data->Integral() > 0)htemp_rel_data->Scale(1/htemp_rel_data->Integral());
      h->GetXaxis()->SetTitle("A");
      h->GetYaxis()->SetTitle("Normalized entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(-1.2,1.2);
      h->SetMaximum(0.1);
      if(j<9) htemp_rel_data->SetLineColor(j+1);
      else    htemp_rel_data->SetLineColor(j+31);
      htemp_rel_data->SetLineWidth(3);
      if(n_ev>100) htemp_rel_data->Draw("HIST SAME");
      leg4.AddEntry(htemp_rel_data, legname);
    }
    leg4.Draw();
    tex->DrawLatex(0.52,0.85,"Data, " + text);
    //tex_lumi->DrawLatex(0.50,0.91,CorrectionObject::_lumitag+"(13TeV)");
    c4->SaveAs(CorrectionObject::_outpath+"plots/control/A_NormDistribution_DATA_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + ".pdf");    

    ///MET over sum pt
    TCanvas* c5 = new TCanvas();
    tdrCanvas(c5,"c5",h,4,10,kSquare,"MC");
    TLegend leg5 = tdrLeg(0.62,0.46,0.85,0.81);
    TH1D* htemp_met_mc;
    for(int j=0; j<n_pt-1; j++){
    //    for(int j=0; j<5; j++){ //TEST
    //    for(int j=5; j<n_pt-1; j++){//TEST
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_met_mc = "hist_mc_METoverJetsPt_"+eta_name+"_"+pt_name;
      htemp_met_mc = (TH1D*)f_mpf_mc->Get(name_met_mc);
      //      htemp_met_mc->Print();
      int n_ev =  htemp_met_mc->GetEntries();
      if(htemp_met_mc->Integral() > 0)htemp_met_mc->Scale(1/htemp_met_mc->Integral());
      h->GetXaxis()->SetTitle("MET/#sum p_{T}");
      h->GetYaxis()->SetTitle("Norm. Entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      // h->SetMaximum(0.3);
      h->GetXaxis()->SetLimits(0,1.2);
      //      h->GetYaxis()->SetLimits(0,0.8);
      h->SetMaximum(0.2);
      if(j<9) htemp_met_mc->SetLineColor(j+1);
      else    htemp_met_mc->SetLineColor(j+31);
      htemp_met_mc->SetLineWidth(3);
      if(n_ev>100) htemp_met_mc->Draw("HIST SAME");
      leg5.AddEntry(htemp_met_mc, legname);
    }

    leg5.Draw();
    tex->DrawLatex(0.52,0.85,"MC, " + text);
    //tex_lumi->DrawLatex(0.6,0.91,"MC");
    c5->SaveAs(CorrectionObject::_outpath+"plots/control/METoverPt_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + ".pdf");
    TCanvas* c6 = new TCanvas();
    tdrCanvas(c6,"c6",h,4,10,kSquare,CorrectionObject::_lumitag);
    TLegend leg6 = tdrLeg(0.62,0.46,0.85,0.81);
    TH1D* htemp_met_data;
    for(int j=0; j<n_pt-1; j++){
    //    for(int j=0; j<5; j++){ //TEST
    //    for(int j=5; j<n_pt-1; j++){//TEST
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_met_data = "hist_data_METoverJetsPt_"+eta_name+"_"+pt_name;
      htemp_met_data = (TH1D*)f_mpf_data->Get(name_met_data);
      int n_ev = htemp_met_data->GetEntries();
      if(htemp_met_data->Integral() > 0)htemp_met_data->Scale(1/htemp_met_data->Integral());
      h->GetXaxis()->SetTitle("MET/#sum p_{T}");
      h->GetYaxis()->SetTitle("Norm. Entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(0,1.2);
      h->GetYaxis()->SetLimits(0,0.8);
      h->SetMaximum(0.2);
      if(j<9) htemp_met_data->SetLineColor(j+1);
      else    htemp_met_data->SetLineColor(j+31);
      htemp_met_data->SetLineWidth(3);
      if(n_ev>100) htemp_met_data->Draw("HIST SAME");
      leg6.AddEntry(htemp_met_data, legname);
    }
    leg6.Draw();
    tex->DrawLatex(0.52,0.85,"Data, " + text);
    //tex_lumi->DrawLatex(0.50,0.91,CorrectionObject::_lumitag+"(13TeV)");
    c6->SaveAs(CorrectionObject::_outpath+"plots/control/METoverPt_DATA_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + ".pdf");    

    ///END MET over sum pt

    delete tex;
    delete htemp_rel_data;
    delete htemp_met_data;
    delete c4;
    delete htemp_rel_mc;
    delete htemp_met_mc;
    delete c3;
    delete htemp_mpf_data;
    delete c2;
    delete htemp_mpf_mc;
    delete c1;
  }





  //pT_ave for MC and data in bins of |eta|

  for(int i=0; i<n_eta-1; i++){
    TCanvas* c1 = new TCanvas();
    tdrCanvas(c1,"c1",h,4,10,kSquare,CorrectionObject::_lumitag);
    TLegend leg1 = tdrLeg(0.62,0.66,0.85,0.81);
    h->GetXaxis()->SetTitle("p_{T}^{ave} [GeV]");
    h->GetXaxis()->SetLimits(0,2000);
    h->GetYaxis()->SetTitle("entries");
    h->GetYaxis()->SetTitleOffset(1.5);
    double maximum = std::max(hdata_pt_ave[i]->GetMaximum(), hmc_pt_ave[i]->GetMaximum());
    h->GetYaxis()->SetRangeUser(0,1.2*maximum);
    hdata_pt_ave[i]->SetLineColor(kBlack);
    hdata_pt_ave[i]->Draw("SAME");
    hmc_pt_ave[i]->SetLineColor(kBlue);
    hmc_pt_ave[i]->Draw("HIST SAME");
    leg1.AddEntry(hdata_pt_ave[i], "DATA");
    leg1.AddEntry(hmc_pt_ave[i], "MC");
    leg1.Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    TString text = eta_range[i] + " < |#eta| < " + eta_range[i+1];
    tex->DrawLatex(0.52,0.85, text);

    c1->SaveAs(CorrectionObject::_outpath+"plots/control/Pt_ave_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + ".pdf");    


    delete tex;
    delete c1;
  }

















  delete c_0;
  delete h;
  
  for(int i=0; i<n_pt-1; i++){
    for(int j=0; j<n_eta-1; j++){
      delete hdata_asymmetry[i][j];
      delete hmc_asymmetry[i][j];
      delete hdata_B[i][j];
      delete hmc_B[i][j];
    }
  }

  for(int i=0; i<n_eta-1; i++){
    delete hdata_pt_ave[i];
    delete hmc_pt_ave[i];
  }


  delete f_mpf_mc;
  delete f_mpf_data;
  delete f_rel_mc;
  delete f_rel_data;








}
