#include "../include/parameters.h"
#include "../include/tdrstyle_mod15.h"
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
#include "THStack.h"



using namespace std;


//void FlavorCorrection_TTree(TString path="/nfs/dust/cms/user/karavdia/JERC/FlavorStudy/PFfractions_2016legacy/", TString _MCfile="/nfs/dust/cms/user/multh/JEC/2016Legacy/ClosureTest/Summer16_07Aug2017_V6/AK4CHS/MC_NoReweighted_NewHF_Binning/uhh2.AnalysisModuleRunner.MC.QCDPt50toInf_pythia8_AK4CHS_RunBCDEFGH.root", TString jettag="AK4PFchs", TString txttag="test", double al_cut=0.3){
void CorrectionObject::FlavorCorrection_TTree(){

  TString path="/nfs/dust/cms/user/karavdia/JERC/FlavorStudy/PFfractions_2016legacy/";
  //  TString _MCfile="/nfs/dust/cms/user/multh/JEC/2016Legacy/ClosureTest/Summer16_07Aug2017_V6/AK4CHS/MC_NoReweighted_NewHF_Binning/uhh2.AnalysisModuleRunner.MC.QCDPt50toInf_pythia8_AK4CHS_RunBCDEFGH.root";
  TString jettag=CorrectionObject::_jettag;
  TString txttag=CorrectionObject::_generator_tag; 
  double al_cut=alpha_cut;

  TString JetDescrib;                                                                                                                            
  if (jettag=="AK4PFchs") JetDescrib = "Anti-k_{t} R = 0.4, PF+CHS";                                                                             
  if (jettag=="AK4PFpuppi") JetDescrib = "Anti-k_{t} R = 0.4, PF+PUPPI";                                                                         
  if (jettag=="AK8PFchs") JetDescrib = "Anti-k_{t} R = 0.8, PF+CHS";                                                                             
  if (jettag=="AK8PFpuppi") JetDescrib = "Anti-k_{t} R = 0.8, PF+PUPPI"; 

  //set up needed histograms
  // TH1D* h_Flavor_Leadingjet = new TH1D("h_Flavor_Leadingjet", "Flavor (physics) of leading jet", 23, -1.5, 21.5);
  // TH1D* h_Flavor_Subleadingjet = new TH1D("h_Flavor_Subleadingjet", "Flavor (physics) of sub-leading jet", 23, -1.5, 21.5);
  // TH1D* h_Flavor_Barreljet = new TH1D("h_Flavor_Barreljet", "Flavor (physics) of barrel jet", 23, -1.5, 21.5);
  TH1D* h_Flavor_Probejet = new TH1D("h_Flavor_Probejet", "Flavor (physics) of probe jet", 23, -1.5, 21.5);

  // //log-equidistant bin-edges
  // double bin_edges[13];
  // double bins_start = 50;
  // double bins_end = 2100;
  // for(int i=0; i<13; i++){
  //   double it = i;
  //   bin_edges[i] = exp((1-it/12)*log(bins_start) + (it/12)*log(bins_end));
  //   cout << bin_edges[i] << endl;
  // }

  // TH1D* h_Gluon_pt = new TH1D("h_Gluon_pt", "p_{T} of gluon jets", 12, bin_edges);
  // TH1D* h_Quark_pt = new TH1D("h_Quark_pt", "p_{T} of quark jets", 12, bin_edges);
  // TH1D* h_Charm_pt = new TH1D("h_Charm_pt", "p_{T} of charm jets", 12, bin_edges);
  // TH1D* h_Bottom_pt = new TH1D("h_Bottom_pt", "p_{T} of bottom jets", 12, bin_edges);
  // TH1D* h_Unmatched_pt = new TH1D("h_Unmatched_pt", "p_{T} of unmatched jets", 12, bin_edges);

  //  TH1D* h_Jet_response = new TH1D("h_Jet_response", "Response of leading jet;Flavor (Physics);Leading Jet Response p_{T,jet}/p_{T,gen-jet}", 7, 0, 7);
  /*h_Jet_response->Fill("u",0.0001);
  h_Jet_response->Fill("d",0.0001);
  h_Jet_response->Fill("s",0.0001);
  h_Jet_response->Fill("c",0.0001);
  h_Jet_response->Fill("b",0.0001);
  h_Jet_response->Fill("g",0.0001);
  h_Jet_response->Fill("unm.",0.0001);*/

  // TH1D* h_u_response = new TH1D("h_u_response", "Response of leading jet (u) ;Leading Jet Response p_{T,jet}/p_{T,gen-jet};", 20, 0.5, 1.5);
  // TH1D* h_d_response = new TH1D("h_d_response", "Response of leading jet (d) ;Leading Jet Response p_{T,jet}/p_{T,gen-jet};", 20, 0.5, 1.5);
  // TH1D* h_s_response = new TH1D("h_s_response", "Response of leading jet (s) ;Leading Jet Response p_{T,jet}/p_{T,gen-jet};", 20, 0.5, 1.5);
  // TH1D* h_c_response = new TH1D("h_c_response", "Response of leading jet (c) ;Leading Jet Response p_{T,jet}/p_{T,gen-jet};", 20, 0.5, 1.5);
  // TH1D* h_b_response = new TH1D("h_b_response", "Response of leading jet (b) ;Leading Jet Response p_{T,jet}/p_{T,gen-jet};", 20, 0.5, 1.5);
  // TH1D* h_g_response = new TH1D("h_g_response", "Response of leading jet (g) ;Leading Jet Response p_{T,jet}/p_{T,gen-jet};", 20, 0.5, 1.5);
  // TH1D* h_unm_response = new TH1D("h_unm_response", "Response of leading jet (unm.) ;Leading Jet Response p_{T,jet}/p_{T,gen-jet};", 20, 0.5, 1.5);

 

  //read out needed values from TTree
  //  TFile* MCfile = new TFile(_MCfile,"READ");
  TTreeReader myReader_MC("AnalysisTree", CorrectionObject::_MCFile);

  TTreeReaderValue<Float_t> barreljet_eta(myReader_MC, "barreljet_eta");
  TTreeReaderValue<Float_t> probejet_eta(myReader_MC, "probejet_eta");
  TTreeReaderValue<Float_t> barreljet_pt(myReader_MC, "barreljet_pt");
  TTreeReaderValue<Float_t> probejet_pt(myReader_MC, "probejet_pt");
  //  TTreeReaderValue<Int_t> flavorBarreljet(myReader_MC, "flavorBarreljet");
  TTreeReaderValue<Int_t> flavorProbejet(myReader_MC, "flavorProbejet");
  //  TTreeReaderValue<Float_t> leadingjet_response(myReader_MC, "leadingjet_response");
  TTreeReaderValue<Float_t> had_n_Efrac(myReader_MC, "probejet_neutHadEF");
  TTreeReaderValue<Float_t> had_ch_Efrac(myReader_MC, "probejet_chHadEF");
  TTreeReaderValue<Float_t> em_n_Efrac(myReader_MC, "probejet_neutEmEF");
  TTreeReaderValue<Float_t> em_ch_Efrac(myReader_MC, "probejet_chEmEF");
  TTreeReaderValue<Float_t> mu_Efrac(myReader_MC, "probejet_muonEF");
  TTreeReaderValue<Float_t> ph_Efrac(myReader_MC, "probejet_photonEF");

  TTreeReaderValue<Float_t> alpha(myReader_MC, "alpha");
  TTreeReaderValue<Float_t> weight(myReader_MC, "weight");



  int n_leadingjet_u = 0, n_leadingjet_d = 0, n_leadingjet_s = 0, n_leadingjet_c = 0, n_leadingjet_b = 0, n_leadingjet_g = 0, n_leadingjet_unm = 0;
  // double energy_fractions_u[6][n_eta], energy_fractions_d[6][n_eta], energy_fractions_s[6][n_eta], energy_fractions_c[6][n_eta], energy_fractions_b[6][n_eta], energy_fractions_g[6][n_eta], energy_fractions_unm[6][n_eta];
  // for(int i=0; i<6; i++){
  //   for(int j=0; j<n_eta-1; j++){
  //   energy_fractions_u[i][j]=0;
  //   energy_fractions_d[i][j]=0;
  //   energy_fractions_s[i][j]=0;
  //   energy_fractions_c[i][j]=0;
  //   energy_fractions_b[i][j]=0;
  //   energy_fractions_g[i][j]=0;
  //   energy_fractions_unm[i][j]=0;
  //   }
  // }

  TH1D* energy_fractions_u[6][n_eta]; TH1D* energy_fractions_d[6][n_eta]; TH1D* energy_fractions_s[6][n_eta]; 
  TH1D*  energy_fractions_c[6][n_eta]; TH1D* energy_fractions_b[6][n_eta]; TH1D* energy_fractions_g[6][n_eta]; TH1D* energy_fractions_unm[6][n_eta];
 for(int i=0; i<6; i++){
    for(int j=0; j<n_eta-1; j++){
      TString add = "_"+i;      add+="_";      add+=j;
      energy_fractions_u[i][j] = new TH1D("energy_fractions_u_"+add,"",100,0,1);
      energy_fractions_d[i][j] = new TH1D("energy_fractions_d_"+add,"",100,0,1);
      energy_fractions_s[i][j] = new TH1D("energy_fractions_s_"+add,"",100,0,1);
      energy_fractions_c[i][j] = new TH1D("energy_fractions_c_"+add,"",100,0,1);
      energy_fractions_b[i][j] = new TH1D("energy_fractions_b_"+add,"",100,0,1);
      energy_fractions_g[i][j] = new TH1D("energy_fractions_g_"+add,"",100,0,1);
      energy_fractions_unm[i][j] = new TH1D("energy_fractions_unm_"+add,"",100,0,1);
    }
 }
 while (myReader_MC.Next()) {
    if(*alpha>al_cut) continue;
        
    double flavor_leadingjet = 0;
    double eta_leadingjet = -999;

    //fill flavor as cross-check, if jet is in barrel region
    
    //1: barreljet -> no eta check needed, but for safety do it nontheless
    // if(fabs(*barreljet_eta) < eta_bins[n_etabarr]){
    //   //  h_Flavor_Barreljet->Fill(*flavorBarreljet);
    //   /*if(*flavorBarreljet == -1) h_Unmatched_pt->Fill(*barreljet_pt);
    //   if(*flavorBarreljet >= 1 && *flavorBarreljet <= 3 ) h_Quark_pt->Fill(*barreljet_pt);
    //   if(*flavorBarreljet == 4) h_Charm_pt->Fill(*barreljet_pt);
    //   if(*flavorBarreljet == 5) h_Bottom_pt->Fill(*barreljet_pt);
    //   if(*flavorBarreljet == 21) h_Gluon_pt->Fill(*barreljet_pt);*/
    // } 
    // if(fabs(*probejet_eta) < eta_bins[n_etabarr]){
    //   h_Flavor_Probejet->Fill(*flavorProbejet);
    //   /*if(*flavorProbejet == -1) h_Unmatched_pt->Fill(*probejet_pt);
    //   if(*flavorProbejet >= 1 && *flavorProbejet <= 3 ) h_Quark_pt->Fill(*probejet_pt);
    //   if(*flavorProbejet == 4) h_Charm_pt->Fill(*probejet_pt);
    //   if(*flavorProbejet == 5) h_Bottom_pt->Fill(*probejet_pt);
    //   if(*flavorProbejet == 21) h_Gluon_pt->Fill(*probejet_pt);*/
    // }

    // //search for leading jet
    // if(*barreljet_pt > *probejet_pt) {
    //   flavor_leadingjet = *flavorBarreljet;
    //   eta_leadingjet = *barreljet_eta;
    //   if(fabs(*barreljet_eta) < eta_bins[n_etabarr]){
    // 	h_Flavor_Leadingjet->Fill(*flavorBarreljet);
    // 	if(*flavorBarreljet == -1) h_Unmatched_pt->Fill(*barreljet_pt);
    // 	if(*flavorBarreljet >= 1 && *flavorBarreljet <= 3 ) h_Quark_pt->Fill(*barreljet_pt);
    // 	if(*flavorBarreljet == 4) h_Charm_pt->Fill(*barreljet_pt);
    // 	if(*flavorBarreljet == 5) h_Bottom_pt->Fill(*barreljet_pt);
    // 	if(*flavorBarreljet == 21) h_Gluon_pt->Fill(*barreljet_pt);
    //   }
    //   if(fabs(*probejet_eta) < eta_bins[n_etabarr])  h_Flavor_Subleadingjet->Fill(*flavorProbejet);
    // }
    // else{
    // 	flavor_leadingjet = *flavorProbejet;
    // 	eta_leadingjet = *probejet_eta;
    //   if(fabs(*probejet_eta) < eta_bins[n_etabarr]){
    // 	h_Flavor_Leadingjet->Fill(*flavorProbejet);
    // 	if(*flavorProbejet == -1) h_Unmatched_pt->Fill(*probejet_pt);
    // 	if(*flavorProbejet >= 1 && *flavorProbejet <= 3 ) h_Quark_pt->Fill(*probejet_pt);
    // 	if(*flavorProbejet == 4) h_Charm_pt->Fill(*probejet_pt);
    // 	if(*flavorProbejet == 5) h_Bottom_pt->Fill(*probejet_pt);
    // 	if(*flavorProbejet == 21) h_Gluon_pt->Fill(*probejet_pt);
    //   }      
    //   if(*barreljet_eta < eta_bins[n_etabarr]) h_Flavor_Subleadingjet->Fill(*flavorBarreljet);  
    // }



  //response plots
  //for each flavor, calculate mean response of the leading jet
    //ToDo: take Mean of PF fraction distribution and not some strange sum as bellow!!!
    //check for eta -> barrel
    for(int j=0; j<n_eta-1; j++){
      //      if((*probejet_pt)>50) continue;
      //      if((*probejet_pt)>150 || (*probejet_pt)<50) continue;
      if((*probejet_pt)<150) continue;
      if(fabs(*probejet_eta)>eta_bins[j+1] || fabs(*probejet_eta)<eta_bins[j]) continue;
      //      cout<<"probejet_pt = "<<(*probejet_pt)<<endl;
      if(*flavorProbejet == 1) {
	//	h_u_response->Fill(*leadingjet_response);
	//	n_leadingjet_u++;
	energy_fractions_u[0][j] ->Fill( *had_n_Efrac,*weight);
	energy_fractions_u[1][j] ->Fill(  *had_ch_Efrac,*weight);
	energy_fractions_u[2][j] ->Fill(  *mu_Efrac,*weight);
	energy_fractions_u[3][j] ->Fill(  *ph_Efrac,*weight);
	energy_fractions_u[4][j] ->Fill(  *em_n_Efrac,*weight);
	energy_fractions_u[5][j] ->Fill(  *em_ch_Efrac,*weight);
      }
      else if(*flavorProbejet == 2){
	//	h_d_response->Fill(*leadingjet_response);
	//	n_leadingjet_d++;
	energy_fractions_d[0][j] ->Fill(  *had_n_Efrac,*weight);
	energy_fractions_d[1][j] ->Fill(  *had_ch_Efrac,*weight);
	energy_fractions_d[2][j] ->Fill(  *mu_Efrac,*weight);
	energy_fractions_d[3][j] ->Fill(  *ph_Efrac,*weight);
	energy_fractions_d[4][j] ->Fill(  *em_n_Efrac,*weight);
	energy_fractions_d[5][j] ->Fill(  *em_ch_Efrac,*weight);

      }
      else if(*flavorProbejet == 3){
	//	h_s_response->Fill(*leadingjet_response);
	//	n_leadingjet_s++;
	energy_fractions_s[0][j] ->Fill(  *had_n_Efrac,*weight);
	energy_fractions_s[1][j] ->Fill(  *had_ch_Efrac,*weight);
	energy_fractions_s[2][j] ->Fill(  *mu_Efrac,*weight);
	energy_fractions_s[3][j] ->Fill(  *ph_Efrac,*weight);
	energy_fractions_s[4][j] ->Fill(  *em_n_Efrac,*weight);
	energy_fractions_s[5][j] ->Fill(  *em_ch_Efrac,*weight);

      }
      else if(*flavorProbejet == 4){
	//	h_c_response->Fill(*leadingjet_response);
	//	n_leadingjet_c++;
	energy_fractions_c[0][j] ->Fill(  *had_n_Efrac,*weight);
	energy_fractions_c[1][j] ->Fill(  *had_ch_Efrac,*weight);
	energy_fractions_c[2][j] ->Fill(  *mu_Efrac,*weight);
	energy_fractions_c[3][j] ->Fill(  *ph_Efrac,*weight);
	energy_fractions_c[4][j] ->Fill(  *em_n_Efrac,*weight);
	energy_fractions_c[5][j] ->Fill(  *em_ch_Efrac,*weight);

      }
      else if(*flavorProbejet == 5){
	//	h_b_response->Fill(*leadingjet_response);
	//	n_leadingjet_b++;
	energy_fractions_b[0][j] ->Fill(  *had_n_Efrac,*weight);
	energy_fractions_b[1][j] ->Fill(  *had_ch_Efrac,*weight);
	energy_fractions_b[2][j] ->Fill(  *mu_Efrac,*weight);
	energy_fractions_b[3][j] ->Fill(  *ph_Efrac,*weight);
	energy_fractions_b[4][j] ->Fill(  *em_n_Efrac,*weight);
	energy_fractions_b[5][j] ->Fill(  *em_ch_Efrac,*weight);

      }
      else if(*flavorProbejet == 21){
	//	h_g_response->Fill(*leadingjet_response);
	//	n_leadingjet_g++;
	energy_fractions_g[0][j] ->Fill(  *had_n_Efrac,*weight);
	energy_fractions_g[1][j] ->Fill(  *had_ch_Efrac,*weight);
	energy_fractions_g[2][j] ->Fill(  *mu_Efrac,*weight);
	energy_fractions_g[3][j] ->Fill(  *ph_Efrac,*weight);
	energy_fractions_g[4][j] ->Fill(  *em_n_Efrac,*weight);
	energy_fractions_g[5][j] ->Fill(  *em_ch_Efrac,*weight);

      }
      else if(*flavorProbejet == -1){
	//	h_unm_response->Fill(*leadingjet_response);
	//	n_leadingjet_unm++;
	energy_fractions_unm[0][j] ->Fill(  *had_n_Efrac,*weight);
	energy_fractions_unm[1][j] ->Fill(  *had_ch_Efrac,*weight);
	energy_fractions_unm[2][j] ->Fill(  *mu_Efrac,*weight);
	energy_fractions_unm[3][j] ->Fill(  *ph_Efrac,*weight);
	energy_fractions_unm[4][j] ->Fill(  *em_n_Efrac,*weight);
	energy_fractions_unm[5][j] ->Fill(  *em_ch_Efrac,*weight);

      }
    // //    if(fabs(eta_leadingjet) < eta_bins[n_etabarr]) {
    //   //    if(fabs(*probejet_eta)<1.3){
    //   if(*flavorProbejet == 1) {
    // 	//	h_u_response->Fill(*leadingjet_response);
    // 	//	n_leadingjet_u++;
    // 	energy_fractions_u[0][j] +=  *had_n_Efrac*(*weight);
    // 	energy_fractions_u[1][j] +=  *had_ch_Efrac*(*weight);
    // 	energy_fractions_u[2][j] +=  *mu_Efrac*(*weight);
    // 	energy_fractions_u[3][j] +=  *ph_Efrac*(*weight);
    // 	energy_fractions_u[4][j] +=  *em_n_Efrac*(*weight);
    // 	energy_fractions_u[5][j] +=  *em_ch_Efrac*(*weight);
    //   }
    //   else if(*flavorProbejet == 2){
    // 	//	h_d_response->Fill(*leadingjet_response);
    // 	//	n_leadingjet_d++;
    // 	energy_fractions_d[0][j] +=  *had_n_Efrac*(*weight);
    // 	energy_fractions_d[1][j] +=  *had_ch_Efrac*(*weight);
    // 	energy_fractions_d[2][j] +=  *mu_Efrac*(*weight);
    // 	energy_fractions_d[3][j] +=  *ph_Efrac*(*weight);
    // 	energy_fractions_d[4][j] +=  *em_n_Efrac*(*weight);
    // 	energy_fractions_d[5][j] +=  *em_ch_Efrac*(*weight);

    //   }
    //   else if(*flavorProbejet == 3){
    // 	//	h_s_response->Fill(*leadingjet_response);
    // 	//	n_leadingjet_s++;
    // 	energy_fractions_s[0][j] +=  *had_n_Efrac*(*weight);
    // 	energy_fractions_s[1][j] +=  *had_ch_Efrac*(*weight);
    // 	energy_fractions_s[2][j] +=  *mu_Efrac*(*weight);
    // 	energy_fractions_s[3][j] +=  *ph_Efrac*(*weight);
    // 	energy_fractions_s[4][j] +=  *em_n_Efrac*(*weight);
    // 	energy_fractions_s[5][j] +=  *em_ch_Efrac*(*weight);

    //   }
    //   else if(*flavorProbejet == 4){
    // 	//	h_c_response->Fill(*leadingjet_response);
    // 	//	n_leadingjet_c++;
    // 	energy_fractions_c[0][j] +=  *had_n_Efrac*(*weight);
    // 	energy_fractions_c[1][j] +=  *had_ch_Efrac*(*weight);
    // 	energy_fractions_c[2][j] +=  *mu_Efrac*(*weight);
    // 	energy_fractions_c[3][j] +=  *ph_Efrac*(*weight);
    // 	energy_fractions_c[4][j] +=  *em_n_Efrac*(*weight);
    // 	energy_fractions_c[5][j] +=  *em_ch_Efrac*(*weight);

    //   }
    //   else if(*flavorProbejet == 5){
    // 	//	h_b_response->Fill(*leadingjet_response);
    // 	//	n_leadingjet_b++;
    // 	energy_fractions_b[0][j] +=  *had_n_Efrac*(*weight);
    // 	energy_fractions_b[1][j] +=  *had_ch_Efrac*(*weight);
    // 	energy_fractions_b[2][j] +=  *mu_Efrac*(*weight);
    // 	energy_fractions_b[3][j] +=  *ph_Efrac*(*weight);
    // 	energy_fractions_b[4][j] +=  *em_n_Efrac*(*weight);
    // 	energy_fractions_b[5][j] +=  *em_ch_Efrac*(*weight);

    //   }
    //   else if(*flavorProbejet == 21){
    // 	//	h_g_response->Fill(*leadingjet_response);
    // 	//	n_leadingjet_g++;
    // 	energy_fractions_g[0][j] +=  *had_n_Efrac*(*weight);
    // 	energy_fractions_g[1][j] +=  *had_ch_Efrac*(*weight);
    // 	energy_fractions_g[2][j] +=  *mu_Efrac*(*weight);
    // 	energy_fractions_g[3][j] +=  *ph_Efrac*(*weight);
    // 	energy_fractions_g[4][j] +=  *em_n_Efrac*(*weight);
    // 	energy_fractions_g[5][j] +=  *em_ch_Efrac*(*weight);

    //   }
    //   else if(*flavorProbejet == -1){
    // 	//	h_unm_response->Fill(*leadingjet_response);
    // 	//	n_leadingjet_unm++;
    // 	energy_fractions_unm[0][j] +=  *had_n_Efrac*(*weight);
    // 	energy_fractions_unm[1][j] +=  *had_ch_Efrac*(*weight);
    // 	energy_fractions_unm[2][j] +=  *mu_Efrac*(*weight);
    // 	energy_fractions_unm[3][j] +=  *ph_Efrac*(*weight);
    // 	energy_fractions_unm[4][j] +=  *em_n_Efrac*(*weight);
    // 	energy_fractions_unm[5][j] +=  *em_ch_Efrac*(*weight);

    //   }
      else {
	cout << "jetflavor: " << *flavorProbejet << endl;
	throw runtime_error("flavor of jet not known.");
      }
    }


  
  } //event loop

  
  // //scale each bin of pt-dependent plots to 1

  // for(int i=1; i<h_Gluon_pt->GetNbinsX()+1; i++){
  //   int n_entries = h_Gluon_pt->GetBinContent(i) + h_Quark_pt->GetBinContent(i) + h_Charm_pt->GetBinContent(i) + h_Bottom_pt->GetBinContent(i) + h_Unmatched_pt->GetBinContent(i);
  //   h_Gluon_pt->SetBinContent(i, h_Gluon_pt->GetBinContent(i) / n_entries);
  //   h_Gluon_pt->SetBinError(i, h_Gluon_pt->GetBinError(i) / n_entries);

  //   h_Quark_pt->SetBinContent(i, h_Quark_pt->GetBinContent(i) / n_entries);
  //   h_Quark_pt->SetBinError(i, h_Quark_pt->GetBinError(i) / n_entries);

  //   h_Charm_pt->SetBinContent(i, h_Charm_pt->GetBinContent(i) / n_entries);
  //   h_Charm_pt->SetBinError(i, h_Charm_pt->GetBinError(i) / n_entries);

  //   h_Bottom_pt->SetBinContent(i, h_Bottom_pt->GetBinContent(i) / n_entries);
  //   h_Bottom_pt->SetBinError(i, h_Bottom_pt->GetBinError(i) / n_entries);

  //   h_Unmatched_pt->SetBinContent(i, h_Unmatched_pt->GetBinContent(i) / n_entries);
  //   h_Unmatched_pt->SetBinError(i, h_Unmatched_pt->GetBinError(i) / n_entries);

  //   //crosscheck
  //   double sum_bincontents = h_Gluon_pt->GetBinContent(i) + h_Quark_pt->GetBinContent(i) + h_Charm_pt->GetBinContent(i) + h_Bottom_pt->GetBinContent(i) + h_Unmatched_pt->GetBinContent(i);
  //   cout << "After scaling -- Sum of bin contents in bin " << i << ": " << sum_bincontents << endl;
  // }



  /*
  //scale Energy_fractions by 1/n_entries(bin) to obtain mean response 
  for(int i=0; i<4; i++) { //0: neutral had, 1: charged had, 2: muon, 3: photon
    energy_fractions_u[i] = energy_fractions_u[i]/n_leadingjet_u;
    energy_fractions_d[i] = energy_fractions_d[i]/n_leadingjet_d;
    energy_fractions_s[i] = energy_fractions_s[i]/n_leadingjet_s;
    energy_fractions_c[i] = energy_fractions_c[i]/n_leadingjet_c;
    energy_fractions_b[i] = energy_fractions_b[i]/n_leadingjet_b;
    energy_fractions_g[i] = energy_fractions_g[i]/n_leadingjet_g;
    energy_fractions_unm[i] = energy_fractions_unm[i]/n_leadingjet_unm;

    if(i==0){ 
      h_neutral_hadron_Efraction->Fill("u",energy_fractions_u[i]);
      h_neutral_hadron_Efraction->Fill("d",energy_fractions_d[i]);
      h_neutral_hadron_Efraction->Fill("s",energy_fractions_s[i]);
      h_neutral_hadron_Efraction->Fill("c",energy_fractions_c[i]);
      h_neutral_hadron_Efraction->Fill("b",energy_fractions_b[i]);
      h_neutral_hadron_Efraction->Fill("g",energy_fractions_g[i]);
      h_neutral_hadron_Efraction->Fill("unm",energy_fractions_unm[i]); 
    }
    else if(i==1){ 
      h_charged_hadron_Efraction->Fill("u",energy_fractions_u[i]);
      h_charged_hadron_Efraction->Fill("d",energy_fractions_d[i]);
      h_charged_hadron_Efraction->Fill("s",energy_fractions_s[i]);
      h_charged_hadron_Efraction->Fill("c",energy_fractions_c[i]);
      h_charged_hadron_Efraction->Fill("b",energy_fractions_b[i]);
      h_charged_hadron_Efraction->Fill("g",energy_fractions_g[i]);
      h_charged_hadron_Efraction->Fill("unm",energy_fractions_unm[i]); 
    }
    else if(i==2){ 
      h_muon_Efraction->Fill("u",energy_fractions_u[i]);
      h_muon_Efraction->Fill("d",energy_fractions_d[i]);
      h_muon_Efraction->Fill("s",energy_fractions_s[i]);
      h_muon_Efraction->Fill("c",energy_fractions_c[i]);
      h_muon_Efraction->Fill("b",energy_fractions_b[i]);
      h_muon_Efraction->Fill("g",energy_fractions_g[i]);
      h_muon_Efraction->Fill("unm",energy_fractions_unm[i]); 
    }
    else if(i==3){ 
      h_photon_Efraction->Fill("u",energy_fractions_u[i]);
      h_photon_Efraction->Fill("d",energy_fractions_d[i]);
      h_photon_Efraction->Fill("s",energy_fractions_s[i]);
      h_photon_Efraction->Fill("c",energy_fractions_c[i]);
      h_photon_Efraction->Fill("b",energy_fractions_b[i]);
      h_photon_Efraction->Fill("g",energy_fractions_g[i]);
      h_photon_Efraction->Fill("unm",energy_fractions_unm[i]); 
    }
    }
*/



  //scale Energy_fractions to stack to unity 
for(int j=0; j<n_eta-1; j++){
  TH1D* h_neutral_hadron_Efraction = new TH1D("h_neutral_hadron_Efraction","Neutral Hadron energy fraction;Flavor (Physics);Probe jet PF energy fraction",7,0,7);
  TH1D* h_charged_hadron_Efraction = new TH1D("h_charged_hadron_Efraction","Charged Hadron energy fraction;Flavor (Physics);Probe jet PF energy fraction",7,0,7);
  TH1D* h_neutral_em_Efraction = new TH1D("h_neutral_em_Efraction","Neutral electromagnetic energy fraction;Flavor (Physics);Probe jet PF energy fraction",7,0,7);
  TH1D* h_charged_em_Efraction = new TH1D("h_charged_em_Efraction","Charged electromagnetic energy fraction;Flavor (Physics);Probe jet PF energy fraction",7,0,7);
  TH1D* h_photon_Efraction = new TH1D("h_photon_Efraction","Photon energy fraction;Flavor (Physics);Probe jet PF energy fraction",7,0,7);
  TH1D* h_muon_Efraction = new TH1D("h_muon_Efraction","Muon energy fraction;Flavor (Physics);Probe jet PF energy fraction",7,0,7);
  double sum_fractions_u = 0, sum_fractions_d = 0, sum_fractions_s = 0, sum_fractions_c = 0, sum_fractions_b = 0, sum_fractions_g = 0, sum_fractions_unm = 0;
  for(int i=0; i<6; i++) {
    sum_fractions_u += energy_fractions_u[i][j]->GetMean();
    sum_fractions_d += energy_fractions_d[i][j]->GetMean();
    sum_fractions_s += energy_fractions_s[i][j]->GetMean();
    sum_fractions_c += energy_fractions_c[i][j]->GetMean();
    sum_fractions_b += energy_fractions_b[i][j]->GetMean();
    sum_fractions_g += energy_fractions_g[i][j]->GetMean();
    sum_fractions_unm += energy_fractions_unm[i][j]->GetMean();
    // cout<<"energy_fractions_u[i][j]->GetMean() = "<<energy_fractions_u[i][j]->GetMean()<<" iPFfrac = "<<i<<" iEta = "<<j<<endl;
    // cout<<"energy_fractions_d[i][j]->GetMean() = "<<energy_fractions_d[i][j]->GetMean()<<" iPFfrac = "<<i<<" iEta = "<<j<<endl;
    // cout<<"energy_fractions_s[i][j]->GetMean() = "<<energy_fractions_s[i][j]->GetMean()<<" iPFfrac = "<<i<<" iEta = "<<j<<endl;
    // cout<<"energy_fractions_c[i][j]->GetMean() = "<<energy_fractions_c[i][j]->GetMean()<<" iPFfrac = "<<i<<" iEta = "<<j<<endl;
    // cout<<"energy_fractions_b[i][j]->GetMean() = "<<energy_fractions_b[i][j]->GetMean()<<" iPFfrac = "<<i<<" iEta = "<<j<<endl;
    // cout<<"energy_fractions_g[i][j]->GetMean() = "<<energy_fractions_g[i][j]->GetMean()<<" iPFfrac = "<<i<<" iEta = "<<j<<endl;
    // cout<<"energy_fractions_unm[i][j]->GetMean() = "<<energy_fractions_unm[i][j]->GetMean()<<" iPFfrac = "<<i<<" iEta = "<<j<<endl;
  }


  for(int i=0; i<6; i++) { //0: neutral had, 1: charged had, 2: muon, 3: photon, 4: neutral em, 5: charged em
    // energy_fractions_u[i][j] = energy_fractions_u[i][j]/sum_fractions_u;
    // energy_fractions_d[i][j] = energy_fractions_d[i][j]/sum_fractions_d;
    // energy_fractions_s[i][j] = energy_fractions_s[i][j]/sum_fractions_s;
    // energy_fractions_c[i][j] = energy_fractions_c[i][j]/sum_fractions_c;
    // energy_fractions_b[i][j] = energy_fractions_b[i][j]/sum_fractions_b;
    // energy_fractions_g[i][j] = energy_fractions_g[i][j]/sum_fractions_g;
    // energy_fractions_unm[i][j] = energy_fractions_unm[i][j]/sum_fractions_unm;
    // cout<<"sum_fractions_u = "<<sum_fractions_u<<endl;
    // cout<<"sum_fractions_d = "<<sum_fractions_d<<endl;
    // cout<<"sum_fractions_s = "<<sum_fractions_s<<endl;
    // cout<<"sum_fractions_c = "<<sum_fractions_c<<endl;
    // cout<<"sum_fractions_b = "<<sum_fractions_b<<endl;
    // cout<<"sum_fractions_g = "<<sum_fractions_g<<endl;
    // cout<<"sum_fractions_unm = "<<sum_fractions_unm<<endl;

    if(i==0){ 
      h_neutral_hadron_Efraction->Fill("u",(energy_fractions_u[i][j]->GetMean())/sum_fractions_u);
      h_neutral_hadron_Efraction->Fill("d",(energy_fractions_d[i][j]->GetMean())/sum_fractions_d);
      h_neutral_hadron_Efraction->Fill("g",(energy_fractions_g[i][j]->GetMean())/sum_fractions_g);
      h_neutral_hadron_Efraction->Fill("c",(energy_fractions_c[i][j]->GetMean())/sum_fractions_c);
      h_neutral_hadron_Efraction->Fill("b",(energy_fractions_b[i][j]->GetMean())/sum_fractions_b);
      h_neutral_hadron_Efraction->Fill("s",(energy_fractions_s[i][j]->GetMean())/sum_fractions_s);
      h_neutral_hadron_Efraction->Fill("unmatch",(energy_fractions_unm[i][j]->GetMean())/sum_fractions_unm); 
      // cout<<"(energy_fractions_u[i][j]->GetMean())/sum_fractions_u = "<<(energy_fractions_u[i][j]->GetMean())/sum_fractions_u<<endl;
      // cout<<"(energy_fractions_d[i][j]->GetMean())/sum_fractions_d = "<<(energy_fractions_d[i][j]->GetMean())/sum_fractions_d<<endl;
      // cout<<"(energy_fractions_s[i][j]->GetMean())/sum_fractions_s = "<<(energy_fractions_s[i][j]->GetMean())/sum_fractions_s<<endl;
      // cout<<"(energy_fractions_c[i][j]->GetMean())/sum_fractions_c = "<<(energy_fractions_c[i][j]->GetMean())/sum_fractions_c<<endl;
      // cout<<"(energy_fractions_b[i][j]->GetMean())/sum_fractions_b = "<<(energy_fractions_b[i][j]->GetMean())/sum_fractions_b<<endl;
      // cout<<"(energy_fractions_g[i][j]->GetMean())/sum_fractions_g = "<<(energy_fractions_g[i][j]->GetMean())/sum_fractions_g<<endl;
    }
    else if(i==1){ 
      h_charged_hadron_Efraction->Fill("u",(energy_fractions_u[i][j]->GetMean())/sum_fractions_u);
      h_charged_hadron_Efraction->Fill("d",(energy_fractions_d[i][j]->GetMean())/sum_fractions_d);
      h_charged_hadron_Efraction->Fill("g",(energy_fractions_g[i][j]->GetMean())/sum_fractions_g);
      h_charged_hadron_Efraction->Fill("c",(energy_fractions_c[i][j]->GetMean())/sum_fractions_c);
      h_charged_hadron_Efraction->Fill("b",(energy_fractions_b[i][j]->GetMean())/sum_fractions_b);
      h_charged_hadron_Efraction->Fill("s",(energy_fractions_s[i][j]->GetMean())/sum_fractions_s);
      h_charged_hadron_Efraction->Fill("unmatch",(energy_fractions_unm[i][j]->GetMean())/sum_fractions_unm); 
      // cout<<"(energy_fractions_u[i][j]->GetMean())/sum_fractions_u = "<<(energy_fractions_u[i][j]->GetMean())/sum_fractions_u<<endl;
      // cout<<"(energy_fractions_d[i][j]->GetMean())/sum_fractions_d = "<<(energy_fractions_d[i][j]->GetMean())/sum_fractions_d<<endl;
      // cout<<"(energy_fractions_s[i][j]->GetMean())/sum_fractions_s = "<<(energy_fractions_s[i][j]->GetMean())/sum_fractions_s<<endl;
      // cout<<"(energy_fractions_c[i][j]->GetMean())/sum_fractions_c = "<<(energy_fractions_c[i][j]->GetMean())/sum_fractions_c<<endl;
      // cout<<"(energy_fractions_b[i][j]->GetMean())/sum_fractions_b = "<<(energy_fractions_b[i][j]->GetMean())/sum_fractions_b<<endl;
      // cout<<"(energy_fractions_g[i][j]->GetMean())/sum_fractions_g = "<<(energy_fractions_g[i][j]->GetMean())/sum_fractions_g<<endl;
    }
    else if(i==2){ 
      h_muon_Efraction->Fill("u",(energy_fractions_u[i][j]->GetMean())/sum_fractions_u);
      h_muon_Efraction->Fill("d",(energy_fractions_d[i][j]->GetMean())/sum_fractions_d);
      h_muon_Efraction->Fill("g",(energy_fractions_g[i][j]->GetMean())/sum_fractions_g);
      h_muon_Efraction->Fill("c",(energy_fractions_c[i][j]->GetMean())/sum_fractions_c);
      h_muon_Efraction->Fill("b",(energy_fractions_b[i][j]->GetMean())/sum_fractions_b);
      h_muon_Efraction->Fill("s",(energy_fractions_s[i][j]->GetMean())/sum_fractions_s);
      h_muon_Efraction->Fill("unmatch",(energy_fractions_unm[i][j]->GetMean())/sum_fractions_unm); 
      // cout<<"(energy_fractions_u[i][j]->GetMean())/sum_fractions_u = "<<(energy_fractions_u[i][j]->GetMean())/sum_fractions_u<<endl;
      // cout<<"(energy_fractions_d[i][j]->GetMean())/sum_fractions_d = "<<(energy_fractions_d[i][j]->GetMean())/sum_fractions_d<<endl;
      // cout<<"(energy_fractions_s[i][j]->GetMean())/sum_fractions_s = "<<(energy_fractions_s[i][j]->GetMean())/sum_fractions_s<<endl;
      // cout<<"(energy_fractions_c[i][j]->GetMean())/sum_fractions_c = "<<(energy_fractions_c[i][j]->GetMean())/sum_fractions_c<<endl;
      // cout<<"(energy_fractions_b[i][j]->GetMean())/sum_fractions_b = "<<(energy_fractions_b[i][j]->GetMean())/sum_fractions_b<<endl;
      // cout<<"(energy_fractions_g[i][j]->GetMean())/sum_fractions_g = "<<(energy_fractions_g[i][j]->GetMean())/sum_fractions_g<<endl;
    }
    else if(i==3){ 
      h_photon_Efraction->Fill("u",(energy_fractions_u[i][j]->GetMean())/sum_fractions_u);
      h_photon_Efraction->Fill("d",(energy_fractions_d[i][j]->GetMean())/sum_fractions_d);
      h_photon_Efraction->Fill("g",(energy_fractions_g[i][j]->GetMean())/sum_fractions_g);
      h_photon_Efraction->Fill("c",(energy_fractions_c[i][j]->GetMean())/sum_fractions_c);
      h_photon_Efraction->Fill("b",(energy_fractions_b[i][j]->GetMean())/sum_fractions_b);
      h_photon_Efraction->Fill("s",(energy_fractions_s[i][j]->GetMean())/sum_fractions_s);
      h_photon_Efraction->Fill("unmatch",(energy_fractions_unm[i][j]->GetMean())/sum_fractions_unm); 
      // cout<<"Photon (energy_fractions_u[i][j]->GetMean())/sum_fractions_u = "<<(energy_fractions_u[i][j]->GetMean())/sum_fractions_u<<endl;
      // cout<<"Photon (energy_fractions_d[i][j]->GetMean())/sum_fractions_d = "<<(energy_fractions_d[i][j]->GetMean())/sum_fractions_d<<endl;
      // cout<<"Photon (energy_fractions_s[i][j]->GetMean())/sum_fractions_s = "<<(energy_fractions_s[i][j]->GetMean())/sum_fractions_s<<endl;
      // cout<<"Photon (energy_fractions_c[i][j]->GetMean())/sum_fractions_c = "<<(energy_fractions_c[i][j]->GetMean())/sum_fractions_c<<endl;
      // cout<<"Photon (energy_fractions_b[i][j]->GetMean())/sum_fractions_b = "<<(energy_fractions_b[i][j]->GetMean())/sum_fractions_b<<endl;
      // cout<<"Photon (energy_fractions_g[i][j]->GetMean())/sum_fractions_g = "<<(energy_fractions_g[i][j]->GetMean())/sum_fractions_g<<endl;
    }
    else if(i==4){ 
      h_neutral_em_Efraction->Fill("u",(energy_fractions_u[i][j]->GetMean())/sum_fractions_u);
      h_neutral_em_Efraction->Fill("d",(energy_fractions_d[i][j]->GetMean())/sum_fractions_d);
      h_neutral_em_Efraction->Fill("g",(energy_fractions_g[i][j]->GetMean())/sum_fractions_g);
      h_neutral_em_Efraction->Fill("c",(energy_fractions_c[i][j]->GetMean())/sum_fractions_c);
      h_neutral_em_Efraction->Fill("b",(energy_fractions_b[i][j]->GetMean())/sum_fractions_b);
      h_neutral_em_Efraction->Fill("s",(energy_fractions_s[i][j]->GetMean())/sum_fractions_s);
      h_neutral_em_Efraction->Fill("unmatch",(energy_fractions_unm[i][j]->GetMean())/sum_fractions_unm); 
      // cout<<"(energy_fractions_u[i][j]->GetMean())/sum_fractions_u = "<<(energy_fractions_u[i][j]->GetMean())/sum_fractions_u<<endl;
      // cout<<"(energy_fractions_d[i][j]->GetMean())/sum_fractions_d = "<<(energy_fractions_d[i][j]->GetMean())/sum_fractions_d<<endl;
      // cout<<"(energy_fractions_s[i][j]->GetMean())/sum_fractions_s = "<<(energy_fractions_s[i][j]->GetMean())/sum_fractions_s<<endl;
      // cout<<"(energy_fractions_c[i][j]->GetMean())/sum_fractions_c = "<<(energy_fractions_c[i][j]->GetMean())/sum_fractions_c<<endl;
      // cout<<"(energy_fractions_b[i][j]->GetMean())/sum_fractions_b = "<<(energy_fractions_b[i][j]->GetMean())/sum_fractions_b<<endl;
      // cout<<"(energy_fractions_g[i][j]->GetMean())/sum_fractions_g = "<<(energy_fractions_g[i][j]->GetMean())/sum_fractions_g<<endl;
    }

    else if(i==5){ 
      h_charged_em_Efraction->Fill("u",(energy_fractions_u[i][j]->GetMean())/sum_fractions_u);
      h_charged_em_Efraction->Fill("d",(energy_fractions_d[i][j]->GetMean())/sum_fractions_d);
      h_charged_em_Efraction->Fill("g",(energy_fractions_g[i][j]->GetMean())/sum_fractions_g);
      h_charged_em_Efraction->Fill("c",(energy_fractions_c[i][j]->GetMean())/sum_fractions_c);
      h_charged_em_Efraction->Fill("b",(energy_fractions_b[i][j]->GetMean())/sum_fractions_b);
      h_charged_em_Efraction->Fill("s",(energy_fractions_s[i][j]->GetMean())/sum_fractions_s);
      h_charged_em_Efraction->Fill("unmatch",(energy_fractions_unm[i][j]->GetMean())/sum_fractions_unm); 
      // cout<<"(energy_fractions_u[i][j]->GetMean())/sum_fractions_u = "<<(energy_fractions_u[i][j]->GetMean())/sum_fractions_u<<endl;
      // cout<<"(energy_fractions_d[i][j]->GetMean())/sum_fractions_d = "<<(energy_fractions_d[i][j]->GetMean())/sum_fractions_d<<endl;
      // cout<<"(energy_fractions_s[i][j]->GetMean())/sum_fractions_s = "<<(energy_fractions_s[i][j]->GetMean())/sum_fractions_s<<endl;
      // cout<<"(energy_fractions_c[i][j]->GetMean())/sum_fractions_c = "<<(energy_fractions_c[i][j]->GetMean())/sum_fractions_c<<endl;
      // cout<<"(energy_fractions_b[i][j]->GetMean())/sum_fractions_b = "<<(energy_fractions_b[i][j]->GetMean())/sum_fractions_b<<endl;
      // cout<<"(energy_fractions_g[i][j]->GetMean())/sum_fractions_g = "<<(energy_fractions_g[i][j]->GetMean())/sum_fractions_g<<endl;
    }
  }




  
  // TCanvas* cfit1 = new TCanvas("cfit1", "cfit1",1);
  // h_u_response->Fit("gaus","","",0.8,1.2);
  // h_Jet_response->Fill("u", h_u_response->GetFunction("gaus")->GetParameter(1));
  // h_Jet_response->SetBinError(1, h_u_response->GetFunction("gaus")->GetParameter(2));
  // h_u_response->Draw("E1 SAME");

  // TCanvas* cfit2 = new TCanvas("cfit2", "cfit2",1);
  // h_d_response->Fit("gaus","","",0.8,1.2);
  // h_Jet_response->Fill("d", h_d_response->GetFunction("gaus")->GetParameter(1));
  // h_Jet_response->SetBinError(2, h_d_response->GetFunction("gaus")->GetParameter(2));
  // h_d_response->Draw("E1 SAME");

  // TCanvas* cfit3 = new TCanvas("cfit3", "cfit3",1);
  // h_s_response->Fit("gaus","","",0.8,1.2);
  // h_Jet_response->Fill("s", h_s_response->GetFunction("gaus")->GetParameter(1));
  // h_Jet_response->SetBinError(3, h_s_response->GetFunction("gaus")->GetParameter(2));
  // h_s_response->Draw("E1 SAME");

  // TCanvas* cfit4 = new TCanvas("cfit4", "cfit4",1);
  // h_c_response->Fit("gaus","","",0.8,1.2);
  // h_Jet_response->Fill("c", h_c_response->GetFunction("gaus")->GetParameter(1));
  // h_Jet_response->SetBinError(4, h_c_response->GetFunction("gaus")->GetParameter(2));
  // h_c_response->Draw("E1 SAME");

  // TCanvas* cfit5 = new TCanvas("cfit5", "cfit5",1);
  // h_b_response->Fit("gaus","","",0.8,1.2);
  // h_Jet_response->Fill("b", h_b_response->GetFunction("gaus")->GetParameter(1));
  // h_Jet_response->SetBinError(5, h_b_response->GetFunction("gaus")->GetParameter(2));
  // h_b_response->Draw("E1 SAME");

  // TCanvas* cfit6 = new TCanvas("cfit6", "cfit6",1);
  // h_g_response->Fit("gaus","","",0.8,1.2);
  // h_Jet_response->Fill("g", h_g_response->GetFunction("gaus")->GetParameter(1));
  // h_Jet_response->SetBinError(6, h_g_response->GetFunction("gaus")->GetParameter(2));
  // h_g_response->Draw("E1 SAME");

  // TCanvas* cfit7 = new TCanvas("cfit7", "cfit7",1);
  // h_unm_response->Fit("gaus","","",0.8,1.2);
  // h_Jet_response->Fill("unm.", h_unm_response->GetFunction("gaus")->GetParameter(1));
  // h_Jet_response->SetBinError(7, h_unm_response->GetFunction("gaus")->GetParameter(2));
  // h_unm_response->Draw("E1 SAME");

  
  //Flavor plots
  // TCanvas* c1 = new TCanvas("c1", "c1", 1);
  // h_Flavor_Leadingjet->Draw();
  // TCanvas* c2 = new TCanvas("c2", "c2", 1);
  // h_Flavor_Subleadingjet->Draw();
  // TCanvas* c3 = new TCanvas("c3", "c3", 1);
  // h_Flavor_Barreljet->Draw();
  TCanvas* c4 = new TCanvas("c4", "c4", 1);
  h_Flavor_Probejet->Draw();

  // //flavor pt-dependent
  // //cosmetics
  // h_Quark_pt->SetMarkerColor(6);
  // h_Quark_pt->SetLineColor(6);
  // h_Quark_pt->SetMarkerStyle(20);
  // h_Gluon_pt->SetMarkerColor(4);
  // h_Gluon_pt->SetLineColor(4);
  // h_Gluon_pt->SetMarkerStyle(20);
  // h_Charm_pt->SetMarkerColor(8);
  // h_Charm_pt->SetLineColor(8);
  // h_Charm_pt->SetMarkerStyle(21);
  // h_Bottom_pt->SetMarkerColor(2);
  // h_Bottom_pt->SetLineColor(2);
  // h_Bottom_pt->SetMarkerStyle(21);
  // h_Unmatched_pt->SetMarkerColor(15);
  // h_Unmatched_pt->SetLineColor(15);
  // h_Unmatched_pt->SetMarkerStyle(33);

  // //Draw
  // TCanvas* c5 = new TCanvas("c5", "c5", 1);
  // h_Gluon_pt->Draw("E1");
  // TCanvas* c6 = new TCanvas("c6", "c6", 1);
  // h_Quark_pt->Draw("E1");
  // TCanvas* c7 = new TCanvas("c7", "c7", 1);
  // h_Charm_pt->Draw("E1");
  // TCanvas* c8 = new TCanvas("c8", "c8", 1);
  // h_Bottom_pt->Draw("E1");
  // TCanvas* c9 = new TCanvas("c9", "c9", 1);
  // h_Unmatched_pt->Draw("E1");




  // //main plot
  TString alVal;
  alVal.Form("%0.2f\n",al_cut);
  TString altitle = " #alpha<"+alVal;

  TString etaVal;
  etaVal.Form("%0.2f\n",eta_bins[n_etabarr]);
 

  // TH1D* h_Unmatched_ForDraw = (TH1D*)h_Unmatched_pt->Clone();
  // h_Unmatched_ForDraw->SetTitle(";p_{T,jet 1} (GeV);Flavor fraction");
  // h_Unmatched_ForDraw->SetMinimum(0.);
  // h_Unmatched_ForDraw->SetMaximum(1.15);
  // h_Unmatched_ForDraw->GetYaxis()->SetTitleSize(0.05);
  // h_Unmatched_ForDraw->GetXaxis()->SetTitleSize(0.05);



  // //lumi_13TeV = "Run2016  12.9 fb^{-1}";
  // bool kSquare = true;

  // TLegend leg1 = tdrLeg(0.43,0.55,0.75,0.8);
  // leg1.AddEntry(h_Quark_pt, "Quark","P");
  // leg1.AddEntry(h_Gluon_pt, "Gluon","P");
  // leg1.AddEntry(h_Charm_pt, "Charm","P");
  // leg1.AddEntry(h_Bottom_pt, "Bottom","P");
  // leg1.AddEntry(h_Unmatched_ForDraw, "Unmatched","P");
  // leg1.Draw();

  TString Latextext1 = JetDescrib;
  //  TString Latextext2 = "QCD dijet, |#eta|<" + etaVal  + altitle;
  TString Latextext2 = "QCD dijet, " + eta_range[j] + " < |#eta| < " + eta_range[j+1]  + altitle;
 // TString text_eta = eta_range[j] + " < |#eta| < " + eta_range[j+1];
 //  tex->DrawLatex(0.52,0.85, text_eta);

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);  
 

  // TCanvas* c10 = new TCanvas();
  // tdrCanvas(c10,"c10", h_Unmatched_ForDraw, 4, 10, kSquare);
  // gPad->SetLogx();



  // h_Unmatched_ForDraw->Draw("E1 SAME");
  // h_Bottom_pt->Draw("E1 SAME");
  // h_Charm_pt->Draw("E1 SAME");
  // h_Quark_pt->Draw("E1 SAME");
  // h_Gluon_pt->Draw("E1 SAME");
  // leg1.Draw();
  // tex->DrawLatex(0.43,0.87,Latextext1);
  // tex->DrawLatex(0.43,0.83,Latextext2);
 
  // c10->SaveAs(path+"plots/FlavorFractions_"+jettag+"_"+txttag+".pdf");

  // //Response Plot
  // TCanvas* c11 = new TCanvas();
  // tdrCanvas(c11,"c11", h_Jet_response, 4, 10, kSquare);
  // h_Jet_response->Draw("E1 SAME");
  // tex->DrawLatex(0.43,0.87,Latextext1);
  // tex->DrawLatex(0.43,0.83,Latextext2);
  
  
  //TCanvas* c_stack =  tdrCanvas("c_stack", h_neutral_hadron_Efraction, 4, 10, kSquare);
  TCanvas* c_stack =  new TCanvas("c_stack","c_stack",1);
 // draw stack
  THStack* stack = new THStack("stack",";Flavor;Probe jet energy fraction");
  h_neutral_hadron_Efraction->SetLineColor(8);
  h_neutral_hadron_Efraction->SetFillColor(8);
  h_charged_hadron_Efraction->SetLineColor(7);
  h_charged_hadron_Efraction->SetFillColor(7);
  h_muon_Efraction->SetLineColor(46);
  h_muon_Efraction->SetFillColor(46);
  h_photon_Efraction->SetLineColor(41);
  h_photon_Efraction->SetFillColor(41);
  h_neutral_em_Efraction->SetLineColor(15);
  h_neutral_em_Efraction->SetFillColor(15);
  h_charged_em_Efraction->SetLineColor(5);
  h_charged_em_Efraction->SetFillColor(5);

  stack->Add(h_neutral_hadron_Efraction);
  stack->Add(h_charged_hadron_Efraction);
  stack->Add(h_neutral_em_Efraction);
  stack->Add(h_charged_em_Efraction);
  stack->Add(h_photon_Efraction);
  stack->Add(h_muon_Efraction);

  stack->Draw("HIST");

  TLegend* leg2 =  new TLegend(0.70,0.78,0.99,0.99);
  leg2 -> AddEntry(h_muon_Efraction, "#mu");
  leg2 -> AddEntry(h_photon_Efraction, "#gamma");
  leg2 -> AddEntry(h_charged_hadron_Efraction, "Ch. hadron");
  leg2 -> AddEntry(h_neutral_hadron_Efraction, "Neut. hadron");
  leg2 -> AddEntry(h_charged_em_Efraction, "Ch. EM");
  leg2 -> AddEntry(h_neutral_em_Efraction, "Neut. EM");
  leg2 -> SetFillColor(10);
  leg2 -> SetBorderSize(0);
  leg2 -> SetTextSize(0.042);
  leg2 -> SetLineColor(1);
  leg2 -> SetTextFont(42);
  leg2 -> Draw();
  tex->DrawLatex(0.22,0.57,Latextext1);
  tex->DrawLatex(0.22,0.53,Latextext2);
 
  c_stack->SaveAs(path+"plots/PFEnergyFractions_"+jettag+"_"+txttag+ "_eta_" + eta_range2[j] + "_" + eta_range2[j+1]+".pdf");
  delete stack;
  delete c_stack;

 }

  cout << "Hello, I'm done" << endl;





}
