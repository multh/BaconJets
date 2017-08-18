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

void CorrectionObject::GenResponsePlots(){
  cout << "--------------- Starting GenResponsePlots() ---------------" << endl << endl;
  gStyle->SetOptStat(0);

  //Table with number of events in each pT- and eta-bin
  
  //Set up histos for ratios of responses
  // TH1D *hmc_response_probejet[n_pt-1][n_eta-1];   // <pTreco>/<pTgen> for tag&probe jet matched to GEN jets
  // TH1D *hmc_response_probetagjet[n_pt-1][n_eta-1];// <pT,probe,gen>/<pT,tag,gen> for tag&probe jet matched to GEN jets
  // TH1D *hmc_response_probetagreco[n_pt-1][n_eta-1];// <pT,probe,RECO>/<pT,tag,RECO> for tag&probe jet matched to GEN jets
  // TH1D *hmc_response_probetagptcl[n_pt-1][n_eta-1];// <pT,probe,parton>/<pT,tag,parton> for tag&probe jet matched to PARTONS
  // TH1D *hmc_response_tagtagjet[n_pt-1][n_eta-1];// <pT,tag,RECO>/<pT,tag,gen> for tag&probe jet matched to GEN jets

  TH1D *hmc_A[n_pt-1][n_eta-1];   // Assymetry_RECO tag&probe jet matched to GEN jets

  TH1I *hmc_jet1_genID[n_pt-1][n_eta-1];// genID for the 1st jet
  TH1I *hmc_jet2_genID[n_pt-1][n_eta-1];// genID for the 1st jet
  TH1I *hmc_jet3_genID[n_pt-1][n_eta-1];// genID for the 1st jet

  TH1D *hmc_probejetpt_flavor[3][n_eta-1];// probe jet pt separated by flavor, 0 = not matched, 1 = quark, 2 = gluon
  TH1D *hmc_tagjetpt_flavor[3][n_eta-1];// tag jet pt separated by flavor, 0 = not matched, 1 = quark, 2 = gluon



  TH1D *hmc_normjetpt[n_pt-1][n_eta-1];//  binning variable used for normalisation
  TH1D *hmc_probejetpt[n_pt-1][n_eta-1];// RECO probe jet pt devided to binning variable, e.g <pT,probe,RECO> = <pT,probe,RECO/pT,ave,RECO> * <pT,ave,RECO>
  TH1D *hmc_probegenjetpt[n_pt-1][n_eta-1];// GEN probe jet pt devided to binning variable, e.g <pT,probe,GEN> = <pT,probe,GEN/pT,ave,RECO> * <pT,ave,RECO>
  TH1D *hmc_tagjetpt[n_pt-1][n_eta-1];// RECO tag jet pt devided to binning variable, e.g <pT,tag,RECO> = <pT,tag,RECO/pT,ave,RECO> * <pT,ave,RECO>
  TH1D *hmc_taggenjetpt[n_pt-1][n_eta-1];// GEN tag jet pt devided to binning variable, e.g <pT,tag,GEN> = <pT,tag,GEN/pT,ave,RECO> * <pT,ave,RECO>

  int count = 0;
 
  /*  TString name3 = "hist_mc_genresponse_probe_";
  TString name4 = "hist_mc_genresponse_probetag_";
  TString name10 = "hist_mc_genresponse_tagtag_";
  TString name5 = "hist_mc_genresponse_probetag_ptcl_";
  TString name2 = "hist_mc_genresponse_probetag_reco_";*/

  TString name11 = "hist_mc_A_";

  TString name6 = "hist_mc_jet1_genID";
  TString name7 = "hist_mc_jet2_genID";
  TString name8 = "hist_mc_jet3_genID";

  TString name9 = "hist_mc_probejetpt_flavor";
  TString name91 = "hist_mc_tagjetpt_flavor";

  TString name99 = "hist_mc_normpt_";
  TString name100 = "hist_mc_probejetrecopt_";
  TString name101 = "hist_mc_probejetgenpt_";
  TString name102 = "hist_mc_tagjetrecopt_";
  TString name103 = "hist_mc_tagjetgenpt_";
 
  for(int j=0; j<n_eta-1; j++){
      TString eta_name = "eta_"+eta_range2[j]+"_"+eta_range2[j+1];
    for(int k=0; k<n_pt-1; k++){
      TString pt_name = "pt_"+pt_range[k]+"_"+pt_range[k+1];
      TString name;
      /*      TString name = name3 + eta_name + "_" + pt_name;
      hmc_response_probejet[k][j] = new TH1D(name,"",nResponseBins, -2.0, 2.0);
      name = name4 + eta_name + "_" + pt_name;
      hmc_response_probetagjet[k][j] = new TH1D(name,"",nResponseBins, -2.0, 2.0);
      name = name10 + eta_name + "_" + pt_name;
      hmc_response_tagtagjet[k][j] = new TH1D(name,"",nResponseBins, -2.0, 2.0);
      name = name5 + eta_name + "_" + pt_name;
      hmc_response_probetagptcl[k][j] = new TH1D(name,"",nResponseBins, -2.0, 2.0);
      name = name2 + eta_name + "_" + pt_name;
      hmc_response_probetagreco[k][j] = new TH1D(name,"",nResponseBins, -2.0, 2.0);*/

      name = name11 + eta_name + "_" + pt_name;
      hmc_A[k][j] = new TH1D(name,"",nResponseBins, -2.0, 2.0);

      name = name99 + eta_name + "_" + pt_name;
      hmc_normjetpt[k][j] = new TH1D(name,"",250, 0, 1000);
      name = name100 + eta_name + "_" + pt_name;
      hmc_probejetpt[k][j] = new TH1D(name,"",nResponseBins, 0.0, 2.0);
      name = name101 + eta_name + "_" + pt_name;
      hmc_probegenjetpt[k][j] = new TH1D(name,"",nResponseBins, 0.0, 2.0);
      name = name102 + eta_name + "_" + pt_name;
      hmc_tagjetpt[k][j] = new TH1D(name,"",nResponseBins, 0.0, 2.0);
      name = name103 + eta_name + "_" + pt_name;
      hmc_taggenjetpt[k][j] = new TH1D(name,"",nResponseBins, 0.0, 2.0);

      name = name6 + eta_name + "_" + pt_name;
      hmc_jet1_genID[k][j] = new TH1I(name,"",15,-5,10);
      name = name7 + eta_name + "_" + pt_name;
      hmc_jet2_genID[k][j] = new TH1I(name,"",15,-5,10);
      name = name8 + eta_name + "_" + pt_name;
      hmc_jet3_genID[k][j] = new TH1I(name,"",15,-5,10);


      //      count++;
    }

    TString name = name9 + eta_name;
    for(int ifl=0;ifl<3;ifl++){
      if(ifl==0) name +="_notmatched";
      if(ifl==1) name +="_quark";
      if(ifl==2) name +="_gluon";
      hmc_probejetpt_flavor[ifl][j] = new TH1D(name,"",100,0,1000);
    }

    name = name91 + eta_name;
    for(int ifl=0;ifl<3;ifl++){
      if(ifl==0) name +="_notmatched";
      if(ifl==1) name +="_quark";
      if(ifl==2) name +="_gluon";
      hmc_tagjetpt_flavor[ifl][j] = new TH1D(name,"",100,0,1000);
    }

  }



  //Get relevant information from MC, loop over MC events 
  TTreeReader myReader_MC("AnalysisTree", CorrectionObject::_MCFile);
  TTreeReaderValue<Float_t> pt_ave_mc(myReader_MC, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_mc(myReader_MC, "probejet_eta");
  TTreeReaderValue<Float_t> probejet_pt_mc(myReader_MC, "probejet_pt");
  TTreeReaderValue<Float_t> barreljet_pt_mc(myReader_MC, "barreljet_pt");
  TTreeReaderValue<Float_t> alpha_mc(myReader_MC, "alpha");
  TTreeReaderValue<Float_t> weight_mc(myReader_MC, "weight");

  TTreeReaderValue<Float_t> probejet_ptgen_mc(myReader_MC, "probejet_ptgen");
  TTreeReaderValue<Float_t> barreljet_ptgen_mc(myReader_MC, "barreljet_ptgen");

  TTreeReaderValue<Float_t> probejet_ptptcl_mc(myReader_MC, "probejet_ptptcl");
  TTreeReaderValue<Float_t> barreljet_ptptcl_mc(myReader_MC, "barreljet_ptptcl");

  TTreeReaderValue<Int_t> jet1_genID_mc(myReader_MC, "jet1_genID");
  TTreeReaderValue<Int_t> jet2_genID_mc(myReader_MC, "jet2_genID");
  TTreeReaderValue<Int_t> jet3_genID_mc(myReader_MC, "jet3_genID");

  TTreeReaderValue<Int_t> flavorProbejet_mc(myReader_MC, "flavorProbejet");
  TTreeReaderValue<Int_t> flavorTagjet_mc(myReader_MC, "flavorBarreljet");

  int icount=0;

  TString pt_binning_var_str = "#bar{p}^{RECO}_{T} [GeV]";//bin in pt_ave, RECO
  TString pt_binning_var_name = "__pT_ave_RECO__";//bin in pt_ave, RECO

  // TString pt_binning_var_str = "p^{tag,GEN}_{T} [GeV]";//bin in pt_tag, GEN
  // TString pt_binning_var_name = "__pT_tag_GEN__";//bin in pt_tag, GEN

  // TString pt_binning_var_str = "p^{probe,GEN}_{T} [GeV]";//bin in pt_probe, GEN
  // TString pt_binning_var_name = "__pT_probe_GEN__";//bin in pt_probe, GEN

  while (myReader_MC.Next()) {
  //  while (myReader_MC.Next() && icount<1000) {
    double pt_binning_var = *pt_ave_mc;//bin in pt_ave, RECO
    //    double pt_binning_var = *barreljet_ptgen_mc;//bin in pt_tag, GEN
    //    double pt_binning_var = *probejet_ptgen_mc;//bin in pt_probe, GEN
  
    if(*alpha_mc>alpha_cut) continue;
    //fill histos in bins of pt and eta
    for(int k=0; k<n_pt-1; k++){
      if(pt_binning_var<pt_bins[k] || pt_binning_var>pt_bins[k+1]) continue;
      for(int j=0; j<n_eta-1; j++){
	if(fabs(*probejet_eta_mc)>eta_bins[j+1] || fabs(*probejet_eta_mc)<eta_bins[j]) continue;
	else{
	  //	  bool matched;
	  if(*probejet_ptgen_mc<0 || *barreljet_ptgen_mc<0){ //not matched
	    //	    matched = false;
	  }
	  else{ //matched

	    ///[BEGIN] Selection according to flavor of tag&probe jets---------------
	    //	    bool flavor_sel=false;
	    // if(*flavorProbejet_mc>0 && *flavorProbejet_mc<6 && *flavorTagjet_mc>0 && *flavorTagjet_mc<6) 
	    //   flavor_sel=true;//QQ
	    // if(*flavorTagjet_mc==21 && *flavorProbejet_mc>0 && *flavorProbejet_mc<6 ) 
	    //   flavor_sel=true;//GQ
	    // if(*flavorTagjet_mc==21 && *flavorProbejet_mc==21) 
	    //   flavor_sel=true;//GG

	    // if(*flavorTagjet_mc>0 && *flavorTagjet_mc<6 && *flavorProbejet_mc==21) 
	    //   flavor_sel=true;//QG
	    // if(!flavor_sel) continue;

	    ///[END] Selection according to flavor of tag&probe jets---------------

	    ///// flavor ////
	  if(*flavorProbejet_mc<0){//not matched
	    hmc_probejetpt_flavor[0][j]->Fill(*probejet_pt_mc,*weight_mc);
	  }
	  else{
	    if(*flavorProbejet_mc>0 && *flavorProbejet_mc<6){
	      hmc_probejetpt_flavor[1][j]->Fill(*probejet_pt_mc,*weight_mc);
	    }
	    else{
	      hmc_probejetpt_flavor[2][j]->Fill(*probejet_pt_mc,*weight_mc);
	    }
	  }
	  if(*flavorTagjet_mc<0){
	    hmc_tagjetpt_flavor[0][j]->Fill(*barreljet_pt_mc,*weight_mc);
	  }
	  else{
	    if(*flavorTagjet_mc>0 && *flavorTagjet_mc<6){
	      hmc_tagjetpt_flavor[1][j]->Fill(*barreljet_pt_mc,*weight_mc);
	    }
	    else{
	      hmc_tagjetpt_flavor[2][j]->Fill(*barreljet_pt_mc,*weight_mc);
	    }
	  }
	  ///// [END] flavor /////
	  hmc_normjetpt[k][j]->Fill(pt_binning_var,*weight_mc);

	  double probejetpt_norm = (*probejet_pt_mc)/pt_binning_var;
	  hmc_probejetpt[k][j]->Fill(probejetpt_norm,*weight_mc);
	  double probegenjetpt_norm = (*probejet_ptgen_mc)/pt_binning_var;
	  hmc_probegenjetpt[k][j]->Fill(probegenjetpt_norm,*weight_mc);

	  double tagjetpt_norm = (*barreljet_pt_mc)/pt_binning_var;
	  hmc_tagjetpt[k][j]->Fill(tagjetpt_norm,*weight_mc);
	  double taggenjetpt_norm = (*barreljet_ptgen_mc)/pt_binning_var;
	  hmc_taggenjetpt[k][j]->Fill(taggenjetpt_norm,*weight_mc);

	  /*	  double genresponse = (*probejet_pt_mc)/(*probejet_ptgen_mc);
	  hmc_response_probejet[k][j]->Fill(genresponse,*weight_mc);
	  double gentagproberesponse = (*probejet_ptgen_mc)/(*barreljet_ptgen_mc);
	  hmc_response_probetagjet[k][j]->Fill(gentagproberesponse,*weight_mc);
	  double gentagprobeptclresponse = (*probejet_ptptcl_mc)/(*barreljet_ptptcl_mc);
	  hmc_response_probetagptcl[k][j]->Fill(gentagprobeptclresponse,*weight_mc);
	  double recoresponse = (*probejet_pt_mc)/(*barreljet_pt_mc);
	  hmc_response_probetagreco[k][j]->Fill(recoresponse,*weight_mc);
	  double tagtagresponse = (*barreljet_pt_mc)/(*barreljet_ptgen_mc);
	  //	    cout<<"barreljet_pt = "<<*barreljet_pt_mc<<" barreljet_ptgen = "<<*barreljet_ptgen_mc<<" tagtagresponse = "<<tagtagresponse<<endl;
	  hmc_response_tagtagjet[k][j]->Fill(tagtagresponse,*weight_mc);*/
	  double assymetry = ((*probejet_pt_mc)-(*barreljet_pt_mc))/((*probejet_pt_mc)+(*barreljet_pt_mc));
	  hmc_A[k][j]->Fill(assymetry,*weight_mc);
	  //	    cout<<"probe_gen = "<<(*probejet_ptgen_mc)<<" barrel_gen = "<<(*barreljet_ptgen_mc)<<endl;
	  //	    cout<<"genresponse = "<<genresponse<<" gentagproberesponse = "<<gentagproberesponse<<endl;
	  }
	  int jet1_genID_mc_val = *jet1_genID_mc;
	  if(jet1_genID_mc_val>-1){
	    jet1_genID_mc_val++;
	    hmc_jet1_genID[k][j]->Fill(jet1_genID_mc_val-1,*weight_mc);
	  }
	  int jet2_genID_mc_val = *jet2_genID_mc;
	  if(jet2_genID_mc_val>-1){
	    jet2_genID_mc_val++;
	    hmc_jet2_genID[k][j]->Fill(jet2_genID_mc_val-2,*weight_mc);
	  }

	  int jet3_genID_mc_val = *jet3_genID_mc;
	  if(jet3_genID_mc_val>-1){
	    jet3_genID_mc_val++;
	    hmc_jet3_genID[k][j]->Fill(jet3_genID_mc_val-3,*weight_mc);
	  }
	  //	  cout<<""<<jet1_genID_mc_val<<" "<<jet2_genID_mc_val<<" "<<jet3_genID_mc_val<<endl;
	 
	}
      }
    }
    icount++;
  } 

  /*  ofstream output;
  output.open(CorrectionObject::_outpath+"plots/control/GenResponse_Number_Events_Pt_Eta_bins_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt");
    

  output << "Number of matched events in each bin for MC" << endl;
  output << "|Eta|:          ";
  double n_tot_MC = 0;
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
	if(hmc_response_probejet[i][j]->GetEntries()/1000 < 0.01)     output << hmc_response_probejet[i][j]->GetEntries() << "      - "; //<1000
	else if(hmc_response_probejet[i][j]->GetEntries()/1000 < 0.1) output << hmc_response_probejet[i][j]->GetEntries() << "     - "; //<1000
	else if(hmc_response_probejet[i][j]->GetEntries()/1000 < 1)   output << hmc_response_probejet[i][j]->GetEntries() << "    - "; //<1000
	else if(hmc_response_probejet[i][j]->GetEntries()/1000 <10)   output << hmc_response_probejet[i][j]->GetEntries() << "   - "; //<10000
	else if(hmc_response_probejet[i][j]->GetEntries()/1000 <100)  output << hmc_response_probejet[i][j]->GetEntries() << "  - ";
	else                                              output << hmc_response_probejet[i][j]->GetEntries() << " - ";
      }
      else output << hmc_response_probejet[i][j]->GetEntries() << endl;
      n_tot_MC+= hmc_response_probejet[i][j]->GetEntries();
    }

  }
  output << endl << endl << "Total number of matched events in MC: " << n_tot_MC << endl;
  */
 
 

  // Dump 1-d distributions of A and B in bins of pT, eta

  TFile* test_out_mc_B = new TFile(CorrectionObject::_outpath+"plots/control/GenResponse_1d_mc_matched.root","RECREATE");
  for(int j=0; j<n_eta-1; j++){
    for(int k=0; k<n_pt-1; k++){     ///k=0 n_pt-1 
      /*      hmc_response_probejet[k][j]->Write();
      hmc_response_probetagjet[k][j]->Write();
      hmc_response_probetagptcl[k][j]->Write();
      hmc_response_probetagreco[k][j]->Write();
      hmc_response_tagtagjet[k][j]->Write();*/
      hmc_jet1_genID[k][j]->Write();
      hmc_jet2_genID[k][j]->Write();
      hmc_jet3_genID[k][j]->Write();
      hmc_normjetpt[k][j]->Write();
      hmc_probejetpt[k][j]->Write();
      hmc_probegenjetpt[k][j]->Write();
      hmc_tagjetpt[k][j]->Write();
      hmc_taggenjetpt[k][j]->Write();
      hmc_A[k][j]->Write();
      for(int ifl=0;ifl<3;ifl++){
	hmc_probejetpt_flavor[ifl][j]->Write();
	hmc_tagjetpt_flavor[ifl][j]->Write();
      }
    }
  }
  test_out_mc_B->Close();
  delete test_out_mc_B;


  //R_MC as a function of pT, in bins of |eta|
  /* double val_rel_mc[n_eta-1][n_pt-1]; //value at pt,eta
  double err_rel_mc[n_eta-1][n_pt-1]; //error of ratio at pt,eta
  double val_tagprobe_mc[n_eta-1][n_pt-1]; //value at pt,eta
  double err_tagprobe_mc[n_eta-1][n_pt-1]; //error of ratio at pt,eta
  double val_tagprobe_ptcl_mc[n_eta-1][n_pt-1]; //value at pt,eta
  double err_tagprobe_ptcl_mc[n_eta-1][n_pt-1]; //error of ratio at pt,eta
  double val_tagprobe_reco_mc[n_eta-1][n_pt-1]; //value at pt,eta
  double err_tagprobe_reco_mc[n_eta-1][n_pt-1]; //error of ratio at pt,eta
  double val_tagtag_mc[n_eta-1][n_pt-1]; //value at pt,eta
  double err_tagtag_mc[n_eta-1][n_pt-1]; //error of ratio at pt,eta

  double val_rel_A_mc[n_eta-1][n_pt-1]; //value at pt,eta
  double err_rel_A_mc[n_eta-1][n_pt-1]; //error of ratio at pt,eta*/

  double val_rel_A_mc[n_eta-1][n_pt-1]; //value at pt,eta
  double err_rel_A_mc[n_eta-1][n_pt-1]; //error of ratio at pt,eta

  double val_probejet_pt[n_eta-1][n_pt-1]; //value at pt,eta
  double err_probejet_pt[n_eta-1][n_pt-1]; //value at pt,eta
  double val_probegenjet_pt[n_eta-1][n_pt-1]; //value at pt,eta
  double err_probegenjet_pt[n_eta-1][n_pt-1]; //value at pt,eta
  double val_tagjet_pt[n_eta-1][n_pt-1]; //value at pt,eta
  double err_tagjet_pt[n_eta-1][n_pt-1]; //value at pt,eta
  double val_taggenjet_pt[n_eta-1][n_pt-1]; //value at pt,eta
  double err_taggenjet_pt[n_eta-1][n_pt-1]; //value at pt,eta

  double val_probeRECO_probeGEN[n_eta-1][n_pt-1];
  double err_probeRECO_probeGEN[n_eta-1][n_pt-1];
  double val_tagRECO_tagGEN[n_eta-1][n_pt-1];
  double err_tagRECO_tagGEN[n_eta-1][n_pt-1];

  for(int i=0; i<n_eta-1; i++){
    for(int j=0; j<n_pt-1; j++){
      //get <response> and error on <response>
      pair <double,double> A_mc = GetValueAndError(hmc_A[j][i]);
      pair<double,double> res_mc_rel_r;
      res_mc_rel_r.first = (1+A_mc.first)/(1-A_mc.first);
      res_mc_rel_r.second = 2/(pow((1-A_mc.first),2)) * A_mc.second;
      val_rel_A_mc[i][j] = res_mc_rel_r.first;
      err_rel_A_mc[i][j] = res_mc_rel_r.second;


      pair <double,double> normpt_mc = GetValueAndError(hmc_normjetpt[j][i]); //<pt_bin> value used for normalisation
      pair <double,double> probejetpt_mc = GetValueAndError(hmc_probejetpt[j][i]);
      val_probejet_pt[j][i] = probejetpt_mc.first*normpt_mc.first;
      err_probejet_pt[j][i] =  ErrorPropagation_AB(probejetpt_mc,normpt_mc);
      pair <double,double> probegenjetpt_mc = GetValueAndError(hmc_probegenjetpt[j][i]);
      val_probegenjet_pt[j][i] = probegenjetpt_mc.first*normpt_mc.first;
      err_probegenjet_pt[j][i] = ErrorPropagation_AB(probegenjetpt_mc,normpt_mc);
      pair <double,double> tagjetpt_mc = GetValueAndError(hmc_tagjetpt[j][i]);
      val_tagjet_pt[j][i] = tagjetpt_mc.first*normpt_mc.first;
      err_tagjet_pt[j][i] = ErrorPropagation_AB(tagjetpt_mc,normpt_mc);
      pair <double,double> taggenjetpt_mc = GetValueAndError(hmc_taggenjetpt[j][i]);
      val_taggenjet_pt[j][i] = taggenjetpt_mc.first*normpt_mc.first;
      err_taggenjet_pt[j][i] = ErrorPropagation_AB(taggenjetpt_mc,normpt_mc);

      if(val_probegenjet_pt[j][i]>0){
      val_probeRECO_probeGEN[j][i] = val_probejet_pt[j][i]/val_probegenjet_pt[j][i];
      pair<double,double> tmp1; tmp1.first = val_probejet_pt[j][i]; tmp1.second = err_probejet_pt[j][i];
      pair<double,double> tmp2; tmp1.first = val_probegenjet_pt[j][i]; tmp1.second = err_probegenjet_pt[j][i];
      err_probeRECO_probeGEN[j][i] = ErrorPropagation_AoverB(tmp1,tmp2);
      }
      else{
	val_probeRECO_probeGEN[j][i] = 0;  err_probeRECO_probeGEN[j][i] =0;
      }
      if(val_taggenjet_pt[j][i]>0){
	val_tagRECO_tagGEN[j][i] = val_tagjet_pt[j][i]/val_taggenjet_pt[j][i];
	pair<double,double> tmp1; tmp1.first = val_tagjet_pt[j][i]; tmp1.second = err_tagjet_pt[j][i];
	pair<double,double> tmp2; tmp1.first = val_taggenjet_pt[j][i]; tmp1.second = err_taggenjet_pt[j][i];
	err_tagRECO_tagGEN[j][i] = ErrorPropagation_AoverB(tmp1,tmp2);

      }
      else{
	val_tagRECO_tagGEN[j][i] = 0;  err_tagRECO_tagGEN[j][i] =0;
      }
      
      cout<<" val_probeRECO_probeGEN[j][i] = "<<val_probeRECO_probeGEN[j][i]<<" "<<probejetpt_mc.first<<" "<<normpt_mc.first<<endl;
      cout<<" err_probeRECO_probeGEN[j][i] = "<<err_probeRECO_probeGEN[j][i]<<" "<<probejetpt_mc.second<<" "<<normpt_mc.second<<endl;
      cout<<""<<endl;
      //      cout<<" val_tagRECO_tagGEN[j][i] = "<<val_tagRECO_tagGEN[j][i]<<" "<<tagjetpt_mc.first<<" "<<normpt_mc.first<<endl;
      
      // //get <response> and error on <response>
      // pair <double,double> ResGEN_mc = GetValueAndError(hmc_response_probejet[j][i]);
      // val_rel_mc[i][j] = ResGEN_mc.first;
      // err_rel_mc[i][j] = ResGEN_mc.second;

      // //get <response> and error on <response>
      // pair <double,double> Res_mc = GetValueAndError(hmc_response_probetagjet[j][i]);
      // //      cout<<" Res_mc.first = "<<Res_mc.first<<" Res_mc.second = "<<Res_mc.second<<endl;
      // val_tagprobe_mc[i][j] = Res_mc.first;
      // err_tagprobe_mc[i][j] = Res_mc.second;

      // //get <response> and error on <response>
      // pair <double,double> ResPTCL_mc = GetValueAndError(hmc_response_probetagptcl[j][i]);
      // val_tagprobe_ptcl_mc[i][j] = ResPTCL_mc.first;
      // err_tagprobe_ptcl_mc[i][j] = ResPTCL_mc.second;

      // pair <double,double> ResRECO_mc = GetValueAndError(hmc_response_probetagreco[j][i]);
      // val_tagprobe_reco_mc[i][j] = ResRECO_mc.first;
      // err_tagprobe_reco_mc[i][j] = ResRECO_mc.second;

      // pair <double,double> ResTAG_mc = GetValueAndError(hmc_response_tagtagjet[j][i]);
      // //      cout<<"ResTAG_mc.first = "<<ResTAG_mc.first<<" ResTAG_mc.second = "<<ResTAG_mc.second<<endl;
      // val_tagtag_mc[i][j] = ResTAG_mc.first;
      // err_tagtag_mc[i][j] = ResTAG_mc.second;
      // //      cout<<"val_tagtag_A_mc[i][j] = "<<val_tagtag_A_mc[i][j]<<" err_tagtag_A_mc[i][j] = "<<err_tagtag_A_mc[i][j]<<endl;


      
    }
  }

  //dummy for tdrCanvas
  TH1D *h = new TH1D("h",";dummy;",41,0,5.191);
  h->SetMaximum(1.2);
  h->SetMinimum(0.8);

  TH1D *hEF = new TH1D("hEF",";dummy;",1000,0,5.191);

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

    TString alVal;
    alVal.Form("%0.2f\n",alpha_cut);
    TString altitle = "{#alpha<"+alVal+"}";

    TGraphErrors *graph_rel_A_mc   = new TGraphErrors(n_pt-1, xbin_tgraph, val_rel_A_mc[i], zero, err_rel_A_mc[i]);
    graph_rel_A_mc   = (TGraphErrors*)CleanEmptyPoints(graph_rel_A_mc);
    //    cout<<"graph_rel_A_mc"<<endl;
    //    graph_rel_A_mc->Print();
    graph_rel_A_mc->SetTitle("");
    // graph_rel_A_mc->GetXaxis()->SetTitle(pt_binning_var_str);
    // graph_rel_A_mc->GetXaxis()->SetTitleSize(0.05);
    // graph_rel_A_mc->GetXaxis()->SetTitleOffset(0.80);
    // graph_rel_A_mc->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    // graph_rel_A_mc->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_rel_A_mc->SetMarkerColor(kOrange+7);
    graph_rel_A_mc->SetMarkerStyle(29);
    graph_rel_A_mc->SetMarkerSize(1.7);
    graph_rel_A_mc->SetLineColor(kOrange+7);
    TString axistitle_A_mc = "(1+<A>)/(1-<A>)";

    TGraphErrors *graph_probeRECO_probeGEN   = new TGraphErrors(n_pt-1, xbin_tgraph, val_probeRECO_probeGEN[i], zero, err_probeRECO_probeGEN[i]);
    //    graph_probeRECO_probeGEN   = (TGraphErrors*)CleanEmptyPoints(graph_probeRECO_probeGEN);
    cout<<"graph_probeRECO_probeGEN"<<endl;
    graph_probeRECO_probeGEN->Print();
    graph_probeRECO_probeGEN->SetTitle("");
    graph_probeRECO_probeGEN->SetMarkerColor(kRed);
    graph_probeRECO_probeGEN->SetMarkerStyle(20);
    graph_probeRECO_probeGEN->SetLineColor(kRed);
    TString axistitle_mc_probeprobe = "<p^{probe,RECO}_{T}>/<p^{probe,GEN}_{T}>";

    TGraphErrors *graph_tagRECO_tagGEN   = new TGraphErrors(n_pt-1, xbin_tgraph, val_tagRECO_tagGEN[i], zero, err_tagRECO_tagGEN[i]);
    //    graph_tagRECO_tagGEN   = (TGraphErrors*)CleanEmptyPoints(graph_tagRECO_tagGEN);
    cout<<"graph_tagRECO_tagGEN"<<endl;
    graph_tagRECO_tagGEN->Print();
    graph_tagRECO_tagGEN->SetTitle("");
    graph_tagRECO_tagGEN->SetMarkerColor(kBlue);
    graph_tagRECO_tagGEN->SetMarkerStyle(20);
    graph_tagRECO_tagGEN->SetLineColor(kBlue);
    TString axistitle_mc_tagtag = "<p^{tag,RECO}_{T}>/<p^{tag,GEN}_{T}>";

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    TString text = eta_range[i] + " < |#eta| < " + eta_range[i+1];

    TCanvas* c_rel = new TCanvas();
    tdrCanvas(c_rel,"c_rel",h,4,10,true,CorrectionObject::_lumitag);
    h->GetXaxis()->SetTitle(pt_binning_var_str);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    h->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_rel_A_mc->Draw("P SAME");
    graph_probeRECO_probeGEN->Draw("P SAME");
    graph_tagRECO_tagGEN->Draw("P SAME");

    gPad->SetLogx();
    TLegend *leg_rel;
    leg_rel = new TLegend(0.35,0.15,0.91,0.49,"","brNDC");//x+0.1
    leg_rel->SetBorderSize(0);
    leg_rel->SetTextSize(0.030);
    leg_rel->SetFillColor(10);
    leg_rel->SetFillStyle(0);
    leg_rel->SetLineColor(1);
    leg_rel->SetTextFont(42);
    leg_rel->SetHeader("R^{MC}_"+altitle+", "+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]); 
    leg_rel->AddEntry(graph_rel_A_mc, axistitle_A_mc,"P");
    leg_rel->AddEntry(graph_probeRECO_probeGEN, axistitle_mc_probeprobe,"P");
    leg_rel->AddEntry(graph_tagRECO_tagGEN, axistitle_mc_tagtag,"P");
    leg_rel->Draw();
    c_rel->SaveAs(CorrectionObject::_outpath+"plots/control/GenResponse_RatioOfAverages_matchedProbeJet_"+pt_binning_var_name+ CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + ".pdf");

    delete graph_probeRECO_probeGEN;
  }

    /* TGraphErrors *graph_rel_A_mc   = new TGraphErrors(n_pt-1, xbin_tgraph, val_rel_A_mc[i], zero, err_rel_A_mc[i]);
    graph_rel_A_mc   = (TGraphErrors*)CleanEmptyPoints(graph_rel_A_mc);

    graph_rel_A_mc->SetTitle("");
    // graph_rel_A_mc->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph_rel_A_mc->GetXaxis()->SetTitle(pt_binning_var_str);
    graph_rel_A_mc->GetXaxis()->SetTitleSize(0.05);
    graph_rel_A_mc->GetXaxis()->SetTitleOffset(0.80);
    graph_rel_A_mc->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    graph_rel_A_mc->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_rel_A_mc->SetMarkerColor(kOrange+7);
    graph_rel_A_mc->SetMarkerStyle(29);
    graph_rel_A_mc->SetMarkerSize(1.7);
    graph_rel_A_mc->SetLineColor(kOrange+7);

    TString alVal;
    alVal.Form("%0.2f\n",alpha_cut);
    TString altitle = "{#alpha<"+alVal+"}";
    // TString axistitle_A_mc = "R^{MC}_"; 
    // axistitle_A_mc += altitle;
    // axistitle_A_mc += "((1+<A>)/(1-<A>))";
    TString axistitle_A_mc = "(1+<A>)/(1-<A>)";

    TGraphErrors *graph_rel_mc   = new TGraphErrors(n_pt-1, xbin_tgraph, val_rel_mc[i], zero, err_rel_mc[i]);
    graph_rel_mc   = (TGraphErrors*)CleanEmptyPoints(graph_rel_mc);

    graph_rel_mc->SetTitle("");
    // //    graph_rel_mc->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    // graph_rel_mc->GetXaxis()->SetTitle(pt_binning_var_str);
    // graph_rel_mc->GetXaxis()->SetTitleSize(0.05);
    // graph_rel_mc->GetXaxis()->SetTitleOffset(0.80);
    // graph_rel_mc->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    // graph_rel_mc->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_rel_mc->SetMarkerColor(kRed);
    graph_rel_mc->SetMarkerStyle(20);
    graph_rel_mc->SetLineColor(kRed);

    // TString axistitle_mc = "R^{MC}_"; 
    // axistitle_mc   += altitle;
    // axistitle_mc += "(<p^{probe,RECO}_{T}/p^{probe,GEN}_{T}>)";
    TString axistitle_mc = "<p^{probe,RECO}_{T}/p^{probe,GEN}_{T}>";

    TGraphErrors *graph_tagprobe_mc   = new TGraphErrors(n_pt-1, xbin_tgraph, val_tagprobe_mc[i], zero, err_tagprobe_mc[i]);
    graph_tagprobe_mc   = (TGraphErrors*)CleanEmptyPoints(graph_tagprobe_mc);

    graph_tagprobe_mc->SetTitle("");
    //    graph_tagprobe_mc->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    // graph_tagprobe_mc->GetXaxis()->SetTitle(pt_binning_var_str);
    // graph_tagprobe_mc->GetXaxis()->SetTitleSize(0.05);
    // graph_tagprobe_mc->GetXaxis()->SetTitleOffset(0.80);
    // graph_tagprobe_mc->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    //    graph_tagprobe_mc->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_tagprobe_mc->SetMarkerColor(kBlack);
    graph_tagprobe_mc->SetMarkerSize(1.1);
    graph_tagprobe_mc->SetMarkerStyle(20);
    graph_tagprobe_mc->SetLineColor(kBlack);

    // TString axistitle_tagprobe = "R^{MC}_"; 
    // axistitle_tagprobe   += altitle;
    // axistitle_tagprobe += "(<p^{probe,GEN}_{T}/p^{tag,GEN}_{T}>)";
    TString axistitle_tagprobe = "<p^{probe,GEN}_{T}/p^{tag,GEN}_{T}>";
    TGraphErrors *graph_tagprobe_ptcl_mc   = new TGraphErrors(n_pt-1, xbin_tgraph, val_tagprobe_ptcl_mc[i], zero, err_tagprobe_ptcl_mc[i]);
    //    graph_tagprobe_ptcl_mc   = (TGraphErrors*)CleanEmptyPoints(graph_tagprobe_ptcl_mc);

    graph_tagprobe_ptcl_mc->SetTitle("");
    //    graph_tagprobe_ptcl_mc->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    // graph_tagprobe_ptcl_mc->GetXaxis()->SetTitle(pt_binning_var_str);
    // graph_tagprobe_ptcl_mc->GetXaxis()->SetTitleSize(0.05);
    // graph_tagprobe_ptcl_mc->GetXaxis()->SetTitleOffset(0.80);
    // graph_tagprobe_ptcl_mc->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    // graph_tagprobe_ptcl_mc->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_tagprobe_ptcl_mc->SetMarkerColor(kGreen);
    graph_tagprobe_ptcl_mc->SetMarkerStyle(20);
    graph_tagprobe_ptcl_mc->SetLineColor(kGreen);

    // TString axistitle_tagprobe_ptcl = "R^{MC}_"; 
    // axistitle_tagprobe_ptcl   += altitle;
    // axistitle_tagprobe_ptcl += "(<p^{probe,PAR}_{T}/p^{tag,PAR}_{T}>)";
    TString axistitle_tagprobe_ptcl = "<p^{probe,PAR}_{T}/p^{tag,PAR}_{T}>";
    TGraphErrors *graph_tagprobe_reco_mc   = new TGraphErrors(n_pt-1, xbin_tgraph, val_tagprobe_reco_mc[i], zero, err_tagprobe_reco_mc[i]);
    graph_tagprobe_reco_mc   = (TGraphErrors*)CleanEmptyPoints(graph_tagprobe_reco_mc);

    graph_tagprobe_reco_mc->SetTitle("");
    //    graph_tagprobe_reco_mc->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    // graph_tagprobe_reco_mc->GetXaxis()->SetTitle(pt_binning_var_str);
    // graph_tagprobe_reco_mc->GetXaxis()->SetTitleSize(0.05);
    // graph_tagprobe_reco_mc->GetXaxis()->SetTitleOffset(0.80);
    // graph_tagprobe_reco_mc->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    // graph_tagprobe_reco_mc->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_tagprobe_reco_mc->SetMarkerColor(kMagenta);
    graph_tagprobe_reco_mc->SetMarkerStyle(20);
    graph_tagprobe_reco_mc->SetLineColor(kMagenta);

    // TString axistitle_tagprobe_reco = "R^{MC}_"; 
    // axistitle_tagprobe_reco   += altitle;
    // axistitle_tagprobe_reco += "(<p^{probe,RECO}_{T}/p^{tag,RECO}_{T}>)";
    TString axistitle_tagprobe_reco = "<p^{probe,RECO}_{T}/p^{tag,RECO}_{T}>";

    TGraphErrors *graph_tagtag_mc   = new TGraphErrors(n_pt-1, xbin_tgraph, val_tagtag_mc[i], zero, err_tagtag_mc[i]);
    graph_tagtag_mc   = (TGraphErrors*)CleanEmptyPoints(graph_tagtag_mc);
    // cout<<"TAGTAG graph"<<endl;
    // graph_tagtag_mc->Print();
    graph_tagtag_mc->SetTitle("");

    //    graph_tagtag_mc->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    // graph_tagtag_mc->GetXaxis()->SetTitle(pt_binning_var_str);
    // graph_tagtag_mc->GetXaxis()->SetTitleSize(0.05);
    // graph_tagtag_mc->GetXaxis()->SetTitleOffset(0.80);
    // graph_tagtag_mc->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    // graph_tagtag_mc->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_tagtag_mc->SetMarkerColor(kBlue);
    graph_tagtag_mc->SetMarkerStyle(29);
    graph_tagtag_mc->SetMarkerSize(1.7);
    graph_tagtag_mc->SetLineColor(kBlue);

    // TString axistitle_tagtag = "R^{MC}_"; 
    // axistitle_tagtag   += altitle;
    // axistitle_tagtag += "(<p^{tag,RECO}_{T}/p^{tag,GEN}_{T}>)";
    TString axistitle_tagtag = "<p^{tag,RECO}_{T}/p^{tag,GEN}_{T}>";

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    TString text = eta_range[i] + " < |#eta| < " + eta_range[i+1];

    TCanvas* c_rel = new TCanvas();
    tdrCanvas(c_rel,"c_rel",h,4,10,true,CorrectionObject::_lumitag);
    //    h->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    h->GetXaxis()->SetTitle(pt_binning_var_str);
    h->GetXaxis()->SetTitleSize(0.05);
    //    h->GetXaxis()->SetTitleOffset(0.90);
    h->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    h->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_rel_A_mc->Draw("P SAME");
    graph_tagtag_mc->Draw("P SAME");
    graph_rel_mc->Draw("P SAME");
    graph_tagprobe_mc->Draw("P SAME");
    graph_tagprobe_ptcl_mc->Draw("P SAME");
    graph_tagprobe_reco_mc->Draw("P SAME");


    // graph_tagtag_mc->Print();
    // graph_rel_A_mc->Print();

    gPad->SetLogx();

    TLegend *leg_rel;
    //    leg_rel = new TLegend(0.35,0.50,0.91,0.89,"","brNDC");//x+0.1
    leg_rel = new TLegend(0.35,0.15,0.91,0.49,"","brNDC");//x+0.1
    //    leg_rel->SetNColumns(2);
    leg_rel->SetBorderSize(0);
    leg_rel->SetTextSize(0.030);
    leg_rel->SetFillColor(10);
    leg_rel->SetFillStyle(0);
    leg_rel->SetLineColor(1);
    leg_rel->SetTextFont(42);
    leg_rel->SetHeader("R^{MC}_"+altitle+", "+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]); 
    leg_rel->AddEntry(graph_rel_A_mc, axistitle_A_mc,"P");
    leg_rel->AddEntry(graph_rel_mc, axistitle_mc,"P");
    leg_rel->AddEntry(graph_tagprobe_mc, axistitle_tagprobe,"P");
    leg_rel->AddEntry(graph_tagprobe_ptcl_mc, axistitle_tagprobe_ptcl,"P");
    leg_rel->AddEntry(graph_tagprobe_reco_mc, axistitle_tagprobe_reco,"P");
    leg_rel->AddEntry(graph_tagtag_mc, axistitle_tagtag,"P");
    leg_rel->Draw();
    //tex->DrawLatex(0.53,0.91,CorrectionObject::_lumitag+"(13TeV)");

    c_rel->SaveAs(CorrectionObject::_outpath+"plots/control/GenResponse_matchedProbeJet_"+pt_binning_var_name+ CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + ".pdf");


    //delete leg_rel;
    delete c_rel;
    //delete leg_mpf;
    delete tex;
    delete graph_rel_mc;
   
  }
   */
  // //Plot 1d genID->OFFLINE jet distributions in a particular eta-bin for different pt-bins on different canvases
  // for(int i=0; i<n_eta-1; i++){
  //   TString eta_name = "eta_"+eta_range2[i]+"_"+eta_range2[i+1];
    
  //   TLatex *tex = new TLatex();
  //   tex->SetNDC();
  //   tex->SetTextSize(0.045); 
  //   TString text = eta_range[i] + " < |#eta| < " + eta_range[i+1];

  //   TLatex *tex_lumi = new TLatex();
  //   tex_lumi->SetNDC();
  //   tex_lumi->SetTextSize(0.045); 
  //   for(int j=0; j<n_pt-1; j++){   ///j=0 j<pt_n-1
  //     TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
  //     TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
  //     int n_ev1 = hmc_jet1_genID[j][i]->GetEntries();
  //     int n_ev2 = hmc_jet2_genID[j][i]->GetEntries();
  //     int n_ev3 = hmc_jet3_genID[j][i]->GetEntries();
  //     // cout<<"Scale = "<<1/hmc_jet1_genID[j][i]->Integral()<<endl;
  //     // if(hmc_jet1_genID[j][i]->Integral() > 0) hmc_jet1_genID[j][i]->Scale(1/hmc_jet1_genID[j][i]->Integral());
  //     // if(hmc_jet2_genID[j][i]->Integral() > 0) hmc_jet2_genID[j][i]->Scale(1/hmc_jet2_genID[j][i]->Integral());
  //     // if(hmc_jet3_genID[j][i]->Integral() > 0) hmc_jet3_genID[j][i]->Scale(1/hmc_jet3_genID[j][i]->Integral());

  //     hmc_jet1_genID[j][i]->Print();
  //     //      if(htemp_mpf_mc->Integral() > 0)htemp_mpf_mc->Scale(1/htemp_mpf_mc->Integral());

  //     if(hmc_jet1_genID[j][i]->Integral() > 0)
  //     h->GetXaxis()->SetTitle("genID-recoID");
  //     //      h->GetYaxis()->SetTitle("Normalized entries");
  //     h->GetYaxis()->SetTitleOffset(1.5);
  //     h->GetXaxis()->SetLimits(-5,10);
  //     h->SetMinimum(1e-1);
  //     h->SetMaximum(1e15);
  //     hmc_jet1_genID[j][i]->SetLineColor(kRed);
  //     hmc_jet1_genID[j][i]->SetMarkerColor(kRed);
  //     hmc_jet1_genID[j][i]->SetLineWidth(3);
  //     hmc_jet1_genID[j][i]->SetMarkerStyle(20);
  //     hmc_jet1_genID[j][i]->SetFillColorAlpha(kRed,0.99);

  //     hmc_jet2_genID[j][i]->SetLineColor(kGreen);
  //     hmc_jet2_genID[j][i]->SetMarkerColor(kGreen);
  //     hmc_jet2_genID[j][i]->SetLineWidth(3);
  //     hmc_jet2_genID[j][i]->SetMarkerStyle(20);
  //     hmc_jet2_genID[j][i]->SetFillColorAlpha(kGreen,0.5);
  //     hmc_jet2_genID[j][i]->SetFillStyle(3005);
  //     //      hmc_jet2_genID[j][i]->SetFillColor(3005kGreen);

  //     hmc_jet3_genID[j][i]->SetLineColor(kBlack);
  //     hmc_jet3_genID[j][i]->SetMarkerColor(kBlack);
  //     hmc_jet3_genID[j][i]->SetLineWidth(3);
  //     hmc_jet3_genID[j][i]->SetMarkerStyle(20);
  //     hmc_jet3_genID[j][i]->SetFillColorAlpha(kBlack,0.25);
      
  //     TCanvas* ctmp = new TCanvas();
  //     tdrCanvas(ctmp,"ctmp",h,4,10,kSquare,CorrectionObject::_lumitag);
  //     TLegend leg2 = tdrLeg(0.35,0.6,0.90,0.89);
  //     leg2.SetHeader(""+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]+", "+legname); 
  //     leg2.AddEntry(hmc_jet1_genID[j][i], "RECO jet1", "l");
  //     leg2.AddEntry(hmc_jet2_genID[j][i], "RECO jet2", "l");
  //     leg2.AddEntry(hmc_jet3_genID[j][i], "RECO jet3", "l");
  //     if(n_ev1>0) hmc_jet1_genID[j][i]->Draw("HIST SAME");
  //     if(n_ev2>0) hmc_jet2_genID[j][i]->Draw("HIST SAME");
  //     if(n_ev3>0) hmc_jet3_genID[j][i]->Draw("HIST SAME");
  //     leg2.Draw();
  //     gPad->SetLogy();
  //     //      tex->DrawLatex(0.47,0.85,"MC, " + text);
  //     ctmp->SaveAs(CorrectionObject::_outpath+"plots/control/Matched_OFFLINEJet_genID_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] 
  // 		   +"_" +pt_name +".pdf");
  //   }
  // }
  
  //Plot Probe jet pt for different flavors on one canvas per eta bin
  for(int i=0; i<n_eta-1; i++){

    TString eta_name = "eta_"+eta_range2[i]+"_"+eta_range2[i+1];
    
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    TString text = eta_range[i] + " < |#eta| < " + eta_range[i+1];

    TLatex *tex_lumi = new TLatex();
    tex_lumi->SetNDC();
    tex_lumi->SetTextSize(0.045); 
    TCanvas* ctmp = new TCanvas();
    tdrCanvas(ctmp,"ctmp",h,4,10,kSquare,CorrectionObject::_lumitag);
    TLegend leg2 = tdrLeg(0.35,0.6,0.90,0.89);
    leg2.SetHeader(""+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]); 
    // TH1D *h_sum = hmc_probejetpt_flavor[0][i];
    // h_sum->Add(hmc_probejetpt_flavor[1][i]);
    // h_sum->Add(hmc_probejetpt_flavor[2][i]);
    for(int j=1; j<3; j++){   //jet flavor
    //    for(int j=2; j>-1; j--){   //jet flavor
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      int n_ev = hmc_probejetpt_flavor[j][i]->GetEntries();
      //      if(hmc_probejetpt_flavor[j][i]->Integral() > 0)hmc_probejetpt_flavor[j][i]->Scale(1/hmc_probejetpt_flavor[j][i]->Integral());
      //      double scale = h_sum->Integral()/hmc_probejetpt_flavor[j][i]->Integral();
      //      cout<<"scale = "<<scale<<endl;
      //      if(h_sum->Integral() > 0 && hmc_probejetpt_flavor[j][i]->Integral() > 0) hmc_probejetpt_flavor[j][i]->Scale(h_sum->Integral()/hmc_probejetpt_flavor[j][i]->Integral());
      h->GetXaxis()->SetTitle("RECO probe jet p_{T}, GeV");
      //      h->GetYaxis()->SetTitle("Normalized entries");
      //      h->GetYaxis()->SetTitle("Normalized Fraction");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(0,500);
      h->SetMinimum(1e-2);
      //      h->SetMaximum(1.0);
      h->SetMaximum(1e14);

      hmc_probejetpt_flavor[j][i]->SetLineColor(kRed+2*j);
      hmc_probejetpt_flavor[j][i]->SetMarkerColor(kRed+2*j);
      // hmc_probejetpt_flavor[j][i]->SetLineColor(1+2*j);
      // hmc_probejetpt_flavor[j][i]->SetMarkerColor(1+2*j);
      hmc_probejetpt_flavor[j][i]->SetLineWidth(3);
      hmc_probejetpt_flavor[j][i]->SetMarkerStyle(20);
      hmc_probejetpt_flavor[j][i]->SetFillColorAlpha(kRed+2*j,1-(j+1)/3.);
      //      hmc_probejetpt_flavor[j][i]->SetFillColorAlpha(1+2*j,1-0.1*(j+1)/3.);
      if(j==0){
	hmc_probejetpt_flavor[j][i]->SetLineColor(kBlue+2*j);
	hmc_probejetpt_flavor[j][i]->SetMarkerColor(kBlue+2*j);
	hmc_probejetpt_flavor[j][i]->SetFillColorAlpha(kBlue+2*j,1-(j+1)/3.);
      }

      TString legtext;
      if(j==0) legtext = "unmatched";
      if(j==1) legtext = "Quarks";
      if(j==2) legtext = "Gluons";
      leg2.AddEntry(hmc_probejetpt_flavor[j][i], legtext, "fl");
      if(n_ev>0) hmc_probejetpt_flavor[j][i]->Draw("HIST SAME");
    }
      leg2.Draw();
      gPad->SetLogy();
      //      tex->DrawLatex(0.47,0.85,"MC, " + text);
      ctmp->SaveAs(CorrectionObject::_outpath+"plots/control/Matched_ProbeJetpt_flavor_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1]+".pdf");
  }

//Plot Tag jet pt for different flavors on one canvas per eta bin
  for(int i=0; i<n_eta-1; i++){

    TString eta_name = "eta_"+eta_range2[i]+"_"+eta_range2[i+1];
    
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    TString text = eta_range[i] + " < |#eta| < " + eta_range[i+1];

    TLatex *tex_lumi = new TLatex();
    tex_lumi->SetNDC();
    tex_lumi->SetTextSize(0.045); 
    TCanvas* ctmp = new TCanvas();
    tdrCanvas(ctmp,"ctmp",h,4,10,kSquare,CorrectionObject::_lumitag);
    TLegend leg2 = tdrLeg(0.35,0.6,0.90,0.89);
    //    TLegend leg2 = tdrLeg(0.35,0.1,0.90,0.49);
    leg2.SetHeader(""+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]); 
    TH1D *h_sum = hmc_tagjetpt_flavor[0][i];
    h_sum->Add(hmc_tagjetpt_flavor[1][i]);
    h_sum->Add(hmc_tagjetpt_flavor[2][i]);
    for(int j=1; j<3; j++){   //jet flavor
    //    for(int j=2; j>-1; j--){   //jet flavor
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      int n_ev = hmc_tagjetpt_flavor[j][i]->GetEntries();
      //      if(h_sum->Integral() > 0) hmc_tagjetpt_flavor[j][i]->Scale(1/h_sum->Integral());
      h->GetXaxis()->SetTitle("RECO tag jet p_{T}, GeV");
      //      h->GetYaxis()->SetTitle("Normalized Fraction");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(0,500);
      h->SetMinimum(1e-2);
      //      h->SetMaximum(100.0);
      h->SetMaximum(1e14);

      hmc_tagjetpt_flavor[j][i]->SetLineColor(kRed+2*j);
      hmc_tagjetpt_flavor[j][i]->SetMarkerColor(kRed+2*j);
      // hmc_tagjetpt_flavor[j][i]->SetLineColor(1+2*j);
      // hmc_tagjetpt_flavor[j][i]->SetMarkerColor(1+2*j);
      hmc_tagjetpt_flavor[j][i]->SetLineWidth(3);
      hmc_tagjetpt_flavor[j][i]->SetMarkerStyle(20);
      hmc_tagjetpt_flavor[j][i]->SetFillColorAlpha(kRed+2*j,1-(j+1)/3.);
      if(j==0){
	hmc_tagjetpt_flavor[j][i]->SetLineColor(kBlue+2*j);
	hmc_tagjetpt_flavor[j][i]->SetMarkerColor(kBlue+2*j);
	hmc_tagjetpt_flavor[j][i]->SetFillColorAlpha(kBlue+2*j,1-(j+1)/3.);
      }
      //      hmc_tagjetpt_flavor[j][i]->SetFillColorAlpha(1+2*j,1-0.1*(j+1)/3.);


      TString legtext;
      if(j==0) legtext = "unmatched";
      if(j==1) legtext = "Quarks";
      if(j==2) legtext = "Gluons";
      leg2.AddEntry(hmc_tagjetpt_flavor[j][i], legtext, "fl");
      if(n_ev>0) hmc_tagjetpt_flavor[j][i]->Draw("HIST SAME");
    }
      leg2.Draw();
      gPad->SetLogy();
      //      tex->DrawLatex(0.47,0.85,"MC, " + text);
      ctmp->SaveAs(CorrectionObject::_outpath+"plots/control/Matched_TagJetpt_flavor_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1]+".pdf");
  }

  
  
}
