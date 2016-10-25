#include "header.h"
#include "tdrstyle_mod15.C"


void steer(){


  gROOT->ProcessLine("bool lumi_binned;");
  gROOT->ProcessLine("bool with_sys;");
  gROOT->ProcessLine("bool is_sys_MC;"); //Systematic only relevant for MC samples, no systematic data files need to be read in

  gROOT->ProcessLine("TString s_jets;");
  gROOT->ProcessLine("TString jettag;");
  gROOT->ProcessLine("TString s_postfix;");
  gROOT->ProcessLine("TString runnr;");
  gROOT->ProcessLine("TString txttag;");
  gROOT->ProcessLine("TString lumitag;");
  gROOT->ProcessLine("TString Generator;");
  gROOT->ProcessLine("TString s_JetID;");
  gROOT->ProcessLine("bool b_JetID;");

  gROOT->ProcessLine("TString dir;");


  gROOT->ProcessLine("TString path;");
  gROOT->ProcessLine("TString path_up;");
  gROOT->ProcessLine("TString path_down;");
  gROOT->ProcessLine("TString path_general;");

  // input data files

  gROOT->ProcessLine("TString realdatfile;");
  gROOT->ProcessLine("TString datfile;");
  gROOT->ProcessLine("TString mcfile;");
  

  //only for HIP mitigation plots
  gROOT->ProcessLine("TString lumis[3];");
  gROOT->ProcessLine("TString path_HIP;");
  gROOT->ProcessLine("TString HIP;");
  gROOT->ProcessLine("TString HIP_mtOff;");
  gROOT->ProcessLine("TString DATA_HIP;");
  gROOT->ProcessLine("TString DATA_mtOff;");


  gROOT->ProcessLine("TString tag;");

  gROOT->ProcessLine("TFile* Realdatafile;"); 
  gROOT->ProcessLine("TFile* datafile;"); 
  gROOT->ProcessLine("TFile* MCfile;"); 
  gROOT->ProcessLine("TFile* MCfile_up, *MCfile_down, *datafile_up, *datafile_down, *Realdatafile_up, *Realdatafile_down;");



  gROOT->ProcessLine("double al_cut;"); 
  gROOT->ProcessLine("int nResponseBins;");
  gROOT->ProcessLine("TString MCname;");

  gROOT->ProcessLine("TString variation;");
  gROOT->ProcessLine(".L HIP_Mitigation_Plots.C");
  gROOT->ProcessLine(".L Control_Plots.C");
  gROOT->ProcessLine(".L kFSR_pT_TTree.C+");
  gROOT->ProcessLine(".L PTextrapolation_TTree_kFSRfit.C+");
  gROOT->ProcessLine(".L L2ResOutput.C");  
  gROOT->ProcessLine(".L JetID_ComparisonPlots.C"); 
  
  gROOT->ProcessLine(".L run.C");



  
  //PYTHIA
  
  //gROOT->ProcessLine("run(\"pythia\", \"AK4CHS\", \"BCD\", false, false)");
  /*gROOT->ProcessLine("run(\"pythia\", \"AK4CHS\", \"EF\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK4CHS\", \"G\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK4CHS\", \"B\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK4CHS\", \"C\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK4CHS\", \"D\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK4CHS\", \"E\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK4CHS\", \"F\", false, true)");
  
  
  gROOT->ProcessLine("run(\"pythia\", \"AK4Puppi\", \"BCD\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK4Puppi\", \"EF\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK4Puppi\", \"G\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK4Puppi\", \"B\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK4Puppi\", \"C\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK4Puppi\", \"D\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK4Puppi\", \"E\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK4Puppi\", \"F\", false, true)");
  
  
  gROOT->ProcessLine("run(\"pythia\", \"AK8CHS\", \"BCD\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK8CHS\", \"EF\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK8CHS\", \"G\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK8CHS\", \"B\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK8CHS\", \"C\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK8CHS\", \"D\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK8CHS\", \"E\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK8CHS\", \"F\", false, true)");
  */
  
  //gROOT->ProcessLine("run(\"pythia\", \"AK8Puppi\", \"BCD\", false, true)");
  //gROOT->ProcessLine("run(\"pythia\", \"AK8Puppi\", \"EF\", false, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK8Puppi\", \"G\", false, true)");
  //gROOT->ProcessLine("run(\"pythia\", \"AK8Puppi\", \"B\", false, true)");
  //gROOT->ProcessLine("run(\"pythia\", \"AK8Puppi\", \"C\", false, true)");
  //gROOT->ProcessLine("run(\"pythia\", \"AK8Puppi\", \"D\", false, true)");
  //gROOT->ProcessLine("run(\"pythia\", \"AK8Puppi\", \"E\", false, true)");
  //gROOT->ProcessLine("run(\"pythia\", \"AK8Puppi\", \"F\", false, true)");
  
  //HERWIG
  /*
  gROOT->ProcessLine("run(\"herwig\", \"AK4CHS\", \"BCD\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK4CHS\", \"EF\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK4CHS\", \"G\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK4CHS\", \"B\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK4CHS\", \"C\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK4CHS\", \"D\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK4CHS\", \"E\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK4CHS\", \"F\", false, true)");

  gROOT->ProcessLine("run(\"herwig\", \"AK4Puppi\", \"BCD\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK4Puppi\", \"EF\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK4Puppi\", \"G\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK4Puppi\", \"B\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK4Puppi\", \"C\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK4Puppi\", \"D\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK4Puppi\", \"E\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK4Puppi\", \"F\", false, true)");

  gROOT->ProcessLine("run(\"herwig\", \"AK8CHS\", \"BCD\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8CHS\", \"EF\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8CHS\", \"G\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8CHS\", \"B\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8CHS\", \"C\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8CHS\", \"D\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8CHS\", \"E\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8CHS\", \"F\", false, true)");
  
  gROOT->ProcessLine("run(\"herwig\", \"AK8Puppi\", \"BCD\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8Puppi\", \"EF\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8Puppi\", \"G\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8Puppi\", \"B\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8Puppi\", \"C\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8Puppi\", \"D\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8Puppi\", \"E\", false, true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8Puppi\", \"F\", false, true)");
  */


  //Only general plots, the collection argument has no meaning.
  /*
  gROOT->ProcessLine("run(\"pythia\", \"AK8Puppi\", \"BCD\", true, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK8Puppi\", \"EF\", true, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK8Puppi\", \"G\", true, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK8Puppi\", \"B\", true, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK8Puppi\", \"C\", true, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK8Puppi\", \"D\", true, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK8Puppi\", \"E\", true, true)");
  gROOT->ProcessLine("run(\"pythia\", \"AK8Puppi\", \"F\", true, true)");
  

  
  gROOT->ProcessLine("run(\"herwig\", \"AK8Puppi\", \"BCD\", true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8Puppi\", \"EF\", true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8Puppi\", \"G\", true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8Puppi\", \"B\", true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8Puppi\", \"C\", true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8Puppi\", \"D\", true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8Puppi\", \"E\", true)");
  gROOT->ProcessLine("run(\"herwig\", \"AK8Puppi\", \"F\", true)");
  */








  // HIP studies: Generator, runnr, 1st bool is irrelevant, JetColl needs to be accurate, 2nd bool = false. Be sure to comment everything else except HIP macro in run.C

  //gROOT->ProcessLine("run(\"pythia\", \"AK4CHS\", \"BCD\", false, false)");
  //gROOT->ProcessLine("run(\"pythia\", \"AK4Puppi\", \"BCD\", false, false)");



  //JetID Comparison Plots
  //you may switch generator, collection, and runnr tags
  gROOT->ProcessLine("run(\"pythia\", \"AK4CHS\", \"BCD\", false, true)");



  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 
}


