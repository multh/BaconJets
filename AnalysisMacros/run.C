//Set of macros, which reads TTree from UHH2/BaconJet output files

void run(){

  // change your absolute path here:
  
    gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_noReweight_AllTriggers_TTree/\";");

  // input data files
  gROOT->ProcessLine("TString datfile = \"uhh2.AnalysisModuleRunner.DATA.RunD_AK4CHS.root\";");
  //  gROOT->ProcessLine("TString mcfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_pythia8_AK4CHS.root\";");
  gROOT->ProcessLine("TString mcfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_herwigpp_AK4CHS.root\";");

  // choose a tag for the txt output files
  gROOT->ProcessLine("TString txttag = \"v2\";");
  // tag for the jet cone size
  gROOT->ProcessLine("TString jettag = \"AK4PFchs\";");
  // choose a systematic uncertainty in pt average ("central" (pt=120GeV), "up" (pt=240GeV), "down" (pt=60GeV), "doubleup" (pt=480GeV) or "nominal" (pt -> mean value of pt ave))
  gROOT->ProcessLine("TString variation= \"nominal\";");
  // Tag for time dependence plots
  //  gROOT->ProcessLine("TString tag = \"_9\";");
  gROOT->ProcessLine("TString tag = \"\";");


  gROOT->ProcessLine("TFile* datafile = new TFile(path+datfile,\"READ\");"); 
  gROOT->ProcessLine("TFile* MCfile = new TFile(path+mcfile,\"READ\");"); 

  //Load common functions for data treatment
  gROOT->ProcessLine(".L UsefulFunctions.C+");

  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  gROOT->ProcessLine("TString MCname = \"herwigpp\";");
  //  gROOT->ProcessLine("TString MCname = \"pythia8\";");
  gROOT->ProcessLine(".L Control_Plots.C");
  gROOT->ProcessLine("Control_Plots(path,datafile,MCfile,MCname);");

  gROOT->ProcessLine(".L kFSR_pT.C+");
  gROOT->ProcessLine("kFSR_pT(true,path,datafile,MCfile);");
  gROOT->ProcessLine(".L kFSR_pT.C+"); //recomplie to have correct labels on plots
  gROOT->ProcessLine("kFSR_pT(false,path,datafile,MCfile);");

  gROOT->ProcessLine(".L PTextrapolation.C");
  gROOT->ProcessLine("PTextrapolation(true,path,datafile,MCfile,txttag,jettag,variation,tag);");
  gROOT->ProcessLine("PTextrapolation(false,path,datafile,MCfile,txttag,jettag,variation,tag);");

  gROOT->ProcessLine(".L AllResPlots.C");
  gROOT->ProcessLine("AllResPlots(path)");

  //W/t kFSR ---------------------------------
  gROOT->ProcessLine(".L PTextrapolation_COMB_noKFSR.C");
  gROOT->ProcessLine("PTextrapolation_COMB_noKFSR(path,datafile,MCfile,txttag,jettag,variation,tag);");
  
  gROOT->ProcessLine(".L AllResPlots_noKFSR.C");
  gROOT->ProcessLine("AllResPlots_noKFSR(path)");
  //-----------------------------------------------

  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  // // Input for Mikko's Global Fit ----------------------
  // gROOT->ProcessLine(".L InputForGlobalFit.C");
  // gROOT->ProcessLine("InputForGlobalFit(path,datafile,MCfile);");
  // //----------------------------------------------------
}
