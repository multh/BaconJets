
void run(){

  // change your absolute path here:
  //  gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_noReweight_18EtaBins_AllTriggers/\";");
  //  gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_noReweight/\";");
    gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_noReweight_17EtaBins_AllTriggers/\";");

  //  gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_noReweight_18EtaBins_HFTriggers/\";");
  //  ROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_noReweight_18EtaBins_AllTriggers/\";");

  //gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_PtReweight/\";");
  //  gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_noReweight_alpha0_1_Selection/\";");
  //gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_Reweight69mb/\";");

  // input data files
  gROOT->ProcessLine("TString datfile = \"uhh2.AnalysisModuleRunner.DATA.RunD_AK4CHS.root\";");
  gROOT->ProcessLine("TString mcfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_pythia8_AK4CHS.root\";");
  //  gROOT->ProcessLine("TString mcfile = \"uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_herwigpp_AK4CHS.root\";");
  //  gROOT->ProcessLine("TString mcfile = \"MC_QCD_Pt-15to7000_Flat_Fall15_25nsV1_MC_Reweight_58mb/uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_pythia8_AK4CHS.root\";");
  //gROOT->ProcessLine("TString mcfile = \"MC_QCD_Pt-15to7000_Flat_Fall15_25nsV1_MC_Reweight_69mb/uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_pythia8_AK4CHS.root\";");
  //  gROOT->ProcessLine("TString mcfile = \"MC_QCD_Pt-15to7000_Flat_Fall15_25nsV1_MC_Reweight_80mb/uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_pythia8_AK4CHS.root\";");

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


  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // // gROOT->ProcessLine("TString MCname = \"herwigpp\";");
  // gROOT->ProcessLine("TString MCname = \"pythia8\";");
  // gROOT->ProcessLine(".L Control_Plots.C");
  // gROOT->ProcessLine("Control_Plots(path,datafile,MCfile,MCname);");
  
  // // gROOT->ProcessLine(".L kFSR.C");
  // // gROOT->ProcessLine("kFSR(true,path,datafile,MCfile);");
  // // gROOT->ProcessLine("kFSR(false,path,datafile,MCfile);");

  // // gROOT->ProcessLine(".L kFSR_pT.C");
  // // gROOT->ProcessLine("kFSR_pT(true,path,datafile,MCfile);");
  // // gROOT->ProcessLine("kFSR_pT(false,path,datafile,MCfile);");
  
  // gROOT->ProcessLine(".L PTextrapolation.C");
  // gROOT->ProcessLine("PTextrapolation(true,path,datafile,MCfile,txttag,jettag,variation,tag);");
  // gROOT->ProcessLine("PTextrapolation(false,path,datafile,MCfile,txttag,jettag,variation,tag);");

  // // gROOT->ProcessLine(".L AllResPlots.C");
  // // gROOT->ProcessLine("AllResPlots(path)");

  // // // // // // // gROOT->ProcessLine(".L PTextrapolation_COMB.C");
  // // // // // // //gROOT->ProcessLine("PTextrapolation_COMB(path,datafile,MCfile,txttag,jettag,variation,tag);");

  

  // // // // //COMB -----------------
  gROOT->ProcessLine(".L PTextrapolation_COMB_noKFSR.C");
  gROOT->ProcessLine("PTextrapolation_COMB_noKFSR(path,datafile,MCfile,txttag,jettag,variation,tag);");

  // gROOT->ProcessLine(".L AllResPlots_noKFSR.C");
  // gROOT->ProcessLine("AllResPlots_noKFSR(path)");
  // // // // //COMB[END] -----------------

  // // // // gROOT->ProcessLine(".L ResPlots_addCOMB.C");
  // // // // gROOT->ProcessLine("ResPlots_addCOMB(path)");
 
  
  // //Input for Global Fit by Mikko ----------------
  // gROOT->ProcessLine(".L InputForGlobalFit.C");
  // gROOT->ProcessLine("InputForGlobalFit(path,datafile,MCfile);");
  // //----------------------------------------------

  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

 

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 

  // gROOT->ProcessLine(".L ResPlots.C");
  // gROOT->ProcessLine("ResPlots(path)");

  // gROOT->ProcessLine(".L ResPlots_OldNew.C");
  // gROOT->ProcessLine("ResPlots_OldNew(path)");

  /*
  gROOT->ProcessLine(".L VarPtave.C");
  gROOT->ProcessLine("VarPtave(path);");
  */
  /*
  gROOT->ProcessLine(".L VarPtaveDijet.C");
  gROOT->ProcessLine("VarPtaveDijet(path);");
  */
  /*
  gROOT->ProcessLine(".L alpha.C");
  gROOT->ProcessLine("alpha(path);");
  */
  /*
  gROOT->ProcessLine(".L Rrel.C");
  gROOT->ProcessLine("Rrel(datafile,MCfile);");
  */
  /*
  gROOT->ProcessLine(".L timedep.C");
  gROOT->ProcessLine("timedep(path);");
  */


}
