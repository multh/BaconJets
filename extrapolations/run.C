
void run(){

  // change your absolute path here:
  gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Summer15_25nsV6_nataliia/\";");
  // input data files
  gROOT->ProcessLine("TString datfile = \"uhh2.AnalysisModuleRunner.DATA.DATA_data_V6_cmssw76X_test2.root\";");
  gROOT->ProcessLine("TString mcfile = \"uhh2.AnalysisModuleRunner.MC.MC_cmssw76x_MC_Flat_no_weight.root\";");
  // choose a tag for the txt output files
  gROOT->ProcessLine("TString txttag = \"v7UHH2\";");
  // tag for the jet cone size
  gROOT->ProcessLine("TString jettag = \"AK4PFchs\";");
  // choose a systematic uncertainty in pt average ("central" (pt=120GeV), "up" (pt=240GeV), "down" (pt=60GeV), "doubleup" (pt=480GeV) or "nominal" (pt -> mean value of pt ave))
  gROOT->ProcessLine("TString variation= \"nominal\";");
  // Tag for time dependence plots
  //  gROOT->ProcessLine("TString tag = \"_9\";");
  gROOT->ProcessLine("TString tag = \"\";");


  gROOT->ProcessLine("TFile* datafile = new TFile(path+datfile,\"READ\");"); 
  gROOT->ProcessLine("TFile* MCfile = new TFile(path+mcfile,\"READ\");"); 


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  
  // gROOT->ProcessLine(".L kFSR.C");
  // gROOT->ProcessLine("kFSR(true,path,datafile,MCfile);");
  // gROOT->ProcessLine("kFSR(false,path,datafile,MCfile);");
  
  
  // gROOT->ProcessLine(".L PTextrapolation.C");
  // gROOT->ProcessLine("PTextrapolation(true,path,datafile,MCfile,txttag,jettag,variation,tag);");
  // gROOT->ProcessLine("PTextrapolation(false,path,datafile,MCfile,txttag,jettag,variation,tag);");
  
  
  // gROOT->ProcessLine(".L ResPlots.C");
  // gROOT->ProcessLine("ResPlots(path)");

  gROOT->ProcessLine(".L ResPlots_OldNew.C");
  gROOT->ProcessLine("ResPlots_OldNew(path)");

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
