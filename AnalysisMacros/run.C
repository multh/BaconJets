//Set of macros, which reads TTree from UHH2/BaconJet output files

void run(){

  // change your absolute path here:
  //  gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV2_noReweight_AllTriggers/\";");
  //gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV2_noReweight_AllTriggers_assym02/\";");
  //gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV2_noReweight_AllTriggers_assym02_nJets_le25/\";");
  //  gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV2_L2L3Res_noReweight/\";");
  //  gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_noReweight_AllTriggers_TTree/\";");
  //  gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_ReweightPU_AllTriggers_TTree/\";");
  //  gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_ReweightPU_HFTriggers_TTree/\";");
  //  gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_ReweightPU_StandTriggers_TTree/\";");
  //gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV2_noReweight/\";");

  //  gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_noReweight/\";");
  //  gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_noReweight/3dJet_50GeV/\";");
  //  gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV2_L2L3Res_noReweight/\";");

  //  gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_80X/Spring16_25nsV1/noHFtriggers_with_HLT_DiPFJetAve80_140/\";");
  gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/karavdia/JEC_80X/Spring16_25nsV1/\";");


  // input data files
  gROOT->ProcessLine("TString realdatfile = \"uhh2.AnalysisModuleRunner.DATA.RunB_AK4CHS.root\";");//TEST
  //  gROOT->ProcessLine("TString realdatfile = \"uhh2.AnalysisModuleRunner.DATA.RunB_AK4PUPPI.root\";");
  //  gROOT->ProcessLine("TString realdatfile = \"uhh2.AnalysisModuleRunner.DATA.RunB_AK8CHS.root\";");
  //  gROOT->ProcessLine("TString realdatfile = \"uhh2.AnalysisModuleRunner.DATA.RunB_AK8PUPPI.root\";");

  gROOT->ProcessLine("TString datfile = \"uhh2.AnalysisModuleRunner.DATA.RunB_AK4CHS.root\";");
  //  gROOT->ProcessLine("TString datfile = \"uhh2.AnalysisModuleRunner.DATA.RunB_AK4PUPPI.root\";");
  //  gROOT->ProcessLine("TString datfile = \"uhh2.AnalysisModuleRunner.DATA.RunB_AK8CHS.root\";");
  //  gROOT->ProcessLine("TString datfile = \"uhh2.AnalysisModuleRunner.DATA.RunB_AK8PUPPI.root\";");


  //  gROOT->ProcessLine("TString datfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_pythia8_AK4CHS.root\";"); 
  //  gROOT->ProcessLine("TString datfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_pythia8_AK4PUPPI.root\";");
  //  gROOT->ProcessLine("TString datfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_pythia8_AK8CHS.root\";");
  //  gROOT->ProcessLine("TString datfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_pythia8_AK8PUPPI.root\";");

  //
  // gROOT->ProcessLine("TString mcfile = \"uhh2.AnalysisModuleRunner.DATA.RunD_AK4CHS.root\";");
  // gROOT->ProcessLine("TString datfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_pythia8_AK4CHS.root\";");

  //  gROOT->ProcessLine("TString mcfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_herwigpp_AK4CHS.root\";");
  //  gROOT->ProcessLine("TString mcfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_herwigpp_AK4PUPPI.root\";");
//  gROOT->ProcessLine("TString mcfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_herwigpp_AK8CHS.root\";");
//  gROOT->ProcessLine("TString mcfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_herwigpp_AK8PUPPI.root\";");

  //  gROOT->ProcessLine("TString mcfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_pythia8_AK8PUPPI.root\";");
//    gROOT->ProcessLine("TString mcfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_pythia8_AK8CHS.root\";");
//gROOT->ProcessLine("TString mcfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_pythia8_AK4PUPPI.root\";");
gROOT->ProcessLine("TString mcfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_pythia8_AK4CHS.root\";");
//  
  // gROOT->ProcessLine("TString datfile = \"uhh2.AnalysisModuleRunner.DATA.RunD_AK8CHS.root\";");
    //  gROOT->ProcessLine("TString mcfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_pythia8_AK8CHS.root\";");
  //  
  
  //  gROOT->ProcessLine("TString mcfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Sum_Pt_pythia8_AK4CHS.root\";");

  // choose a tag for the txt output files
  gROOT->ProcessLine("TString txttag = \"v2\";");
  // tag for the jet cone size
  gROOT->ProcessLine("TString jettag = \"AK4PFchs\";");
  //  gROOT->ProcessLine("TString jettag = \"AK4PFpuppi\";");
  //  gROOT->ProcessLine("TString jettag = \"AK8PFchs\";");
  //  gROOT->ProcessLine("TString jettag = \"AK8PFpuppi\";");
  // choose a systematic uncertainty in pt average ("central" (pt=120GeV), "up" (pt=240GeV), "down" (pt=60GeV), "doubleup" (pt=480GeV) or "nominal" (pt -> mean value of pt ave))
  //  gROOT->ProcessLine("TString variation= \"nominal\";");
  // Tag for time dependence plots
  //  gROOT->ProcessLine("TString tag = \"_9\";");
  gROOT->ProcessLine("TString tag = \"\";");

  gROOT->ProcessLine("TFile* Realdatafile = new TFile(path+realdatfile,\"READ\");"); 

  gROOT->ProcessLine("TFile* datafile = new TFile(path+datfile,\"READ\");"); 
  gROOT->ProcessLine("TFile* MCfile = new TFile(path+mcfile,\"READ\");"); 

  gROOT->ProcessLine("int nResponseBins = 100;");
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  //  gROOT->ProcessLine("TString MCname = \"herwigpp\";");
  gROOT->ProcessLine("TString MCname = \"pythia8\";");
  gROOT->ProcessLine(".L Control_Plots.C");
  gROOT->ProcessLine("Control_Plots(path,datafile,MCfile,MCname,Realdatafile);");

  gROOT->ProcessLine("double al_cut=0.3;"); 
  // // // // // // // // // // // // // // // // // // // //  gROOT->ProcessLine("double al_cut=0.7;"); 
  gROOT->ProcessLine(".L kFSR_pT_TTree.C+");
  gROOT->ProcessLine("kFSR_pT_TTree(true,path,datafile,MCfile,al_cut,nResponseBins);");
  // // // // // // // // gROOT->ProcessLine("kFSR_pT_TTree(false,path,datafile,MCfile,al_cut,nResponseBins);");
  // // // // // // // // // // gROOT->ProcessLine(".L kFSR_TTree.C+");
  // // // // // // // // // // gROOT->ProcessLine("kFSR_TTree(true,path,datafile,MCfile,al_cut,nResponseBins);");
  // // // // // // // // // // gROOT->ProcessLine("kFSR_TTree(false,path,datafile,MCfile,al_cut,nResponseBins);");

  gROOT->ProcessLine("TString variation= \"nominal\";");
  // // //   gROOT->ProcessLine("TString variation= \"central\";");
  // // //  gROOT->ProcessLine("TString variation= \"up\";");
  // // //  gROOT->ProcessLine("TString variation= \"down\";");
  // // // gROOT->ProcessLine("TString variation= \"doubleup\";");
  gROOT->ProcessLine(".L PTextrapolation_TTree.C+");
  gROOT->ProcessLine("PTextrapolation_TTree(true,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);");
  gROOT->ProcessLine("PTextrapolation_TTree(false,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);");
  gROOT->ProcessLine(".L AllResPlots.C");
  gROOT->ProcessLine("AllResPlots(path,al_cut,variation)");
  
  // //  gROOT->ProcessLine("TString variation= \"central\";");
  // //  gROOT->ProcessLine("TString variation= \"up\";");
  // //  gROOT->ProcessLine("TString variation= \"down\";");
  // // gROOT->ProcessLine("TString variation= \"doubleup\";");
  // // gROOT->ProcessLine(".L PTextrapolation_TTree.C+");
  // // gROOT->ProcessLine("PTextrapolation_TTree(true,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut);");
  // // gROOT->ProcessLine("PTextrapolation_TTree(false,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut);"); 
  

  // // // // // // // // // //W/t kFSR ---------------------------------
  // gROOT->ProcessLine("TString variation= \"nominal\";");
  // // // //  gROOT->ProcessLine("TString variation= \"central\";");
  // // // //gROOT->ProcessLine("TString variation= \"up\";"); 
  // // // //gROOT->ProcessLine("TString variation= \"down\";");
  // // // // gROOT->ProcessLine("TString variation= \"doubleup\";");
  // gROOT->ProcessLine(".L PTextrapolation_COMB_noKFSR_TTree.C");
  // gROOT->ProcessLine("PTextrapolation_COMB_noKFSR_TTree(path,datafile,MCfile,txttag,jettag,variation,tag,al_cut);");
  // gROOT->ProcessLine(".L AllResPlots_noKFSR.C");
  // gROOT->ProcessLine("AllResPlots_noKFSR(path,al_cut)");
  // // // // gROOT->ProcessLine(".L ResPtDependence.C");
  // // // // gROOT->ProcessLine("TString method= \"COMB\";");
  // // // // // //  gROOT->ProcessLine("TString method= \"MPF\";");
  // // // // gROOT->ProcessLine("ResPtDependence(path,al_cut,method)");
  // // // // // // // // //-----------------------------------------------

  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  // // Input for Mikko's Global Fit ----------------------
  // gROOT->ProcessLine(".L InputForGlobalFit_TTree.C+");
  // gROOT->ProcessLine("InputForGlobalFit_TTree(path,datafile,MCfile);");
 //  //----------------------------------------------------
 //  gROOT->ProcessLine(".L PTextrapolation.C");
 // // choose a systematic uncertainty in pt average ("central" (pt=120GeV), "up" (pt=240GeV), "down" (pt=60GeV), "doubleup" (pt=480GeV) or "nominal" (pt -> mean value of pt ave))
 //  gROOT->ProcessLine("TString variation= \"central\";");
 //  gROOT->ProcessLine("PTextrapolation(true,path,datafile,MCfile,txttag,jettag,variation,tag);");
 //  gROOT->ProcessLine("PTextrapolation(false,path,datafile,MCfile,txttag,jettag,variation,tag);");
 //  gROOT->ProcessLine("TString variation= \"up\";");
 //  gROOT->ProcessLine("PTextrapolation(true,path,datafile,MCfile,txttag,jettag,variation,tag);");
 //  gROOT->ProcessLine("PTextrapolation(false,path,datafile,MCfile,txttag,jettag,variation,tag);");
 //  gROOT->ProcessLine("TString variation= \"down\";");
 //  gROOT->ProcessLine("PTextrapolation(true,path,datafile,MCfile,txttag,jettag,variation,tag);");
 //  gROOT->ProcessLine("PTextrapolation(false,path,datafile,MCfile,txttag,jettag,variation,tag);");
 //  gROOT->ProcessLine("TString variation= \"doubleup\";");
 //  gROOT->ProcessLine("PTextrapolation(true,path,datafile,MCfile,txttag,jettag,variation,tag);");
 //  gROOT->ProcessLine("PTextrapolation(false,path,datafile,MCfile,txttag,jettag,variation,tag);");


 // gROOT->ProcessLine(".L kFSR_pT.C+");
  // gROOT->ProcessLine("kFSR_pT(true,path,datafile,MCfile);");
  // gROOT->ProcessLine(".L kFSR_pT.C+"); //recomplie to have correct labels on plots
  // gROOT->ProcessLine("kFSR_pT(false,path,datafile,MCfile);");
  // gROOT->ProcessLine(".L PTextrapolation.C");
  // gROOT->ProcessLine("PTextrapolation(true,path,datafile,MCfile,txttag,jettag,variation,tag);");
  // gROOT->ProcessLine("PTextrapolation(false,path,datafile,MCfile,txttag,jettag,variation,tag);");
// // gROOT->ProcessLine(".L PTextrapolation_COMB_noKFSR.C");
  // // gROOT->ProcessLine("PTextrapolation_COMB_noKFSR(path,datafile,MCfile,txttag,jettag,variation,tag);");
}
