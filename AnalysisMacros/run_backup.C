//Set of macros, which reads TTree from UHH2/BaconJet output files

using namespace std;

void run(TString Generator, TString Collection, TString Runs){

  gROOT->ProcessLine("bool lumi_binned = false;");
  gROOT->ProcessLine("bool with_sys = false;");
  gROOT->ProcessLine("bool is_sys_MC = true;"); //Systematic only relevant for MC samples, no systematic data files need to be read in

  gROOT->ProcessLine("TString s_jets = \"AK4CHS\";");
  gROOT->ProcessLine("TString jettag = \"AK4PFchs\";");
  gROOT->ProcessLine("TString s_postfix = \"\";");
  //gROOT->ProcessLine("TString s_postfix = \"_NoJERSmearing\";");
  if(Runs == "BCD")  gROOT->ProcessLine("TString runnr = \"BCD\";");
  else if(Runs == "EF")  gROOT->ProcessLine("TString runnr = \"EF\";");
  else if(Runs == "G")  gROOT->ProcessLine("TString runnr = \"G\";");
  else if(Runs == "B")  gROOT->ProcessLine("TString runnr = \"B\";");
  else if(Runs == "C")  gROOT->ProcessLine("TString runnr = \"C\";");
  else if(Runs == "D")  gROOT->ProcessLine("TString runnr = \"D\";");
  else if(Runs == "E")  gROOT->ProcessLine("TString runnr = \"E\";");
  else if(Runs == "F")  gROOT->ProcessLine("TString runnr = \"F\";");
  else throw runtime_error("run.C: Invalid RunNr. specified.");

  // choose a tag for the txt output files
  
  if(Generator == "pythia")gROOT->ProcessLine("TString txttag = \"pythia8_v4\";");
  else if (Generator == "herwig")gROOT->ProcessLine("TString txttag = \"herwigpp_v4\";");
  else throw runtime_error("run.C: Invalid Generator specified.");





  // change your absolute path here:
  gROOT->ProcessLine("TString dir =\"RunBCDEFG_noRes\" ");

   gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/reimersa/JEC/JEC_80X_standalone/\" + dir + \"/\" + s_jets + s_postfix + \"/Run\" + runnr +\"/\";");
  gROOT->ProcessLine("TString path_up = \"/nfs/dust/cms/user/reimersa/JEC/JEC_80X_standalone/\" + dir + \"_Sys_JER_Up/\" + s_jets + s_postfix + \"/Run\" + runnr +\"/\";");
  gROOT->ProcessLine("TString path_down = \"/nfs/dust/cms/user/reimersa/JEC/JEC_80X_standalone/\" + dir + \"_Sys_JER_Down/\" + s_jets + s_postfix + \"/Run\" + runnr +\"/\";");
  //gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/reimersa/JEC/JEC_80X_standalone/FlavorCorrections/\";");
  //gROOT->ProcessLine("TString path = \"/nfs/dust/cms/user/reimersa/JEC/JEC_80X_standalone/test_1/\";");

  // input data files

  gROOT->ProcessLine("TString realdatfile = \"../uhh2.AnalysisModuleRunner.DATA.DATA_Run\"+runnr+\"_\"+s_jets+\".root\";");
  gROOT->ProcessLine("TString datfile = \"../uhh2.AnalysisModuleRunner.DATA.DATA_Run\"+runnr+\"_\"+s_jets+\".root\";");
  //gROOT->ProcessLine("TString realdatfile = \"uhh2.AnalysisModuleRunner.DATA.DATA_RunB_\"+s_jets+\".root\";");
  //gROOT->ProcessLine("TString datfile = \"uhh2.AnalysisModuleRunner.DATA.DATA_RunB_\"+s_jets+\".root\";");
  gROOT->ProcessLine("TString mcfile = \"../uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_herwigpp_\"+s_jets+\".root\";");
  
  //For HIP mitigation studies --> also change datafile if needed.
  //gROOT->ProcessLine("TString mcfile = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt-15to7000_Flat_HIP_mtOff_0p6.root\";");
 

  //only for HIP mitigation plots
  gROOT->ProcessLine("TString lumis[3] = {\"0p6\",\"0p8\",\"1p2\"};");
  gROOT->ProcessLine("TString HIP = \"../uhh2.AnalysisModuleRunner.MC.QCD_Pt-15to7000_Flat_HIP_\";");
  gROOT->ProcessLine("TString HIP_mtOff = \"../uhh2.AnalysisModuleRunner.MC.QCD_Pt-15to7000_Flat_HIP_mtOff_\";");
  gROOT->ProcessLine("TString DATA_HIP = \"../uhh2.AnalysisModuleRunner.DATA.DATA_Run2016E_29Jul_HIP\";");
  gROOT->ProcessLine("TString DATA_mtOff = \"../uhh2.AnalysisModuleRunner.DATA.DATA_RunBCD_\" + s_jets;");








  // choose a systematic uncertainty in pt average ("central" (pt=120GeV), "up" (pt=240GeV), "down" (pt=60GeV), "doubleup" (pt=480GeV) or "nominal" (pt -> mean value of pt ave))

  // Tag for time dependence plots
  //gROOT->ProcessLine("TString tag = \"_9\";");
  //gROOT->ProcessLine("TString tag = \"_pythia8\";");
  gROOT->ProcessLine("TString tag = \"\";");

  gROOT->ProcessLine("TFile* Realdatafile = new TFile(path+realdatfile,\"READ\");"); 
  gROOT->ProcessLine("TFile* datafile = new TFile(path+datfile,\"READ\");"); 
  gROOT->ProcessLine("TFile* MCfile = new TFile(path+mcfile,\"READ\");"); 
  gROOT->ProcessLine("TFile* MCfile_up, *MCfile_down, *datafile_up, *datafile_down, *Realdatafile_up, *Realdatafile_down;");
  gROOT->ProcessLine("if(with_sys) {MCfile_up = new TFile(path_up+mcfile,\"READ\"); MCfile_down = new TFile(path_down+mcfile,\"READ\");}");
  gROOT->ProcessLine("if(with_sys && is_sys_MC) {datafile_up = datafile; datafile_down = datafile; Realdatafile_up = Realdatafile; Realdatafile_down = Realdatafile;}");
  gROOT->ProcessLine("if(with_sys && !is_sys_MC){datafile_up = new TFile(path_up+datfile,\"READ\"); datafile_down = new TFile(path_down+datfile,\"READ\"); Realdatafile_up = new TFile(path_up+realdatfile,\"READ\"); Realdatafile_down = new TFile(path_down+realdatfile,\"READ\");}");



  gROOT->ProcessLine("double al_cut=0.3;"); 
  gROOT->ProcessLine("int nResponseBins = 100;");

 


  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  /*
  gROOT->ProcessLine(".L HIP_Mitigation_Plots.C");
  gROOT->ProcessLine("HIP_Mitigation_Plots(path,HIP,HIP_mtOff,DATA_HIP,DATA_mtOff,lumis, al_cut);");
  */
  
  /*
  gROOT->ProcessLine("TString MCname = \"herwigpp\";");
  //gROOT->ProcessLine("TString MCname = \"pythia8\";");
  gROOT->ProcessLine(".L Control_Plots.C");
  gROOT->ProcessLine("Control_Plots(lumi_binned,path,datafile,MCfile,MCname,Realdatafile);");
  
  
  
  
  gROOT->ProcessLine(".L kFSR_pT_TTree.C+");
  gROOT->ProcessLine("kFSR_pT_TTree(true,false,1,path,datafile,MCfile,txttag,al_cut,nResponseBins);");
  gROOT->ProcessLine("if(lumi_binned){for(int i=4; i<n_lumi;i++){kFSR_pT_TTree(true,lumi_binned,i,path,datafile,MCfile,txttag,al_cut,nResponseBins);}}");
  gROOT->ProcessLine("if(with_sys){kFSR_pT_TTree(true,false,1,path_up,datafile_up,MCfile_up,txttag,al_cut,nResponseBins); kFSR_pT_TTree(true,false,1,path_down,datafile_down,MCfile_down,txttag,al_cut,nResponseBins);}");

 
  
  
  gROOT->ProcessLine("TString variation= \"nominal\";");
  gROOT->ProcessLine(".L PTextrapolation_TTree_kFSRfit.C+");
  gROOT->ProcessLine("PTextrapolation_TTree_kFSRfit(true,false,1,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);");
  gROOT->ProcessLine("PTextrapolation_TTree_kFSRfit(false,false,1,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);");//WATCH OUT for kFSR fit results here!!!
  gROOT->ProcessLine("if(lumi_binned){for(int i=4; i<n_lumi;i++){PTextrapolation_TTree_kFSRfit(true,lumi_binned,i,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);}}");
  gROOT->ProcessLine("if(lumi_binned){for(int i=4; i<n_lumi;i++){PTextrapolation_TTree_kFSRfit(false,lumi_binned,i,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);}}");
  gROOT->ProcessLine("if(with_sys){PTextrapolation_TTree_kFSRfit(true,false,1,path_up,datafile_up,MCfile_up,txttag,jettag,variation,tag,al_cut,nResponseBins); PTextrapolation_TTree_kFSRfit(true,false,1,path_down,datafile_down,MCfile_down,txttag,jettag,variation,tag,al_cut,nResponseBins);}");
  gROOT->ProcessLine("if(with_sys){PTextrapolation_TTree_kFSRfit(false,false,1,path_up,datafile_up,MCfile_up,txttag,jettag,variation,tag,al_cut,nResponseBins); PTextrapolation_TTree_kFSRfit(false,false,1,path_down,datafile_down,MCfile_down,txttag,jettag,variation,tag,al_cut,nResponseBins);}");//WATCH OUT for kFSR fit results here!!!
  
  
  gROOT->ProcessLine("variation= \"central\";");
  gROOT->ProcessLine(".L PTextrapolation_TTree_kFSRfit.C+");
  gROOT->ProcessLine("PTextrapolation_TTree_kFSRfit(true,false,1,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);");
  gROOT->ProcessLine("PTextrapolation_TTree_kFSRfit(false,false,1,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);");//WATCH OUT for kFSR fit results here!!!
  gROOT->ProcessLine("if(lumi_binned){for(int i=4; i<n_lumi;i++){PTextrapolation_TTree_kFSRfit(true,lumi_binned,i,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);}}");
  gROOT->ProcessLine("if(lumi_binned){for(int i=4; i<n_lumi;i++){PTextrapolation_TTree_kFSRfit(false,lumi_binned,i,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);}}");

  gROOT->ProcessLine("variation= \"up\";");
  gROOT->ProcessLine(".L PTextrapolation_TTree_kFSRfit.C+");
  gROOT->ProcessLine("PTextrapolation_TTree_kFSRfit(true,false,1,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);");
  gROOT->ProcessLine("PTextrapolation_TTree_kFSRfit(false,false,1,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);");//WATCH OUT for kFSR fit results here!!!
  gROOT->ProcessLine("if(lumi_binned){for(int i=4; i<n_lumi;i++){PTextrapolation_TTree_kFSRfit(true,lumi_binned,i,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);}}");
  gROOT->ProcessLine("if(lumi_binned){for(int i=4; i<n_lumi;i++){PTextrapolation_TTree_kFSRfit(false,lumi_binned,i,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);}}");

  gROOT->ProcessLine("variation= \"doubleup\";");
  gROOT->ProcessLine(".L PTextrapolation_TTree_kFSRfit.C+");
  gROOT->ProcessLine("PTextrapolation_TTree_kFSRfit(true,false,1,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);");
  gROOT->ProcessLine("PTextrapolation_TTree_kFSRfit(false,false,1,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);");//WATCH OUT for kFSR fit results here!!!
  gROOT->ProcessLine("if(lumi_binned){for(int i=4; i<n_lumi;i++){PTextrapolation_TTree_kFSRfit(true,lumi_binned,i,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);}}");
  gROOT->ProcessLine("if(lumi_binned){for(int i=4; i<n_lumi;i++){PTextrapolation_TTree_kFSRfit(false,lumi_binned,i,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);}}");

  gROOT->ProcessLine("variation= \"down\";"); 
  gROOT->ProcessLine(".L PTextrapolation_TTree_kFSRfit.C+");
  gROOT->ProcessLine("PTextrapolation_TTree_kFSRfit(true,false,1,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);");
  gROOT->ProcessLine("PTextrapolation_TTree_kFSRfit(false,false,1,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);");//WATCH OUT for kFSR fit results here!!!
  gROOT->ProcessLine("if(lumi_binned){for(int i=4; i<n_lumi;i++){PTextrapolation_TTree_kFSRfit(true,lumi_binned,i,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);}}");
  gROOT->ProcessLine("if(lumi_binned){for(int i=4; i<n_lumi;i++){PTextrapolation_TTree_kFSRfit(false,lumi_binned,i,path,datafile,MCfile,txttag,jettag,variation,tag,al_cut,nResponseBins);}}");
  
  
  gROOT->ProcessLine(".L L2ResOutput.C");                                          
  gROOT->ProcessLine("L2ResOutput(false,with_sys,1,path,path_up,path_down,txttag,jettag,tag,al_cut)");
  gROOT->ProcessLine("if(lumi_binned){for(int i=4; i<n_lumi;i++){L2ResOutput(lumi_binned,false,i,path,path_up,path_down,txttag,jettag,tag,al_cut);}}");
  */
  /*
  gROOT->ProcessLine(".L FlavorCorrection_TTree.C");                                          
  gROOT->ProcessLine("FlavorCorrection_TTree(path,MCfile,jettag, txttag,al_cut)");
  */

  //gROOT->ProcessLine(".L L2ResTxtTest.C");                                          
  //gROOT->ProcessLine("L2ResTxtTest(path,txttag,jettag,tag,al_cut)"); 
  //gROOT->ProcessLine(".L L2ResALLJetColl.C");  
  //gROOT->ProcessLine("L2ResALLJetColl(path,txttag,jettag,tag)");

  //gROOT->ProcessLine(".L AllResPlots.C");
  //gROOT->ProcessLine("AllResPlots(path,txttag,jettag,variation,tag,al_cut)");

  //gROOT->ProcessLine(".q");
}
