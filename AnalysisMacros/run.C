//Set of macros, which reads TTree from UHH2/BaconJet output files

//include "tdrstyle_mod15.C"



void run(TString Generator, TString Collection, TString Runs, bool only_AllCollections, bool WithJetID){

  gROOT->ProcessLine(" lumi_binned = false;");
  gROOT->ProcessLine(" with_sys = false;");
  gROOT->ProcessLine(" is_sys_MC = true;"); //Systematic only relevant for MC samples, no systematic data files need to be read in

  //Jet Collection
  if(Collection == "AK4CHS"){
    gROOT->ProcessLine(" s_jets = \"AK4CHS\";");
    gROOT->ProcessLine(" jettag = \"AK4PFchs\";");
  }
  else if(Collection == "AK8CHS"){
    gROOT->ProcessLine(" s_jets = \"AK8CHS\";");
    gROOT->ProcessLine(" jettag = \"AK8PFchs\";");
  }
  else if(Collection == "AK4Puppi"){
    gROOT->ProcessLine(" s_jets = \"AK4Puppi\";");
    gROOT->ProcessLine(" jettag = \"AK4PFpuppi\";");
  }
  else if(Collection == "AK8Puppi"){
    gROOT->ProcessLine(" s_jets = \"AK8Puppi\";");
    gROOT->ProcessLine(" jettag = \"AK8PFpuppi\";");
  }
  else throw runtime_error("run.C: Invalid Jet-collection specified.");


  if(Collection == "AK4CHS" || Collection == "AK4Puppi" ) gROOT->ProcessLine(" s_postfix = \"\";");
  else gROOT->ProcessLine(" s_postfix = \"_NoJERSmearing\";");

  //Runs
  if(Runs == "BCD")  gROOT->ProcessLine(" runnr = \"BCD\";");
  else if(Runs == "EF")  gROOT->ProcessLine(" runnr = \"EF\";");
  else if(Runs == "G")  gROOT->ProcessLine(" runnr = \"G\";");
  else if(Runs == "B")  gROOT->ProcessLine(" runnr = \"B\";");
  else if(Runs == "C")  gROOT->ProcessLine(" runnr = \"C\";");
  else if(Runs == "D")  gROOT->ProcessLine(" runnr = \"D\";");
  else if(Runs == "E")  gROOT->ProcessLine(" runnr = \"E\";");
  else if(Runs == "F")  gROOT->ProcessLine(" runnr = \"F\";");
  else throw runtime_error("run.C: Invalid RunNr. specified.");

  // choose a tag for the txt output files
  //Generator
  if(Generator == "pythia") {gROOT->ProcessLine(" txttag = \"pythia8_v4\";"); gROOT->ProcessLine(" Generator = \"pythia\";");}
  else if (Generator == "herwig") {gROOT->ProcessLine(" txttag = \"herwigpp_v4\";"); gROOT->ProcessLine(" Generator = \"herwig\";");}
  else throw runtime_error("run.C: Invalid Generator specified.");

  
  //Runs
  if(Runs == "BCD")  gROOT->ProcessLine(" lumitag=\"RunBCD  12.9 fb^{-1}\";");
  else if(Runs == "EF")  gROOT->ProcessLine(" lumitag=\"RunEF  7.3 fb^{-1}\";");
  else if(Runs == "G")  gROOT->ProcessLine(" lumitag=\"RunG  1.8 fb^{-1}\";");
  else if(Runs == "B")  gROOT->ProcessLine(" lumitag=\"RunB  5.8 fb^{-1}\";");
  else if(Runs == "C")  gROOT->ProcessLine(" lumitag=\"RunC  2.6 fb^{-1}\";");
  else if(Runs == "D")  gROOT->ProcessLine(" lumitag=\"RunD  4.3 fb^{-1}\";");
  else if(Runs == "E")  gROOT->ProcessLine(" lumitag=\"RunE  4.1 fb^{-1}\";");
  else if(Runs == "F")  gROOT->ProcessLine(" lumitag=\"RunF  3.2 fb^{-1}\";");
  else throw runtime_error("run.C: Invalid RunNr. specified.");

  //WithJetID?
  if(WithJetID) {
    gROOT->ProcessLine("s_JetID = \"_Loose_JetPFID\";");
    gROOT->ProcessLine("b_JetID = true;");
  }
  else{
    gROOT->ProcessLine("s_JetID = \"\";");
    gROOT->ProcessLine("b_JetID = false;");
  }


  // change your absolute path here:
  gROOT->ProcessLine(" dir =\"RunBCDEFG_noRes\" ");

  gROOT->ProcessLine(" path = \"/nfs/dust/cms/user/reimersa/JEC/JEC_80X_standalone/\" + dir + \"/\" + s_jets + s_JetID + s_postfix + \"/Run\" + runnr +\"/\";"); 
  gROOT->ProcessLine(" path_up = \"/nfs/dust/cms/user/reimersa/JEC/JEC_80X_standalone/\" + dir + \"_Sys_JER_Up/\" + s_jets + s_JetID + s_postfix + \"/Run\" + runnr +\"/\";");
  gROOT->ProcessLine(" path_down = \"/nfs/dust/cms/user/reimersa/JEC/JEC_80X_standalone/\" + dir + \"_Sys_JER_Down/\" + s_jets + s_JetID + s_postfix + \"/Run\" + runnr +\"/\";");
  gROOT->ProcessLine(" path_general = \"/nfs/dust/cms/user/reimersa/JEC/JEC_80X_standalone/\" + dir +\"/\";"); //for general plots in the very end
 
  // input files

  gROOT->ProcessLine(" realdatfile = \"../uhh2.AnalysisModuleRunner.DATA.DATA_Run\"+runnr+\"_\"+s_jets+\".root\";");
  gROOT->ProcessLine(" datfile = \"../uhh2.AnalysisModuleRunner.DATA.DATA_Run\"+runnr+\"_\"+s_jets+\".root\";");
  if(Generator == "herwig") gROOT->ProcessLine(" mcfile = \"../uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_herwigpp_\"+s_jets+\".root\";");
  else if(Generator == "pythia") gROOT->ProcessLine(" mcfile = \"../uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_\"+s_jets+\".root\";");
  
 

  //only for HIP mitigation plots
  gROOT->ProcessLine(" lumis[0] = \"0p6\";");
  gROOT->ProcessLine(" lumis[1] = \"0p8\";");
  gROOT->ProcessLine(" lumis[2] = \"1p2\";");
  if(Collection == "AK4CHS") {
    cout << "AK4CHS collection specified." << endl;
    gROOT->ProcessLine(" path_HIP = \"/nfs/dust/cms/user/reimersa/JEC/JEC_80X_standalone/\" + dir + \"/\" + s_jets + s_JetID + s_postfix +\"/\";"); 
  }
  else if(Collection == "AK4Puppi") {
    cout << "AK4Puppi collection specified." << endl; 
    gROOT->ProcessLine(" path_HIP = \"/nfs/dust/cms/user/reimersa/JEC/JEC_80X_standalone/\" + dir + \"/\" + s_jets + s_JetID + s_postfix +\"/\";"); 
  }
  else { 
    cout << "Collection different from AK4CHS or AK4Puppi specified! HIP_path will be generic " << endl; 
    gROOT->ProcessLine(" path_HIP = path");
  }
  gROOT->ProcessLine(" HIP = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt-15to7000_Flat_HIP_\";");
  gROOT->ProcessLine(" HIP_mtOff = \"uhh2.AnalysisModuleRunner.MC.QCD_Pt-15to7000_Flat_HIP_mtOff_\";");
  gROOT->ProcessLine(" DATA_HIP = \"uhh2.AnalysisModuleRunner.DATA.DATA_Run2016E_29Jul_HIP\";");
  gROOT->ProcessLine(" DATA_mtOff = \"uhh2.AnalysisModuleRunner.DATA.DATA_RunBCD_\" + s_jets;");








  // choose a systematic uncertainty in pt average ("central" (pt=120GeV), "up" (pt=240GeV), "down" (pt=60GeV), "doubleup" (pt=480GeV) or "nominal" (pt -> mean value of pt ave))

  // Tag for time dependence plots
  //gROOT->ProcessLine("TString tag = \"_9\";");
  //gROOT->ProcessLine("TString tag = \"_pythia8\";");
  gROOT->ProcessLine(" tag = \"\";");

  gROOT->ProcessLine("Realdatafile = new TFile(path+realdatfile,\"READ\");"); 
  gROOT->ProcessLine(" datafile = new TFile(path+datfile,\"READ\");"); 
  gROOT->ProcessLine(" MCfile = new TFile(path+mcfile,\"READ\");"); 
  gROOT->ProcessLine("if(with_sys) {MCfile_up = new TFile(path_up+mcfile,\"READ\"); MCfile_down = new TFile(path_down+mcfile,\"READ\");}");
  gROOT->ProcessLine("if(with_sys && is_sys_MC) {datafile_up = datafile; datafile_down = datafile; Realdatafile_up = Realdatafile; Realdatafile_down = Realdatafile;}");
  gROOT->ProcessLine("if(with_sys && !is_sys_MC){datafile_up = new TFile(path_up+datfile,\"READ\"); datafile_down = new TFile(path_down+datfile,\"READ\"); Realdatafile_up = new TFile(path_up+realdatfile,\"READ\"); Realdatafile_down = new TFile(path_down+realdatfile,\"READ\");}");



  gROOT->ProcessLine("al_cut=0.3;"); 
  gROOT->ProcessLine("nResponseBins = 100;");



  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

  
  if(Generator == "herwig") gROOT->ProcessLine(" MCname = \"herwigpp\";");
  else if (Generator == "pythia") gROOT->ProcessLine(" MCname = \"pythia8\";");
  else throw runtime_error("run.C: Invalid Generator specified.");

  gROOT->ProcessLine("cout << \"Opening MC file: \" << path+mcfile << endl;");
  


  // ********************************************* starting analysis *********************************** //
  
  //gROOT->ProcessLine("HIP_Mitigation_Plots(path_HIP,HIP,HIP_mtOff,DATA_HIP,DATA_mtOff,s_jets,lumis, al_cut);");
  
  if(!only_AllCollections){
    /*
      gROOT->ProcessLine("Control_Plots(lumi_binned,path,datafile,MCfile,MCname,Realdatafile);");
      
    
    gROOT->ProcessLine("kFSR_pT_TTree(true,false,1,path,datafile,MCfile,txttag,lumitag,al_cut,nResponseBins);");
    gROOT->ProcessLine("if(lumi_binned){for(int i=4; i<n_lumi;i++){kFSR_pT_TTree(true,lumi_binned,i,path,datafile,MCfile,txttag,lumitag,al_cut,nResponseBins);}}");
    gROOT->ProcessLine("if(with_sys){kFSR_pT_TTree(true,false,1,path_up,datafile_up,MCfile_up,txttag,lumitag,al_cut,nResponseBins); kFSR_pT_TTree(true,false,1,path_down,datafile_down,MCfile_down,txttag,lumitag,al_cut,nResponseBins);}");
      
 
    
    gROOT->ProcessLine("PTextrapolation_TTree_kFSRfit(true,false,1,path,datafile,MCfile,txttag,jettag,variation,tag,Generator,runnr,lumitag,al_cut,nResponseBins);");
    gROOT->ProcessLine("PTextrapolation_TTree_kFSRfit(false,false,1,path,datafile,MCfile,txttag,jettag,variation,tag,Generator,runnr,lumitag,al_cut,nResponseBins);");//WATCH OUT for kFSR fit results here!!!
    gROOT->ProcessLine("if(lumi_binned){for(int i=4; i<n_lumi;i++){PTextrapolation_TTree_kFSRfit(true,lumi_binned,i,path,datafile,MCfile,txttag,jettag,variation,tag,Generator,runnr,lumitag,al_cut,nResponseBins);}}");
    gROOT->ProcessLine("if(lumi_binned){for(int i=4; i<n_lumi;i++){PTextrapolation_TTree_kFSRfit(false,lumi_binned,i,path,datafile,MCfile,txttag,jettag,variation,tag,Generator,runnr,lumitag,al_cut,nResponseBins);}}");
    gROOT->ProcessLine("if(with_sys){PTextrapolation_TTree_kFSRfit(true,false,1,path_up,datafile_up,MCfile_up,txttag,jettag,variation,tag,Generator,runnr,lumitag,al_cut,nResponseBins); PTextrapolation_TTree_kFSRfit(true,false,1,path_down,datafile_down,MCfile_down,txttag,jettag,variation,tag,Generator,runnr,lumitag,al_cut,nResponseBins);}");
    gROOT->ProcessLine("if(with_sys){PTextrapolation_TTree_kFSRfit(false,false,1,path_up,datafile_up,MCfile_up,txttag,jettag,variation,tag,Generator,runnr,lumitag,al_cut,nResponseBins); PTextrapolation_TTree_kFSRfit(false,false,1,path_down,datafile_down,MCfile_down,txttag,jettag,variation,tag,Generator,runnr,lumitag,al_cut,nResponseBins);}");//WATCH OUT for kFSR fit results here!!!
    
 
                      
    gROOT->ProcessLine("L2ResOutput(false,with_sys,1,path,path_up,path_down,txttag,jettag,tag,lumitag,al_cut)");
    gROOT->ProcessLine("if(lumi_binned){for(int i=4; i<n_lumi;i++){L2ResOutput(lumi_binned,false,i,path,path_up,path_down,txttag,jettag,tag,lumitag,al_cut);}}");
    */
    /*
    if(Generator == "pythia"){
      gROOT->ProcessLine(".L FlavorCorrection_TTree.C");                                          
      gROOT->ProcessLine("FlavorCorrection_TTree(path,MCfile,jettag, txttag,al_cut)");
    }
    else cout << "Flavor plots not created for herwig samples." << endl;
    */
    if(WithJetID) gROOT->ProcessLine("JetID_ComparisonPlots(path, mcfile, datfile, txttag, s_jets, runnr)");
  }
  else{
    gROOT->ProcessLine(".L L2ResALLJetColl.C");  
    gROOT->ProcessLine("L2ResALLJetColl(path_general,txttag,tag,lumitag,runnr,b_JetID)");
  }

  //gROOT->ProcessLine(".L L2ResTxtTest.C");                                          
  //gROOT->ProcessLine("L2ResTxtTest(path,txttag,jettag,tag,al_cut)"); 
  //gROOT->ProcessLine(".L AllResPlots.C");
  //gROOT->ProcessLine("AllResPlots(path,txttag,jettag,variation,tag,al_cut)");

}
