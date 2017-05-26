#include <cmath>
#include <iostream>
#include "../include/CorrectionObject.h"
#include "../include/parameters.h"
#include <TString.h>


using namespace std;



int main(){
  cout << "Hello from main(). What am I going to do?" << endl << endl;


  TString generator    = "pythia";
  bool    closure_test    = false;
  bool    trigger_fwd     = true;    //Use for Weight Calc
  bool    trigger_central = false;     //Use for Weight Calc
  TString collection    = "AK4CHS";

    TString input_path   = "/nfs/dust/cms/user/multh/JEC/2016ReReco/Residuals/Summer16_03Feb2017_V3/AK4CHS/MC_Reweighted_chsMET_ForWeights_Down/FWD/";
    TString weight_path  = "/nfs/dust/cms/user/multh/JEC/2016ReReco/Residuals/Summer16_03Feb2017_V3/AK4CHS/MC_Reweighted_chsMET_ForWeights_Down/";


  //  TString input_path   = "/nfs/dust/cms/user/multh/JEC/2016ReReco/Residuals/Summer16_03Feb2017_V0/AK4CHS/MC_Reweighted_chsMET/";
  // TString weight_path  = "/nfs/dust/cms/user/multh/JEC/2016ReReco/Residuals/Summer16_03Feb2017_V0/AK4CHS/MC_Reweighted_chsMET_ForWeights/";

  // TString input_path   = "/nfs/dust/cms/user/multh/JEC/2016Legacy/Residuals/Summer16_18Apr_EGReg_chsMET/AK4CHS/MC_Reweighted_chsMET/";
  // TString weight_path  = "/nfs/dust/cms/user/multh/JEC/2016Legacy/Residuals/Summer16_18Apr_EGReg_chsMET/AK4CHS/MC_Reweighted_chsMET_ForWeights/";  //Check Path for Weights 


   //eine Klasse: enthaelt Info ueber runnr, Generator, collection, Strings zu MC/DATA-files, memberfunctions: controlPlots, kFSR etc.
  vector<CorrectionObject> Objects;

   Objects.emplace_back(CorrectionObject("BCD", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
   Objects.emplace_back(CorrectionObject("EFearly", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
   Objects.emplace_back(CorrectionObject("FlateG", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
   Objects.emplace_back(CorrectionObject("H", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
  
  //Objects.emplace_back(CorrectionObject("BCDEFGH", generator,collection, input_path, weight_path, closure_test));

  cout << "testobject is " << Objects[0] << endl;

  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CalculateMCWeights();

if(trigger_central && !trigger_fwd){
  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CalculateMCWeights_TriggerThresholds(true);  //Central Triggers
}
else if(!trigger_central && trigger_fwd){
  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CalculateMCWeights_TriggerThresholds(false);  //FWD Triggers 
}
else {
	cout<<"No Weight Calculation"<<endl;
}
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].ControlPlots();

  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].kFSR();
// for(unsigned int i=0; i<Objects.size(); i++) Objects[i].kFSR_CorrectFormulae();
// for(unsigned int i=0; i<Objects.size(); i++) Objects[i].kFSR_CorrectFormulae_eta();

  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation(true);
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation(false);
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative(true);
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative(false);
  //  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative_CorrectFormulae(true);
  //  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative_CorrectFormulae(false);
 
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative_CorrectFormulae_eta(true);
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative_CorrectFormulae_eta(false);

// for(unsigned int i=0; i<Objects.size(); i++) Objects[i].L2ResOutput();
// for(unsigned int i=0; i<Objects.size(); i++) Objects[i].L2ResOutput_eta();
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].InputForGlobalFit();
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].InputForGlobalFit_eta_0_13();

  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FinalControlPlots();
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FinalControlPlots_CorrectFormulae();
// for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FinalControlPlots_CorrectFormulae_eta();

//for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FullCycle_CorrectFormulae();
  
  //  Objects[0].L2ResAllRuns();
  //  Objects[0].L2ResOverlay(true);
  //  Objects[0].L2ResOverlay(false);
  //  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FinalControlPlots_CorrectFormulae();

// for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FullCycle_CorrectFormulae();
  
//  Objects[0].L2ResAllRuns();
//  Objects[0].L2ResOverlay();

  cout << endl << "Closing MC and DATA files." << endl;
  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CloseFiles();
  cout << "Going to return 0 now, cya." << endl << endl;
  return 0;
}
