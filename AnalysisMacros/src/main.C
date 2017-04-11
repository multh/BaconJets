//#includes
#include <cmath>
#include <iostream>
#include "../include/CorrectionObject.h"
#include "../include/parameters.h"
#include <TString.h>


using namespace std;



int main(){
  cout << "Hello from main(). What am I going to do?" << endl << endl;


  TString generator      = "pythia";
  bool closure_test      = true;      //false
  TString collection     = "AK4CHS";
  
  
  TString input_path  = "/nfs/dust/cms/user/multh/JEC/2016ReReco/ClosureTest/Summer16_23Sep2016_V6/" + collection + "/MC_Reweighted_PU_pt_hat_cut_ForWeights/CENTRAL/";
  TString weight_path = "/nfs/dust/cms/user/multh/JEC/2016ReReco/ClosureTest/Summer16_23Sep2016_V6/" + collection + "/MC_Reweighted_PU_pt_hat_cut_ForWeights/";  //Check Path for Weights 

 
   //eine Klasse: enthaelt Info ueber runnr, Generator, collection, Strings zu MC/DATA-files, memberfunctions: controlPlots, kFSR etc.
  vector<CorrectionObject> Objects;

  //Objects
  Objects.emplace_back(CorrectionObject("BCD", generator,collection, input_path, weight_path, closure_test));
  Objects.emplace_back(CorrectionObject("EFearly", generator,collection, input_path, weight_path, closure_test));
  //Objects.emplace_back(CorrectionObject("FlateG", generator, "AK4CHS",version, closure_test));
  //Objects.emplace_back(CorrectionObject("H", generator, collection, input_path, weight_path, closure_test));
  //Objects.emplace_back(CorrectionObject("BCDEFGH", generator, "AK4CHS",version, closure_test));

  //Objects for Closure Test
 
 
  //Objects.emplace_back(CorrectionObject("BCD", generator,collection, input_path, weight_path, closure_test));
  //Objects.emplace_back(CorrectionObject("EFearly", generator,collection, input_path, weight_path, closure_test));
  //Objects.emplace_back(CorrectionObject("FlateG", generator,collection, input_path, weight_path, closure_test));
  //Objects.emplace_back(CorrectionObject("H", generator, collection, input_path, weight_path, closure_test));
  //Objects.emplace_back(CorrectionObject("BCDEFGH", generator,collection, input_path, weight_path, closure_test));


  cout << "testobject is " << Objects[0] << endl;

  //Reweight Samples 
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CalculateMCWeights();
  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CalculateMCWeights_TriggerThresholds(true);    //Central Triggers
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CalculateMCWeights_TriggerThresholds(false);   //FWD Triggers 


  // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].ControlPlots();

  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].kFSR();
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].kFSR_CorrectFormulae();

  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation(true);
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation(false);
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative(true);
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative(false);
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative_CorrectFormulae(true);
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative_CorrectFormulae(false);

  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].L2ResOutput();
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].InputForGlobalFit();
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].InputForGlobalFit_eta_0_13();

  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FinalControlPlots();
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FinalControlPlots_CorrectFormulae();

  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FullCycle_CorrectFormulae();
  
  //Objects[0].L2ResAllRuns();
  //Objects[0].L2ResOverlay(true);
  //Objects[0].L2ResOverlay(false);





  cout << endl << "Closing MC and DATA files." << endl;
  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CloseFiles();
  cout << "Going to return 0 now, cya." << endl << endl;
  return 0;
}
