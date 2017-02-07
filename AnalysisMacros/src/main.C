//#includes
#include <cmath>
#include <iostream>
#include "../include/CorrectionObject.h"
#include "../include/parameters.h"
#include <TString.h>


using namespace std;



int main(){
  cout << "Hello from main(). What am I going to do?" << endl << endl;


  TString generator = "pythia";
  bool closure_test = false;


   //eine Klasse: enthaelt Info ueber runnr, Generator, collection, Strings zu MC/DATA-files, memberfunctions: controlPlots, kFSR etc.
  vector<CorrectionObject> Objects;

  //Objects.emplace_back(CorrectionObject("BCD", generator, "AK4CHS", true, closure_test));
  Objects.emplace_back(CorrectionObject("EFearly", generator, "AK4CHS", true, closure_test));
  //Objects.emplace_back(CorrectionObject("FlateG", generator, "AK4CHS", true, closure_test));
  //Objects.emplace_back(CorrectionObject("H", generator, "AK4CHS", true, closure_test));
  //Objects.emplace_back(CorrectionObject("BCDEFGH", generator, "AK4CHS", true, closure_test));

  cout << "testobject is " << Objects[0] << endl;

  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CalculateMCWeights();


  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].ControlPlots();

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
  
  //Objects[0].L2ResAllRuns();
  //Objects[0].L2ResOverlay();





  cout << endl << "Closing MC and DATA files." << endl;
  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CloseFiles();
  cout << "Going to return 0 now, cya." << endl << endl;
  return 0;
}
