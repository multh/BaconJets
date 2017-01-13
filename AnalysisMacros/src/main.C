//#includes
#include <cmath>
#include <iostream>
#include "../include/CorrectionObject.h"
#include "../include/parameters.h"
#include <TString.h>


using namespace std;



int main(){
  cout << "Hello from main(). What am I going to do?" << endl << endl;

   //eine Klasse: enthaelt Info ueber runnr, Generator, collection, Strings zu MC/DATA-files, memberfunctions: controlPlots, kFSR etc.
  vector<CorrectionObject> Objects;
  //Closure tests
  //Objects.emplace_back(CorrectionObject("BCD", "pythia", "AK4CHS", true, true));
  //Objects.emplace_back(CorrectionObject("EFearly", "pythia", "AK4CHS", true, true));
  //Objects.emplace_back(CorrectionObject("FlateG", "pythia", "AK4CHS", true, true));
  //Objects.emplace_back(CorrectionObject("H", "pythia", "AK4CHS", true, true));

  //New residuals
  Objects.emplace_back(CorrectionObject("BCD", "pythia", "AK4CHS", true, false));
  Objects.emplace_back(CorrectionObject("EFearly", "pythia", "AK4CHS", true, false));
  Objects.emplace_back(CorrectionObject("FlateG", "pythia", "AK4CHS", true, false));
  Objects.emplace_back(CorrectionObject("H", "pythia", "AK4CHS", true, false));
  cout << "testobject is " << Objects[0] << endl;


  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].ControlPlots();
  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].kFSR();
  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation(true);
  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation(false);
  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].L2ResOutput();
  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].InputForGlobalFit();
  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FinalControlPlots();
  Objects[0].L2ResAllRuns();
  


  cout << endl << "Closing MC and DATA files." << endl;
  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CloseFiles();
  cout << "Going to return 0 now, cya." << endl << endl;
  return 0;
}
