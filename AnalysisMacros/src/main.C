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
  //Objects.emplace_back(CorrectionObject("BCD", "pythia", "AK4CHS", true));
  //Objects.emplace_back(CorrectionObject("E", "pythia", "AK4CHS", true));
  Objects.emplace_back(CorrectionObject("Fearly", "pythia", "AK4CHS", true));
  Objects.emplace_back(CorrectionObject("FlateG", "pythia", "AK4CHS", true));
  cout << "testobject is " << Objects[0] << endl;

  for(unsigned int i=0; i<Objects.size(); i++){
    Objects[i].ControlPlots();
    Objects[i].kFSR();
    Objects[i].Pt_Extrapolation(true);
    Objects[i].Pt_Extrapolation(false);
    Objects[i].L2ResOutput();
  }


  cout << endl << "Closing MC and DATA files." << endl;
  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CloseFiles();
  cout << "Going to return 0 now, cya." << endl << endl;
  return 0;
}
