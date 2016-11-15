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
  //Objects.emplace_back(CorrectionObject("G", "pythia", "AK4CHS"));
  Objects.emplace_back(CorrectionObject("BCD", "pythia", "AK4CHS", true));
  cout << "testobject is " << Objects[0] << endl;

  Objects[0].ControlPlots();
  Objects[0].kFSR();
  Objects[0].Pt_Extrapolation(true);
  Objects[0].Pt_Extrapolation(false);
  Objects[0].L2ResOutput();

  cout << endl << "Closing MC and DATA files." << endl;
  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CloseFiles();
  cout << "Going to return 0 now, cya." << endl;
  return 0;
}
