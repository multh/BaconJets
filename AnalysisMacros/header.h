#ifndef  HEADER_H
#define  HEADER_H
#include "TString.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include <iostream>
#include <sstream>
#include <cmath>
#include "TString.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

// define number of bins and bin ranges
const int n_etabarr=5; // needed for the normalization to 1 in the barrel
//bins for Fall15_V2!
const int n_alpha = 9;
TString alpha_range[n_alpha] = {"a005", "a010", "a015", "a020", "a025", "a030", "a035", "a040", "a045"};//tmp solution! should be changed to alpha_range everywhere in the code
double alpha_bins[n_alpha] = {0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350,  0.400, 0.450};
const int n_pt = 9;
TString pt_range[n_pt] = {"56.000", "78.000", "100.000", "168.000", "232.000", "300.000", "366.000", "453.000", "562.000"};
double pt_bins[n_pt] = {56, 78,  100, 168, 232, 300, 366, 453, 562};




/* AK4CHS */
/* /\* 2015 *\/ */
const int n_eta = 19;
double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "3.839", "5.191"};
TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "3489", "3839", "5191"};

/* const int n_eta = 17; */
/* double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.191}; */
/* TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "5.191"}; */
/* TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "5191"}; */


/* //AK4PUPPI */
/* const int n_eta = 19; */
/* double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191}; */
/* TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "3.839", "5.191"}; */
/* TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "3489", "3839", "5191"}; */


/* //AK8CHS */
/* const int n_eta = 19; */
/* double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191}; */
/* TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "3.839", "5.191"}; */
/* TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139",  "3489", "3839", "5191"}; */

/* /\* //AK8PUPPI *\/ */
/* const int n_eta = 19; */
/* double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 4.013, 5.191}; */
/* TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "4.013", "5.191"}; */
/* TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139",  "3489", "4013", "5191"}; */


/* /\* /\\* /\\\* /\\\\* /\\\\\* //bins for Fall15_V1! *\\\\\/ *\\\\/ *\\\/ *\\/ *\/ */
/* const int n_eta = 17; */
/* double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.191}; */
/* TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "5.191"}; */
/* TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "5191"}; */

/* const int n_eta = 19; */
/* double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.664, 5.191}; */
/* TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.314", "3.664", "5.191"}; */
/* TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "3314", "3664", "5191"}; */

/* const int n_eta = 19; */
/* double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191}; */
/* TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "3.839", "5.191"}; */
/* TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "3489", "3839", "5191"}; */


/* const int n_eta = 18; */
/* double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 4.363, 5.191}; */
/* TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "4.363", "5.191"}; */
/* TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139",  "4363", "5191"}; */


/* /\* /\\* /\\\* //bins for Fall15_V2! *\\\/ *\\/ *\/ */
/* const int n_eta = 18; */
/* double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 5.191}; */
/* TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "5.191"}; */
/* TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139",  "3489", "5191"}; */


/* const int n_eta = 19; */
/* double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.664, 5.191}; */
/* TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "3.664", "5.191"}; */
/* TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139",  "3489", "3664", "5191"}; */


/* const int n_eta = 20; */
/* double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139,  3.314, 3.489, 3.664,  5.191}; */
/* TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.314", "3.489", "3.664","5.191"}; */
/* TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "3314", "3489", "3664", "5191"}; */

/* const int n_eta = 21; */
/* double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139,  3.314, 3.489, 3.664, 3.839, 5.191}; */
/* TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.314", "3.489", "3.664", "3.839", "5.191"}; */
/* TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "3314", "3489", "3664",  "3839", "5191"}; */

/* const int n_eta = 25; */
/* double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139,  3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 5.191}; */
/* TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.314", "3.489", "3.664", "3.839", "4.013", "4.191", "4.363", "4.538", "5.191"}; */
/* TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "3314", "3489", "3664",  "3839", "4013", "4191", "4363", "4538", "5191"}; */

/* const int n_alpha = 7; */
/* TString alpha_range[n_alpha] = {"a040", "a050", "a060","a070", "a080", "a090", "a100"};//tmp solution! should be changed to alpha_range everywhere in the code */
/* double alpha_bins[n_alpha] = {0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000 }; */
/* const int n_eta = 16; */
/* double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 3.139, 5.191}; */
/* TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "3.139", "5.191"}; */
/* TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "3139", "5191"}; */

/* const int n_pt = 11; */
/* TString pt_range[n_pt] = {"56.000", "78.000","89.000", "100.000", "168.000", "232.000", "266.000", "300.000", "366.000", "453.000", "562.000"}; */
/* double pt_bins[n_pt] = {56, 78, 89, 100, 168, 232, 266, 300, 366, 453, 562}; */

/* const int n_pt = 18; */
/* TString pt_range[n_pt] = {"56.000","67.000", "78.000","83.000", "89.000", "100.000","117.000","134.000", "168.000","176.000","184.000","192.000", "200.000", "232.000", "300.000", "366.000", "453.000", "562.000"}; */
/* double pt_bins[n_pt] = {56, 67, 78, 83, 89, 100, 117, 134, 168, 176, 184, 192, 200, 232, 300, 366, 453, 562}; */



/* const int n_eta = 21; */
/* double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 4.013,  5.191}; */
/* TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.314", "3.489", "3.664", "4.013", "5.191"}; */
/* TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "3314", "3489", "3664", "4013",  "5191"}; */


/* const int n_eta = 22; */
/* double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 4.013, 4.363, 5.191}; */
/* TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.314", "3.489", "3.664", "4.013", "4.363", "5.191"}; */
/* TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "3314", "3489", "3664", "4013", "4363", "5191"}; */
/* const int n_eta = 12; */
/* double eta_bins[n_eta] = {3.139, 3.314,3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191}; */
/* TString eta_range[n_eta] = {"3.139", "3.489", "3.664", "3.839", "4.013", "4.191", "4.363", "4.538", "4.716", "4.889","5.191"}; */
/* TString eta_range2[n_eta] = {"3139",  "3489", "3664", "3839", "4013", "4191", "4363", "4538", "4716", "4889", "5191"}; */

/* const int n_eta = 16; */
/* double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.83, 2.043, 2.322, 2.65, 2.853, 3.139, 3.489,  5.191}; */
/* TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.830", "2.043", "2.322", "2.650", "2.853", "3.139", "3.489", "5.191"}; */
/* TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "1830", "2043", "2322", "2650", "2853", "3139",  "3489",  "5191"}; */



/* //bin used for test of macros --------------- */
/* const int n_alpha = 4; */
/* TString alpha_range[n_alpha] = {"a015", "a020", "a030", "a040"}; */
/* double alpha_bins[n_alpha] = {0.150, 0.200, 0.300, 0.400}; */
/* const int n_pt = 5; */
/* TString pt_range[n_pt] = {"56.000", "100.000", "232.000", "366.000", "562.000"}; */
/* double pt_bins[n_pt] = {56, 100, 232, 366, 562}; */
/* //const int n_eta = 5; */
/* /\* double eta_bins[n_eta] = {0, 1.305, 2.5,  3.139, 5.191}; *\/ */
/* /\* TString eta_range[n_eta] = {"0.000",  "1.305", "2.500", "3.139","5.191"}; *\/ */
/* /\* TString eta_range2[n_eta] = {"00", "1305", "25","3139","5191"}; *\/ */
/* const int n_eta = 3; */
/* double eta_bins[n_eta] = {2.5,  3.139, 5.191}; */
/* TString eta_range[n_eta] = {"2.500", "3.139","5.191"}; */
/* TString eta_range2[n_eta] = {"25","3139","5191"}; */
/* //TString eta_range3[n_eta] = {"0.0", "1.305", "2.5", "3.139", "5.191"}; */
/* //-------------------------------------------- */

//bins used to produce input for global fit
const int n_alpha_common = 4;
TString alpha_range_common[n_alpha_common] = {"a010","a015", "a020", "a030"};//for Global fit we produce result only in few alpha values
double alpha_bins_common[n_alpha_common] = {0.100, 0.150, 0.200, 0.300};

const int n_eta_common = 8;
double eta_common_bins[n_eta_common] ={0, 0.783, 1.305, 1.93, 2.5, 2.964, 3.2, 5.191};
TString eta_common_range[n_eta_common] = {"0.000", "0.783", "1.305", "1.930", "2.500", "2.964", "3.200", "5.191"};
TString eta_output[n_eta_common-1] = {"eta0000-0783", "eta0783-1305","eta1305-1930","eta1930-2500","eta2500-2964","eta2964-3200","eta3200-5191"};



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Some useful functions --------------------------------------------------
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//Calculate ratio between MC and DATA
pair<double,double> Rmc_to_Rdata(pair<double,double> mc, pair<double,double> data){
  pair<double,double> out;
  out.first =0;   out.second =0;
  if(abs(data.first)>1e-3){ 
    double ratio = mc.first/data.first;
    double ratioError = sqrt(pow(mc.second,2)+pow(data.second,2));
    /* double ratio = mc.first;//TEST */
    /* double ratioError = sqrt(pow(mc.second,2)); //TEST */
  out.first =ratio;   out.second =ratioError;
  }
  return out;
}

//Clean points not filled due to low statistic
TGraphErrors* CleanEmptyPoints(TGraphErrors* input){
  //  input->Print();
  double *Yval = input->GetY();
  double *YvalError = input->GetEY();
  double *Xval = input->GetX();
  double *XvalError = input->GetEX();
  int count=0;
  vector<double> Xnew,Ynew,Xerrornew,Yerrornew;
  for(int i=0;i<input->GetN();i++){
    if(Yval[i]!=0){
      count++;
      Xnew.push_back(Xval[i]);       
      Ynew.push_back(Yval[i]);
      Xerrornew.push_back(XvalError[i]);       
      Yerrornew.push_back(YvalError[i]);
    }
  }

  const int NnewSize =  count;
  double Xnew_m[NnewSize],Ynew_m[NnewSize],Xerrornew_m[NnewSize],Yerrornew_m[NnewSize]; //because silly ROOT doesn't know how to treat vectors
  for(int i=0;i<NnewSize;i++){
    Xnew_m[i] = Xnew[i];     Ynew_m[i] = Ynew[i];
    Xerrornew_m[i] = Xerrornew[i];     Yerrornew_m[i] = Yerrornew[i];
  }

  TGraphErrors* output = new TGraphErrors(count,Xnew_m,Ynew_m,Xerrornew_m,Yerrornew_m);
  if(input->GetN()!=output->GetN()) cout<<"Number of points in input: "<<input->GetN()<<" in output: "<<output->GetN()<<endl;
  return output;
}

TString ToString(int num) {
  ostringstream start;
  start<<num;
  TString start1=start.str();
  return start1;
}

//Get hist for variable from tree with particular selection 
//NB: Name of the tree should be <OutputTree Name=""/> parameter in UHH2/BaconJet xml file! ("Events" in current version)
TH1D* GetHist(TFile *rootfile, TString selection, TString varName, int nbins, double low, double up){
  TH1D* hist = new TH1D("hist","",nbins,low,up);
  TTree *tree = (TTree*)rootfile->Get("Events");
  int Nev = tree->Project("hist",varName,selection);
  return hist;
}


//Get response values for MPF or pt-balance methods for particular alpha cut value in eta and pT ranges
// returns mean and meanError of the response
//NB: if number of events below NevMin returns 0!
pair<double,double> Response(TFile *rootfile,double alpha, double eta_low, double eta_up, double pT_low, double pT_up, bool isMPF, int NevMin=100){
  pair<double,double> out;
  out.first = 0; out.second = 0;
  TString varName = "rel_r";
  if(isMPF) varName = "mpf_r";
  TString selection = "weight*(alpha<";
  selection +=alpha;
  selection += " && fabs(probejet_eta)<";
  selection +=eta_up;
  selection += " && fabs(probejet_eta)>=";
  selection +=eta_low;
  selection += " && pt_ave<";
  selection +=pT_up;
  selection +=" && pt_ave>=";
  selection +=pT_low;
  selection +=")";
  TH1D* hist = GetHist(rootfile, selection, varName, 100, 0, 2.5);
  int Nev = hist->GetEntries();
  if(Nev>NevMin){
  std::cout<<"Number of events = "<<Nev<<std::endl;
  double mean = hist->GetMean();
  double meanError = hist->GetMeanError();
  out.first = mean;
  out.second = meanError;
  }
  delete hist;
  return out;
}

pair<double,double> GetValueAndError(TH1D *hin){
  pair<double,double> res;
  res.first = 0; res.second = 0;
  //if(hin->GetEntries()>50){
  if(hin->GetEntries()>100){
  //  if(hin->GetEntries()>200){
  //if(hin->GetEntries()>400){
    res.first = hin->GetMean();
    res.second = hin->GetMeanError();
    //res.second = hin->GetRMS()/sqrt(hin->Integral());//TEST: exctly the same as GetMeanError()!
    //    cout<<"Entries = "<<hin->GetEntries()<<" MeanError() = "<<hin->GetMeanError()<<endl;
  }
  return res;
}

#endif
