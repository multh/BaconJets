#ifndef  HEADER_H
#define  HEADER_H
#include "TString.h"

// define number of bins and bin ranges

const int n_etabarr=5; // needed for the normalization to 1 in the barrel

const int n_alpha = 9;
TString alpha_range[n_alpha] = {"a005", "a010", "a015", "a020", "a025", "a030", "a035", "a040", "a045"};//tmp solution! should be changed to alpha_range everywhere in the code
double alpha_bins[n_alpha] = {0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350,  0.400, 0.450};


const int n_pt = 9;
TString pt_range[n_pt] = {"56.000", "78.000", "100.000", "168.000", "232.000", "300.000", "366.000", "453.000", "562.000"};
double pt_bins[n_pt] = {56, 78, 100, 168, 232, 300, 366, 453, 562};
const int n_eta = 17;
double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 5.191};
TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "5.191"};
TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139",  "3489", "5191"};



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
#endif
