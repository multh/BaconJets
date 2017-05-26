#pragma once

#include <TString.h>

using namespace std;


const double alpha_cut = 0.3;
const TString s_alpha_cut = "0.3";



const int n_alpha = 9;
const int n_eta = 19;

//Negative Eta Range
//const int n_eta = 37;
const int n_pt = 11;
const int nResponseBins = 100;// 100
const int n_etabarr=5; // needed for the normalization to 1 in the barrel



const TString alpha_range[n_alpha] = {"a005", "a010", "a015", "a020", "a025", "a030", "a035", "a040", "a045"};
const TString pt_range[n_pt]= {"51", "73", "95", "129", "163", "230", "299", "365", "453", "566", "1000"};

const TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "3.839", "5.191"};
const TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "3489", "3839", "5191"};

//Negative Eta Range
//const TString eta_range[n_eta] = {"-5.191","-3.839","-3.489","-3.139","-2.964","-2.853", "-2.65", "-2.5", "-2.322", "-2.172", "-1.93", "-1.653", "-1.479", "-1.305", "-1.044", "-0.783", "-0.522", "-0.261"," 0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "3.839", "5.191"};

//const TString eta_range2[n_eta] = {"-5191","-3839","-3489","-3139","-2964","-2853", "-265", "-25", "-2322", "-2172", "-193", "-1653", "-1479", "-1305", "-1044", "-0783", "-0522", "-0261","00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "3489", "3839", "5191"};



const double alpha_bins[n_alpha] = {0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350,  0.400, 0.450};
const double pt_bins[n_pt]       = {51, 73, 95, 129, 163, 230, 299, 365, 453, 566, 1000};
const double eta_bins[n_eta]     = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};

//Negative Eta Range
//const double eta_bins[n_eta]     = {-5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -1.93, -1.653, -1.479, -1.305, -1.044, -0.783, -0.522, -0.261, 0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};



const int n_alpha_common = 4;
const TString alpha_range_common[n_alpha_common] = {"a10","a15", "a20", "a30"};//for Global fit we produce result only in few alpha values
const double alpha_bins_common[n_alpha_common] = {0.100, 0.150, 0.200, 0.300};

const int n_eta_common = 8;
const double eta_common_bins[n_eta_common] ={0, 0.783, 1.305, 1.93, 2.5, 2.964, 3.2, 5.191};
const TString eta_common_range[n_eta_common] = {"0.000", "0.783", "1.305", "1.930", "2.500", "2.964", "3.200", "5.191"};
const TString eta_output[n_eta_common-1] = {"eta00-08", "eta08-13","eta13-19","eta19-25","eta25-30","eta30-32","eta32-52"};

const int n_eta_common_2 = 7;
const double eta_common_bins_2[n_eta_common] ={0, 1.305, 1.93, 2.5, 2.964, 3.2, 5.191};
const TString eta_common_range_2[n_eta_common] = {"0.000", "1.305", "1.930", "2.500", "2.964", "3.200", "5.191"};
const TString eta_output_2[n_eta_common-1] = {"eta00-13","eta13-19","eta19-25","eta25-30","eta30-32","eta32-52"};
