#ifndef  HEADER_H
#define  HEADER_H

// define number of bins
const int n_alpha = 5;
const int n_eta = 16;
const int n_pt = 9;

const int n_etabarr=5; // needed for the normalization to 1 in the barrel

// define the bin ranges
TString alpha_range[n_alpha-1] = {"a01", "a02", "a03", "a04"};
double alpha_bins[n_alpha] = {0.000, 0.100, 0.200, 0.300, 0.400};

TString pt_range[n_pt] = {"55.000", "76.000", "93.000", "172.000", "232.000", "300.000", "366.000", "452.000", "558.000"};
double pt_bins[n_pt] = {55, 76, 93, 172, 232, 300, 366, 452, 558};

double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.191};
TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "5.191"};
TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139",  "5191"};
TString eta_range3[n_eta] = {"0.0", "0.261", "0.522", "0.783", "1.044", "1.305", "1.653", "1.93", "2.172", "2.322", "2.5", "2.65", "2.853", "2.964", "3.139", "5.191"};
 
