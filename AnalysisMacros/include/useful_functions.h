#pragma once

#include <iostream>
#include <utility>
#include <cmath>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TFile.h>

std::pair<double,double> Rmc_to_Rdata(std::pair<double,double> mc, std::pair<double,double> data);

std::pair<double,double> GetValueAndError(TH1D *hin);

TGraphErrors* CleanEmptyPoints(TGraphErrors* input);

TH1D* GetHist(TFile *rootfile, TString selection, TString varName, int nbins, double low, double up);
