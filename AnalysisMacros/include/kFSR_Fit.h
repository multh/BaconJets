#pragma once

#include <iostream>
#include <utility>
#include <cmath>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>

TF1* kFSR_Fit(TGraphErrors* mean, TGraphErrors* width, int i, int j);

