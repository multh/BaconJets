#include "header.h"
#include "TGraph.h"
#include <iostream>
#include <sstream>
#include <cmath>
#include "TString.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

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
  TString selection = "alpha<";
  selection +=alpha;
  selection += " && probejet_eta<";
  selection +=eta_up;
  selection += " && probejet_eta>=";
  selection +=eta_low;
  selection += " && pt_ave<";
  selection +=pT_up;
  selection +=" && pt_ave>=";
  selection +=pT_low;
  TH1D* hist = GetHist(rootfile, selection, varName, 200, 0, 2.5);
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

//Calculate ratio between MC and DATA
pair<double,double> Rmc_to_Rdata(pair<double,double> mc, pair<double,double> data){
  pair<double,double> out;
  out.first =0;   out.second =0;
  if(abs(data.first)>1e-3){ 
  double ratio = mc.first/data.first;
  double ratioError = sqrt(pow(mc.second,2)+pow(data.second,2));
  out.first =ratio;   out.second =ratioError;
  }
  return out;
}

//Clean points not filled due to low statistic
TGraphErrors* CleanEmptyPoints(TGraphErrors* input){
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
