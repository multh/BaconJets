#include "../include/useful_functions.h"

using namespace std;

//Calculate ratio between MC and DATA
pair<double,double> Rmc_to_Rdata(pair<double,double> mc, pair<double,double> data){
  pair<double,double> out;
  out.first =0;   out.second =0;
  if(abs(data.first)>1e-3){ 
    double ratio = mc.first/data.first;
    double ratioError = sqrt(pow(mc.second/data.first,2) + pow(mc.first*data.second/data.first/data.first,2));

    out.first =ratio;   out.second =ratioError;
  }
  return out;
}

//Get mean value and RMS of a TH1D
pair<double,double> GetValueAndError(TH1D *hin){
  pair<double,double> res;
  res.first = 0; res.second = 0;
  if(hin->GetEntries()>30){
    res.first = hin->GetMean();
    // GetMeanError calculates the uncertainty on the mean value, arising due to limited statistics in the sample. We dont care for the width itself, only the uncertainty on the predicted mean is relevant.
    res.second = hin->GetMeanError();
    /*
    //Using median
    hin->ComputeIntegral();
    double median = 0;
    double q_med = 0.5;
    double q_res_low = 0.16;
    double q_res_high = 0.84;
    double res_low = 0;
    double res_high = 0;
    double res_symm = 0;
    hin->GetQuantiles(1, &median, &q_med);    
    hin->GetQuantiles(1, &res_low, &q_res_low);
    hin->GetQuantiles(1, &res_high, &q_res_high);
    res_symm = ((res_high - median) + (median - res_low))/2 / sqrt(hin->GetEntries());
    //Using gauss fit to core
    TF1 *f1 = new TF1("f1", "gaus", 0.75, 1.25);
    hin->Fit("f1", "R");
    double gaus_mean = f1->GetParameter(1);
    double gaus_sigma = f1->GetParameter(2) / sqrt(hin->GetEntries());
    cout << "Mean: " << res.first << ", error: " << res.second << ". " << endl;
    cout << "If using median -- value: " << median  << " + "<< res_high-median  << " - " << median - res_low << ", symm error: " << res_symm << endl;
    cout << "If using gauss  -- value: " << gaus_mean  << " +- "<< gaus_sigma << endl;
    delete f1;
    */
  }
  return res;
}


//Clean points not filled due to low statistics
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
    Xnew_m[i] = Xnew[i];
    Ynew_m[i] = Ynew[i];
    Xerrornew_m[i] = Xerrornew[i];
    Yerrornew_m[i] = Yerrornew[i];
  }

  TGraphErrors* output = new TGraphErrors(count,Xnew_m,Ynew_m,Xerrornew_m,Yerrornew_m);
  if(input->GetN()!=output->GetN()) cout<<"Number of points in input: "<<input->GetN()<<" in output: "<<output->GetN()<<endl;
  return output;
}



//Get hist for variable from tree with particular selection 
//NB: Name of the tree should be <OutputTree Name=""/> parameter in UHH2/BaconJet xml file! ("AnalysisTree" in current version)
TH1D* GetHist(TFile *rootfile, TString selection, TString varName, int nbins, double low, double up){
  TH1D* hist = new TH1D("hist","",nbins,low,up);
  TTree *tree = (TTree*)rootfile->Get("AnalysisTree");
  int Nev = tree->Project("hist",varName,selection);
  return hist;
}
