#include "../include/parameters.h"
#include "../include/useful_functions.h"

#include <TStyle.h>
#include <TROOT.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TMatrixD.h>
#include <TTree.h>
#include <TLegend.h>
#include <TMath.h>
#include <TMinuit.h>

#include <assert.h> 
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>


// --------------------------------- //
// any data to be used by chi2_linear for evaluating the function:
struct fit_data{
   // x, y values:
   std::vector<double> x_val, y_val;
   // variance-covariance matrix for y values:
   TMatrixD y_cov;
   // inverted cov matrix; calculated by chi2_linear "on demand".
   TMatrixD y_cov_inv;
   
   void reset(){
      x_val.clear();
      y_val.clear();
      y_cov.ResizeTo(0,0);
      y_cov_inv.ResizeTo(0,0);
   }

   void CheckPoints(){
      std::vector<int> RemovedPoints;
      TMatrixD y_cov_new;
      int j = 0;
    
      for(int i = 0; i < y_val.size(); i++) {
        std::cout << "i: " << i << "   j: " << j << std::endl;
         if( y_val.at(i) == 0) {
            x_val.erase(x_val.begin()+i);
            y_val.erase(y_val.begin()+i);
            RemovedPoints.push_back(j);
            i = i-1;
         }
         j++;
      }    
      for(int i = 0; i < x_val.size(); i++) {
          std::cout << "x: " << x_val.at(i) << std::endl;

      }

        std::cout<< "Removed Points: " << RemovedPoints.size() << std::endl;
        std::cout<< "Remaining Points: " << x_val.size() << std::endl;

      y_cov_new.ResizeTo(x_val.size(),x_val.size());
      for( int i=0; i < x_val.size(); i++) {
         for(int k= 0; k < x_val.size(); k++) {
            y_cov_new(i,k) = y_cov(i+RemovedPoints.size(),k+RemovedPoints.size());
         }
      }
      y_cov.ResizeTo(0,0);
      y_cov.ResizeTo(x_val.size(),x_val.size());
      y_cov = y_cov_new;  
   }
};

fit_data data;

// --------------------------------- //
// the chi^2 to minimize for fitting a linear function
//   y = p[0]*x + p[1]
// with fit parameters p[0], p[1] to data with known x and y and covariance
// matrix for y.
void chi2_linear(Int_t& npar, Double_t* grad, Double_t& fval, Double_t* p, Int_t status){
    if(data.y_cov_inv.GetNcols()==0){
        double dummy;
        int ncols = data.y_cov.GetNcols();
        data.y_cov_inv.ResizeTo(ncols, ncols);
        data.y_cov_inv = data.y_cov.Invert(&dummy);
    }
    const size_t ndata = data.x_val.size(); // number of data points in x,y graph to fit to
    std::vector<double> delta_y(ndata);
    for(size_t i=0; i<ndata; ++i){
        delta_y[i] = data.x_val[i]*p[0] + p[1] - data.y_val[i];
    }
    // now calculate the chi2, i.e.
    //  dy^T * C^{-1} * dy
    // where C is the variance--covariance matrix and dy = (y_data - y_pred)
    // This could probably be implemented in ROOT, but it's so simple, we just do it here:
    fval = 0.0;
    for(size_t i=0; i<ndata; ++i){
        for(size_t j=0; j<ndata; ++j){
            fval += delta_y[i] * delta_y[j] * data.y_cov_inv(i,j);
        }
    }
}

// --------------------------------- //
void make_lin_fit(double & slope, double & d_slope, double & offset, double & d_offset){
  TMinuit *min = new TMinuit();
  //min->SetPrintLevel(-1);
    min->SetPrintLevel(0);
    int err = min->DefineParameter(0, "slope", slope, d_slope, -1, 1);
    cout<<"Step size: "<<d_slope<<endl;
    assert(err==0);
    err = min->DefineParameter(1, "offset", offset, d_offset, 0.5,1.5);
    assert(err==0);
    min->SetFCN(chi2_linear);
    min->mnmigr();
    min->GetParameter(0, slope, d_slope);
    min->GetParameter(1, offset, d_offset);
}


// --------------------------------- //
TF1* kFSR_Fit(TGraphErrors* h_mean, TGraphErrors* h_width, int i, int j)
{

   std::vector<double> alpha;
   alpha.push_back(0.05); 
   alpha.push_back(0.1); 
   alpha.push_back(0.15); 
   alpha.push_back(0.2); 
   alpha.push_back(0.25); 
   alpha.push_back(0.3);
   alpha.push_back(0.35); 
   alpha.push_back(0.4);
   alpha.push_back(0.45);

   double *mean_old  = h_mean->GetY();
   double *err_old   = h_mean->GetEY();;
   double *width_old = h_width->GetEY();;
       
  int count=0;
  vector<double> mean, err, width;
  for(int i=0;i<h_mean->GetN();i++){
  
      count++;
      mean .push_back(mean_old[i]);       
      err  .push_back(err_old[i]);
      width.push_back(width_old[i]);       
  }

         // Covariance matrices needed for fitting 
         TMatrixD Mean_cov;

         Mean_cov.ResizeTo(n_alpha, n_alpha);
   

         // fill covariance matrix for data and mc
         for(int ialpha=0; ialpha < n_alpha; ++ialpha){
	   for (Int_t jalpha =0; jalpha < n_alpha; jalpha++){
               if( ialpha <= jalpha ) {
                      
		 Mean_cov(ialpha, jalpha) = pow(width[ialpha],2)*(pow(err[jalpha],2)/pow(width[jalpha],2));
    
               }
               else {
        	 Mean_cov(ialpha, jalpha) = pow(width[jalpha],2)*(pow(err[ialpha],2)/pow(width[ialpha],2));
               }
            }
         }        
  
 
         // fit linear extrapolation function
         TF1 *lin_extrapol = new TF1("lin_extrapol","[0]+[1]*x",0,alpha.back()+0.05); 
   
         //fit extrapolation function to the TGraphErrors
                   
         // fit
         data.reset();
         data.x_val = alpha;
         data.y_val = mean;
         data.y_cov.ResizeTo(n_alpha,n_alpha);
         data.y_cov = Mean_cov;
         data.CheckPoints();
         
         // choose start values for MC fit
         double slope = (mean[n_alpha-3] - mean[n_alpha-5])/(alpha_bins[n_alpha-3] - alpha_bins[n_alpha-5]);
         double d_slope = abs(slope*0.2);
         double offset = mean[n_alpha-3] - (slope*alpha_bins[n_alpha-5]);
         double d_offset = abs(offset*0.2);
             
         std::cout << "MC start values: " << "slope: " << slope <<" d_slope: "<<d_slope << "  offset: " << offset <<" d_offset: "<<d_offset<< std::endl; 
         make_lin_fit(slope, d_slope, offset, d_offset);
         std::cout << "MC fit values: " << "slope: " << slope << " offset: " << offset << std::endl; 

         lin_extrapol->SetParameter(0, offset);
         lin_extrapol->SetParError(0, d_offset);
         lin_extrapol->SetParameter(1, slope);
         lin_extrapol->SetParError(1, d_slope);

	 cout<<"Get function"<<endl;
	 
         h_mean->GetListOfFunctions()->Add(lin_extrapol);
         
         data.reset();

         // draw extrapolations data + mc
         TCanvas *c = new TCanvas("c","",600,600);
         c->DrawFrame(0,0.05,alpha.back()+0.05,1.47,(";Threshold #alpha_{max};#mu_{A}"));
	 h_mean->GetYaxis()->SetRangeUser(0.8, 1.2);
         h_mean->SetMarkerStyle(20);
         h_mean->SetMarkerColor(kRed+1);
         h_mean->SetLineColor(kRed+1);
	 cout<<"Draw Hist"<<endl;
         h_mean->Draw("P");
	 cout<<"After Hist"<<endl;
         TF1* Temp = new TF1();
         h_mean->GetFunction("lin_extrapol")->SetLineColor(kRed+1);
         h_mean->GetFunction("lin_extrapol")->SetLineStyle(2);
	 cout<<"Set Function Style"<<endl;
	   
	 Temp=(TF1*) h_mean->GetFunction("lin_extrapol")->Clone();
         Temp->SetRange(0.1,1);
	 Temp->SetLineStyle(1);
	 cout<<"Draw"<<endl;
         Temp->Draw("same");
	 cout<<"Print"<<endl;
  	 c->Print("/nfs/dust/cms/user/multh/JEC/2016ReReco/Residuals/Summer16_03Feb2017_V3/AK4CHS/MC_Reweighted_chsMET_NewTriggerSetup/RunBCDEFGH/plots/kFSR_Pt_eta_"+eta_range2[i]+"_"+eta_range2[i+1]+"_pT_"+pt_range[j]+"_"+pt_range[j+1]+".pdf");




	 return lin_extrapol;
}


