// This script stores ratio between MC and DATA responses as function of pT in different eta bins and for different alpha values

#include "header.h"


//Read histrogram from input file; method = mpf for MPF, = r_rel for pT-balance
TH1D* GetHistPtBin(TFile *rootfile, int alphabin, int etabin,int pTbin, TString method){
  TString histname = alpha_range_common[alphabin]+"/eta_"+eta_common_range[etabin]+"_"+eta_common_range[etabin+1]+"/pt_"+pt_range[pTbin]+"_"+pt_range[pTbin+1]+"/"+method;
  //  TString histname = "/eta_"+eta_common_range[etabin]+"_"+eta_common_range[etabin+1]+"/pt_"+pt_range[pTbin]+"_"+pt_range[pTbin+1]+"/"+method;//TEST
  TH1D* hist = (TH1D*)rootfile->Get(histname);
  return hist;
}

//Calculate ratio between MC and DATA responses
//return ratio and error
pair<double,double> Rmc_to_Rdata(TH1D* mc, TH1D* data){
  double Ratio = mc->GetMean() / data->GetMean();
  //double error = sqrt( pow(data->GetRMS(),2)/data->Integral() + pow(mc->GetRMS(),2)/mc->Integral());
  double error = sqrt( pow(data->GetMeanError(),2) + pow(mc->GetMeanError(),2)); //TEST
  // cout<<"MC_mean = "<<mc->GetMean()<<" ("<<mc->GetEntries()<<") DATA_mean = "<<data->GetMean()<<" ("<<data->GetEntries()<<") ratio = "<<Ratio<<endl;
  // cout<<"MC MeanError = "<<data->GetMeanError()<<" MC tricky RMS = "<<sqrt(pow(mc->GetRMS(),2)/mc->Integral())<<endl;
  // cout<<" "<<endl;
  pair<double,double> out;
  out.first = Ratio;
  out.second = error;
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
    if(Yval[i]!=0 && YvalError[i]!=0){
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

void InputForGlobalFit(TString path, TFile* datafile, TFile* MCfile){
  const int min_number_events = 100;//required number of events in response histogram
  gStyle->SetOptFit(0);
  //datafile->Print();
  //MCfile->Print();


  // get the histos and calculate ration in MC to DATA responses
  //TH1D* data,mc;
  double mc_al_mpf[n_alpha_common-1][n_eta_common-1][n_pt-1]; 
  double err_mc_al_mpf[n_alpha_common-1][n_eta_common-1][n_pt-1]; 
  double mc_al_pTbal[n_alpha_common-1][n_eta_common-1][n_pt-1]; 
  double err_mc_al_pTbal[n_alpha_common-1][n_eta_common-1][n_pt-1]; 
  double data_al_mpf[n_alpha_common-1][n_eta_common-1][n_pt-1]; 
  double err_data_al_mpf[n_alpha_common-1][n_eta_common-1][n_pt-1]; 
  double data_al_pTbal[n_alpha_common-1][n_eta_common-1][n_pt-1]; 
  double err_data_al_pTbal[n_alpha_common-1][n_eta_common-1][n_pt-1]; 

  double ratio_al_mpf[n_alpha_common-1][n_eta_common-1][n_pt-1]; //ratio at pt,eta,alpha bins
  double err_ratio_al_mpf[n_alpha_common-1][n_eta_common-1][n_pt-1]; //error of ratio at pt,eta,alpha bins
  double ratio_al_pTbal[n_alpha_common-1][n_eta_common-1][n_pt-1]; //ratio at pt,eta,alpha bins
  double err_ratio_al_pTbal[n_alpha_common-1][n_eta_common-1][n_pt-1]; //error of ratio at pt,eta,alpha bins


  //Initialize with 0 values
  for(int i=0; i<n_alpha_common-1; i++){
    for(int j=0; j<n_eta_common-1; j++){
      for(int k=0; k<n_pt-1; k++){
	ratio_al_mpf[i][j][k] = 0;
	err_ratio_al_mpf[i][j][k] = 0;
	ratio_al_pTbal[i][j][k] = 0;
	err_ratio_al_pTbal[i][j][k] = 0;
	mc_al_mpf[i][j][k] = 0;
	err_mc_al_mpf[i][j][k] = 0; 
	mc_al_pTbal[i][j][k] = 0;  
	err_mc_al_pTbal[i][j][k] = 0;
	data_al_mpf[i][j][k] = 0; 
	err_data_al_mpf[i][j][k] = 0; 
	data_al_pTbal[i][j][k] = 0; 
	err_data_al_pTbal[i][j][k] = 0; 
      }
    }
  }

  for(int i=0; i<n_alpha_common-1; i++){
    for(int j=0; j<n_eta_common-1; j++){
      for(int k=0; k<n_pt-1; k++){
	//	cout<<"For alpha = "<<alpha_range_common[i]<<" eta = "<<eta_common_range[j]<<" pT = "<<pt_range[k]<<endl;
	TH1D* mc_mpf  = GetHistPtBin(MCfile,i,j,k,"mpf");
	TH1D* data_mpf  = GetHistPtBin(datafile,i,j,k,"mpf");
	if(mc_mpf->GetEntries()<min_number_events || data_mpf->GetEntries()<min_number_events) continue;
	pair<double,double> ratio_pair_mpf = Rmc_to_Rdata(mc_mpf,data_mpf);
	ratio_al_mpf[i][j][k+1] =  ratio_pair_mpf.first; //TEST
	err_ratio_al_mpf[i][j][k+1] =  ratio_pair_mpf.second; //TEST
	TH1D* mc = GetHistPtBin(MCfile,i,j,k,"r_rel");
	TH1D* data = GetHistPtBin(datafile,i,j,k,"r_rel");
	if(mc->GetEntries()<min_number_events || data->GetEntries()<min_number_events) continue;
	pair<double,double> ratio_pair = Rmc_to_Rdata(mc,data);
	ratio_al_pTbal[i][j][k] =  ratio_pair.first;
	err_ratio_al_pTbal[i][j][k] =  ratio_pair.second;
	mc_al_mpf[i][j][k] = mc_mpf->GetMean(); 
	err_mc_al_mpf[i][j][k] = mc_mpf->GetMeanError(); 
	mc_al_pTbal[i][j][k] = mc->GetMean();  
	err_mc_al_pTbal[i][j][k] = mc->GetMeanError();
	data_al_mpf[i][j][k] = data_mpf->GetMean(); 
	err_data_al_mpf[i][j][k] = data_mpf->GetMeanError(); 
	data_al_pTbal[i][j][k] = data->GetMean(); 
	err_data_al_pTbal[i][j][k] = data->GetMeanError(); 
      }
    }
  }


  //Store results in TGraphErrors
  double res_xbin_tgraph[n_pt-1];// = {(pt_bins[0]+pt_bins[1])/2
  double res_zero[n_pt-1];
  for(int i=0;i<n_pt-1;i++){
    res_xbin_tgraph[i]=(pt_bins[i]+pt_bins[i+1])/2;
    res_zero[i]=(pt_bins[i+1]-pt_bins[i])/2 ;
  }

  TGraphErrors *mpf_data[n_alpha_common-1][n_eta_common-1];
  TGraphErrors *rrel_data[n_alpha_common-1][n_eta_common-1];
  TGraphErrors *mpf_mc[n_alpha_common-1][n_eta_common-1];
  TGraphErrors *rrel_mc[n_alpha_common-1][n_eta_common-1];
  TGraphErrors *mpf_ratio[n_alpha_common-1][n_eta_common-1];
  TGraphErrors *rrel_ratio[n_alpha_common-1][n_eta_common-1];

  for(int i=0; i<n_alpha_common-1; i++){
    for(int j=0; j<n_eta_common-1; j++){
      mpf_data[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph,data_al_mpf[i][j],res_zero,err_data_al_mpf[i][j]);
      mpf_mc[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph,mc_al_mpf[i][j],res_zero,err_mc_al_mpf[i][j]);
      rrel_data[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph,data_al_pTbal[i][j],res_zero,err_data_al_pTbal[i][j]);
      rrel_mc[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph,mc_al_pTbal[i][j],res_zero,err_mc_al_pTbal[i][j]);
      mpf_ratio[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph,ratio_al_mpf[i][j],res_zero,err_ratio_al_mpf[i][j]);
      rrel_ratio[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph,ratio_al_pTbal[i][j],res_zero,err_ratio_al_pTbal[i][j]);

    }
  }

  // Cleaning for empty points (with low statistic)
  for(int i=0; i<n_alpha_common-1; i++){
    for(int j=0; j<n_eta_common-1; j++){
      mpf_data[i][j] = CleanEmptyPoints(mpf_data[i][j]);
      mpf_mc[i][j] = CleanEmptyPoints(mpf_mc[i][j]);
      rrel_data[i][j] = CleanEmptyPoints(rrel_data[i][j]);
      rrel_mc[i][j] = CleanEmptyPoints(rrel_mc[i][j]);
      mpf_ratio[i][j] = CleanEmptyPoints(mpf_ratio[i][j]);
      rrel_ratio[i][j] = CleanEmptyPoints(rrel_ratio[i][j]);

      mpf_data[i][j]->SetName("mpfchs_dijet_"+alpha_range_common[i]);
      mpf_mc[i][j]->SetName("mpfchs_dijet_"+alpha_range_common[i]);
      rrel_data[i][j]->SetName("ptchs_dijet_"+alpha_range_common[i]);
      rrel_mc[i][j]->SetName("ptchs_dijet_"+alpha_range_common[i]);
      mpf_ratio[i][j]->SetName("mpfchs_dijet_"+alpha_range_common[i]);
      rrel_ratio[i][j]->SetName("ptchs_dijet_"+alpha_range_common[i]);
      }
    }
 


  //Save results in root file
  TFile* outputfile = new TFile(path+"JECcombifile_Dijet.root","RECREATE");
  outputfile->Print();
 // TString eta_output[n_eta_common-1] = {"eta0000-0261", "eta0216-0522","eta0522-0783","eta0783-1044","eta1044-1305","eta1305-1653","eta1653-1930","eta1930-2172","eta2172-2322","eta2322-2500","eta2500-2650","eta2650-2853","eta2853-2964","eta2964-3139",
 // 				"eta3139-3489","eta3489-5191"};//TMP
  for(int i=0; i<n_alpha_common-1; i++){
    for(int j=0; j<n_eta_common-1; j++){
      if(i==0){
	outputfile->mkdir("ratio/"+eta_output[j]);
	outputfile->mkdir("data/"+eta_output[j]);
	outputfile->mkdir("mc/"+eta_output[j]);
      }
      outputfile->cd("ratio/"+eta_output[j]);
      mpf_ratio[i][j]->Write();
      rrel_ratio[i][j]->Write();
    
      outputfile->cd("data/"+eta_output[j]);
      mpf_data[i][j]->Write();
      rrel_data[i][j]->Write();
      
      outputfile->cd("mc/"+eta_output[j]);
      mpf_mc[i][j]->Write();
      rrel_mc[i][j]->Write();
    }
 }

  cout<<"Draw result for alpha = "<<alpha_range_common[3]<<" eta = "<<eta_common_range[n_eta_common-2]<<" "<<eta_common_range[n_eta_common-1]<<endl;
  mpf_ratio[3][n_eta_common-2]->Draw();

  outputfile->Write();
  outputfile->Close();

}
