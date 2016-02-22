// with this script kFSR extrapolations over alpha in bins of eta will be calculated for different pT bins
// works for MPF and pt-balance methods separetly

#include "header.h"


// TColor ColorCode(int k){
//   TColor colorName = "kBlack";
//   switch(k)
//     {
//     case 0:
//      colorName = "kBlack";
//     case 1:
//      colorName = "kPink";
//     case 2:
//       colorName = "kOrange+7";
//     case 3:
//       colorName = "kBlue";
//     case 4:
//      colorName = "kAzure";
//     case 5:
//      colorName = "kOrange-1";
//     case 6:
//      colorName = "kBlue";
//     case 7:
//       colorName = "kGreen";
//     case 8:
//       colorName = "kMagenta";
//     case 9:
//       colorName = "kSpring";
//     case 10:
//       colorName = "kRed";
//     default:
//       cout<<"Error, bad input, quitting\n";
//       break;
//     }
//   return colorName;
// }

TString ToString(int num) {
  ostringstream start;
  start<<num;
  TString start1=start.str();
  return start1;
}

//Read histrogram from input file; method = mpf for MPF, = r_rel for pT-balance
TH1D* GetHistSumPt(TFile *rootfile, int alphabin, int etabin, TString method){
TH1D* hist = (TH1D*)rootfile->Get(alpha_range[alphabin]+"/eta_"+eta_range[etabin]+"_"+eta_range[etabin+1]+"/"+method);
 return hist;
}

//Read histrogram from input file; method = mpf for MPF, = r_rel for pT-balance
TH1D* GetHistPtBin(TFile *rootfile, int alphabin, int etabin,int pTbin, TString method){
  TString histname = alpha_range[alphabin]+"/eta_"+eta_range[etabin]+"_"+eta_range[etabin+1]+"/pt_"+pt_range[pTbin]+"_"+pt_range[pTbin+1]+"/"+method;
  TH1D* hist = (TH1D*)rootfile->Get(histname);
  return hist;
}

//Calculate ratio between MC and DATA responses
//return ratio and error
pair<double,double> Rmc_to_Rdata(TH1D* mc, TH1D* data){
  double Ratio = mc->GetMean() / data->GetMean();
  //double error = sqrt( pow(data->GetRMS(),2)/data->Integral() + pow(mc->GetRMS(),2)/mc->Integral());
  double error = sqrt( pow(data->GetMeanError(),2) + pow(mc->GetMeanError(),2)); //TEST
  cout<<"MC_mean = "<<mc->GetMean()<<" ("<<mc->GetEntries()<<") DATA_mean = "<<data->GetMean()<<" ("<<data->GetEntries()<<") ratio = "<<Ratio<<endl;
  cout<<"MC MeanError = "<<data->GetMeanError()<<" MC tricky RMS = "<<sqrt(pow(mc->GetRMS(),2)/mc->Integral())<<endl;
  cout<<" "<<endl;
  pair<double,double> out;
  out.first = Ratio;
  out.second = error;
  return out;
}


//void kFSR(bool mpfMethod(false), TString path, TFile* datafile, TFile* MCfile){
void kFSR_pT(bool mpfMethod, TString path, TFile* datafile, TFile* MCfile){

  gStyle->SetOptFit(0);
  //datafile->Print();
  //MCfile->Print();

  // get the histos and calculate ration in MC to DATA responses
  //TH1D* data,mc;
  double ratio_al[n_pt-1][n_eta-1][n_alpha-1]; //ratio at pt,eta,alpha bins
  double err_ratio_al[n_pt-1][n_eta-1][n_alpha-1]; //error of ratio at pt,eta,alpha bins

  for(int i=0; i<n_alpha-1; i++){
    for(int j=0; j<n_eta-1; j++){
      for(int k=0; k<n_pt-1; k++){
	cout<<"For alpha = "<<alpha_range[i]<<" eta = "<<eta_range[j]<<" pT = "<<pt_range[k]<<endl;
	if(mpfMethod){
	  TH1D* mc  = GetHistPtBin(MCfile,i,j,k,"mpf");
	  TH1D* data  = GetHistPtBin(datafile,i,j,k,"mpf");
	  pair<double,double> ratio_pair = Rmc_to_Rdata(mc,data);
	  ratio_al[k][j][i] =  ratio_pair.first;
	  err_ratio_al[k][j][i] =  ratio_pair.second;
	}
	else{
	  TH1D* mc = GetHistPtBin(MCfile,i,j,k,"r_rel");
	  TH1D* data = GetHistPtBin(datafile,i,j,k,"r_rel");
	  pair<double,double> ratio_pair = Rmc_to_Rdata(mc,data);
	  ratio_al[k][j][i] =  ratio_pair.first;
	  err_ratio_al[k][j][i] =  ratio_pair.second;
	}
      }
    }
  }

  //Divide reponses by value at alpha=0.2
  int alpha_02_bin = 3;//position of alpha=0.2 bin
  for(int j=0; j<n_eta-1; j++){
    for(int k=0; k<n_pt-1; k++){
      double norm_al02 = ratio_al[k][j][alpha_02_bin];
      double err_norm_al02 = err_ratio_al[k][j][alpha_02_bin];
      for(int i=0; i<n_alpha-1; i++){
  	ratio_al[k][j][i] =   ratio_al[k][j][i]/norm_al02;
  	err_ratio_al[k][j][i] = sqrt(abs(pow(err_ratio_al[k][j][i],2)-pow(err_norm_al02,2)));
      }
    }
  }
  
  TLegend *leg1;
  leg1 = new TLegend(0.53,0.54,0.65,0.88,"","brNDC");//x+0.1
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.035);
  leg1->SetFillColor(10);
  leg1->SetLineColor(1);
  leg1->SetTextFont(42);

  TGraph *graph[n_pt-1][n_eta-1];//set of points vs alpha
  TMultiGraph *pTgraph[n_eta-1];//set of different pT bins in on eta bin
  double xbin_tgraph[n_alpha-1],zero[n_alpha-1];
  for(int i=1;i<n_alpha;i++){
    xbin_tgraph[i-1] = alpha_bins[i];
    zero[i-1] = 0;
  }
  for(int j=0; j<n_eta-1; j++){
    pTgraph[j] = new TMultiGraph();
    for(int k=0; k<n_pt-1; k++){
      graph[k][j] = new TGraphErrors(n_alpha-1,xbin_tgraph,ratio_al[k][j],zero,err_ratio_al[k][j]);
      //      graph[k][j]->SetMarkerStyle(20);
      // graph[k][j]->SetMarkerColor(kBlue-3+k);
      // graph[k][j]->SetLineColor(kBlue-3+k);
      graph[k][j]->SetMarkerSize(1.3);
      graph[k][j]->SetMarkerStyle(20+k);
      if(k<4){ //skip yellow color
      graph[k][j]->SetMarkerColor(k+1);
      graph[k][j]->SetLineColor(k+1);
      }
      else{
	graph[k][j]->SetMarkerColor(k+2);
	graph[k][j]->SetLineColor(k+2);
      }
      pTgraph[j]->Add(graph[k][j]);
      TString pTbin_label = "";
      pTbin_label+=pt_bins[k];
      pTbin_label+=" < p_{T} < ";
      pTbin_label+=pt_bins[k+1];
      //      cout<<j<<" "<<k<<" "<<pTbin_label<<endl;
      if(j==0) leg1->AddEntry(graph[k][j],pTbin_label,"epl");
    }
  }

// //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // // Plots

  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.03,1,0.47,1);

 
  // create output .dat file, including the kFSR extrapolation (alpha->0)
  FILE *fp; TH1D* kFSR_MPF; TH1D* kFSR_DiJet;
  if(mpfMethod){
    fp = fopen("KFSR_MPF_extrapolation.dat","w");
    kFSR_MPF = new TH1D("kfsr_mpf","kfsr_mpf", n_eta-1,eta_bins);
  }
  else{
    fp = fopen("KFSR_DiJet_extrapolation.dat","w");
    kFSR_DiJet = new TH1D("kfsr_dijet","kfsr_dijet", n_eta-1,eta_bins);
  }

  TH1D* plotkfsr = new TH1D("kfsr","kfsr", n_eta-1,eta_bins);


  


  //create plots
  TCanvas* a[n_eta-1];
  TString plotname[n_eta-1];
  TF1 *pol1[n_eta-1];
  for (int j=0; j<n_eta-1; j++){
  //  for (int j=0; j<1; j++){//TEST
    if(mpfMethod){
      plotname[j]="mpf_kfsr_diffPt_eta_"+eta_range[j]+"_"+eta_range[j+1];
    }
    else{
      plotname[j]="dijet_kfsr_diffPt_eta_"+eta_range[j]+"_"+eta_range[j+1];
    }
    a[j] = new TCanvas(plotname[j], plotname[j], 800,700);
    gStyle->SetOptTitle(0);
    pTgraph[j]->Draw("AP");
    //    pTgraph[j]->GetYaxis()->SetRangeUser(0.8,1.4);
    //    pTgraph[j]->GetYaxis()->SetRangeUser(0.9,1.1);
    pTgraph[j]->GetYaxis()->SetRangeUser(0.92,1.08);
    pTgraph[j]->GetYaxis()->SetTitle("(R_{MC}/R_{DATA})/(R_{MC}/R_{DATA})_{#alpha<0.2}");
    pTgraph[j]->GetXaxis()->SetTitle("cut on #alpha");

    pol1[j] = new TF1("pol1","pol1",0.14,0.46);  //TEST
    pol1[j]->SetParameters(0,0);
    //    TF1 *pol1 = new TF1("pol1","pol1",0.19,0.46); 
    //TF1 *pol1 = new TF1("pol1","pol1",0.09,0.46); 
    pTgraph[j]->Fit(pol1[j],"R");
    // //TF1 *pol1 = new TF1("pol1","pol1",0.14,0.46);  //TEST alpha range in fit 
    // TF1 *pol1 = new TF1("pol1","pol1",0.14,0.31);  //TEST alpha range in fit 
    // graph1[j]->Fit(pol1,"R");
    // //    graph1[j]->Fit("pol1","R");
    // graph1[j]->GetXaxis()->SetTitle("#alpha");
    // graph1[j]->GetXaxis()->SetTitleSize(0.05);
    // //graph1[j]->GetYaxis()->SetTitle("kFSR");
    // graph1[j]->GetXaxis()->SetLimits(0.,0.5);
    // graph1[j]->GetYaxis()->SetRangeUser(0.92,1.08);
    line->SetLineStyle(2);
    line->Draw("SAME");
    // fill the output.dat file
    if (fp!=NULL) {
      Float_t value = pol1[j]->GetParameter(0);
      Float_t uncert = pol1[j]->GetParError(0);
      fprintf(fp, "%f %f\n",value,uncert);
    }
    plotkfsr->SetBinContent(j+1,pol1[j]->GetParameter(0));
    plotkfsr->SetBinError(j+1,pol1[j]->GetParError(0));
    if(mpfMethod){
      kFSR_MPF->SetBinContent(j+1,pol1[j]->GetParameter(0));
      kFSR_MPF->SetBinError(j+1,pol1[j]->GetParError(0));
    }
    else{
      kFSR_DiJet->SetBinContent(j+1,pol1[j]->GetParameter(0));
      kFSR_DiJet->SetBinError(j+1,pol1[j]->GetParError(0));
    }
   
    if(mpfMethod){
      leg1->SetHeader("MPF, "+eta_range3[j]+"#leq|#eta|<"+eta_range3[j+1]);
    }
    else{
      leg1->SetHeader("p_{T} balance, "+eta_range3[j]+"#leq|#eta|<"+eta_range3[j+1]);
    }
    // leg1->AddEntry(graph1[j], "R(MC)/R(DATA)","P");
    // leg1->AddEntry(pol1, "linear fit","L");
    leg1->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    tex->DrawLatex(0.64,0.91,"2.11fb^{-1} (13TeV)");

    //save the plots
    if(mpfMethod){
      a[j]->Print("plots/kFSR_MPF_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
    }
    else{
      a[j]->Print("plots/kFSR_Pt_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
    }

  }
  fclose(fp);


  // create output file including the kFSR plot
  if(mpfMethod){
    TFile* outputfile = new TFile(path+"Histo_KFSR_MPF_L1.root","RECREATE");
    kFSR_MPF->Write();
    outputfile->Write();
    outputfile->Close();
  }
  else{
    TFile* outputfile = new TFile(path+"Histo_KFSR_DiJet_L1.root","RECREATE");
    kFSR_DiJet->Write();
    outputfile->Write();
    outputfile->Close();
  }


}
