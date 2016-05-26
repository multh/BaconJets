// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// with this script pT extrapolations in bins of eta will be calculated
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "header.h"



void PTextrapolation_COMB_noKFSR_TTree(TString path, TFile* datafile, TFile* MCfile, TString txttag, TString jettag, TString variation, TString tag, double al_cut=0.2, int nResponseBins=100){

  gStyle->SetOptFit(000);


  // Tag for time dependence plots
  //TString tag = "";

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




  // systematic uncertainty: "down" means that 0.5*ptave will be used in the loglinear fit. "up" means that 2*ptave will be used in the loglinear fit
  int syst = 0;
  if(variation=="central") syst=0;
  if(variation=="down") syst=1;
  if(variation=="up") syst=2;
  if(variation=="doubleup") syst=3;
  if(variation=="nominal") syst=5;

  // get the histos for pt average
  TH1D* ptave_data[n_eta-1];
  for(int i=0; i<n_eta-1; i++){
    TString selection = "alpha<";
    selection+=al_cut;
    selection+=" && fabs(probejet_eta)<";
    selection+=eta_range[i+1];
    selection+=" && fabs(probejet_eta)>=";
    selection+=eta_range[i];
    TString var1 = "pt_ave";   
    ptave_data[i] = (TH1D*)GetHist(datafile, selection, var1, 60,0,600)->Clone();
  }

 
 // get ratio for MC to DATA responses
  double ratio_mpf[n_eta-1][n_pt-1]; //ratio at pt,eta bins for alpha = al_cut
  double err_ratio_mpf[n_eta-1][n_pt-1]; //error of ratio at pt,eta bins for alpha = al_cut
  double ratio_dijet[n_eta-1][n_pt-1]; //ratio at pt,eta bins for alpha = al_cut
  double err_ratio_dijet[n_eta-1][n_pt-1]; //error of ratio at pt,eta bins for alpha = al_cut

  TH1D *hdata_rel_r[n_pt-1][n_eta-1];// pT-balance responce for data
  TH1D *hdata_mpf_r[n_pt-1][n_eta-1];//MPF responce for data
  TH1D *hmc_rel_r[n_pt-1][n_eta-1];// pT-balance responce for MC
  TH1D *hmc_mpf_r[n_pt-1][n_eta-1];//MPF responce for MC
  int count = 0;
  TString name1 = "hist_data_rel_r_";
  TString name2 = "hist_data_mpf_r_";
  TString name3 = "hist_mc_rel_r_";
  TString name4 = "hist_mc_mpf_r_";
    for(int j=0; j<n_eta-1; j++){
      for(int k=0; k<n_pt-1; k++){
	ratio_dijet[j][k] = 0;
	err_ratio_dijet[j][k] = 0;
	ratio_mpf[j][k] = 0;
	err_ratio_mpf[j][k] = 0;
	TString name = name1; name+=count;
	hdata_rel_r[k][j] = new TH1D(name,"",nResponseBins, 0, 2.5);
	name = name2;name+=count;
	hdata_mpf_r[k][j] = new TH1D(name,"",nResponseBins, 0, 2.5);
	name = name3; name+=count;
	hmc_rel_r[k][j] = new TH1D(name,"",nResponseBins, 0, 2.5);
	name = name4; name+=count;
	hmc_mpf_r[k][j] = new TH1D(name,"",nResponseBins, 0, 2.5);
	count++;
      }
    }



// Create the tree reader and its data containers
   TTreeReader myReader_DATA("Events", datafile);
   TTreeReaderValue<Float_t> pt_ave_data(myReader_DATA, "pt_ave");
   TTreeReaderValue<Float_t> probejet_eta_data(myReader_DATA, "probejet_eta");
   TTreeReaderValue<Float_t> alpha_data(myReader_DATA, "alpha");
   TTreeReaderValue<Float_t> rel_r_data(myReader_DATA, "rel_r");
   TTreeReaderValue<Float_t> mpf_r_data(myReader_DATA, "mpf_r");
   TTreeReaderValue<Float_t> weight_data(myReader_DATA, "weight");
   
   while (myReader_DATA.Next()) {
     //  cout<<"DATA point: alpha = "<<*alpha_data<<" eta_probe = "<<*probejet_eta_data<<" pT_ave = "<<*pt_ave_data<<endl;
     if(*alpha_data>al_cut) continue;
     for(int k=0; k<n_pt-1; k++){
   	   if(*pt_ave_data<pt_bins[k] || *pt_ave_data>pt_bins[k+1]) continue;
	   for(int j=0; j<n_eta-1; j++){
	     if(fabs(*probejet_eta_data)>eta_bins[j+1] || fabs(*probejet_eta_data)<eta_bins[j]) continue;
	     else{
	       // std::cout<<"DATA: in eta range = "<<eta_bins[j]<<", "<<eta_bins[j+1]
	       // 		<<" and pT range = "<<pt_bins[k]<<", "<<pt_bins[k+1]<<std::endl;
	       // std::cout<<"probe_eta = "<<*probejet_eta_data<<" pT_ave = "<<*pt_ave_data<<std::endl;

	       hdata_rel_r[k][j]->Fill(*rel_r_data,*weight_data);
	       hdata_mpf_r[k][j]->Fill(*mpf_r_data,*weight_data);
	     }
	   }
     }
   }

   TTreeReader myReader_MC("Events", MCfile);
   TTreeReaderValue<Float_t> pt_ave_mc(myReader_MC, "pt_ave");
   TTreeReaderValue<Float_t> probejet_eta_mc(myReader_MC, "probejet_eta");
   TTreeReaderValue<Float_t> alpha_mc(myReader_MC, "alpha");
   TTreeReaderValue<Float_t> rel_r_mc(myReader_MC, "rel_r");
   TTreeReaderValue<Float_t> mpf_r_mc(myReader_MC, "mpf_r");
   TTreeReaderValue<Float_t> weight_mc(myReader_MC, "weight");
   while (myReader_MC.Next()) {
     if(*alpha_mc>al_cut) continue;
     for(int k=0; k<n_pt-1; k++){
       if(*pt_ave_mc<pt_bins[k] || *pt_ave_mc>pt_bins[k+1]) continue;
       for(int j=0; j<n_eta-1; j++){
   	 if(fabs(*probejet_eta_mc)>eta_bins[j+1] || fabs(*probejet_eta_mc)<eta_bins[j]) continue;
   	   else{
   	     hmc_rel_r[k][j]->Fill(*rel_r_mc,*weight_mc);
   	     hmc_mpf_r[k][j]->Fill(*mpf_r_mc,*weight_mc);
   	   }
       }
     }
   }


     for(int j=0; j<n_eta-1; j++){
      for(int k=0; k<n_pt-1; k++){
	cout<<" "<<endl;
	std::cout<<"Ratio: For eta range = "<<eta_bins[j]<<", "<<eta_bins[j+1]<<" and pT range = "<<pt_bins[k]<<", "<<pt_bins[k+1]<<std::endl;
	pair<double,double> res_mc_rel_r,res_data_rel_r;
	pair<double,double> res_mc_mpf_r,res_data_mpf_r;
	res_mc_rel_r = GetValueAndError(hmc_rel_r[k][j]);
	res_data_rel_r = GetValueAndError(hdata_rel_r[k][j]);
	res_mc_mpf_r = GetValueAndError(hmc_mpf_r[k][j]);
	res_data_mpf_r = GetValueAndError(hdata_mpf_r[k][j]);

	
	pair<double,double> ratio_res_rel_r;
	if(res_mc_rel_r.first>0 && res_data_rel_r.first>0)
	  ratio_res_rel_r = Rmc_to_Rdata(res_mc_rel_r,res_data_rel_r);
	else 
	  ratio_res_rel_r.first = 0;

	pair<double,double> ratio_res_mpf_r;
	if(res_mc_mpf_r.first>0 && res_data_mpf_r.first>0)
	  ratio_res_mpf_r = Rmc_to_Rdata(res_mc_mpf_r,res_data_mpf_r);
	else 
	  ratio_res_mpf_r.first = 0;

	ratio_dijet[j][k] = ratio_res_rel_r.first;
	err_ratio_dijet[j][k] = ratio_res_rel_r.second;
	ratio_mpf[j][k] = ratio_res_mpf_r.first;
	err_ratio_mpf[j][k] = ratio_res_mpf_r.second;
	//	cout<<"ratio_dijet[k][j] = "<<ratio_dijet[k][j]<<" ratio_mpf[k][j] = "<<ratio_mpf[k][j]<<endl;
      }
    }
 
  // create and fill tgrapherrors
  double xbin_tgraph[n_pt-1];
  double zero[n_pt-1];
  for(int i=0;i<n_pt-1;i++){
    xbin_tgraph[i]=(pt_bins[i]+pt_bins[i+1])/2;
    zero[i]=(pt_bins[i+1]-pt_bins[i])/2 ;
  }
  TGraphErrors *graph1_mpf[n_eta-1];
  TGraphErrors *graph1_dijet[n_eta-1];
  TMultiGraph *mg1[n_eta-1]; 
  // double Vcov[3][n_eta-1];//covarance matrix for log lin fit result
  for(int j=0; j<n_eta-1; j++){
    std::cout<<"Cleaning For eta range = "<<eta_bins[j]<<", "<<eta_bins[j+1]<<std::endl;
    graph1_mpf[j] = new TGraphErrors(n_pt-1, xbin_tgraph, ratio_mpf[j], zero, err_ratio_mpf[j]);
    graph1_mpf[j] = (TGraphErrors*)CleanEmptyPoints(graph1_mpf[j]);
    graph1_dijet[j] = new TGraphErrors(n_pt-1, xbin_tgraph, ratio_dijet[j], zero, err_ratio_dijet[j]);
    graph1_dijet[j] = (TGraphErrors*)CleanEmptyPoints(graph1_dijet[j]);
    mg1[j] = new TMultiGraph();
  }

  
  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,pt_bins[n_pt-1]+10,1);

  // create a function for the loglinear fit
  TF1 * f1[n_eta-1];
  // create a function for the constant fit
  TF1 * f2[n_eta-1];

  // create output .dat file, including the kFSR extrapolation (alpha->0)
  FILE *fp; FILE *fp2; FILE *l2resfile;
  fp = fopen(path+"output/pT_COMB_extrapolations.dat","w");
  fp2 = fopen(path+"output/pT_COMB_constantExtrapolation.dat","w");
  l2resfile = fopen(path+"output/L2Res_COMB.dat","w");

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // Plots

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  TCanvas* asd[n_eta-1];
  TString plotname[n_eta-1];
  double Vcov[3][n_eta-1];//covarance matrix for log lin fit results
  for (int j=0; j<n_eta-1; j++){


  //for (int j=n_eta-2; j<n_eta-1; j++){//TEST
    plotname[j]="comb_ptextra_eta_"+eta_range[j]+"_"+eta_range[j+1];
    asd[j] = new TCanvas(plotname[j], plotname[j], 800,700);
    gStyle->SetOptTitle(0);
    gPad->SetLogx();
   
    //    f1[j] = new TF1(plotname[j]+"f1","[0]+[1]*TMath::Log(x)", 76 , 558);
    f1[j] = new TF1(plotname[j]+"f1","[0]+[1]*TMath::Log(x)", 50 , 570);
    //    f1[j] = new TF1(plotname[j]+"f1","[0]+[1]*TMath::Log(x)", 0 , 600);
    f1[j]->SetParameters(1,0);
    f2[j] = new TF1(plotname[j]+"f2","pol0",0,600);
    f1[j]->SetParameter(0,1);
    f2[j]->SetLineColor(kGreen+1);
    f1[j]->SetLineColor(kGreen+1);
    f2[j]->SetLineStyle(3);

    graph1_mpf[j]->SetMarkerColor(kRed+1);
    graph1_mpf[j]->SetLineColor(kRed+1);

    graph1_dijet[j]->SetMarkerColor(kBlue+1);
    graph1_dijet[j]->SetLineColor(kBlue+1);

    graph1_mpf[j]->SetMarkerStyle(20);
    graph1_dijet[j]->SetMarkerStyle(20);
    bool is_graph_points = false;
    if(graph1_mpf[j]->GetN()>0){ 
      mg1[j]->Add(graph1_mpf[j]);
      is_graph_points = true;
    }
    if(graph1_dijet[j]->GetN()>0){
      mg1[j]->Add(graph1_dijet[j]);
      is_graph_points = true;
    }
    if(is_graph_points){
    TFitResultPtr fitloglin = mg1[j]->Fit(plotname[j]+"f1","SR");
    TMatrixDSym cov = fitloglin->GetCovarianceMatrix();
    //    cov.Print();
    Vcov[0][j] = cov(0,0);     Vcov[1][j] = cov(1,1);     Vcov[2][j] = cov(0,1);
    //    cout<<"00 = "<<Vcov[0]<<" 11 = "<<Vcov[1]<<" 01 = "<<Vcov[2]<<endl;
    TFitResultPtr fitconst = mg1[j]->Fit(plotname[j]+"f2","R");
    mg1[j]->Draw("AP");
    f1[j]->Draw("same");
    f2[j]->Draw("same");
    //  mg2[j]->Draw("P SAME");
    mg1[j]->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    mg1[j]->GetXaxis()->SetTitleSize(0.05);
    mg1[j]->GetXaxis()->SetTitleOffset(0.80);
    mg1[j]->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+10);
    line->SetLineStyle(2);
    line->Draw("SAME");
 

    if (fp2!=NULL) {
      // getting the p0 parameter from the constant fit
      Float_t value = f2[j]->GetParameter(0);
      fprintf(fp2, "%f\n",value);
    }
    if (l2resfile!=NULL) {
      fprintf(l2resfile, "%f\n", eta_bins[j]);
    }
    asd[j]->Modified(); asd[j]->Update();
    TLegend *leg1;
    leg1 = new TLegend(0.12,0.68,0.35,0.88,"","brNDC");//x+0.1
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.038);
    leg1->SetFillColor(10);
    leg1->SetLineColor(1);
    leg1->SetTextFont(42);
    leg1->SetHeader("COMB #bar{p}_{T} extrapolation, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
   
    leg1->AddEntry(graph1_mpf[j], "R(MC)/R(DATA) MPF","P");
    leg1->AddEntry(graph1_dijet[j], "R(MC)/R(DATA) p_{T} balance","P");
    leg1->AddEntry(f1[j], "loglinear fit","L");
    leg1->AddEntry(f2[j], "constant fit","L");
    leg1->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    tex->DrawLatex(0.64,0.91,"2.11fb^{-1} (13TeV)");
    TString chi2_loglin = "loglinear fit #chi^{2}/n.d.f = ";
    chi2_loglin += trunc(f1[j]->GetChisquare());
    chi2_loglin +="/";
    chi2_loglin +=trunc(f1[j]->GetNDF());
    TString chi2_const = "constant fit #chi^{2}/n.d.f = ";
    chi2_const+=trunc(f2[j]->GetChisquare());
    chi2_const+="/";
    chi2_const+=trunc(f2[j]->GetNDF());

    TLatex *tex2 = new TLatex();
    tex2->SetNDC();
    tex2->SetTextSize(0.035); 
    tex2->DrawLatex(0.51,0.73,chi2_loglin);
    tex2->DrawLatex(0.51,0.69,chi2_const);


    mg1[j]->GetYaxis()->SetRangeUser(0.9,1.15);
    // if(j<11){
    //   mg1[j]->GetYaxis()->SetRangeUser(0.95,1.10);
    // }
    // if(j>=11){
    //   mg1[j]->GetYaxis()->SetRangeUser(0.95,1.30);
    // }
    //save plots
    asd[j]->Print(path+"plots/pTextrapolation_COMB_pT_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
    }
  }
  fclose(fp);
  fclose(fp2);


  
  //==============================================================================================================================
  double flat[n_eta-1];
  double loglin[n_eta-1];
  for (int j=0; j<n_eta-1; j++){
    flat[j] = f2[j]->GetParameter(0);
    loglin[j] = f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean());
    cout<<"ptave_data[j]->GetMean() = "<<ptave_data[j]->GetMean()<<" flat = "<<flat[j]<<" loglin = "<<loglin[j]<<endl;
  }
  double flat_norm = 0; double loglin_norm = 0;
  for (int j=0; j<n_etabarr; j++){
    flat_norm += flat[j];
    loglin_norm += loglin[j];
  }
  flat_norm = flat_norm/n_etabarr;
  loglin_norm = loglin_norm/n_etabarr;
  for (int j=0; j<n_eta-1; j++){
    flat[j] = flat[j]/flat_norm;
    loglin[j] = loglin[j]/loglin_norm;
  }

  cout << "normalization factor (flat): " << flat_norm << endl;
  cout << "normalization factor (logl): " << loglin_norm << endl;
  
  TH1D* Residual_logpt_COMB = new TH1D("res_logpt_comb","res_logpt_comb", n_eta-1,eta_bins);
  TH1D* Residual_const_COMB = new TH1D("res_const_comb","res_const_comb", n_eta-1,eta_bins);
  TH1D* ptave_const_COMB = new TH1D("ptave_const_comb","ptave_const_comb", n_eta-1,eta_bins);
  TH1D* ptave_logpt_COMB = new TH1D("ptave_logpt_comb","ptave_logpt_comb", n_eta-1,eta_bins);
  ofstream output, output_loglin, uncerts, uncerts_loglin;
  output.open(path+"output/Fall15_25ns_COMB_FLAT_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt");
  output_loglin.open(path+"output/Fall15_25ns_COMB_LOGLIN_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt");
  uncerts.open(path+"output/Fall15_25ns_COMB_FLAT_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt.STAT");
  uncerts_loglin.open(path+"output/Fall15_25ns_COMB_LOGLIN_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt.STAT");
  //output  << "JetEta (abs)      " << "p[0]" << "       " <<  "kFSR * p[0] " << "  statistical unc" << endl;
  //output_loglin  << "JetEta (abs)      "  << "p[0]        "  << "p[1]"  << "       " <<  "kFSR*[p[0]+p[1]*Log(max(ptmin,min(ptmax, x)) )] " << "  statistical unc" << endl;
  output  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(208))/208))) Correction L2Relative}" << endl;
  output_loglin  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(208))/208))) Correction L2Relative}" << endl;
  uncerts << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}" << endl;
  uncerts_loglin << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}" << endl;
  
  for (int j=n_eta-1; j>0; --j){
    output << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   55   " << ptave_data[j-1]->FindLastBinAbove(10.)*10 << "   " << 1/flat_norm << " " << 1.0 << "   " << f2[j-1]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;
    output_loglin << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   55   " << ptave_data[j-1]->FindLastBinAbove(10.)*10 << "   " << 1/loglin_norm << " " << 1.0 << "   " << f1[j-1]->GetParameter(0) << " " << f1[j-1]->GetParameter(1) << "   1 0.0000 0.0" << endl;
    uncerts << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " << ptave_data[j-1]->FindLastBinAbove(10.)*10 << "       " << sqrt( f2[j-1]->GetParError(0) ) / flat_norm  << endl;
    uncerts_loglin << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " << ptave_data[j-1]->FindLastBinAbove(10.)*10 << "       " << sqrt(Vcov[0][j-1]+Vcov[1][j-1]*pow(TMath::Log(ptave_data[j-1]->GetMean()),2)+2*Vcov[2][j-1]*TMath::Log(ptave_data[j-1]->GetMean())) / loglin_norm << endl;
  }
  
  for (int j=0; j<n_eta-1; j++){
    output << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   55   " << ptave_data[j]->FindLastBinAbove(10.)*10 << "   " << 1/flat_norm << " " << 1.0 << "   " << f2[j]->GetParameter(0) << " 0   1 0.0000 0.0"  << endl;
      output_loglin << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   55   " << ptave_data[j]->FindLastBinAbove(10.)*10 << "   " << 1/loglin_norm << " " << 1.0 << "   " << f1[j]->GetParameter(0) << " " << f1[j]->GetParameter(1) << "   1 0.0000 0.0" << endl;
      uncerts << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " << ptave_data[j]->FindLastBinAbove(10.)*10 << "       " << sqrt( f2[j]->GetParError(0)) / flat_norm  << endl;
      uncerts_loglin << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " << ptave_data[j]->FindLastBinAbove(10.)*10 << "       " << sqrt(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(ptave_data[j]->GetMean()),2)+2*Vcov[2][j]*TMath::Log(ptave_data[j]->GetMean())) / loglin_norm << endl;

      cout<<"uncert for log lin["<<j<<"] "<<f1[j]->GetParError(0)<<" "<<f1[j]->GetParError(1)<<" together = "<<sqrt(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(ptave_data[j]->GetMean()),2)+2*Vcov[2][j]*TMath::Log(ptave_data[j]->GetMean()))<<endl;//TEST

      if(syst==5){
      	Residual_logpt_COMB->SetBinContent(j+1, (f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean())));
      	Residual_logpt_COMB->SetBinError(j+1, sqrt(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(ptave_data[j]->GetMean()),2)+2*Vcov[2][j]*TMath::Log(ptave_data[j]->GetMean())));
      }
      if(syst==0){
  	Residual_logpt_COMB->SetBinContent(j+1,(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(120)));
  	Residual_logpt_COMB->SetBinError(j+1, sqrt(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(120),2)+2*Vcov[2][j]*TMath::Log(120)));
      }
      if(syst==1){
      	Residual_logpt_COMB->SetBinContent(j+1,1.0*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(60) ) );
	Residual_logpt_COMB->SetBinError(j+1, sqrt(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(60),2)+2*Vcov[2][j]*TMath::Log(60)));
      }
      if(syst==2){
      	Residual_logpt_COMB->SetBinContent(j+1,1.0*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(240) ) );
	Residual_logpt_COMB->SetBinError(j+1, sqrt(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(240),2)+2*Vcov[2][j]*TMath::Log(240)));
      }
      if(syst==3){
      	Residual_logpt_COMB->SetBinContent(j+1,1.0*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(480) ) );
	Residual_logpt_COMB->SetBinError(j+1, sqrt(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(480),2)+2*Vcov[2][j]*TMath::Log(480)));
      }

      Residual_const_COMB->SetBinContent(j+1,f2[j]->GetParameter(0));
      Residual_const_COMB->SetBinError(j+1,f2[j]->GetParError(0));
      ptave_const_COMB->SetBinContent(j+1,f2[j]->GetParameter(0));
      ptave_const_COMB->SetBinError(j+1,f2[j]->GetParError(0));
      ptave_logpt_COMB->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()));
      ptave_logpt_COMB->SetBinError(j+1,sqrt(Vcov[0][j]+Vcov[1][j]*pow(TMath::Log(ptave_data[j]->GetMean()),2)+2*Vcov[2][j]*TMath::Log(ptave_data[j]->GetMean())));
    }




    TFile* outputfile;
    if(syst==5){
      outputfile = new TFile(path+"Histo_Res_COMB_L1"+tag+".root","RECREATE");
    }
    if(syst==0){
     outputfile = new TFile(path+"Histo_Res_COMB_L1_"+variation+".root","RECREATE");
    }
    if(syst==1){
     outputfile = new TFile(path+"Histo_Res_COMB_L1_"+variation+".root","RECREATE");
    }
    if(syst==2){
      outputfile = new TFile(path+"Histo_Res_COMB_L1_"+variation+".root","RECREATE");
    }
    if(syst==3){
      outputfile = new TFile(path+"Histo_Res_COMB_L1_"+variation+".root","RECREATE");
    }

    Residual_logpt_COMB->Write();
    Residual_const_COMB->Write();
    outputfile->Write();
    outputfile->Close();
    TFile* outputfile2 = new TFile(path+"Histo_ptave_COMB_L1.root","RECREATE");
    ptave_const_COMB->Write();
    ptave_logpt_COMB->Write();
    outputfile2->Write();
    outputfile2->Close();
  




}
