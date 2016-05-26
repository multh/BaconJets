// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// with this script pT extrapolations in bins of eta will be calculated
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "header.h"


void PTextrapolation(bool mpfMethod, TString path, TFile* datafile, TFile* MCfile, TString txttag, TString jettag, TString variation, TString tag){
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
    TString selection = "alpha<0.3 && probejet_eta<";
    selection+=eta_range[i+1];
    selection+=" && probejet_eta>=";
    selection+=eta_range[i];
    TString var1 = "pt_ave";   
    ptave_data[i] = (TH1D*)GetHist(datafile, selection, var1, 300,0,3000)->Clone();
  }



  // get ratio for MC to DATA responses
  double ratio_mpf[n_eta-1][n_pt-1]; //ratio at pt,eta bins for alpha = 0.3
  double err_ratio_mpf[n_eta-1][n_pt-1]; //error of ratio at pt,eta bins for alpha = 0.3
  for(int j=0; j<n_eta-1; j++){
    for(int k=0; k<n_pt-1; k++){
      std::cout<<"For eta range = "<<eta_bins[j]<<", "<<eta_bins[j+1]<<" and pT range = "<<pt_bins[k]<<", "<<pt_bins[k+1]<<std::endl;
      pair<double,double> res_mc_mpf =  Response(MCfile,0.3,eta_bins[j],eta_bins[j+1],pt_bins[k],pt_bins[k+1],mpfMethod);
      pair<double,double> res_data_mpf =  Response(datafile,0.3,eta_bins[j],eta_bins[j+1],pt_bins[k],pt_bins[k+1],mpfMethod);
      pair<double,double> ratio_res_mpf = Rmc_to_Rdata(res_mc_mpf,res_data_mpf);
      ratio_mpf[j][k] = ratio_res_mpf.first;
      err_ratio_mpf[j][k] = ratio_res_mpf.second;
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
 
 
  for(int j=0; j<n_eta-1; j++){
    graph1_mpf[j] = new TGraphErrors(n_pt-1, xbin_tgraph, ratio_mpf[j], zero, err_ratio_mpf[j]);
    graph1_mpf[j] = (TGraphErrors*)CleanEmptyPoints(graph1_mpf[j]);
  }

  

  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,pt_bins[n_pt-1]+10,1);

  // create a function for the loglinear fit
  TF1 * f1[n_eta-1];
  // create a function for the constant fit
  TF1 * f2[n_eta-1];

  // create output .dat file, including the kFSR extrapolation (alpha->0)
  FILE *fp; FILE *fp2; FILE *l2resfile;
  if(mpfMethod){
    fp = fopen(path+"output/pT_MPF_extrapolations.dat","w");
    fp2 = fopen(path+"output/pT_MPF_constantExtrapolation.dat","w");
    l2resfile = fopen(path+"output/L2Res_MPF.dat","w");
  }
  else{
    fp = fopen(path+"output/pT_DiJet_extrapolations.dat","w");
    fp2 = fopen(path+"output/pT_DiJet_constantExtrapolation.dat","w");
    l2resfile = fopen(path+"output/L2Res_DiJet.dat","w");
  }


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // Plots

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  TCanvas* asd[n_eta-1];
  TString plotname[n_eta-1];

  for (int j=0; j<n_eta-1; j++){
    if(mpfMethod){
      plotname[j]="mpf_ptextra_eta_"+eta_range[j]+"_"+eta_range[j+1];
    }
    else{
      plotname[j]="dijet_ptextra_eta_"+eta_range[j]+"_"+eta_range[j+1];
    }
    asd[j] = new TCanvas(plotname[j], plotname[j], 800,700);
    gStyle->SetOptTitle(0);
    gPad->SetLogx();
    graph1_mpf[j]->SetMarkerColor(kBlue);
    graph1_mpf[j]->SetMarkerStyle(20);
    graph1_mpf[j]->SetLineColor(kBlue);
   
 
    f1[j] = new TF1(plotname[j]+"f1","[0]+[1]*TMath::Log(x)", 76 , 558);
    f1[j]->SetParameters(1,0);
    f2[j] = new TF1(plotname[j]+"f2","pol0");
    f2[j]->SetLineColor(kBlue);
    f2[j]->SetLineStyle(3);

    graph1_mpf[j]->Fit(plotname[j]+"f1","");
    graph1_mpf[j]->Fit(plotname[j]+"f2","+ SAME");
    graph1_mpf[j]->Draw("AP");
    graph1_mpf[j]->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph1_mpf[j]->GetXaxis()->SetTitleSize(0.05);
    graph1_mpf[j]->GetXaxis()->SetTitleOffset(0.80);
    graph1_mpf[j]->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+10);
    graph1_mpf[j]->GetYaxis()->SetRangeUser(0.95,1.10);
 
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
    TPaveStats *st = ((TPaveStats*)(graph1_mpf[j]->GetListOfFunctions()->FindObject("stats")));
    if (st) {
      st->SetTextColor(kRed);
      st->SetX1NDC(0.69); st->SetX2NDC(0.99);
      st->SetY1NDC(0.65); st->SetY2NDC(0.8);
    }
    st = ((TPaveStats*)(graph1_mpf[j]->GetListOfFunctions()->FindObject("stats")));
    if (st) {
      st->SetTextColor(graph1_mpf[j]->GetLineColor());
      st->SetX1NDC(0.69); st->SetX2NDC(0.99);
      st->SetY1NDC(0.85); st->SetY2NDC(0.95);
    }
    asd[j]->Modified(); asd[j]->Update();

    TLegend *leg1;
    leg1 = new TLegend(0.12,0.68,0.35,0.88,"","brNDC");//x+0.1
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.038);
    leg1->SetFillColor(10);
    leg1->SetLineColor(1);
    leg1->SetTextFont(42);
    if(mpfMethod){
      leg1->SetHeader("MPF #bar{p}_{T} extrapolation, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
    }
    else{
      leg1->SetHeader("p_{T} balance #bar{p}_{T} extrapolation, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
    }
    leg1->AddEntry(graph1_mpf[j], "R(MC)/R(DATA)","P");
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

    //save plots
    if(mpfMethod){
      asd[j]->Print(path+"plots/pTextrapolation_MPF_pT_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
    }
    else{
      asd[j]->Print(path+"plots/pTextrapolation_Pt_pT_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
    }
  }
  fclose(fp);
  fclose(fp2);


  //==============================================================================================================================
  // some couts for cross check
  for(int i=0; i<n_eta-1; i++){
    cout << "mean value, eta " << eta_range[i] << " to " << eta_range[i+1] << ": " << ptave_data[i]->GetMean() << endl;
    //cout << "max value, eta " << eta_range[i] << " to " << eta_range[i+1] << ": " << ptave_data[i]->FindLastBinAbove() << endl;
  }
  for(int i=0; i<n_eta-1; i++){
    cout << "loglin fit value " << eta_range[i] << " to " << eta_range[i+1] << ": " << f1[i]->GetParameter(1) << endl;
  }
  for(int i=0; i<n_eta-1; i++){
    cout << "max value of pt " << eta_range[i] << " to " << eta_range[i+1] << ": " << ptave_data[i]->FindLastBinAbove(0.) << endl;
  }

  for(int i=0; i<n_eta-1; i++){
    for(int j=0; j<graph1_mpf[i]->GetN(); j++) {
      cout << "uncert, eta " << eta_range[i] << " to " << eta_range[i+1] << ", pt range " << pt_range[j] << " to " << pt_range[j+1] << ": "  << graph1_mpf[i]->GetEY()[j] << endl;
    }
  }


  //==============================================================================================================================
  // get the kFSR plots and calculate residuals
  if(mpfMethod){
    TFile* kfsr_mpf = new TFile(path+"Histo_KFSR_MPF_L1.root","READ");
    TH1D* hist_kfsr_mpf = (TH1D*)kfsr_mpf->Get("kfsr_mpf");
    /*
    for(int i=0; i<50; i++){
      hist_kfsr_mpf->SetBinContent(i,1);
    }
    */
    double flat[n_eta-1];
    double loglin[n_eta-1];
    for (int j=0; j<n_eta-1; j++){
      flat[j] = hist_kfsr_mpf->GetBinContent(j+1)*f2[j]->GetParameter(0);
      loglin[j] = hist_kfsr_mpf->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()) );
    }
    double flat_norm, loglin_norm;
    for (int j=0; j<n_etabarr; j++){
      flat_norm += flat[j];
      loglin_norm += loglin[j];
    }
    flat_norm = flat_norm/n_etabarr;
    loglin_norm = loglin_norm/n_etabarr;
    for (int j=0; j<n_eta-1; j++){
      flat[j] = flat[j]/flat_norm;
      loglin[j] = loglin[j]/loglin_norm;
      //cout << flat[j] << endl;
    }

    cout << "normalization factor (flat): " << flat_norm << endl;
    cout << "normalization factor (logl): " << loglin_norm << endl;

    TH1D* Residual_logpt_MPF = new TH1D("res_logpt_mpf","res_logpt_mpf", n_eta-1,eta_bins);
    TH1D* Residual_const_MPF = new TH1D("res_const_mpf","res_const_mpf", n_eta-1,eta_bins);
    TH1D* ptave_const_MPF = new TH1D("ptave_const_mpf","ptave_const_mpf", n_eta-1,eta_bins);
    TH1D* ptave_logpt_MPF = new TH1D("ptave_logpt_mpf","ptave_logpt_mpf", n_eta-1,eta_bins);
    ofstream output, output_loglin, uncerts, uncerts_loglin;
    output.open(path+"output/Fall15_25ns_MPF_FLAT_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt");
    output_loglin.open(path+"output/Fall15_25ns_MPF_LOGLIN_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt");
    uncerts.open(path+"output/Fall15_25ns_MPF_FLAT_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt.STAT");
    uncerts_loglin.open(path+"output/Fall15_25ns_MPF_LOGLIN_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt.STAT");
    output  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(208))/208))) Correction L2Relative}" << endl;
    output_loglin  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(208))/208))) Correction L2Relative}" << endl;
    uncerts << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}" << endl;
    uncerts_loglin << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}" << endl;
    for (int j=n_eta-1; j>0; --j){
      output << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   55   " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "   " << 1/flat_norm << " " << hist_kfsr_mpf->GetBinContent(j) << "   " << f2[j-1]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;
      output_loglin << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   55   " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "   " << 1/loglin_norm << " " << hist_kfsr_mpf->GetBinContent(j) << "   " << f1[j-1]->GetParameter(0) << " " << f1[j-1]->GetParameter(1) << "   1 0.0000 0.0" << endl;
      uncerts << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "       " << sqrt( pow(f2[j-1]->GetParameter(0)*hist_kfsr_mpf->GetBinError(j),2)+pow( hist_kfsr_mpf->GetBinContent(j)*f2[j-1]->GetParError(0),2)  ) / flat_norm  << endl;
      uncerts_loglin << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "       " << sqrt(  pow(f1[j-1]->GetParameter(0)+f1[j-1]->GetParameter(1)*TMath::Log(ptave_data[j-1]->GetMean()),2)*pow(hist_kfsr_mpf->GetBinError(j),2)   + pow(hist_kfsr_mpf->GetBinContent(j),2)*(pow(f2[j-1]->GetParError(0),2)+pow(f1[j-1]->GetParError(1),2)) ) / loglin_norm << endl;
    }
  
    for (int j=0; j<n_eta-1; j++){
      output << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   55   " << ptave_data[j]->FindLastBinAbove(0.)*10 << "   " << 1/flat_norm << " " << hist_kfsr_mpf->GetBinContent(j+1) << "   " << f2[j]->GetParameter(0) << " 0   1 0.0000 0.0"  << endl;
      output_loglin << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   55   " << ptave_data[j]->FindLastBinAbove(0.)*10 << "   " << 1/loglin_norm << " " << hist_kfsr_mpf->GetBinContent(j+1) << "   " << f1[j]->GetParameter(0) << " " << f1[j]->GetParameter(1) << "   1 0.0000 0.0" << endl;
      uncerts << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " << ptave_data[j]->FindLastBinAbove(0.)*10 << "       " << sqrt( pow(f2[j]->GetParameter(0)*hist_kfsr_mpf->GetBinError(j+1),2)+pow( hist_kfsr_mpf->GetBinContent(j+1)*f2[j]->GetParError(0),2)  ) / flat_norm  << endl;
      uncerts_loglin << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " << ptave_data[j]->FindLastBinAbove(0.)*10 << "       " << sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()),2)*pow(hist_kfsr_mpf->GetBinError(j+1),2)   + pow(hist_kfsr_mpf->GetBinContent(j+1),2)*(pow(f2[j]->GetParError(0),2)+pow(f1[j]->GetParError(1),2)) ) / loglin_norm << endl;

      if(syst==5){
	Residual_logpt_MPF->SetBinContent(j+1,hist_kfsr_mpf->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()) ) );
	Residual_logpt_MPF->SetBinError(j+1, sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()),2)*pow(hist_kfsr_mpf->GetBinError(j+1),2)   + pow(hist_kfsr_mpf->GetBinContent(j+1),2)*(pow(f2[j]->GetParError(0),2)+pow(f1[j]->GetParError(1),2)) )  );
      }
      if(syst==0){
	Residual_logpt_MPF->SetBinContent(j+1,hist_kfsr_mpf->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(120) ) );
	Residual_logpt_MPF->SetBinError(j+1, sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(120),2)*pow(hist_kfsr_mpf->GetBinError(j+1),2)   + pow(hist_kfsr_mpf->GetBinContent(j+1),2)*(pow(f2[j]->GetParError(0),2)+pow(f1[j]->GetParError(1),2)) )  );
      }
      if(syst==1){
	Residual_logpt_MPF->SetBinContent(j+1,hist_kfsr_mpf->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(60) ) );
      }
      if(syst==2){
	Residual_logpt_MPF->SetBinContent(j+1,hist_kfsr_mpf->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(240) ) );
      }
      if(syst==3){
	Residual_logpt_MPF->SetBinContent(j+1,hist_kfsr_mpf->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(480) ) );
      }

      Residual_const_MPF->SetBinContent(j+1,hist_kfsr_mpf->GetBinContent(j+1)*f2[j]->GetParameter(0));
      Residual_const_MPF->SetBinError(j+1, sqrt( pow(f2[j]->GetParameter(0)*hist_kfsr_mpf->GetBinError(j+1),2)+pow( hist_kfsr_mpf->GetBinContent(j+1)*f2[j]->GetParError(0),2)  )  );
      ptave_const_MPF->SetBinContent(j+1,f2[j]->GetParameter(0));
      ptave_const_MPF->SetBinError(j+1,f2[j]->GetParError(0));
      ptave_logpt_MPF->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()));
      ptave_logpt_MPF->SetBinError(j+1,f2[j]->GetParError(0));
    }




    TFile* outputfile;
    if(syst==5){
      outputfile = new TFile(path+"Histo_Res_MPF_L1"+tag+".root","RECREATE");
    }
    if(syst==0){
     outputfile = new TFile(path+"Histo_Res_MPF_L1_"+variation+".root","RECREATE");
    }
    if(syst==1){
     outputfile = new TFile(path+"Histo_Res_MPF_L1_"+variation+".root","RECREATE");
    }
    if(syst==2){
      outputfile = new TFile(path+"Histo_Res_MPF_L1_"+variation+".root","RECREATE");
    }
    if(syst==3){
      outputfile = new TFile(path+"Histo_Res_MPF_L1_"+variation+".root","RECREATE");
    }

    Residual_logpt_MPF->Write();
    Residual_const_MPF->Write();
    outputfile->Write();
    outputfile->Close();
    TFile* outputfile2 = new TFile(path+"Histo_ptave_MPF_L1.root","RECREATE");
    ptave_const_MPF->Write();
    ptave_logpt_MPF->Write();
    outputfile2->Write();
    outputfile2->Close();
  }



  //==============================================================================================================================
  // DIJET balance
  else{
    TFile* kfsr_dijet = new TFile(path+"Histo_KFSR_DiJet_L1.root","READ");
    TH1D* hist_kfsr_dijet = (TH1D*)kfsr_dijet->Get("kfsr_dijet");
    /*
    for(int i=0; i<50; i++){
      hist_kfsr_dijet->SetBinContent(i,1);
    }
    */
    double flat[n_eta-1];
    double loglin[n_eta-1];
    for (int j=0; j<n_eta-1; j++){
      flat[j] = hist_kfsr_dijet->GetBinContent(j+1)*f2[j]->GetParameter(0);
      loglin[j] = hist_kfsr_dijet->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()) );
    }
    double flat_norm, loglin_norm;
    for (int j=0; j<n_etabarr; j++){
      flat_norm += flat[j];
      loglin_norm += loglin[j];
    }
    flat_norm = flat_norm/n_etabarr;
    loglin_norm = loglin_norm/n_etabarr;
    for (int j=0; j<n_eta-1; j++){
      flat[j] = flat[j]/flat_norm;
      loglin[j] = loglin[j]/loglin_norm;
      //cout << flat[j] << endl;
    }

    cout << "normalization factor (flat): " << flat_norm << endl;
    cout << "normalization factor (logl): " << loglin_norm << endl;


    TH1D* Residual_logpt_DiJet = new TH1D("res_logpt_dijet","res_logpt_dijet", n_eta-1,eta_bins);
    TH1D* Residual_const_DiJet = new TH1D("res_const_dijet","res_const_dijet", n_eta-1,eta_bins);
    TH1D* ptave_const_DiJet = new TH1D("ptave_const_dijet","ptave_const_dijet", n_eta-1,eta_bins);
    TH1D* ptave_logpt_DiJet = new TH1D("ptave_logpt_dijet","ptave_logpt_dijet", n_eta-1,eta_bins);
    ofstream output, output_loglin, uncerts, uncerts_loglin;
    output.open(path+"output/Fall15_25ns_pT_FLAT_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt");
    output_loglin.open(path+"output/Fall15_25ns_pT_LOGLIN_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt");
    uncerts.open(path+"output/Fall15_25ns_pT_FLAT_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt.STAT");
    uncerts_loglin.open(path+"output/Fall15_25ns_pT_LOGLIN_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt.STAT");
    output  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(208))/208))) Correction L2Relative}" << endl;
    output_loglin  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(208))/208))) Correction L2Relative}" << endl;
    uncerts << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}" << endl;
    uncerts_loglin << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}" << endl;

    for (int j=n_eta-1; j>0; --j){
      output << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   55   " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "   " << 1/flat_norm << " " << hist_kfsr_dijet->GetBinContent(j) << "   " << f2[j-1]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;
      output_loglin << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   55   " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "   " << 1/loglin_norm << " " << hist_kfsr_dijet->GetBinContent(j) << "   " << f1[j-1]->GetParameter(0) << " " << f1[j-1]->GetParameter(1) <<  "   1 0.0000 0.0" << endl;
      uncerts << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "       " << sqrt( pow(f2[j-1]->GetParameter(0)*hist_kfsr_dijet->GetBinError(j),2)+pow( hist_kfsr_dijet->GetBinContent(j)*f2[j-1]->GetParError(0),2)  ) / flat_norm  << endl;
      uncerts_loglin << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "       " << sqrt(  pow(f1[j-1]->GetParameter(0)+f1[j-1]->GetParameter(1)*TMath::Log(ptave_data[j-1]->GetMean()),2)*pow(hist_kfsr_dijet->GetBinError(j),2)   + pow(hist_kfsr_dijet->GetBinContent(j),2)*(pow(f2[j-1]->GetParError(0),2)+pow(f1[j-1]->GetParError(1),2)) ) / loglin_norm << endl;
    }


    for (int j=0; j<n_eta-1; j++){
      
      output << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   55   " << ptave_data[j]->FindLastBinAbove(0.)*10 << "   " << 1/flat_norm << " " << hist_kfsr_dijet->GetBinContent(j+1) << "   " << f2[j]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;
      output_loglin << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   55   " << ptave_data[j]->FindLastBinAbove(0.)*10 << "   " << 1/loglin_norm << " " << hist_kfsr_dijet->GetBinContent(j+1) << "   " << f1[j]->GetParameter(0) << " " << f1[j]->GetParameter(1) <<  "   1 0.0000 0.0" << endl;
      uncerts << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " << ptave_data[j]->FindLastBinAbove(0.)*10 << "       " << sqrt( pow(f2[j]->GetParameter(0)*hist_kfsr_dijet->GetBinError(j+1),2)+pow( hist_kfsr_dijet->GetBinContent(j+1)*f2[j]->GetParError(0),2)  ) / flat_norm  << endl;
      uncerts_loglin << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " << ptave_data[j]->FindLastBinAbove(0.)*10 << "       " << sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()),2)*pow(hist_kfsr_dijet->GetBinError(j+1),2)   + pow(hist_kfsr_dijet->GetBinContent(j+1),2)*(pow(f2[j]->GetParError(0),2)+pow(f1[j]->GetParError(1),2)) ) / loglin_norm  << endl;


      if(syst==5){
	Residual_logpt_DiJet->SetBinContent(j+1,hist_kfsr_dijet->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()) ) );
	Residual_logpt_DiJet->SetBinError(j+1, sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()),2)*pow(hist_kfsr_dijet->GetBinError(j+1),2)   + pow(hist_kfsr_dijet->GetBinContent(j+1),2)*(pow(f2[j]->GetParError(0),2)+pow(f1[j]->GetParError(1),2)) )  );
      }
      if(syst==0){
	Residual_logpt_DiJet->SetBinContent(j+1,hist_kfsr_dijet->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(120) ) );
	Residual_logpt_DiJet->SetBinError(j+1, sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(120),2)*pow(hist_kfsr_dijet->GetBinError(j+1),2)   + pow(hist_kfsr_dijet->GetBinContent(j+1),2)*(pow(f2[j]->GetParError(0),2)+pow(f1[j]->GetParError(1),2)) )  );
      }
      if(syst==1){
	Residual_logpt_DiJet->SetBinContent(j+1,hist_kfsr_dijet->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(60) ) );
      }
      if(syst==2){
	Residual_logpt_DiJet->SetBinContent(j+1,hist_kfsr_dijet->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(240) ) );
      }
      if(syst==3){
	Residual_logpt_DiJet->SetBinContent(j+1,hist_kfsr_dijet->GetBinContent(j+1)*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(480) ) );
      }


      Residual_const_DiJet->SetBinContent(j+1,hist_kfsr_dijet->GetBinContent(j+1)*f2[j]->GetParameter(0));
      Residual_const_DiJet->SetBinError(j+1, sqrt( pow(f2[j]->GetParameter(0)*hist_kfsr_dijet->GetBinError(j+1),2)+pow( hist_kfsr_dijet->GetBinContent(j+1)*f2[j]->GetParError(0),2)  )  );
      ptave_const_DiJet->SetBinContent(j+1,f2[j]->GetParameter(0));
      ptave_const_DiJet->SetBinError(j+1,f2[j]->GetParError(0));
      ptave_logpt_DiJet->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()));
      ptave_logpt_DiJet->SetBinError(j+1,f2[j]->GetParError(0));
    }





    TFile* outputfile;
    if(syst==5){
      outputfile = new TFile(path+"Histo_Res_DiJet_L1"+tag+".root","RECREATE");
    }
    if(syst==0){
     outputfile = new TFile(path+"Histo_Res_DiJet_L1_"+variation+".root","RECREATE");
    }
    if(syst==1){
      outputfile = new TFile(path+"Histo_Res_DiJet_L1_"+variation+".root","RECREATE");
    }
    if(syst==2){
      outputfile = new TFile(path+"Histo_Res_DiJet_L1_"+variation+".root","RECREATE");
    }
    if(syst==3){
      outputfile = new TFile(path+"Histo_Res_DiJet_L1_"+variation+".root","RECREATE");
    }

    Residual_logpt_DiJet->Write();
    Residual_const_DiJet->Write();
    outputfile->Write();
    outputfile->Close();
    TFile* outputfile2 = new TFile(path+"Histo_ptave_DiJet_L1.root","RECREATE");
    ptave_const_DiJet->Write();
    ptave_logpt_DiJet->Write();
    outputfile2->Write();
    outputfile2->Close();
  }






}
