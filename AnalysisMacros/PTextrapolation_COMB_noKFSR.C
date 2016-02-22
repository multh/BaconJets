// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// with this script pT extrapolations in bins of eta will be calculated
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "header.h"



void PTextrapolation_COMB_noKFSR(TString path, TFile* datafile, TFile* MCfile, TString txttag, TString jettag, TString variation, TString tag){

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
  syst=0;//ToDo: apply variations

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
  double ratio_dijet[n_eta-1][n_pt-1]; //ratio at pt,eta bins for alpha = 0.3
  double err_ratio_dijet[n_eta-1][n_pt-1]; //error of ratio at pt,eta bins for alpha = 0.3
  for(int j=0; j<n_eta-1; j++){
    for(int k=0; k<n_pt-1; k++){
      std::cout<<"For eta range = "<<eta_bins[j]<<", "<<eta_bins[j+1]<<" and pT range = "<<pt_bins[k]<<", "<<pt_bins[k+1]<<std::endl;
      pair<double,double> res_mc_mpf =  Response(MCfile,0.3,eta_bins[j],eta_bins[j+1],pt_bins[k],pt_bins[k+1],true);
      pair<double,double> res_data_mpf =  Response(datafile,0.3,eta_bins[j],eta_bins[j+1],pt_bins[k],pt_bins[k+1],true);
      pair<double,double> ratio_res_mpf = Rmc_to_Rdata(res_mc_mpf,res_data_mpf);
      ratio_mpf[j][k] = ratio_res_mpf.first;
      err_ratio_mpf[j][k] = ratio_res_mpf.second;
      pair<double,double> res_mc_dijet =  Response(MCfile,0.3,eta_bins[j],eta_bins[j+1],pt_bins[k],pt_bins[k+1],false);
      pair<double,double> res_data_dijet =  Response(datafile,0.3,eta_bins[j],eta_bins[j+1],pt_bins[k],pt_bins[k+1],false);
      pair<double,double> ratio_res_dijet = Rmc_to_Rdata(res_mc_dijet,res_data_dijet);
      ratio_dijet[j][k] = ratio_res_dijet.first;
      err_ratio_dijet[j][k] = ratio_res_dijet.second;
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
 
  for(int j=0; j<n_eta-1; j++){
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
  fp = fopen("output/pT_COMB_extrapolations.dat","w");
  fp2 = fopen("output/pT_COMB_constantExtrapolation.dat","w");
  l2resfile = fopen("output/L2Res_COMB.dat","w");

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
  //for (int j=n_eta-2; j<n_eta-1; j++){//TEST
    plotname[j]="comb_ptextra_eta_"+eta_range[j]+"_"+eta_range[j+1];
    asd[j] = new TCanvas(plotname[j], plotname[j], 800,700);
    gStyle->SetOptTitle(0);
    gPad->SetLogx();
   
    //    f1[j] = new TF1(plotname[j]+"f1","[0]+[1]*TMath::Log(x)", 76 , 558);
    f1[j] = new TF1(plotname[j]+"f1","[0]+[1]*TMath::Log(x)", 0 , 600);
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
    if(graph1_mpf[j]->GetN()>0) mg1[j]->Add(graph1_mpf[j]);
    if(graph1_dijet[j]->GetN()>0) mg1[j]->Add(graph1_dijet[j]);

    mg1[j]->Fit(plotname[j]+"f1","R");
    mg1[j]->Fit(plotname[j]+"f2","R");
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
    // TPaveStats *st = ((TPaveStats*)(graph1[j]->GetListOfFunctions()->FindObject("stats")));
    // if (st) {
    //   st->SetTextColor(kRed);
    //   st->SetX1NDC(0.69); st->SetX2NDC(0.99);
    //   st->SetY1NDC(0.65); st->SetY2NDC(0.8);
    // }
    // st = ((TPaveStats*)(graph2[j]->GetListOfFunctions()->FindObject("stats")));
    // if (st) {
    //   st->SetTextColor(graph2[j]->GetLineColor());
    //   st->SetX1NDC(0.69); st->SetX2NDC(0.99);
    //   st->SetY1NDC(0.85); st->SetY2NDC(0.95);
    // }
    // asd[j]->Modified(); asd[j]->Update();



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


    mg1[j]->GetYaxis()->SetRangeUser(0.95,1.10);
    // if(j<11){
    //   mg1[j]->GetYaxis()->SetRangeUser(0.95,1.10);
    // }
    // if(j>=11){
    //   mg1[j]->GetYaxis()->SetRangeUser(0.95,1.30);
    // }
    //save plots
    asd[j]->Print("plots/pTextrapolation_COMB_pT_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
  }
  fclose(fp);
  fclose(fp2);


  
  //==============================================================================================================================
  // get the kFSR plots for mpfMethod and calculate residuals
  // TFile* kfsr_comb = new TFile(path+"Histo_KFSR_MPF_L1.root","READ");
  // TH1D* hist_kfsr_comb = (TH1D*)kfsr_comb->Get("kfsr_mpf");
  // TFile* kfsr_comb = new TFile(path+"Histo_KFSR_DiJet_L1.root","READ");
  // TH1D* hist_kfsr_comb = (TH1D*)kfsr_comb->Get("kfsr_dijet");
  double flat[n_eta-1];
  double loglin[n_eta-1];
  for (int j=0; j<n_eta-1; j++){
    flat[j] = f2[j]->GetParameter(0);
    loglin[j] = f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean());
    cout<<"ptave_data[j]->GetMean() = "<<ptave_data[j]->GetMean()<<endl;
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
  }

  cout << "normalization factor (flat): " << flat_norm << endl;
  cout << "normalization factor (logl): " << loglin_norm << endl;
  
  TH1D* Residual_logpt_COMB = new TH1D("res_logpt_comb","res_logpt_comb", n_eta-1,eta_bins);
  TH1D* Residual_const_COMB = new TH1D("res_const_comb","res_const_comb", n_eta-1,eta_bins);
  TH1D* ptave_const_COMB = new TH1D("ptave_const_comb","ptave_const_comb", n_eta-1,eta_bins);
  TH1D* ptave_logpt_COMB = new TH1D("ptave_logpt_comb","ptave_logpt_comb", n_eta-1,eta_bins);
  ofstream output, output_loglin, uncerts, uncerts_loglin;
  output.open("output/Fall15_25ns_COMB_FLAT_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt");
  output_loglin.open("output/Fall15_25ns_COMB_LOGLIN_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt");
  uncerts.open("output/Fall15_25ns_COMB_FLAT_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt.STAT");
  uncerts_loglin.open("output/Fall15_25ns_COMB_LOGLIN_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt.STAT");
  //output  << "JetEta (abs)      " << "p[0]" << "       " <<  "kFSR * p[0] " << "  statistical unc" << endl;
  //output_loglin  << "JetEta (abs)      "  << "p[0]        "  << "p[1]"  << "       " <<  "kFSR*[p[0]+p[1]*Log(max(ptmin,min(ptmax, x)) )] " << "  statistical unc" << endl;
  output  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(208))/208))) Correction L2Relative}" << endl;
  output_loglin  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(208))/208))) Correction L2Relative}" << endl;
  uncerts << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}" << endl;
  uncerts_loglin << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}" << endl;
  
  for (int j=n_eta-1; j>0; --j){
    output << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   55   " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "   " << 1/flat_norm << " " << 1.0 << "   " << f2[j-1]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;
    output_loglin << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   55   " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "   " << 1/loglin_norm << " " << 1.0 << "   " << f1[j-1]->GetParameter(0) << " " << f1[j-1]->GetParameter(1) << "   1 0.0000 0.0" << endl;
    uncerts << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "       " << sqrt( pow(f2[j-1]->GetParameter(0)*1,2)+pow( 1.0*f2[j-1]->GetParError(0),2)  ) / flat_norm  << endl;
    uncerts_loglin << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " << ptave_data[j-1]->FindLastBinAbove(0.)*10 << "       " << sqrt(  pow(f1[j-1]->GetParameter(0)+f1[j-1]->GetParameter(1)*TMath::Log(ptave_data[j-1]->GetMean()),2)*pow(1,2)   + pow(1.0,2)*(pow(f2[j-1]->GetParError(0),2)+pow(f1[j-1]->GetParError(1),2)) ) / loglin_norm << endl;
  }
  
  for (int j=0; j<n_eta-1; j++){
    output << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   55   " << ptave_data[j]->FindLastBinAbove(0.)*10 << "   " << 1/flat_norm << " " << 1.0 << "   " << f2[j]->GetParameter(0) << " 0   1 0.0000 0.0"  << endl;
      output_loglin << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   55   " << ptave_data[j]->FindLastBinAbove(0.)*10 << "   " << 1/loglin_norm << " " << 1.0 << "   " << f1[j]->GetParameter(0) << " " << f1[j]->GetParameter(1) << "   1 0.0000 0.0" << endl;
      uncerts << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " << ptave_data[j]->FindLastBinAbove(0.)*10 << "       " << sqrt( pow(f2[j]->GetParameter(0)*1,2)+pow( 1.0*f2[j]->GetParError(0),2)  ) / flat_norm  << endl;
      uncerts_loglin << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " << ptave_data[j]->FindLastBinAbove(0.)*10 << "       " << sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()),2)*pow(1,2)   + pow(1.0,2)*(pow(f2[j]->GetParError(0),2)+pow(f1[j]->GetParError(1),2)) ) / loglin_norm << endl;


      // if(syst==5){
      // 	Residual_logpt_COMB->SetBinContent(j+1,1.0*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()) ) );
      // 	Residual_logpt_COMB->SetBinError(j+1, sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()),2)*pow(1,2)   + pow(1.0,2)*(pow(f2[j]->GetParError(0),2)+pow(f1[j]->GetParError(1),2)) )  );
      // }
      if(syst==0){
  	Residual_logpt_COMB->SetBinContent(j+1,(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean())));
       	cout<<" LogLin Residual in bin#"<<j+1<<" "<<(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()))<<" +/- "<<f2[j]->GetParError(0)<<endl;
  	Residual_logpt_COMB->SetBinError(j+1,f2[j]->GetParError(0));
      }
      else{
	cout<<"This is not properly implemented yet!"<<endl;
      }

      // if(syst==1){
      // 	Residual_logpt_COMB->SetBinContent(j+1,1.0*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(60) ) );
      // }
      // if(syst==2){
      // 	Residual_logpt_COMB->SetBinContent(j+1,1.0*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(240) ) );
      // }
      // if(syst==3){
      // 	Residual_logpt_COMB->SetBinContent(j+1,1.0*(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(480) ) );
      // }

      Residual_const_COMB->SetBinContent(j+1,f2[j]->GetParameter(0));
      cout<<"Const Residual in bin#"<<j+1<<" "<<f2[j]->GetParameter(0)<<endl;
      Residual_const_COMB->SetBinError(j+1,f2[j]->GetParError(0));
      ptave_const_COMB->SetBinContent(j+1,f2[j]->GetParameter(0));
      ptave_const_COMB->SetBinError(j+1,f2[j]->GetParError(0));
      ptave_logpt_COMB->SetBinContent(j+1,f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()));
      ptave_logpt_COMB->SetBinError(j+1,f2[j]->GetParError(0));
    }




    TFile* outputfile;
    outputfile = new TFile(path+"Histo_Res_COMB_L1.root","RECREATE");

    // if(syst==5){
    //   outputfile = new TFile(path+"Histo_Res_COMB_L1"+tag+".root","RECREATE");
    // }
    // if(syst==0){
    //  outputfile = new TFile(path+"Histo_Res_COMB_L1_"+variation+".root","RECREATE");
    // }
    // if(syst==1){
    //  outputfile = new TFile(path+"Histo_Res_COMB_L1_"+variation+".root","RECREATE");
    // }
    // if(syst==2){
    //   outputfile = new TFile(path+"Histo_Res_COMB_L1_"+variation+".root","RECREATE");
    // }
    // if(syst==3){
    //   outputfile = new TFile(path+"Histo_Res_COMB_L1_"+variation+".root","RECREATE");
    // }

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
