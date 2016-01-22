// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// with this script pT extrapolations in bins of eta will be calculated
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "header.h"


TString ToString(int num) {
  ostringstream start;
  start<<num;
  TString start1=start.str();
  return start1;
}

void PTextrapolation(bool mpfMethod(false), TString path, TFile* datafile, TFile* MCfile, TString txttag, TString jettag, TString variation, TString tag){
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


  // get the histos 
  TH1D* data[n_eta-1][n_pt-1];
  TH1D* mc[n_eta-1][n_pt-1];
  TH1D* ptave_data[n_eta-1];
  TH1D* ptaverebin_data[n_eta-1];
  
  for(int i=0; i<n_eta-1; i++){
    ptave_data[i] = (TH1D*)datafile->Get("/a02/eta_"+eta_range[i]+"_"+eta_range[i+1]+"/pt_ave");
    ptaverebin_data[i] = (TH1D*)datafile->Get("/a02/eta_"+eta_range[i]+"_"+eta_range[i+1]+"/pt_ave_rebin");
    for(int j=0; j<n_pt-1; j++){
      if(mpfMethod){
	data[i][j] = (TH1D*)datafile->Get("/eta_"+eta_range[i]+"_"+eta_range[i+1]+"/pt_"+pt_range[j]+"_"+pt_range[j+1]+"/mpf");
	mc[i][j] = (TH1D*)MCfile->Get("/eta_"+eta_range[i]+"_"+eta_range[i+1]+"/pt_"+pt_range[j]+"_"+pt_range[j+1]+"/mpf");
      }
      else{
	data[i][j] = (TH1D*)datafile->Get("/eta_"+eta_range[i]+"_"+eta_range[i+1]+"/pt_"+pt_range[j]+"_"+pt_range[j+1]+"/r_rel");
	mc[i][j] = (TH1D*)MCfile->Get("/eta_"+eta_range[i]+"_"+eta_range[i+1]+"/pt_"+pt_range[j]+"_"+pt_range[j+1]+"/r_rel");
      }
    }
  }

 
  
  // create the ratio histogram
  TH1D* ratio[n_eta-1];
  //TF1* fit_func[n_eta-1];
  for (int j=0; j<n_eta-1; j++){
    TString numstr=ToString(j);
    if(mpfMethod){
      TString histoname="mpf_ptextra_histo"+numstr;
      ratio[j] = new TH1D(histoname,histoname,n_pt-1,pt_bins);
    }
    else{
      TString histoname="kfsr_ptextra_histo"+numstr;
      ratio[j] = new TH1D(histoname,histoname,n_pt-1,pt_bins);
    }
    //fit_func[j] = new TF1("lin_fit[j]", "pol1", 0., n_alpha);
  }


  
  // fill ratio histograms
  for(int j=0; j<n_eta-1; j++){
    for(int i=0; i<n_pt-1; i++){
      ratio[j]->SetBinContent(i+1, (mc[j][i]->GetMean() / data[j][i]->GetMean()));
      ratio[j]->SetBinError(i+1, sqrt( pow(data[j][i]->GetRMS(),2)/data[j][i]->Integral())); // + pow(mc[j][i]->GetRMS(),2)/mc[j][i]->Integral() ) );
      //ratio[j]->SetBinError(i+1, sqrt( (mc[j][i]->GetMean()/pow(data[j][i]->GetMean(),2))* pow(data[j][i]->GetRMS(),2)/data[j][i]->Integral() + (1/data[j][i]->GetMean())* pow(mc[j][i]->GetRMS(),2)/mc[j][i]->Integral() ) );
      //ratio[j]->SetBinError(i+1, 1/sqrt(data[j][i]->Integral()));
      //ratio[j]->SetBinError(i+1, 1/data[j][i]->Integral());
    }
  }
 

  // create and fill tgrapherrors
  double xbin_tgraph[n_pt-1];
  double zero[n_pt-1];
  for(int i=0;i<n_pt-1;i++){
    xbin_tgraph[i]=(pt_bins[i]+pt_bins[i+1])/2;
    zero[i]=(pt_bins[i+1]-pt_bins[i])/2 ;
  }
  TGraphErrors *graph1[n_eta-1];
  TGraphErrors *graph2[n_eta-1];
  double content[n_pt-1];
  double error[n_pt-1];
  for(int j=0; j<n_eta-1; j++){
    for(int i=0; i<n_pt-1; i++){
      content[i] = ratio[j]->GetBinContent(i+1);
      error[i] = ratio[j]->GetBinError(i+1);
    }
    graph1[j] = new TGraphErrors(n_pt-1, xbin_tgraph, content , zero, error);
    graph2[j] = new TGraphErrors(n_pt-1, xbin_tgraph, content , zero, error);
  }

  
  // clean up nonphysical values
  for(int i=0; i<n_eta-1; i++){
    for(int j= graph1[i]->GetN()-1; j != -1; --j) {
      if(graph1[i]->GetY()[j]!=graph1[i]->GetY()[j] || graph1[i]->GetEY()[j]!=graph1[i]->GetEY()[j]) graph1[i]->RemovePoint(j);
    }
  }
  

  // remove all data points with less than 100 entries  
  for(int j=0; j<n_eta-1; j++){
    for(int i=n_pt-2; i!=-1; --i){
      if(data[j][i]->GetEntries()<100){
	graph1[j]->RemovePoint(i);
      }
    }
  }
  // make sure that the cut was working
  for(int i=0; i<n_eta-1; i++){
      cout << graph1[i]->GetN() << endl;
  }
  
  for(int j=0; j<n_eta-1; j++){
    graph2[j] =  (TGraphErrors*)graph1[j]->Clone();	
  }


  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,pt_bins[8]+20,1);

  // create a function for the loglinear fit
  TF1 * f1[n_eta-1];
  // create a function for the constant fit
  TF1 * f2[n_eta-1];

  // create output .dat file, including the kFSR extrapolation (alpha->0)
  if(mpfMethod){
    FILE *fp = fopen("pT_MPF_extrapolations.dat","w");
    FILE *fp2 = fopen("pT_MPF_constantExtrapolation.dat","w");
    FILE *l2resfile = fopen("L2Res_MPF.dat","w");
  }
  else{
    FILE *fp = fopen("pT_DiJet_extrapolations.dat","w");
    FILE *fp2 = fopen("pT_DiJet_constantExtrapolation.dat","w");
    FILE *l2resfile = fopen("L2Res_DiJet.dat","w");
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
    graph1[j]->SetMarkerColor(kBlue);
    graph1[j]->SetMarkerStyle(20);
    graph1[j]->SetLineColor(kBlue);
    graph2[j]->SetMarkerColor(kBlue);
    graph2[j]->SetMarkerStyle(20);
    graph2[j]->SetLineColor(kBlue);
 
    f1[j] = new TF1(plotname[j]+"f1","[0]+[1]*TMath::Log(x)", 76 , 558);
    f1[j]->SetParameters(1,0);
    f2[j] = new TF1(plotname[j]+"f2","pol0");
    f2[j]->SetLineColor(kBlue);
    f2[j]->SetLineStyle(3);

    graph1[j]->Fit(plotname[j]+"f1","");
    graph1[j]->Fit(plotname[j]+"f2","+ SAME");
    graph2[j]->Fit(plotname[j]+"f2","");
    graph1[j]->Draw("AP");
    graph2[j]->Draw("P SAME");
    graph1[j]->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph1[j]->GetXaxis()->SetTitleSize(0.05);
    graph1[j]->GetXaxis()->SetTitleOffset(0.80);
    //graph1[j]->GetYaxis()->SetTitle("kFSR");
    graph1[j]->GetXaxis()->SetLimits(30,pt_bins[8]+20);
    if(j<11){
      graph1[j]->GetYaxis()->SetRangeUser(0.95,1.10);
    }
    if(j>=11){
      graph1[j]->GetYaxis()->SetRangeUser(0.95,1.30);
    }

 
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
    TPaveStats *st = ((TPaveStats*)(graph1[j]->GetListOfFunctions()->FindObject("stats")));
    if (st) {
      st->SetTextColor(kRed);
      st->SetX1NDC(0.69); st->SetX2NDC(0.99);
      st->SetY1NDC(0.65); st->SetY2NDC(0.8);
    }
    st = ((TPaveStats*)(graph2[j]->GetListOfFunctions()->FindObject("stats")));
    if (st) {
      st->SetTextColor(graph2[j]->GetLineColor());
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
      leg1->SetHeader("MPF #bar{p}_{T} extrapolation, "+eta_range3[j]+"#leq|#eta|<"+eta_range3[j+1]);
    }
    else{
      leg1->SetHeader("p_{T} balance #bar{p}_{T} extrapolation, "+eta_range3[j]+"#leq|#eta|<"+eta_range3[j+1]);
    }
    leg1->AddEntry(graph1[j], "R(MC)/R(DATA)","P");
    leg1->AddEntry(f1[j], "loglinear fit","L");
    leg1->AddEntry(f2[j], "constant fit","L");
    leg1->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    tex->DrawLatex(0.64,0.91,"2.11fb^{-1} (13TeV)");

    //save plots
    if(mpfMethod){
      asd[j]->Print("plots/pTextrapolation_MPF_pT_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
    }
    else{
      asd[j]->Print("plots/pTextrapolation_Pt_pT_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
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
    cout << "max value of pt " << eta_range[i] << " to " << eta_range[i+1] << ": " << ptaverebin_data[i]->FindLastBinAbove(0.) << endl;
  }

  for(int i=0; i<n_eta-1; i++){
    for(int j=0; j<graph1[i]->GetN(); j++) {
      cout << "uncert, eta " << eta_range[i] << " to " << eta_range[i+1] << ", pt range " << pt_range[j] << " to " << pt_range[j+1] << ": "  << graph1[i]->GetEY()[j] << endl;
    }
  }


  //==============================================================================================================================
  // get the kFSR plots and calculate residuals
  if(mpfMethod){
    TFile* kfsr_mpf = new TFile(path+"Histo_KFSR_MPF_L1.root","READ");
    TH1D* hist_kfsr_mpf = kfsr_mpf->Get("kfsr_mpf");
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
    ofstream output, output_loglin, uncerts, uncerts_loglin;
    output.open("Summer15_25ns_MPF_FLAT_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt");
    output_loglin.open("Summer15_25ns_MPF_LOGLIN_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt");
    uncerts.open("Summer15_25ns_MPF_FLAT_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt.STAT");
    uncerts_loglin.open("Summer15_25ns_MPF_LOGLIN_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt.STAT");
    //output  << "JetEta (abs)      " << "p[0]" << "       " <<  "kFSR * p[0] " << "  statistical unc" << endl;
    //output_loglin  << "JetEta (abs)      "  << "p[0]        "  << "p[1]"  << "       " <<  "kFSR*[p[0]+p[1]*Log(max(ptmin,min(ptmax, x)) )] " << "  statistical unc" << endl;
    output  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(208))/208))) Correction L2Relative}" << endl;
    output_loglin  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(208))/208))) Correction L2Relative}" << endl;
    uncerts << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}" << endl;
    uncerts_loglin << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}" << endl;

    //TString ptmax[n_eta-1] = {"2160","2060","1950","1810","1970","1590","1380","1260"," 920"," 850"," 790"," 710"," 560"," 510"," 510"};
    for (int j=n_eta-1; j>0; --j){
      output << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   55   " << ptaverebin_data[j-1]->FindLastBinAbove(0.)*10 << "   " << 1/flat_norm << " " << hist_kfsr_mpf->GetBinContent(j) << "   " << f2[j-1]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;
      output_loglin << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   55   " << ptaverebin_data[j-1]->FindLastBinAbove(0.)*10 << "   " << 1/loglin_norm << " " << hist_kfsr_mpf->GetBinContent(j) << "   " << f1[j-1]->GetParameter(0) << " " << f1[j-1]->GetParameter(1) << "   1 0.0000 0.0" << endl;
      uncerts << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " << ptaverebin_data[j-1]->FindLastBinAbove(0.)*10 << "       " << sqrt( pow(f2[j-1]->GetParameter(0)*hist_kfsr_mpf->GetBinError(j),2)+pow( hist_kfsr_mpf->GetBinContent(j)*f2[j-1]->GetParError(0),2)  ) / flat_norm  << endl;
      uncerts_loglin << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " << ptaverebin_data[j-1]->FindLastBinAbove(0.)*10 << "       " << sqrt(  pow(f1[j-1]->GetParameter(0)+f1[j-1]->GetParameter(1)*TMath::Log(ptave_data[j-1]->GetMean()),2)*pow(hist_kfsr_mpf->GetBinError(j),2)   + pow(hist_kfsr_mpf->GetBinContent(j),2)*(pow(f2[j-1]->GetParError(0),2)+pow(f1[j-1]->GetParError(1),2)) ) / loglin_norm << endl;
    }
  
    for (int j=0; j<n_eta-1; j++){
      //output << std::setprecision(5)  << eta_range[j]<< "    " << eta_range[j+1] << "    " <<  f2[j]->GetParameter(0)  << "     " << flat[j] << "        " << sqrt( pow(f2[j]->GetParameter(0)*hist_kfsr_mpf->GetBinError(j+1),2)+pow( hist_kfsr_mpf->GetBinContent(j+1)*f2[j]->GetParError(0),2)  ) / flat_norm  << endl;
      //output_loglin << fixed << std::setprecision(5)  << eta_range[j]<< "    " << eta_range[j+1] << "    " <<  f1[j]->GetParameter(0)  << "     " << f1[j]->GetParameter(1) << "            " << fixed << std::setprecision(10) << loglin[j] << "        " << sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()),2)*pow(hist_kfsr_mpf->GetBinError(j+1),2)   + pow(hist_kfsr_mpf->GetBinContent(j+1),2)*(pow(f2[j]->GetParError(0),2)+pow(f1[j]->GetParError(1),2)) ) / loglin_norm  << endl;
      output << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   55   " << ptaverebin_data[j]->FindLastBinAbove(0.)*10 << "   " << 1/flat_norm << " " << hist_kfsr_mpf->GetBinContent(j+1) << "   " << f2[j]->GetParameter(0) << " 0   1 0.0000 0.0"  << endl;
      output_loglin << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   55   " << ptaverebin_data[j]->FindLastBinAbove(0.)*10 << "   " << 1/loglin_norm << " " << hist_kfsr_mpf->GetBinContent(j+1) << "   " << f1[j]->GetParameter(0) << " " << f1[j]->GetParameter(1) << "   1 0.0000 0.0" << endl;
      uncerts << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " << ptaverebin_data[j]->FindLastBinAbove(0.)*10 << "       " << sqrt( pow(f2[j]->GetParameter(0)*hist_kfsr_mpf->GetBinError(j+1),2)+pow( hist_kfsr_mpf->GetBinContent(j+1)*f2[j]->GetParError(0),2)  ) / flat_norm  << endl;
      uncerts_loglin << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " << ptaverebin_data[j]->FindLastBinAbove(0.)*10 << "       " << sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()),2)*pow(hist_kfsr_mpf->GetBinError(j+1),2)   + pow(hist_kfsr_mpf->GetBinContent(j+1),2)*(pow(f2[j]->GetParError(0),2)+pow(f1[j]->GetParError(1),2)) ) / loglin_norm << endl;


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
    }





    if(syst==5){
      TFile* outputfile = new TFile(path+"Histo_Res_MPF_L1"+tag+".root","RECREATE");
    }
    if(syst==0){
      TFile* outputfile = new TFile(path+"Histo_Res_MPF_L1_"+variation+".root","RECREATE");
    }
    if(syst==1){
      TFile* outputfile = new TFile(path+"Histo_Res_MPF_L1_"+variation+".root","RECREATE");
    }
    if(syst==2){
      TFile* outputfile = new TFile(path+"Histo_Res_MPF_L1_"+variation+".root","RECREATE");
    }
    if(syst==3){
      TFile* outputfile = new TFile(path+"Histo_Res_MPF_L1_"+variation+".root","RECREATE");
    }

    Residual_logpt_MPF->Write();
    Residual_const_MPF->Write();
    outputfile->Write();
    outputfile->Close();
    TFile* outputfile2 = new TFile(path+"Histo_ptaveConst_MPF_L1.root","RECREATE");
    ptave_const_MPF->Write();
    outputfile2->Write();
    outputfile2->Close();
  }



  //==============================================================================================================================
  // DIJET balance
  else{
    TFile* kfsr_dijet = new TFile(path+"Histo_KFSR_DiJet_L1.root","READ");
    TH1D* hist_kfsr_dijet = kfsr_dijet->Get("kfsr_dijet");
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
    ofstream output, output_loglin, uncerts, uncerts_loglin;
    output.open("Summer15_25ns_pT_FLAT_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt");
    output_loglin.open("Summer15_25ns_pT_LOGLIN_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt");
    uncerts.open("Summer15_25ns_pT_FLAT_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt.STAT");
    uncerts_loglin.open("Summer15_25ns_pT_LOGLIN_L2Residual_"+txttag+"_"+jettag+"_"+variation+".txt.STAT");
    //output  << "JetEta (abs)      " << "p[0]" << "       " <<  "kFSR * p[0] " << "  statistical unc" << endl;
    //output_loglin  << "JetEta (abs)      "  << "p[0]        "  << "p[1]"  << "       " <<  "kFSR*[p[0]+p[1]*Log(max(ptmin,min(ptmax, x)) )] " << "  statistical unc" << endl;

    output  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(208))/208))) Correction L2Relative}" << endl;
    output_loglin  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(208))/208))) Correction L2Relative}" << endl;
    uncerts << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}" << endl;
    uncerts_loglin << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}" << endl;

    //TString ptmax[n_eta-1] = {"2160","2060","1950","1810","1970","1590","1380","1260"," 920"," 850"," 790"," 710"," 560"," 510"," 510"};
    for (int j=n_eta-1; j>0; --j){
      output << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   55   " << ptaverebin_data[j-1]->FindLastBinAbove(0.)*10 << "   " << 1/flat_norm << " " << hist_kfsr_dijet->GetBinContent(j) << "   " << f2[j-1]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;
      output_loglin << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   55   " << ptaverebin_data[j-1]->FindLastBinAbove(0.)*10 << "   " << 1/loglin_norm << " " << hist_kfsr_dijet->GetBinContent(j) << "   " << f1[j-1]->GetParameter(0) << " " << f1[j-1]->GetParameter(1) <<  "   1 0.0000 0.0" << endl;
      uncerts << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " << ptaverebin_data[j-1]->FindLastBinAbove(0.)*10 << "       " << sqrt( pow(f2[j-1]->GetParameter(0)*hist_kfsr_dijet->GetBinError(j),2)+pow( hist_kfsr_dijet->GetBinContent(j)*f2[j-1]->GetParError(0),2)  ) / flat_norm  << endl;
      uncerts_loglin << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " << ptaverebin_data[j-1]->FindLastBinAbove(0.)*10 << "       " << sqrt(  pow(f1[j-1]->GetParameter(0)+f1[j-1]->GetParameter(1)*TMath::Log(ptave_data[j-1]->GetMean()),2)*pow(hist_kfsr_dijet->GetBinError(j),2)   + pow(hist_kfsr_dijet->GetBinContent(j),2)*(pow(f2[j-1]->GetParError(0),2)+pow(f1[j-1]->GetParError(1),2)) ) / loglin_norm << endl;
    }


    for (int j=0; j<n_eta-1; j++){
      //output << std::setprecision(5)  << eta_range[j]<< "    " << eta_range[j+1] << "    " <<  f2[j]->GetParameter(0)  << "     " << flat[j] << "        " << sqrt( pow(f2[j]->GetParameter(0)*hist_kfsr_dijet->GetBinError(j+1),2)+pow( hist_kfsr_dijet->GetBinContent(j+1)*f2[j]->GetParError(0),2)  ) / flat_norm  << endl;
      //output_loglin << std::setprecision(5)  << eta_range[j]<< "    " << eta_range[j+1] << "    " <<  f1[j]->GetParameter(0)  << "     " << f1[j]->GetParameter(1) << "            " << loglin[j] << "        " << sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()),2)*pow(hist_kfsr_dijet->GetBinError(j+1),2)   + pow(hist_kfsr_dijet->GetBinContent(j+1),2)*(pow(f2[j]->GetParError(0),2)+pow(f1[j]->GetParError(1),2)) ) / loglin_norm  << endl;

      output << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   55   " << ptaverebin_data[j]->FindLastBinAbove(0.)*10 << "   " << 1/flat_norm << " " << hist_kfsr_dijet->GetBinContent(j+1) << "   " << f2[j]->GetParameter(0) << " 0   1 0.0000 0.0" << endl;
      output_loglin << fixed << std::setprecision(6)  << "   " << eta_range[j]<< "  " << eta_range[j+1] << "   11   10 6500   55   " << ptaverebin_data[j]->FindLastBinAbove(0.)*10 << "   " << 1/loglin_norm << " " << hist_kfsr_dijet->GetBinContent(j+1) << "   " << f1[j]->GetParameter(0) << " " << f1[j]->GetParameter(1) <<  "   1 0.0000 0.0" << endl;
      uncerts << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " << ptaverebin_data[j]->FindLastBinAbove(0.)*10 << "       " << sqrt( pow(f2[j]->GetParameter(0)*hist_kfsr_dijet->GetBinError(j+1),2)+pow( hist_kfsr_dijet->GetBinContent(j+1)*f2[j]->GetParError(0),2)  ) / flat_norm  << endl;
      uncerts_loglin << fixed << std::setprecision(6) << " " << eta_range[j] << "  " << eta_range[j+1] << "  3    55    " << ptaverebin_data[j]->FindLastBinAbove(0.)*10 << "       " << sqrt(  pow(f1[j]->GetParameter(0)+f1[j]->GetParameter(1)*TMath::Log(ptave_data[j]->GetMean()),2)*pow(hist_kfsr_dijet->GetBinError(j+1),2)   + pow(hist_kfsr_dijet->GetBinContent(j+1),2)*(pow(f2[j]->GetParError(0),2)+pow(f1[j]->GetParError(1),2)) ) / loglin_norm  << endl;


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
    }






    if(syst==5){
      TFile* outputfile = new TFile(path+"Histo_Res_DiJet_L1"+tag+".root","RECREATE");
    }
    if(syst==0){
      TFile* outputfile = new TFile(path+"Histo_Res_DiJet_L1_"+variation+".root","RECREATE");
    }
    if(syst==1){
      TFile* outputfile = new TFile(path+"Histo_Res_DiJet_L1_"+variation+".root","RECREATE");
    }
    if(syst==2){
      TFile* outputfile = new TFile(path+"Histo_Res_DiJet_L1_"+variation+".root","RECREATE");
    }
    if(syst==3){
      TFile* outputfile = new TFile(path+"Histo_Res_DiJet_L1_"+variation+".root","RECREATE");
    }


    Residual_logpt_DiJet->Write();
    Residual_const_DiJet->Write();
    outputfile->Write();
    outputfile->Close();
    TFile* outputfile2 = new TFile(path+"Histo_ptaveConst_DiJet_L1.root","RECREATE");
    ptave_const_DiJet->Write();
    outputfile2->Write();
    outputfile2->Close();
  }






}
