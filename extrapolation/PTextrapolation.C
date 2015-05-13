
// with this script kFSR extrapolations over alpha in bins of eta will be calculated

TString ToString(int num) {
  ostringstream start;
  start<<num;
  TString start1=start.str();
  return start1;
}

void PTextrapolation(){

  // input DATA file
  TFile* data2012 = new TFile("/nfs/dust/cms/user/mstoev/sFrame_new/JEC/output/uhh2.AnalysisModuleRunner.DATA.DATA_test1234.root","READ");

  // input MC file
  TFile* mc8tev = new TFile("/nfs/dust/cms/user/mstoev/sFrame_new/JEC/output/uhh2.AnalysisModuleRunner.MC.QCDPt15to3000_8tev_test1234.root","READ");

  // define number of bins
  const int n_alpha = 5;
  const int n_eta = 30;
  const int n_pt = 8;

  // define the bin ranges
  TString alpha_range[n_alpha] = {"0.000", "0.100", "0.200", "0.300", "0.400"};
  TString eta_range[n_eta] = {"0.000", "0.087", "0.174", "0.261", "0.348", "0.435", "0.522", "0.609", "0.696", "0.783", "0.879", "0.957", "1.044", "1.131", "1.218", "1.305", "1.392", "1.479", "1.566", "1.653", "1.830", "1.930", "2.043", "2.172", "2.322", "2.500", "2.853", "2.964", "3.139", "5.232"};
  TString pt_range[n_pt] = {"66.000", "107.000", "191.000", "240.000", "306.000", "379.000", "468.000", "900.000"};
  double alpha_bins[n_alpha] = {0.000, 0.100, 0.200, 0.300, 0.400};
  double eta_bins[n_eta] = {0.000, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.830, 1.930, 2.043, 2.172, 2.322, 2.500, 2.853, 2.964, 3.139, 5.232};
  double pt_bins[n_pt] = {66, 107, 191, 240, 306, 379, 468, 900};


  TString eta_range2[n_eta] = {"0000", "0087", "0174", "0261", "0348", "0435", "0522", "0609", "0696", "0783", "0879", "0957", "1044", "1131", "1218", "1305", "1392", "1479", "1566", "1653", "1830", "1930", "2043", "2172", "2322", "2500", "2853", "2964", "3139", "5232"};


  // get the histos 
  TH1D* data[n_eta-1][n_pt-1];
  TH1D* mc[n_eta-1][n_pt-1];
  
  for(int i=0; i<n_eta-1; i++){
    for(int j=0; j<n_pt-1; j++){
      data[i][j] = (TH1D*)data2012->Get("/eta_"+eta_range[i]+"_"+eta_range[i+1]+"/pt_"+pt_range[j]+"_"+pt_range[j+1]+"/mpf");
      mc[i][j] = (TH1D*)mc8tev->Get("/eta_"+eta_range[i]+"_"+eta_range[i+1]+"/pt_"+pt_range[j]+"_"+pt_range[j+1]+"/mpf");
    }
  }

  
  double xbin_tgraph[n_pt-1] = {107,191,240,306,379,468,900};
  double zero[n_pt-1] = {0, 0, 0, 0, 0, 0, 0};
  TGraphErrors *graph1[n_eta-1];
  

  

  // create the ratio histogram
  TH1D* ratio[n_eta-1];
  //TF1* fit_func[n_eta-1];
  for (int j=0; j<n_eta-1; j++){
    TString numstr=ToString(j);
    TString histoname="histo"+numstr;
    ratio[j] = new TH1D(histoname,histoname,n_pt-1,pt_bins);
    //fit_func[j] = new TF1("lin_fit[j]", "pol1", 0., n_alpha);
  }

  /*
  double kfsr[n_eta-1];
  ifstream inp;
  inp.open("extrapolation.dat");
  for (int i=0; i<n_eta-1; i++){
    inp >> kfsr[i];
  }
  inp.close();
  */

  
  // fill ratio histograms
  for(int j=0; j<n_eta-1; j++){
    for(int i=0; i<n_pt-1; i++){
      ratio[j]->SetBinContent(i+1, (mc[i][j]->GetMean() / data[i][j]->GetMean()));
      ratio[j]->SetBinError(i+1, 1/data[i][j]->Integral());
    }
  }
  


  for(int j=0; j<n_eta-1; j++){
    double content[n_pt-1] = {ratio[j]->GetBinContent(1),ratio[j]->GetBinContent(2),ratio[j]->GetBinContent(3),ratio[j]->GetBinContent(4),ratio[j]->GetBinContent(5),ratio[j]->GetBinContent(6),ratio[j]->GetBinContent(7)};
    double error[n_pt-1] = {ratio[j]->GetBinError(1),ratio[j]->GetBinError(2),ratio[j]->GetBinError(3),ratio[j]->GetBinError(4),ratio[j]->GetBinError(5),ratio[j]->GetBinError(6),ratio[j]->GetBinError(7)};
    graph1[j] = new TGraphErrors(n_pt-1, xbin_tgraph, content , zero, error);
  }
  


  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,1000,1);

  // create a function for the loglinear fit
  TF1 * f1[n_eta-1];
  // create a function for the constant fit
  TF1 * f2[n_eta-1];

  // create output .dat file, including the kFSR extrapolation (alpha->0)
  FILE *fp = fopen("pTextrapolations.dat","w");
  FILE *fp2 = fopen("pTconstantExtrapolation.dat","w");
  TCanvas* a[n_eta-1];
  for (int j=0; j<n_eta-1; j++){
    TString numstr=ToString(j);
    TString plotname="eta_"+eta_range[j]+"_"+eta_range[j+1];
    a[j] = new TCanvas(plotname, plotname, 800,600);
    gStyle->SetOptTitle(0);
    gPad->SetLogx();
    graph1[j]->SetMarkerColor(kBlue);
    graph1[j]->SetMarkerStyle(20);
    graph1[j]->SetLineColor(kBlue);
    graph1[j]->Draw("AP");
    f1[j] = new TF1(plotname+"f1","[0]*TMath::Power(x,[1])");
    f2[j] = new TF1(plotname+"f2","[0]");
    f2[j]->SetLineColor(kBlue);
    f2[j]->SetLineStyle(3);
    graph1[j]->Fit(plotname+"f1","W");
    graph1[j]->Fit(plotname+"f2","+");
    graph1[j]->GetXaxis()->SetTitle("#bar{p}_{T}");
    //graph1[j]->GetYaxis()->SetTitle("kFSR");
    graph1[j]->GetXaxis()->SetLimits(30,1000);
    graph1[j]->GetYaxis()->SetRangeUser(0.95,1.05);
    line->SetLineStyle(2);
    line->Draw("SAME");
    // fill the output.dat file
    if (fp!=NULL) {
      // getting the slope parameter from the loglinear fit
      Float_t value = f1[j]->GetParameter(1);
      fprintf(fp, "%f\n",value);
    }
    if (fp2!=NULL) {
      // getting the slope parameter from the loglinear fit
      Float_t value = f2[j]->GetParameter(0);
      fprintf(fp2, "%f\n",value);
    }
    TLegend *leg1;
    leg1 = new TLegend(0.15,0.65,0.60,0.85,"","brNDC");//x+0.1
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.035);
    leg1->SetFillColor(10);
    leg1->SetLineColor(1);
    leg1->SetTextFont(42);
    leg1->SetHeader("eta_"+eta_range[j]+"_"+eta_range[j+1]);
    leg1->AddEntry(graph1[j], "R(MC)/R(DATA)","P");
    leg1->AddEntry(f1[j], "loglinear fit","L");
    leg1->AddEntry(f2[j], "constant fit","L");
    leg1->Draw();

    //a[j]->Print("kFSR_eta_"+eta_range[j]+"_"+eta_range[j+1]+".pdf");
    a[j]->Print("pTextrapol_pT_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");

  }
  fclose(fp);
  fclose(fp2);


}
