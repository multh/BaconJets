
// with this script kFSR extrapolations over alpha in bins of eta will be calculated

TString ToString(int num) {
  ostringstream start;
  start<<num;
  TString start1=start.str();
  return start1;
}

void kFSR(){

  // input DATA file
  TFile* data2012 = new TFile("/nfs/dust/cms/user/mstoev/sFrame_new/JEC/output/uhh2.AnalysisModuleRunner.DATA.DATA_test1234.root","READ"); //uhh2.AnalysisModuleRunner.DATA.QCDPt15to3000_2012data_test.root","READ");

  // input MC file
  TFile* mc8tev = new TFile("/nfs/dust/cms/user/mstoev/sFrame_new/JEC/output/uhh2.AnalysisModuleRunner.MC.QCDPt15to3000_8tev_test1234.root","READ");

  // define number of bins
  const int n_alpha = 5;
  const int n_eta = 30;

  // define the bin ranges
  TString alpha_range[n_alpha] = {"0.000", "0.100", "0.200", "0.300", "0.400"};
  TString eta_range[n_eta] = {"0.000", "0.087", "0.174", "0.261", "0.348", "0.435", "0.522", "0.609", "0.696", "0.783", "0.879", "0.957", "1.044", "1.131", "1.218", "1.305", "1.392", "1.479", "1.566", "1.653", "1.830", "1.930", "2.043", "2.172", "2.322", "2.500", "2.853", "2.964", "3.139", "5.232"};
  double alpha_bins[n_alpha] = {0.000, 0.100, 0.200, 0.300, 0.400};
  double eta_bins[n_eta] = {0.000, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.830, 1.930, 2.043, 2.172, 2.322, 2.500, 2.853, 2.964, 3.139, 5.232};


  TString eta_range2[n_eta] = {"0000", "0087", "0174", "0261", "0348", "0435", "0522", "0609", "0696", "0783", "0879", "0957", "1044", "1131", "1218", "1305", "1392", "1479", "1566", "1653", "1830", "1930", "2043", "2172", "2322", "2500", "2853", "2964", "3139", "5232"};


  // get the histos 
  TH1D* data[n_alpha-1][n_eta-1];
  TH1D* mc[n_alpha-1][n_eta-1];
  // for the dijet balance
  TH1D* data_dijet[n_alpha-1][n_eta-1];
  TH1D* mc_dijet[n_alpha-1][n_eta-1];
  
  for(int i=0; i<n_alpha-1; i++){
    for(int j=0; j<n_eta-1; j++){
      data[i][j] = (TH1D*)data2012->Get("alpha_"+alpha_range[i]+"_"+alpha_range[i+1]+"/eta_"+eta_range[j]+"_"+eta_range[j+1]+"/mpf");
      mc[i][j] = (TH1D*)mc8tev->Get("alpha_"+alpha_range[i]+"_"+alpha_range[i+1]+"/eta_"+eta_range[j]+"_"+eta_range[j+1]+"/mpf");
      data_dijet[i][j] = (TH1D*)data2012->Get("alpha_"+alpha_range[i]+"_"+alpha_range[i+1]+"/eta_"+eta_range[j]+"_"+eta_range[j+1]+"/asym");
      mc_dijet[i][j] = (TH1D*)mc8tev->Get("alpha_"+alpha_range[i]+"_"+alpha_range[i+1]+"/eta_"+eta_range[j]+"_"+eta_range[j+1]+"/asym");
    }
  }


  double xbin_tgraph[n_alpha-1] = {0.1,0.2,0.3,0.4};
  double zero[n_alpha-1] = {0, 0, 0, 0};
  TGraphErrors *graph1[n_eta-1];
  TGraphErrors *graph_dijet[n_eta-1];


  // create the ratio histograms
  TH1D* ratio[n_eta-1];
  TH1D* ratio_dijet[n_eta-1];
  //TF1* fit_func[n_eta-1];
  for (int j=0; j<n_eta-1; j++){
    TString numstr=ToString(j);
    TString histoname="histo"+numstr;
    ratio[j] = new TH1D(histoname,histoname,n_alpha-1,alpha_bins);
    ratio_dijet[j] = new TH1D("dijet_"+histoname,"dijet_"+histoname,n_alpha-1,alpha_bins);
  }


  // fill ratio histograms
  for(int j=0; j<n_eta-1; j++){
    for(int i=0; i<n_alpha-1; i++){
      ratio[j]->SetBinContent(i+1, (mc[i][j]->GetMean() / data[i][j]->GetMean()) / (mc[1][j]->GetMean() / data[1][j]->GetMean()) );
      ratio[j]->SetBinError(i+1, 1/data[i][j]->Integral());
      ratio_dijet[j]->SetBinContent(i+1, (  ( ((1+mc_dijet[i][j]->GetMean())/(1-mc_dijet[i][j]->GetMean())) / ((1+data_dijet[i][j]->GetMean())/(1-data_dijet[i][j]->GetMean())) ) / ( ((1+mc_dijet[1][j]->GetMean())/(1-mc_dijet[1][j]->GetMean())) / ((1+data_dijet[1][j]->GetMean())/(1-data_dijet[1][j]->GetMean())) )  ) );
      ratio_dijet[j]->SetBinError(i+1, 1/data_dijet[i][j]->Integral());
    }
  }


  for(int j=0; j<n_eta-1; j++){
    double content[n_alpha-1] = {ratio[j]->GetBinContent(1),ratio[j]->GetBinContent(2),ratio[j]->GetBinContent(3),ratio[j]->GetBinContent(4)};
    double error[n_alpha-1] = {ratio[j]->GetBinError(1),ratio[j]->GetBinError(2),ratio[j]->GetBinError(3),ratio[j]->GetBinError(4)};
    graph1[j] = new TGraphErrors(n_alpha-1, xbin_tgraph, content , zero, error);
    double content_dijet[n_alpha-1] = {ratio_dijet[j]->GetBinContent(1),ratio_dijet[j]->GetBinContent(2),ratio_dijet[j]->GetBinContent(3),ratio_dijet[j]->GetBinContent(4)};
    double error_dijet[n_alpha-1] = {ratio_dijet[j]->GetBinError(1),ratio_dijet[j]->GetBinError(2),ratio_dijet[j]->GetBinError(3),ratio_dijet[j]->GetBinError(4)};
    graph_dijet[j] = new TGraphErrors(n_alpha-1, xbin_tgraph, content_dijet , zero, error_dijet);
  }



  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,0.45,1);

  // choose mpf or dijet balance
  bool mpfMethod(true);


  if(mpfMethod){
  // create output .dat file, including the kFSR extrapolation (alpha->0)
  FILE *fp = fopen("extrapolation.dat","w");
  TCanvas* a[n_eta-1];
  for (int j=0; j<n_eta-1; j++){
    TString numstr=ToString(j);
    TString plotname="eta_"+eta_range[j]+"_"+eta_range[j+1];
    a[j] = new TCanvas(plotname, plotname, 800,600);
    gStyle->SetOptTitle(0);
    graph1[j]->SetMarkerColor(kBlue);
    graph1[j]->SetMarkerStyle(20);
    graph1[j]->SetLineColor(kBlue);
    graph1[j]->Draw("AP");
    graph1[j]->Fit("pol1","R");
    graph1[j]->GetXaxis()->SetTitle("#alpha");
    //graph1[j]->GetYaxis()->SetTitle("kFSR");
    graph1[j]->GetXaxis()->SetLimits(0.,0.45);
    graph1[j]->GetYaxis()->SetRangeUser(0.997,1.003);
    line->SetLineStyle(2);
    line->Draw("SAME");
    // fill the output.dat file
    if (fp!=NULL) {
      Float_t value = pol1->GetParameter(0);
      fprintf(fp, "%f\n",value);
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
    leg1->AddEntry(pol1, "linear fit","L");
    leg1->Draw();

    a[j]->Print("kFSR_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");

  }
  fclose(fp);
  }

  // if dijet balance..
  else{
  // create output .dat file, including the kFSR extrapolation (alpha->0)
  FILE *fp = fopen("extrapolation.dat","w");
  TCanvas* a[n_eta-1];
  for (int j=0; j<n_eta-1; j++){
    TString numstr=ToString(j);
    TString plotname="eta_"+eta_range[j]+"_"+eta_range[j+1];
    a[j] = new TCanvas(plotname, plotname, 800,600);
    gStyle->SetOptTitle(0);
    graph_dijet[j]->SetMarkerColor(kBlue);
    graph_dijet[j]->SetMarkerStyle(20);
    graph_dijet[j]->SetLineColor(kBlue);
    graph_dijet[j]->Draw("AP");
    graph_dijet[j]->Fit("pol1","R");
    graph_dijet[j]->GetXaxis()->SetTitle("#alpha");
    //graph1[j]->GetYaxis()->SetTitle("kFSR");
    graph_dijet[j]->GetXaxis()->SetLimits(0.,0.45);
    graph_dijet[j]->GetYaxis()->SetRangeUser(0.997,1.003);
    line->SetLineStyle(2);
    line->Draw("SAME");
    // fill the output.dat file
    if (fp!=NULL) {
      Float_t value = pol1->GetParameter(0);
      fprintf(fp, "%f\n",value);
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
    leg1->AddEntry(pol1, "linear fit","L");
    leg1->Draw();

    a[j]->Print("kFSR_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");

  }
  fclose(fp);
  }

}
