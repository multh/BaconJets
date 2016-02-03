// with this script kFSR extrapolations over alpha in bins of eta will be calculated

#include "header.h"


TString ToString(int num) {
  ostringstream start;
  start<<num;
  TString start1=start.str();
  return start1;
}

//void kFSR(bool mpfMethod(false), TString path, TFile* datafile, TFile* MCfile){
void kFSR(bool mpfMethod, TString path, TFile* datafile, TFile* MCfile){

  gStyle->SetOptFit(0);
  // datafile->Print();
  // MCfile->Print();

  // get the histos 
  TH1D* data[n_alpha-1][n_eta-1];
  TH1D* mc[n_alpha-1][n_eta-1];
 
 

  for(int i=0; i<n_alpha-1; i++){
    for(int j=0; j<n_eta-1; j++){
      if(mpfMethod){
	//	cout<<alpha_range[i]<<" "<<eta_range[j]<<" "<<eta_range[j+1]<<endl;
	mc[i][j] = (TH1D*)MCfile->Get(alpha_range[i]+"/eta_"+eta_range[j]+"_"+eta_range[j+1]+"/mpf");
	//	mc[i][j]->Print();
	data[i][j] = (TH1D*)datafile->Get(alpha_range[i]+"/eta_"+eta_range[j]+"_"+eta_range[j+1]+"/mpf");
	//	data[i][j]->Print();

      }
      else{
	data[i][j] = (TH1D*)datafile->Get(alpha_range[i]+"/eta_"+eta_range[j]+"_"+eta_range[j+1]+"/r_rel");
	mc[i][j] = (TH1D*)MCfile->Get(alpha_range[i]+"/eta_"+eta_range[j]+"_"+eta_range[j+1]+"/r_rel");
	//	data[i][j]->Print();
	//	mc[i][j]->Print();
      }
    }
  }

 


  // create the ratio histograms
  TH1D* ratio[n_eta-1];
  //TF1* fit_func[n_eta-1];
  for (int j=0; j<n_eta-1; j++){
    TString numstr=ToString(j);
    TString histoname="histo"+numstr;
    ratio[j] = new TH1D(histoname,histoname,n_alpha-1,alpha_bins);
    //    ratio[j]->Print();
  }


  // fill ratio histograms
  for(int j=0; j<n_eta-1; j++){
    for(int i=0; i<n_alpha-1; i++){
      ratio[j]->SetBinContent(i+1, (mc[i][j]->GetMean() / data[i][j]->GetMean()) / (mc[1][j]->GetMean() / data[1][j]->GetMean()) );
      //ratio[j]->SetBinError(i+1, 1/sqrt(data[i][j]->Integral()));
      ratio[j]->SetBinError(i+1, sqrt( pow(data[i][j]->GetRMS(),2)/data[i][j]->Integral() + pow(mc[i][j]->GetRMS(),2)/mc[i][j]->Integral() ));
      //ratio[j]->SetBinError(i+1, sqrt( pow(data[i][j]->GetRMS(),2)/data[i][j]->Integral())); // + pow(mc[i][j]->GetRMS(),2)/mc[i][j]->Integral() ) );
      //ratio[j]->SetBinError(i+1, sqrt( pow(data[i][j]->GetRMS(),2)/data[i][j]->Integral() + pow(mc[i][j]->GetRMS(),2)/mc[i][j]->Integral() ) );
    }
  }


  // create and fill tgrapherrors
  double xbin_tgraph[n_alpha-1] = {0.1,0.2,0.3,0.4};
  double zero[n_alpha-1] = {0, 0, 0, 0};
  TGraphErrors *graph1[n_eta-1];

  double content[n_alpha-1];
  double error[n_alpha-1];
  double content_dijet[n_alpha-1];
  double error_dijet[n_alpha-1];
  for(int j=0; j<n_eta-1; j++){
    for(int i=0; i<n_alpha-1; i++){
      content[i] = (ratio[j]->GetBinContent(i+1));
      error[i] = (ratio[j]->GetBinError(i+1));
    }
    graph1[j] = new TGraphErrors(n_alpha-1, xbin_tgraph, content , zero, error);
  }



  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,0.45,1);

 
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


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // Plots

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  //create plots
  TCanvas* a[n_eta-1];
  TString plotname[n_eta-1];

  for (int j=0; j<n_eta-1; j++){
    if(mpfMethod){
      plotname[j]="mpf_kfsr_eta_"+eta_range[j]+"_"+eta_range[j+1];
    }
    else{
      plotname[j]="dijet_kfsr_eta_"+eta_range[j]+"_"+eta_range[j+1];
    }
    a[j] = new TCanvas(plotname[j], plotname[j], 800,700);
    gStyle->SetOptTitle(0);
    graph1[j]->SetMarkerColor(kBlue);
    graph1[j]->SetMarkerStyle(20);
    graph1[j]->SetLineColor(kBlue);
    graph1[j]->Draw("AP");
    TF1 *pol1 = new TF1("pol1","pol1"); 
    graph1[j]->Fit(pol1,"R");
    //    graph1[j]->Fit("pol1","R");
    graph1[j]->GetXaxis()->SetTitle("#alpha");
    graph1[j]->GetXaxis()->SetTitleSize(0.05);
    //graph1[j]->GetYaxis()->SetTitle("kFSR");
    graph1[j]->GetXaxis()->SetLimits(0.,0.45);
    graph1[j]->GetYaxis()->SetRangeUser(0.92,1.08);
    line->SetLineStyle(2);
    line->Draw("SAME");
    // fill the output.dat file
    if (fp!=NULL) {
      Float_t value = pol1->GetParameter(0);
      Float_t uncert = pol1->GetParError(0);
      fprintf(fp, "%f %f\n",value,uncert);
    }
    plotkfsr->SetBinContent(j+1,pol1->GetParameter(0));
    plotkfsr->SetBinError(j+1,pol1->GetParError(0));
    if(mpfMethod){
      kFSR_MPF->SetBinContent(j+1,pol1->GetParameter(0));
      kFSR_MPF->SetBinError(j+1,pol1->GetParError(0));
    }
    else{
      kFSR_DiJet->SetBinContent(j+1,pol1->GetParameter(0));
      kFSR_DiJet->SetBinError(j+1,pol1->GetParError(0));
    }
    TLegend *leg1;
    leg1 = new TLegend(0.13,0.68,0.35,0.88,"","brNDC");//x+0.1
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.045);
    leg1->SetFillColor(10);
    leg1->SetLineColor(1);
    leg1->SetTextFont(42);
    if(mpfMethod){
      leg1->SetHeader("MPF, "+eta_range3[j]+"#leq|#eta|<"+eta_range3[j+1]);
    }
    else{
      leg1->SetHeader("p_{T} balance, "+eta_range3[j]+"#leq|#eta|<"+eta_range3[j+1]);
    }
    leg1->AddEntry(graph1[j], "R(MC)/R(DATA)","P");
    leg1->AddEntry(pol1, "linear fit","L");
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
