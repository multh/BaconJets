{


  // change the alpha value here:
  //  TString alpha[1] = {"a010"};
  //  TString alpha[1] = {"a015"};
  //  TString alpha[1] = {"a020"};
  //  TString alpha[1] = {"a025"};
  TString alpha[1] = {"a02"};
  TString alphaName[1] = {"a020"};//TMP due to tricky naming

  // choose the path for the input (uhh2) files and for the output file
  TString path = "/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_noReweight/";

  // input files
  TFile* DATA = new TFile(path+"uhh2.AnalysisModuleRunner.DATA.RunD_AK4CHS.root","READ");
  TFile* QCD = new TFile(path+"uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_pythia8_AK4CHS.root","READ");
  DATA->Print();
  QCD->Print();

  const int n_pt = 9;
  const int n_eta = 17;
  const int n_eta_barrel = 3;



  // required eta binning by mikko
  // TString eta_range[n_eta] = {"0.000", "1.300","1.900","2.500","3.000","3.200","5.000"};
  // double etabins[n_eta] = {0,1.3,1.9,2.5,3.0,3.2,5.0};
  double eta_bins[n_eta] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 5.191};
  TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "5.191"};


  double etabins_barrel[n_eta_barrel] = {0,0.8,1.3};
  // pt ranges from trigger studies
  TString pt_range[n_pt] = {"56.000", "78.000", "100.000", "168.000", "232.000", "300.000", "366.000", "453.000", "562.000"};
  double pt_bins[n_pt] = {56, 78, 100, 168, 232, 300, 366, 453, 562};
  TString eta_range_barrel[n_eta_barrel] = {"0.000", "0.800","1.300"};





  //******************************************************************************************************************************
  //******************************************************************************************************************************
  //******************************************************************************************************************************
  //******************************************************************************************************************************
  //******************************************************************************************************************************



  // get the histos 
  TH1D* data[n_eta-1][n_pt-1];
  TH1D* mc[n_eta-1][n_pt-1];
  // TH1D* data_barrel[n_eta_barrel-1][n_pt-1];
  // TH1D* mc_barrel[n_eta_barrel-1][n_pt-1];
  

  for(int i=0; i<n_eta-1; i++){
    for(int j=0; j<n_pt-1; j++){
      data[i][j] = (TH1D*)DATA->Get(alpha[0]+"/eta_"+eta_range[i]+"_"+eta_range[i+1]+"/pt_"+pt_range[j]+"_"+pt_range[j+1]+"/mpf");
      mc[i][j] = (TH1D*)QCD->Get(alpha[0]+"/eta_"+eta_range[i]+"_"+eta_range[i+1]+"/pt_"+pt_range[j]+"_"+pt_range[j+1]+"/mpf");
      data[i][j]->Print();
      mc[i][j]->Print();
    }
  }

  // for(int i=0; i<n_eta_barrel-1; i++){
  //   for(int j=0; j<n_pt-1; j++){
  //     data_barrel[i][j] = (TH1D*)DATA->Get(alpha[0]+"/eta_"+eta_range_barrel[i]+"_"+eta_range_barrel[i+1]+"/pt_"+pt_range[j]+"_"+pt_range[j+1]+"/mpf");
  //     mc_barrel[i][j] = (TH1D*)QCD->Get(alpha[0]+"/eta_"+eta_range_barrel[i]+"_"+eta_range_barrel[i+1]+"/pt_"+pt_range[j]+"_"+pt_range[j+1]+"/mpf");
  //   }
  // }


  TH1D* dataR[n_eta-1][n_pt-1];
  TH1D* mcR[n_eta-1][n_pt-1];
  //  TH1D* dataR_barrel[n_eta_barrel-1][n_pt-1];
  //  TH1D* mcR_barrel[n_eta_barrel-1][n_pt-1];

  for(int i=0; i<n_eta-1; i++){
    for(int j=0; j<n_pt-1; j++){
      dataR[i][j] = (TH1D*)DATA->Get(alpha[0]+"/eta_"+eta_range[i]+"_"+eta_range[i+1]+"/pt_"+pt_range[j]+"_"+pt_range[j+1]+"/r_rel");
      mcR[i][j] = (TH1D*)QCD->Get(alpha[0]+"/eta_"+eta_range[i]+"_"+eta_range[i+1]+"/pt_"+pt_range[j]+"_"+pt_range[j+1]+"/r_rel");
    }
  }

  // for(int i=0; i<n_eta_barrel-1; i++){
  //   for(int j=0; j<n_pt-1; j++){
  //     dataR_barrel[i][j] = (TH1D*)DATA->Get(alpha[0]+"/eta_"+eta_range_barrel[i]+"_"+eta_range_barrel[i+1]+"/pt_"+pt_range[j]+"_"+pt_range[j+1]+"/r_rel");
  //     mcR_barrel[i][j] = (TH1D*)QCD->Get(alpha[0]+"/eta_"+eta_range_barrel[i]+"_"+eta_range_barrel[i+1]+"/pt_"+pt_range[j]+"_"+pt_range[j+1]+"/r_rel");
  //   }
  // }
 
  // Tgraph to plot result
  double res_xbin_tgraph[n_pt-1];// = {(pt_bins[0]+pt_bins[1])/2, 149, 215.5, 273, 342.5, 423.5, 684};//, 900}; // 86.5, 149, 215.5, 273, 342.5, 423.5, 684};//, 900}; 
  double res_zero[n_pt-1];
  for(int i=0;i<n_pt-1;i++){
    res_xbin_tgraph[i]=(pt_bins[i]+pt_bins[i+1])/2;
    res_zero[i]=(pt_bins[i+1]-pt_bins[i])/2 ;
  }

  //  double res_xbin_tgraph[7] = {66, 107, 191, 240, 306, 379, 468};//, 900};
  //double res_zero[n_pt-1] = {20.5, 42, 24.5, 33, 36.5, 44.5, 216};
  TGraphErrors *mpf_data[n_eta-1];
  TGraphErrors *rrel_data[n_eta-1];
  TGraphErrors *mpf_mc[n_eta-1];
  TGraphErrors *rrel_mc[n_eta-1];
  TGraphErrors *mpf_ratio[n_eta-1];
  TGraphErrors *rrel_ratio[n_eta-1];
  // TGraphErrors *mpf_data_barrel[n_eta-1];
  // TGraphErrors *rrel_data_barrel[n_eta-1];
  // TGraphErrors *mpf_mc_barrel[n_eta-1];
  // TGraphErrors *rrel_mc_barrel[n_eta-1];
  // TGraphErrors *mpf_ratio_barrel[n_eta-1];
  // TGraphErrors *rrel_ratio_barrel[n_eta-1];




  double res_y1[n_pt-1];
  double res_ey_stat1[n_pt-1];
  double res_y2[n_pt-1];
  double res_ey_stat2[n_pt-1];
  double res_y3[n_pt-1];
  double res_ey_stat3[n_pt-1];
  double res_y4[n_pt-1];
  double res_ey_stat4[n_pt-1];
  double res_y5[n_pt-1];
  double res_ey_stat5[n_pt-1];
  double res_y6[n_pt-1];
  double res_ey_stat6[n_pt-1];

  for(int i=0; i<n_eta-1; i++){
    for(int j=0; j<n_pt-1; j++){
      res_y1[j]=data[i][j]->GetMean();
      res_ey_stat1[j] = (data[i][j]->GetRMS() /sqrt(data[i][j]->Integral()));
      res_y2[j] = (dataR[i][j]->GetMean());
      res_ey_stat2[j] = (dataR[i][j]->GetRMS()/sqrt(dataR[i][j]->Integral()));
      res_y3[j]       = (mc[i][j]->GetMean());
      res_ey_stat3[j] = (mc[i][j]->GetRMS()/sqrt(mc[i][j]->Integral()));
      res_y4[j]       = (mcR[i][j]->GetMean());
      res_ey_stat4[j] = (mcR[i][j]->GetRMS()/sqrt(mcR[i][j]->Integral()));
      res_y5[j] = (data[i][j]->GetMean()/mc[i][j]->GetMean());
      res_ey_stat5[j] = (sqrt(pow(data[i][j]->GetRMS(),2)/data[i][j]->Integral()));
      //res_ey_stat5[j] = (sqrt(pow(data[i][j]->GetRMS(),2)/data[i][j]->Integral() + pow(mc[i][j]->GetRMS(),2)/mc[i][j]->Integral() ) );
      res_y6[j] = (dataR[i][j]->GetMean()/mcR[i][j]->GetMean());
      res_ey_stat6[j] = (sqrt(pow(dataR[i][j]->GetRMS(),2)/dataR[i][j]->Integral()));
      //res_ey_stat6[j] = (sqrt(pow(dataR[i][j]->GetRMS(),2)/dataR[i][j]->Integral() + pow(mcR[i][j]->GetRMS(),2)/mcR[i][j]->Integral() ) );
    }
    mpf_data[i] = new TGraphErrors(n_pt-1, res_xbin_tgraph, res_y1, res_zero, res_ey_stat1);
    mpf_data[i]->SetName("mpfchs_dijet_"+alphaName[0]);
    for(int j= mpf_data[i]->GetN()-1; j != -1; --j) {
      if(mpf_data[i]->GetY()[j]==0) mpf_data[i]->RemovePoint(j);
      if(mpf_data[i]->GetEY()[j]!=mpf_data[i]->GetEY()[j] || mpf_data[i]->GetY()[j]!=mpf_data[i]->GetY()[j]  || mpf_data[i]->GetEY()[j]==0) mpf_data[i]->RemovePoint(j);
    }
    for(int j= mpf_data[i]->GetN()-1; j != -1; --j) {
      //cout << mpf_data[i]->GetEY()[j] << endl;
    }
    cout << endl;
    rrel_data[i] = new TGraphErrors(n_pt-1, res_xbin_tgraph, res_y2, res_zero, res_ey_stat2);
    rrel_data[i]->SetName("ptchs_dijet_"+alphaName[0]);
    for(int j= rrel_data[i]->GetN()-1; j != -1; --j) {
      if(rrel_data[i]->GetY()[j]==0) rrel_data[i]->RemovePoint(j);
      if(rrel_data[i]->GetY()[j]!=rrel_data[i]->GetY()[j] || rrel_data[i]->GetY()[j]!=rrel_data[i]->GetY()[j] || rrel_data[i]->GetEY()[j]==0) rrel_data[i]->RemovePoint(j);
    }
    mpf_mc[i] = new TGraphErrors(n_pt-1, res_xbin_tgraph, res_y3, res_zero, res_ey_stat3);
    mpf_mc[i]->SetName("mpfchs_dijet_"+alphaName[0]);
    for(int j= mpf_mc[i]->GetN()-1; j != -1; --j) {
      if(mpf_mc[i]->GetY()[j]==0) mpf_mc[i]->RemovePoint(j);
    }
    rrel_mc[i] = new TGraphErrors(n_pt-1, res_xbin_tgraph, res_y4, res_zero, res_ey_stat4);
    rrel_mc[i]->SetName("ptchs_dijet_"+alphaName[0]);
    for(int j= rrel_mc[i]->GetN()-1; j != -1; --j) {
      if(rrel_mc[i]->GetY()[j]==0) rrel_mc[i]->RemovePoint(j);
    }
    mpf_ratio[i] = new TGraphErrors(n_pt-1, res_xbin_tgraph, res_y5, res_zero, res_ey_stat5);
    mpf_ratio[i]->SetName("mpfchs_dijet_"+alphaName[0]);
    for(int j= mpf_ratio[i]->GetN()-1; j != -1; --j) {
      if(mpf_ratio[i]->GetY()[j]==0  || mpf_ratio[i]->GetY()[j]!=mpf_ratio[i]->GetY()[j] || mpf_ratio[i]->GetEY()[j]!=mpf_ratio[i]->GetEY()[j] || mpf_ratio[i]->GetEY()[j]==0) mpf_ratio[i]->RemovePoint(j);
    }
    for(int j= mpf_ratio[i]->GetN()-1; j != -1; --j) {
      cout << mpf_ratio[i]->GetY()[j] << endl;
    }
    cout << endl;
    rrel_ratio[i] = new TGraphErrors(n_pt-1, res_xbin_tgraph, res_y6, res_zero, res_ey_stat6);
    rrel_ratio[i]->SetName("ptchs_dijet_"+alphaName[0]);
    for(int j= rrel_ratio[i]->GetN()-1; j != -1; --j) {
      if(rrel_ratio[i]->GetY()[j]==0 || rrel_ratio[i]->GetY()[j]!=rrel_ratio[i]->GetY()[j] || rrel_ratio[i]->GetEY()[j]!=rrel_ratio[i]->GetEY()[j] || rrel_ratio[i]->GetEY()[j]==0) rrel_ratio[i]->RemovePoint(j);
    }
  }

  cout<<"Before cleaning"<<endl;
  // remove all data point with less than 100 entries  
  for(int j=0; j<n_eta-1; j++){
    for(int i=n_pt-2; i!=-1; --i){
      if(data[j][i]->GetEntries()<100){
	mpf_data[j]->RemovePoint(i);
	mpf_ratio[j]->RemovePoint(i);
	cout<<"removing points "<<j<<" bin"<<endl;
      }
      if(dataR[j][i]->GetEntries()<20){
	rrel_data[j]->RemovePoint(i);
	rrel_ratio[j]->RemovePoint(i);
	cout<<"removing points "<<j<<" bin"<<endl;
    }
    }
  }
  cout<<"After cleaning"<<endl;


  // double barr_res_y1[n_pt-1];
  // double barr_res_ey_stat1[n_pt-1];
  // double barr_res_y2[n_pt-1];
  // double barr_res_ey_stat2[n_pt-1];
  // double barr_res_y3[n_pt-1];
  // double barr_res_ey_stat3[n_pt-1];
  // double barr_res_y4[n_pt-1];
  // double barr_res_ey_stat4[n_pt-1];
  // double barr_res_y5[n_pt-1];
  // double barr_res_ey_stat5[n_pt-1];
  // double barr_res_y6[n_pt-1];
  // double barr_res_ey_stat6[n_pt-1];

  // for(int i=0; i<2; i++){
  //   for(int j=0; j<n_pt-1; j++){
  //     cout<<"j="<<j<<endl;
  //     data_barrel[i][j]->Print();
  //     barr_res_y1[j]       = (data_barrel[i][j]->GetMean());
  //     cout<<" barr_res_y1[j] = "<< barr_res_y1[j] <<endl;
  //     barr_res_ey_stat1[j] = (1/sqrt(data_barrel[i][j]->Integral()));
  //     barr_res_y2[j] = (dataR_barrel[i][j]->GetMean());
  //     barr_res_ey_stat2[j] = (1/sqrt(dataR_barrel[i][j]->Integral()));
  //     barr_res_y3[j]       = (mc_barrel[i][j]->GetMean());
  //     barr_res_ey_stat3[j] = (1/sqrt(mc_barrel[i][j]->Integral()));
  //     barr_res_y4[j]       = (mcR_barrel[i][j]->GetMean());
  //     barr_res_ey_stat4[j] = (1/sqrt(mcR_barrel[i][j]->Integral()));
  //     barr_res_y5[j] = (data_barrel[i][j]->GetMean()/mc_barrel[i][j]->GetMean());
  //     barr_res_ey_stat5[j] = (sqrt(pow(data_barrel[i][j]->GetRMS(),2)/data_barrel[i][j]->Integral()) );
  //     //barr_res_ey_stat5[j] = (sqrt(pow(data_barrel[i][j]->GetRMS(),2)/data_barrel[i][j]->Integral() + pow(mc_barrel[i][j]->GetRMS(),2)/mc_barrel[i][j]->Integral() ) );
  //     barr_res_y6[j] = (dataR_barrel[i][j]->GetMean()/mcR_barrel[i][j]->GetMean());
  //     barr_res_ey_stat6[j] = (sqrt(pow(dataR_barrel[i][j]->GetRMS(),2)/dataR_barrel[i][j]->Integral() ) );
  //     //barr_res_ey_stat6[j] = (sqrt(pow(dataR_barrel[i][j]->GetRMS(),2)/dataR_barrel[i][j]->Integral() + pow(mcR_barrel[i][j]->GetRMS(),2)/mcR_barrel[i][j]->Integral() ) );
  //   }
  //   cout<<"i="<<i<<endl;
  //   mpf_data_barrel[i] = new TGraphErrors(n_pt-1, res_xbin_tgraph, barr_res_y1, res_zero, barr_res_ey_stat1);
  //   mpf_data_barrel[i]->SetName("mpfchs_dijet_"+alpha[0]);
  //   for(int j= mpf_data_barrel[i]->GetN()-1; j != -1; --j) {
  //     if(mpf_data_barrel[i]->GetY()[j]==0) mpf_data_barrel[i]->RemovePoint(j);
  //   }
  //   rrel_data_barrel[i] = new TGraphErrors(n_pt-1, res_xbin_tgraph, barr_res_y2, res_zero, barr_res_ey_stat2);
  //   rrel_data_barrel[i]->SetName("ptchs_dijet_"+alpha[0]);
  //   for(int j= rrel_data_barrel[i]->GetN()-1; j != -1; --j) {
  //     if(rrel_data_barrel[i]->GetY()[j]==0) rrel_data_barrel[i]->RemovePoint(j);
  //   }
  //   mpf_mc_barrel[i] = new TGraphErrors(n_pt-1, res_xbin_tgraph, barr_res_y3, res_zero, barr_res_ey_stat3);
  //   mpf_mc_barrel[i]->SetName("mpfchs_dijet_"+alpha[0]);
  //   for(int j= mpf_mc_barrel[i]->GetN()-1; j != -1; --j) {
  //     if(mpf_mc_barrel[i]->GetY()[j]==0) mpf_mc_barrel[i]->RemovePoint(j);
  //   }
  //   rrel_mc_barrel[i] = new TGraphErrors(n_pt-1, res_xbin_tgraph, barr_res_y4, res_zero, barr_res_ey_stat4);
  //   rrel_mc_barrel[i]->SetName("ptchs_dijet_"+alpha[0]);
  //   for(int j= rrel_mc_barrel[i]->GetN()-1; j != -1; --j) {
  //     if(rrel_mc_barrel[i]->GetY()[j]==0) rrel_mc_barrel[i]->RemovePoint(j);
  //   }
  //   mpf_ratio_barrel[i] = new TGraphErrors(n_pt-1, res_xbin_tgraph, barr_res_y5, res_zero, barr_res_ey_stat5);
  //   mpf_ratio_barrel[i]->SetName("mpfchs_dijet_"+alpha[0]);
  //   for(int j= mpf_ratio_barrel[i]->GetN()-1; j != -1; --j) {
  //     if(mpf_ratio_barrel[i]->GetY()[j]==0 || mpf_ratio_barrel[i]->GetEY()[j]>10) mpf_ratio_barrel[i]->RemovePoint(j);
  //   }
  //   rrel_ratio_barrel[i] = new TGraphErrors(n_pt-1, res_xbin_tgraph, barr_res_y6, res_zero, barr_res_ey_stat6);
  //   rrel_ratio_barrel[i]->SetName("ptchs_dijet_"+alpha[0]);
  //   for(int j= rrel_ratio_barrel[i]->GetN()-1; j != -1; --j) {
  //     if(rrel_ratio_barrel[i]->GetY()[j]==0 || rrel_ratio_barrel[i]->GetEY()[j]>10) rrel_ratio_barrel[i]->RemovePoint(j);
  //   }
  // }
  


  //mpf_ratio[5]->RemovePoint(6);
  //rrel_ratio[5]->RemovePoint(6);
  mpf_ratio[n_eta-2]->Draw();


  // create output file
  TFile* outputfile = new TFile(path+"JECcombifile_"+alpha[0]+".root","RECREATE");

  cout<<"Util here was fine!"<<endl;
  //  TString eta_output[n_eta-1] = {"eta00-13", "eta13-19","eta19-25","eta25-30","eta30-32","eta32-50"};
  TString eta_output[n_eta-1] = {"eta0000-0261", "eta0216-0522","eta0522-0783","eta0783-1044","eta1044-1305","eta1305-1653","eta1653-1930","eta1930-2172","eta2172-2322","eta2322-2500","eta2500-2650","eta2650-2853","eta2853-2964","eta2964-3139",
				 "eta3139-3489","eta3489-5191"};

  // TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "5.191"};

  //TString eta_output[n_eta-1] = {"eta000-025", "eta025-050","eta050-075","eta075-100","eta100-130","eta130-160","eta160-190","eta190-220","eta220-250","eta250-270","eta270-285","eta285-300","eta300-320","eta320-500"};
  //  TString eta_output_barrel[2] = {"eta00-08", "eta08-13"};

  for(int i=0; i<n_eta-1; i++){
    outputfile->mkdir("ratio/"+eta_output[i]);
    outputfile->cd("ratio/"+eta_output[i]);
    mpf_ratio[i]->Write();
    rrel_ratio[i]->Write();
  }

  for(int i=0; i<n_eta-1; i++){
    outputfile->mkdir("data/"+eta_output[i]);
    outputfile->cd("data/"+eta_output[i]);
    mpf_data[i]->Write();
    rrel_data[i]->Write();
  }

  for(int i=0; i<n_eta-1; i++){
    outputfile->mkdir("mc/"+eta_output[i]);
    outputfile->cd("mc/"+eta_output[i]);
    mpf_mc[i]->Write();
    rrel_mc[i]->Write();
  }

  // for(int i=0; i<2; i++){
  //   outputfile->mkdir("ratio/"+eta_output_barrel[i]);
  //   outputfile->cd("ratio/"+eta_output_barrel[i]);
  //   mpf_ratio_barrel[i]->Write();
  //   rrel_ratio_barrel[i]->Write();
  // }

  // for(int i=0; i<2; i++){
  //   outputfile->mkdir("data/"+eta_output_barrel[i]);
  //   outputfile->cd("data/"+eta_output_barrel[i]);
  //   mpf_data_barrel[i]->Write();
  //   rrel_data_barrel[i]->Write();
  // }

  // for(int i=0; i<2; i++){
  //   outputfile->mkdir("mc/"+eta_output_barrel[i]);
  //   outputfile->cd("mc/"+eta_output_barrel[i]);
  //   mpf_mc_barrel[i]->Write();
  //   rrel_mc_barrel[i]->Write();
  // }



  outputfile->Write();
  outputfile->Close();
  






}
