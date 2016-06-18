//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// This macro reads L2 residuals for LOGLIN pt extrapolation 
// from txt files
// 2D plots are created for visual check
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "header.h"
#include "tdrstyle_mod15.C"

void L2ResTxtTest(TString path, TString txttag, TString jettag, TString tag, double al_cut=0.2){
  //lumi_13TeV = "2.11 fb^{-1}";
  lumi_13TeV = "0.8 fb^{-1}";
      //lumi_13TeV = "589.3 pb^{-1}";
  gStyle->SetPalette(55);  
gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);  
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.21);
  // For the axis titles:
  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.06, "XYZ");
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.1);
  gStyle->SetTitleOffset(1.45,"Z");
  // For the axis labels:
  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");
  // For the axis:
  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);



  TString JetDescrib;
  if (jettag=="AK4PFchs") JetDescrib = "Anti-k_{t} R = 0.4, PF+CHS";
  if (jettag=="AK4PFpuppi") JetDescrib = "Anti-k_{t} R = 0.4, PF+PUPPI"; 
  if (jettag=="AK8PFchs") JetDescrib = "Anti-k_{t} R = 0.8, PF+CHS";                                                                             
  if (jettag=="AK8PFpuppi") JetDescrib = "Anti-k_{t} R = 0.8, PF+PUPPI"; 

  TString campain_prefix = "Spring16_25ns_";
  TString method[2] = {"MPF", "pT"};
  ifstream input_file_val[2];
  ifstream input_file_err[2];
  TH2D *hcorr_eta_pT[2];   TH2D *hcorr_eta_pT_extr[2];
  TH2D *herrcorr_eta_pT[2];   TH2D *herrcorr_eta_pT_extr[2];
  double pT_min = 1000;     double pT_max = 0.;
  for(int i=0;i<2;i++){
  //  for(int i=0;i<1;i++){
    TString histname = "hcorr_eta_pT_"+method[i];
    TString histerrname = "herrcorr_eta_pT_"+method[i];
    TString name = path+"output/"+campain_prefix+method[i]+"_LOGLIN_L2Residual_"+txttag+"_"+jettag+".txt";
    TString name_err = path+"output/"+campain_prefix+method[i]+"_LOGLIN_L2Residual_"+txttag+"_"+jettag+".txt.STAT";
    //    cout<<name.Data()<<endl;
    input_file_val[i].open(name.Data(),ios::in);
    input_file_err[i].open(name_err.Data(),ios::in);
    string line;
    getline (input_file_val[i],line);    
    //    cout<<line<<endl;  
    double eta1, eta2, x,y,z, pT1, pT2, norm, kFSR, f0, f1, kFSR_err, cov00, cov11, cov01;
    Int_t nlines = 0;
    vector<double> eta; vector<double> pT1_val; vector<double> pT2_val;
    vector<double> norm_val; 
    vector<double> kFSR_val;     
    vector<double> f0_val;    vector<double> f1_val;
    while (1) {
      input_file_val[i] >> eta1 >> eta2 >> x >> y >> z >> pT1 >> pT2 >> norm >> kFSR >> f0 >> f1 >> x >> y >> z;
      eta.push_back(eta1);
      pT1_val.push_back(pT1); pT2_val.push_back(pT2);
      norm_val.push_back(norm);
      kFSR_val.push_back(kFSR);
      f0_val.push_back(f0); f1_val.push_back(f1);
      if(pT_min>pT1) pT_min = pT1;
      if(pT_max<pT2) pT_max = pT2;
      if (!input_file_val[i].good()) break;
      if (nlines < 3) 
      //      if (nlines < 22) 
	printf("eta1=%8f, pT2=%8f, f0=%8f, f1=%8f\n",eta1,pT2,f0,f1);
      nlines++;
      //      cout<<"last eta "<<eta1<<" "<<eta2<<endl;
    }
    //    cout<<"last eta "<<eta1<<" "<<eta2<<endl;
    //    cout<<eta.size()<<" "<<eta[eta.size()-1]<<endl;
    eta[eta.size()-1] = eta2;
    double eta_vals[eta.size()];
    for(int j=0;j<eta.size();j++){
      //      cout<<" "<<eta[j];
      eta_vals[j] = eta[j];

    }
    printf(" found %d points\n",nlines);
    input_file_val[i].close();
    //    int nPt = 100;
    int nPt = 40;
    hcorr_eta_pT[i] = new TH2D(histname,method[i]+"; #eta; #bar{p}_{T}, GeV",eta.size()-1,eta_vals,nPt,pT_min,pT_max+200.);
    hcorr_eta_pT_extr[i] = new TH2D(histname+"_extr",method[i]+"; #eta; #bar{p}_{T}, GeV",eta.size()-1,eta_vals,nPt,pT_min,pT_max+200.);
    herrcorr_eta_pT[i] = new TH2D(histerrname,method[i]+"; #eta; #bar{p}_{T}, GeV",eta.size()-1,eta_vals,nPt,pT_min,pT_max+200.);
    herrcorr_eta_pT_extr[i] = new TH2D(histerrname+"_extr",method[i]+"; #eta; #bar{p}_{T}, GeV",eta.size()-1,eta_vals,nPt,pT_min,pT_max+200.);
    double DpT = (pT_max+200. - (pT_min))/nPt;
    //  cout<<"Dpt = "<<DpT<<endl;
    vector<TLine*> pT_max_levels;
    for(int xi=0;xi<eta.size()-1;xi++){
      //      cout<<eta[xi+1]<<endl;
      TLine *pT_max_levels_curr = new TLine(eta[xi],pT2_val[xi],eta[xi+1],pT2_val[xi]);
      pT_max_levels_curr->SetLineColor(kBlack);
      pT_max_levels_curr->SetLineWidth(2);
      pT_max_levels_curr->SetLineStyle(2);
      pT_max_levels.push_back(pT_max_levels_curr);
      for(int yi=0;yi<nPt;yi++){
      //  for(int yi=70;yi<nPt;yi++){
	double pt_curr= pT_min + yi*DpT;
	//	cout<<pt_curr<<" "<< pT1_val[xi] <<" "<<pT2_val[xi] <<endl;
	double corr = norm_val[xi]*kFSR_val[xi]*(f0_val[xi]+log(pt_curr)*f1_val[xi]);
	hcorr_eta_pT_extr[i]->Fill(0.5*(eta[xi+1]+eta[xi]),pt_curr,corr);
	if(pt_curr>pT1_val[xi] && pt_curr<pT2_val[xi]){
	  //	  if(corr<0.7)	  cout<<0.5*(eta[xi+1]+eta[xi])<<" "<<f0_val[xi]<<" "<<f1_val[xi]<<" "<<kFSR_val[xi]<<" "<<pt_curr<<" "<<corr<<endl;
	  //	  cout<<0.5*(eta[xi+1]+eta[xi])<<" "<<pt_curr<<" "<<corr<<endl;
	  hcorr_eta_pT[i]->Fill(0.5*(eta[xi+1]+eta[xi]),pt_curr,corr);
	}
      }
    }

    // //   // Draw results                             
    //    bool kSquare = true;            
    hcorr_eta_pT_extr[i]->GetZaxis()->SetTitle("L2Residuals");
    hcorr_eta_pT_extr[i]->GetZaxis()->SetTitleSize(0.05);
    hcorr_eta_pT_extr[i]->GetZaxis()->SetLabelSize(0.05);
    //    hcorr_eta_pT_extr[i]->GetZaxis()->SetRangeUser(0.0,3.0);
    hcorr_eta_pT_extr[i]->GetZaxis()->SetRangeUser(0.55,1.25);
    TCanvas *c = new TCanvas("L2Res_LOGLIN_"+method[i],"",800,600);
    hcorr_eta_pT_extr[i]->Draw("COLZ");
    c->SetLogy(1);
    for(int ieta=0;ieta<eta.size()-1;ieta++)
      pT_max_levels[ieta]->Draw("SAME");
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045);
    //    tex->DrawLatex(0.8,0.94,lumi_13TeV);
    tex->DrawLatex(0.05,0.94,JetDescrib);
    CMS_lumi(c,4,10);
    c->SaveAs(path+"plots/L2Res_eta_pT_"+method[i]+"_"+jettag+"_"+txttag+".pdf");  
    // //    pT_max_levels[0]->Draw();
 
    vector<double> kFSR_err_val, cov00_val, cov11_val, cov01_val;
    nlines = 0;
    getline (input_file_err[i],line);
    //    cout<<line<<endl;  
    while (1) {  
      input_file_err[i] >> eta1 >> eta2 >> x >> pT1 >> pT2 >> kFSR_err >> cov00 >> cov11 >> cov01;        
      if (!input_file_err[i].good()) break;
      if (nlines < 3) printf("eta1=%8f, pT2=%8f, kFSR_err=%8f, cov00=%8f\n",eta1,pT2,kFSR_err,cov00);
      nlines++;
      kFSR_err_val.push_back(kFSR_err);
      cov00_val.push_back(cov00); cov11_val.push_back(cov11); cov01_val.push_back(cov01);
    }     
    printf(" found %d points\n",nlines);

    for(int xi=0;xi<eta.size()-1;xi++){
      for(int yi=0;yi<nPt;yi++){
	double pt_curr= pT_min + yi*DpT;
	// sqrt(pow((f0_val[xi]+log(pt_curr)*f1_val[xi]),2)*pow(kFSR_err_val[xi],2)   
	//      +pow(kFSR_val[xi],2)*(cov00_val[xi]+
	// 			   cov11_val[xi]*pow(TMath::Log(pt_curr),2)+2*cov22_val[xi]*TMath::Log(pt_curr)));


	double corr_err = sqrt(pow((f0_val[xi]+TMath::Log(pt_curr)*f1_val[xi]),2)*pow(kFSR_err_val[xi],2)                           
					     +pow(kFSR_val[xi],2)*fabs(cov00_val[xi]+
	 							   cov11_val[xi]*pow(TMath::Log(pt_curr),2)
								   +2*cov01_val[xi]*TMath::Log(pt_curr)));
	herrcorr_eta_pT_extr[i]->Fill(0.5*(eta[xi+1]+eta[xi]),pt_curr,corr_err);
	if(pt_curr>pT1_val[xi] && pt_curr<pT2_val[xi]){
	  // cout<<0.5*(eta[xi+1]+eta[xi])<<" "<<pt_curr<<" "<<corr_err<<" "
	  //     <<pow((f0_val[xi]+log(pt_curr)*f1_val[xi]),2)*pow(kFSR_err_val[xi],2)<<" "
	  //     <<pow(kFSR_val[xi],2)*fabs(cov00_val[xi]+                                                                 
	  // 	 cov11_val[xi]*pow(TMath::Log(pt_curr),2)                                       
	  // 	 +2*cov01_val[xi]*TMath::Log(pt_curr))
	  //     <<endl;
	  herrcorr_eta_pT[i]->Fill(0.5*(eta[xi+1]+eta[xi]),pt_curr,corr_err);
	}
      }
    }
    herrcorr_eta_pT_extr[i]->GetZaxis()->SetTitle("Stat. uncertainty");
    herrcorr_eta_pT_extr[i]->GetZaxis()->SetTitleSize(0.05);
    herrcorr_eta_pT_extr[i]->GetZaxis()->SetLabelSize(0.05);
    //    herrcorr_eta_pT_extr[i]->GetZaxis()->SetRangeUser(0.,0.05);
    herrcorr_eta_pT_extr[i]->GetZaxis()->SetRangeUser(0.,0.1);
    TCanvas *c2 = new TCanvas("L2Res_STAT_LOGLIN_"+method[i],"",800,600);
    herrcorr_eta_pT_extr[i]->Draw("COLZ");
    c2->SetLogy(1);
    for(int ieta=0;ieta<eta.size()-1;ieta++)
      pT_max_levels[ieta]->Draw("SAME");
    //    TLatex *tex = new TLatex();
    //    tex->SetNDC();
    //    tex->SetTextSize(0.045);
    //    tex->DrawLatex(0.8,0.94,lumi_13TeV);
    tex->DrawLatex(0.05,0.94,JetDescrib);
    CMS_lumi(c2,4,10);
    c2->SaveAs(path+"plots/L2Res_eta_pT_STAT_"+method[i]+"_"+jettag+"_"+txttag+".pdf");  
    // //    pT_max_levels[0]->Draw();
 

    // input_file_err[i].close();
    // herrcorr_eta_pT[i] = new TH2D(histerrname,method[i]+"; #eta; #bar{p}_{T}, GeV",eta.size(),eta_vals,20,5,2005.);
    // 
  }

}

