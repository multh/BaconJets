#include "header.h"
//#include "tdrstyle_mod14.C"
#include "tdrstyle_mod15.C"

void L2ResOutput(TString path, TString txttag, TString jettag, TString tag, double al_cut=0.2){

  TFile* f_Res_mpf = new TFile(path+"Histo_Res_MPF_L1_"+txttag+"_"+jettag+tag+".root","READ");   
  TFile* f_Res_dijet = new TFile(path+"Histo_Res_DiJet_L1_"+txttag+"_"+jettag+tag+".root","READ");  

  //    TString JetDescrib = "Anti-k_{t} R = 0.4, PF+CHS";
  //  TString JetDescrib = "Anti-k_{t} R = 0.4, PF+PUPPI";
  //  TString JetDescrib = "Anti-k_{t} R = 0.8, PF+CHS";
  //  TString JetDescrib = "Anti-k_{t} R = 0.8, PF+PUPPI";
  TString JetDescrib;                                                                                                                            
  if (jettag=="AK4PFchs") JetDescrib = "Anti-k_{t} R = 0.4, PF+CHS";                                                                             
  if (jettag=="AK4PFpuppi") JetDescrib = "Anti-k_{t} R = 0.4, PF+PUPPI";                                                                         
  if (jettag=="AK8PFchs") JetDescrib = "Anti-k_{t} R = 0.8, PF+CHS";                                                                             
  if (jettag=="AK8PFpuppi") JetDescrib = "Anti-k_{t} R = 0.8, PF+PUPPI"; 

  //plot results for "nominal" variation
 

  // get the (R_{MC}/R_{DATA}) hists for MPF and pt balance
  TH1D* pt_depend_const_mpf = (TH1D*)f_Res_mpf->Get("ptave_const_mpf");
  TH1D* pt_depend_logpt_mpf = (TH1D*)f_Res_mpf->Get("ptave_logpt_mpf");

  TH1D* pt_depend_const_dijet = (TH1D*)f_Res_dijet->Get("ptave_const_dijet");
  TH1D* pt_depend_logpt_dijet = (TH1D*)f_Res_dijet->Get("ptave_logpt_dijet");

  // get the kFSR hists for MPF and pt balance
  TH1D* kfsr_mpf = (TH1D*)f_Res_mpf->Get("kfsr_mpf");
  TH1D* kfsr_mpf_fit = (TH1D*)f_Res_mpf->Get("hist_kfsr_fit_mpf");
  TH1D* kfsr_dijet = (TH1D*)f_Res_dijet->Get("kfsr_dijet");
  TH1D* kfsr_dijet_fit = (TH1D*)f_Res_dijet->Get("hist_kfsr_fit_dijet");

  //get L2Res hists for MPF and pt balance
  TH1D* res_const_mpf_kfsrfit = (TH1D*)f_Res_mpf->Get("res_const_mpf");
  TH1D* res_const_dijet_kfsrfit = (TH1D*)f_Res_dijet->Get("res_const_dijet");
  TH1D* res_logpt_mpf_kfsrfit = (TH1D*)f_Res_mpf->Get("res_logpt_mpf");
  TH1D* res_logpt_dijet_kfsrfit = (TH1D*)f_Res_dijet->Get("res_logpt_dijet");

  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,5.191,1);
  res_const_mpf_kfsrfit->SetLineWidth(2);
  res_const_mpf_kfsrfit->SetLineColor(kRed+1);
  res_const_dijet_kfsrfit->SetLineWidth(2);
  res_const_dijet_kfsrfit->SetLineColor(kBlue+1);
  res_const_mpf_kfsrfit->SetLineStyle(2);
  res_const_dijet_kfsrfit->SetLineStyle(2);
  res_logpt_mpf_kfsrfit->SetLineWidth(2);
  res_logpt_mpf_kfsrfit->SetLineColor(kRed+1);
  res_logpt_dijet_kfsrfit->SetLineWidth(2);
  res_logpt_dijet_kfsrfit->SetLineColor(kBlue+1);

  pt_depend_const_mpf->SetLineWidth(2);
  pt_depend_const_mpf->SetLineColor(kRed+1);
  pt_depend_const_dijet->SetLineWidth(2);
  pt_depend_const_dijet->SetLineColor(kBlue+1);
  pt_depend_const_mpf->SetLineStyle(2);
  pt_depend_const_dijet->SetLineStyle(2);
  pt_depend_logpt_mpf->SetLineWidth(2);
  pt_depend_logpt_mpf->SetLineColor(kRed+1);
  pt_depend_logpt_dijet->SetLineWidth(2);
  pt_depend_logpt_dijet->SetLineColor(kBlue+1);


  kfsr_mpf->SetLineWidth(2);
  kfsr_mpf->SetLineColor(kRed+1);
  kfsr_dijet->SetLineWidth(2);
  kfsr_dijet->SetLineColor(kBlue+1);
  // kfsr_mpf_fit->SetLineWidth(2);
  // kfsr_mpf_fit->SetLineColor(kRed+1);
  // kfsr_dijet_fit->SetLineWidth(2);
  // kfsr_dijet_fit->SetLineColor(kBlue+1);


//   // Draw results
  TH1D *h = new TH1D("h",";|#eta|;Relative correction",41,0,5.191);
  // h->SetMaximum(1.35);
  // h->SetMinimum(0.95);
  h->SetMaximum(1.2);
  h->SetMinimum(0.8);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  //  lumi_13TeV = "Run2015  2.1 fb^{-1}";
  lumi_13TeV = "Run2016  0.8 fb^{-1}";
  //lumi_13TeV = "589.3 pb^{-1}";
  bool kSquare = true;


  TLegend *leg1 = tdrLeg(0.17,0.19,0.40,0.40);
  leg1 -> AddEntry(pt_depend_const_mpf, "MPF Flat","L");
  leg1 -> AddEntry(pt_depend_const_dijet, "Pt Flat","L");
  leg1 -> AddEntry(pt_depend_logpt_mpf, "MPF Loglin","L");
  leg1 -> AddEntry(pt_depend_logpt_dijet, "Pt Loglin","L");
  leg1->Draw();

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);  
  tex->DrawLatex(0.45,0.87,JetDescrib);


  TCanvas *c2 = tdrCanvas("c2",h,4,10,kSquare);
  //  pt_depend_const_mpf->GetYaxis()->SetTitle("(R^{MC}/R^{data})_{#alpha<0.3}");
  TString alVal;
  alVal.Form("%0.2f\n",al_cut);
  TString altitle = "{#alpha<"+alVal+"}";
  TString axistitle = "(R^{MC}/R^{data})_";
  axistitle +=altitle;
  h->GetYaxis()->SetTitle(axistitle);
  h->GetYaxis()->SetRangeUser(0.81,1.15);
  pt_depend_const_mpf->Draw("E1 SAME");
  pt_depend_const_dijet->Draw("E1 SAME");
  pt_depend_logpt_mpf->Draw("E1 SAME");
  pt_depend_logpt_dijet->Draw("E1 SAME");
  line->Draw("SAME");
  leg1->Draw();
  tex->DrawLatex(0.45,0.87,JetDescrib);
  c2->SaveAs(path+"plots/Ratio_"+jettag+"_"+txttag+".pdf");


  //todo: continue here!
  TCanvas *c3 = tdrCanvas("c3",h,4,10,kSquare);
  h->GetYaxis()->SetTitle("k_{FSR}");
  h->GetYaxis()->SetRangeUser(0.81,1.15);
  kfsr_dijet_fit->Draw("E3 SAME");
  kfsr_mpf_fit->Draw("E3 SAME");
  kfsr_mpf->Draw("E1 SAME");
  kfsr_dijet->Draw("E1 SAME");
  line->Draw("SAME");

  //TLegend *leg2 = tdrLeg(0.17,0.49,0.40,0.80);
  TLegend *leg2 = tdrLeg(0.17,0.19,0.40,0.30);
  leg2 -> AddEntry(kfsr_mpf, "MPF","L");
  leg2 -> AddEntry(kfsr_dijet, "Pt","L");
  leg2->Draw();
  tex->DrawLatex(0.45,0.87,JetDescrib);
  c3->SaveAs(path+"plots/kFSR_"+jettag+"_"+txttag+".pdf");


  TCanvas *c4 = tdrCanvas("L2res_kFSRfit",h,4,10,kSquare);
  h->GetYaxis()->SetTitle("Relative correction");
  h->GetYaxis()->SetRangeUser(0.8,1.2);
  res_const_mpf_kfsrfit->Draw("E1 SAME");
  res_const_dijet_kfsrfit->Draw("E1 SAME");
  res_logpt_mpf_kfsrfit->Draw("E1 SAME");
  res_logpt_dijet_kfsrfit->Draw("E1 SAME");
  leg1->Draw();
  tex->DrawLatex(0.45,0.87,JetDescrib);
  line->Draw();
  c4->SaveAs(path+"plots/L2Res_kFSRfit_"+jettag+"_"+txttag+".pdf");

  //pt-dependence of L2Res (Money plot)

  TString var="";
  TH1D* res_logpt_mpf_kfsrfit_var[4];
  TH1D* res_logpt_dijet_kfsrfit_var[4];
  for(int i=0;i<4;i++){
    if(i==0) var="central"; 
    if(i==1) var="down";
    if(i==2) var="up";
    if(i==3) var="doubleup";
    //    cout<<var<<endl; 
    TFile* f_Res_mpf_var = new TFile(path+"Histo_Res_MPF_L1_"+txttag+"_"+jettag+"_"+var+".root","READ");   
    TFile* f_Res_dijet_var = new TFile(path+"Histo_Res_DiJet_L1_"+txttag+"_"+jettag+"_"+var+".root","READ");  
    res_logpt_mpf_kfsrfit_var[i] = (TH1D*)f_Res_mpf_var->Get("res_logpt_mpf");
    res_logpt_dijet_kfsrfit_var[i] = (TH1D*)f_Res_dijet_var->Get("res_logpt_dijet");
  }

  TCanvas *c5 = tdrCanvas("L2res_logpt_MPF_kFSRfit_ptDepend",h,4,10,kSquare);
  res_logpt_mpf_kfsrfit->SetLineColor(kBlack);
  res_logpt_mpf_kfsrfit->Draw("E1 SAME");
  for(int i=0;i<4;i++){
    res_logpt_mpf_kfsrfit_var[i]->SetLineColor(kRed-3*i);
    res_logpt_mpf_kfsrfit_var[i]->SetMarkerColor(kRed-3*i);
    res_logpt_mpf_kfsrfit_var[i]->SetMarkerStyle(20+i);
    res_logpt_mpf_kfsrfit_var[i]->Draw("E1 SAME"); 
  }

  leg2 = tdrLeg(0.17,0.19,0.40,0.42); 
  leg2 -> AddEntry(res_logpt_mpf_kfsrfit_var[1] , "60 GeV","LP");  
  leg2 -> AddEntry(res_logpt_mpf_kfsrfit_var[0] , "120 GeV","LP");
  leg2 -> AddEntry(res_logpt_mpf_kfsrfit_var[2] , "240 GeV","LP"); 
  leg2 -> AddEntry(res_logpt_mpf_kfsrfit_var[3] , "480 GeV","LP");
  leg2 -> AddEntry(res_logpt_mpf_kfsrfit, "Nominal","LP");
  leg2->Draw();              
  tex->DrawLatex(0.45,0.87,JetDescrib);      
  c5->SaveAs(path+"plots/L2Res_logpt_MPF_kFSRfit_"+jettag+"_"+txttag+".pdf");
  TCanvas *c6 = tdrCanvas("L2res_logpt_DiJet_kFSRfit_ptDepend",h,4,10,kSquare);
  res_logpt_dijet_kfsrfit->SetLineColor(kBlack);
  res_logpt_dijet_kfsrfit->Draw("E1 SAME");
  for(int i=0;i<4;i++){
    res_logpt_dijet_kfsrfit_var[i]->SetLineColor(kBlue+3*i);
    res_logpt_dijet_kfsrfit_var[i]->SetMarkerColor(kBlue+3*i);
    res_logpt_dijet_kfsrfit_var[i]->SetMarkerStyle(20+i);
    res_logpt_dijet_kfsrfit_var[i]->Draw("E1 SAME");  

  }
  leg2 = tdrLeg(0.17,0.19,0.40,0.42); 
  leg2 -> AddEntry(res_logpt_dijet_kfsrfit_var[1] , "60 GeV","LP");  
  leg2 -> AddEntry(res_logpt_dijet_kfsrfit_var[0] , "120 GeV","LP");
  leg2 -> AddEntry(res_logpt_dijet_kfsrfit_var[2] , "240 GeV","LP"); 
  leg2 -> AddEntry(res_logpt_dijet_kfsrfit_var[3] , "480 GeV","LP");
  leg2 -> AddEntry(res_logpt_dijet_kfsrfit, "Nominal","LP");
  leg2->Draw();
  tex->DrawLatex(0.45,0.87,JetDescrib);                    
  c6->SaveAs(path+"plots/L2Res_logpt_DiJet_kFSRfit_"+jettag+"_"+txttag+".pdf");

  /*
  //use kFSR values in bins
  TCanvas *c1 = tdrCanvas("L2res_kFSRbins",h,4,10,kSquare);
  TH1D* pt_depend_const_mpf_kFSRbins = (TH1D*)pt_depend_const_mpf->Clone();
  TH1D* pt_depend_const_dijet_kFSRbins = (TH1D*)pt_depend_const_dijet->Clone();
  TH1D* pt_depend_logpt_mpf_kFSRbins = (TH1D*)pt_depend_logpt_mpf->Clone();
  TH1D* pt_depend_logpt_dijet_kFSRbins = (TH1D*)pt_depend_logpt_dijet->Clone();

  pt_depend_const_mpf_kFSRbins->Multiply(kfsr_mpf);
  pt_depend_logpt_mpf_kFSRbins->Multiply(kfsr_mpf);
  pt_depend_const_dijet_kFSRbins->Multiply(kfsr_dijet);
  pt_depend_logpt_dijet_kFSRbins->Multiply(kfsr_dijet);
  pt_depend_const_mpf_kFSRbins->SetName("res_const_mpf");
  pt_depend_logpt_mpf_kFSRbins->SetName("res_logpt_mpf");
  pt_depend_const_dijet_kFSRbins->SetName("res_const_dijet");
  pt_depend_logpt_dijet_kFSRbins->SetName("res_logpt_dijet");

  h->GetYaxis()->SetTitle("Relative correction");
  h->GetYaxis()->SetRangeUser(0.81,1.15);
  pt_depend_const_mpf_kFSRbins->Draw("E1 SAME");
  pt_depend_const_dijet_kFSRbins->Draw("E1 SAME");
  pt_depend_logpt_mpf_kFSRbins->Draw("E1 SAME");
  pt_depend_logpt_dijet_kFSRbins->Draw("E1 SAME");
  leg1->Draw();
  tex->DrawLatex(0.47,0.87,JetDescrib);
  line->Draw();
  c1->SaveAs(path+"plots/L2Res_kFSRbins_"+jettag+"_"+txttag+"_"+variation+".pdf");

    // use kFSR fit function instead of values in bins
  TCanvas *c4 = tdrCanvas("c4",h,4,10,kSquare);
  TH1D* pt_depend_const_mpf_kFSR = (TH1D*)pt_depend_const_mpf->Clone();
  TH1D* pt_depend_const_dijet_kFSR = (TH1D*)pt_depend_const_dijet->Clone();
  TH1D* pt_depend_logpt_mpf_kFSR = (TH1D*)pt_depend_logpt_mpf->Clone();
  TH1D* pt_depend_logpt_dijet_kFSR = (TH1D*)pt_depend_logpt_dijet->Clone();
  pt_depend_const_mpf_kFSR->SetName("res_const_mpf_kFSRfit");                                                            
  pt_depend_logpt_mpf_kFSR->SetName("res_logpt_mpf_kFSRfit");
  pt_depend_const_dijet_kFSR->SetName("res_const_dijet_kFSRfit");                                                                     
  pt_depend_logpt_dijet_kFSR->SetName("res_logpt_dijet_kFSRfit"); 
  pt_depend_const_mpf_kFSR->Multiply(hint_mpf);
  pt_depend_logpt_mpf_kFSR->Multiply(hint_mpf);
  pt_depend_const_dijet_kFSR->Multiply(hint_dijet);
  pt_depend_logpt_dijet_kFSR->Multiply(hint_dijet);

  pt_depend_const_mpf_kFSR->GetYaxis()->SetTitle("Relative correction");
  pt_depend_const_mpf_kFSR->GetYaxis()->SetTitleSize(0.05);
  pt_depend_const_mpf_kFSR->GetYaxis()->SetRangeUser(0.81,1.15);
  pt_depend_const_mpf_kFSR->GetXaxis()->SetTitle("|#eta|");
  pt_depend_const_mpf_kFSR->GetXaxis()->SetTitleSize(0.05);
  pt_depend_const_mpf_kFSR->Draw("E1 SAME");
  pt_depend_const_dijet_kFSR->Draw("E1 SAME");
  pt_depend_logpt_mpf_kFSR->Draw("E1 SAME");
  pt_depend_logpt_dijet_kFSR->Draw("E1 SAME");
  leg1->Draw();
  tex->DrawLatex(0.47,0.87,JetDescrib);
  line->Draw();
  c4->SaveAs(path+"plots/L2Res_kFSRfit_"+jettag+"_"+txttag+"_"+variation+".pdf");

  //ROOT output file with all histrograms 
  TFile* outputfilenorm = new TFile(path+"/Histos_Res_L1_"+jettag+"_"+txttag+"_"+variation+".root","RECREATE");
  //  consmpf->Write();
  //  consdijet->Write();
  pt_depend_const_mpf->Write();
  pt_depend_const_dijet->Write();
  //  logptmpf->Write();
  //  logptdijet->Write();
  pt_depend_logpt_mpf->Write();
  pt_depend_logpt_dijet->Write();
  kfsr_mpf->Write();
  kfsr_dijet->Write();
  pt_depend_const_mpf_kFSRbins->Write();
  pt_depend_logpt_mpf_kFSRbins->Write();
  pt_depend_const_dijet_kFSRbins->Write();
  pt_depend_logpt_dijet_kFSRbins->Write();
  pt_depend_const_mpf_kFSR->Write();
  pt_depend_logpt_mpf_kFSR->Write();
  pt_depend_const_dijet_kFSR->Write();
  pt_depend_logpt_dijet_kFSR->Write(); 
  outputfilenorm->Write();
  outputfilenorm->Print();
  outputfilenorm->Close();

  //txt file to be used in JEC
  ofstream output, output_loglin, uncerts, uncerts_loglin;                                                                                      
  output.open(path+"output/Spring16_25ns_MPF_FLAT_L2Residual_"+txttag+"_"+jettag+".txt");                                         
  output_loglin.open(path+"output/Spring16_25ns_MPF_LOGLIN_L2Residual_"+txttag+"_"+jettag+".txt");                                
  uncerts.open(path+"output/Spring16_25ns_MPF_FLAT_L2Residual_"+txttag+"_"+jettag+".txt.STAT");                                   
  uncerts_loglin.open(path+"output/Spring16_25ns_MPF_LOGLIN_L2Residual_"+txttag+"_"+jettag+".txt.STAT");                          
  output  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-\
0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(2\
08))/208))) Correction L2Relative}" << endl;                                                                                                      
  output_loglin  << "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*\
pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath\
::Log(208))/208))) Correction L2Relative}" << endl;                                                                                               
  uncerts << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}" << endl;                                                                           
  uncerts_loglin << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}" << endl; 

  for (int j=n_eta-1; j>0; --j){                                                                                                                
    output << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   55   " << ptave_data[j-1]-\
      >FindLastBinAbove(0.)*10 << "   " << 1/flat << " " << hist_kfsr_mpf->GetBinContent(j) << "   " << f2[j-1]->GetParameter(0) << " 0   1 0.0000\
 0.0" << endl;                                                                                                                                    
    output_loglin << fixed << std::setprecision(6)  << "  -" << eta_range[j]<< " -" << eta_range[j-1] << "   11   10 6500   55   " << ptave_dat\
      a[j-1]->FindLastBinAbove(0.)*10 << "   " << 1/loglin << " " << hist_kfsr_mpf->GetBinContent(j) << "   " << f1[j-1]->GetParameter(0) << " " <\
      < f1[j-1]->GetParameter(1) << "   1 0.0000 0.0" << endl;                                                                                          
    uncerts << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " << ptave_data[j-1]->FindLastBi\
      nAbove(0.)*10 << "       " << sqrt( pow(f2[j-1]->GetParameter(0)*hist_kfsr_mpf->GetBinError(j),2)+pow( hist_kfsr_mpf->GetBinContent(j)*f2[j-1]->G\
													     etParError(0),2)  ) / flat  << endl;                                                                                                         
    uncerts_loglin << fixed << std::setprecision(6) << "-" << eta_range[j] << " -" << eta_range[j-1] << "  3    55    " << ptave_data[j-1]->Fin\
      dLastBinAbove(0.)*10 << "       " << sqrt(  pow(f1[j-1]->GetParameter(0)+f1[j-1]->GetParameter(1)*TMath::Log(ptave_data[j-1]->GetMean()),2)*pow(h\
																		      ist_kfsr_mpf->GetBinError(j),2)  + pow(hist_kfsr_mpf->GetBinContent(j),2)*(Vcov[0][j-1]+Vcov[1][j-1]*pow(TMath::Log(ptave_data[j-1]->GetMean()),2\
																															       )+2*Vcov[2][j-1]*TMath::Log(ptave_data[j-1]->GetMean()))) / loglin << endl;                                                                  
  }                

  */
}

//  LocalWords:  pt_depend_const_
