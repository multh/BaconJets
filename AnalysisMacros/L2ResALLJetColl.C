#include "header.h"
//#include "tdrstyle_mod14.C"
#include "tdrstyle_mod15.C"

void L2ResALLJetColl(TString path, TString txttag, TString jettag, TString tag){
//   // Draw results
  TH1D *h = new TH1D("h",";|#eta|;Relative correction",41,0,5.191);
  h->SetMaximum(1.2);
  h->SetMinimum(0.8);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  //  lumi_13TeV = "Run2015  2.1 fb^{-1}";
  lumi_13TeV = "Run2016  0.8 fb^{-1}";
  //lumi_13TeV = "589.3 pb^{-1}";
  bool kSquare = true;
  TLegend *leg1 = tdrLeg(0.17,0.19,0.35,0.4);
  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,5.191,1);
  // tag for the jet cone size                                             
  TString jettag_v[4]={"AK4PFchs","AK4PFpuppi","AK8PFchs","AK8PFpuppi"};
  TCanvas *c2 = tdrCanvas("c2",h,4,10,kSquare);
  for(int j=0;j<4;j++){
    TString jettag=jettag_v[j];
    TString JetDescrib;                                                                                                                          
  if (jettag=="AK4PFchs") JetDescrib = "AK4CHS";                                                                             
  if (jettag=="AK4PFpuppi") JetDescrib = "AK4PUPPI";
  if (jettag=="AK8PFchs") JetDescrib = "AK8CHS"; 
  if (jettag=="AK8PFpuppi") JetDescrib = "AK8PUPPI"; 
  TFile* f_Res_mpf = new TFile(path+"Histo_Res_MPF_L1_"+txttag+"_"+jettag+tag+".root","READ");   
  //  TFile* f_Res_dijet = new TFile(path+"Histo_Res_DiJet_L1_"+txttag+"_"+jettag+tag+".root","READ");  

  //plot results for "nominal" variation

  //get L2Res hists for MPF and pt balance
  TH1D* res_logpt_mpf_kfsrfit = (TH1D*)f_Res_mpf->Get("res_logpt_mpf");
  //  TH1D* res_logpt_dijet_kfsrfit = (TH1D*)f_Res_dijet->Get("res_logpt_dijet");
  res_logpt_mpf_kfsrfit->SetLineWidth(2);
  res_logpt_mpf_kfsrfit->SetLineColor(1+j);
  //  res_logpt_dijet_kfsrfit->SetLineWidth(2);
  //  res_logpt_dijet_kfsrfit->SetLineColor(kBlue+1+j);
  leg1 -> AddEntry(res_logpt_mpf_kfsrfit,JetDescrib,"L");
  res_logpt_mpf_kfsrfit->Draw("E SAME");  
}
  leg1->Draw();

  c2->SaveAs(path+"plots/L2_MPF_LOGLIN_nominal_ALL_jets_collections.pdf");
  c2->SaveAs(path+"plots/L2_MPF_LOGLIN_nominal_ALL_jets_collections.root");

}
