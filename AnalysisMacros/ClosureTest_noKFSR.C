#include "header.h"
#include "tdrstyle_mod14.C"

void ClosureTest_noKFSR(TString path_1="/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_noReweight_AllTriggers_TTree/",TString path_2="/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV2_noReweight_AllTriggers/", double al_cut=0.2){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  

  // get the (R_{MC}/R_{DATA}) root file for MPF and pt balance
  TFile* rrcomb_1 = new TFile(path_1+"Histos_Res_norm_L1.root","READ"); 
  TH1D* rrconstcomb_1 = (TH1D*)rrcomb_1->Get("ptave_const_comb");
  TH1D* rrlogcomb_1 = (TH1D*)rrcomb_1->Get("ptave_logpt_comb");
  TFile* rrcomb_2 = new TFile(path_2+"Histos_Res_norm_L1.root","READ"); 
  TH1D* rrconstcomb_2 = (TH1D*)rrcomb_2->Get("ptave_const_comb");
  TH1D* rrlogcomb_2 = (TH1D*)rrcomb_2->Get("ptave_logpt_comb");


  
  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0.,1,5.191,1);
  rrconstcomb_1->SetLineWidth(2);
  rrconstcomb_2->SetLineWidth(2);
  rrconstcomb_1->SetLineStyle(2);
  rrconstcomb_2->SetLineStyle(2);
  rrconstcomb_1->SetLineColor(kGreen+1);
  rrconstcomb_1->SetMarkerSize(1.1);
  rrconstcomb_2->SetMarkerSize(1.2);
  rrconstcomb_1->SetMarkerStyle(25);
  rrconstcomb_1->SetMarkerColor(kGreen+1);
  rrconstcomb_2->SetLineColor(kGreen+4);
  rrconstcomb_2->SetMarkerStyle(26);
  rrconstcomb_2->SetMarkerColor(kGreen+4);
  rrlogcomb_1->SetLineWidth(2);
  rrlogcomb_2->SetLineWidth(2);
  rrlogcomb_1->SetLineColor(kGreen+1);
  rrlogcomb_1->SetMarkerSize(1.1);
  rrlogcomb_2->SetMarkerSize(1.2);
  rrlogcomb_1->SetMarkerStyle(21);
  rrlogcomb_1->SetMarkerColor(kGreen+1);
  rrlogcomb_2->SetLineColor(kGreen+4);
  rrlogcomb_2->SetMarkerStyle(22);
  rrlogcomb_2->SetMarkerColor(kGreen+4);
  // rrconstcomb_2->Draw();
  // rrconstcomb_1->Draw("same");

  // Draw results
  TH1D *h = new TH1D("h",";|#eta|;Relative correction",41,0,5.191);
  // h->SetMaximum(1.35);
  // h->SetMinimum(0.95);
  h->SetMaximum(1.2);
  h->SetMinimum(0.8);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  lumi_13TeV = "2.11 fb^{-1}";
  bool kSquare = true;
  //  bool kSquare = false;
  TCanvas *c2 = tdrCanvas("c2",h,4,0,kSquare);
  // TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  // pad1->SetBottomMargin(0); // Upper and lower plot are joined
  // pad1->SetGridx();         // Vertical grid
  // pad1->Draw();             // Draw the upper pad: pad1
  // pad1->cd();               // pad1 becomes the current pad
  TString alVal;
  alVal.Form("%0.2f\n",al_cut);
  TString altitle = "{#alpha<"+alVal+"}";
  TString axistitle = "(R^{MC}/R^{data})_";
  axistitle +=altitle;
  rrconstcomb_2->GetYaxis()->SetTitle(axistitle);
  rrconstcomb_2->GetYaxis()->SetTitleSize(0.05);
  rrconstcomb_2->GetYaxis()->SetRangeUser(0.95,1.15);
  rrconstcomb_2->GetXaxis()->SetTitle("|#eta|");
  rrconstcomb_2->GetXaxis()->SetTitleSize(0.05);
  rrconstcomb_2->Draw("E1");
  rrconstcomb_1->Draw("E1 SAME");
  rrlogcomb_2->Draw("E1 SAME");
  rrlogcomb_1->Draw("E1 SAME");
  line->Draw("SAME");
  TLegend *leg2 = tdrLeg(0.17,0.62,0.40,0.82);
  leg2 -> AddEntry(rrlogcomb_1, "Fall15_25nsV1 COMB LOGLIN","LP");
  leg2 -> AddEntry(rrconstcomb_1, "Fall15_25nsV1 COMB FLAT","LP");
  leg2 -> AddEntry(rrlogcomb_2, "Fall15_25nsV2 COMB LOGLIN","LP");
  leg2 -> AddEntry(rrconstcomb_2, "Fall15_25nsV2 COMB FLAT","LP");
  leg2->Draw();

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  tex->DrawLatex(0.47,0.87,"Anti-k_{t} R = 0.4, PF+CHS");

  // c2->cd();          // Go back to the main canvas before defining pad2
  // TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  // pad2->SetTopMargin(0);
  // pad2->SetBottomMargin(0.2);
  // //  pad2->SetGridx(); // vertical grid
  // pad2->Draw();
  // pad2->cd();       // pad2 becomes the current pad
  // TH1D *hratio = (TH1D*)rrconstcomb_2->Clone("hratio");
  // hratio->Divide(rrconstcomb_1);
  // hratio->SetMarkerStyle(20);
  // hratio->SetMarkerColor(1);
  // hratio->SetLineColor(1);
  // hratio->GetYaxis()->SetRangeUser(0.995,1.005);
  // hratio->GetYaxis()->SetTitleSize(0.05);
  // //  hratio->GetYaxis()->SetLabelSize(15);
  // hratio->Draw("ep");
  // hratio->GetYaxis()->SetTitle("ratio V2/V1");
  c2->SaveAs(path_2+"plots/ClosureTest_25ns_2p11fb_nokFSR.pdf");
}
