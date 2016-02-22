// plot control plots of important parameters for L2Res determination

#include "header.h"


TString ToStringC(int num) {
  ostringstream start;
  start<<num;
  TString start1=start.str();
  return start1;
}

void Control_Plots(TString path, TFile* datafile, TFile* MCfile, TString GenName="pythia8"){
  TString dirName = "Selection";
  //  TString dirName = "JetMatching";
  gStyle->SetOptFit(0);
  TH1F *pt_jet1_DATA,*pt_jet2_DATA,*pt_jet3_DATA;
  TH1F *pt_jet1_MC,*pt_jet2_MC,*pt_jet3_MC;
  TH1F *eta_jet1_DATA,*eta_jet2_DATA,*eta_jet3_DATA;
  TH1F *eta_jet1_MC,*eta_jet2_MC,*eta_jet3_MC;
  pt_jet1_DATA = (TH1F*)datafile->Get(dirName+"/pt_1");
  pt_jet2_DATA = (TH1F*)datafile->Get(dirName+"/pt_2");
  pt_jet3_DATA = (TH1F*)datafile->Get(dirName+"/pt_3");
  pt_jet1_MC = (TH1F*)MCfile->Get(dirName+"/pt_1");
  pt_jet2_MC = (TH1F*)MCfile->Get(dirName+"/pt_2");
  pt_jet3_MC = (TH1F*)MCfile->Get(dirName+"/pt_3");
  eta_jet1_DATA = (TH1F*)datafile->Get(dirName+"/eta_1");
  eta_jet2_DATA = (TH1F*)datafile->Get(dirName+"/eta_2");
  eta_jet3_DATA = (TH1F*)datafile->Get(dirName+"/eta_3");
  eta_jet1_MC = (TH1F*)MCfile->Get(dirName+"/eta_1");
  eta_jet2_MC = (TH1F*)MCfile->Get(dirName+"/eta_2");
  eta_jet3_MC = (TH1F*)MCfile->Get(dirName+"/eta_3");

  TH1F *pt_ave_DATA = (TH1F*)datafile->Get(dirName+"/pt_ave");
  TH1F *pt_rel_DATA = (TH1F*)datafile->Get(dirName+"/pt_rel");
  TH1F *pt_ave_MC = (TH1F*)MCfile->Get(dirName+"/pt_ave");
  TH1F *pt_rel_MC = (TH1F*)MCfile->Get(dirName+"/pt_rel");
  TH1F *asym_DATA = (TH1F*)datafile->Get(dirName+"/asym");
  TH1F *asym_MC = (TH1F*)MCfile->Get(dirName+"/asym");
  TH1F *r_dijet_DATA = (TH1F*)datafile->Get(dirName+"/r_rel");
  TH1F *r_dijet_MC = (TH1F*)MCfile->Get(dirName+"/r_rel");
  TH1F *r_mpf_DATA = (TH1F*)datafile->Get(dirName+"/mpf");
  TH1F *r_mpf_MC = (TH1F*)MCfile->Get(dirName+"/mpf");
  TH1F *N_PV_DATA = (TH1F*)datafile->Get(dirName+"/N_PV");
  TH1F *N_PV_MC = (TH1F*)MCfile->Get(dirName+"/N_PV");

  TH2F *ptjet3_vs_alpha_DATA = (TH2F*)datafile->Get(dirName+"/ptjet3_vs_alpha");
  TH2F *ptjet3_vs_alpha_MC = (TH2F*)MCfile->Get(dirName+"/ptjet3_vs_alpha");
  TH2F *pt_ave_vs_alpha_DATA = (TH2F*)datafile->Get(dirName+"/pt_ave_vs_alpha");
  TH2F *pt_ave_vs_alpha_MC = (TH2F*)MCfile->Get(dirName+"/pt_ave_vs_alpha");
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // Plots

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  //create plots
  gStyle->SetOptStat(0);
  TCanvas* a = new TCanvas();
  a->Divide(3,2);
  a->cd(1);
  pt_jet1_DATA->SetMarkerStyle(20);
  pt_jet1_DATA->SetMarkerSize(0.5);
  pt_jet1_DATA->SetMarkerColor(1);
  pt_jet1_DATA->SetLineColor(1);
  pt_jet1_DATA->GetXaxis()->SetTitle("p_{T}, GeV");
  pt_jet1_DATA->Scale(1/pt_jet1_DATA->Integral());
  pt_jet1_DATA->GetYaxis()->SetRangeUser(0,0.1);
  pt_jet1_DATA->Draw();
  pt_jet1_MC->SetMarkerStyle(22);
  pt_jet1_MC->SetMarkerSize(0.5);
  pt_jet1_MC->SetMarkerColor(2);
  pt_jet1_MC->SetLineColor(2);
  pt_jet1_MC->Scale(1/pt_jet1_MC->Integral());
  pt_jet1_MC->Draw("same");
  TLegend *leg1;
  leg1 = new TLegend(0.63,0.68,0.9,0.85,"","brNDC");//x+0.1
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.045);
  leg1->SetFillColor(10);
  leg1->SetLineColor(1);
  leg1->SetTextFont(42);
  leg1->AddEntry(eta_jet3_DATA,"DATA","lp");
  leg1->AddEntry(eta_jet3_MC,"MC","lp");
  leg1->Draw();
  a->cd(2);
  pt_jet2_DATA->SetMarkerStyle(20);
  pt_jet2_DATA->SetMarkerSize(0.5);
  pt_jet2_DATA->SetMarkerColor(1);
  pt_jet2_DATA->SetLineColor(1);
  pt_jet2_DATA->GetXaxis()->SetTitle("p_{T}, GeV");
  pt_jet2_DATA->Scale(1/pt_jet2_DATA->Integral());
  pt_jet2_DATA->GetYaxis()->SetRangeUser(0,0.1);
  pt_jet2_DATA->Draw();
  pt_jet2_MC->SetMarkerStyle(22);
  pt_jet2_MC->SetMarkerSize(0.5);
  pt_jet2_MC->SetMarkerColor(2);
  pt_jet2_MC->SetLineColor(2);
  pt_jet2_MC->Scale(1/pt_jet2_MC->Integral());
  pt_jet2_MC->Draw("same");
  a->cd(3);
  // gPad->SetLogx();   gPad->SetLogy();
  pt_jet3_DATA->SetMarkerStyle(20);
  pt_jet3_DATA->SetMarkerSize(0.5);
  pt_jet3_DATA->SetMarkerColor(1);
  pt_jet3_DATA->SetLineColor(1);
  pt_jet3_DATA->GetXaxis()->SetTitle("p_{T}, GeV");
  pt_jet3_DATA->Scale(1/pt_jet3_DATA->Integral());
  pt_jet3_DATA->GetYaxis()->SetRangeUser(0,0.1);
  pt_jet3_DATA->Draw();
  pt_jet3_MC->SetMarkerStyle(22);
  pt_jet3_MC->SetMarkerSize(0.5);
  pt_jet3_MC->SetMarkerColor(2);
  pt_jet3_MC->SetLineColor(2);
  pt_jet3_MC->Scale(1/pt_jet3_MC->Integral());
  pt_jet3_MC->Draw("same");
 
  a->cd(4);
  eta_jet1_DATA->SetMarkerStyle(20);
  eta_jet1_DATA->SetMarkerSize(0.5);
  eta_jet1_DATA->SetMarkerColor(1);
  eta_jet1_DATA->SetLineColor(1);
  eta_jet1_DATA->GetXaxis()->SetTitle("#eta");
  eta_jet1_DATA->Scale(1./eta_jet1_DATA->Integral());
  eta_jet1_DATA->GetYaxis()->SetRangeUser(0,0.04);
  eta_jet1_DATA->Draw();
  eta_jet1_MC->SetMarkerStyle(22);
  eta_jet1_MC->SetMarkerSize(0.5);
  eta_jet1_MC->SetMarkerColor(2);
  eta_jet1_MC->SetLineColor(2);
  eta_jet1_MC->Scale(1./eta_jet1_MC->Integral());
  eta_jet1_MC->Draw("same");
  a->cd(5);
  eta_jet2_DATA->SetMarkerStyle(20);
  eta_jet2_DATA->SetMarkerSize(0.5);
  eta_jet2_DATA->SetMarkerColor(1);
  eta_jet2_DATA->SetLineColor(1);
  eta_jet2_DATA->GetXaxis()->SetTitle("#eta");
  eta_jet2_DATA->Scale(1./eta_jet2_DATA->Integral());
  eta_jet2_DATA->GetYaxis()->SetRangeUser(0,0.04);
  eta_jet2_DATA->Draw();
  eta_jet2_MC->SetMarkerStyle(22);
  eta_jet2_MC->SetMarkerSize(0.5);
  eta_jet2_MC->SetMarkerColor(2);
  eta_jet2_MC->SetLineColor(2);
  eta_jet2_MC->Scale(1./eta_jet2_MC->Integral());
  eta_jet2_MC->Draw("same");
  a->cd(6);
  eta_jet3_DATA->SetMarkerStyle(20);
  eta_jet3_DATA->SetMarkerSize(0.5);
  eta_jet3_DATA->SetMarkerColor(1);
  eta_jet3_DATA->SetLineColor(1);
  eta_jet3_DATA->GetXaxis()->SetTitle("#eta");
  eta_jet3_DATA->Scale(1./eta_jet3_DATA->Integral());
  eta_jet3_DATA->GetYaxis()->SetRangeUser(0,0.04);
  eta_jet3_DATA->Draw();
  eta_jet3_MC->SetMarkerStyle(22);
  eta_jet3_MC->SetMarkerSize(0.5);
  eta_jet3_MC->SetMarkerColor(2);
  eta_jet3_MC->SetLineColor(2);  
  eta_jet3_MC->Scale(1./eta_jet3_MC->Integral());
  eta_jet3_MC->Draw("same");
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045); 
  tex->DrawLatex(0.44,0.81,"MC ["+GenName+"]");
  tex->DrawLatex(0.44,0.71,"DATA [2.11fb^{-1} (13TeV)]");
  a->Print("plots/Control_Plots_"+dirName+"_"+GenName+"_jets.pdf");
  
  TCanvas* b = new TCanvas();
  b->Divide(3,2);
  b->cd(1);
  N_PV_DATA->SetMarkerStyle(20);
  N_PV_DATA->SetMarkerSize(0.5);
  N_PV_DATA->SetMarkerColor(1);
  N_PV_DATA->SetLineColor(1);
  //  N_PV_DATA->GetXaxis()->SetTitle("N PVtx");
  N_PV_DATA->Scale(1./N_PV_DATA->Integral());
  N_PV_DATA->GetYaxis()->SetRangeUser(0,0.15);
  N_PV_DATA->Draw();
  N_PV_MC->SetMarkerStyle(22);
  N_PV_MC->SetMarkerSize(0.5);
  N_PV_MC->SetMarkerColor(2);
  N_PV_MC->SetLineColor(2);
  N_PV_MC->Scale(1./N_PV_MC->Integral());
  N_PV_MC->Draw("same");
  TLegend *leg2;
  leg2 = new TLegend(0.63,0.68,0.9,0.86,"","brNDC");//x+0.1
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.045);
  leg2->SetFillColor(10);
  leg2->SetLineColor(1);
  leg2->SetTextFont(42);
  leg2->AddEntry(eta_jet3_DATA,"DATA","lp");
  leg2->AddEntry(eta_jet3_MC,"MC","lp");
  leg2->Draw();
  b->cd(2);
  pt_ave_DATA->SetMarkerStyle(20);
  pt_ave_DATA->SetMarkerSize(0.5);
  pt_ave_DATA->SetMarkerColor(1);
  pt_ave_DATA->SetLineColor(1);
  // pt_ave_DATA->Scale(1./pt_ave_DATA->Integral());
  //  pt_ave_DATA->GetYaxis()->SetRangeUser(0,0.1);
  pt_ave_DATA->Draw();
  pt_ave_MC->SetMarkerStyle(22);
  pt_ave_MC->SetMarkerSize(0.5);
  pt_ave_MC->SetMarkerColor(2);
  pt_ave_MC->SetLineColor(2);
  //  pt_ave_MC->Scale(1./pt_ave_MC->Integral());
  pt_ave_MC->Draw("same");
  b->cd(3);
  pt_rel_DATA->SetMarkerStyle(20);
  pt_rel_DATA->SetMarkerSize(0.5);
  pt_rel_DATA->SetMarkerColor(1);
  pt_rel_DATA->SetLineColor(1);
  pt_rel_DATA->Scale(1./pt_rel_DATA->Integral());
  pt_rel_DATA->GetYaxis()->SetRangeUser(0,0.15);
  pt_rel_DATA->Draw();
  pt_rel_MC->SetMarkerStyle(22);
  pt_rel_MC->SetMarkerSize(0.5);
  pt_rel_MC->SetMarkerColor(2);
  pt_rel_MC->SetLineColor(2);
  pt_rel_MC->Scale(1./pt_rel_MC->Integral());
  pt_rel_MC->Draw("same");
  //  gPad->SetLogy();
  b->cd(4);
  asym_DATA->SetMarkerStyle(20);
  asym_DATA->SetMarkerSize(0.5);
  asym_DATA->SetMarkerColor(1);
  asym_DATA->SetLineColor(1);
  asym_DATA->Scale(1./asym_DATA->Integral());
  asym_DATA->GetYaxis()->SetRangeUser(0,0.1);
  asym_DATA->Draw();
  asym_MC->SetMarkerStyle(22);
  asym_MC->SetMarkerSize(0.5);
  asym_MC->SetMarkerColor(2);
  asym_MC->SetLineColor(2);
  asym_MC->Scale(1./asym_MC->Integral());
  asym_MC->Draw("same");
  b->cd(5);
  r_dijet_DATA->SetMarkerStyle(20);
  r_dijet_DATA->SetMarkerSize(0.5);
  r_dijet_DATA->SetMarkerColor(1);
  r_dijet_DATA->SetLineColor(1);
  r_dijet_DATA->Scale(1./r_dijet_DATA->Integral());
  r_dijet_DATA->GetYaxis()->SetRangeUser(0,0.04);
  r_dijet_DATA->Draw();
  r_dijet_MC->SetMarkerStyle(22);
  r_dijet_MC->SetMarkerSize(0.5);
  r_dijet_MC->SetMarkerColor(2);
  r_dijet_MC->SetLineColor(2);
  r_dijet_MC->Scale(1./r_dijet_MC->Integral());
  r_dijet_MC->Draw("same");
  leg1->Draw();
  b->cd(6);
  r_mpf_DATA->SetMarkerStyle(20);
  r_mpf_DATA->SetMarkerSize(0.5);
  r_mpf_DATA->SetMarkerColor(1);
  r_mpf_DATA->SetLineColor(1);
  r_mpf_DATA->Scale(1./r_mpf_DATA->Integral());
  r_mpf_DATA->GetYaxis()->SetRangeUser(0,0.05);
  r_mpf_DATA->Draw();
  r_mpf_MC->SetMarkerStyle(22);
  r_mpf_MC->SetMarkerSize(0.5);
  r_mpf_MC->SetMarkerColor(2);
  r_mpf_MC->SetLineColor(2);
  r_mpf_MC->Scale(1./r_mpf_MC->Integral());
  r_mpf_MC->Draw("same");
  TLatex *tex2 = new TLatex();
  tex2->SetNDC();
  tex2->SetTextSize(0.045); 
  tex2->DrawLatex(0.44,0.81,"MC ["+GenName+"]");
  tex2->DrawLatex(0.44,0.71,"DATA [2.11fb^{-1} (13TeV)]");
  b->Print("plots/Control_Plots_"+dirName+"_"+GenName+"_dijet.pdf");
  
  TCanvas* c = new TCanvas();
  c->Divide(2,2);
  c->cd(1);
  ptjet3_vs_alpha_DATA->SetTitle("DATA");
  ptjet3_vs_alpha_DATA->GetXaxis()->SetTitle("#alpha");
  ptjet3_vs_alpha_DATA->GetYaxis()->SetTitle("p_{T} jet3, GeV");
  ptjet3_vs_alpha_DATA->GetXaxis()->SetTitleSize(0.05);
  ptjet3_vs_alpha_DATA->GetYaxis()->SetTitleSize(0.05);
  ptjet3_vs_alpha_DATA->Draw("colz");
  c->cd(2);
  ptjet3_vs_alpha_MC->SetTitle("DATA");
  ptjet3_vs_alpha_MC->GetXaxis()->SetTitle("#alpha");
  ptjet3_vs_alpha_MC->GetYaxis()->SetTitle("p_{T} jet3, GeV");
  ptjet3_vs_alpha_MC->GetXaxis()->SetTitleSize(0.05);
  ptjet3_vs_alpha_MC->GetYaxis()->SetTitleSize(0.05);
  ptjet3_vs_alpha_MC->SetTitle("MC");
  ptjet3_vs_alpha_MC->Draw("colz");
  c->cd(3);
  pt_ave_vs_alpha_DATA->SetTitle("DATA");
  pt_ave_vs_alpha_DATA->GetXaxis()->SetTitle("#alpha");
  pt_ave_vs_alpha_DATA->GetYaxis()->SetTitle("p^{ave}_{T}, GeV");
  pt_ave_vs_alpha_DATA->GetXaxis()->SetTitleSize(0.05);
  pt_ave_vs_alpha_DATA->GetYaxis()->SetTitleSize(0.05);
  pt_ave_vs_alpha_DATA->Draw("colz");
  c->cd(4);
  pt_ave_vs_alpha_MC->GetXaxis()->SetTitle("#alpha");
  pt_ave_vs_alpha_MC->GetYaxis()->SetTitle("p^{ave}_{T}, GeV");
  pt_ave_vs_alpha_MC->GetXaxis()->SetTitleSize(0.05);
  pt_ave_vs_alpha_MC->GetYaxis()->SetTitleSize(0.05);
  pt_ave_vs_alpha_MC->SetTitle("MC");
  pt_ave_vs_alpha_MC->Draw("colz");
  c->Print("plots/Control_Plots_"+dirName+"_"+GenName+"_2Dalpha.pdf");



  // Plot trigger efficiency plots
  //  TString     sTriggerThreshold[] = {"40","60", "80", "140", "200", "260", "320", "400","500"};
  TString     sTriggerThreshold[] = {"40","60", "200", "260", "320", "400","500"};
  Float_t scale[] = {4.66251,53.6801,1.43334,6.2315,3.5902,1.2033}; //normalisation: obtained from const fit of the middle of the trigger efficiency curves
  Float_t     fTriggerThreshold[10] = {56,78,100,168,232,300,366,453,562,600};
  TH1F *histo_DATA_A[9]; 
  TH1F *histo_DATA_B[9]; 
  TH1F *histo_DATA_eff[9]; 
  TCanvas* d = new TCanvas();
  TLegend *leg3;
  leg3 = new TLegend(0.13,0.58,0.5,0.88,"","brNDC");//x+0.1
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.045);
  leg3->SetFillColor(10);
  leg3->SetLineColor(1);
  leg3->SetTextFont(42);
  TString trigName = "HLT_DiPFJetAve";
  d->Divide(3,2);
  for(int i=0; i < 6; i++){
    d->cd(i+1);
    histo_DATA_B[i] = (TH1F*)datafile-> Get(dirName+"/pt_ave_hltDiPFJetAve"+sTriggerThreshold[i]);
    histo_DATA_A[i] = (TH1F*)datafile-> Get(dirName+"/pt_ave_hltDiPFJetAve"+sTriggerThreshold[i+1]);
    histo_DATA_eff[i] = (TH1F*)histo_DATA_A[i]->Clone();
    histo_DATA_eff[i]->Divide(histo_DATA_B[i]);
    //    TF1 *fa1 = new TF1("fa1","pol0",300,400); 
    //    TF1 *fa1 = new TF1("fa1","pol0",400,500); 
    //    TF1 *fa1 = new TF1("fa1","pol0",500,600); 
    // TF1 *fa1 = new TF1("fa1","pol0",550,600); 
    // histo_DATA_eff[i]->Fit(fa1,"R");
    histo_DATA_eff[i]->Scale(1./scale[i]);
    histo_DATA_eff[i]->SetMarkerColor(kRed-1+i);
    histo_DATA_eff[i]->SetLineColor(kRed-1+i);
    histo_DATA_eff[i]->SetMarkerStyle(20);
    //    TF1 *fstep = new TF1("fstep","0.5*(1+TMath::Erf((x-[0])/(TMath::Sqrt(2)*[1])))",50,200); 
    //    TF1 *fstep = new TF1("fstep","TMath::Erf(x-[0])+[1]",50,200); 
    //    TF1 *fstep = new TF1("fstep","TMath::Erf(x)",0,600); 
    //  histo_DATA_eff[i]->Fit(fstep,"R");
    histo_DATA_eff[i]->GetYaxis()->SetRangeUser(-0.1,3.);
    histo_DATA_eff[i]->Draw();
  }
  d->Print("plots/Control_Plots_"+dirName+"_"+GenName+"_TriggerEff.pdf");
  TCanvas* e = new TCanvas();
  for(int i=0; i < 6; i++){
    if(i==0){
      histo_DATA_eff[i]->GetYaxis()->SetRangeUser(-0.1,3.);
      histo_DATA_eff[i]->SetTitle("");
      histo_DATA_eff[i]->Draw();
    }
    else{
      histo_DATA_eff[i]->Draw("same");
    }
    leg3->AddEntry(histo_DATA_eff[i],trigName+sTriggerThreshold[i+1],"lp");
  }
  leg3->Draw();
  e->Print("plots/Control_Plots_"+dirName+"_"+GenName+"_TriggerEffOnePlot.pdf");


 TCanvas* f = new TCanvas();
 TH2F *mpf_vs_etaProbe_DATA = (TH2F*)datafile->Get(dirName+"/mpf_vs_etaProbe");
 TH2F *mpf_vs_etaProbe_MC = (TH2F*)MCfile->Get(dirName+"/mpf_vs_etaProbe");
 TH2F *r_rel_vs_etaProbe_DATA = (TH2F*)datafile->Get(dirName+"/r_rel_vs_etaProbe");
 TH2F *r_rel_vs_etaProbe_MC = (TH2F*)MCfile->Get(dirName+"/r_rel_vs_etaProbe");
 f->Divide(2,2);
 f->cd(1);
 r_rel_vs_etaProbe_DATA->SetTitle("DATA");
 r_rel_vs_etaProbe_DATA->Draw("colz");
 f->cd(2);
 r_rel_vs_etaProbe_MC->SetTitle("MC");
 r_rel_vs_etaProbe_MC->Draw("colz");
 f->cd(3);
 mpf_vs_etaProbe_DATA->SetTitle("DATA");
 mpf_vs_etaProbe_DATA->Draw("colz");
 f->cd(4);
 mpf_vs_etaProbe_MC->SetTitle("MC");
 mpf_vs_etaProbe_MC->Draw("colz");
 f->Print("plots/Control_Plots_"+dirName+"_"+GenName+"_ResponsesVsEta.pdf");

 TCanvas* g = new TCanvas();
 TH2F *pt_ave_vs_etaProbe_DATA = (TH2F*)datafile->Get(dirName+"/pt_ave_vs_etaProbe");
 TH2F *pt_ave_vs_etaProbe_MC = (TH2F*)MCfile->Get(dirName+"/pt_ave_vs_etaProbe");


 TLine *lineEta[2*n_eta];
for(int i=0;i<n_eta;i++){
  lineEta[i] = new TLine(-1*eta_bins[i],0,-1*eta_bins[i],1000);
   lineEta[i]->SetLineColor(2);
   lineEta[i]->SetLineStyle(2);
 }
 for(int i=n_eta;i<2*n_eta;i++){
   lineEta[i] = new TLine(eta_bins[i-n_eta],0,eta_bins[i-n_eta],1000);
   lineEta[i]->SetLineColor(2);
   lineEta[i]->SetLineStyle(2);
   lineEta[i]->SetLineWidth(2);
 }

 TLine *linePt[n_pt];
 for(int i=0;i<n_pt;i++){
   linePt[i] = new TLine(-5.2,pt_bins[i],5.2,pt_bins[i]);
   linePt[i]->SetLineColor(2);
   linePt[i]->SetLineStyle(2);
   linePt[i]->SetLineWidth(2);
 }

 g->Divide(2,1);
 g->cd(1);
 pt_ave_vs_etaProbe_DATA->SetTitle("DATA");
 pt_ave_vs_etaProbe_DATA->Draw("colz");
for(int i=0;i<2*n_eta;i++)
  lineEta[i]->Draw("same");
for(int i=0;i<n_pt;i++)
   linePt[i]->Draw("same");
 g->cd(2);
 pt_ave_vs_etaProbe_MC->SetTitle("MC");
 pt_ave_vs_etaProbe_MC->Draw("colz");
for(int i=0;i<2*n_eta;i++)
  lineEta[i]->Draw("same");
for(int i=0;i<n_pt;i++)
   linePt[i]->Draw("same");

  g->Print("plots/Control_Plots_"+dirName+"_"+GenName+"_Pt_aveVsEtaProbe.pdf");

}
