// plot control plots of important parameters for L2Res determination

#include "header.h"


TString ToStringC(int num) {
  ostringstream start;
  start<<num;
  TString start1=start.str();
  return start1;
}

void Control_Plots(TString path, TFile* datafile, TFile* MCfile, TString GenName="pythia8",TFile* Realdatafile=NULL){
  // TString DATAtitle = "DATA";
  // TString MCtitle = "MC";
  //  TString MCtitle = "Herwigpp";
  TString MCtitle = "Pythia8";
 TString DATAtitle = "Pythia8";

  // // TString DATAtitle = "Herwigpp";
  // // TString MCtitle = "Pythia8"; 
 //  TString dirName = "Selection";
  //  TString dirName = "JetMatching";
  //  TString dirName = "diJet";
 TString dirName = "noCuts";
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

  TH1F *Njets_DATA = (TH1F*)datafile->Get(dirName+"/N_jets");
  TH1F *Njets_MC = (TH1F*)MCfile->Get(dirName+"/N_jets");

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
  TH1F *nPU_DATA = (TH1F*)datafile->Get(dirName+"/nPu");
  TH1F *nPU_MC = (TH1F*)MCfile->Get(dirName+"/nPu");

  TH2F *ptjet3_vs_alpha_DATA = (TH2F*)datafile->Get(dirName+"/ptjet3_vs_alpha");
  TH2F *ptjet3_vs_alpha_MC = (TH2F*)MCfile->Get(dirName+"/ptjet3_vs_alpha");
  TH2F *pt_ave_vs_alpha_DATA = (TH2F*)datafile->Get(dirName+"/pt_ave_vs_alpha");
  TH2F *pt_ave_vs_alpha_MC = (TH2F*)MCfile->Get(dirName+"/pt_ave_vs_alpha");

 TH1F *Eta_pos_DATA = (TH1F*)datafile->Get(dirName+"/eta_probe_pos");
 TH1F *Eta_neg_DATA = (TH1F*)datafile->Get(dirName+"/eta_probe_neg");
 TH1F *Eta_assym_top_DATA = (TH1F*)Eta_pos_DATA->Clone();
 Eta_assym_top_DATA->Add(Eta_neg_DATA,-1);
 TH1F *Eta_assym_bot_DATA = (TH1F*)Eta_pos_DATA->Clone();
 Eta_assym_bot_DATA->Add(Eta_neg_DATA,+1);
 Eta_assym_top_DATA->Divide(Eta_assym_bot_DATA);

 // Eta_neg_DATA->Add(Eta_pos_DATA,-1);
 // Eta_neg_DATA->Divide(Eta_pos_DATA);
 TH1F *Eta_pos_MC = (TH1F*)MCfile->Get(dirName+"/eta_probe_pos");
 TH1F *Eta_neg_MC = (TH1F*)MCfile->Get(dirName+"/eta_probe_neg");
 TH1F *Eta_assym_top_MC = (TH1F*)Eta_pos_MC->Clone();
 Eta_assym_top_MC->Add(Eta_neg_MC,-1);
 TH1F *Eta_assym_bot_MC = (TH1F*)Eta_pos_MC->Clone();
 Eta_assym_bot_MC->Add(Eta_neg_MC,+1);
 Eta_assym_top_MC->Divide(Eta_assym_bot_MC);
 // Eta_neg_MC->Add(Eta_pos_MC,-1);
 // Eta_neg_MC->Divide(Eta_pos_MC);
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
  // pt_jet1_DATA->Scale(1/pt_jet1_DATA->Integral());
  //  pt_jet1_DATA->GetYaxis()->SetRangeUser(0,0.1);
  pt_jet1_DATA->Draw();
  pt_jet1_MC->SetMarkerStyle(22);
  pt_jet1_MC->SetMarkerSize(0.5);
  pt_jet1_MC->SetMarkerColor(2);
  pt_jet1_MC->SetLineColor(2);
  //  pt_jet1_MC->Scale(1/pt_jet1_MC->Integral());
  pt_jet1_MC->Draw("same");
  TLegend *leg1;
  leg1 = new TLegend(0.63,0.68,0.9,0.85,"","brNDC");//x+0.1
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.045);
  leg1->SetFillColor(10);
  leg1->SetLineColor(1);
  leg1->SetTextFont(42);
  leg1->AddEntry(pt_jet1_DATA,DATAtitle,"lp");
  leg1->AddEntry(pt_jet1_MC,MCtitle,"lp");
  leg1->Draw();
  a->cd(2);
  pt_jet2_DATA->SetMarkerStyle(20);
  pt_jet2_DATA->SetMarkerSize(0.5);
  pt_jet2_DATA->SetMarkerColor(1);
  pt_jet2_DATA->SetLineColor(1);
  pt_jet2_DATA->GetXaxis()->SetTitle("p_{T}, GeV");
  //  pt_jet2_DATA->Scale(1/pt_jet2_DATA->Integral());
  //  pt_jet2_DATA->GetYaxis()->SetRangeUser(0,0.1);
  pt_jet2_DATA->Draw();
  pt_jet2_MC->SetMarkerStyle(22);
  pt_jet2_MC->SetMarkerSize(0.5);
  pt_jet2_MC->SetMarkerColor(2);
  pt_jet2_MC->SetLineColor(2);
  //  pt_jet2_MC->Scale(1/pt_jet2_MC->Integral());
  pt_jet2_MC->Draw("same");
  a->cd(3);
  // gPad->SetLogx();   gPad->SetLogy();
  pt_jet3_DATA->SetMarkerStyle(20);
  pt_jet3_DATA->SetMarkerSize(0.5);
  pt_jet3_DATA->SetMarkerColor(1);
  pt_jet3_DATA->SetLineColor(1);
  pt_jet3_DATA->GetXaxis()->SetTitle("p_{T}, GeV");
  //  pt_jet3_DATA->Scale(1/pt_jet3_DATA->Integral());
  // pt_jet3_DATA->GetYaxis()->SetRangeUser(0,0.1);
  pt_jet3_DATA->Draw();
  pt_jet3_MC->SetMarkerStyle(22);
  pt_jet3_MC->SetMarkerSize(0.5);
  pt_jet3_MC->SetMarkerColor(2);
  pt_jet3_MC->SetLineColor(2);
  //  pt_jet3_MC->Scale(1/pt_jet3_MC->Integral());
  pt_jet3_MC->Draw("same");
  // TLatex *tex = new TLatex();
  // tex->SetNDC();
  // tex->SetTextSize(0.045); 
  // tex->DrawLatex(0.44,0.81,"MC ["+GenName+"]");
  //  tex->DrawLatex(0.44,0.71,"DATA [2.11fb^{-1} (13TeV)]");
  a->cd(4);
  eta_jet1_DATA->SetMarkerStyle(20);
  eta_jet1_DATA->SetMarkerSize(0.5);
  eta_jet1_DATA->SetMarkerColor(1);
  eta_jet1_DATA->SetLineColor(1);
  eta_jet1_DATA->GetXaxis()->SetTitle("#eta");
  //  eta_jet1_DATA->Scale(1./eta_jet1_DATA->Integral());
  //  eta_jet1_DATA->GetYaxis()->SetRangeUser(0,0.04);
  eta_jet1_DATA->Draw();
  eta_jet1_MC->SetMarkerStyle(22);
  eta_jet1_MC->SetMarkerSize(0.5);
  eta_jet1_MC->SetMarkerColor(2);
  eta_jet1_MC->SetLineColor(2);
  //  eta_jet1_MC->Scale(1./eta_jet1_MC->Integral());
  eta_jet1_MC->Draw("same");
  a->cd(5);
  eta_jet2_DATA->SetMarkerStyle(20);
  eta_jet2_DATA->SetMarkerSize(0.5);
  eta_jet2_DATA->SetMarkerColor(1);
  eta_jet2_DATA->SetLineColor(1);
  eta_jet2_DATA->GetXaxis()->SetTitle("#eta");
  //  eta_jet2_DATA->Scale(1./eta_jet2_DATA->Integral());
  //  eta_jet2_DATA->GetYaxis()->SetRangeUser(0,0.04);
  eta_jet2_DATA->Draw();
  eta_jet2_MC->SetMarkerStyle(22);
  eta_jet2_MC->SetMarkerSize(0.5);
  eta_jet2_MC->SetMarkerColor(2);
  eta_jet2_MC->SetLineColor(2);
  //  eta_jet2_MC->Scale(1./eta_jet2_MC->Integral());
  eta_jet2_MC->Draw("same");
 
  a->cd(6);
  // eta_jet3_DATA->SetMarkerStyle(20);
  // eta_jet3_DATA->SetMarkerSize(0.5);
  // eta_jet3_DATA->SetMarkerColor(1);
  // eta_jet3_DATA->SetLineColor(1);
  // eta_jet3_DATA->GetXaxis()->SetTitle("#eta");
  // //  eta_jet3_DATA->Scale(1./eta_jet3_DATA->Integral());
  // //  eta_jet3_DATA->GetYaxis()->SetRangeUser(0,0.04);
  // eta_jet3_DATA->Draw();
  // eta_jet3_MC->SetMarkerStyle(22);
  // eta_jet3_MC->SetMarkerSize(0.5);
  // eta_jet3_MC->SetMarkerColor(2);
  // eta_jet3_MC->SetLineColor(2);  
  // //  eta_jet3_MC->Scale(1./eta_jet3_MC->Integral());
  // eta_jet3_MC->Draw("same");

  Njets_DATA->SetMarkerStyle(20);
  Njets_DATA->SetMarkerSize(0.5);
  Njets_DATA->SetMarkerColor(1);
  Njets_DATA->SetLineColor(1);
  Njets_DATA->GetXaxis()->SetTitle("Number of jets");
  //  Njets_DATA->Scale(1./Njets_DATA->Integral());
  //  Njets_DATA->GetYaxis()->SetRangeUser(0,0.04);
  Njets_DATA->Draw();
  Njets_MC->SetMarkerStyle(22);
  Njets_MC->SetMarkerSize(0.5);
  Njets_MC->SetMarkerColor(2);
  Njets_MC->SetLineColor(2);  
  //  Njets_MC->Scale(1./Njets_MC->Integral());
  Njets_MC->Draw("same");

  leg1->Draw();

 
  a->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_jets.pdf");
  
  TCanvas* b = new TCanvas();
  b->Divide(3,2);
  b->cd(1);
  N_PV_DATA->SetMarkerStyle(20);
  N_PV_DATA->SetMarkerSize(0.5);
  N_PV_DATA->SetMarkerColor(1);
  N_PV_DATA->SetLineColor(1);
  //  N_PV_DATA->GetXaxis()->SetTitle("N PVtx");
  //  N_PV_DATA->Scale(1./N_PV_DATA->Integral());
  //  N_PV_DATA->GetYaxis()->SetRangeUser(0,0.15);
  N_PV_DATA->Draw();
  N_PV_MC->SetMarkerStyle(22);
  N_PV_MC->SetMarkerSize(0.5);
  N_PV_MC->SetMarkerColor(2);
  N_PV_MC->SetLineColor(2);
  //  N_PV_MC->Scale(1./N_PV_MC->Integral());
  N_PV_MC->Draw("same");
  TLegend *leg2;
  leg2 = new TLegend(0.63,0.68,0.9,0.86,"","brNDC");//x+0.1
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.045);
  leg2->SetFillColor(10);
  leg2->SetLineColor(1);
  leg2->SetTextFont(42);
  leg2->AddEntry(Njets_DATA,DATAtitle,"lp");
  leg2->AddEntry(Njets_MC,MCtitle,"lp");
  leg2->Draw();
  b->cd(2);
  // nPU_DATA->SetMarkerStyle(20);
  // nPU_DATA->SetMarkerSize(0.5);
  // nPU_DATA->SetMarkerColor(1);
  // nPU_DATA->SetLineColor(1);
  // //  nPU_DATA->GetXaxis()->SetTitle("N PVtx");
  // //  nPU_DATA->Scale(1./nPU_DATA->Integral());
  // //  nPU_DATA->GetYaxis()->SetRangeUser(0,0.15);
  // nPU_DATA->Draw();
  // nPU_MC->SetMarkerStyle(22);
  // nPU_MC->SetMarkerSize(0.5);
  // nPU_MC->SetMarkerColor(2);
  // nPU_MC->SetLineColor(2);
  // //  nPU_MC->Scale(1./nPU_MC->Integral());
  // nPU_MC->Draw("same");
  Eta_assym_top_DATA->SetMarkerStyle(20);
  Eta_assym_top_DATA->SetMarkerSize(0.5);
  Eta_assym_top_DATA->SetMarkerColor(1);
  Eta_assym_top_DATA->SetLineColor(1);
  Eta_assym_top_DATA->SetTitle("Asymmetry (eta)");
  //  Eta_assym_top_DATA->GetYaxis()->SetTitle("(#eta^{probe}_{negative} - #eta^{probe}_{positive})/#eta^{probe}_{positive}");
  Eta_assym_top_DATA->GetYaxis()->SetTitle("(N^{+}-N^{-})/(N^{+}+N^{-})");
  Eta_assym_top_DATA->GetYaxis()->SetRangeUser(-0.35,0.35);
  Eta_assym_top_DATA->Draw();
  Eta_assym_top_MC->SetMarkerStyle(22);
  Eta_assym_top_MC->SetMarkerSize(0.5);
  Eta_assym_top_MC->SetMarkerColor(2);
  Eta_assym_top_MC->SetLineColor(2);
  Eta_assym_top_MC->Draw("same");
  leg2->Draw();
  b->cd(3);
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
  // b->cd(3);
  // pt_rel_DATA->SetMarkerStyle(20);
  // pt_rel_DATA->SetMarkerSize(0.5);
  // pt_rel_DATA->SetMarkerColor(1);
  // pt_rel_DATA->SetLineColor(1);
  // //  pt_rel_DATA->Scale(1./pt_rel_DATA->Integral());
  // //  pt_rel_DATA->GetYaxis()->SetRangeUser(0,0.15);
  // pt_rel_DATA->Draw();
  // pt_rel_MC->SetMarkerStyle(22);
  // pt_rel_MC->SetMarkerSize(0.5);
  // pt_rel_MC->SetMarkerColor(2);
  // pt_rel_MC->SetLineColor(2);
  // //  pt_rel_MC->Scale(1./pt_rel_MC->Integral());
  // pt_rel_MC->Draw("same");
  //
  b->cd(4);
  gPad->SetLogy();
  asym_DATA->SetMarkerStyle(20);
  asym_DATA->SetMarkerSize(0.5);
  asym_DATA->SetMarkerColor(1);
  asym_DATA->SetLineColor(1);
  //  asym_DATA->Scale(1./asym_DATA->Integral());
  //  asym_DATA->GetYaxis()->SetRangeUser(0,0.1);
  asym_DATA->Draw();
  asym_MC->SetMarkerStyle(22);
  asym_MC->SetMarkerSize(0.5);
  asym_MC->SetMarkerColor(2);
  asym_MC->SetLineColor(2);
  //  asym_MC->Scale(1./asym_MC->Integral());
  asym_MC->Draw("same");
  leg1->Draw();
  b->cd(5);
  r_dijet_DATA->SetMarkerStyle(20);
  r_dijet_DATA->SetMarkerSize(0.5);
  r_dijet_DATA->SetMarkerColor(1);
  r_dijet_DATA->SetLineColor(1);
  //  r_dijet_DATA->Scale(1./r_dijet_DATA->Integral());
  //  r_dijet_DATA->GetYaxis()->SetRangeUser(0,0.04);
  r_dijet_DATA->Draw();
  r_dijet_MC->SetMarkerStyle(22);
  r_dijet_MC->SetMarkerSize(0.5);
  r_dijet_MC->SetMarkerColor(2);
  r_dijet_MC->SetLineColor(2);
  //  r_dijet_MC->Scale(1./r_dijet_MC->Integral());
  r_dijet_MC->Draw("same");
  leg1->Draw();
  b->cd(6);
  r_mpf_DATA->SetMarkerStyle(20);
  r_mpf_DATA->SetMarkerSize(0.5);
  r_mpf_DATA->SetMarkerColor(1);
  r_mpf_DATA->SetLineColor(1);
  //  r_mpf_DATA->Scale(1./r_mpf_DATA->Integral());
  //  r_mpf_DATA->GetYaxis()->SetRangeUser(0,0.05);
  r_mpf_DATA->Draw();
  r_mpf_MC->SetMarkerStyle(22);
  r_mpf_MC->SetMarkerSize(0.5);
  r_mpf_MC->SetMarkerColor(2);
  r_mpf_MC->SetLineColor(2);
  //  r_mpf_MC->Scale(1./r_mpf_MC->Integral());
  r_mpf_MC->Draw("same");
  // TLatex *tex2 = new TLatex();
  // tex2->SetNDC();
  // tex2->SetTextSize(0.045); 
  // tex2->DrawLatex(0.44,0.81,"MC ["+GenName+"]");
  // tex2->DrawLatex(0.44,0.71,"DATA [2.11fb^{-1} (13TeV)]");
  b->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_dijet.pdf");
  
 //  TCanvas* c = new TCanvas();
 //  c->Divide(2,2);
 //  c->cd(1);
 //  ptjet3_vs_alpha_DATA->SetTitle(DATAtitle);
 //  ptjet3_vs_alpha_DATA->GetXaxis()->SetTitle("#alpha");
 //  ptjet3_vs_alpha_DATA->GetYaxis()->SetTitle("p_{T} jet3, GeV");
 //  ptjet3_vs_alpha_DATA->GetXaxis()->SetTitleSize(0.05);
 //  ptjet3_vs_alpha_DATA->GetYaxis()->SetTitleSize(0.05);
 //  ptjet3_vs_alpha_DATA->Draw("colz");
 //  c->cd(2);
 //  ptjet3_vs_alpha_MC->SetTitle(DATAtitle);
 //  ptjet3_vs_alpha_MC->GetXaxis()->SetTitle("#alpha");
 //  ptjet3_vs_alpha_MC->GetYaxis()->SetTitle("p_{T} jet3, GeV");
 //  ptjet3_vs_alpha_MC->GetXaxis()->SetTitleSize(0.05);
 //  ptjet3_vs_alpha_MC->GetYaxis()->SetTitleSize(0.05);
 //  ptjet3_vs_alpha_MC->SetTitle(MCtitle);
 //  ptjet3_vs_alpha_MC->Draw("colz");
 //  c->cd(3);
 //  pt_ave_vs_alpha_DATA->SetTitle(DATAtitle);
 //  pt_ave_vs_alpha_DATA->GetXaxis()->SetTitle("#alpha");
 //  pt_ave_vs_alpha_DATA->GetYaxis()->SetTitle("p^{ave}_{T}, GeV");
 //  pt_ave_vs_alpha_DATA->GetXaxis()->SetTitleSize(0.05);
 //  pt_ave_vs_alpha_DATA->GetYaxis()->SetTitleSize(0.05);
 //  pt_ave_vs_alpha_DATA->Draw("colz");
 //  c->cd(4);
 //  pt_ave_vs_alpha_MC->GetXaxis()->SetTitle("#alpha");
 //  pt_ave_vs_alpha_MC->GetYaxis()->SetTitle("p^{ave}_{T}, GeV");
 //  pt_ave_vs_alpha_MC->GetXaxis()->SetTitleSize(0.05);
 //  pt_ave_vs_alpha_MC->GetYaxis()->SetTitleSize(0.05);
 //  pt_ave_vs_alpha_MC->SetTitle(MCtitle);
 //  pt_ave_vs_alpha_MC->Draw("colz");
 //  c->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_2Dalpha.pdf");



 //  // Plot trigger efficiency plots
 //  //  TString     sTriggerThreshold[] = {"40","60", "80", "140", "200", "260", "320", "400","500"};
 //  TString     sTriggerThreshold[] = {"40","60", "200", "260", "320", "400","500"};
 //  Float_t scale[] = {4.66251,53.6801,1.43334,6.2315,3.5902,1.2033}; //normalisation: obtained from const fit of the middle of the trigger efficiency curves
 //  Float_t     fTriggerThreshold[10] = {56,78,100,168,232,300,366,453,562,600};
 //  TH1F *histo_DATA_A[9]; 
 //  TH1F *histo_DATA_B[9]; 
 //  TH1F *histo_DATA_eff[9]; 
 //  TCanvas* d = new TCanvas();
 //  TLegend *leg3;
 //  leg3 = new TLegend(0.13,0.58,0.5,0.88,"","brNDC");//x+0.1
 //  leg3->SetBorderSize(0);
 //  leg3->SetTextSize(0.045);
 //  leg3->SetFillColor(10);
 //  leg3->SetLineColor(1);
 //  leg3->SetTextFont(42);
 //  TString trigName = "HLT_DiPFJetAve";
 //  d->Divide(3,2);
 //  for(int i=0; i < 6; i++){
 //    d->cd(i+1);
 //    histo_DATA_B[i] = (TH1F*)datafile-> Get(dirName+"/pt_ave_hltDiPFJetAve"+sTriggerThreshold[i]);
 //    histo_DATA_A[i] = (TH1F*)datafile-> Get(dirName+"/pt_ave_hltDiPFJetAve"+sTriggerThreshold[i+1]);
 //    histo_DATA_eff[i] = (TH1F*)histo_DATA_A[i]->Clone();
 //    histo_DATA_eff[i]->Divide(histo_DATA_B[i]);
 //    //    TF1 *fa1 = new TF1("fa1","pol0",300,400); 
 //    //    TF1 *fa1 = new TF1("fa1","pol0",400,500); 
 //    //    TF1 *fa1 = new TF1("fa1","pol0",500,600); 
 //    // TF1 *fa1 = new TF1("fa1","pol0",550,600); 
 //    // histo_DATA_eff[i]->Fit(fa1,"R");
 //    //    histo_DATA_eff[i]->Scale(1./scale[i]);
 //    histo_DATA_eff[i]->SetMarkerColor(kRed-1+i);
 //    histo_DATA_eff[i]->SetLineColor(kRed-1+i);
 //    histo_DATA_eff[i]->SetMarkerStyle(20);
 //    //    TF1 *fstep = new TF1("fstep","0.5*(1+TMath::Erf((x-[0])/(TMath::Sqrt(2)*[1])))",50,200); 
 //    //    TF1 *fstep = new TF1("fstep","TMath::Erf(x-[0])+[1]",50,200); 
 //    //    TF1 *fstep = new TF1("fstep","TMath::Erf(x)",0,600); 
 //    //  histo_DATA_eff[i]->Fit(fstep,"R");
 //    histo_DATA_eff[i]->GetYaxis()->SetRangeUser(-0.1,3.);
 //    histo_DATA_eff[i]->Draw();
 //  }
 //  d->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_TriggerEff.pdf");
 //  TCanvas* e = new TCanvas();
 //  for(int i=0; i < 6; i++){
 //    if(i==0){
 //      histo_DATA_eff[i]->GetYaxis()->SetRangeUser(-0.1,3.);
 //      histo_DATA_eff[i]->SetTitle("");
 //      histo_DATA_eff[i]->Draw();
 //    }
 //    else{
 //      histo_DATA_eff[i]->Draw("same");
 //    }
 //    leg3->AddEntry(histo_DATA_eff[i],trigName+sTriggerThreshold[i+1],"lp");
 //  }
 //  leg3->Draw();
 //  e->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_TriggerEffOnePlot.pdf");


 TCanvas* f = new TCanvas();
 TH2F *mpf_vs_etaProbe_DATA = (TH2F*)datafile->Get(dirName+"/mpf_vs_etaProbe");
 TH2F *mpf_vs_etaProbe_MC = (TH2F*)MCfile->Get(dirName+"/mpf_vs_etaProbe");
 TH2F *r_rel_vs_etaProbe_DATA = (TH2F*)datafile->Get(dirName+"/r_rel_vs_etaProbe");
 TH2F *r_rel_vs_etaProbe_MC = (TH2F*)MCfile->Get(dirName+"/r_rel_vs_etaProbe");
 f->Divide(2,2);
 f->cd(1);
 r_rel_vs_etaProbe_DATA->SetTitle(DATAtitle);
 r_rel_vs_etaProbe_DATA->Draw("colz");
 f->cd(2);
 r_rel_vs_etaProbe_MC->SetTitle(MCtitle);
 r_rel_vs_etaProbe_MC->Draw("colz");
 f->cd(3);
 mpf_vs_etaProbe_DATA->SetTitle(DATAtitle);
 mpf_vs_etaProbe_DATA->Draw("colz");
 f->cd(4);
 mpf_vs_etaProbe_MC->SetTitle(MCtitle);
 mpf_vs_etaProbe_MC->Draw("colz");
 f->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_ResponsesVsEta.pdf");

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
 pt_ave_vs_etaProbe_DATA->SetTitle(DATAtitle);
 pt_ave_vs_etaProbe_DATA->Draw("colz");
for(int i=0;i<2*n_eta;i++)
  lineEta[i]->Draw("same");
for(int i=0;i<n_pt;i++)
   linePt[i]->Draw("same");
 g->cd(2);
 pt_ave_vs_etaProbe_MC->SetTitle(MCtitle);
 pt_ave_vs_etaProbe_MC->Draw("colz");
for(int i=0;i<2*n_eta;i++)
  lineEta[i]->Draw("same");
for(int i=0;i<n_pt;i++)
   linePt[i]->Draw("same");

  g->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_Pt_aveVsEtaProbe.pdf");

 TCanvas* h1 = new TCanvas();
 TH2F *Rrel_vs_assym_DATA = (TH2F*)datafile->Get(dirName+"/Rrel_vs_assym");
 TH2F *Rrel_vs_assym_MC = (TH2F*)MCfile->Get(dirName+"/Rrel_vs_assym");

 h1->Divide(2,1);
 h1->cd(1);
 Rrel_vs_assym_DATA->SetTitle(DATAtitle);
 Rrel_vs_assym_DATA->Draw("colz");

 h1->cd(2);
 Rrel_vs_assym_MC->SetTitle(MCtitle);
 Rrel_vs_assym_MC->Draw("colz");
 h1->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_Rrel_vs_assym.pdf");

 TCanvas* h2 = new TCanvas();
 TH2F *Rmpf_vs_assym_DATA = (TH2F*)datafile->Get(dirName+"/Rmpf_vs_assym");
 TH2F *Rmpf_vs_assym_MC = (TH2F*)MCfile->Get(dirName+"/Rmpf_vs_assym");

 h2->Divide(2,1);
 h2->cd(1);
 Rmpf_vs_assym_DATA->SetTitle(DATAtitle);
 Rmpf_vs_assym_DATA->Draw("colz");

 h2->cd(2);
 Rmpf_vs_assym_MC->SetTitle(MCtitle);
 Rmpf_vs_assym_MC->Draw("colz");
 h2->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_Rmpf_vs_assym.pdf");


 TCanvas* h3 = new TCanvas();
 TH2F *alpha_vs_assym_DATA = (TH2F*)datafile->Get(dirName+"/alpha_vs_assym");
 TH2F *alpha_vs_assym_MC = (TH2F*)MCfile->Get(dirName+"/alpha_vs_assym");

 h3->Divide(2,1);
 h3->cd(1);
 alpha_vs_assym_DATA->SetTitle(DATAtitle);
 alpha_vs_assym_DATA->Draw("colz");

 h3->cd(2);
 alpha_vs_assym_MC->SetTitle(MCtitle);
 alpha_vs_assym_MC->Draw("colz");
 h3->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_alpha_vs_assym.pdf");


 TCanvas* h4 = new TCanvas();
 TH2F *Njets_alpha_DATA = (TH2F*)datafile->Get(dirName+"/Njets_vs_alpha");
 TH2F *Njets_alpha_MC = (TH2F*)MCfile->Get(dirName+"/Njets_vs_alpha");

 h4->Divide(2,1);
 h4->cd(1);
 Njets_alpha_DATA->SetTitle(DATAtitle);
 Njets_alpha_DATA->Draw("colz");

 h4->cd(2);
 Njets_alpha_MC->SetTitle(MCtitle);
 Njets_alpha_MC->Draw("colz");
 h4->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_N_jets_alpha.pdf");

 TCanvas* h5 = new TCanvas();
 TH2F *alphaSum_vs_alpha3_DATA = (TH2F*)datafile->Get(dirName+"/alphaSum_vs_alpha3");
 TH2F *alphaSum_vs_alpha3_MC = (TH2F*)MCfile->Get(dirName+"/alphaSum_vs_alpha3");

 h5->Divide(2,1);
 h5->cd(1);
 alphaSum_vs_alpha3_DATA->SetTitle(DATAtitle);
 alphaSum_vs_alpha3_DATA->Draw("colz");

 h5->cd(2);
 alphaSum_vs_alpha3_MC->SetTitle(MCtitle);
 alphaSum_vs_alpha3_MC->Draw("colz");
 h5->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_alphaSum_vs_alpha3.pdf");


 TCanvas* h6 = new TCanvas();
 TH2F *Pt_vs_JetN_DATA = (TH2F*)datafile->Get(dirName+"/Pt_vs_JetN");
 TH2F *Pt_vs_JetN_MC = (TH2F*)MCfile->Get(dirName+"/Pt_vs_JetN");

 h6->Divide(2,1);
 h6->cd(1);
 Pt_vs_JetN_DATA->SetTitle(DATAtitle);
 Pt_vs_JetN_DATA->Draw("colz");

 h6->cd(2);
 Pt_vs_JetN_MC->SetTitle(MCtitle);
 Pt_vs_JetN_MC->Draw("colz");
 h6->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_Pt_vs_JetN.pdf");

 TCanvas* h7 = new TCanvas();
 TH2F *eta_vs_JetN_DATA = (TH2F*)datafile->Get(dirName+"/eta_vs_JetN");
 TH2F *eta_vs_JetN_MC = (TH2F*)MCfile->Get(dirName+"/eta_vs_JetN");

 h7->Divide(2,1);
 h7->cd(1);
 eta_vs_JetN_DATA->SetTitle(DATAtitle);
 eta_vs_JetN_DATA->Draw("colz");

 h7->cd(2);
 eta_vs_JetN_MC->SetTitle(MCtitle);
 eta_vs_JetN_MC->Draw("colz");
 h7->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_eta_vs_JetN.pdf");

 // TCanvas* h8 = new TCanvas();
 // //TH1F *diffPt1_vs_alpha3_DATA = (TH1F*)datafile->Get(dirName+"/diffPt1_vs_alpha3");
 // TH2F *diffPt1_vs_alpha3_MC = (TH2F*)MCfile->Get(dirName+"/diffPt1_vs_alpha3");
 // TH2F *diffPt2_vs_alpha3_MC = (TH2F*)MCfile->Get(dirName+"/diffPt2_vs_alpha3");

 // h8->Divide(2,1);
 // h8->cd(1);
 // diffPt1_vs_alpha3_MC->SetTitle("1st jet");
 // diffPt1_vs_alpha3_MC->Draw("colz");

 // h8->cd(2);
 // diffPt2_vs_alpha3_MC->SetTitle("2nd jet");
 // diffPt2_vs_alpha3_MC->Draw("colz");
 // h8->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_diffPt_vs_alpha3.pdf");


 // TCanvas* h9 = new TCanvas();
 // //TH1F *diffPt1_vs_alpha3_DATA = (TH1F*)datafile->Get(dirName+"/diffPt1_vs_alpha3");
 // TH2F *diffPt1_vs_alphaSum_MC = (TH2F*)MCfile->Get(dirName+"/diffPt1_vs_alphaSum");
 // TH2F *diffPt2_vs_alphaSum_MC = (TH2F*)MCfile->Get(dirName+"/diffPt2_vs_alphaSum");

 // h9->Divide(2,1);
 // h9->cd(1);
 // diffPt1_vs_alphaSum_MC->SetTitle("1st jet");
 // diffPt1_vs_alphaSum_MC->Draw("colz");

 // h9->cd(2);
 // diffPt2_vs_alphaSum_MC->SetTitle("2nd jet");
 // diffPt2_vs_alphaSum_MC->Draw("colz");
 // h9->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_diffPt_vs_alphaSum.pdf");


 TCanvas* h10 = new TCanvas();
 //TH1F *diffPt1_vs_alpha3_DATA = (TH1F*)datafile->Get(dirName+"/diffPt1_vs_alpha3");
 TH2F *diffPtGen_vs_alpha3_MC = (TH2F*)MCfile->Get(dirName+"/diffPtGen_vs_alpha3");
 TH2F *diffPtGen_vs_alphaSum_MC = (TH2F*)MCfile->Get(dirName+"/diffPtGen_vs_alphaSum");

 h10->Divide(2,1);
 h10->cd(1);
 diffPtGen_vs_alpha3_MC->SetTitle("#alpha_{3}");
 diffPtGen_vs_alpha3_MC->Draw("colz");

 h10->cd(2);
 diffPtGen_vs_alphaSum_MC->SetTitle("#alpha_{Sum}");
 diffPtGen_vs_alphaSum_MC->Draw("colz");
 h10->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_diffPtGEN_vs_alpha.pdf");


 
 
 TH2F *diffPtGen_vs_alpha3_DATA = (TH2F*)datafile->Get(dirName+"/diffPtGen_vs_alpha3");
 if(DATAtitle=="Pythia8" ||  DATAtitle== "Herwigpp"){
   TCanvas* h11 = new TCanvas();
   TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
   pad1->SetBottomMargin(0); // Upper and lower plot are joined
   //   pad1->SetGridx();         // Vertical grid
   pad1->Draw();             // Draw the upper pad: pad1
   pad1->cd();               // pad1 becomes the current pad
   // h11->Divide(2,1);
   // h11->cd(1);
   TProfile *diffPtGen_vs_alpha3_MC_Profile = diffPtGen_vs_alpha3_MC->ProfileX();
   diffPtGen_vs_alpha3_MC_Profile->GetYaxis()->SetTitle("<(pT_{1}-pT_{2})_{gen}/pT^{gen}_{ave}>");
   diffPtGen_vs_alpha3_MC_Profile->SetTitle("<(pT_{1}-pT_{2})_{gen}/pT^{gen}_{ave}>");
   // diffPtGen_vs_alpha3_MC_Profile->Print();
   diffPtGen_vs_alpha3_MC_Profile->SetName("diffPtGen_vs_alpha3_pfx_MC");
   TF1 *fit_lin = new TF1("diffPtGen_vs_alpha3_MC_1_fit","pol1",0.19,0.36); //AK4
   // TF1 *fit_lin = new TF1("diffPtGen_vs_alpha3_MC_1_fit","pol1",0.19,0.41); //AK8Puppi
   // TF1 *fit_lin = new TF1("diffPtGen_vs_alpha3_MC_1_fit","pol1",0.19,0.46); //AK8CHS
   diffPtGen_vs_alpha3_MC_Profile->Fit("diffPtGen_vs_alpha3_MC_1_fit","SR");
   diffPtGen_vs_alpha3_MC_Profile->SetMarkerStyle(20);
   diffPtGen_vs_alpha3_MC_Profile->SetMarkerColor(kOrange+7);
   diffPtGen_vs_alpha3_MC_Profile->GetXaxis()->SetRangeUser(0,0.5);
   TString chi2 = "#chi^{2}/n.d.f = ";
   chi2 += trunc(fit_lin->GetChisquare());
   chi2 +="/";
   chi2 +=trunc(fit_lin->GetNDF());
   TLatex *tex3 = new TLatex();
   tex3->SetNDC();
   tex3->SetTextSize(0.035); 
   TString p0 = Form("p0 = %g", fit_lin->GetParameter(0));
   p0 += Form(" +/- %g",fit_lin->GetParError(0));
   TString p1 = Form("p1 = %g", fit_lin->GetParameter(1));
   p1 += Form(" +/- %g",fit_lin->GetParError(1));
   diffPtGen_vs_alpha3_MC_Profile->Draw();
   tex3->SetTextColor(kOrange+7);
   
   TProfile *diffPtGen_vs_alpha3_DATA_Profile = diffPtGen_vs_alpha3_DATA->ProfileX();
   diffPtGen_vs_alpha3_DATA_Profile->SetName("diffPtGen_vs_alpha3_pfx_DATA");
   diffPtGen_vs_alpha3_DATA_Profile->GetYaxis()->SetTitle("<(pT_{1}-pT_{2})_{gen}/pT^{gen}_{ave}>");
   diffPtGen_vs_alpha3_DATA_Profile->SetTitle("");
   TH1F *ratio = (TH1F*)diffPtGen_vs_alpha3_DATA_Profile->Clone("ratio_DATA_MC");
   //   diffPtGen_vs_alpha3_DATA_Profile->Print();
   TF1 *fit_lin_DATA = new TF1("diffPtGen_vs_alpha3_DATA_1_fit","pol1",0.19,0.36); //AK4
   diffPtGen_vs_alpha3_DATA_Profile->Fit("diffPtGen_vs_alpha3_DATA_1_fit","SR");
   diffPtGen_vs_alpha3_DATA_Profile->SetMarkerStyle(21);
   diffPtGen_vs_alpha3_DATA_Profile->SetMarkerColor(kGreen+3);
   diffPtGen_vs_alpha3_DATA_Profile->GetXaxis()->SetRangeUser(0,0.5);
   TString chi2_DATA = "#chi^{2}/n.d.f = ";
   chi2_DATA += trunc(fit_lin_DATA->GetChisquare());
   chi2_DATA +="/";
   chi2_DATA +=trunc(fit_lin_DATA->GetNDF());
   TLatex *tex3_DATA = new TLatex();
   tex3_DATA->SetNDC();
   tex3_DATA->SetTextSize(0.035); 
   TString p0_DATA = Form("p0 = %g", fit_lin_DATA->GetParameter(0));
   p0_DATA += Form(" +/- %g",fit_lin_DATA->GetParError(0));
   TString p1_DATA = Form("p1 = %g", fit_lin_DATA->GetParameter(1));
   p1_DATA += Form(" +/- %g",fit_lin_DATA->GetParError(1));
   //  
   tex3_DATA->SetTextColor(kGreen+3);
   tex3_DATA->DrawLatex(0.55,0.45,DATAtitle);
   tex3_DATA->DrawLatex(0.55,0.4,chi2_DATA);
   tex3_DATA->DrawLatex(0.55,0.35,p0_DATA);
   tex3_DATA->DrawLatex(0.55,0.3,p1_DATA);
   //   diffPtGen_vs_alpha3_MC_Profile->Draw("hist SAME");
   diffPtGen_vs_alpha3_DATA_Profile->Draw("SAME");
   diffPtGen_vs_alpha3_DATA_Profile->GetYaxis()->SetTitleSize(15);
   diffPtGen_vs_alpha3_DATA_Profile->GetYaxis()->SetTitleFont(43);
   diffPtGen_vs_alpha3_DATA_Profile->GetYaxis()->SetTitleOffset(1.2);
   diffPtGen_vs_alpha3_DATA_Profile->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   diffPtGen_vs_alpha3_DATA_Profile->GetYaxis()->SetLabelSize(15);
   diffPtGen_vs_alpha3_MC_Profile->Draw("SAME");
   //   diffPtGen_vs_alpha3_DATA_Profile->GetYaxis()->SetLabelSize(0.);
   // TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
   // axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   // axis->SetLabelSize(15);
   // axis->Draw();

   tex3->DrawLatex(0.2,0.85,MCtitle);
   tex3->DrawLatex(0.2,0.8,chi2);
   tex3->DrawLatex(0.2,0.75,p0);
   tex3->DrawLatex(0.2,0.7,p1);
   h11->cd();          // Go back to the main canvas before defining pad2
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.2);
   //   pad2->SetGridx(); // vertical grid
   pad2->Draw();
   pad2->cd();       // pad2 becomes the current pad
   ratio->GetXaxis()->SetRangeUser(0,0.5);
   ratio->Divide(diffPtGen_vs_alpha3_MC_Profile);
   ratio->SetMarkerStyle(21);
   ratio->SetFillColor(14);
   ratio->SetLineColor(14);
   ratio->Draw("hist");
   TString rtitle = DATAtitle+"/"+MCtitle;
   ratio->GetYaxis()->SetRangeUser(0.85,1.15);
   ratio->GetYaxis()->SetTitle(rtitle);
   ratio->GetYaxis()->SetNdivisions(505);
   ratio->GetYaxis()->SetTitleSize(15);
   ratio->GetYaxis()->SetTitleFont(43);
   ratio->GetYaxis()->SetTitleOffset(1.2);
   ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   ratio->GetYaxis()->SetLabelSize(15);

   // X axis ratio plot settings
   ratio->GetXaxis()->SetTitleSize(15);
   ratio->GetXaxis()->SetTitleFont(43);
   ratio->GetXaxis()->SetTitleOffset(4.);
   ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   ratio->GetXaxis()->SetLabelSize(15);
   h11->Print(path+"plots/Control_Plots_"+dirName+"_AsData_"+DATAtitle+"_AsMC_"+MCtitle+"_diffPtGEN_vs_alpha_Profile.pdf");
 }


 

// TF1 *fit_lin = new TF1("diffPtGen_vs_alpha3_MC_1_fit","pol1",0,0.35);
 // TF1 *fit_lin = new TF1("diffPtGen_vs_alpha3_MC_1_fit","pol1",0.04,0.46);
 
 // h11->cd(2);
 // diffPtGen_vs_alpha3_MC_Profile->Draw("E");
// TProfile *diffPtGen_vs_alphaSum_MC_Profile = diffPtGen_vs_alphaSum_MC->ProfileX();
//  diffPtGen_vs_alphaSum_MC_Profile->GetYaxis()->SetTitle("<(pT_{1}-pT_{2})_{gen}/pT^{gen}_{ave}>");
//  diffPtGen_vs_alphaSum_MC_Profile->SetTitle("<(pT_{1}-pT_{2})_{gen}/pT^{gen}_{ave}>");
//  TF1 *fit_lin_2 = new TF1("diffPtGen_vs_alphaSum_MC_1_fit","pol1",0,0.7);
//  diffPtGen_vs_alphaSum_MC_Profile->Fit("diffPtGen_vs_alphaSum_MC_1_fit","SR");
//  diffPtGen_vs_alphaSum_MC_Profile->SetMarkerStyle(20);
//  diffPtGen_vs_alphaSum_MC_Profile->GetXaxis()->SetRangeUser(0,0.9);
//  TString chi2_2 = "#chi^{2}/n.d.f = ";
//  chi2_2 += trunc(fit_lin_2->GetChisquare());
//  chi2_2 +="/";
//  chi2_2 +=trunc(fit_lin_2->GetNDF());
//  TString p0_2 = Form("p0 = %g", fit_lin_2->GetParameter(0));
//  p0_2 += Form(" +/- %g",fit_lin_2->GetParError(0));
//  TString p1_2 = Form("p1 = %g", fit_lin_2->GetParameter(1));
//  p1_2 += Form(" +/- %g",fit_lin_2->GetParError(1));
//  diffPtGen_vs_alphaSum_MC_Profile->Draw();
//  tex3->DrawLatex(0.2,0.85,chi2_2);
//  tex3->DrawLatex(0.2,0.8,p0_2);
//  tex3->DrawLatex(0.2,0.75,p1_2);

 // diffPtGen_vs_alphaSum_MC->FitSlicesY();
 // TH1F*  diffPtGen_vs_alphaSum_MC_1 = (TH1F*)gDirectory->Get("diffPtGen_vs_alphaSum_1");
 // diffPtGen_vs_alphaSum_MC_1->Draw();


 TCanvas* h12 = new TCanvas();
 TH1D *alpha3_MC = ptjet3_vs_alpha_MC->ProjectionX("pxMC");
 alpha3_MC->SetMarkerStyle(20);
 alpha3_MC->SetMarkerColor(kOrange+7);
 alpha3_MC->GetXaxis()->SetRangeUser(0,1.1);
 double integ_alpha3_MC = alpha3_MC->Integral();
 alpha3_MC->Scale(1./integ_alpha3_MC);
 // TH2F *ptjet3_vs_alpha_DATA = (TH2F*)datafile->Get(dirName+"/ptjet3_vs_alpha");
 TH1D *alpha3_DATA = ptjet3_vs_alpha_DATA->ProjectionX("pxDATA");
 alpha3_DATA->SetMarkerStyle(21);
 alpha3_DATA->SetMarkerColor(kGreen+3);
 alpha3_DATA->GetXaxis()->SetRangeUser(0,1.1);
 double integ_alpha3_DATA = alpha3_DATA->Integral();
 alpha3_DATA->Scale(1./integ_alpha3_DATA);
 alpha3_DATA->GetXaxis()->SetTitle("#alpha");
 alpha3_DATA->SetTitle("");
 alpha3_DATA->Draw();
 alpha3_MC->Draw("same");
 TLegend *leg12;
  leg12 = new TLegend(0.63,0.68,0.9,0.85,"","brNDC");//x+0.1
  leg12->SetBorderSize(0);
  leg12->SetTextSize(0.045);
  leg12->SetFillColor(10);
  leg12->SetLineColor(1);
  leg12->SetTextFont(42);
  leg12->AddEntry(alpha3_DATA,DATAtitle,"lp");
  leg12->AddEntry(alpha3_MC,MCtitle,"lp");
 if(Realdatafile){
   TH2F *ptjet3_vs_alpha_RDATA = (TH2F*)Realdatafile->Get(dirName+"/ptjet3_vs_alpha");
   TH1D *alpha3_RDATA = ptjet3_vs_alpha_RDATA->ProjectionX("pxRDATA");
   alpha3_RDATA->SetMarkerStyle(20);
   alpha3_RDATA->GetXaxis()->SetRangeUser(0,1.1);
   double integ_alpha3_RDATA = alpha3_RDATA->Integral();
   alpha3_RDATA->Scale(1./integ_alpha3_RDATA);
   alpha3_RDATA->Draw("same");
   leg12->AddEntry(alpha3_RDATA,"DATA","lp");
 }
  leg12->Draw();
  h12->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_alpha3.pdf");
 // TCanvas* h12 = new TCanvas();
 // TH2F *AbsdiffPtGen_vs_alpha3_MC = (TH2F*)MCfile->Get(dirName+"/AbsdiffPtGen_vs_alpha3");
 // TH2F *AbsdiffPtGen_vs_alphaSum_MC = (TH2F*)MCfile->Get(dirName+"/AbsdiffPtGen_vs_alphaSum");

 // h12->Divide(2,1);
 // h12->cd(1);
 // AbsdiffPtGen_vs_alpha3_MC->SetTitle("#alpha_{3}");
 // AbsdiffPtGen_vs_alpha3_MC->Draw("colz");

 // h12->cd(2);
 // AbsdiffPtGen_vs_alphaSum_MC->SetTitle("#alpha_{Sum}");
 // AbsdiffPtGen_vs_alphaSum_MC->Draw("colz");
 // h12->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_AbsdiffPtGEN_vs_alpha.pdf");

 // TCanvas* h13 = new TCanvas();
 // h13->Divide(2,1);
 // h13->cd(1);
 // AbsdiffPtGen_vs_alpha3_MC->FitSlicesY();
 // TH1F* AbsdiffPtGen_vs_alpha3_MC_1 = (TH1F*)gDirectory->Get("AbsdiffPtGen_vs_alpha3_1");
 // AbsdiffPtGen_vs_alpha3_MC_1->Draw();
 // h13->cd(2);
 // AbsdiffPtGen_vs_alphaSum_MC->FitSlicesY();
 // TH1F*  AbsdiffPtGen_vs_alphaSum_MC_1 = (TH1F*)gDirectory->Get("AbsdiffPtGen_vs_alphaSum_1");
 // AbsdiffPtGen_vs_alphaSum_MC_1->Draw();
 // h13->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_AbsdiffPtGEN_vs_alpha_Profile.pdf");
}
