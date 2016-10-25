// plot control plots of important parameters for L2Res determination

#include "header.h"


void Control_Plots(bool divide_by_lumi, TString path, TFile* datafile, TFile* MCfile, TString GenName="pythia8",TFile* Realdatafile=NULL){
 
  TString DATAtitle = "DATA";
  TString MCtitle = "MC";
  //  TString MCtitle = "Herwigpp";
  //  TString MCtitle = "Pythia8";
  //  TString DATAtitle = "Pythia8";

  // // TString DATAtitle = "Herwigpp";
  // // TString MCtitle = "Pythia8"; 
  TString dirName = "Selection";
  //  TString dirName = "JetMatching";
  //    TString dirName = "diJet";
  //  TString dirName = "noCuts";
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
  leg1->AddEntry(pt_jet1_DATA,DATAtitle,"lp");
  leg1->AddEntry(pt_jet1_MC,MCtitle,"lp");
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

  Njets_DATA->SetMarkerStyle(20);
  Njets_DATA->SetMarkerSize(0.5);
  Njets_DATA->SetMarkerColor(1);
  Njets_DATA->SetLineColor(1);
  Njets_DATA->GetXaxis()->SetTitle("Number of jets");
  Njets_DATA->Scale(1./Njets_DATA->Integral());
  Njets_DATA->GetYaxis()->SetRangeUser(0,0.1);
  Njets_DATA->Draw();
  Njets_MC->SetMarkerStyle(22);
  Njets_MC->SetMarkerSize(0.5);
  Njets_MC->SetMarkerColor(2);
  Njets_MC->SetLineColor(2);  
  Njets_MC->Scale(1./Njets_MC->Integral());
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
  N_PV_DATA->GetXaxis()->SetTitle("N good PV");
  N_PV_DATA->Scale(1./N_PV_DATA->Integral());
  N_PV_DATA->GetYaxis()->SetRangeUser(0,0.1);
  N_PV_DATA->Draw();
  N_PV_MC->SetMarkerStyle(22);
  N_PV_MC->SetMarkerSize(0.5);
  N_PV_MC->SetMarkerColor(2);
  N_PV_MC->SetLineColor(2);
  N_PV_MC->Scale(1./N_PV_MC->Integral());
  N_PV_MC->Draw("same");

  TLegend *leg2;
  leg2 = new TLegend(0.63,0.68,0.87,0.86,"","brNDC");//x+0.1
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.045);
  leg2->SetFillColor(10);
  leg2->SetLineColor(1);
  leg2->SetTextFont(42);
  leg2->AddEntry(Njets_DATA,DATAtitle,"lp");
  leg2->AddEntry(Njets_MC,MCtitle,"lp");
  leg2->Draw();
  b->cd(2);
 
  leg2->Draw();
  b->cd(3);
  pt_ave_DATA->SetMarkerStyle(20);
  pt_ave_DATA->SetMarkerSize(0.5);
  pt_ave_DATA->SetMarkerColor(1);
  pt_ave_DATA->SetLineColor(1);
  pt_ave_DATA->Scale(1./pt_ave_DATA->Integral());
  pt_ave_DATA->GetYaxis()->SetRangeUser(0,0.03);
  pt_ave_DATA->Draw();
  pt_ave_MC->SetMarkerStyle(22);
  pt_ave_MC->SetMarkerSize(0.5);
  pt_ave_MC->SetMarkerColor(2);
  pt_ave_MC->SetLineColor(2);
  pt_ave_MC->Scale(1./pt_ave_MC->Integral());
  pt_ave_MC->Draw("same");

  b->cd(4);
  gPad->SetLogy();
  asym_DATA->SetMarkerStyle(20);
  asym_DATA->SetMarkerSize(0.5);
  asym_DATA->SetMarkerColor(1);
  asym_DATA->SetLineColor(1);
  asym_DATA->Scale(1./asym_DATA->Integral());
  //  asym_DATA->GetYaxis()->SetRangeUser(0,0.1);
  asym_DATA->Draw();
  asym_MC->SetMarkerStyle(22);
  asym_MC->SetMarkerSize(0.5);
  asym_MC->SetMarkerColor(2);
  asym_MC->SetLineColor(2);
  asym_MC->Scale(1./asym_MC->Integral());
  asym_MC->Draw("same");
  leg1->Draw();
  b->cd(5);
  r_dijet_DATA->SetMarkerStyle(20);
  r_dijet_DATA->SetMarkerSize(0.5);
  r_dijet_DATA->SetMarkerColor(1);
  r_dijet_DATA->SetLineColor(1);
  r_dijet_DATA->Scale(1./r_dijet_DATA->Integral());
  r_dijet_DATA->GetYaxis()->SetRangeUser(0,0.03);
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
  r_mpf_DATA->GetYaxis()->SetRangeUser(0,0.04);
  r_mpf_DATA->Draw();
  r_mpf_MC->SetMarkerStyle(22);
  r_mpf_MC->SetMarkerSize(0.5);
  r_mpf_MC->SetMarkerColor(2);
  r_mpf_MC->SetLineColor(2);
  r_mpf_MC->Scale(1./r_mpf_MC->Integral());
  r_mpf_MC->Draw("same");
  b->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_dijet.pdf");
  



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
  TH2F *mpf_vs_etaProbe_DATA_norm = (TH2F*)mpf_vs_etaProbe_DATA->Clone("mpf_vs_etaProbe_DATA_norm");
  TH2F *mpf_vs_etaProbe_MC_norm = (TH2F*)mpf_vs_etaProbe_MC->Clone("mpf_vs_etaProbe_MC_norm");
  TH2F *r_rel_vs_etaProbe_DATA_norm = (TH2F*)r_rel_vs_etaProbe_DATA->Clone("r_rel_vs_etaProbe_DATA_norm");
  TH2F *r_rel_vs_etaProbe_MC_norm = (TH2F*)r_rel_vs_etaProbe_MC->Clone("r_rel_vs_etaProbe_MC_norm");




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
  /*TCanvas* f_norm = new TCanvas();
  f_norm->Divide(2,2);
  f_norm->cd(1);
  r_rel_vs_etaProbe_DATA_norm->SetTitle(DATAtitle+" normalised");				    
  r_rel_vs_etaProbe_DATA_norm->Draw("colz");
  f_norm->cd(2);	
  r_rel_vs_etaProbe_MC_norm->SetTitle(MCtitle+" normalised");					    
  r_rel_vs_etaProbe_MC_norm->Draw("colz");
  f_norm->cd(3);
  mpf_vs_etaProbe_DATA_norm->SetTitle(DATAtitle+" normalised");					    
  mpf_vs_etaProbe_DATA_norm->Draw("colz");						    
  f_norm->cd(4);
  mpf_vs_etaProbe_MC_norm->SetTitle(MCtitle+" normalised");					    
  mpf_vs_etaProbe_MC_norm->Draw("colz");
  f->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_ResponsesVsEta_norm.pdf");*/

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



  //If lumi_bins are activated: plot 2d-plots of jet-response vs inst lumi for all eta and different eta bins-->use AnalysisTree observables
  //For DATA only, of course
  if(divide_by_lumi){

    //Get values from TTree
    TTreeReader myReader_DATA("AnalysisTree", datafile);
    TTreeReaderValue<Float_t> probejet_eta_data(myReader_DATA, "probejet_eta");
    TTreeReaderValue<Float_t> alpha_data(myReader_DATA, "alpha");
    TTreeReaderValue<Float_t> rel_r_data(myReader_DATA, "rel_r");
    TTreeReaderValue<Float_t> mpf_r_data(myReader_DATA, "mpf_r");
    TTreeReaderValue<Float_t> inst_lumi_data(myReader_DATA, "instantaneous_lumi");
   
    //Set up 2d-histos and one lumi-hist
    TH1D* h_lumi = new TH1D("h_lumi",";instantaneous lumi;events",50,0,0.01);
    TH2D* h_rel_r_lumi_alleta = new TH2D("h_rel_r_lumi_alleta","0 < |#eta| < "+eta_range[n_eta-1]+";instantaneous lumi;p_{T} balance",10, 0, 0.01, 100, 0, 6);
    TH2D* h_rel_r_lumi_etabins[n_eta-1];
    TH2D* h_mpf_r_lumi_alleta = new TH2D("h_mpf_r_lumi_alleta","0 < |#eta| < "+eta_range[n_eta-1]+";instantaneous lumi;MPF response", 10, 0, 0.01, 100, 0, 6);
    TH2D* h_mpf_r_lumi_etabins[n_eta-1];
    for(int i=0; i<n_eta-1; i++){
      stringstream ss;
      ss << i;
      string str = ss.str();
      TString rel_r_name = "h_rel_r_lumi_etabins"+str;
      TString mpf_r_name = "h_mpf_r_lumi_etabins"+str;
      TString rel_r_title = eta_range[i]+" < |#eta| < "+eta_range[i+1]+";instantaneous lumi;p_{T} balance";
      TString mpf_r_title = eta_range[i]+" < |#eta| < "+eta_range[i+1]+";instantaneous lumi;MPF response";
      h_rel_r_lumi_etabins[i] = new TH2D(rel_r_name,rel_r_title, 10, 0, 0.01, 100,0,6);
      h_mpf_r_lumi_etabins[i] = new TH2D(mpf_r_name,mpf_r_title,10, 0, 0.01, 100,0,6);
    }

    //loop over all events to fill histograms
    while (myReader_DATA.Next()) {
      //regardless of eta values
      h_rel_r_lumi_alleta->Fill(*inst_lumi_data,*rel_r_data);
      h_mpf_r_lumi_alleta->Fill(*inst_lumi_data,*mpf_r_data);
      h_lumi->Fill(*inst_lumi_data);

      //in bins of |eta|
      for(int i=0; i<n_eta-1; i++){
	if(*probejet_eta_data>eta_bins[i] && *probejet_eta_data<eta_bins[i+1]){
	  h_rel_r_lumi_etabins[i]->Fill(*inst_lumi_data,*rel_r_data);
	  h_mpf_r_lumi_etabins[i]->Fill(*inst_lumi_data,*mpf_r_data);
	}
      }
    }

    //Draw plots
    TCanvas* c_lumi = new TCanvas("c_lumi","c_lumi",1);
    h_lumi->Draw();
    c_lumi->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_Instantaneous_Lumi.pdf");
    TCanvas* c_mpf_all = new TCanvas();
    h_mpf_r_lumi_alleta->Draw("colz");
    c_mpf_all->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_mpf_r_vs_lumi_eta_full.pdf");
    TCanvas* c_rel_all = new TCanvas();
    h_rel_r_lumi_alleta->Draw("colz");
    c_rel_all->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_rel_r_vs_lumi_eta_full.pdf");
    TCanvas *c_mpf_eta[n_eta-1], *c_rel_eta[n_eta-1];
    for(int i=0; i<n_eta-1; i++){
      TString s_eta_range = eta_range[i]+"_"+eta_range[i+1];
      c_mpf_eta[i] = new TCanvas();
      h_mpf_r_lumi_etabins[i]->Draw("colz");
      c_mpf_eta[i]->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_mpf_r_vs_lumi_eta_"+s_eta_range+".pdf");
      c_rel_eta[i] = new TCanvas();
      h_rel_r_lumi_etabins[i]->Draw("colz");
      c_rel_eta[i]->Print(path+"plots/Control_Plots_"+dirName+"_"+GenName+"_rel_r_vs_lumi_eta_"+s_eta_range+".pdf");
    }

  }//divide by lumi --> for 2d plots vs inst lumi



  datafile->Close();
  MCfile->Close();
}
