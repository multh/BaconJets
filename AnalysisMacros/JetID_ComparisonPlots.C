#include "header.h"

void JetID_ComparisonPlots(TString path_with, TString MCFile_name, TString datfile_name, TString Generator, TString Collection, TString runnr){
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  //pretty unflexible - adjust paths and filenames by hand! Structure is assumed to be kept the same in the future
  //only do sth if initial collection is WITH jetID and always compare to 'without JetID'

  //get nominal hists from both, with and without JetID
  //get files
  TString path_without = path_with;
  path_without.ReplaceAll(Collection+"_Loose_JetPFID",Collection);
  TFile* MC_without = new TFile(path_without+MCFile_name,"READ");
  TFile* DATA_without = new TFile(path_without+datfile_name,"READ");
  TFile* MC_with = new TFile(path+MCFile_name,"READ");
  TFile* DATA_with = new TFile(path+datfile_name,"READ");

  //get hists without ID
  TH2F* mpf_vs_eta_MC_without = (TH2F*)MC_without->Get("Selection/mpf_vs_etaProbe");
  TH2F* mpf_vs_eta_DATA_without = (TH2F*)DATA_without->Get("Selection/mpf_vs_etaProbe");
  TH2F* r_rel_vs_eta_MC_without = (TH2F*)MC_without->Get("Selection/r_rel_vs_etaProbe");
  TH2F* r_rel_vs_eta_DATA_without = (TH2F*)DATA_without->Get("Selection/r_rel_vs_etaProbe");
  mpf_vs_eta_MC_without->RebinX(3);
  mpf_vs_eta_DATA_without->RebinX(3);
  r_rel_vs_eta_MC_without->RebinX(3);
  r_rel_vs_eta_DATA_without->RebinX(3);

  //now with
  TH2F* mpf_vs_eta_MC_with = (TH2F*)MC_with->Get("Selection/mpf_vs_etaProbe");
  TH2F* mpf_vs_eta_DATA_with = (TH2F*)DATA_with->Get("Selection/mpf_vs_etaProbe");
  TH2F* r_rel_vs_eta_MC_with = (TH2F*)MC_with->Get("Selection/r_rel_vs_etaProbe");
  TH2F* r_rel_vs_eta_DATA_with = (TH2F*)DATA_with->Get("Selection/r_rel_vs_etaProbe");
  mpf_vs_eta_MC_with->RebinX(3);
  mpf_vs_eta_DATA_with->RebinX(3);
  r_rel_vs_eta_MC_with->RebinX(3);
  r_rel_vs_eta_DATA_with->RebinX(3);


  //2d-hists normalize to integral ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  mpf_vs_eta_MC_without->Scale(1/mpf_vs_eta_MC_without->Integral());
  mpf_vs_eta_DATA_without->Scale(1/mpf_vs_eta_DATA_without->Integral());
  r_rel_vs_eta_MC_without->Scale(1/r_rel_vs_eta_MC_without->Integral());
  r_rel_vs_eta_DATA_without->Scale(1/r_rel_vs_eta_DATA_without->Integral());
  mpf_vs_eta_MC_with->Scale(1/mpf_vs_eta_MC_with->Integral());
  mpf_vs_eta_DATA_with->Scale(1/mpf_vs_eta_DATA_with->Integral());
  r_rel_vs_eta_MC_with->Scale(1/r_rel_vs_eta_MC_with->Integral());
  r_rel_vs_eta_DATA_with->Scale(1/r_rel_vs_eta_DATA_with->Integral());

  //clone the ones to norm
  TH2F* mpf_vs_eta_MC_with_norm = (TH2F*)mpf_vs_eta_MC_with->Clone("mpf_vs_eta_MC_with_norm");
  TH2F* mpf_vs_eta_DATA_with_norm = (TH2F*)mpf_vs_eta_DATA_with->Clone("mpf_vs_eta_DATA_with_norm");
  TH2F* r_rel_vs_eta_MC_with_norm = (TH2F*)r_rel_vs_eta_MC_with->Clone("r_rel_vs_eta_MC_with_norm");
  TH2F* r_rel_vs_eta_DATA_with_norm = (TH2F*)r_rel_vs_eta_DATA_with->Clone("r_rel_vs_eta_DATA_with_norm");

  //normalize the ones with to the ones without
  mpf_vs_eta_MC_with_norm->Divide(mpf_vs_eta_MC_without);	 
  mpf_vs_eta_DATA_with_norm->Divide(mpf_vs_eta_DATA_without);	 
  r_rel_vs_eta_MC_with_norm->Divide(r_rel_vs_eta_MC_without);	 
  r_rel_vs_eta_DATA_with_norm->Divide(r_rel_vs_eta_DATA_without);

  
  




  TCanvas *a = new TCanvas();
  a->Divide(2,2);
  a->cd(1);
  r_rel_vs_eta_DATA_with_norm->GetZaxis()->SetRangeUser(0.8,1.2);
  r_rel_vs_eta_DATA_with_norm->SetTitle("DATA");
  r_rel_vs_eta_DATA_with_norm->Draw("colz");
  a->cd(2);
  r_rel_vs_eta_MC_with_norm->GetZaxis()->SetRangeUser(0.8,1.2);
  r_rel_vs_eta_MC_with_norm->SetTitle("MC");
  r_rel_vs_eta_MC_with_norm->Draw("colz");
  a->cd(3);
  mpf_vs_eta_DATA_with_norm->GetZaxis()->SetRangeUser(0.8,1.2);
  mpf_vs_eta_DATA_with_norm->SetTitle("DATA");
  mpf_vs_eta_DATA_with_norm->Draw("colz");
  a->cd(4);
  mpf_vs_eta_MC_with_norm->GetZaxis()->SetRangeUser(0.8,1.2);
  mpf_vs_eta_MC_with_norm->SetTitle("MC");
  mpf_vs_eta_MC_with_norm->Draw("colz");
  a->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_Responses_vs_Eta_norm.pdf");






 TCanvas *b = new TCanvas();
  b->Divide(2,2);
  b->cd(1);
  r_rel_vs_eta_DATA_with->SetTitle("DATA");
  r_rel_vs_eta_DATA_with->Draw("colz");
  b->cd(2);
  r_rel_vs_eta_MC_with->SetTitle("MC");
  r_rel_vs_eta_MC_with->Draw("colz");
  b->cd(3);
  mpf_vs_eta_DATA_with->SetTitle("DATA");
  mpf_vs_eta_DATA_with->Draw("colz");
  b->cd(4);
  mpf_vs_eta_MC_with->SetTitle("MC");
  mpf_vs_eta_MC_with->Draw("colz");
  b->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_Responses_vs_Eta_WithID.pdf");

 TCanvas *c = new TCanvas();
  c->Divide(2,2);
  c->cd(1);
  r_rel_vs_eta_DATA_without->SetTitle("DATA");
  r_rel_vs_eta_DATA_without->Draw("colz");
  c->cd(2);
  r_rel_vs_eta_MC_without->SetTitle("MC");
  r_rel_vs_eta_MC_without->Draw("colz");
  c->cd(3);
  mpf_vs_eta_DATA_without->SetTitle("DATA");
  mpf_vs_eta_DATA_without->Draw("colz");
  c->cd(4);
  mpf_vs_eta_MC_without->SetTitle("MC");
  mpf_vs_eta_MC_without->Draw("colz");
  c->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_Responses_vs_Eta_WithoutID.pdf");

  
  TCanvas *tmp = new TCanvas();
  r_rel_vs_eta_DATA_with_norm->Draw("colz");
  tmp->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_DATA_R_rel_vs_Eta_norm.pdf");
  r_rel_vs_eta_MC_with_norm->Draw("colz");
  tmp->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_MC_R_rel_vs_Eta_norm.pdf");
  mpf_vs_eta_DATA_with_norm->Draw("colz");
  tmp->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_DATA_Mpf_vs_Eta_norm.pdf");
  mpf_vs_eta_MC_with_norm->Draw("colz");
  tmp->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_MC_Mpf_vs_Eta_norm.pdf");

  r_rel_vs_eta_DATA_with->Draw("colz");
  tmp->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_DATA_R_rel_vs_Eta_WithID.pdf");
  r_rel_vs_eta_MC_with->Draw("colz");
  tmp->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_MC_R_rel_vs_Eta_WithID.pdf");
  mpf_vs_eta_DATA_with->Draw("colz");
  tmp->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_DATA_Mpf_vs_Eta_WithID.pdf");
  mpf_vs_eta_MC_with->Draw("colz");
  tmp->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_MC_Mpf_vs_Eta_WithID.pdf");

  r_rel_vs_eta_DATA_without->Draw("colz");
  tmp->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_DATA_R_rel_vs_Eta_WithoutID.pdf");
  r_rel_vs_eta_MC_without->Draw("colz");
  tmp->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_MC_R_rel_vs_Eta_WithoutID.pdf");
  mpf_vs_eta_DATA_without->Draw("colz");
  tmp->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_DATA_Mpf_vs_Eta_WithoutID.pdf");
  mpf_vs_eta_MC_without->Draw("colz");
  tmp->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_MC_Mpf_vs_Eta_WithoutID.pdf");


  /*
  //profiles
  TCanvas *c_prof_r_rel_data_norm = new TCanvas();
  TProfile *prof_r_rel_DATA_with_norm = (TProfile*)r_rel_vs_eta_DATA_with_norm->ProfileX();
  prof_r_rel_DATA_with_norm->SetTitle("DATA, with/without JetID, profile of normalized");
  prof_r_rel_DATA_with_norm->GetXaxis()->SetTitle("#eta probe");
  prof_r_rel_DATA_with_norm->GetYaxis()->SetTitle("p_{T} bal. response");
  prof_r_rel_DATA_with_norm->GetYaxis()->SetRangeUser(0.7,1.3);
  prof_r_rel_DATA_with_norm->Draw();
  c_prof_r_rel_data_norm->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_DATA_R_rel_vs_Eta_norm_profile.pdf");
  TCanvas *c_prof_r_rel_mc_norm = new TCanvas();
  TProfile *prof_r_rel_MC_with_norm = (TProfile*)r_rel_vs_eta_MC_with_norm->ProfileX();
  prof_r_rel_MC_with_norm->SetTitle("MC, with/without JetID, profile of normalized");
  prof_r_rel_MC_with_norm->GetXaxis()->SetTitle("#eta probe");
  prof_r_rel_MC_with_norm->GetYaxis()->SetTitle("p_{T} bal. response");
  prof_r_rel_MC_with_norm->GetYaxis()->SetRangeUser(0.7,1.3);
  prof_r_rel_MC_with_norm->Draw();
  c_prof_r_rel_mc_norm->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_MC_R_rel_vs_Eta_norm_profile.pdf");
  TCanvas *c_prof_mpf_data_norm = new TCanvas();
  TProfile *prof_mpf_DATA_with_norm = (TProfile*)mpf_vs_eta_DATA_with_norm->ProfileX();
  prof_mpf_DATA_with_norm->SetTitle("DATA, with/without JetID, profile of normalized");
  prof_mpf_DATA_with_norm->GetXaxis()->SetTitle("#eta probe");
  prof_mpf_DATA_with_norm->GetYaxis()->SetTitle("MPF response");
  prof_mpf_DATA_with_norm->GetYaxis()->SetRangeUser(0.7,1.3);
  prof_mpf_DATA_with_norm->Draw();
  c_prof_r_rel_data_norm->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_DATA_Mpf_vs_Eta_norm_profile.pdf");
  TCanvas *c_prof_mpf_mc_norm = new TCanvas();
  TProfile *prof_mpf_MC_with_norm = (TProfile*)mpf_vs_eta_MC_with_norm->ProfileX();
  prof_mpf_MC_with_norm->SetTitle("MC, with/without JetID, profile of normalized");
  prof_mpf_MC_with_norm->GetXaxis()->SetTitle("#eta probe");
  prof_mpf_MC_with_norm->GetYaxis()->SetTitle("MPF response");
  prof_mpf_MC_with_norm->GetYaxis()->SetRangeUser(0.7,1.3);
  prof_mpf_MC_with_norm->Draw();
  c_prof_r_rel_data_norm->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_MC_Mpf_vs_Eta_norm_profile.pdf");
  */

  TCanvas *c_prof_r_rel_data_with = new TCanvas("c_prof_r_rel_data_with", "Data, withID, rrel", 1);
  TProfile *prof_r_rel_DATA_with = (TProfile*)r_rel_vs_eta_DATA_with->ProfileX();
  prof_r_rel_DATA_with->SetTitle("DATA, with JetID, profile");
  prof_r_rel_DATA_with->GetXaxis()->SetTitle("#eta probe");
  prof_r_rel_DATA_with->GetYaxis()->SetTitle("p_{T} bal. response");
  prof_r_rel_DATA_with->GetYaxis()->SetRangeUser(0.7,1.3);
  prof_r_rel_DATA_with->Draw();
  c_prof_r_rel_data_with->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_DATA_R_rel_vs_Eta_WithID_profile.pdf");
  TCanvas *c_prof_r_rel_mc_with = new TCanvas("c_prof_r_rel_mc_with", "MC, withID, rrel", 1);
  TProfile *prof_r_rel_MC_with = (TProfile*)r_rel_vs_eta_MC_with->ProfileX();
  prof_r_rel_MC_with->SetTitle("MC, with JetID, profile");
  prof_r_rel_MC_with->GetXaxis()->SetTitle("#eta probe");
  prof_r_rel_MC_with->GetYaxis()->SetTitle("p_{T} bal. response");
  prof_r_rel_MC_with->GetYaxis()->SetRangeUser(0.7,1.3);
  prof_r_rel_MC_with->Draw();
  c_prof_r_rel_mc_with->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_MC_R_rel_vs_Eta_WithID_profile.pdf");
  TCanvas *c_prof_mpf_data_with = new TCanvas("c_prof_mpf_data_with", "Data, withID, mpf", 1);
  TProfile *prof_mpf_DATA_with = (TProfile*)mpf_vs_eta_DATA_with->ProfileX();
  prof_mpf_DATA_with->SetTitle("DATA, with JetID, profile");
  prof_mpf_DATA_with->GetXaxis()->SetTitle("#eta probe");
  prof_mpf_DATA_with->GetYaxis()->SetTitle("MPF response");
  prof_mpf_DATA_with->GetYaxis()->SetRangeUser(0.7,1.3);
  prof_mpf_DATA_with->Draw();
  c_prof_mpf_data_with->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_DATA_Mpf_vs_Eta_WithID_profile.pdf");
  TCanvas *c_prof_mpf_mc_with = new TCanvas("c_prof_mpf_mc_with", "MC, withID, mpf", 1);
  TProfile *prof_mpf_MC_with = (TProfile*)mpf_vs_eta_MC_with->ProfileX();
  prof_mpf_MC_with->SetTitle("MC, with JetID, profile");
  prof_mpf_MC_with->GetXaxis()->SetTitle("#eta probe");
  prof_mpf_MC_with->GetYaxis()->SetTitle("MPF response");
  prof_mpf_MC_with->GetYaxis()->SetRangeUser(0.7,1.3);
  prof_mpf_MC_with->Draw();
  c_prof_mpf_mc_with->SaveAs(path_with+"plots/JetIDPlots_"+Generator+"_"+Collection+"_MC_Mpf_vs_Eta_WithID_profile.pdf");



  TCanvas  *c_prof_r_rel_data_without = new TCanvas();
  TProfile *prof_r_rel_DATA_without   = (TProfile*)r_rel_vs_eta_DATA_without->ProfileX();
  prof_r_rel_DATA_without->SetTitle("DATA, without JetID, profile");
  prof_r_rel_DATA_without->GetXaxis()->SetTitle("#eta probe");
  prof_r_rel_DATA_without->GetYaxis()->SetTitle("p_{T} bal. response");
  prof_r_rel_DATA_without->GetYaxis()->SetRangeUser(0.7,1.3);
  prof_r_rel_DATA_without->Draw();
  c_prof_r_rel_data_without->SaveAs(path_without+"plots/JetIDPlots_"+Generator+"_"+Collection+"_DATA_R_rel_vs_Eta_WithoutID_profile.pdf");
  TCanvas  *c_prof_r_rel_mc_without = new TCanvas();
  TProfile *prof_r_rel_MC_without   = (TProfile*)r_rel_vs_eta_MC_without->ProfileX();
  prof_r_rel_MC_without->SetTitle("MC, without JetID, profile");
  prof_r_rel_MC_without->GetXaxis()->SetTitle("#eta probe");
  prof_r_rel_MC_without->GetYaxis()->SetTitle("p_{T} bal. response");
  prof_r_rel_MC_without->GetYaxis()->SetRangeUser(0.7,1.3);
  prof_r_rel_MC_without->Draw();
  c_prof_r_rel_mc_without->SaveAs(path_without+"plots/JetIDPlots_"+Generator+"_"+Collection+"_MC_R_rel_vs_Eta_WithoutID_profile.pdf");
  TCanvas  *c_prof_mpf_data_without = new TCanvas();
  TProfile *prof_mpf_DATA_without   = (TProfile*)mpf_vs_eta_DATA_without->ProfileX();
  prof_mpf_DATA_without->SetTitle("DATA, without JetID, profile");
  prof_mpf_DATA_without->GetXaxis()->SetTitle("#eta probe");
  prof_mpf_DATA_without->GetYaxis()->SetTitle("MPF response");
  prof_mpf_DATA_without->GetYaxis()->SetRangeUser(0.7,1.3);
  prof_mpf_DATA_without->Draw();
  c_prof_mpf_data_without->SaveAs(path_without+"plots/JetIDPlots_"+Generator+"_"+Collection+"_DATA_Mpf_vs_Eta_WithoutID_profile.pdf");
  TCanvas  *c_prof_mpf_mc_without = new TCanvas();
  TProfile *prof_mpf_MC_without   = (TProfile*)mpf_vs_eta_MC_without->ProfileX();
  prof_mpf_MC_without->SetTitle("MC, without JetID, profile");
  prof_mpf_MC_without->GetXaxis()->SetTitle("#eta probe");
  prof_mpf_MC_without->GetYaxis()->SetTitle("MPF response");
  prof_mpf_MC_without->GetYaxis()->SetRangeUser(0.7,1.3);
  prof_mpf_MC_without->Draw();
  c_prof_mpf_mc_without->SaveAs(path_without+"plots/JetIDPlots_"+Generator+"_"+Collection+"_MC_Mpf_vs_Eta_WithoutID_profile.pdf");


}
