//CompMacro.C
#include "../include/parameters.h"

//***********************************************************************
//
//Macro for Background Composition after every Selectionstep
//
//***********************************************************************


void CompMacro() {

  //  c1->SetOptTitle(0);
 
   cout<<"Before files"<<endl;
 // open a file and get a histogram
 TFile *mpf_new = new TFile("/nfs/dust/cms/user/multh/JEC/2016Legacy/Residuals/Summer16_07Aug2017_V6/AK4CHS/MC_NoReweighted_NewHF_Binning/RunBCD/plots/pT_extrapolation_mpf.root");
 TFile *pt_new  = new TFile("/nfs/dust/cms/user/multh/JEC/2016Legacy/Residuals/Summer16_07Aug2017_V6/AK4CHS/MC_NoReweighted_NewHF_Binning/RunBCD/plots/pT_extrapolation_pt.root");
 TFile *mpf_old = new TFile("/nfs/dust/cms/user/multh/JEC/2016Legacy/Residuals/Summer16_07Aug2017_V3/AK4CHS/MC_NoReweighted_NewTriggerThresholds/RunBCD/plots/pT_extrapolation_mpf.root");
 TFile *pt_old  = new TFile("/nfs/dust/cms/user/multh/JEC/2016Legacy/Residuals/Summer16_07Aug2017_V3/AK4CHS/MC_NoReweighted_NewTriggerThresholds/RunBCD/plots/pT_extrapolation_pt.root");
    cout<<"After files"<<endl;

 TH1F *h_pt_extrapolation_mpf_old[n_eta-1];
 TH1F *h_pt_extrapolation_mpf_new[n_eta-1];
 TH1F *h_pt_extrapolation_pT_old[n_eta-1];
 TH1F *h_pt_extrapolation_pT_new[n_eta-1];

   cout<<"Before Loop"<<endl;
 for(int i=0; i<n_eta-1; i++){
 h_pt_extrapolation_mpf_old[i] = new TH1F("pt_extrapolation_old"+eta_range[i]+"_"+eta_range[i+1],"R_MC/R_Data",10,0,1);
 h_pt_extrapolation_mpf_new[i] = new TH1F("pt_extrapolation_new"+eta_range[i]+"_"+eta_range[i+1],"R_MC/R_Data",10,0,1);

 cout<<"create Hists"<<endl;
TCanvas *c1 = new TCanvas("c1","c1",650,500);
   h_pt_extrapolation_mpf_old[i]    = (TH1F*)mpf_old->Get("pTextrapolation_MPF_pythia8_pT_"+eta_range2[i]+"_"+eta_range2[i+1]);

   h_pt_extrapolation_mpf_new[i]    = (TH1F*)mpf_new->Get("pTextrapolation_MPF_pythia8_pT_"+eta_range2[i]+"_"+eta_range2[i+1]);
    cout<<"After Hists files"<<endl;

  c1->SetLogx();
 
 // h_pt_extrapolation_mpf_old[i]->GetYaxis()->SetTitle("R");
 h_pt_extrapolation_mpf_old[i]->SetMarkerColor(kRed);
 h_pt_extrapolation_mpf_old[i]->SetLineColor(kRed);
 h_pt_extrapolation_mpf_old[i]->SetMarkerStyle(20);

 h_pt_extrapolation_mpf_old[i]->Draw("AP");
 cout<<"Draw First Hist"<<endl;
 // h_pt_extrapolation_mpf_new->GetYaxis()->SetTitle("R_MC/R_Data");
 h_pt_extrapolation_mpf_new[i]->SetMarkerColor(kBlue);
 h_pt_extrapolation_mpf_new[i]->SetLineColor(kBlue);
 h_pt_extrapolation_mpf_new[i]->SetMarkerStyle(20);
 //h_pt_extrapolation_mpf_new[i]->GetFunction("dijet_ptextra_pythia8_eta_"+eta_range[i]+"_"+eta_range[i+1]+"f1")->SetLineColor(kBlue);
 h_pt_extrapolation_mpf_new[i]->Draw("same P");

cout<<"Draw Second Hist"<<endl;

 if(i==13){ TLegend* leg1 = new TLegend(0.12,0.17,0.47,0.3);
  leg1 -> AddEntry( h_pt_extrapolation_mpf_new[i] ,"MPF New","LP");  
  leg1 -> AddEntry( h_pt_extrapolation_mpf_old[i] ,"MPF previous","LP"); 
  leg1->SetTextSize(0.05);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
 
  leg1->Draw();
}
 else{ TLegend* leg1 = new TLegend(0.12,0.8,0.4,0.7);
  leg1 -> AddEntry( h_pt_extrapolation_mpf_new[i] ,"MPF New","LP");  
  leg1 -> AddEntry( h_pt_extrapolation_mpf_old[i] ,"MPF previous","LP"); 
  leg1->SetTextSize(0.05);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
 
  leg1->Draw();
 }


 TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    TString text = eta_range[i] + " < |#eta| < " + eta_range[i+1];
    tex->DrawLatex(0.52,0.85, text);
 
  c1->SaveAs("pTextrapolation_MPF_eta_"+ eta_range2[i] + "_" + eta_range2[i+1] + ".pdf");
cout<<"Save"<<endl;
 }

   cout<<"Before Pt balance"<<endl;

 for(int i=0; i<n_eta-1; i++){
 h_pt_extrapolation_pT_old[i] = new TH1F("pt_extrapolation_old"+eta_range2[i]+"_"+eta_range2[i+1],"R_MC/R_Data",10,0,1);
 h_pt_extrapolation_pT_new[i] = new TH1F("pt_extrapolation_new"+eta_range2[i]+"_"+eta_range2[i+1],"R_MC/R_Data",10,0,1);

 cout<<"create Hists"<<endl;
TCanvas *c1 = new TCanvas("c1","c1",650,500);


   h_pt_extrapolation_pT_old[i]    = (TH1F*)pt_old->Get("pTextrapolation_Pt_pythia8_pT_"+eta_range2[i]+"_"+eta_range2[i+1]);

   h_pt_extrapolation_pT_new[i]    = (TH1F*)pt_new->Get("pTextrapolation_Pt_pythia8_pT_"+eta_range2[i]+"_"+eta_range2[i+1]);
 

 c1->SetLogx();
 
 // h_pt_extrapolation_pT_old[i]->GetYaxis()->SetTitle("R");
 h_pt_extrapolation_pT_old[i]->SetMarkerColor(kRed);
 h_pt_extrapolation_pT_old[i]->SetLineColor(kRed);
 h_pt_extrapolation_pT_old[i]->SetMarkerStyle(20);

 h_pt_extrapolation_pT_old[i]->Draw("AP");
 cout<<"Draw First Hist"<<endl;
 // h_pt_extrapolation_pT_new->GetYaxis()->SetTitle("R_MC/R_Data");
 h_pt_extrapolation_pT_new[i]->SetMarkerColor(kBlue);
 h_pt_extrapolation_pT_new[i]->SetLineColor(kBlue);
 h_pt_extrapolation_pT_new[i]->SetMarkerStyle(20);
 h_pt_extrapolation_pT_new[i]->Draw("same P");
cout<<"Draw Second Hist"<<endl;

 if(i==13){ TLegend* leg1 = new TLegend(0.12,0.17,0.47,0.3);
  leg1 -> AddEntry( h_pt_extrapolation_pT_new[i] ,"Pt balance New","LP");  
  leg1 -> AddEntry( h_pt_extrapolation_pT_old[i] ,"Pt balance previous","LP"); 
  leg1->SetTextSize(0.05);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
 
  leg1->Draw();
}
 else{ TLegend* leg1 = new TLegend(0.12,0.8,0.4,0.7);
  leg1 -> AddEntry( h_pt_extrapolation_pT_new[i] ,"Pt balance New","LP");  
  leg1 -> AddEntry( h_pt_extrapolation_pT_old[i] ,"Pt balance previous","LP"); 
  leg1->SetTextSize(0.05);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
 
  leg1->Draw();
 }


 TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    TString text = eta_range[i] + " < |#eta| < " + eta_range[i+1];
    tex->DrawLatex(0.52,0.85, text);
 
  c1->SaveAs("pTextrapolation_Pt_eta_"+ eta_range2[i] + "_" + eta_range2[i+1] + ".pdf");
cout<<"Save"<<endl;
 }



}
