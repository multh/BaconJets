#include "header.h"


TString ToString(int num) {
  ostringstream start;
  start<<num;
  TString start1=start.str();
  return start1;
}

void Rrel(TFile* datafile, TFile* MCfile){

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //gStyle->SetErrorX(0.0);


  TH2D* databins[n_eta-1];
  TH2D* qcdbins[n_eta-1];
  TH2D* data50bins[n_eta-1];

  for (int i=0; i<n_eta-1; i++){
    databins[i] = (TH2D*)datafile->Get("/a02/eta_"+eta_range[i]+"_"+eta_range[i+1]+"/Rrel_vs_Npv");
    qcdbins[i] = (TH2D*)MCfile->Get("/a02/eta_"+eta_range[i]+"_"+eta_range[i+1]+"/Rrel_vs_Npv");
  }
  
  TProfile* tpf[n_eta-1];
  TProfile* tpf_MC[n_eta-1];
  for (int i=0; i<n_eta-1; i++){
    TString numstr=ToString(i);
    TString plotname="eta_"+eta_range[i]+"_"+eta_range[i+1];
    tpf[i] = databins[i]->ProfileX(plotname,0,1000);
    tpf_MC[i] = qcdbins[i]->ProfileX(plotname+"MC",0,100);
  }
  

  for (int j=0; j<n_eta-1; j++){
    for (int i=0; i<41; i++){
      cout << "content " << tpf[j]->GetBinContent(i) << endl;
      cout << "error   " << tpf[j]->GetBinError(i) << endl;
      if(tpf[j]->GetBinError(i)==0){
	tpf[j]->SetBinContent(i,0);
      }
      if(tpf_MC[j]->GetBinError(i)==0){
	tpf_MC[j]->SetBinContent(i,0);
      }
    }
  }


  TLine *line = new TLine(0.,1,40.0,1);

  TCanvas* a[n_eta-1];
  for (int j=0; j<n_eta-1; j++){
    TString numstr=ToString(j);
    TString histname="eta_"+eta_range[j]+"_"+eta_range[j+1];
    a[j] = new TCanvas(histname, histname, 800,700);
    gStyle->SetOptTitle(0);
    tpf[j]->GetXaxis()->SetTitle("NPV");
    tpf[j]->GetYaxis()->SetTitle("R_{rel}");
    tpf[j]->GetXaxis()->SetTitleOffset(0.8);
    tpf[j]->GetYaxis()->SetTitleOffset(0.8);
    tpf[j]->GetXaxis()->SetTitleSize(0.050);
    tpf[j]->GetYaxis()->SetTitleSize(0.055);
    tpf[j]->GetXaxis()->SetRangeUser(0,40);
    tpf[j]->GetYaxis()->SetRangeUser(0.50,1.50);
    tpf[j]->SetMarkerStyle(20);
    tpf[j]->SetMarkerSize(1.3);
    tpf[j]->SetLineColor(1);
    tpf[j]->Draw();

    tpf_MC[j]->SetMarkerStyle(23);
    tpf_MC[j]->SetMarkerSize(1.3);
    tpf_MC[j]->SetMarkerColor(kRed);
    tpf_MC[j]->SetLineColor(kRed);
    tpf_MC[j]->Draw("SAME");

    line->SetLineStyle(2);
    line->Draw("SAME");

    TLegend *leg1;
    leg1 = new TLegend(0.15,0.65,0.35,0.85,"","brNDC");//x+0.1
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.045);
    leg1->SetFillColor(10);
    leg1->SetLineColor(1);
    leg1->SetTextFont(42);
    leg1->SetHeader(eta_range3[j]+"#leq|#eta|<"+eta_range3[j+1]);
    leg1->AddEntry(tpf[j], "DATA","P");
    leg1->AddEntry(tpf_MC[j], "QCD","P");
    leg1->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    tex->DrawLatex(0.64,0.91,"2.11fb^{-1} (13TeV)");

    a[j]->Print("plots/Rrel_vs_NPV_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");
  }



  /*
  TCanvas * c = new TCanvas("c", "c", 800,700);
  tpf->GetXaxis()->SetTitle("NPV");
  tpf->GetYaxis()->SetTitle("R_{rel}");
  tpf->GetXaxis()->SetTitleSize(0.045);
  tpf->GetYaxis()->SetTitleSize(0.045);
  tpf->GetYaxis()->SetRangeUser(0.80,1.20);
  tpf->SetMarkerStyle(20);
  tpf->SetMarkerSize(1.3);
  tpf->SetLineColor(1);
  tpf->Draw();
  tpf_MC->SetMarkerStyle(23);
  tpf_MC->SetMarkerSize(1.3);
  tpf_MC->SetMarkerColor(kRed);
  tpf_MC->SetLineColor(kRed);
  tpf_MC->Draw("SAME");
 
  line->SetLineStyle(2);
  line->Draw("SAME");

  TLegend *leg1;
  leg1 = new TLegend(0.20,0.70,0.45,0.85);//x+0.1
  leg1 -> SetBorderSize(0);
  leg1 -> SetTextSize(0.042);
  leg1 -> SetFillColor(10);
  leg1 -> SetLineColor(1);
  leg1 -> SetTextFont(42);
  TLegendEntry* entries[2];
  leg1->SetHeader(eta_range3[a]+"<|#eta|<"+eta_range3[a+1]+", p_{T}^{corr}");
  entries[1]=leg1 -> AddEntry(tpf, "DATA","P");
  entries[0]=leg1 -> AddEntry(tpf_MC, "MC","P");
  leg1->Draw();


  c->SaveAs("rrelvsnpv/Rrel_vs_NPV_eta320to500_ptCorr_L1.pdf");
  */


}

