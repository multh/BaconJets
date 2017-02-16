#include "../include/CorrectionObject.h"
#include "../include/parameters.h"
#include "../include/useful_functions.h"

#include <TH1.h>
#include <TH1D.h>
#include <TString.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TFile.h>
#include <TLine.h>

#include "../include/tdrstyle_mod15.h"






using namespace std;

//Function to overlay current results with another set of corrections
void CorrectionObject::L2ResOverlay(){
  cout << "--------------- Starting L2ResOverlay() ---------------" << endl << endl;

  TCanvas* c1 = new TCanvas();
  TH1D *h = new TH1D("h",";|#eta|;Relative correction",41,0,5.191);
  h->SetMaximum(1.2); //1.2
  h->SetMinimum(0.8); //0.8
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  tdrCanvas(c1,"c1",h,4,10,true,CorrectionObject::_lumitag);

  TLegend leg1 = tdrLeg(0.17,0.19,0.35,0.5);

  TLine *line = new TLine(0.,1,5.191,1);              

  const int n_runs = 1;
  //TString runnr_v[n_runs]={"BCD","EFearly","FlateG","H"};
  TString runnr_v[n_runs]={"EFearly"};
  TCanvas* c2 = new TCanvas;
  tdrCanvas(c2,"c2",h,4,10,true,"V3 vs. newest");
  TFile* f_Res_mpf;
  TFile* f_Res_mpf_old;
  TH1D* res_logpt_mpf_kfsrfit;
  TH1D* res_logpt_mpf_kfsrfit_old;
  TH1D* res_diff[n_runs];
  for(int j=0;j<n_runs;j++){
    TString runnr=runnr_v[j];
    TString path = CorrectionObject::_outpath + "../Run" + runnr + "/" + "Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root";
    TString path_old = "/nfs/dust/cms/user/reimersa/JEC/2016ReReco/Residuals/Summer16_23Sep2016_V1/AK4CHS/AveResponse/Run" + runnr + "/" + "Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root";
    f_Res_mpf = new TFile(path,"READ");   
    f_Res_mpf_old = new TFile(path_old,"READ");  
 
    res_logpt_mpf_kfsrfit = (TH1D*)f_Res_mpf->Get("res_logpt_mpf");
    res_logpt_mpf_kfsrfit_old = (TH1D*)f_Res_mpf_old->Get("res_logpt_mpf");
    res_diff[j] = (TH1D*)res_logpt_mpf_kfsrfit_old->Clone();
    res_diff[j]->Add(res_logpt_mpf_kfsrfit,-1);

    res_logpt_mpf_kfsrfit->SetMarkerStyle(1);
    res_logpt_mpf_kfsrfit_old->SetMarkerStyle(1); 
    res_diff[j]->SetMarkerStyle(1); 
    res_logpt_mpf_kfsrfit_old->SetLineStyle(2); // dashed line 
    res_logpt_mpf_kfsrfit->SetLineWidth(2);
    res_logpt_mpf_kfsrfit_old->SetLineWidth(2);
    res_diff[j]->SetLineWidth(2);
    res_logpt_mpf_kfsrfit->SetLineColor(1+j); 
    res_logpt_mpf_kfsrfit_old->SetLineColor(1+j); 
    res_diff[j]->SetLineColor(1+j); 
    if(j == 2){
      res_logpt_mpf_kfsrfit->SetLineColor(kGreen-2); 
      res_logpt_mpf_kfsrfit_old->SetLineColor(kGreen-2);
      res_diff[j]->SetLineColor(kGreen-2);
    }
    leg1.AddEntry(res_logpt_mpf_kfsrfit,runnr+"newest","L");
    leg1.AddEntry(res_logpt_mpf_kfsrfit_old,runnr+"V3","L");
    res_logpt_mpf_kfsrfit->Draw("E SAME");  
    res_logpt_mpf_kfsrfit_old->Draw("E SAME");  
  }
  
  line->Draw("SAME");
  leg1.Draw();

  TString save_as = CorrectionObject::_outpath + "../plots/L2Res_MPF_LOGLIN_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_nominal_VersionComparison_Summer16_23Sep2016_V3_vs_CorrectFormulae_MCReweighting.pdf";
  cout << "save-as = " << save_as << endl;
  c2->SaveAs(save_as); 

  TCanvas* c3 = new TCanvas;
  h->SetMaximum(0.3); 
  h->SetMinimum(-0.3); 
  h->GetYaxis()->SetTitle("V3 - 'newest'");
  tdrCanvas(c3,"c3",h,4,10,true,"V3 - newest");
  TLegend leg2 = tdrLeg(0.17,0.19,0.35,0.5);
  for(int i=0; i<n_runs; i++){
    TString runnr=runnr_v[i];
    leg2.AddEntry(res_diff[i],runnr+": (V3-newest)","L");
    res_diff[i]->Draw("E SAME");
  }
  line->Draw("SAME");
  leg2.Draw();

  TString save_as_diff = CorrectionObject::_outpath + "../plots/L2Res_MPF_LOGLIN_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_nominal_VersionComparison_Summer16_23Sep2016_V3_vs_CorrectFormulae_MCReweighting_diff.pdf";
  cout << "save-as_diff = " << save_as_diff << endl;
  c3->SaveAs(save_as_diff); 




  //delete everything
  for(int i=0; i<n_runs; i++) delete res_diff[i];
  delete res_logpt_mpf_kfsrfit;
  delete res_logpt_mpf_kfsrfit_old;
  delete f_Res_mpf;
  delete f_Res_mpf_old;
  //delete c2;
  delete line;
  //delete leg1;
  delete h;
  delete c1;


}
