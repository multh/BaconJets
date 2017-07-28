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

void CorrectionObject::L2ResAllRuns(){
  cout << "--------------- Starting L2ResAllRuns() ---------------" << endl << endl;

  TCanvas* c1 = new TCanvas();
  TH1D *h = new TH1D("h",";|#eta|;Relative correction",41,0,5.191);
  h->SetMaximum(1.2); //1.2
  h->SetMinimum(0.8); //0.8
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  tdrCanvas(c1,"c1",h,4,10,true,CorrectionObject::_lumitag);

  TLegend leg1 = tdrLeg(0.40,0.19,0.65,0.4);

  TLine *line = new TLine(0,1,5.191,1);              
                            
  TString runnr_v[5]={"BCD","EFearly","FlateG","H", "BCDEFGH"};
  TCanvas* c2 = new TCanvas;
  tdrCanvas(c2,"c2",h,4,10,true,"All Runs");
  TFile* f_Res_mpf;
  TH1D* res_logpt_mpf_kfsrfit;
  for(int j=0;j<5;j++){
    TString runnr=runnr_v[j];
    TString path =CorrectionObject::_input_path+"abs_eta/Run" + runnr + "/" + "Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root";
    f_Res_mpf = new TFile(path,"READ");   
    res_logpt_mpf_kfsrfit = (TH1D*)f_Res_mpf->Get("res_logpt_mpf");
    res_logpt_mpf_kfsrfit->SetMarkerStyle(1);
    res_logpt_mpf_kfsrfit->SetLineWidth(2);
    res_logpt_mpf_kfsrfit->SetLineColor(1+j); 
    if(j == 2) res_logpt_mpf_kfsrfit->SetLineColor(kGreen-2); 
    if(j == 4) res_logpt_mpf_kfsrfit->SetLineColor(kCyan-3); 
    leg1.AddEntry(res_logpt_mpf_kfsrfit,runnr,"L");
    res_logpt_mpf_kfsrfit->Draw("E SAME");  

  }

  line->Draw("SAME");
  leg1.Draw();

  TString save_as = CorrectionObject::_input_path+"/L2Res_MPF_LOGLIN_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_nominal_ALL_runs_BCD_EFearly_FlateG_H.pdf";
  cout << "save-as = " << save_as << endl;
  c2->SaveAs(save_as); 





  //delete everything
  delete res_logpt_mpf_kfsrfit;
  delete f_Res_mpf;
  //delete c2;
  delete line;
  //delete leg1;
  delete h;
  delete c1;


}
