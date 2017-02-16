#include "../include/CorrectionObject.h"
#include "../include/parameters.h"
#include "../include/useful_functions.h"

#include <TStyle.h>
#include <TH1.h>
#include <TH1D.h>
#include <TString.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TFile.h>
#include <TLine.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>



using namespace std;

void CorrectionObject::kFSR(){
  cout << "--------------- Starting kFSR() ---------------" << endl << endl;
  TStyle* m_gStyle = new TStyle();
  m_gStyle->SetOptFit(0);

// get ratio for MC to DATA responses
  double ratio_al_rel_r[n_pt-1][n_eta-1][n_alpha]; //ratio at pt,eta,alpha bins
  double err_ratio_al_rel_r[n_pt-1][n_eta-1][n_alpha]; //error of ratio at pt,eta,alpha bins
  double ratio_al_mpf_r[n_pt-1][n_eta-1][n_alpha]; //ratio at pt,eta,alpha bins
  double err_ratio_al_mpf_r[n_pt-1][n_eta-1][n_alpha]; //error of ratio at pt,eta,alpha bins
  TH1D *hdata_rel_r[n_pt-1][n_eta-1][n_alpha];// pT-balance response for data
  TH1D *hdata_mpf_r[n_pt-1][n_eta-1][n_alpha];//MPF response for data
  TH1D *hmc_rel_r[n_pt-1][n_eta-1][n_alpha];// pT-balanse responce for MC
  TH1D *hmc_mpf_r[n_pt-1][n_eta-1][n_alpha];//MPF response for MC
  int count = 0;
  TString name1 = "hist_data_rel_r_";
  TString name2 = "hist_data_mpf_r_";
  TString name3 = "hist_mc_rel_r_";
  TString name4 = "hist_mc_mpf_r_";
  for(int i=0; i<n_alpha; i++){
    for(int j=0; j<n_eta-1; j++){
      for(int k=0; k<n_pt-1; k++){
	ratio_al_rel_r[k][j][i] = 0;
	err_ratio_al_rel_r[k][j][i] = 0;
	ratio_al_mpf_r[k][j][i] = 0;
	err_ratio_al_mpf_r[k][j][i] = 0;
	TString name = name1; name+=count;
	hdata_rel_r[k][j][i] = new TH1D(name,"",nResponseBins, 0, 2.5);
	name = name2;name+=count;
	hdata_mpf_r[k][j][i] = new TH1D(name,"",nResponseBins, 0, 2.5);
	name = name3; name+=count;
	hmc_rel_r[k][j][i] = new TH1D(name,"",nResponseBins, 0, 2.5);
	name = name4; name+=count;
	hmc_mpf_r[k][j][i] = new TH1D(name,"",nResponseBins, 0, 2.5);
	count++;
      }
    }
  }
  cout << "Set up a total of " << count << " histograms." << endl;
  
  // Get relevant quantities from DATA, loop over events
  TTreeReader myReader_DATA("AnalysisTree", CorrectionObject::_DATAFile);
  TTreeReaderValue<Float_t> pt_ave_data(myReader_DATA, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_data(myReader_DATA, "probejet_eta");
  TTreeReaderValue<Float_t> alpha_data(myReader_DATA, "alpha");
  TTreeReaderValue<Float_t> rel_r_data(myReader_DATA, "rel_r");
  TTreeReaderValue<Float_t> mpf_r_data(myReader_DATA, "mpf_r");   
  TTreeReaderValue<Float_t> weight_data(myReader_DATA, "weight");
  int idx = 0;
  
  cout << "starting to loop over DATA events." << endl;

  while (myReader_DATA.Next()) {
    for(int k=0; k<n_pt-1; k++){
       if(*pt_ave_data<pt_bins[k] || *pt_ave_data>pt_bins[k+1]) continue;
       for(int j=0; j<n_eta-1; j++){
	 if(fabs(*probejet_eta_data)>eta_bins[j+1] || fabs(*probejet_eta_data)<eta_bins[j]) continue;
	     
	 for(int i=0; i<n_alpha; i++){
	   if(*alpha_data>alpha_bins[i]) continue;
	   else{
	     hdata_rel_r[k][j][i]->Fill(*rel_r_data,*weight_data);
	     hdata_mpf_r[k][j][i]->Fill(*mpf_r_data,*weight_data);
	     idx++;
	     if(idx%1000000==0) cout << "looping over data-TTree: Idx = " << idx << endl;
	   }
	 }
       }
     }
   }
  cout << "Finished running over DATA events. Read in a total of " << idx << " events." << endl;

   // Get relevant quantities from MC, loop over events
   TTreeReader myReader_MC("AnalysisTree", CorrectionObject::_MCFile);
   TTreeReaderValue<Float_t> pt_ave_mc(myReader_MC, "pt_ave");
   TTreeReaderValue<Float_t> probejet_eta_mc(myReader_MC, "probejet_eta");
   TTreeReaderValue<Float_t> alpha_mc(myReader_MC, "alpha");
   TTreeReaderValue<Float_t> rel_r_mc(myReader_MC, "rel_r");
   TTreeReaderValue<Float_t> mpf_r_mc(myReader_MC, "mpf_r");
   TTreeReaderValue<Float_t> weight_mc(myReader_MC, "weight");
   idx = 0;

   while (myReader_MC.Next()) {
     for(int k=0; k<n_pt-1; k++){
       if(*pt_ave_mc<pt_bins[k] || *pt_ave_mc>pt_bins[k+1]) continue;
       for(int j=0; j<n_eta-1; j++){
   	 if(fabs(*probejet_eta_mc)>eta_bins[j+1] || fabs(*probejet_eta_mc)<eta_bins[j]) continue;
   	 for(int i=0; i<n_alpha; i++){
   	   if(*alpha_mc>alpha_bins[i]) continue;
   	   else{
   	     hmc_rel_r[k][j][i]->Fill(*rel_r_mc,*weight_mc);
   	     hmc_mpf_r[k][j][i]->Fill(*mpf_r_mc,*weight_mc);
   	     idx++;
   	     if(idx%1000000==0) cout << "looping over MC-TTree: Idx = " << idx << endl;
   	   }
   	 }
       }
     }
   }



   //normalize histograms to their integral
   for(int i=0; i<n_alpha; i++){
     for(int j=0; j<n_eta-1; j++){
       for(int k=0; k<n_pt-1; k++){
	 double integ_rel_mc = hmc_rel_r[k][j][i]->Integral();
	 hmc_rel_r[k][j][i]->Scale(1./integ_rel_mc);
	 double integ_mpf_mc = hmc_mpf_r[k][j][i]->Integral();
	 hmc_mpf_r[k][j][i]->Scale(1./integ_mpf_mc);
	 double integ_rel_data = hdata_rel_r[k][j][i]->Integral();
	 hdata_rel_r[k][j][i]->Scale(1./integ_rel_data);
	 double integ_mpf_data = hdata_mpf_r[k][j][i]->Integral();
	 hdata_mpf_r[k][j][i]->Scale(1./integ_mpf_data);
       }
     }
   }

   //Get mean responses and respective errors from filled histograms
   for(int i=0; i<n_alpha; i++){
     for(int j=0; j<n_eta-1; j++){
       int n_sum = 0;
       for(int k=0; k<n_pt-1; k++){
	 if(k==0 && j==0) cout << "Alpha-bin from " << alpha_bins[i] << " to " << alpha_bins[i+1] << endl;
	 if(k==0) cout<<eta_bins[j]<<" ";
	 cout<<"& "<<hmc_rel_r[k][j][i]->GetEntries();//MC
	 n_sum += hmc_rel_r[k][j][i]->GetEntries(); //MC
	 //cout<<"& "<<hdata_rel_r[k][j][i]->GetEntries();//DATA
	 //n_sum += hdata_rel_r[k][j][i]->GetEntries(); //DATA

	 pair<double,double> res_mc_rel_r,res_data_rel_r;
	 pair<double,double> res_mc_mpf_r,res_data_mpf_r;
	 res_mc_rel_r = GetValueAndError(hmc_rel_r[k][j][i]);
	 res_data_rel_r = GetValueAndError(hdata_rel_r[k][j][i]);
	 res_mc_mpf_r = GetValueAndError(hmc_mpf_r[k][j][i]);
	 res_data_mpf_r = GetValueAndError(hdata_mpf_r[k][j][i]);

	 pair<double,double> ratio_res_rel_r;
	 if(res_mc_rel_r.first>0 && res_data_rel_r.first>0)
	   ratio_res_rel_r = Rmc_to_Rdata(res_mc_rel_r,res_data_rel_r);
	 else 
	   ratio_res_rel_r.first = 0;
	 pair<double,double> ratio_res_mpf_r;
	 if(res_mc_mpf_r.first>0 && res_data_mpf_r.first>0)
	   ratio_res_mpf_r = Rmc_to_Rdata(res_mc_mpf_r,res_data_mpf_r);
	 else 
	   ratio_res_mpf_r.first = 0;

	 ratio_al_rel_r[k][j][i] = ratio_res_rel_r.first;
	 err_ratio_al_rel_r[k][j][i] = ratio_res_rel_r.second;
	 ratio_al_mpf_r[k][j][i] = ratio_res_mpf_r.first;
	 err_ratio_al_mpf_r[k][j][i] = ratio_res_mpf_r.second;
       }
       cout << ",   SUM: " << n_sum <<endl;
     }
   }



   //Normalization of hists to value at alpha = alpha_cut
 
  //1) find bin with alpha = alpha_cut: bin no. al_ref
   int al_ref=0;
   for(int i=0; i<n_alpha; i++){
     if(fabs(alpha_bins[i]-alpha_cut)<1e-4) 
       al_ref=i;
   }

   //2) Normalize values and errors of responses to value at alpha = alpha_cut
   
   for(int j=0; j<n_eta-1; j++){
     for(int k=0; k<n_pt-1; k++){
       double norm_alref_rel_r = ratio_al_rel_r[k][j][al_ref];
       double err_norm_alref_rel_r = err_ratio_al_rel_r[k][j][al_ref];
       double norm_alref_mpf_r = ratio_al_mpf_r[k][j][al_ref];
       double err_norm_alref_mpf_r = err_ratio_al_mpf_r[k][j][al_ref];
       for(int i=0; i<n_alpha; i++){
	 if(norm_alref_rel_r>0){ //WHAT IS HAPPENING HERE? NO PROPER ERROR PROPAGATION !?!
	   ratio_al_rel_r[k][j][i] =   ratio_al_rel_r[k][j][i]/norm_alref_rel_r; //original
	   err_ratio_al_rel_r[k][j][i] = sqrt(abs(pow(err_ratio_al_rel_r[k][j][i],2)-pow(err_norm_alref_rel_r,2)));
	   //err_ratio_al_rel_r[k][j][i] = sqrt(abs(pow(err_ratio_al_rel_r[k][j][i] / (ratio_al_rel_r[k][j][i]) ,2)+pow(err_norm_alref_rel_r / norm_alref_rel_r,2))) * ratio_al_rel_r[k][j][i] / norm_alref_rel_r ; //self
	   //ratio_al_rel_r[k][j][i] =   ratio_al_rel_r[k][j][i]/norm_alref_rel_r;
	   if(i == al_ref) err_ratio_al_rel_r[k][j][i] = 0.;
	 }
	 if(norm_alref_mpf_r>0){
	   ratio_al_mpf_r[k][j][i] =   ratio_al_mpf_r[k][j][i]/norm_alref_mpf_r;
	   err_ratio_al_mpf_r[k][j][i] = sqrt(abs(pow(err_ratio_al_mpf_r[k][j][i],2)-pow(err_norm_alref_mpf_r,2)));
	   //err_ratio_al_mpf_r[k][j][i] = sqrt(abs(pow(err_ratio_al_mpf_r[k][j][i] / (ratio_al_mpf_r[k][j][i]) ,2)+pow(err_norm_alref_mpf_r / norm_alref_mpf_r,2))) * err_ratio_al_mpf_r[k][j][i] / norm_alref_mpf_r;
	   //ratio_al_mpf_r[k][j][i] =   ratio_al_mpf_r[k][j][i]/norm_alref_mpf_r;
	   if(i == al_ref) err_ratio_al_mpf_r[k][j][i] = 0.;
	 }
       }
     }
   }


   // Build the Multigraphs containing the ratio of responses (MC/DATA) as a function of alpha
   TGraphErrors *graph_rel_r[n_pt-1][n_eta-1];  //set of points vs alpha
   TMultiGraph *pTgraph_rel_r[n_eta-1];         //set of different pT bins in on eta bin
   TGraphErrors *graph_mpf_r[n_pt-1][n_eta-1];  //set of points vs alpha
   TMultiGraph *pTgraph_mpf_r[n_eta-1];         //set of different pT bins in on eta bin

   //Define legend
   TLegend *leg1;
   leg1 = new TLegend(0.17,0.68,0.65,0.89,"","brNDC");//x+0.1
   leg1->SetBorderSize(0);
   leg1->SetTextSize(0.03);
   leg1->SetFillColor(10);
   leg1->SetLineColor(1);
   leg1->SetTextFont(42);
   leg1->SetNColumns(2);


   double xbin_tgraph[n_alpha],zero[n_alpha];
   for(int i=0;i<n_alpha;i++){
     xbin_tgraph[i] = alpha_bins[i];
     zero[i] = 0;
   }


   for(int j=0; j<n_eta-1; j++){
     pTgraph_rel_r[j] = new TMultiGraph();
     pTgraph_mpf_r[j] = new TMultiGraph();
     for(int k=0; k<n_pt-1; k++){

       graph_rel_r[k][j] = new TGraphErrors(n_alpha,xbin_tgraph,ratio_al_rel_r[k][j],zero,err_ratio_al_rel_r[k][j]);
       graph_rel_r[k][j] = (TGraphErrors*)CleanEmptyPoints(graph_rel_r[k][j]);

       graph_rel_r[k][j]->SetMarkerSize(1.5);
       graph_rel_r[k][j]->SetMarkerStyle(20);
       if(k<9){
	 graph_rel_r[k][j]->SetMarkerColor(k+1);
	 graph_rel_r[k][j]->SetLineColor(k+1);
       }
       else{
	 graph_rel_r[k][j]->SetMarkerColor(k+19);
	 graph_rel_r[k][j]->SetLineColor(k+19);
       }
       if(graph_rel_r[k][j]->GetN()>0) 
	 {
	   pTgraph_rel_r[j]->Add(graph_rel_r[k][j]); //one multigraph consisting of several TGraphErrors. One Multigraph for each eta bin
	 }

       TString pTbin_label = "";
       pTbin_label+=pt_bins[k];
       pTbin_label+=" < p_{T} < ";
       pTbin_label+=pt_bins[k+1];
       if(j==0) leg1->AddEntry(graph_rel_r[k][j],pTbin_label,"epl");
      
       graph_mpf_r[k][j] = new TGraphErrors(n_alpha,xbin_tgraph,ratio_al_mpf_r[k][j],zero,err_ratio_al_mpf_r[k][j]);
       graph_mpf_r[k][j] = (TGraphErrors*)CleanEmptyPoints(graph_mpf_r[k][j]);
       graph_mpf_r[k][j]->SetMarkerSize(1.3);
       graph_mpf_r[k][j]->SetMarkerStyle(20);
       if(k<9){
	 graph_mpf_r[k][j]->SetMarkerColor(k+1);
	 graph_mpf_r[k][j]->SetLineColor(k+1);
       }
       else{
	 graph_mpf_r[k][j]->SetMarkerColor(k+19);
	 graph_mpf_r[k][j]->SetLineColor(k+19);
       }
       if(graph_mpf_r[k][j]->GetN()>0) 
	 pTgraph_mpf_r[j]->Add(graph_mpf_r[k][j]);
     }
   }

   /* +++++++++++++++++ PLOTS ++++++++++++++++++ */
   //Do the well-known, colorful kFSR linear fit plots
  

   //Create horizontal line for plotting ("ideal value")
   TLine *line = new TLine(alpha_bins[0]-0.01,1,alpha_bins[n_alpha-1]+0.01,1);

   //First set up the output files   
   //Create output .dat file, including the kFSR extrapolation (alpha->0)
   FILE *fp_rel_r, *fp_mpf_r; 
   TH1D *kFSR_MPF, *kFSR_DiJet, *plotkfsr;
   
   cout << "Opening .dat files at: " << CorrectionObject::_outpath+"output/KFSR_MPF_extrapolation.dat" << endl;
   fp_mpf_r = fopen(CorrectionObject::_outpath+"output/KFSR_MPF_extrapolation.dat","w");
   fp_rel_r = fopen(CorrectionObject::_outpath+"output/KFSR_DiJet_extrapolation.dat","w");

   kFSR_MPF = new TH1D("kfsr_mpf","kfsr_mpf", n_eta-1,eta_bins);
   kFSR_DiJet = new TH1D("kfsr_dijet","kfsr_dijet", n_eta-1,eta_bins);
   plotkfsr = new TH1D("kfsr","kfsr", n_eta-1,eta_bins);


   // Start with pT-balance
   //create plots
   TCanvas* a[n_eta-1];
   TString plotname[n_eta-1];
   TF1 *pol1[n_eta-1];
   for (int j=0; j<n_eta-1; j++){ //n_eta-1
     plotname[j]="dijet_kfsr_diffPt_eta_"+eta_range[j]+"_"+eta_range[j+1];
     a[j] = new TCanvas(plotname[j], plotname[j], 800,700);
     m_gStyle->SetOptTitle(0);
     pTgraph_rel_r[j]->Draw("AP");
     pTgraph_rel_r[j]->GetYaxis()->SetRangeUser(0.92,1.08);
     pTgraph_rel_r[j]->GetXaxis()->SetRangeUser(0,alpha_bins[n_alpha-1]+0.01);
     pTgraph_rel_r[j]->GetYaxis()->SetTitle("(R_{MC}/R_{DATA})");
     pTgraph_rel_r[j]->GetXaxis()->SetTitle("cut on #alpha");

     pol1[j] = new TF1("pol1","pol1",0.09,0.36);

     pol1[j]->SetParameters(0,0);
     if(j == 13){
       pol1[j]->SetParameters(0.985,0.05);
     }
     pTgraph_rel_r[j]->Fit(pol1[j],"R");
     line->SetLineStyle(2);
     line->Draw("SAME");

     // fill the output.dat file
     if (fp_rel_r!=NULL) {
       Float_t value = pol1[j]->GetParameter(0);
       Float_t uncert = pol1[j]->GetParError(0);
       fprintf(fp_rel_r, "%f %f\n",value,uncert);
     }
     plotkfsr->SetBinContent(j+1,pol1[j]->GetParameter(0));
     plotkfsr->SetBinError(j+1,pol1[j]->GetParError(0));
     kFSR_DiJet->SetBinContent(j+1,pol1[j]->GetParameter(0));
     kFSR_DiJet->SetBinError(j+1,pol1[j]->GetParError(0));

     leg1->SetHeader("p_{T} balance, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
     leg1->Draw();

     TLatex *tex = new TLatex();
     tex->SetNDC();
     tex->SetTextSize(0.045); 
     tex->DrawLatex(0.64,0.91,CorrectionObject::_lumitag);


     TString chi2_loglin = "#chi^{2}/n.d.f = ";
     chi2_loglin += trunc(pol1[j]->GetChisquare());
     chi2_loglin +="/";
     chi2_loglin +=trunc(pol1[j]->GetNDF());

     TLatex *tex2 = new TLatex();
     tex2->SetNDC();
     tex2->SetTextSize(0.035); 
     tex2->DrawLatex(0.64,0.35,chi2_loglin);
     cout << "Printing kFSR plots to " << CorrectionObject::_outpath+"plots/kFSR_Pt_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf" << endl;
     a[j]->Print(CorrectionObject::_outpath+"plots/kFSR_Pt_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");


     //delete stuff

     delete tex2;
     delete tex;
   }
   
   cout << endl << endl << "finished all fits for rel" << endl  << endl;
   fclose(fp_rel_r);
   cout << "closed some file, now opening output file" << endl;
   
   // create output file including the kFSR plot
   TFile* outputfile_rel_r;
   cout << "Creating output-rootfile:" << CorrectionObject::_outpath+"Histo_KFSR_DiJet_"+CorrectionObject::_generator_tag+"_L1.root" << endl;
   outputfile_rel_r = new TFile(CorrectionObject::_outpath+"Histo_KFSR_DiJet_"+CorrectionObject::_generator_tag+"_L1.root","RECREATE");
   
   cout << "now writing kFSR values" << endl;
   kFSR_DiJet->Write();
   outputfile_rel_r->Write();
   cout << "closing output file" << endl;
   outputfile_rel_r->Close();
   cout << "closed output file" << endl;



   // And now MPF results
   //create plots
   TCanvas* b[n_eta-1];
   for (int j=0; j<n_eta-1; j++){
     plotname[j]="mpf_kfsr_diffPt_eta_"+eta_range[j]+"_"+eta_range[j+1];
     b[j] = new TCanvas(plotname[j], plotname[j], 800,700);
     m_gStyle->SetOptTitle(0);

     pTgraph_mpf_r[j]->Draw("AP");
     pTgraph_mpf_r[j]->GetYaxis()->SetRangeUser(0.92,1.08);
     pTgraph_mpf_r[j]->GetXaxis()->SetRangeUser(0,alpha_bins[n_alpha-1]+0.01);
     pTgraph_mpf_r[j]->GetYaxis()->SetTitle("(R_{MC}/R_{DATA})");
     pTgraph_mpf_r[j]->GetXaxis()->SetTitle("cut on #alpha");
     pol1[j] = new TF1("pol1","pol1",0.09,0.36);  //TEST AK4
     pol1[j]->SetParameters(0,0);
     pTgraph_mpf_r[j]->Fit(pol1[j],"R");
     line->SetLineStyle(2);
     line->Draw("SAME");


     // fill the output.dat file
     if (fp_mpf_r!=NULL) {
       Float_t value = pol1[j]->GetParameter(0);
       Float_t uncert = pol1[j]->GetParError(0);
       fprintf(fp_mpf_r, "%f %f\n",value,uncert);
     }
     plotkfsr->SetBinContent(j+1,pol1[j]->GetParameter(0));
     plotkfsr->SetBinError(j+1,pol1[j]->GetParError(0));
     kFSR_MPF->SetBinContent(j+1,pol1[j]->GetParameter(0));
     kFSR_MPF->SetBinError(j+1,pol1[j]->GetParError(0));

     leg1->SetHeader("MPF, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
     leg1->Draw();

     TLatex *tex = new TLatex();
     tex->SetNDC();
     tex->SetTextSize(0.045); 
     tex->DrawLatex(0.64,0.91,CorrectionObject::_lumitag);

     TString chi2_loglin = "#chi^{2}/n.d.f = ";
     chi2_loglin += trunc(pol1[j]->GetChisquare());
     chi2_loglin +="/";
     chi2_loglin +=trunc(pol1[j]->GetNDF());

     TLatex *tex2 = new TLatex();
     tex2->SetNDC();
     tex2->SetTextSize(0.035); 
     tex2->DrawLatex(0.64,0.35,chi2_loglin);

     //save the plots
     cout << "Saving the MPF plots to " << CorrectionObject::_outpath + "plots/kFSR_MPF_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf" << endl;
     b[j]->Print(CorrectionObject::_outpath + "plots/kFSR_MPF_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");

     //delete stuff
     delete tex2;
     delete tex;
   }
   fclose(fp_mpf_r);


   // create output file including the kFSR plot
   TFile* outputfile_mpf_r;
   cout << "Creating output-rootfile: " << CorrectionObject::_outpath+"Histo_KFSR_MPF_"+CorrectionObject::_generator_tag+"_L1.root" << endl;
   outputfile_mpf_r = new TFile(CorrectionObject::_outpath+"Histo_KFSR_MPF_"+CorrectionObject::_generator_tag+"_L1.root","RECREATE");
   kFSR_MPF->Write();
   outputfile_mpf_r->Write();
   outputfile_mpf_r->Close();




   cout << "+++++++++++++++++ Finished kFSR() +++++++++++++++++++" << endl;

   //delete everything
   




   for (int j=0; j<n_eta-1; j++){ //n_eta-1
     delete pol1[j];
   }

   delete outputfile_mpf_r;
   delete outputfile_rel_r;
   
   


   delete plotkfsr;
   delete kFSR_DiJet;
   delete kFSR_MPF;
   delete line;
   

   for(int j=0; j<n_eta-1; j++){
     for(int k=0; k<n_pt-1; k++){
       delete graph_rel_r[k][j];
       delete graph_mpf_r[k][j];
     }
   }
   for(int j=0; j<n_eta-1; j++){
     delete pTgraph_rel_r[j];
     delete pTgraph_mpf_r[j];
   }
   


   delete leg1;
   for(int i=0; i<n_alpha; i++){
     for(int j=0; j<n_eta-1; j++){
       for(int k=0; k<n_pt-1; k++){
	 delete hdata_rel_r[k][j][i];
	 delete hdata_mpf_r[k][j][i];
	 delete hmc_rel_r[k][j][i];
	 delete hmc_mpf_r[k][j][i];
       }
     }
   }
   delete  m_gStyle;    
   cout << "++++++++++++ Deleted everything in kFSR(), exiting ++++++++++++++" << endl;
}
