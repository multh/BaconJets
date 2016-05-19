// a small script to reweight nPu distribution

// ROOT includes
#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>

// weighting factors will be stored;

TString     nPu_reweighting_dbase = "nPu_reweighting";
TString     sminBiasXsec = "58000";
//TString     sminBiasXsec = "69000";
//TString     sminBiasXsec = "80000";
TString     sMCsample = "Flat";
//TString MCgen = "herwigpp";
TString MCgen = "pythia8";

// main function
int reweight_nPu_run2() {
//     TFile *file_DATA = new TFile ("/nfs/dust/cms/user/kovalch/DataPileup/PuWeights/MyDataPileup_V6_minBiasXsec58000_pileupJSON_151102.root");
//     TFile *file_DATA = new TFile ("/nfs/dust/cms/user/kovalch/DataPileup/PuWeights/MyDataPileup_V6_minBiasXsec69000_pileupJSON_151102.root");
 TFile *file_DATA = new TFile ("/nfs/dust/cms/user/kovalch/DataPileup/PuWeights/MyDataPileup_V6_minBiasXsec"+sminBiasXsec+"_pileupJSON_151102.root");
//     TFile *file_DATA = new TFile ("/nfs/dust/cms/user/kovalch/DataPileup/PuWeights/MyDataPileup_V6_minBiasXsec80000_pileupJSON_151102.root");

    //my files
  //  TFile *file_DATA = new TFile ("/nfs/dust/cms/user/kovalch/DataPileup/PuWeights/MyDataPileup_V6_minBiasXsec58000_pileupJSON_151102.root");
  //  TFile *file_DATA = new TFile ("/afs/desy.de/user/k/karavdia/CMSSW_7_6_3/src/UHH2/BaconJets/reweighting/pileup_76Xcampaign/MyDataPileup_V6_minBiasXsec"+sminBiasXsec+"_pileupJSON_160207.root");
  //  TFile *file_DATA = new TFile ("/nfs/dust/cms/user/kovalch/DataPileup/PuWeights/MyDataPileup_V6_minBiasXsec80000_pileupJSON_151102.root");

    //    TFile *file_MC = new TFile ("/nfs/dust/cms/user/kovalch/sFrame/JEC/V5/uhh2.AnalysisModuleRunner.MC.MC_no_weight_Flat_V6.root");
//     TFile *file_MC = new TFile ("/nfs/dust/cms/user/kovalch/sFrame/JEC/V5/uhh2.AnalysisModuleRunner.MC.MC_no_weight_new_mc_Asympt_V6.root");

//  TFile *file_DATA = new TFile ("/afs/desy.de/user/k/karavdia/CMSSW_7_6_3/src/UHH2/BaconJets/reweighting/pileup_76Xcampaign/MyDataPileup_V6_minBiasXsec58000_pileupJSON_160207.root");
// TFile *file_DATA = new TFile ("/afs/desy.de/user/k/karavdia/CMSSW_7_6_3/src/UHH2/BaconJets/reweighting/pileup_76Xcampaign/MyDataPileup_V6_minBiasXsec69000_pileupJSON_160207.root");
//  TFile *file_DATA = new TFile ("/afs/desy.de/user/k/karavdia/CMSSW_7_6_3/src/UHH2/BaconJets/reweighting/pileup_76Xcampaign/MyDataPileup_V6_minBiasXsec80000_pileupJSON_160207.root");
  // TFile *file_MC = new TFile ("/nfs/dust/cms/user/karavdia/JEC_76X/MC_QCD_Pt-15to7000_Flat_pythia8_Fall15_25nsV1_MC_noReweight/uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_pythia8.root");
  TFile *file_MC = new TFile ("/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_noReweight_AllTriggers_TTree/uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_"+MCgen+"_AK4CHS.root");
    TH1F * histo_DATA = (TH1F*) file_DATA -> Get("pileup");
    TH1F * histo_MC = (TH1F*) file_MC -> Get("Selection/nPu");

    // normalize DATA and MC to one
    cout<<"histo_DATA->Integral() = "<<histo_DATA->Integral()<<" histo_MC->Integral() = "<<histo_MC->Integral()<<endl;
    histo_DATA -> Scale(1./histo_DATA->Integral());
    histo_MC -> Scale(1./histo_MC->Integral());


    histo_DATA  -> SetMarkerStyle(20);
    histo_MC    -> SetLineWidth(2);
    histo_MC    -> SetLineColor(4);

    // normalize MC to DATA
//     histo_MC -> Scale(histo_DATA->Integral()/histo_MC->Integral());

    histo_DATA -> Draw();
    histo_MC -> Draw("sameHIST");
    histo_DATA -> Draw("sameE1");
    gROOT -> GetListOfCanvases()->Print("01_before_reweighting_run2_Selection.png");

    // histogram for data/mc ratio
    // to cross-check the weight
    //    const int npoints = 60;
    const int npoints = 50;
    std::vector<double> result(npoints);
    double s = 0.0;
    for(int npu=0; npu<npoints; ++npu) {
        double npu_estimated = histo_DATA->GetBinContent(histo_DATA->GetXaxis()->FindBin(npu));
        double npu_probs_Summer2012 = histo_MC->GetBinContent(histo_MC->GetXaxis()->FindBin(npu));
	cout<<"npu_estimated = "<<npu_estimated<<" npu_probs_Summer2012 = "<<npu_probs_Summer2012<<endl;
        result[npu] = npu_estimated / npu_probs_Summer2012;
    //    cout <<" q "<<result[npu]<<endl;
        s += npu_estimated;
    }
    // normalize weights such that the total sum of weights over the whole sample is 1.0, i.e., sum_i  result[i] * npu_probs_Summer2012[i] should be 1.0 (!)
    for(int n=0; n<npoints; ++n) {
        result[n] /= s;
        //cout << "results = "<< result[n] << endl;
    }

    TH1F * histo_substr = (TH1F*) histo_DATA -> Clone("histo_substr");
    // histo_substr->Print();
    // histo_MC->Print();
    histo_substr -> Scale(1./histo_substr->Integral()); //TEST
    histo_MC -> Scale(1./histo_MC->Integral()); //TEST
    histo_substr -> Divide(histo_MC);
    histo_substr -> Scale(1/s);
    // to cross-check the weight
//     for(int i=0; i<60; ++i) {
//         cout << "histo = "<< histo_substr->GetBinContent(histo_substr->GetXaxis()->FindBin(i)) <<endl;
//     }
    histo_substr -> Draw();
    gROOT -> GetListOfCanvases()->Print("02_weighting_function_run2_Selection.png");


    TH1F *  histo_nPu_reweighted = (TH1F*) histo_MC -> Clone("histo_nPu_reweighted");
    histo_nPu_reweighted -> Multiply(histo_substr);
    histo_nPu_reweighted -> Draw();
    gROOT -> GetListOfCanvases()->Print("03_nPu_reweighted_run2_Selection.png");


    histo_nPu_reweighted -> SetLineColor(kRed);
    histo_DATA -> Draw("E1");
    histo_nPu_reweighted -> Draw("sameHIST");
    gROOT -> GetListOfCanvases()->Print("04_after_reweighting_run2_Selection.png");

    // now store the reweighting function in order to not bothering with parametrization, doing simple
    // bin-by-bin correction, hence storing simply a weighting factor for each bin 

    // TH1F variant

    //    TFile   * out_file = new TFile ("PUweight_V6_minBiasXsec"+sminBiasXsec+"_pileupJSON_151102_new"+sMCsample+"MCSel.root", "recreate");
    TFile   * out_file = new TFile ("PUweight_RunII_76X_V1_minBiasXsec"+sminBiasXsec+"_pileupJSON_"+sMCsample+"MCSel_"+MCgen+"_.root", "recreate");
    out_file -> cd();
    histo_substr -> SetDirectory(gDirectory);
    histo_substr -> Write();
    out_file -> Close();

    // finished successfully
    return 0;
}

