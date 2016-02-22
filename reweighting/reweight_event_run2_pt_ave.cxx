// a small script to reweight nPu distribution

// ROOT includes
#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>
using namespace std;

// weighting factors will be stored;

TString     nPu_reweighting_dbase = "reweight_event";


// main function
int reweight_event_run2_pt_ave() {


  //    TString     sTriggerThreshold[] = {"40","60", "80", "140", "200", "260", "320", "400","500"};
  TString     sTriggerThreshold[] = {"40","60", "60","60","200", "260", "320", "400","500"};
    //    Float_t     fTriggerThreshold[10] = {55,76,93,172,232,300,366,453,558,600};
    Float_t     fTriggerThreshold[] = {56,78,93,172,232,300,366,453,562,600};

    float scale_factor = 0;
    float mc_integral = 0;
    float data_integral = 0;
    float data_all_integral = 0;
    int bmin = 1;
    int bmax = 1;
    int size = sizeof(sTriggerThreshold)/sizeof(double);

    TFile *file_DATA = new TFile ("/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_ReweightPU_AllTriggers_TTree/uhh2.AnalysisModuleRunner.DATA.RunD_AK4CHS.root");
    TFile *file_MC = new TFile ("/nfs/dust/cms/user/karavdia/JEC_76X/Fall15_25nsV1_ReweightPU_AllTriggers_TTree/uhh2.AnalysisModuleRunner.MC.QCD_Pt15to7000_pythia8_AK4CHS.root");


//  TFile *file_DATA = new TFile ("/nfs/dust/cms/user/kovalch/sFrame/JEC/V6/uhh2.AnalysisModuleRunner.DATA.DATAdata_1200pt_ave_V6_CDv3Dv4.root");//finite bining
//  TFile *file_MC = new TFile ("/nfs/dust/cms/user/kovalch/sFrame/JEC/V5/uhh2.AnalysisModuleRunner.MC.MC_with_MCandPU69_new_central_smearing_Asympt_V6_1200.root");//finite bining
//     TFile *file_MC = new TFile ("/nfs/dust/cms/user/kovalch/sFrame/JEC/V5/uhh2.AnalysisModuleRunner.MC.MC_weight_pt_ave_pt_hat_1_PUweight_69000_Asympt_V6_pt_ave1200.root");
//     TFile *file_MC = new TFile ("/nfs/dust/cms/user/kovalch/sFrame/JEC/V5/uhh2.AnalysisModuleRunner.MC.MC_with_MCandPU69_scale_up_smearing_Asympt_V6_1200.root");
//     TFile *file_MC = new TFile ("/nfs/dust/cms/user/kovalch/sFrame/JEC/V5/uhh2.AnalysisModuleRunner.MC.MC_with_MCandPU69_scale_down_smearing_Asympt_V6_1200.root");
//     TFile *file_MC = new TFile ("/nfs/dust/cms/user/kovalch/sFrame/JEC/V5/uhh2.AnalysisModuleRunner.MC.MC_with_MCandPU80_central_smearing_Asympt_V6_1200.root");

    TH1F * histo_MC = (TH1F*) file_MC -> Get("Selection/pt_ave");
    TH1F * histo_DATA_all = (TH1F*) file_DATA -> Get("Selection/pt_ave");
    for(int i=0; i < 9; i++){
    //    for(int i=0; i < 7; i++){

        TH1F * histo_DATA = (TH1F*) file_DATA -> Get("Selection/pt_ave_hltDiPFJetAve"+sTriggerThreshold[i]);

        bmin = histo_MC->FindBin(fTriggerThreshold[i]); 
        bmax = histo_MC->FindBin(fTriggerThreshold[i+1]-0.5); 
//         cout<<"bmin: "<<bmin<<endl;
//         cout<<"bmax: "<<bmax<<endl;

        mc_integral = histo_MC->Integral(bmin,bmax); //cout<<"mc INt: "<<mc_integral<<endl;
        data_integral = histo_DATA->Integral();//cout<<"data INt: "<<data_integral<<endl;
        data_all_integral = histo_DATA_all->Integral(bmin,bmax); //cout<<"data_all INt: "<<data_all_integral<<endl;

        scale_factor = data_all_integral / mc_integral; cout<<"scale : "<<scale_factor<<endl;
    }
    histo_DATA_all -> SetLineColor(2);
    histo_DATA_all ->Draw("histo");
    histo_MC -> Draw("same histo");
    // finished successfully
    return 0;
}

