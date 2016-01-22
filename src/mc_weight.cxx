// system headers
#include<iostream>
#include<fstream>
#include <sys/stat.h>

// ROOT headers
#include<TCollection.h>
#include<TMath.h>
#include<TFile.h>
#include<TRandom3.h>
#include<TSystem.h>
#include<TVector3.h>
#include <vector>

#include "UHH2/BaconJets/include/mc_weight.h"
#include "UHH2/BaconJets/include/constants.h"
#include "UHH2/bacondataformats/interface/TJet.hh"
using namespace std;

namespace uhh2bacon {


McWeight::McWeight(uhh2::Context & ctx) :
    context(ctx),
    event(0),
    fPuReweighting_histo(NULL)
{
    h_jets = context.declare_event_input<TClonesArray>("AK4PFCHS");
    h_eventInfo = context.declare_event_input<baconhep::TEventInfo>("Info");

    fPuReweighting_histoname.push_back(hPuReweighting_histo69A);
    fPuReweighting_histoname.push_back(hPuReweighting_histo80A);
    fPuReweighting_histoname.push_back(hPuReweighting_histo69F);
    fPuReweighting_histoname.push_back(hPuReweighting_histo80F);

    TString DATABASE_PATH = "/nfs/dust/cms/user/kovalch/DataPileup/PuWeights";

    TString fPuReweighting_filename[] = {"PUweight_V6_minBiasXsec69000_pileupJSON_151102_newAsymptMCSel.root","PUweight_V6_minBiasXsec80000_pileupJSON_151102_newAsymptMCSel.root", "PUweight_V6_minBiasXsec69000_pileupJSON_151102_newFlatMCSel.root", "PUweight_V6_minBiasXsec80000_pileupJSON_151102_newFlatMCSel.root"};

    for (int j = 0; j < fPuReweighting_histoname.size(); j++) {
        file = new TFile (DATABASE_PATH+"/"+fPuReweighting_filename[j]);
        fPuReweighting_histoname[j] = (TH1F*) file -> Get("histo_substr");

        if (fPuReweighting_histoname[j] == NULL) {
            cout << "was not possible to retrieve a histogram from "<< fPuReweighting_filename[j]<< endl;
            abort();
        } else {
            cout << "successfully got the PU reweighting histogram from "<< fPuReweighting_filename[j]<< endl;
        }
    }
}

void McWeight::SetEvent(uhh2::Event& evt)
{
    event = &evt;
    assert(event);

}

//gets a weighting factor for true-PU reweighting 
float  McWeight::getPuReweighting(TString MC_option, int minBiasXsec) {

    assert(event);
    const baconhep::TEventInfo & info = event->get(h_eventInfo);
    baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);
    assert(eventInfo);

    Float_t   Pu_true = eventInfo->nPU;

    // declare variable to store a weighting factor
    Double_t    weighting_factor = 1;
    unsigned    bin = 0;
    // sanity check: a pointer to the histogram should be valid
    for (int i = 0; i < fPuReweighting_histoname.size(); i++) {
        if (fPuReweighting_histoname[i] == NULL) {
            cout << "was not possible to retrieve a histogram " << endl;
            abort();
        }
    }
    if(MC_option == "Asympt"){
        if(minBiasXsec == 69){
            bin = fPuReweighting_histoname[0]->FindBin(Pu_true);
            weighting_factor = fPuReweighting_histoname[0]->GetBinContent(bin);
        } else if(minBiasXsec == 80){
            bin = fPuReweighting_histoname[1]->FindBin(Pu_true);
            weighting_factor = fPuReweighting_histoname[1]->GetBinContent(bin);
        }
    } else if(MC_option == "Flat"){
        if(minBiasXsec == 69){
            bin = fPuReweighting_histoname[2]->FindBin(Pu_true);
            weighting_factor = fPuReweighting_histoname[2]->GetBinContent(bin);
        } else if(minBiasXsec == 80){
            bin = fPuReweighting_histoname[3]->FindBin(Pu_true);
            weighting_factor = fPuReweighting_histoname[3]->GetBinContent(bin);
        }
    }

    if (weighting_factor!=0) return      weighting_factor;
    else return 1;
}
//gets a weighting factor for event reweighting 

float  McWeight::getEvReweighting(int  direction, TString MC_option, int minBiasXsec, TString TriggerType) {

    assert(event);

    const TClonesArray & js = event->get(h_jets);
    baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
    baconhep::TJet* jet2 = (baconhep::TJet*)js[1];

    int njets = js.GetEntries();
    Double_t    ev_weighting_factor = 1;
    unsigned    ev_bin = 0;
    float scale_factor1[n_pt_bins];
    for (int i = 0; i < n_pt_bins; i++){
        if(MC_option == "Asympt" && minBiasXsec == 69){
            if(direction == 0){
                scale_factor1[i] = scale_factor_centrA[i]; //central
            } else if(direction == 1){
                scale_factor1[i] = scale_factor_upA[i]; //scale up
            } else if(direction == -1){
                scale_factor1[i] = scale_factor_downA[i]; //scale down
            } else if(direction == 99){
                scale_factor1[i] = scale_factor_noJERA[i]; //no scale
            } else if(direction == 11){
                scale_factor1[i] = scale_factor_centrA_PU69_set1[i]; //set1
            } else if(direction == 12){
                scale_factor1[i] = scale_factor_centrA_PU69_set2[i]; //set2
            } else if(direction == 13){
                scale_factor1[i] = scale_factor_centrA_PU69_set3[i]; //set3
            } else if(direction == 14){
                scale_factor1[i] = scale_factor_centrA_PU69_set4[i]; //set4
            } else if(direction == 15){
                scale_factor1[i] = scale_factor_centrA_PU69_set5[i]; //set5
            } else if(direction == 16){
                scale_factor1[i] = scale_factor_centrA_PU69_set6[i]; //set6
            } else if(direction == 17){
                scale_factor1[i] = scale_factor_centrA_PU69_set7[i]; //set7
            } else if(direction == 18){
                scale_factor1[i] = scale_factor_centrA_PU69_set8[i]; //set8
            } else if(direction == 19){
                scale_factor1[i] = scale_factor_centrA_PU69_set9[i]; //set9
            }
        }
        if(direction == 0 && MC_option == "Asympt" && minBiasXsec == 80){
            scale_factor1[i] = scale_factor_centrA_PU80[i]; //central 80
        }
        if(MC_option == "Flat" && minBiasXsec == 69){
            if(direction == 0){
                scale_factor1[i] = scale_factor_centrF[i]; //central
            } else if(direction == 99){
                scale_factor1[i] = scale_factor_noJERF[i]; //no scale
            }
        }
    }
    // njets >= 2
    if (njets>=2) {
    double pt_ave = (jet1->pt + jet2->pt)/2;
        if (TriggerType == "Nominal"){
            // declare variable to store a weighting factor
            if ((pt_ave >= s_Pt_Ave40_cut) && (pt_ave < s_Pt_Ave60_cut)) {
                ev_weighting_factor = scale_factor1[0];
            } else if ((pt_ave >= s_Pt_Ave60_cut) && (pt_ave < s_Pt_Ave80_cut)) {
                ev_weighting_factor = scale_factor1[1];
            } else if ((pt_ave >= s_Pt_Ave80_cut) && (pt_ave < s_Pt_Ave140_cut)) {
                ev_weighting_factor = scale_factor1[2];
            } else if ((pt_ave >= s_Pt_Ave140_cut) && (pt_ave < s_Pt_Ave200_cut)) {
                ev_weighting_factor = scale_factor1[3];
            } if ((pt_ave >= s_Pt_Ave200_cut) && (pt_ave < s_Pt_Ave260_cut)) {
                ev_weighting_factor = scale_factor1[4];
            } else if ((pt_ave >= s_Pt_Ave260_cut) && (pt_ave < s_Pt_Ave320_cut)) {
                ev_weighting_factor = scale_factor1[5];
            } else if ((pt_ave >= s_Pt_Ave320_cut) && (pt_ave < s_Pt_Ave400_cut)) {
                ev_weighting_factor = scale_factor1[6];
            } else if ((pt_ave >= s_Pt_Ave400_cut) && (pt_ave < s_Pt_Ave500_cut)) {
                ev_weighting_factor = scale_factor1[7];
            } else if ((pt_ave >= s_Pt_Ave500_cut)) {
                ev_weighting_factor = scale_factor1[8];
            }
        } else if (TriggerType == "HF"){
            if ((pt_ave >= s_Pt_Ave60HF_cut) && (pt_ave < s_Pt_Ave80HF_cut)) {
                ev_weighting_factor = scale_factor1[0];
            } else if ((pt_ave >= s_Pt_Ave80HF_cut) && (pt_ave < s_Pt_Ave100HF_cut)) {
                ev_weighting_factor = scale_factor1[1];
            } else if ((pt_ave >= s_Pt_Ave100HF_cut) && (pt_ave < s_Pt_Ave160HF_cut)) {
                ev_weighting_factor = scale_factor1[2];
            } else if ((pt_ave >= s_Pt_Ave160HF_cut) && (pt_ave < s_Pt_Ave220HF_cut)) {
                ev_weighting_factor = scale_factor1[3];
            } if ((pt_ave >= s_Pt_Ave220HF_cut) && (pt_ave < s_Pt_Ave300HF_cut)) {
                ev_weighting_factor = scale_factor1[4];
            } else if ((pt_ave >= s_Pt_Ave300HF_cut)) {
                ev_weighting_factor = scale_factor1[5];
            }
        }

   //cout << "for event with pt_ave = "<< pt_ave << " apply scale_factor = "<<ev_weighting_factor<<endl;
    }
    if (ev_weighting_factor!=0) return      ev_weighting_factor;
    else return 1;
}

McWeight::~McWeight()
{
}

}