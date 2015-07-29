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

    TString DATABASE_PATH = "/nfs/dust/cms/user/kovalch/DataPileup/PuWeights_run2";

    file = new TFile (DATABASE_PATH+"/lowPu_golden_json.root");
    fPuReweighting_histo = (TH1F*) file -> Get("histo_substr");

}

void McWeight::SetEvent(uhh2::Event& evt)
{
    event = &evt;
    assert(event);

}

//gets a weighting factor for true-PU reweighting 
float  McWeight::getPuReweighting() {

    assert(event);
    const baconhep::TEventInfo & info = event->get(h_eventInfo);
    baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);
    assert(eventInfo);

    Float_t   Pu_true = eventInfo->nPU;
    // to check which triggers were fired

    // declare variable to store a weighting factor
    Double_t    weighting_factor = 1;
    unsigned    bin = 0;

    bin = fPuReweighting_histo -> FindBin(Pu_true);
    weighting_factor = fPuReweighting_histo -> GetBinContent(bin);

    // finished successfully, return the weighting factor
//       cout <<" ev nu: "<<eventInfo->evtNum << " nPU = "<< Pu_true<< " weighting_factor = "<< weighting_factor<< " bin # = "<< bin <<endl;

    if (weighting_factor!=0) return      weighting_factor;
    else return 1;
}
//gets a weighting factor for event reweighting 

float  McWeight::getEvReweighting() {

    assert(event);

    const TClonesArray & js = event->get(h_jets);
    baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
    baconhep::TJet* jet2 = (baconhep::TJet*)js[1];

    int njets = js.GetEntries();
    Double_t    ev_weighting_factor = 1;
    unsigned    ev_bin = 0;
    // njets >= 2
    if (njets>=2) {
    double pt_ave = (jet1->pt + jet2->pt)/2;

    // declare variable to store a weighting factor
    if ((pt_ave >= s_Pt_Ave40_cut) && (pt_ave < s_Pt_Ave80_cut)) {
        ev_weighting_factor = scale_factor1[0]*scale_factor2[0];

    } else if ((pt_ave >= s_Pt_Ave80_cut) && (pt_ave < s_Pt_Ave140_cut)) {
        ev_weighting_factor = scale_factor1[1]*scale_factor2[1];

    } else if ((pt_ave >= s_Pt_Ave140_cut) && (pt_ave < s_Pt_Ave200_cut)) {
        ev_weighting_factor = scale_factor1[2]*scale_factor2[2];

    } if ((pt_ave >= s_Pt_Ave200_cut) && (pt_ave < s_Pt_Ave260_cut)) {
        ev_weighting_factor = scale_factor1[3]*scale_factor2[3];

    } else if ((pt_ave >= s_Pt_Ave260_cut) && (pt_ave < s_Pt_Ave320_cut)) {
        ev_weighting_factor = scale_factor1[4]*scale_factor2[4];

    } else if ((pt_ave >= s_Pt_Ave320_cut) && (pt_ave < s_Pt_Ave400_cut)) {
        ev_weighting_factor = scale_factor1[5]*scale_factor2[5];

    } else if ((pt_ave >= s_Pt_Ave400_cut)) {
        ev_weighting_factor = scale_factor1[6]*scale_factor2[6];

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