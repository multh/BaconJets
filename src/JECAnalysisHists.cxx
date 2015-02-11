#include "UHH2/BaconJets/include/JECAnalysisHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/bacondataformats/interface/TJet.hh"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <stdlib.h>

using namespace std;
using namespace uhh2;

JECAnalysisHists::JECAnalysisHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // jets
  book<TH1F>("N_jets", "N_{jets}", 20, -0.5, 19.5);  
  book<TH1F>("pt","p_{T} all jets",100,0,1500);
  book<TH1F>("eta","#eta all jets",100,-5,5);
  book<TH1F>("phi","#phi all jets",50,-M_PI,M_PI);
 
  book<TH1F>("pt_1","p_{T} jet 1",100,0,1500);
  book<TH1F>("pt_2","p_{T} jet 2",100,0,1500);
  book<TH1F>("pt_3","p_{T} jet 3",100,0,1500);

  book<TH1F>("pt_barrel","p_{T} barrel jet",100,0,1500);
  book<TH1F>("pt_probe","p_{T} probe jet",100,0,1500);

  book<TH1F>("eta_1","#eta jet 1",100,-5,5);
  book<TH1F>("eta_2","#eta jet 2",100,-5,5);
  book<TH1F>("eta_3","#eta jet 3",100,-5,5);

  book<TH1F>("asym","asymmetrie jet 1 and jet 2",100,-1,1);
  book<TH1F>("generic_asym","generic asymmetrie jet 1 and jet 2",100,-1,1);

  book<TH1F>("DeltaPhi_Jet1_Jet2", "#Delta#Phi(first jet, second jet)", 70, 0, 7);
  book<TH1F>("pt_rel","p_{T} jet 3 / #overline{P}_{T}", 50, 0, 1);
  book<TH1F>("generic_pt_rel","generic p_{T} jet 3 / #overline{P}_{T}", 50, 0, 1);

  book<TH2F>("ptrel_vs_deltaphi","delta phi vs pt_rel", 50, 0, 1 ,64, 0, 3.20);

  h_jets = ctx.get_handle<TClonesArray>("Jet04");

}


void JECAnalysisHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.

  const TClonesArray & js = event.get(h_jets);

  double weight = event.weight;

  Int_t njets = js.GetEntries();
  hist("N_jets")->Fill(njets, weight);

  for (int i=0; i<njets; i++){
    baconhep::TJet* jets = (baconhep::TJet*)js[i];
    hist("pt")->Fill(jets->pt, weight);
    hist("eta")->Fill(jets->eta, weight);
    hist("phi")->Fill(jets->phi, weight);
  }
  
  baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
  baconhep::TJet* jet2 = (baconhep::TJet*)js[1];
  hist("pt_1")->Fill(jet1->pt, weight);
  hist("eta_1")->Fill(jet1->eta, weight);


  double barreljet = 0.0;
  double probejet = 0.0;
  int numb = rand() % 2 + 1; 
  if(numb==1){
    if(fabs(jet1->eta)<1.3){
      barreljet += jet1->pt;
      probejet += jet2->pt;
      hist("pt_barrel")->Fill(jet1->pt, weight);
      hist("pt_probe")->Fill(jet2->pt, weight);
      hist("asym")->Fill((jet2->pt - jet1->pt) / (jet2->pt + jet1->pt), weight);
    }
  }
  if(numb==2){
    if(fabs(jet2->eta)<1.3){
      barreljet += jet2->pt;
      probejet += jet1->pt;
      hist("pt_barrel")->Fill(jet2->pt, weight);
      hist("pt_probe")->Fill(jet1->pt, weight);
      hist("asym")->Fill((jet1->pt - jet2->pt) / (jet2->pt + jet1->pt), weight);
    }
  }
 

 

  hist("pt_2")->Fill(jet2->pt, weight);
  hist("eta_2")->Fill(jet2->eta, weight);

  if(jet1->eta<1.3){
    hist("generic_asym")->Fill((jet2->pt - jet1->pt) / (jet2->pt + jet1->pt), weight);
  }
  if(jet2->eta<1.3){
    hist("generic_asym")->Fill((jet1->pt - jet2->pt) / (jet2->pt + jet1->pt), weight);
  }

  double deltaPhi = abs(jet1->phi - jet2->phi);
  hist("DeltaPhi_Jet1_Jet2")->Fill(deltaPhi, weight);

  
  baconhep::TJet* jet3 = (baconhep::TJet*)js[2];
  if (njets>2){
  hist("pt_3")->Fill(jet3->pt, weight);
  hist("eta_3")->Fill(jet3->eta, weight);
  hist("pt_rel")->Fill(jet3->pt/(0.5*(barreljet + probejet)),weight);
  hist("generic_pt_rel")->Fill(jet3->pt/(0.5*(jet1->pt + jet2->pt)),weight);
  hist("ptrel_vs_deltaphi")->Fill(jet3->pt/(0.5*(barreljet + probejet)),deltaPhi);
  }

}

JECAnalysisHists::~JECAnalysisHists(){}

