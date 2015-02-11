#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "../include/JECAnalysisHists.h"

#include "UHH2/bacondataformats/interface/TJet.hh"

#include "TClonesArray.h"

using namespace std;
using namespace uhh2;

namespace uhh2BaconJets {

class DijetSel: public AnalysisModule {
public:
    
    explicit DijetSel(Context & ctx);
    virtual bool process(Event & event) override;

private:
    
  Event::Handle<TClonesArray> h_jets;

  std::unique_ptr<Hists> h_nocuts, h_sel;

};


DijetSel::DijetSel(Context & ctx){
  h_jets = ctx.declare_event_input<TClonesArray>("Jet04");
  h_nocuts.reset(new JECAnalysisHists(ctx,"noCuts"));
  h_sel.reset(new JECAnalysisHists(ctx,"Selection"));
}


bool DijetSel::process(Event & event) {
  const TClonesArray & js = event.get(h_jets);
  //cout << "read n_pv = " << js.GetEntries() << endl;
  /*if(js.GetEntries() > 0){
    // note: pvs[0] returns a "TObject *"; have to cast to "TVertex *":
    const baconhep::TJet * vtx = dynamic_cast<const baconhep::TJet*>(js[0]);
    assert(vtx);
    //cout << "jet pt: " << vtx->pt << ", " << vtx->eta << ", " << vtx->phi << endl;
    }*/
  const baconhep::TJet * jet = dynamic_cast<const baconhep::TJet*>(js[0]);
  assert(jet);
  

  int njets = js.GetEntries();

  baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
  baconhep::TJet* jet2 = (baconhep::TJet*)js[1];
  baconhep::TJet* jet3 = (baconhep::TJet*)js[2];

  double deltaPhi = abs(jet1->phi - jet2->phi);
  double barreljet = 0.0;
  double probejet = 0.0;
  int numb = rand() % 2 + 1; 
  if(numb==1){
    if(fabs(jet1->eta)<1.3){
      barreljet += jet1->pt;
      probejet += jet2->pt;
    }
  }
  if(numb==2){
    if(fabs(jet2->eta)<1.3){
      barreljet += jet2->pt;
      probejet += jet1->pt;
    }
  }



  // njets >= 2
  if (njets<2) return false;
  // eta<1.3 for jet1 or jet2

  
  // asymm < 0.7
  //if ((fabs(probejet - barreljet) / (probejet + barreljet)) > 0.7) return false;
  /*
  if ((fabs(jet2->pt - jet1->pt) / (jet2->pt + jet1->pt)) > 0.7) return false;
  if ((fabs(jet1->pt - jet2->pt) / (jet2->pt + jet1->pt)) > 0.7) return false;
  */

  // fill histos without any cuts
  h_nocuts->fill(event);

  // delta phi < 2.7
  if (deltaPhi < 2.7) return false;


  // p_t,rel < 0.2
  if (njets>2){
    if ((2*(jet3->pt))/(barreljet + probejet) > 0.2) return false;
  }
  /*
  if (njets>2){
    if ((2*(jet3->pt))/(jet1->pt + jet2->pt) > 0.2) return false;
  }
  */

  // fill histos after dijet event selection
  h_sel->fill(event);


  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the ExampleModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(DijetSel)

}
