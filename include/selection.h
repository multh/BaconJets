//#include <iostream>

#include "TClonesArray.h"

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"


#include "UHH2/bacondataformats/interface/TEventInfo.hh"
//#include "UHH2/bacondataformats/interface/BaconAnaDefs.hh"

namespace uhh2bacon {

class Selection {

    private:
    uhh2::Context& context;
    uhh2::Event* event;

    uhh2::Event::Handle<TClonesArray> h_jets;
    uhh2::Event::Handle<baconhep::TEventInfo> h_eventInfo;
    uhh2::Event::Handle<TClonesArray> h_pv;
    public:
    Selection(uhh2::Context & ctx);
    ~Selection();

    void SetEvent(uhh2::Event& evt);
    bool Trigger();
    bool DiJet();
    bool DiJetAdvanced();
    bool goodPVertex();
    bool FullSelection();


};

}