//#include <iostream>
// ROOT includes
#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TString.h>
#include <TList.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TMinuit.h>
#include <TLorentzVector.h>
#include <TProfile.h>

#include "TClonesArray.h"

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/bacondataformats/interface/TEventInfo.hh"
//#include "UHH2/bacondataformats/interface/BaconAnaDefs.hh"
using namespace std;

namespace uhh2bacon {

class PileupData {

    private:
    uhh2::Context& context;
    uhh2::Event* event;

    uhh2::Event::Handle<TClonesArray> h_jets;
    uhh2::Event::Handle<baconhep::TEventInfo> h_eventInfo;

    bool        fPuReweighting;

    public:
    PileupData(uhh2::Context & ctx);
    ~PileupData();



    void        SetEvent(uhh2::Event& evt);
    float       getDataPU(int Run, int Ls);

};

}