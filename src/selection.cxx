#include "UHH2/BaconJets/include/selection.h"

#include "UHH2/bacondataformats/interface/TJet.hh"
#include "UHH2/bacondataformats/interface/TVertex.hh"
#include "UHH2/BaconJets/include/constants.h"
namespace uhh2bacon {

Selection::Selection(uhh2::Context & ctx) :
    context(ctx),
    event(0)
{
  h_jets = context.declare_event_input<TClonesArray>("nt_AK4PFCluster");
  h_eventInfo = context.declare_event_input<baconhep::TEventInfo>("Info");
  h_pv = context.declare_event_input<TClonesArray>("PV");

//   h_jets = context.declare_event_autput<TClonesArray>("Jet05");
//   h_eventInfo = context.declare_event_autput<baconhep::TEventInfo>("Info");
}

void Selection::SetEvent(uhh2::Event& evt)
{
   event = &evt;
   assert(event);
}


bool Selection::Trigger()
{
    assert(event);

    const TClonesArray & js = event->get(h_jets);
    const baconhep::TEventInfo & info = event->get(h_eventInfo);
//   const baconhep::TJet * jet = dynamic_cast<const baconhep::TJet*>(js[0]);
//   assert(jet);

    baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);
    assert(eventInfo);

    baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
     baconhep::TJet* jet2 = (baconhep::TJet*)js[1];

    double avePt = (jet1->pt + jet2->pt)/2;

    bool trigger40fired = false;
    bool trigger80fired = false;
    bool trigger140fired = false;
    bool trigger200fired = false;
    bool trigger260fired = false;
    bool trigger320fired = false;
    bool trigger400fired = false;


    if(eventInfo->triggerBits[0]==1)  trigger40fired = true;
    if(eventInfo->triggerBits[1]==1)  trigger80fired = true;
    if(eventInfo->triggerBits[2]==1)  trigger140fired = true;
    if(eventInfo->triggerBits[3]==1)  trigger200fired = true;
    if(eventInfo->triggerBits[4]==1)  trigger260fired = true;
    if(eventInfo->triggerBits[5]==1)  trigger320fired = true;
    if(eventInfo->triggerBits[6]==1)  trigger400fired = true;

    if (avePt < s_Pt_Ave40_cut) return false;

    if (avePt >= s_Pt_Ave40_cut  && avePt < s_Pt_Ave80_cut  && trigger40fired) return true;
    if (avePt >= s_Pt_Ave80_cut  && avePt < s_Pt_Ave140_cut && trigger80fired) return true;
    if (avePt >= s_Pt_Ave140_cut && avePt < s_Pt_Ave200_cut && trigger140fired) return true;
    if (avePt >= s_Pt_Ave200_cut && avePt < s_Pt_Ave260_cut && trigger200fired) return true;
    if (avePt >= s_Pt_Ave260_cut && avePt < s_Pt_Ave320_cut && trigger260fired) return true;
    if (avePt >= s_Pt_Ave320_cut && avePt < s_Pt_Ave400_cut && trigger320fired) return true;
    if (avePt >= s_Pt_Ave400_cut && trigger400fired) return true;

 return false;
}

bool Selection::DiJet()
{
    assert(event);

    const TClonesArray & js = event->get(h_jets);
    const baconhep::TEventInfo & info = event->get(h_eventInfo);
    baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);
    assert(eventInfo);

    int njets = js.GetEntries();

    // njets >= 2
    if (njets>=2) return true;

 return false;
}
bool Selection::DiJetAdvanced()
{
    assert(event);

    const TClonesArray & js = event->get(h_jets);

//   const baconhep::TJet * jet = dynamic_cast<const baconhep::TJet*>(js[0]);
//   assert(jet);

    int njets = js.GetEntries();

     // njets >= 2
    if (njets < 2) return false;

    baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
    baconhep::TJet* jet2 = (baconhep::TJet*)js[1];

    // at least one barrel jet
    if((fabs(jet1->eta) >= s_eta_barr) && (fabs(jet2->eta) >= s_eta_barr)) return false; 

    // delta phi < 2.7
    double deltaPhi = std::min(std::abs(double(jet1->phi) - double(jet2->phi)),2*M_PI-std::abs(jet2->phi - jet1->phi));
    if (deltaPhi < s_delta_phi) return false;

  // |asymm| < 0.7
  if (fabs((jet2->pt - jet1->pt) / (jet2->pt + jet1->pt)) > s_asymm) return false;

  // p_t,rel < 0.2
  //if (njets>2){
  //baconhep::TJet* jet3 = (baconhep::TJet*)js[2];
  // if ((2*(jet3->pt))/(jet1->pt + jet2->pt) > s_pt_rel) return false;
  //}

 return true;
}

bool Selection::AlphaCut()
{
  assert(event);
  
  const TClonesArray & js = event->get(h_jets);
  int njets = js.GetEntries();
  if (njets < 2) return false;
  
  baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
  baconhep::TJet* jet2 = (baconhep::TJet*)js[1];
  
  // p_t,rel < 0.2
  if (njets>2){
    baconhep::TJet* jet3 = (baconhep::TJet*)js[2];
    if ((2*(jet3->pt))/(jet1->pt + jet2->pt) > s_pt_rel) return false;
  }

 return true;
}



bool Selection::goodPVertex()
{
    assert(event);

    const TClonesArray & pvs = event->get(h_pv);

    Int_t nvertices = pvs.GetEntries();
    // require in the event that there is at least one reconstructed vertex
    if(nvertices<=0) return false;
    float nPrVer = 0;
    // pick the first (i.e. highest sum pt) verte
    for (int i=0; i<nvertices; i++){
        baconhep::TVertex* vertices = (baconhep::TVertex*)pvs[i];
        // require that the vertex meets certain criteria

        if(/*(vertices->nTracksFit > s_n_PvTracks) && */(fabs(vertices->z) < s_n_Pv_z) && (fabs(vertices->y) < s_n_Pv_xy) && (fabs(vertices->x) < s_n_Pv_xy)){
            nPrVer++;
        }
    }

 return true;
}

bool Selection::jetIds(float csv_threshold)
{
    assert(event);
    const TClonesArray & js = event->get(h_jets);
    Int_t njets = js.GetEntries();

    for (int i=0; i<njets; i++){
        baconhep::TJet* jets = (baconhep::TJet*)js[i];
        //if (jets->csv < csv_threshold) return false;
    }
 return true;
}



bool Selection::FullSelection()
{
    return Trigger()&&DiJet()&&DiJetAdvanced()&&goodPVertex();

}

Selection::~Selection()
{
}

}
