#include "UHH2/BaconJets/include/selection.h"

#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/PrimaryVertex.h"
// #include "UHH2/BaconTrans/baconheaders/TJet.hh"
//#include "UHH2/BaconTrans/baconheaders/TVertex.hh"
#include "UHH2/BaconJets/include/constants.h"

#include "TVector2.h"

namespace uhh2bacon {

Selection::Selection(uhh2::Context & ctx) :
    context(ctx),
    event(0)
{
  // auto jetCollection = ctx.get("jetCollection");
  // h_jets = ctx.declare_event_input<TClonesArray>(jetCollection);
  // //    h_jets = context.declare_event_input<TClonesArray>("AK4PFCHS");
  // // h_jets = context.declare_event_input<TClonesArray>("AK4PFPUPPI");
  // h_eventInfo = context.declare_event_input<baconhep::TEventInfo>("Info");
  // h_pv = context.declare_event_input<TClonesArray>("PV");

  tt_gen_pthat = ctx.declare_event_output<float>("gen_pthat");
  tt_gen_weight = ctx.declare_event_output<float>("gen_weight");
  tt_jet1_pt = ctx.declare_event_output<float>("jet1_pt");
  tt_jet2_pt = ctx.declare_event_output<float>("jet2_pt");
  tt_jet3_pt = ctx.declare_event_output<float>("jet3_pt");
  tt_jet1_ptRaw = ctx.declare_event_output<float>("jet1_ptRaw");
  tt_jet2_ptRaw = ctx.declare_event_output<float>("jet2_ptRaw");
  tt_jet3_ptRaw = ctx.declare_event_output<float>("jet3_ptRaw");
  tt_nvertices = ctx.declare_event_output<int>("nvertices");
  tt_probejet_eta = ctx.declare_event_output<float>("probejet_eta");
  tt_probejet_phi = ctx.declare_event_output<float>("probejet_phi");
  tt_probejet_pt = ctx.declare_event_output<float>("probejet_pt");
  tt_probejet_ptRaw = ctx.declare_event_output<float>("probejet_ptRaw");
  tt_barreljet_eta = ctx.declare_event_output<float>("barreljet_eta");
  tt_barreljet_phi = ctx.declare_event_output<float>("barreljet_phi");
  tt_barreljet_pt = ctx.declare_event_output<float>("barreljet_pt");
  tt_barreljet_ptRaw = ctx.declare_event_output<float>("barreljet_ptRaw");
  tt_pt_ave = ctx.declare_event_output<float>("pt_ave");
  tt_alpha = ctx.declare_event_output<float>("alpha");
  tt_rel_r = ctx.declare_event_output<float>("rel_r");
  tt_mpf_r = ctx.declare_event_output<float>("mpf_r");
  tt_asymmetry = ctx.declare_event_output<float>("asymmetry");
  tt_nPU = ctx.declare_event_output<int>("nPU");

}

void Selection::SetEvent(uhh2::Event& evt)
{
   event = &evt;
   assert(event);
}

  bool Selection::PUpthat(uhh2::Event& evt)
  {
    assert(event);
  
   double  pt_hat = event->genInfo->binningValues()[0];
   double  PU_pt_hat = event->genInfo->PU_pT_hat_max();
  
    double Ratio = PU_pt_hat/pt_hat;

    if(Ratio < 1) return true;

    return false;
  }



bool Selection::PtMC(uhh2::Event& evt)
{
  assert(event);
  //  std::cout<<"evt.get(tt_pt_ave) = "<<evt.get(tt_pt_ave)<<" s_Pt_Ave40_cut = "<<s_Pt_Ave40_cut<<std::endl;
  if (evt.get(tt_pt_ave) < s_Pt_AveMC_cut) 
    return false;
  return true;
}


bool Selection::DiJet()
{
    assert(event);
    const int njets = event->jets->size();
    if (njets>=2) return true;

    return false;
}

bool Selection::DiJetAdvanced(uhh2::Event& evt)
{
    assert(event);

    const int njets = event->jets->size();
    if (njets < 2) return false;
    Jet* jet1 = &event->jets->at(0);// leading jet
    Jet* jet2 = &event->jets->at(1);// sub-leading jet

    // at least one barrel jet
    if((fabs(jet1->eta()) >= s_eta_barr) && (fabs(jet2->eta()) >= s_eta_barr)) return false; 

    // delta phi > 2.7
    double deltaPhi = std::abs(TVector2::Phi_mpi_pi(jet1->phi() - jet2->phi()));
    if (deltaPhi < s_delta_phi) return false;

    // |asymm| < 0.7
    if (fabs((event->get(tt_jet2_pt) - event->get(tt_jet1_pt)) / (event->get(tt_jet2_pt) + event->get(tt_jet1_pt))) > s_asymm) return false;

    //(pTgen1 < 1.5*pThat || pTreco1 < 1.5* pTgen1)
    if(!event->isRealData){
      if(event->genjets->size() < 1) return false;
      if(!(event->genjets->at(0).pt() < 1.5*event->genInfo->binningValues()[0] || event->jets->at(0).pt() < 1.5*event->genjets->at(0).pt())) return false;
    }
    return true;
}

  int Selection::goodPVertex()
  {
    assert(event);
    Int_t nvertices = event->pvs->size();
    // require in the event that there is at least one reconstructed vertex
    if(nvertices<=0) return 0;//false;
    float nPrVer = 0;
    // pick the first (i.e. highest sum pt) verte
    for (int i=0; i<nvertices; i++){
      PrimaryVertex* vertices = &event->pvs->at(i);
      // require that the vertex meets certain criteria
      //if(vertices->nTracksFit)
      //std::cout<<" vertices->nTracksFit = "<<vertices->nTracksFit<<" "<<fabs(vertices->z)<<" cut at "<<s_n_PvTracks<<std::endl;
      //        if((vertices->nTracksFit > s_n_PvTracks) && (fabs(vertices->z) < s_n_Pv_z) && (fabs(vertices->y) < s_n_Pv_xy) && (fabs(vertices->x) < s_n_Pv_xy)){
      //SHOULD BE LIKE
      //cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2.")
      // std::cout<<fabs(vertices->z())<<" "<<fabs(vertices->rho())<<" "<<vertices->ndof()<<" "<<vertices->chi2()<<" "<<std::endl;
      if((fabs(vertices->z()) < 24.) && (fabs(vertices->rho()) < 2.) && (vertices->ndof() >= 4) 
	 && (vertices->chi2()) > 0){
	nPrVer++;
      }
    }
    return nPrVer;
    //    std::cout<<" nPrVer = "<<nPrVer<<" all vtxs = "<<nvertices<<std::endl;
    //    event->set(tt_nGoodvertices,nPrVer); 
    // goodVtx = nPrVer; 
    // if(nPrVer<=0) return false;
    // else
    //   return true;
  }

// bool Selection::goodPVertex()
// {
//     assert(event);

//     const TClonesArray & pvs = event->get(h_pv);

//     Int_t nvertices = pvs.GetEntries();
//     // require in the event that there is at least one reconstructed vertex
//     if(nvertices<=0) return false;
//     float nPrVer = 0;
//     // pick the first (i.e. highest sum pt) verte
//     for (int i=0; i<nvertices; i++){
//         baconhep::TVertex* vertices = (baconhep::TVertex*)pvs[i];
//         // require that the vertex meets certain criteria

//         if((vertices->nTracksFit > s_n_PvTracks) && (fabs(vertices->z) < s_n_Pv_z) && (fabs(vertices->y) < s_n_Pv_xy) && (fabs(vertices->x) < s_n_Pv_xy)){
//             nPrVer++;
//         }
//     }

//  return true;
// }





// bool Selection::FullSelection()
// {
//     return DiJet()&&goodPVertex();

// }

Selection::~Selection()
{
}

}
