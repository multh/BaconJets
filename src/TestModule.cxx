#include <iostream>
#include <memory>
#include <stdlib.h>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "../include/JECAnalysisHists.h"

#include "UHH2/BaconJets/include/selection.h"
#include "UHH2/BaconJets/include/jet_corrections.h"
#include "UHH2/BaconJets/include/mc_weight.h"
#include "UHH2/BaconJets/include/constants.h"
#include "UHH2/BaconJets/include/TSetTree.h"

#include "UHH2/BaconTrans/baconheaders/TGenEventInfo.hh"
#include "UHH2/BaconTrans/baconheaders/TJet.hh"
#include "UHH2/BaconTrans/baconheaders/TEventInfo.hh"
#include "UHH2/BaconTrans/baconheaders/BaconAnaDefs.hh"

#include "UHH2/BaconJets/include/pileup_data.h"
#include "UHH2/BaconJets/include/data_corrections.h"
#include "TClonesArray.h"
#include "TString.h"


TTree   *fCurrentTree;
Int_t   Runnr;
Int_t   Eventnr;
TFile   *fCurrentTreeFile;
using namespace std;
using namespace uhh2;
using namespace baconhep;

namespace uhh2bacon {

  class TestModule: public uhh2::AnalysisModule {
  public:
    explicit TestModule(Context & ctx);
    virtual bool process(Event & event) override;
    ~TestModule();

  private:

    Event::Handle<TClonesArray> h_jets  ;
    //    Event::Handle<std::vector<uhh2bacon::TJet>> h_jets;
    Event::Handle<baconhep::TEventInfo> h_eventInfo;
    Event::Handle<baconhep::TGenEventInfo> h_genInfo;
    Event::Handle<TClonesArray> h_pv;
    std::unique_ptr<Hists> h_nocuts, h_sel, h_dijet, h_match;
    std::vector<JECAnalysisHists> h_pt_bins, h_eta_bins, h_pt_bins_a01, h_eta_bins_a01, h_eta_bins_a03, h_eta_bins_a04, h_noalpha_bins,  h_eta_bins_mikko_a10, h_eta_bins_mikko_a15, h_eta_bins_mikko_a20, h_eta_bins_mikko_a30;

    Selection sel;
    //    std::unique_ptr<Selection> sel;
    JetCorrections jetcorr;
    McWeight mcweight;
    bool is_mc;
    bool is_data;
    TSetTree cSetTree;
    uhh2bacon::PileupData  pileupData;

    //Additional vars in Event, specific for dijet
    Event::Handle<float> tt_gen_pthat; Event::Handle<float> tt_gen_weight;
    Event::Handle<float> tt_jet1_pt;     Event::Handle<float> tt_jet2_pt;     Event::Handle<float> tt_jet3_pt;
    Event::Handle<float> tt_jet1_ptRaw;  Event::Handle<float> tt_jet2_ptRaw;  Event::Handle<float> tt_jet3_ptRaw;
    Event::Handle<int> tt_nvertices;
    Event::Handle<float> tt_probejet_eta;  Event::Handle<float> tt_probejet_phi; Event::Handle<float> tt_probejet_pt; Event::Handle<float> tt_probejet_ptRaw;
    Event::Handle<float> tt_barreljet_eta;  Event::Handle<float> tt_barreljet_phi; Event::Handle<float> tt_barreljet_pt; Event::Handle<float> tt_barreljet_ptRaw;
    Event::Handle<float> tt_pt_ave;
    Event::Handle<float> tt_alpha;
    Event::Handle<float> tt_rel_r; Event::Handle<float> tt_mpf_r; Event::Handle<float> tt_asymmetry; Event::Handle<float> tt_nPU;

  };


  TestModule::TestModule(Context & ctx) :
    sel(ctx),
    jetcorr(ctx),
    mcweight(ctx),
    pileupData(ctx)
// // //     datacorr(ctx),
    //  cSetTree()
  {

    cSetTree = TSetTree();

    auto dataset_type = ctx.get("dataset_type");
    is_mc = dataset_type  == "MC";
    is_data = dataset_type  == "DATA";

    h_jets = ctx.declare_event_input<TClonesArray>("AK4PFCHS");
    //    h_jets = ctx.declare_event_input<std::vector<baconhep::TJet>>("AK4PFCHS");
    h_eventInfo = ctx.declare_event_input<baconhep::TEventInfo>("Info");
    if(is_mc){ /// apply for MC only
      h_genInfo = ctx.declare_event_input<baconhep::TGenEventInfo>("GenEvtInfo");
    }
//     h_pv = ctx.get_handle<TClonesArray>("PV");
    h_pv = ctx.declare_event_input<TClonesArray>("PV");


    h_nocuts.reset(new JECAnalysisHists(ctx,"noCuts"));
    h_dijet.reset(new JECAnalysisHists(ctx,"diJet"));
    h_match.reset(new JECAnalysisHists(ctx,"JetMatching"));
    h_sel.reset(new JECAnalysisHists(ctx,"Selection"));


    // int size = sizeof(eta_range)/sizeof(double); // to get size of string object
    // eta_range.size() // to get size of vector

    std::vector<std::string> pt_range_name;
    for( unsigned int i=0; i < pt_range.size(); ++i ){
      char pt_buffer [50];
      sprintf (pt_buffer, "%5.3f", pt_range[i]);
      pt_range_name.push_back(pt_buffer);
    }
    std::vector<std::string> eta_range_name;
    for( unsigned int j=0; j < eta_range.size(); ++j ){
      char eta_buffer [50];
      sprintf (eta_buffer, "%5.3f", eta_range[j]);
      eta_range_name.push_back(eta_buffer);
    }
    // for Mikkos combination root file
    std::vector<std::string> eta_range_mikko_name;
    for( unsigned int j=0; j < eta_range_mikko.size(); ++j ){
      char eta_buffer_mikko [50];
      sprintf (eta_buffer_mikko, "%5.3f", eta_range_mikko[j]);
      eta_range_mikko_name.push_back(eta_buffer_mikko);
    }

    // for( unsigned int k=0; k < alpha_range.size()-1; ++k ){
    //     for( unsigned int j=0; j < eta_range.size()-1; ++j ){
    for( unsigned int i=0; i < pt_range.size()-1; ++i ){
      h_pt_bins.push_back(JECAnalysisHists(ctx,(std::string)("a02/pt_"+pt_range_name[i]+"_"+pt_range_name[i+1])));
      h_pt_bins_a01.push_back(JECAnalysisHists(ctx,(std::string)("a01/pt_"+pt_range_name[i]+"_"+pt_range_name[i+1])));

    }
    //     }
    //   }
    for( unsigned int i=0; i < eta_range.size()-1; ++i ){
      h_eta_bins.push_back(JECAnalysisHists(ctx,(std::string)("a02/eta_"+eta_range_name[i]+"_"+eta_range_name[i+1])));
      h_eta_bins_a01.push_back(JECAnalysisHists(ctx,(std::string)("a01/eta_"+eta_range_name[i]+"_"+eta_range_name[i+1])));
      h_eta_bins_a03.push_back(JECAnalysisHists(ctx,(std::string)("a03/eta_"+eta_range_name[i]+"_"+eta_range_name[i+1])));
      h_eta_bins_a04.push_back(JECAnalysisHists(ctx,(std::string)("a04/eta_"+eta_range_name[i]+"_"+eta_range_name[i+1])));

    }

    for( unsigned int j=0; j < eta_range_mikko.size()-1; ++j ){
      for( unsigned int i=0; i < pt_range.size()-1; ++i ){

        h_eta_bins_mikko_a10.push_back(JECAnalysisHists(ctx,(std::string)("a10/eta_"+eta_range_mikko_name[j]+"_"+eta_range_mikko_name[j+1]+"/pt_"+pt_range_name[i]+"_"+pt_range_name[i+1])));
        h_eta_bins_mikko_a15.push_back(JECAnalysisHists(ctx,(std::string)("a15/eta_"+eta_range_mikko_name[j]+"_"+eta_range_mikko_name[j+1]+"/pt_"+pt_range_name[i]+"_"+pt_range_name[i+1])));
        h_eta_bins_mikko_a20.push_back(JECAnalysisHists(ctx,(std::string)("a20/eta_"+eta_range_mikko_name[j]+"_"+eta_range_mikko_name[j+1]+"/pt_"+pt_range_name[i]+"_"+pt_range_name[i+1])));
        h_eta_bins_mikko_a30.push_back(JECAnalysisHists(ctx,(std::string)("a30/eta_"+eta_range_mikko_name[j]+"_"+eta_range_mikko_name[j+1]+"/pt_"+pt_range_name[i]+"_"+pt_range_name[i+1])));
      }
    }

    // histos for the pT extrapolations
    for( unsigned int j=0; j < eta_range.size()-1; ++j ){
      for( unsigned int i=0; i < pt_range.size()-1; ++i ){
        h_noalpha_bins.push_back(JECAnalysisHists(ctx,(std::string)("/eta_"+eta_range_name[j]+"_"+eta_range_name[j+1]+"/pt_"+pt_range_name[i]+"_"+pt_range_name[i+1])));
      }
    }

    //Additional vars in event
    //

    // gen_pthat = 0; gen_weight = 0;
    // jet1_pt = 0; jet2_pt = 0; jet3_pt = 0; 
    // jet1_ptRaw = 0; jet2_ptRaw = 0; jet3_ptRaw = 0; 
    // nvertices = 0;
    // probejet_eta = 0; probejet_phi = 0; probejet_pt = 0; probejet_ptRaw = 0;
    // barreljet_eta = 0; barreljet_phi = 0; barreljet_pt = 0; barreljet_ptRaw = 0;
    // pt_ave = 0; alpha = 0; rel_r = 0; mpf_r = 0;
    // asymmetry = 0; nPU = 0;
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
    tt_nPU = ctx.declare_event_output<float>("nPU");
  }

  TestModule::~TestModule() {
    cSetTree.general();
    //      cPuData.general();

  }

  bool TestModule::process(Event & event) {

    sel.SetEvent(event);
    jetcorr.SetEvent(event);
    mcweight.SetEvent(event);
    pileupData.SetEvent(event);
//     datacorr.SetEvent(event);
    // float j3L1corr =1.;
    // float j1L1corr =1.;
    // float j2L1corr =1.;
    // const baconhep::TEventInfo & info = event.get(h_eventInfo);
    // baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);

    const baconhep::TEventInfo & eventInfo = event.get(h_eventInfo);
    const TClonesArray & js = event.get(h_jets);
    //    std::vector<baconhep::TJet> js = event.get(h_jets);
    // //! JER smearing
    // if(is_mc){ /// apply for MC only

    //     //! matching from GEN to RECO
    //     if(!jetcorr.JetMatching()) return false;
    //     //! JER smearing
    //     //0 = central; 1 = scale up; other = scale down
    //     if(!jetcorr.JetResolutionSmearer(0)) return false;

    // }

    // baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
    // baconhep::TJet* jet2 = (baconhep::TJet*)js[1];

    baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
    baconhep::TJet* jet2 = (baconhep::TJet*)js[1];
    Int_t njets = js.GetEntries();
    std::cout<<"Number of jets = "<<js.GetEntries()<<std::endl;
    //    std::cout<<"evtNum = "<<eventInfo->evtNum<<std::endl;
    std::cout<<"evtNum = "<<eventInfo.evtNum<<std::endl;
    baconhep::TJet* jet3;

    const TClonesArray & pvs = event.get(h_pv);
    event.set(tt_jet1_pt,jet1->pt);
    event.set(tt_jet2_pt,jet2->pt);
    event.set(tt_jet1_ptRaw,jet1->ptRaw);
    event.set(tt_jet2_ptRaw,jet2->ptRaw);
    // event.jet1_pt = jet1->pt;
    // event.jet2_pt = jet2->pt;
    // event.jet1_ptRaw = jet1->ptRaw;
    // event.jet2_ptRaw = jet2->ptRaw;
    std::cout<<"event.jet1_ptRaw = "<<event.get(tt_jet1_ptRaw)<<" event.jet1_pt = "<<event.get(tt_jet1_pt)<<std::endl;
    std::cout<<"event.jet2_ptRaw = "<<event.get(tt_jet2_ptRaw)<<" event.jet2_pt = "<<event.get(tt_jet2_pt)<<std::endl;
    if (njets > 2) {
        jet3 = (baconhep::TJet*)js[2];
	event.set(tt_jet3_pt, jet3->pt);
	event.set(tt_jet3_ptRaw, jet3->ptRaw);
	//        event.jet3_pt = jet3->pt;
	//        event.jet3_ptRaw = jet3->ptRaw;
    }
    float pt_ave = (event.get(tt_jet1_pt) + event.get(tt_jet2_pt))/2;
    event.set(tt_pt_ave,pt_ave);
    //    event.pt_ave = pt_ave;
    std::cout<<"pt_ave = "<<pt_ave<<endl;
    event.set(tt_gen_pthat,0);//set default values for DATA
    event.set(tt_gen_weight, 0);//set default values for DATA
    if(is_mc){ /// apply for MC only
        const baconhep::TGenEventInfo & geninfo = event.get(h_genInfo);
        baconhep::TGenEventInfo* genInfo= new baconhep::TGenEventInfo(geninfo);
	event.set(tt_gen_pthat,genInfo->pthat);
	event.set(tt_gen_weight, genInfo->weight);
        // event.gen_pthat    = genInfo->pthat;
        // event.gen_weight   = genInfo->weight;

        //! Reweighting
        //event.weight = event.weight * genInfo->weight * mcweight.getPuReweighting("Asympt", 69) * mcweight.getEvReweighting(0, "Asympt", 69);
	if ((event.get(tt_pt_ave) - event.get(tt_gen_pthat))/event.get(tt_gen_pthat) > 1) return false;
         //0 = central; 1 = scale up; -1 = scale down; 99 = no scale
         //MC option:  Asympt; or Flat
         //minBiasXsec for PU:  69; or 80 mb
        //11 = set 1; 12 = set2 ...
	event.weight = event.weight * event.get(tt_gen_weight) * mcweight.getPuReweighting("Asympt", 69)* mcweight.getEvReweighting(0, "Asympt", 69);
	 std::cout<<"event.weight = "<<event.weight<<endl;
    }


    float probejet_eta  = -99.;
    float probejet_pt   = 0;
    float probejet_phi  = -99.;
    float probejet_ptRaw = -99.;

    float barrel_eta    = -99.;
    float barrel_pt     = 0;
    float barrel_phi    = -99.;
    float barrel_ptRaw  = -99.;

    float rel_r     = -99.;
    float mpf_r     = -99.;
    float asymmetry = -99.;

    TVector2 pt, met;
    TVector2* MET = new TVector2(1,1);
    //    MET->SetMagPhi(eventInfo->pfMET ,eventInfo->pfMETphi);
    MET->SetMagPhi(eventInfo.pfMET ,eventInfo.pfMETphi);

    //    met.Set(eventInfo->pfMET * cos(eventInfo->pfMETphi),eventInfo->pfMET * sin(eventInfo->pfMETphi));
    met.Set(eventInfo.pfMET * cos(eventInfo.pfMETphi),eventInfo.pfMET * sin(eventInfo.pfMETphi));

    int ran = rand();
    int numb = ran % 2 + 1;
    if ((fabs(jet1->eta)<s_eta_barr)&&(fabs(jet2->eta)<s_eta_barr)) {
        if(numb==1){
            probejet_eta = jet2->eta;
            probejet_pt = event.get(tt_jet2_pt);
            probejet_phi = jet2->phi;
            probejet_ptRaw = event.get(tt_jet2_ptRaw);

            barrel_eta = jet1->eta;
            barrel_pt = event.get(tt_jet1_pt);
            barrel_phi = jet1->phi;
            barrel_ptRaw = event.get(tt_jet1_ptRaw);

            asymmetry = (event.get(tt_jet2_pt) - event.get(tt_jet1_pt))/(event.get(tt_jet2_pt) + event.get(tt_jet1_pt));
            rel_r = event.get(tt_jet2_pt) / event.get(tt_jet1_pt);

            pt.Set(event.get(tt_jet1_pt) * cos(jet1->phi),event.get(tt_jet1_pt) * sin(jet1->phi));
            mpf_r = 1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py());
        }
        if(numb==2){
            probejet_eta = jet1->eta;
            probejet_pt = event.get(tt_jet1_pt);
            probejet_phi = jet1->phi;
            probejet_ptRaw = event.get(tt_jet1_ptRaw);

            barrel_eta = jet2->eta;
            barrel_pt = event.get(tt_jet2_pt);
            barrel_phi = jet2->phi;
            barrel_ptRaw = event.get(tt_jet2_ptRaw);

            asymmetry = (event.get(tt_jet1_pt) - event.get(tt_jet2_pt))/(event.get(tt_jet1_pt) + event.get(tt_jet2_pt));
            rel_r = event.get(tt_jet1_pt) / event.get(tt_jet2_pt);

            pt.Set(event.get(tt_jet2_pt) * cos(jet2->phi),event.get(tt_jet2_pt) * sin(jet2->phi));
            mpf_r = 1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py());

        }
    } else if ((fabs(jet1->eta)<s_eta_barr)||(fabs(jet2->eta)<s_eta_barr)){
        if(fabs(jet1->eta)<s_eta_barr){
            probejet_eta = jet2->eta;
            probejet_pt = event.get(tt_jet2_pt);
            probejet_phi = jet2->phi;
            probejet_ptRaw = event.get(tt_jet2_ptRaw);

            barrel_eta = jet1->eta;
            barrel_pt = event.get(tt_jet1_pt);
            barrel_phi = jet1->phi;
            barrel_ptRaw = event.get(tt_jet1_ptRaw);

            asymmetry = (event.get(tt_jet2_pt) - event.get(tt_jet1_pt))/(event.get(tt_jet2_pt) + event.get(tt_jet1_pt));
            rel_r = event.get(tt_jet2_pt) / event.get(tt_jet1_pt);

            pt.Set(event.get(tt_jet1_pt) * cos(jet1->phi),event.get(tt_jet1_pt) * sin(jet1->phi));
            mpf_r = 1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py());

        }else{
            probejet_eta = jet1->eta;
            probejet_pt = event.get(tt_jet1_pt);
            probejet_phi = jet1->phi;
            probejet_ptRaw = event.get(tt_jet1_ptRaw);

            barrel_eta = jet2->eta;
            barrel_pt = event.get(tt_jet2_pt);
            barrel_phi = jet2->phi;
            barrel_ptRaw = event.get(tt_jet2_ptRaw);

            asymmetry = (event.get(tt_jet1_pt) - event.get(tt_jet2_pt))/(event.get(tt_jet1_pt) + event.get(tt_jet2_pt));
            rel_r = event.get(tt_jet1_pt) / event.get(tt_jet2_pt);

            pt.Set(event.get(tt_jet2_pt) * cos(jet2->phi),event.get(tt_jet2_pt) * sin(jet2->phi));
            mpf_r = 1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py());
        }

    }

    event.set(tt_probejet_eta,probejet_eta);
    event.set(tt_probejet_phi,probejet_phi);
    event.set(tt_probejet_pt,probejet_pt);
    event.set(tt_probejet_ptRaw,probejet_ptRaw);
    event.set(tt_barreljet_eta,barrel_eta);
    event.set(tt_barreljet_phi,barrel_phi);
    event.set(tt_barreljet_pt,barrel_pt);
    event.set(tt_barreljet_ptRaw,barrel_ptRaw);
    event.set(tt_asymmetry,asymmetry);
    event.set(tt_rel_r,rel_r);
    event.set(tt_mpf_r,mpf_r);

    event.set(tt_nvertices,pvs.GetEntries());
    //    event.nvertices = pvs.GetEntries();
//     if ((event.nvertices < 14.) || (event.nvertices >= 16.) ) return false;

    float alpha = 0.;
    if (njets > 2) {
      alpha = (2*(event.get(tt_jet3_pt)))/(event.get(tt_jet1_pt) + event.get(tt_jet2_pt));
    }
    // event.alpha = alpha;
    event.set(tt_alpha,alpha);

    if(!sel.DiJet()) return false;

    h_nocuts->fill(event);


    if(!sel.DiJetAdvanced(event)) return false;

    h_dijet->fill(event);


    h_match->fill(event);

    if(is_data){
        if(!sel.Trigger(event)) return false;

    } 
    if (event.get(tt_alpha) < 0.2) {
      h_sel->fill(event);
    }


    //! fill histos after dijet event selection


    //   if(!sel.goodPVertex()) return false;
    //! for alpha < 0.4
    if (event.get(tt_alpha) < 0.4) {
        for( unsigned int i=0; i < eta_range.size()-1; ++i ){
	  if ((fabs(event.get(tt_probejet_eta))>=eta_range[i])&&(fabs(event.get(tt_probejet_eta))<eta_range[i+1])) h_eta_bins_a04[i].fill(event, ran);
        }
    }
    //! for alpha < 0.3
    if (event.get(tt_alpha) < 0.3) {
        for( unsigned int i=0; i < eta_range.size()-1; ++i ){
	  if ((fabs(event.get(tt_probejet_eta))>=eta_range[i])&&(fabs(event.get(tt_probejet_eta))<eta_range[i+1])) h_eta_bins_a03[i].fill(event, ran);
        }
        for( unsigned int j=0; j < eta_range_mikko.size()-1; ++j ){
	  if ((fabs(event.get(tt_probejet_eta))>=eta_range_mikko[j])&&(fabs(event.get(tt_probejet_eta))<eta_range_mikko[j+1])) {
                for( unsigned int i=0; i < pt_range.size()-1; ++i ){
		  if ((event.get(tt_pt_ave)>=pt_range[i])&&(event.get(tt_pt_ave)<pt_range[i+1])) {
                        h_eta_bins_mikko_a30[j*(pt_range.size()-1)+i].fill(event, ran);
                    }
                }
            }
        }
    }
    //! for alpha < 0.2
    if (event.get(tt_alpha) < 0.2) {
        for( unsigned int i=0; i < eta_range.size()-1; ++i ){
	  if ((fabs(event.get(tt_probejet_eta))>=eta_range[i])&&(fabs(event.get(tt_probejet_eta))<eta_range[i+1])) h_eta_bins[i].fill(event, ran);
        }

        for( unsigned int i=0; i < pt_range.size()-1; ++i ){
	  if ((event.get(tt_pt_ave) >= pt_range[i])&&(event.get(tt_pt_ave) < pt_range[i+1])) h_pt_bins[i].fill(event, ran);//j*pt_range.size()+i
        }
        for( unsigned int j=0; j < eta_range_mikko.size()-1; ++j ){
	  if ((fabs(event.get(tt_probejet_eta))>=eta_range_mikko[j])&&(fabs(event.get(tt_probejet_eta))<eta_range_mikko[j+1])) {
                for( unsigned int i=0; i < pt_range.size()-1; ++i ){
		  if ((event.get(tt_pt_ave)>=pt_range[i])&&(event.get(tt_pt_ave)<pt_range[i+1])) {
                        h_eta_bins_mikko_a20[j*(pt_range.size()-1)+i].fill(event, ran);
                    }
                }
            }
        }
        for( unsigned int j=0; j < eta_range.size()-1; ++j ){
	  if ((fabs(event.get(tt_probejet_eta))>=eta_range[j])&&(fabs(event.get(tt_probejet_eta))<eta_range[j+1])) {
                for( unsigned int i=0; i < pt_range.size()-1; ++i ){
		  if ((event.get(tt_pt_ave)>=pt_range[i])&&(event.get(tt_pt_ave)<pt_range[i+1])) {
                        h_noalpha_bins[j*(pt_range.size()-1)+i].fill(event, ran);//j*pt_range.size()+i
                    }
                }
            }
        }
    }
    //! for alpha < 0.15
    if (event.get(tt_alpha) < 0.15) {
        for( unsigned int j=0; j < eta_range_mikko.size()-1; ++j ){
	  if ((fabs(event.get(tt_probejet_eta))>=eta_range_mikko[j])&&(fabs(event.get(tt_probejet_eta))<eta_range_mikko[j+1])) {
                for( unsigned int i=0; i < pt_range.size()-1; ++i ){
		  if ((event.get(tt_pt_ave)>=pt_range[i])&&(event.get(tt_pt_ave)<pt_range[i+1])) {
                        h_eta_bins_mikko_a15[j*(pt_range.size()-1)+i].fill(event, ran);
                    }
                }
            }
        }
    }

    //! for alpha < 0.1

    if (event.get(tt_alpha) < 0.1) {
        for( unsigned int i=0; i < pt_range.size()-1; ++i ){
	  if ((event.get(tt_pt_ave) >= pt_range[i])&&(event.get(tt_pt_ave) < pt_range[i+1])) h_pt_bins_a01[i].fill(event, ran);//j*pt_range.size()+i
        }
        for( unsigned int i=0; i < eta_range.size()-1; ++i ){
	  if ((fabs(event.get(tt_probejet_eta))>=eta_range[i])&&(fabs(event.get(tt_probejet_eta))<eta_range[i+1])) h_eta_bins_a01[i].fill(event, ran);
        }
        for( unsigned int j=0; j < eta_range_mikko.size()-1; ++j ){
	  if ((fabs(event.get(tt_probejet_eta))>=eta_range_mikko[j])&&(fabs(event.get(tt_probejet_eta))<eta_range_mikko[j+1])) {
                for( unsigned int i=0; i < pt_range.size()-1; ++i ){
		  if ((event.get(tt_pt_ave)>=pt_range[i])&&(event.get(tt_pt_ave)<pt_range[i+1])) {
                        h_eta_bins_mikko_a10[j*(pt_range.size()-1)+i].fill(event, ran);
                    }
                }
            }
        }
    }


    //cSetTree.fillTree(eventInfo->evtNum, eventInfo->runNum, pt_ave, probejet_eta, barrel_eta, probejet_phi, barrel_phi, probejet_pt, barrel_pt, event.jet1_pt, event.jet2_pt, jet1->phi, jet2->phi, jet1->eta, jet2->eta, event.jet1_ptRaw, event.jet2_ptRaw, probejet_ptRaw, barrel_ptRaw);

    //     }

    return true;
  }

  // as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
  // make sure the ExampleModule is found by class name. This is ensured by this macro:
  UHH2_REGISTER_ANALYSIS_MODULE(TestModule)

}
