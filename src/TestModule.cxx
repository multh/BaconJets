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

#include "UHH2/bacondataformats/interface/TGenEventInfo.hh"
#include "UHH2/bacondataformats/interface/TJet.hh"
#include "UHH2/bacondataformats/interface/TEventInfo.hh"
#include "UHH2/bacondataformats/interface/BaconAnaDefs.hh"

#include "TClonesArray.h"
#include "TString.h"




using namespace std;
using namespace uhh2;

namespace uhh2bacon {

class TestModule: public AnalysisModule {
public:

    explicit TestModule(Context & ctx);
    virtual bool process(Event & event) override;

private:

  Event::Handle<TClonesArray> h_jets  ;
  Event::Handle<baconhep::TEventInfo> h_eventInfo;
  Event::Handle<baconhep::TGenEventInfo> h_genInfo;

  std::unique_ptr<Hists> h_nocuts, h_sel, h_dijet, h_match;
//   std::vector<double> eta_range, pt_range, alpha_range;
  std::vector<JECAnalysisHists> h_pt_bins, h_noalpha_bins;

  Selection sel;
  JetCorrections jetcorr;
  McWeight mcweight;
  bool is_mc;

};


TestModule::TestModule(Context & ctx) :
  sel(ctx),
  jetcorr(ctx),
  mcweight(ctx)
{
  auto dataset_type = ctx.get("dataset_type");
  is_mc = dataset_type  == "MC";
  h_jets = ctx.declare_event_input<TClonesArray>("Jet05");
  h_eventInfo = ctx.declare_event_input<baconhep::TEventInfo>("Info");
  if(is_mc){ /// apply for MC only
    h_genInfo = ctx.declare_event_input<baconhep::TGenEventInfo>("GenEvtInfo");
  }

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
  std::vector<std::string> alpha_range_name;
  for( unsigned int k=0; k < alpha_range.size(); ++k ){
    char alpha_buffer [50];
    sprintf (alpha_buffer, "%5.3f", alpha_range[k]);
    alpha_range_name.push_back(alpha_buffer);
  }
  //cout << "alpha range "<<alpha_range_name[3]<<"eta range "<<eta_range_name[3]<< "pt range "<<pt_range_name[3]<< endl;

  // histos for the kFSR extrapolations
  for( unsigned int k=0; k < alpha_range.size()-1; ++k ){
    for( unsigned int j=0; j < eta_range.size()-1; ++j ){
        for( unsigned int i=0; i < pt_range.size()-1; ++i ){
            h_pt_bins.push_back(JECAnalysisHists(ctx,(std::string)("alpha_"+alpha_range_name[k]+"_"+alpha_range_name[k+1]+"/eta_"+eta_range_name[j]+"_"+eta_range_name[j+1]+"/pt_"+pt_range_name[i]+"_"+pt_range_name[i+1])));
        }
    }
  }

  // histos for the pT extrapolations
  for( unsigned int j=0; j < eta_range.size()-1; ++j ){
    for( unsigned int i=0; i < pt_range.size()-1; ++i ){
      h_noalpha_bins.push_back(JECAnalysisHists(ctx,(std::string)("/eta_"+eta_range_name[j]+"_"+eta_range_name[j+1]+"/pt_"+pt_range_name[i]+"_"+pt_range_name[i+1])));
    }
  }


}


bool TestModule::process(Event & event) {

  sel.SetEvent(event);
  jetcorr.SetEvent(event);
  mcweight.SetEvent(event);

  const TClonesArray & js = event.get(h_jets);
  baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
  baconhep::TJet* jet2 = (baconhep::TJet*)js[1];
  Int_t njets = js.GetEntries();

  const baconhep::TEventInfo & info = event.get(h_eventInfo);
  baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);

  if(is_mc){ /// apply for MC only
    // to reweight MC
//    cout << "w/o mc weight = " << event.weight << endl;
    const baconhep::TGenEventInfo & geninfo = event.get(h_genInfo);
    baconhep::TGenEventInfo* genInfo= new baconhep::TGenEventInfo(geninfo);
    event.weight = event.weight * genInfo->weight * mcweight.getPuReweighting() * mcweight.getEvReweighting();
  //  cout << "with mc weight = " << event.weight << endl;

  // doing the matching from GEN to RECO
    if(!jetcorr.JetMatching()) return false;

    // JER smearing
    if(!jetcorr.JetResolutionSmearer()) return false;
  }

  h_nocuts->fill(event);

  if(!sel.DiJet()) return false;
  h_dijet->fill(event);


  if(!sel.DiJetAdvanced()) return false;
  h_match->fill(event);


  if(!sel.Trigger()) return false;


  if(!sel.goodPVertex()) return false;


  if(!sel.jetIds(s_working_point_csv_threshold)) return false;

  double probejet_eta = -99.;
  double probejet_pt = -99.;

  int ran = rand();
  int numb = ran % 2 + 1;
  if ((fabs(jet1->eta)<s_eta_barr)&&(fabs(jet2->eta)<s_eta_barr)) {
    if(numb==1){
        probejet_eta = jet2->eta;
        probejet_pt = jet2->pt;
    }
    if(numb==2){
        probejet_eta = jet1->eta;
        probejet_pt = jet1->pt;
    }
  } else if ((fabs(jet1->eta)<s_eta_barr)||(fabs(jet2->eta)<s_eta_barr)){
    if(fabs(jet1->eta)<s_eta_barr){
        probejet_eta = jet2->eta;
        probejet_pt = jet2->pt;
    }
    else{
        probejet_eta = jet1->eta;
        probejet_pt = jet1->pt;
    }
  }
  double alpha = 0.;
  if (njets > 2) {
    baconhep::TJet* jet3 = (baconhep::TJet*)js[2];
        alpha = (2*(jet3->pt))/(jet1->pt + jet2->pt);
        //cout << "alpha = "<< alpha << endl;
  }

  // fill histos for the kFSR extrapolations
  // no cut on alpha is required since we want to extrapolate as a function of alpha
  for( unsigned int k=0; k < alpha_range.size()-1; ++k ){
    if ((alpha>=alpha_range[k])&&(alpha<alpha_range[k+1])) {
        for( unsigned int j=0; j < eta_range.size()-1; ++j ){
            if ((fabs(probejet_eta)>=eta_range[j])&&(fabs(probejet_eta)<eta_range[j+1])) {
                for( unsigned int i=0; i < pt_range.size()-1; ++i ){
                    if ((probejet_pt>=pt_range[i])&&(probejet_pt<pt_range[i+1])) {
                        h_pt_bins[k*(eta_range.size()-1)*(pt_range.size()-1)+j*(pt_range.size()-1)+i].fill(event, ran);//j*pt_range.size()+i
                        //cout <<"eta range = "<< eta_range[j]<<" - "<< eta_range[j+1]<< "pt range = "<< pt_range[i]<<" - "<< pt_range[i+1]<<endl;
                        //cout <<"eta value = "<< fabs(probejet_eta) << " pt value = "<< probejet_pt <<endl;
                    }
                }
            }
        }
    }
  }

  // alpha<0.2
  if(!sel.AlphaCut()) return false;
  h_sel->fill(event);

  // fill histos for the pT extrapolations
  // alpha<0.2 needs to be required!
  for( unsigned int j=0; j < eta_range.size()-1; ++j ){
    if ((fabs(probejet_eta)>=eta_range[j])&&(fabs(probejet_eta)<eta_range[j+1])) {
      for( unsigned int i=0; i < pt_range.size()-1; ++i ){
	if ((probejet_pt>=pt_range[i])&&(probejet_pt<pt_range[i+1])) {
	  h_noalpha_bins[j*(pt_range.size()-1)+i].fill(event, ran);//j*pt_range.size()+i
	}
      }
    }
  }



  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the ExampleModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TestModule)

}
