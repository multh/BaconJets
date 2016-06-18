#include <iostream>
#include <memory>
#include <stdlib.h>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
//#include "UHH2/core/include/EventHelper.h"
#include "../include/JECAnalysisHists.h"

#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/LumiSelection.h>
#include <UHH2/common/include/TriggerSelection.h>

#include "UHH2/BaconJets/include/selection.h"
#include "UHH2/BaconJets/include/jet_corrections.h"
#include "UHH2/BaconJets/include/mc_weight.h"
#include "UHH2/BaconJets/include/constants.h"
#include "UHH2/BaconJets/include/TSetTree.h"
//#include "UHH2/BaconJets/include/dijet_event.h"
#include "UHH2/core/include/Jet.h"

// #include "UHH2/BaconTrans/baconheaders/TGenEventInfo.hh"
// #include "UHH2/BaconTrans/baconheaders/TJet.hh"
// #include "UHH2/BaconTrans/baconheaders/TEventInfo.hh"
// #include "UHH2/BaconTrans/baconheaders/BaconAnaDefs.hh"

//#include "UHH2/BaconJets/include/pileup_data.h"
//#include "UHH2/BaconJets/include/data_corrections.h"
#include "TClonesArray.h"
#include "TString.h"

//#include "UHH2/common/include/MCWeight.h"

//TTree   *fCurrentTree;
Int_t   Runnr;
Int_t   Eventnr;
//TFile   *fCurrentTreeFile;
using namespace std;
using namespace uhh2;
//using namespace baconhep;
//using uhh2::detail::EventHelper;


  class TestModule: public uhh2::AnalysisModule {

  public:
    explicit TestModule(uhh2::Context&);
    virtual bool process(uhh2::Event&) override;
    ~TestModule();

  protected:
    // correctors
    std::unique_ptr<JetCorrector> jet_corrector;
    // selections
    std::unique_ptr<uhh2::Selection> lumi_sel;
    //    std::unique_ptr<uhh2::AndSelection> trigger_sel;
    std::unique_ptr<uhh2::Selection> trigger40_sel;
    std::unique_ptr<uhh2::Selection> trigger60_sel;
    std::unique_ptr<uhh2::Selection> trigger80_sel;
    std::unique_ptr<uhh2::Selection> trigger140_sel;
    std::unique_ptr<uhh2::Selection> trigger200_sel;
    std::unique_ptr<uhh2::Selection> trigger260_sel;
    std::unique_ptr<uhh2::Selection> trigger320_sel;
    std::unique_ptr<uhh2::Selection> trigger400_sel;
    std::unique_ptr<uhh2::Selection> trigger500_sel;
    std::unique_ptr<uhh2::Selection> trigger60_HFJEC_sel;
    std::unique_ptr<uhh2::Selection> trigger80_HFJEC_sel;
    std::unique_ptr<uhh2::Selection> trigger100_HFJEC_sel;
    std::unique_ptr<uhh2::Selection> trigger160_HFJEC_sel;
    std::unique_ptr<uhh2::Selection> trigger220_HFJEC_sel;
    std::unique_ptr<uhh2::Selection> trigger300_HFJEC_sel;
    //// Data/MC scale factors
    std::unique_ptr<uhh2::AnalysisModule> pileupSF;

   


    // float gen_pthat; //pt hat (from QCD simulation)
    //  float gen_weight;// weight from MC
    //  float jet1_pt, jet2_pt, jet3_pt; //leading, subleading and the 3rd jet pt (corrected)
    //  float jet1_ptRaw, jet2_ptRaw, jet3_ptRaw;//leading, subleading and the 3rd jet pt (not corrected)
    //  float nvertices;//number of vertices
    //  float probejet_eta, probejet_phi, probejet_pt, probejet_ptRaw;// probe jet parameters
    //  float barreljet_eta, barreljet_phi, barreljet_pt, barreljet_ptRaw;//reference jet parameters 
    //  float pt_ave;//pt average of leading and subleading jets
    //  float alpha;// pt of the 3rd jet divided by  pt_ave
    //  float rel_r, mpf_r; //responces from pT-balance and MPF method
    //  float asymmetry;//asymmetry=(p_{T}^{probe}-p_{T}^{barrel})/(p_{T}^{probe}+p_{T}^{barrel})
    //  float nPU;//number of pile-up vertices; a-ka mu in MC is the poisson mean of pileup that we use for pileup reweighing
    //  float ev_weight;//weight of the event
    //  float jets_pt;//sum of additional jets pT (does _not_ include leading and subleading jets)
    //  int nJets;//number of jets

    Event::Handle<float> tt_gen_pthat; Event::Handle<float> tt_gen_weight;
    Event::Handle<float> tt_jet1_pt;     Event::Handle<float> tt_jet2_pt;     Event::Handle<float> tt_jet3_pt;
    Event::Handle<float> tt_jet1_ptRaw;  Event::Handle<float> tt_jet2_ptRaw;  Event::Handle<float> tt_jet3_ptRaw;
    Event::Handle<int> tt_nvertices;
    Event::Handle<float> tt_probejet_eta;  Event::Handle<float> tt_probejet_phi; Event::Handle<float> tt_probejet_pt; Event::Handle<float> tt_probejet_ptRaw;
    Event::Handle<float> tt_barreljet_eta;  Event::Handle<float> tt_barreljet_phi; Event::Handle<float> tt_barreljet_pt; Event::Handle<float> tt_barreljet_ptRaw;
    Event::Handle<float> tt_pt_ave;
    Event::Handle<float> tt_alpha;
    Event::Handle<float> tt_rel_r; Event::Handle<float> tt_mpf_r; Event::Handle<float> tt_asymmetry; Event::Handle<int> tt_nPU;
    Event::Handle<float> tt_ev_weight;
    Event::Handle<float> tt_jets_pt;//sum of jets pT
    Event::Handle<int> tt_jet_n;//number of jets
    Event::Handle<float> tt_rho;//event energy density
    Event::Handle<int> tt_nGoodvertices;

    //    Event::Handle<dijet_event> tt_dijet_event;
    //  private:
    // //    std::unique_ptr<EventHelper> eh;
    // // Event::Handle<TClonesArray> h_jets  ;
    // //    Event::Handle<std::vector<uhh2bacon::TJet>> h_jets;
    // // Event::Handle<baconhep::TEventInfo> h_eventInfo;
    // // Event::Handle<baconhep::TGenEventInfo> h_genInfo;

    // Event::Handle<TClonesArray> h_pv;
    // //std::unique_ptr<Hists> h_nocuts, h_sel, h_dijet, h_match;
    std::unique_ptr<JECAnalysisHists> h_nocuts, h_sel, h_dijet, h_match;
    std::unique_ptr<JECAnalysisHists> h_trg40, h_trg60, h_trg80, h_trg140, h_trg200,h_trg260,h_trg320,h_trg400,h_trg500;
    std::unique_ptr<JECAnalysisHists> h_trgHF60, h_trgHF80,h_trgHF100, h_trgHF160,h_trgHF220, h_trgHF300;     
    uhh2bacon::Selection sel;
    // //    JetCorrections jetcorr;
    // std::unique_ptr<McWeight> mcweight; //todo: do we need it?
    // bool is_mc;
    // bool is_data;
    // bool is_mc_reweight;
    // //  TSetTree cSetTree;
    // uhh2bacon::PileupData  pileupData;
    // double jets_pt;//sum of jets pT



    // // //Additional vars in Event, specific for dijet
    
    // std::unique_ptr<MCLumiWeight> fMCLumiWeight;
  };

  TestModule::TestModule(uhh2::Context & ctx) :
    sel(ctx)
  {
    const bool isMC = (ctx.get("dataset_type") == "MC");
    //// COMMON MODULES
    if(!isMC) lumi_sel.reset(new LumiSelection(ctx));

    if(!isMC){
    const std::string& trigger40 = ctx.get("trigger40", "NULL");
    const std::string& trigger60 = ctx.get("trigger60", "NULL");
    const std::string& trigger80 = ctx.get("trigger80", "NULL");
    const std::string& trigger140 = ctx.get("trigger140", "NULL");
    const std::string& trigger200 = ctx.get("trigger200", "NULL");
    const std::string& trigger260 = ctx.get("trigger260", "NULL");
    const std::string& trigger320 = ctx.get("trigger320", "NULL");
    const std::string& trigger400 = ctx.get("trigger400", "NULL");
    const std::string& trigger500 = ctx.get("trigger500", "NULL");
    const std::string& trigger60_HFJEC = ctx.get("trigger60_HFJEC", "NULL");
    const std::string& trigger80_HFJEC = ctx.get("trigger80_HFJEC", "NULL");
    const std::string& trigger100_HFJEC = ctx.get("trigger100_HFJEC", "NULL");
    const std::string& trigger160_HFJEC = ctx.get("trigger160_HFJEC", "NULL");
    const std::string& trigger220_HFJEC = ctx.get("trigger220_HFJEC", "NULL");
    const std::string& trigger300_HFJEC = ctx.get("trigger300_HFJEC", "NULL");

      // const std::string& trigger = ctx.get("trigger", "NULL");
      if(trigger40 != "NULL") trigger40_sel.reset(new TriggerSelection(trigger40));
      else trigger40_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger60 != "NULL") trigger60_sel.reset(new TriggerSelection(trigger60));
      else trigger60_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger80 != "NULL") trigger80_sel.reset(new TriggerSelection(trigger80));
      else trigger80_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger140 != "NULL") trigger140_sel.reset(new TriggerSelection(trigger140));
      else trigger140_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger200 != "NULL") trigger200_sel.reset(new TriggerSelection(trigger200));
      else trigger200_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger260 != "NULL") trigger260_sel.reset(new TriggerSelection(trigger260));
      else trigger260_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger320 != "NULL") trigger320_sel.reset(new TriggerSelection(trigger320));
      else trigger320_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger400 != "NULL") trigger400_sel.reset(new TriggerSelection(trigger400));
      else trigger400_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger500 != "NULL") trigger500_sel.reset(new TriggerSelection(trigger500));
      else trigger500_sel.reset(new uhh2::AndSelection(ctx));

      if(trigger60_HFJEC != "NULL") trigger60_HFJEC_sel.reset(new TriggerSelection(trigger60_HFJEC));
      else trigger60_HFJEC_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger80_HFJEC != "NULL") trigger80_HFJEC_sel.reset(new TriggerSelection(trigger80_HFJEC));
      else trigger80_HFJEC_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger100_HFJEC != "NULL") trigger100_HFJEC_sel.reset(new TriggerSelection(trigger100_HFJEC));
      else trigger100_HFJEC_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger160_HFJEC != "NULL") trigger160_HFJEC_sel.reset(new TriggerSelection(trigger160_HFJEC));
      else trigger160_HFJEC_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger220_HFJEC != "NULL") trigger220_HFJEC_sel.reset(new TriggerSelection(trigger220_HFJEC));
      else trigger220_HFJEC_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger300_HFJEC != "NULL") trigger300_HFJEC_sel.reset(new TriggerSelection(trigger300_HFJEC));
      else trigger300_HFJEC_sel.reset(new uhh2::AndSelection(ctx));
    }

    //Jet collection used in the analysis is defined in xml config, parameters: "JetCollection"(input),"JetLabel"(output,label)
    const std::string& jetLabel = ctx.get("JetLabel");
    std::vector<std::string> JEC_corr;
    if(isMC){
      if(jetLabel == "AK4CHS") JEC_corr = JERFiles::Spring16_25ns_L123_AK4PFchs_MC;
      if(jetLabel == "AK8CHS") JEC_corr = JERFiles::Spring16_25ns_L123_AK8PFchs_MC;
      if(jetLabel == "AK4PUPPI") JEC_corr = JERFiles::Spring16_25ns_L123_AK4PFPuppi_MC;
      if(jetLabel == "AK8PUPPI") JEC_corr = JERFiles::Spring16_25ns_L123_AK8PFPuppi_MC;
    }
    else {
      if(jetLabel == "AK4CHS") JEC_corr = JERFiles::Spring16_25ns_L123_noRes_AK4PFchs_DATA;
      if(jetLabel == "AK8CHS") JEC_corr = JERFiles::Spring16_25ns_L123_noRes_AK4PFchs_DATA;//ToDo: change for different jet collections!!!
      if(jetLabel == "AK4PUPPI") JEC_corr = JERFiles::Spring16_25ns_L123_noRes_AK4PFchs_DATA;
      if(jetLabel == "AK8PUPPI") JEC_corr = JERFiles::Spring16_25ns_L123_noRes_AK4PFchs_DATA;
    }
    jet_corrector.reset(new JetCorrector(ctx, JEC_corr));

    //output
    // ctx.undeclare_all_event_output();   
    // //pileup (define it after undeclaring all other variables to keep the weights in the output)
    // pileupSF.reset(new MCPileupReweight(ctx));

    //    tt_dijet_event = ctx.declare_event_output<dijet_event>("dijet");
    //Store only vars needed for the dijet analysis
    tt_gen_pthat = ctx.declare_event_output<float>("gen_pthat");
    tt_gen_weight = ctx.declare_event_output<float>("gen_weight");
    tt_jet1_pt = ctx.declare_event_output<float>("jet1_pt");
    tt_jet2_pt = ctx.declare_event_output<float>("jet2_pt");
    tt_jet3_pt = ctx.declare_event_output<float>("jet3_pt");
    tt_jet1_ptRaw = ctx.declare_event_output<float>("jet1_ptRaw");
    tt_jet2_ptRaw = ctx.declare_event_output<float>("jet2_ptRaw");
    tt_jet3_ptRaw = ctx.declare_event_output<float>("jet3_ptRaw");
    tt_nvertices = ctx.declare_event_output<int>("nvertices");
    tt_nGoodvertices = ctx.declare_event_output<int>("nGoodvertices");
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
    tt_ev_weight = ctx.declare_event_output<float>("weight");
    tt_jets_pt= ctx.declare_event_output<float>("sum_jets_pt");
    tt_jet_n= ctx.declare_event_output<int>("Njet");
    tt_rho = ctx.declare_event_output<float>("rho");

    h_nocuts.reset(new JECAnalysisHists(ctx,"noCuts"));
    h_dijet.reset(new JECAnalysisHists(ctx,"diJet"));
    h_match.reset(new JECAnalysisHists(ctx,"JetMatching"));
    h_sel.reset(new JECAnalysisHists(ctx,"Selection"));

    h_trg40.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve40"));
    h_trg60.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve60"));
    h_trg80.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve80"));
    h_trg140.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve140"));
    h_trg200.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve200"));
    h_trg260.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve260"));
    h_trg320.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve320"));
    h_trg400.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve400"));
    h_trg500.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve500"));

    h_trgHF60.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve60_HFJEC"));
    h_trgHF80.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve80_HFJEC"));
    h_trgHF100.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve100_HFJEC"));
    h_trgHF160.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve160_HFJEC"));
    h_trgHF220.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve220_HFJEC"));
    h_trgHF300.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve300_HFJEC"));

  };
 //  TestModule::TestModule(Context & ctx) :
//     sel(ctx),
//     mcweight(ctx),
//     pileupData(ctx)
//   {



//     h_eventInfo = ctx.declare_event_input<baconhep::TEventInfo>("Info");
//     if(is_mc){ /// apply for MC only
//       h_genInfo = ctx.declare_event_input<baconhep::TGenEventInfo>("GenEvtInfo");
//     }
// //     h_pv = ctx.get_handle<TClonesArray>("PV");
//     h_pv = ctx.declare_event_input<TClonesArray>("PV");

//     //    eh.reset(new EventHelper(ctx));
//     h_nocuts.reset(new JECAnalysisHists(ctx,"noCuts"));
//     h_dijet.reset(new JECAnalysisHists(ctx,"diJet"));
//     h_match.reset(new JECAnalysisHists(ctx,"JetMatching"));
//     h_sel.reset(new JECAnalysisHists(ctx,"Selection"));
//     fMCLumiWeight.reset(new MCLumiWeight(ctx));






//   }

  TestModule::~TestModule() {
    //  cSetTree.general();
    //      cPuData.general();

  }

  bool TestModule::process(Event & event) {
    //    cout<<"NEW EVENT"<<endl;

    /* CMS-certified luminosity sections */
    if(event.isRealData){
      //      if(!lumi_sel->passes(event)) return false;
      if(!lumi_sel->passes(event)){
	//	std::cout<<"Sorry, bad lumi sec"<<std::endl;
	return false;
      }
    }
    //    std::cout<<"JSON OK"<<std::endl;
    //// Data/MC scale factors
    // //pileup
    // pileupSF->process(event);

    sort_by_pt<Jet>(*event.jets);
    const int jet_n = event.jets->size();
    //    std::cout<<"jet_n = "<<jet_n<<std::endl;
    if(jet_n<2) return false;
    jet_corrector->process(event);
    //todo: add GenericJetResolutionSmearer?
    //    std::cout<<"Jet_Corrector OK"<<std::endl;
    Jet* jet1 = &event.jets->at(0);// leading jet
    Jet* jet2 = &event.jets->at(1);// sub-leading jet
    float jet1_pt = jet1->pt(); float jet2_pt = jet2->pt();
    float pt_ave = (jet1_pt + jet2_pt)/2.;
    //// HLT selection
    bool pass_trigger40=false; bool pass_trigger60=false; bool pass_trigger80=false;
    bool pass_trigger140=false; bool pass_trigger200=false; bool pass_trigger260=false;
    bool pass_trigger320=false; bool pass_trigger400=false; bool pass_trigger500=false;
    bool pass_trigger60_HFJEC=false; bool pass_trigger80_HFJEC=false;
    bool pass_trigger100_HFJEC=false; bool pass_trigger160_HFJEC=false;
    bool pass_trigger220_HFJEC=false; bool pass_trigger300_HFJEC=false;
    double trg_thresh[9] = {56,78,100,168,232,300,366,453,562};
    double trgHF_thresh[6] = {77,131,154,244,321,426};
    if(event.isRealData){
      // cout << " =================== " << endl;
      // cout << "Available triggers: " << endl;
      // for(const auto & tname : event.get_current_triggernames()){
      // 	cout << tname  << endl;
      // }
      //      cout << " =================== " << endl;
      //       pass_trigger = trigger_sel->passes(event);
      pass_trigger40 = (trigger40_sel->passes(event) && pt_ave>trg_thresh[0]);
      pass_trigger60 = (trigger60_sel->passes(event) && pt_ave>trg_thresh[1]);
      pass_trigger80 = (trigger80_sel->passes(event) && pt_ave>trg_thresh[2]);
      pass_trigger140 = (trigger140_sel->passes(event) && pt_ave>trg_thresh[3]);
      pass_trigger200 = (trigger200_sel->passes(event) && pt_ave>trg_thresh[4]);
      pass_trigger260 = (trigger260_sel->passes(event) && pt_ave>trg_thresh[5]);
      pass_trigger320 = (trigger320_sel->passes(event) && pt_ave>trg_thresh[6]);
      pass_trigger400 = (trigger400_sel->passes(event) && pt_ave>trg_thresh[7]);
      pass_trigger500 = (trigger500_sel->passes(event) && pt_ave>trg_thresh[8]);
      pass_trigger60_HFJEC = (trigger60_HFJEC_sel->passes(event) && pt_ave>trgHF_thresh[0]);
      pass_trigger80_HFJEC = (trigger80_HFJEC_sel->passes(event) && pt_ave>trgHF_thresh[1]);
      pass_trigger100_HFJEC = (trigger100_HFJEC_sel->passes(event) && pt_ave>trgHF_thresh[2]);
      pass_trigger160_HFJEC = (trigger160_HFJEC_sel->passes(event) && pt_ave>trgHF_thresh[3]);
      pass_trigger220_HFJEC = (trigger220_HFJEC_sel->passes(event) && pt_ave>trgHF_thresh[4]);
      pass_trigger300_HFJEC = (trigger300_HFJEC_sel->passes(event) && pt_ave>trgHF_thresh[5]);
      const bool pass_trigger = (pass_trigger40 || pass_trigger60 || pass_trigger140 || pass_trigger200 
				 || pass_trigger260 || pass_trigger320 || pass_trigger400 || pass_trigger500
				 || pass_trigger60_HFJEC || pass_trigger80_HFJEC || pass_trigger100_HFJEC
				 || pass_trigger160_HFJEC || pass_trigger220_HFJEC || pass_trigger300_HFJEC);
	//const bool pass_trigger = (pass_trigger40 || pass_trigger60 || pass_trigger80 || pass_trigger140);//TEST 
      if(!pass_trigger)
	return false;
    }
    //    std::cout<<"HLT OK"<<std::endl;



    //    std::cout<<"eta1 = "<<fabs(jet1->eta())<<" eta2 = "<<fabs(jet2->eta())<<std::endl;
    //  sel.SetEvent(event);
    Jet* jet_probe = jet1; Jet* jet_barrel = jet2;
    if ((fabs(jet1->eta())<s_eta_barr)&&(fabs(jet2->eta())<s_eta_barr)) {
      int ran = rand();
      int numb = ran % 2 + 1;
        if(numb==1){
	  jet_probe = jet2;
	  jet_barrel = jet1;
	}
        if(numb==2){
	  jet_probe = jet1;
	  jet_barrel = jet2;
        }
    } else if ((fabs(jet1->eta())<s_eta_barr)||(fabs(jet2->eta())<s_eta_barr)){
        if(fabs(jet1->eta())<s_eta_barr){
	  jet_probe = jet2;
	  jet_barrel = jet1;
        }else{
	  jet_probe = jet1;
	  jet_barrel = jet2;
        }
    }

    //read or calculated values for dijet events
    float gen_pthat = 0; //pt hat (from QCD simulation) //todo!
    float gen_weight = 0;
    if(!event.isRealData)
      gen_weight = event.weight;
    float nvertices = event.pvs->size(); 
    float nPU = 0 ;//todo for data?
    if(!event.isRealData)
      nPU = event.genInfo->pileup_TrueNumInteractions();



    float ev_weight = event.weight;

    auto factor_raw1 = jet1->JEC_factor_raw();     auto factor_raw2 = jet2->JEC_factor_raw();
    float jet1_ptRaw = jet1_pt*factor_raw1;  float jet2_ptRaw = jet2_pt*factor_raw2;
    float probejet_eta = jet_probe->eta(); 
    float  probejet_phi = jet_probe->phi(); 
    float  probejet_pt = jet_probe->pt(); 
    auto factor_raw_probe = jet_probe->JEC_factor_raw();
    float probejet_ptRaw = probejet_pt*factor_raw_probe;
    float barreljet_eta = jet_barrel->eta(); 
    float  barreljet_phi = jet_barrel->phi(); 
    float  barreljet_pt = jet_barrel->pt(); 
    auto factor_raw_barrel = jet_barrel->JEC_factor_raw();
    float barreljet_ptRaw = barreljet_pt*factor_raw_barrel;

    float jet3_pt = 0; float jet3_ptRaw = 0;
    if(jet_n>2){
      Jet* jet3 = &event.jets->at(2);
      jet3_pt = jet3->pt();
      auto factor_raw3 = jet3->JEC_factor_raw();
      jet3_ptRaw = jet3_pt*factor_raw3;
    }
    float alpha = jet3_pt/pt_ave;
    float asymmetry = (probejet_pt - barreljet_pt)/(probejet_pt + barreljet_pt);
    float rel_r = (1+asymmetry)/(1-asymmetry);
    TVector2 pt, met;
    met.Set(event.met->pt() * cos(event.met->phi()),event.met->pt() * sin(event.met->phi()));
    pt.Set(barreljet_pt * cos(barreljet_phi),barreljet_pt* sin(barreljet_phi));
    float mpf_r = 1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py());
    float jets_pt = 0;
    for(int i=2;i<jet_n;i++){
      jets_pt += ((Jet*)&event.jets->at(i))->pt();
    }

    //fill the containers
    //    dijet_event t_dijet_event;
    event.set(tt_gen_pthat,gen_pthat);
    event.set(tt_gen_weight,gen_weight);
    event.set(tt_jet1_pt,jet1_pt);
    event.set(tt_jet2_pt,jet2_pt);
    event.set(tt_jet3_pt,jet3_pt);
    event.set(tt_jet1_ptRaw,jet1_ptRaw);
    event.set(tt_jet2_ptRaw,jet2_ptRaw);
    event.set(tt_jet3_ptRaw,jet3_ptRaw);
    event.set(tt_nvertices,nvertices);
    event.set(tt_probejet_eta,probejet_eta);
    event.set(tt_probejet_phi,probejet_phi);
    event.set(tt_probejet_pt,probejet_pt);
    event.set(tt_probejet_ptRaw,probejet_ptRaw);
    event.set(tt_barreljet_eta,barreljet_eta);
    event.set(tt_barreljet_phi,barreljet_phi);
    event.set(tt_barreljet_pt,barreljet_pt);
    event.set(tt_barreljet_ptRaw,barreljet_ptRaw);
    event.set(tt_pt_ave,pt_ave);
    event.set(tt_alpha,alpha);
    event.set(tt_asymmetry,asymmetry);
    event.set(tt_rel_r,rel_r);
    event.set(tt_mpf_r,mpf_r);
    event.set(tt_nPU,nPU);
    event.set(tt_ev_weight,ev_weight);
    event.set(tt_jets_pt,jets_pt);
    event.set(tt_jet_n,jet_n);
    event.set(tt_rho,event.rho);    

    sel.SetEvent(event);
    //good primary vertex
    int nGoodVts = sel.goodPVertex();
    //    std::cout<<"nGoodVts = "<<nGoodVts<<std::endl;
    if(nGoodVts<=0) return false;
    event.set(tt_nGoodvertices, nGoodVts);
    if(!sel.DiJet()) return false;
    h_nocuts->fill(event);
    if(!sel.DiJetAdvanced(event)) return false;
    h_dijet->fill(event);
    h_match->fill(event);
    if(event.isRealData){
     if(pass_trigger40) h_trg40->fill(event); 
     if(pass_trigger60) h_trg60->fill(event); 
     if(pass_trigger80) h_trg80->fill(event); 
     if(pass_trigger140) h_trg140->fill(event); 
     if(pass_trigger200) h_trg200->fill(event); 
     if(pass_trigger260) h_trg260->fill(event);
     if(pass_trigger320) h_trg320->fill(event);  
     if(pass_trigger400) h_trg400->fill(event);  
     if(pass_trigger500) h_trg500->fill(event);  
     if(pass_trigger60_HFJEC) h_trgHF60->fill(event);  
     if(pass_trigger80_HFJEC) h_trgHF80->fill(event);  
     if(pass_trigger100_HFJEC) h_trgHF100->fill(event);  
     if(pass_trigger160_HFJEC) h_trgHF160->fill(event);  
     if(pass_trigger220_HFJEC) h_trgHF220->fill(event);  
     if(pass_trigger300_HFJEC) h_trgHF300->fill(event);  
    }
    else
      if(!sel.PtMC(event)) return false; // For MC only one Pt threshold
    if (event.get(tt_alpha) < 0.3) {
      h_sel->fill(event);
    }

//     // cout<<"Yff! "<<event.weight<<endl;
//     event.set(tt_ev_weight,event.weight);

//   std::cout<<"jet1_pt = "<<jet1_pt<<" jet2_pt = "<<jet2_pt<<" jet3_pt = "<<jet3_pt<<std::endl;
//    std::cout<<"jet1_ptRaw = "<<jet1_ptRaw<<" jet2_ptRaw = "<<jet2_ptRaw<<" jet3_ptRaw = "<<jet3_ptRaw<<std::endl;
//    std::cout<<"mpf_r = "<<mpf_r<<" rel_r = "<<rel_r<<" asymmetry = "<<asymmetry<<" wgt = "<<gen_weight<<" pvN = "<<nvertices<<" wgt = "<<ev_weight<<" true PU = "<<nPU<<std::endl;
    return true;
  }

//     // cout<<"NEW EVENT"<<endl;
//     sel.SetEvent(event);
//     jetcorr.SetEvent(event);
//     mcweight.SetEvent(event);
//     pileupData.SetEvent(event);
// //     datacorr.SetEvent(event);
//     // float j3L1corr =1.;
//     // float j1L1corr =1.;
//     // float j2L1corr =1.;
//     // const baconhep::TEventInfo & info = event.get(h_eventInfo);
//     // baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);

//     const baconhep::TEventInfo & eventInfo = event.get(h_eventInfo);
//     const TClonesArray & js = event.get(h_jets);
//     int nPU_tt = eventInfo.nPU;
//     //    cout<<"nPU = "<<nPU_tt<<endl;
//     event.set(tt_nPU,nPU_tt);

//     //    std::vector<baconhep::TJet> js = event.get(h_jets);
//     // //! JER smearing
//     // if(is_mc){ /// apply for MC only

//     //     //! matching from GEN to RECO
//     //     if(!jetcorr.JetMatching()) return false;
//     //     //! JER smearing
//     //     //0 = central; 1 = scale up; -1 = scale down
//     //     if(!jetcorr.JetResolutionSmearer(0)) return false;

//     // }

//     // baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
//     // baconhep::TJet* jet2 = (baconhep::TJet*)js[1];

//     baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
//     baconhep::TJet* jet2 = (baconhep::TJet*)js[1];
//     Int_t njets = js.GetEntries();
//     //   std::cout<<"Number of jets = "<<js.GetEntries()<<std::endl;
//     //    std::cout<<"evtNum = "<<eventInfo->evtNum<<std::endl;
//     // std::cout<<"evtNum = "<<eventInfo.evtNum<<std::endl;
//     baconhep::TJet* jet3;

//     const TClonesArray & pvs = event.get(h_pv);
//     event.set(tt_jet1_pt,jet1->pt);
//     event.set(tt_jet2_pt,jet2->pt);
//     event.set(tt_jet1_ptRaw,jet1->ptRaw);
//     event.set(tt_jet2_ptRaw,jet2->ptRaw);
//     // event.jet1_pt = jet1->pt;
//     // event.jet2_pt = jet2->pt;
//     // event.jet1_ptRaw = jet1->ptRaw;
//     // event.jet2_ptRaw = jet2->ptRaw;
//     //  std::cout<<"event.jet1_ptRaw = "<<event.get(tt_jet1_ptRaw)<<" event.jet1_pt = "<<event.get(tt_jet1_pt)<<std::endl;
//     //    std::cout<<"event.jet2_ptRaw = "<<event.get(tt_jet2_ptRaw)<<" event.jet2_pt = "<<event.get(tt_jet2_pt)<<std::endl;
//     if (njets > 2) {
//         jet3 = (baconhep::TJet*)js[2];
// 	event.set(tt_jet3_pt, jet3->pt);
// 	event.set(tt_jet3_ptRaw, jet3->ptRaw);
// 	//        event.jet3_pt = jet3->pt;
// 	//        event.jet3_ptRaw = jet3->ptRaw;
//     }
//     else{
//       event.set(tt_jet3_pt, -100);
//       event.set(tt_jet3_ptRaw, -100);
//     }
//     float pt_ave = (event.get(tt_jet1_pt) + event.get(tt_jet2_pt))/2;
//     event.set(tt_pt_ave,pt_ave);
//     //    event.pt_ave = pt_ave;
//     //   std::cout<<"pt_ave = "<<pt_ave<<endl;
//     event.set(tt_gen_pthat,0);//set default values for DATA
//     event.set(tt_gen_weight, 0);//set default values for DATA

//     jets_pt = 0;
//     for(int i=0;i<njets;i++)
//       jets_pt += fabs(((baconhep::TJet*)js[i])->pt);
//     event.set(tt_jets_pt,jets_pt);

//    // //  //!!!NO reweighting for reweighting hists
//     if(is_mc){ /// apply for MC only
//         const baconhep::TGenEventInfo & geninfo = event.get(h_genInfo);
//         baconhep::TGenEventInfo* genInfo= new baconhep::TGenEventInfo(geninfo);
// 	event.set(tt_gen_pthat,genInfo->pthat);
// 	event.set(tt_gen_weight, genInfo->weight);
// 	if ((event.get(tt_pt_ave) - event.get(tt_gen_pthat))/event.get(tt_gen_pthat) > 1) return false;
// 	fMCLumiWeight->process(event);
// 	//	std::cout<<"event.weight "<<event.weight<<std::endl;
// 	if(is_mc_reweight){
//         // event.gen_pthat    = genInfo->pthat;
//         // event.gen_weight   = genInfo->weight;

//         //! Reweighting
//         //event.weight = event.weight * genInfo->weight * mcweight.getPuReweighting("Asympt", 69) * mcweight.getEvReweighting(0, "Asympt", 69);
// 	//	if ((event.get(tt_pt_ave) - event.get(tt_gen_pthat))/event.get(tt_gen_pthat) > 1) return false; //ToDO: include special cut
//          //0 = central; 1 = scale up; -1 = scale down; 99 = no scale(no smearing!)
//          //MC option:  Asympt; or Flat
//          //minBiasXsec for PU:  69; or 80 mb
//         //11 = set 1; 12 = set2 ...
// 	//	cout<<"Before: "<<event.weight<<endl;
// 	//	event.weight = event.weight * mcweight.getPuReweighting("Flat", 58); //TEST
// 	event.weight = event.weight * mcweight.getPuReweighting("Flat", 69); //TEST
// 	//	event.weight = event.weight * mcweight.getEvReweighting(99, "Flat", 69) * mcweight.getPuReweighting("Flat", 69);
// 	//	event.weight = event.weight * mcweight.getPuReweighting("Flat", 80); //TEST

// 	//	event.weight = event.weight * event.get(tt_gen_weight) * mcweight.getPuReweighting("Flat", 58); //ToDo: run it 1st!
// 	//	event.weight = event.weight * event.get(tt_gen_weight) * mcweight.getPuReweighting("Flat", 69); //ToDo: run it 1st!

// 	//	event.weight = event.weight * event.get(tt_gen_weight) * mcweight.getPuReweighting("Flat", 80); //ToDo: run it 1st!
// //	event.weight = event.weight * event.get(tt_gen_weight) * mcweight.getPuReweighting("Flat", 69)* mcweight.getEvReweighting(0, "Flat", 69);
// //	cout<<"After: "<<event.weight<<endl;
// 	// std::cout<<"event.weight = "<<event.weight<<endl;
// 	}
//     }


//     float probejet_eta  = -99.;
//     float probejet_pt   = 0;
//     float probejet_phi  = -99.;
//     float probejet_ptRaw = -99.;

//     float barrel_eta    = -99.;
//     float barrel_pt     = 0;
//     float barrel_phi    = -99.;
//     float barrel_ptRaw  = -99.;

//     float rel_r     = -99.;
//     float mpf_r     = -99.;
//     float asymmetry = -99.;

//     TVector2 pt, met;
//     TVector2* MET = new TVector2(1,1);
//     //    MET->SetMagPhi(eventInfo->pfMET ,eventInfo->pfMETphi);
//     MET->SetMagPhi(eventInfo.pfMET ,eventInfo.pfMETphi);

//     //    met.Set(eventInfo->pfMET * cos(eventInfo->pfMETphi),eventInfo->pfMET * sin(eventInfo->pfMETphi));
//     met.Set(eventInfo.pfMET * cos(eventInfo.pfMETphi),eventInfo.pfMET * sin(eventInfo.pfMETphi));

   

//     if ((fabs(jet1->eta)<s_eta_barr)&&(fabs(jet2->eta)<s_eta_barr)) {
//       int ran = rand();
//       int numb = ran % 2 + 1;
//       //      cout<<"numb = "<<numb<<endl;
//         if(numb==1){
//             probejet_eta = jet2->eta;
//             probejet_pt = event.get(tt_jet2_pt);
//             probejet_phi = jet2->phi;
//             probejet_ptRaw = event.get(tt_jet2_ptRaw);

//             barrel_eta = jet1->eta;
//             barrel_pt = event.get(tt_jet1_pt);
//             barrel_phi = jet1->phi;
//             barrel_ptRaw = event.get(tt_jet1_ptRaw);

//             asymmetry = (event.get(tt_jet2_pt) - event.get(tt_jet1_pt))/(event.get(tt_jet2_pt) + event.get(tt_jet1_pt));
//             rel_r = event.get(tt_jet2_pt) / event.get(tt_jet1_pt);

//             pt.Set(event.get(tt_jet1_pt) * cos(jet1->phi),event.get(tt_jet1_pt) * sin(jet1->phi));
//             mpf_r = 1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py());
//         }
//         if(numb==2){
//             probejet_eta = jet1->eta;
//             probejet_pt = event.get(tt_jet1_pt);
//             probejet_phi = jet1->phi;
//             probejet_ptRaw = event.get(tt_jet1_ptRaw);

//             barrel_eta = jet2->eta;
//             barrel_pt = event.get(tt_jet2_pt);
//             barrel_phi = jet2->phi;
//             barrel_ptRaw = event.get(tt_jet2_ptRaw);

//             asymmetry = (event.get(tt_jet1_pt) - event.get(tt_jet2_pt))/(event.get(tt_jet1_pt) + event.get(tt_jet2_pt));
//             rel_r = event.get(tt_jet1_pt) / event.get(tt_jet2_pt);

//             pt.Set(event.get(tt_jet2_pt) * cos(jet2->phi),event.get(tt_jet2_pt) * sin(jet2->phi));
//             mpf_r = 1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py());

//         }
//     } else if ((fabs(jet1->eta)<s_eta_barr)||(fabs(jet2->eta)<s_eta_barr)){
//         if(fabs(jet1->eta)<s_eta_barr){
//             probejet_eta = jet2->eta;
//             probejet_pt = event.get(tt_jet2_pt);
//             probejet_phi = jet2->phi;
//             probejet_ptRaw = event.get(tt_jet2_ptRaw);

//             barrel_eta = jet1->eta;
//             barrel_pt = event.get(tt_jet1_pt);
//             barrel_phi = jet1->phi;
//             barrel_ptRaw = event.get(tt_jet1_ptRaw);

//             asymmetry = (event.get(tt_jet2_pt) - event.get(tt_jet1_pt))/(event.get(tt_jet2_pt) + event.get(tt_jet1_pt));
//             rel_r = event.get(tt_jet2_pt) / event.get(tt_jet1_pt);

//             pt.Set(event.get(tt_jet1_pt) * cos(jet1->phi),event.get(tt_jet1_pt) * sin(jet1->phi));
//             mpf_r = 1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py());

//         }else{
//             probejet_eta = jet1->eta;
//             probejet_pt = event.get(tt_jet1_pt);
//             probejet_phi = jet1->phi;
//             probejet_ptRaw = event.get(tt_jet1_ptRaw);

//             barrel_eta = jet2->eta;
//             barrel_pt = event.get(tt_jet2_pt);
//             barrel_phi = jet2->phi;
//             barrel_ptRaw = event.get(tt_jet2_ptRaw);

//             asymmetry = (event.get(tt_jet1_pt) - event.get(tt_jet2_pt))/(event.get(tt_jet1_pt) + event.get(tt_jet2_pt));
//             rel_r = event.get(tt_jet1_pt) / event.get(tt_jet2_pt);

//             pt.Set(event.get(tt_jet2_pt) * cos(jet2->phi),event.get(tt_jet2_pt) * sin(jet2->phi));
//             mpf_r = 1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py());
//         }

//     }

//     event.set(tt_probejet_eta,probejet_eta);
//     event.set(tt_probejet_phi,probejet_phi);
//     event.set(tt_probejet_pt,probejet_pt);
//     event.set(tt_probejet_ptRaw,probejet_ptRaw);
//     event.set(tt_barreljet_eta,barrel_eta);
//     event.set(tt_barreljet_phi,barrel_phi);
//     event.set(tt_barreljet_pt,barrel_pt);
//     event.set(tt_barreljet_ptRaw,barrel_ptRaw);
//     event.set(tt_asymmetry,asymmetry);
//     event.set(tt_rel_r,rel_r);
//     event.set(tt_mpf_r,mpf_r);

//     event.set(tt_nvertices,pvs.GetEntries());
//     //    event.nvertices = pvs.GetEntries();
// //     if ((event.nvertices < 14.) || (event.nvertices >= 16.) ) return false;

//     float alpha = 0.;
//     float alpha_sum = 0.;
//     if (njets > 2) {
//       alpha = (2*(event.get(tt_jet3_pt)))/(event.get(tt_jet1_pt) + event.get(tt_jet2_pt));
//       alpha_sum = (2*(event.get(tt_jets_pt)-(event.get(tt_jet1_pt) + event.get(tt_jet2_pt))))/(event.get(tt_jet1_pt) + event.get(tt_jet2_pt));
//     }
//     // event.alpha = alpha;
//     event.set(tt_alpha,alpha);
//     event.set(tt_alpha_sum,alpha_sum);

//     if(!sel.DiJet()) return false;

//     h_nocuts->fill(event);
//     //  if(js.GetEntries()>25) return false; //TEST cut events with too high jet multiplicity

//     if(!sel.DiJetAdvanced(event)) return false;

//     h_dijet->fill(event);


//     h_match->fill(event);

//     if(is_data){
//         if(!sel.Trigger(event)) return false;
//     }
//     else
//       if(!sel.PtMC(event)) return false;

//     // if( event.get(tt_jet3_pt) 50.) return false;//27.04.2016: add cut on the 3rd jet

//     //cout<<"Fill hist for selection"<<endl;
//     //if (event.get(tt_alpha) < 0.2) {
//     if (event.get(tt_alpha) < 0.3) { //18.02.2016: change nominal alpha cut to 0.3
//       h_sel->fill(event);
//     }

//     // cout<<"Yff! "<<event.weight<<endl;
//     event.set(tt_ev_weight,event.weight);
//     return true;
//   }


  // as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
  // make sure the ExampleModule is found by class name. This is ensured by this macro:
  UHH2_REGISTER_ANALYSIS_MODULE(TestModule)


