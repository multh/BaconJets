#include <iostream>
#include <memory>
#include <stdlib.h>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
//#include "UHH2/core/include/EventHelper.h"
#include "../include/JECAnalysisHists.h"
#include "../include/JECRunnumberHists.h"

#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/LumiSelection.h> //includes also LuminosityHists.h
#include <UHH2/common/include/TriggerSelection.h>
#include "UHH2/common/include/CleaningModules.h"

#include "UHH2/BaconJets/include/selection.h"
#include "UHH2/BaconJets/include/jet_corrections.h"
#include "UHH2/BaconJets/include/mc_weight.h"
#include "UHH2/BaconJets/include/constants.h"
#include "UHH2/BaconJets/include/TSetTree.h"
//#include "UHH2/BaconJets/include/dijet_event.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/common/include/PrintingModules.h"


//#include "UHH2/BaconJets/include/pileup_data.h"
//#include "UHH2/BaconJets/include/data_corrections.h"
#include "TClonesArray.h"
#include "TString.h"


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
    std::unique_ptr<JetCorrector> jet_corrector, jet_corrector_BCD, jet_corrector_E, jet_corrector_Fearly, jet_corrector_FlateGH;
    std::unique_ptr<JetCorrector> jet_corrector_noRes, jet_corrector_noRes_BCD, jet_corrector_noRes_E, jet_corrector_noRes_Fearly, jet_corrector_noRes_FlateGH;
    std::unique_ptr<GenericJetResolutionSmearer> jetER_smearer; 
    std::unique_ptr<JetLeptonCleaner> jetleptoncleaner, JLC_BCD, JLC_E, JLC_Fearly, JLC_FlateGH;
    std::unique_ptr<JetLeptonCleaner> jetleptoncleaner_noRes, JLC_noRes_BCD, JLC_noRes_E, JLC_noRes_Fearly, JLC_noRes_FlateGH; 
    std::unique_ptr<JetCleaner> jetcleaner;
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
    unique_ptr<AnalysisModule>  Jet_printer, GenParticles_printer;

   


 

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
    Event::Handle<int> tt_partonFlavor; //only MC
    Event::Handle<int> tt_flavorBarreljet, tt_flavorProbejet, tt_flavorLeadingjet, tt_flavorSubleadingjet; //only MC
    Event::Handle<float> tt_response_leadingjet;
    Event::Handle<float> tt_had_n_Efrac, tt_had_ch_Efrac, tt_mu_Efrac, tt_ph_Efrac;
    Event::Handle<float> tt_inst_lumi, tt_integrated_lumi_in_bin, tt_integrated_lumi;
    Event::Handle<int> tt_lumibin;

 
    std::unique_ptr<JECAnalysisHists> h_nocuts, h_sel, h_dijet, h_match, h_final;
    std::unique_ptr<JECAnalysisHists> h_trg40, h_trg60, h_trg80, h_trg140, h_trg200,h_trg260,h_trg320,h_trg400,h_trg500;
    std::unique_ptr<JECAnalysisHists> h_trgHF60, h_trgHF80,h_trgHF100, h_trgHF160,h_trgHF220, h_trgHF300;   
    std::unique_ptr<LuminosityHists> h_lumi_nocuts, h_lumi_sel, h_lumi_dijet, h_lumi_match, h_lumi_final;    
    std::unique_ptr<LuminosityHists> h_lumi_Trig40, h_lumi_Trig60, h_lumi_Trig80, h_lumi_Trig140, h_lumi_Trig200, h_lumi_Trig260, h_lumi_Trig320, h_lumi_Trig400, h_lumi_Trig500;
    std::unique_ptr<LuminosityHists> h_lumi_TrigHF60, h_lumi_TrigHF80, h_lumi_TrigHF100, h_lumi_TrigHF160, h_lumi_TrigHF220, h_lumi_TrigHF300;
    std::unique_ptr<JECRunnumberHists> h_runnr_input;
    uhh2bacon::Selection sel;

    bool debug;
    bool isMC, split_JEC, ClosureTest;
    string jetLabel;
    TString dataset_version;
    JetId Jet_PFID;
    int n_evt;
    
    std::map<run_lumi, double> rl2lumi;
    std::map<run_lumi, double> rl2intlumi;
    TBranch * brun ;
    TBranch * blumiblock;
    TBranch * bilumi;
    double integrated_lumi;
    vector<run_lumi> upper_binborders_runnrs;
    vector<double> lumi_in_bins;


  };

  TestModule::TestModule(uhh2::Context & ctx) :
    sel(ctx)
  {

    cout << "start" << endl;
    isMC = (ctx.get("dataset_type") == "MC");
    //// COMMON MODULES
    if(!isMC) lumi_sel.reset(new LumiSelection(ctx));
    Jet_PFID = JetPFID(JetPFID::WP_LOOSE);
    jetcleaner.reset(new JetCleaner(ctx, Jet_PFID));

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
    /* //old
    jetLabel = ctx.get("JetLabel");
    std::vector<std::string> JEC_corr;
    if(isMC){
      if(jetLabel == "AK4CHS") JEC_corr = JERFiles::Spring16_25ns_L123_AK4PFchs_MC;
      if(jetLabel == "AK8CHS") JEC_corr = JERFiles::Spring16_25ns_L123_AK4PFchs_MC; //CARE, AK4 corrections applied!
      if(jetLabel == "AK4PUPPI") JEC_corr = JERFiles::Spring16_25ns_L23_AK4PFPuppi_MC;
      if(jetLabel == "AK8PUPPI") JEC_corr = JERFiles::Spring16_25ns_L23_AK4PFPuppi_MC;//CARE, AK4 corrections applied!
    }
    else { 
      if(jetLabel == "AK4CHS") JEC_corr = JERFiles::Spring16_25ns_L123_noRes_AK4PFchs_DATA;
      if(jetLabel == "AK8CHS") JEC_corr = JERFiles::Spring16_25ns_L123_noRes_AK4PFchs_DATA; //CARE, AK4 corrections applied!
      if(jetLabel == "AK4PUPPI") JEC_corr = JERFiles::Spring16_25ns_L23_noRes_AK4PFPuppi_DATA;
      if(jetLabel == "AK8PUPPI") JEC_corr = JERFiles::Spring16_25ns_L23_noRes_AK4PFPuppi_DATA;//CARE, AK4 corrections applied!
    }
    jet_corrector.reset(new JetCorrector(ctx, JEC_corr));
    if(isMC) jetER_smearer.reset(new GenericJetResolutionSmearer(ctx)); 
    jetleptoncleaner.reset(new JetLeptonCleaner(ctx, JEC_corr));
*/

    //new
    jetLabel = ctx.get("JetLabel");
    dataset_version = ctx.get("dataset_version");
    split_JEC = (ctx.get("Split_JEC") == "true");
    ClosureTest = (ctx.get("ClosureTest") == "true");
    std::vector<std::string> JEC_corr_noRes, JEC_corr_noRes_BCD, JEC_corr_noRes_E, JEC_corr_noRes_Fearly, JEC_corr_noRes_FlateGH; 
    std::vector<std::string> JEC_corr,       JEC_corr_BCD, JEC_corr_E, JEC_corr_Fearly, JEC_corr_FlateGH, JEC_corr_MC_FlateGH;
    if(isMC){
      //for MC
      if(jetLabel == "AK4CHS"){
	JEC_corr       = JERFiles::Spring16_25ns_V8_G_L123_AK4PFchs_MC; //noRes only for DATA ;)
	//JEC_corr       = JERFiles::Spring16_25ns_V8_G_L1RC23_AK4PFchs_MC;
	JEC_corr_noRes = JERFiles::Spring16_25ns_V8_G_L123_AK4PFchs_MC; 
      }
      if(jetLabel == "AK8CHS"){
	JEC_corr       = JERFiles::Spring16_25ns_V8_G_L123_AK8PFchs_MC; 
	JEC_corr_noRes = JERFiles::Spring16_25ns_V8_G_L123_AK8PFchs_MC; 
      }
      if(jetLabel == "AK4PUPPI"){
	JEC_corr       = JERFiles::Spring16_25ns_V8_G_L23_AK4PFPuppi_MC;
	JEC_corr_noRes = JERFiles::Spring16_25ns_V8_G_L23_AK4PFPuppi_MC;
      }
      if(jetLabel == "AK8PUPPI"){
	JEC_corr       = JERFiles::Spring16_25ns_V8_G_L23_AK8PFPuppi_MC;
	JEC_corr_noRes = JERFiles::Spring16_25ns_V8_G_L23_AK8PFPuppi_MC;
      }
    }
    else { 
      //for calculating residual corrections DATA
      if(jetLabel == "AK4CHS"){
	//residual
	JEC_corr_noRes         = JERFiles::Spring16_25ns_V8_G_L123_noRes_AK4PFchs_DATA;
	JEC_corr_noRes_BCD     = JERFiles::Spring16_25ns_V8_BCD_L123_noRes_AK4PFchs_DATA;
	JEC_corr_noRes_E       = JERFiles::Spring16_25ns_V8_E_L123_noRes_AK4PFchs_DATA;
	JEC_corr_noRes_Fearly  = JERFiles::Spring16_25ns_V8_F_L123_noRes_AK4PFchs_DATA;
	JEC_corr_noRes_FlateGH = JERFiles::Spring16_25ns_V8_G_L123_noRes_AK4PFchs_DATA;
	//closure
	JEC_corr         = JERFiles::Spring16_25ns_V8_G_L123_AK4PFchs_DATA;
	JEC_corr_BCD     = JERFiles::Spring16_25ns_V8_BCD_L123_AK4PFchs_DATA;
	JEC_corr_E       = JERFiles::Spring16_25ns_V8_E_L123_AK4PFchs_DATA;
	JEC_corr_Fearly  = JERFiles::Spring16_25ns_V8_F_L123_AK4PFchs_DATA;
	JEC_corr_FlateGH = JERFiles::Spring16_25ns_V8_G_L123_AK4PFchs_DATA;
	//JEC_corr_FlateGH = JERFiles::Spring16_25ns_V8_G_L1RC23_AK4PFchs_DATA;
      }
      if(jetLabel == "AK8CHS"){
	//residual
	JEC_corr_noRes         = JERFiles::Spring16_25ns_V8_G_L123_noRes_AK8PFchs_DATA; 
	JEC_corr_noRes_BCD     = JERFiles::Spring16_25ns_V8_BCD_L123_noRes_AK8PFchs_DATA;
 	JEC_corr_noRes_E       = JERFiles::Spring16_25ns_V8_E_L123_noRes_AK8PFchs_DATA; 
	JEC_corr_noRes_Fearly  = JERFiles::Spring16_25ns_V8_F_L123_noRes_AK8PFchs_DATA; 
	JEC_corr_noRes_FlateGH = JERFiles::Spring16_25ns_V8_G_L123_noRes_AK8PFchs_DATA; 
	JEC_corr         = JERFiles::Spring16_25ns_V8_G_L123_AK8PFchs_DATA; 
	JEC_corr_BCD     = JERFiles::Spring16_25ns_V8_BCD_L123_AK8PFchs_DATA;
 	JEC_corr_E       = JERFiles::Spring16_25ns_V8_E_L123_AK8PFchs_DATA; 
	JEC_corr_Fearly  = JERFiles::Spring16_25ns_V8_F_L123_AK8PFchs_DATA; 
	JEC_corr_FlateGH = JERFiles::Spring16_25ns_V8_G_L123_AK8PFchs_DATA; 
      }
      if(jetLabel == "AK4PUPPI"){
	//residual
	JEC_corr_noRes         = JERFiles::Spring16_25ns_V8_G_L23_noRes_AK4PFPuppi_DATA;
	JEC_corr_noRes_BCD     = JERFiles::Spring16_25ns_V8_BCD_L23_noRes_AK4PFPuppi_DATA;
	JEC_corr_noRes_E       = JERFiles::Spring16_25ns_V8_E_L23_noRes_AK4PFPuppi_DATA;
	JEC_corr_noRes_Fearly  = JERFiles::Spring16_25ns_V8_F_L23_noRes_AK4PFPuppi_DATA;
	JEC_corr_noRes_FlateGH = JERFiles::Spring16_25ns_V8_G_L23_noRes_AK4PFPuppi_DATA;
	//closure
	JEC_corr         = JERFiles::Spring16_25ns_V8_G_L23_AK4PFPuppi_DATA;
	JEC_corr_BCD     = JERFiles::Spring16_25ns_V8_BCD_L23_AK4PFPuppi_DATA;
	JEC_corr_E       = JERFiles::Spring16_25ns_V8_E_L23_AK4PFPuppi_DATA;
	JEC_corr_Fearly  = JERFiles::Spring16_25ns_V8_F_L23_AK4PFPuppi_DATA;
	JEC_corr_FlateGH = JERFiles::Spring16_25ns_V8_G_L23_AK4PFPuppi_DATA;
      }
      if(jetLabel == "AK8PUPPI"){
	//residual
	JEC_corr_noRes         = JERFiles::Spring16_25ns_V8_G_L23_noRes_AK8PFPuppi_DATA;
	JEC_corr_noRes_BCD     = JERFiles::Spring16_25ns_V8_BCD_L23_noRes_AK8PFPuppi_DATA;
	JEC_corr_noRes_E       = JERFiles::Spring16_25ns_V8_E_L23_noRes_AK8PFPuppi_DATA;
	JEC_corr_noRes_Fearly  = JERFiles::Spring16_25ns_V8_F_L23_noRes_AK8PFPuppi_DATA;
	JEC_corr_noRes_FlateGH = JERFiles::Spring16_25ns_V8_G_L23_noRes_AK8PFPuppi_DATA;
	//closure
	JEC_corr         = JERFiles::Spring16_25ns_V8_G_L23_AK8PFPuppi_DATA;
	JEC_corr_BCD     = JERFiles::Spring16_25ns_V8_BCD_L23_AK8PFPuppi_DATA;
	JEC_corr_E       = JERFiles::Spring16_25ns_V8_E_L23_AK8PFPuppi_DATA;
	JEC_corr_Fearly  = JERFiles::Spring16_25ns_V8_F_L23_AK8PFPuppi_DATA;
	JEC_corr_FlateGH = JERFiles::Spring16_25ns_V8_G_L23_AK8PFPuppi_DATA;
      }
    }
    





    
  
      //for closure test
      if(ClosureTest){
	if(!isMC && split_JEC){ //these only exist for DATA
	  jet_corrector_BCD.reset(new JetCorrector(ctx, JEC_corr_BCD));
	  jet_corrector_E.reset(new JetCorrector(ctx, JEC_corr_E));
	  jet_corrector_Fearly.reset(new JetCorrector(ctx, JEC_corr_Fearly));
	  jet_corrector_FlateGH.reset(new JetCorrector(ctx, JEC_corr_FlateGH));
	  JLC_BCD.reset(new JetLeptonCleaner(ctx, JEC_corr_BCD));
	  JLC_E.reset(new JetLeptonCleaner(ctx, JEC_corr_E));
	  JLC_Fearly.reset(new JetLeptonCleaner(ctx, JEC_corr_Fearly));
	  JLC_FlateGH.reset(new JetLeptonCleaner(ctx, JEC_corr_FlateGH));
	}
	jet_corrector.reset(new JetCorrector(ctx, JEC_corr));
	jetleptoncleaner.reset(new JetLeptonCleaner(ctx, JEC_corr));
      }
      //for residuals
      else{
	if(!isMC && split_JEC){
	  jet_corrector_noRes_BCD.reset(new JetCorrector(ctx, JEC_corr_noRes_BCD));
	  jet_corrector_noRes_E.reset(new JetCorrector(ctx, JEC_corr_noRes_E));
	  jet_corrector_noRes_Fearly.reset(new JetCorrector(ctx, JEC_corr_noRes_Fearly));
	  jet_corrector_noRes_FlateGH.reset(new JetCorrector(ctx, JEC_corr_noRes_FlateGH));
	  JLC_noRes_BCD.reset(new JetLeptonCleaner(ctx, JEC_corr_noRes_BCD));
	  JLC_noRes_E.reset(new JetLeptonCleaner(ctx, JEC_corr_noRes_E));
	  JLC_noRes_Fearly.reset(new JetLeptonCleaner(ctx, JEC_corr_noRes_Fearly));
	  JLC_noRes_FlateGH.reset(new JetLeptonCleaner(ctx, JEC_corr_noRes_FlateGH));
	}
	jet_corrector_noRes.reset(new JetCorrector(ctx, JEC_corr_noRes));
	jetleptoncleaner_noRes.reset(new JetLeptonCleaner(ctx, JEC_corr_noRes));
      }
    

    if(isMC) jetER_smearer.reset(new GenericJetResolutionSmearer(ctx)); 


    //output
    ctx.undeclare_all_event_output();   
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
    tt_partonFlavor = ctx.declare_event_output<int>("partonFlavor");
    tt_flavorBarreljet = ctx.declare_event_output<int>("flavorBarreljet");
    tt_flavorProbejet = ctx.declare_event_output<int>("flavorProbejet");
    tt_flavorLeadingjet = ctx.declare_event_output<int>("flavorLeadingjet");
    tt_flavorSubleadingjet = ctx.declare_event_output<int>("flavorSubleadingjet");
    tt_response_leadingjet = ctx.declare_event_output<float>("leadingjet_response");
    tt_had_n_Efrac = ctx.declare_event_output<float>("neutralhad_Efraction");
    tt_had_ch_Efrac = ctx.declare_event_output<float>("chargedhad_Efraction");
    tt_mu_Efrac = ctx.declare_event_output<float>("mu_Efraction");
    tt_ph_Efrac = ctx.declare_event_output<float>("photon_Efraction");
    tt_inst_lumi = ctx.declare_event_output<float>("instantaneous_lumi");
    tt_integrated_lumi_in_bin = ctx.declare_event_output<float>("integrated_lumi_in_bin");
    tt_lumibin = ctx.declare_event_output<int>("lumibin");
    tt_integrated_lumi = ctx.declare_event_output<float>("integrated_lumi");




    h_runnr_input.reset(new JECRunnumberHists(ctx,"Runnr_input"));
    h_nocuts.reset(new JECAnalysisHists(ctx,"NoCuts"));
    h_dijet.reset(new JECAnalysisHists(ctx,"diJet"));
    h_match.reset(new JECAnalysisHists(ctx,"JetMatching"));
    h_sel.reset(new JECAnalysisHists(ctx,"Selection"));
    h_final.reset(new JECAnalysisHists(ctx,"Final"));

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

    h_lumi_nocuts.reset(new LuminosityHists(ctx,"Lumi_noCuts"));  
    h_lumi_sel.reset(new LuminosityHists(ctx,"Lumi_Selection"));
    h_lumi_dijet.reset(new LuminosityHists(ctx,"Lumi_diJet"));
    h_lumi_match.reset(new LuminosityHists(ctx,"Lumi_JetMatching"));
    h_lumi_final.reset(new LuminosityHists(ctx,"Lumi_Final"));
    h_lumi_Trig40.reset(new LuminosityHists(ctx,"Lumi_Trig40"));  
    h_lumi_Trig60.reset(new LuminosityHists(ctx,"Lumi_Trig60")); 
    h_lumi_Trig80.reset(new LuminosityHists(ctx,"Lumi_Trig80")); 
    h_lumi_Trig140.reset(new LuminosityHists(ctx,"Lumi_Trig140")); 
    h_lumi_Trig200.reset(new LuminosityHists(ctx,"Lumi_Trig200")); 
    h_lumi_Trig260.reset(new LuminosityHists(ctx,"Lumi_Trig260")); 
    h_lumi_Trig320.reset(new LuminosityHists(ctx,"Lumi_Trig320")); 
    h_lumi_Trig400.reset(new LuminosityHists(ctx,"Lumi_Trig400")); 
    h_lumi_Trig500.reset(new LuminosityHists(ctx,"Lumi_Trig500")); 
    h_lumi_TrigHF60.reset(new LuminosityHists(ctx,"Lumi_TrigHF60")); 
    h_lumi_TrigHF80.reset(new LuminosityHists(ctx,"Lumi_TrigHF80")); 
    h_lumi_TrigHF100.reset(new LuminosityHists(ctx,"Lumi_TrigHF100")); 
    h_lumi_TrigHF160.reset(new LuminosityHists(ctx,"Lumi_TrigHF160")); 
    h_lumi_TrigHF220.reset(new LuminosityHists(ctx,"Lumi_TrigHF220")); 
    h_lumi_TrigHF300.reset(new LuminosityHists(ctx,"Lumi_TrigHF300")); 
    
    Jet_printer.reset(new JetPrinter("Jet-Printer", 0));
    GenParticles_printer.reset(new GenParticlesPrinter(ctx));

    debug = false;
    n_evt = 0;
    
    string lumifile = ctx.get("lumi_file");
    std::unique_ptr<TFile> file(TFile::Open(lumifile.c_str(), "read"));
    TTree * tree = dynamic_cast<TTree*>(file->Get("AnalysisTree"));
    if(!tree){
      throw runtime_error("LuminosityHists: Did not find TTree 'AnalysisTree' in file ;" + lumifile + "'");
    }
    // only fetch branches we really need:
    brun = tree->GetBranch("run");
    blumiblock = tree->GetBranch("luminosityBlock");
    bilumi = tree->GetBranch("intgRecLumi");



    run_lumi rl;
    double ilumi;
    double intlumi_pb = 0;
    brun->SetAddress(&rl.run);
    blumiblock->SetAddress(&rl.lumiblock);
    bilumi->SetAddress(&ilumi);

    //loop over all lumiblocks to save the map between run/lumiblock and stored lumi of the lumiblock (to be divided by 23s)
    auto ientries = tree->GetEntries();
    for(auto ientry = 0l; ientry < ientries; ientry++){
      for(auto b : {brun, blumiblock, bilumi}){
	b->GetEntry(ientry);
      }
      double ilumi_pb = ilumi * 1e-6; // convert units in file (microbarn) to pb.
      intlumi_pb += ilumi_pb;
      rl2lumi.insert(make_pair(rl, ilumi_pb));
      rl2intlumi.insert(make_pair(rl, intlumi_pb));
    }
   

    double ilumi_current_bin = 0.0;
    run_lumi last_entry;
    for(const auto & rl : rl2lumi){
      ilumi_current_bin += rl.second;
      if(ilumi_current_bin >= 2000){
	upper_binborders_runnrs.push_back(rl.first);
	lumi_in_bins.push_back(ilumi_current_bin - rl.second);
	ilumi_current_bin = ilumi_current_bin - 2000;
      }
      last_entry = rl.first;
    }
    upper_binborders_runnrs.push_back(last_entry); //this is not exactly an UPPER limit because it is equal to the highest possible entry, not greater than it...created exception for this case.
    lumi_in_bins.push_back(ilumi_current_bin);
    


  };





  TestModule::~TestModule() {

  }

  bool TestModule::process(Event & event) {
    n_evt++;

    if(!isMC){ //split up RunF into Fearly and Flate (the latter has to be hadd'ed to RunG manually)
      if(dataset_version.Contains("Fearly")){
	if(event.run >= s_runnr_Fearly) return false;
      }
      else if(dataset_version.Contains("Flate")){
	if(event.run < s_runnr_Fearly) return false;
      }
    }

    h_runnr_input->fill(event);

    /* CMS-certified luminosity sections */
    if(event.isRealData){
      if(!lumi_sel->passes(event)){
	return false;
      }
    }

    //h_input->fill(event);

    int event_in_lumibin = -1;
    double fill_event_integrated_lumi = 0;
    double inst_lumi = -1;
    double int_lumi_event = -1;
    if(event.isRealData){
      run_lumi rl_event{event.run, event.luminosityBlock};
      double lumiblock_lumi = rl2lumi[rl_event];
      inst_lumi = lumiblock_lumi/23;
      int_lumi_event = rl2intlumi[rl_event];

      vector<run_lumi>::iterator it;
      if(!(rl_event < upper_binborders_runnrs.back())){
	if(upper_binborders_runnrs.back() < rl_event) throw runtime_error("TestModule: run_lumi of event greater than largest bin-border.");
	else it = prev(upper_binborders_runnrs.end()); //force the entries with the highest run_lumi to enter the last bin instead of overflow.
      }
      else it = upper_bound(upper_binborders_runnrs.begin(), upper_binborders_runnrs.end(), rl_event); //find the first entry in the vector of binborders that is greater than rl_event
      
      event_in_lumibin = distance(upper_binborders_runnrs.begin(), it); //find how many elements of the vector of binborders are smaller than 'it', this is the bin to be filled
      fill_event_integrated_lumi = lumi_in_bins.at(event_in_lumibin);
    }
    
 
    event.weight = 1;

    int n_jets_beforeCleaner = event.jets->size();
    //JetID
    if(jetLabel == "AK4CHS" || jetLabel == "AK8CHS") jetcleaner->process(event);
    int n_jets_afterCleaner = event.jets->size();
    //discard events if not all jets fulfill JetID instead of just discarding single jets
    if (n_jets_beforeCleaner != n_jets_afterCleaner) return false;
    sort_by_pt<Jet>(*event.jets);
    //h_cleaner->fill(event);

    const int jet_n = event.jets->size();
    if(jet_n<2) return false;
    //h_2jets->fill(event);

    /* //old
    jetleptoncleaner->process(event);
    jet_corrector->process(event);
    */
    //semi-new
    /*
    if(ClosureTest){
      if(jetLabel == "AK4PUPPI" && event.run < s_runnr_Fearly && !isMC) throw runtime_error("TestModule.cxx: Closure test for AK4Puppi collection only valid for runG data. Here, runnr < 278802.");
      if(isMC){
	JLC_MC_FlateGH->process(event);
	jet_corrector_MC_FlateGH->process(event);
      }
      else{
	//apply different corrections for events in different run intervals
	if(event.run <= s_runnr_BCD){
	  JLC_BCD->process(event);
	  jet_corrector_BCD->process(event);
	}
	else if(event.run <= s_runnr_E){
	  JLC_E->process(event);
	  jet_corrector_E->process(event);
	}
	else if(event.run < s_runnr_Fearly){ //< is correct, not <=
	  JLC_Fearly->process(event);
	  jet_corrector_Fearly->process(event); 
	}
	else if(event.run >= s_runnr_Fearly){
	  JLC_FlateGH->process(event);
	  jet_corrector_FlateGH->process(event);
	}
	else throw runtime_error("TestModule.cxx: run number not covered by if-statements in process-routine.");
      }
    }
    else{
      jetleptoncleaner_noRes->process(event);
      jet_corrector_noRes->process(event);
    }
*/ 
    if(ClosureTest){
      if(!isMC && split_JEC){ //these only exist for DATA
	if(event.run <= s_runnr_BCD){
	  jet_corrector_BCD->process(event);
	  JLC_BCD->process(event);
	}
	else if(event.run <= s_runnr_E){
	  jet_corrector_E->process(event);
	  JLC_E->process(event);
	}
	else if(event.run < s_runnr_Fearly){ //< is correct, not <=
	  jet_corrector_Fearly->process(event);
	  JLC_Fearly->process(event);
	}
	else if(event.run >= s_runnr_Fearly){
	  jet_corrector_FlateGH->process(event);
	  JLC_FlateGH->process(event);
	}
	else throw runtime_error("TestModule.cxx: run number not covered by if-statements in process-routine.");
      }
      else{
	jet_corrector->process(event);
	jetleptoncleaner->process(event);
      }
    }
    else{
      if(!isMC && split_JEC){ //these only exist for DATA
	if(event.run <= s_runnr_BCD){
	  jet_corrector_noRes_BCD->process(event);
	  JLC_noRes_BCD->process(event);
	}
	else if(event.run <= s_runnr_E){
	  jet_corrector_noRes_E->process(event);
	  JLC_noRes_E->process(event);
	}
	else if(event.run < s_runnr_Fearly){ //< is correct, not <=
	  jet_corrector_noRes_Fearly->process(event);
	  JLC_noRes_Fearly->process(event);
	}
	else if(event.run >= s_runnr_Fearly){
	  jet_corrector_noRes_FlateGH->process(event);
	  JLC_noRes_FlateGH->process(event);
	}
	else throw runtime_error("TestModule.cxx: run number not covered by if-statements in process-routine.");
      }
      else{
	jet_corrector_noRes->process(event);
	jetleptoncleaner_noRes->process(event);
      }
    }

    if(jetLabel == "AK4CHS" ) {
      if(jetER_smearer.get()) jetER_smearer->process(event);  
    }
    //h_corr->fill(event);

    // if(jet_n>2) {
    //   if (event.jets->at(2).pt()<15) return false;
    //   //h_thirdPt15->fill(event);
    // }



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
    //double trg_thresh[9] = {56,78,100,168,232,300,366,453,562}; //2015
    //double trgHF_thresh[6] = {77,131,154,244,321,426}; //2015
    double trg_thresh[9] = {s_Pt_Ave40_cut,s_Pt_Ave60_cut,s_Pt_Ave80_cut,s_Pt_Ave140_cut,s_Pt_Ave200_cut,s_Pt_Ave260_cut,s_Pt_Ave320_cut,s_Pt_Ave400_cut,s_Pt_Ave500_cut}; 
    double trgHF_thresh[6] = {s_Pt_Ave60HF_cut,s_Pt_Ave80HF_cut,s_Pt_Ave100HF_cut,s_Pt_Ave160HF_cut,s_Pt_Ave220HF_cut,s_Pt_Ave300HF_cut}; 

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
      const bool pass_trigger = (pass_trigger40 || pass_trigger60 || pass_trigger80 || pass_trigger140 || pass_trigger200 
				 || pass_trigger260 || pass_trigger320 || pass_trigger400 || pass_trigger500
				 || pass_trigger60_HFJEC || pass_trigger80_HFJEC || pass_trigger100_HFJEC
				 || pass_trigger160_HFJEC || pass_trigger220_HFJEC || pass_trigger300_HFJEC);
      if(debug){
	cout << "before triggers: " << endl;
	cout << " Evt# "<<event.event<<" Run: "<<event.run<<" " << endl;
      }

      if(!pass_trigger)
	return false;
    }

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
    } 
    else if ((fabs(jet1->eta())<s_eta_barr)||(fabs(jet2->eta())<s_eta_barr)){
      if(fabs(jet1->eta())<s_eta_barr){
	jet_probe = jet2;
	jet_barrel = jet1;
      }
      else{
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
    if(!event.isRealData) nPU = event.genInfo->pileup_TrueNumInteractions();

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

    int flavor = 0;

    double had_n_Efrac = event.jets->at(0).neutralHadronEnergyFraction();
    double had_ch_Efrac = event.jets->at(0).chargedHadronEnergyFraction();
    double mu_Efrac = event.jets->at(0).muonEnergyFraction();
    double ph_Efrac = event.jets->at(0).photonEnergyFraction();

    //fill the containers
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
    event.set(tt_partonFlavor,flavor); 
    event.set(tt_had_n_Efrac,had_n_Efrac);
    event.set(tt_had_ch_Efrac,had_ch_Efrac);
    event.set(tt_mu_Efrac,mu_Efrac);    
    event.set(tt_ph_Efrac,ph_Efrac); 
    event.set(tt_inst_lumi,inst_lumi);
    event.set(tt_integrated_lumi_in_bin,fill_event_integrated_lumi);
    event.set(tt_lumibin,event_in_lumibin);
    event.set(tt_integrated_lumi,int_lumi_event);

    sel.SetEvent(event);
    //good primary vertex
    int nGoodVts = sel.goodPVertex();

    if(debug){
      cout << "before good vertex selection : " << endl;
      cout << " Evt# "<<event.event<<" Run: "<<event.run<<" " << endl;
    }

    if(nGoodVts<=0) return false;
    event.set(tt_nGoodvertices, nGoodVts);

    if(debug){
      cout << "before 'dijet selection' : " << endl;
      cout << " Evt# "<<event.event<<" Run: "<<event.run<<" " << endl;
    }

    if(!sel.DiJet()) return false;
    h_nocuts->fill(event);
    h_lumi_nocuts->fill(event);

   if(debug){
     cout << "before 'dijet advanced selection' : " << endl;
     cout << " Evt# "<<event.event<<" Run: "<<event.run<<" " << endl;
   }

    if(!sel.DiJetAdvanced(event)) return false;
    h_dijet->fill(event);
    h_lumi_dijet->fill(event);
    h_match->fill(event);
    h_lumi_match->fill(event);
    if(event.isRealData){
      if(pass_trigger40) {h_trg40->fill(event); h_lumi_Trig40->fill(event);}
      if(pass_trigger60) {h_trg60->fill(event); h_lumi_Trig60->fill(event);} 
      if(pass_trigger80) {h_trg80->fill(event); h_lumi_Trig80->fill(event);}
      if(pass_trigger140) {h_trg140->fill(event); h_lumi_Trig140->fill(event);}
      if(pass_trigger200) {h_trg200->fill(event); h_lumi_Trig200->fill(event);}
      if(pass_trigger260) {h_trg260->fill(event); h_lumi_Trig260->fill(event);}
      if(pass_trigger320) {h_trg320->fill(event); h_lumi_Trig320->fill(event);} 
      if(pass_trigger400) {h_trg400->fill(event); h_lumi_Trig400->fill(event);}
      if(pass_trigger500) {h_trg500->fill(event); h_lumi_Trig500->fill(event);}
      if(pass_trigger60_HFJEC) {h_trgHF60->fill(event); h_lumi_TrigHF60->fill(event);}  
      if(pass_trigger80_HFJEC) {h_trgHF80->fill(event); h_lumi_TrigHF80->fill(event);}
      if(pass_trigger100_HFJEC) {h_trgHF100->fill(event); h_lumi_TrigHF100->fill(event);}
      if(pass_trigger160_HFJEC) {h_trgHF160->fill(event); h_lumi_TrigHF160->fill(event);}
      if(pass_trigger220_HFJEC) {h_trgHF220->fill(event); h_lumi_TrigHF220->fill(event);}
      if(pass_trigger300_HFJEC) {h_trgHF300->fill(event); h_lumi_TrigHF300->fill(event);}
    }
    else{    
      if(debug){
	cout << "before Pt selection (MC only) : " << endl;
	cout << " Evt# "<<event.event<<" Run: "<<event.run<<" " << endl;
      }
      if(!sel.PtMC(event)) return false; // For MC only one Pt threshold
    }
    if (event.get(tt_alpha) < 0.3) {
      h_sel->fill(event);
      h_lumi_sel->fill(event);
    }

    h_final->fill(event);
    h_lumi_final->fill(event);


    if(debug){
      cout<<"-- Event -- "<<endl;
      cout<<" Evt# "<<event.event<<" Run: "<<event.run<<" "<<endl;
      cout<<" Npv = "<<event.get(tt_nvertices)<<" jet_pt_ave = "<<event.get(tt_pt_ave)<<" MET = "<<met.Mod()<<endl;
      cout<<"Probe: "<<event.get(tt_probejet_eta)<<" "<<event.get(tt_probejet_phi)
	  <<" "<<event.get(tt_probejet_pt)<<" "<<event.get(tt_probejet_ptRaw)<<endl;
      cout<<" Barrel: "<<event.get(tt_barreljet_eta)<<" "<<event.get(tt_barreljet_phi)
	  <<" "<<event.get(tt_barreljet_pt)<<" "<<event.get(tt_barreljet_ptRaw)<<endl;
      cout<<" "<<event.get(tt_asymmetry)<<" "<<event.get(tt_rel_r)<<" "<<event.get(tt_mpf_r)<<""<<endl;
      cout<<" "<<endl; 
    }

    if(debug && isMC){
      Jet_printer->process(event);
      GenParticles_printer->process(event);
      cout << "event has " << event.genjets->size() << " GenJets" << endl;
      for(size_t i=0; i< event.genjets->size(); ++i){
        const auto & jet = (*event.genjets)[i];
        cout << " GenJet[" << i << "]: pt=" << jet.pt() << "; eta=" << jet.eta() << "; phi=" << jet.phi() <<  endl;
      }
    }
 
  
    
    if(isMC){    
      double flavor_barreljet = 0;
      double flavor_probejet = 0;
      double flavor_leadingjet = 0;
      double flavor_subleadingjet = 0;
      const unsigned int genjets_n = event.genjets->size();
      int idx_jet_matching_genjet[genjets_n];

      //match genp to gen-jets
      int idx_j=0;
      int idx_genp_min = -1;
      //this array contains one idx for each jet in the event. If -1: unmatched, else: idx of the closest genpart with dR<=0.2
      int idx_matched_jets[jet_n];
      for(int i=0; i<jet_n; i++){
	idx_matched_jets[i] = -1;
      }

      //matching gen- and reco-jets
      for(unsigned int i=0; i<event.genjets->size(); i++){
	double dR_min = 99999; int idx_matching_jet = -1;
	for(unsigned int j=0; j<event.jets->size(); j++){
	  double dR = deltaR(event.jets->at(j), event.genjets->at(i));
	  if(debug) cout << "dR between GenJet " << i << " and RecoJet " << j << ": " << dR << endl;
	  if(dR<dR_min){
	    dR_min = dR; 
	    if(dR_min<0.1) idx_matching_jet = j;
	  }
	}
	idx_jet_matching_genjet[i] = idx_matching_jet;
	if(debug) cout << "the jet matching the genjet no. " << i << " is jet no. " << idx_matching_jet << endl;
      }
      /////////////////////


      for(Particle & genj : *event.genjets){
	double dr_min = 99999;
	double dr_cut = 0;
	if(jetLabel == "AK4CHS" || jetLabel == "AK4PUPPI") dr_cut = 0.2;
	else if (jetLabel == "AK8CHS" || jetLabel == "AK8PUPPI")dr_cut = 0.4;
	else throw runtime_error("TestModule.cxx: Invalid jet-label specified.");

	int idx_g = 0;
	for(GenParticle & genp: *event.genparticles){
	  double dr = deltaR(genj,genp);
	  if(dr < dr_min){
	    dr_min = dr;
	    idx_genp_min = idx_g;	
	  }	
	  //cout << "dr between genjet " << idx_j << " and genp (flavor: " << genp.flavor() << ") " << idx_g << "= " << dr << endl;
	  idx_g++;
	}
	if(dr_min <= dr_cut) {
	  if(debug) cout << "genjet " << idx_j << " is matched to genparticle " << idx_genp_min << " of flavor " << event.genparticles->at(idx_genp_min).flavor() << " within dR = " << dr_min << ". " <<  endl; 
	  if(idx_jet_matching_genjet[idx_j] >= 0) idx_matched_jets[idx_jet_matching_genjet[idx_j]] = idx_genp_min;
	}
	idx_j++;
      }

      //only consider jets that could be matched to a genparticle, these shall take the partons flavor by definition
      //TEST
      if(debug){
	for (int i=0; i<jet_n; i++){
	  if(idx_matched_jets[i] != -1) cout << "Jet no. " << i << " is matching genpart no. " << idx_matched_jets[i] << endl;
	}
      }

      // flavor-quantities

      if(debug && event.genjets->size() <2) cout << "WARNING: GENjets size < 2" << endl;

      //only consider the barreljet, is it leading or sub-leading jet?
      int idx_barreljet = -1;
      if(fabs(jet1->pt() - jet_barrel->pt()) < 0.001) idx_barreljet = 0;
      else if (fabs(jet2->pt() - jet_barrel->pt()) < 0.001) idx_barreljet = 1;
      else throw runtime_error("first two jets are not the barrel jets, how could this happen?");
    
      //obtain flavor of the barreljet
      //-1: unmatched, 0: alpha too large, >0: flavor of matching genparticle 
      if(idx_matched_jets[idx_barreljet] != -1)	flavor_barreljet = fabs(event.genparticles->at(idx_matched_jets[idx_barreljet]).flavor());
      else flavor_barreljet = -1;
      if(debug) cout << "barreljet is jet no. " << idx_barreljet << ", alpha = " << event.get(tt_alpha) << ", flavor of barreljet = " << flavor_barreljet << endl;
    

      //also for probe jets
      int idx_probejet = fabs(idx_barreljet - 1);
      //obtain flavor of the probejet
      //-1: unmatched,  >0: flavor of matching genparticle 
      if(idx_matched_jets[idx_probejet] != -1) flavor_probejet = fabs(event.genparticles->at(idx_matched_jets[idx_probejet]).flavor());
      else flavor_probejet = -1;
      if(debug) cout << "probejet is jet no. " << idx_probejet << ", alpha = " << event.get(tt_alpha) << ", flavor of probejet = " << flavor_probejet << endl;
      
      
      //same for leading jet
      //-1: unmatched, 0: alpha too large, >0: flavor of matching genparticle 
      if(idx_matched_jets[0] != -1) flavor_leadingjet = fabs(event.genparticles->at(idx_matched_jets[0]).flavor());
      else flavor_leadingjet = -1;
      if(debug) cout << "leadingjet is jet no. " << 0 << ", alpha = " << event.get(tt_alpha) << ", flavor of leadingjet = " << flavor_leadingjet << endl;
      

      //same for subleading jet
      //-1: unmatched, 0: alpha too large, >0: flavor of matching genparticle 
      if(idx_matched_jets[1] != -1) flavor_subleadingjet = fabs(event.genparticles->at(idx_matched_jets[1]).flavor());
      else flavor_subleadingjet = -1;
      if(debug) cout << "subleadingjet is jet no. " << 1 << ", alpha = " << event.get(tt_alpha) << ", flavor of subleadingjet = " << flavor_subleadingjet << endl;

      event.set(tt_flavorBarreljet,flavor_barreljet);   
      event.set(tt_flavorProbejet,flavor_probejet);  
      event.set(tt_flavorLeadingjet,flavor_leadingjet);  
      event.set(tt_flavorSubleadingjet,flavor_subleadingjet);  

      //response of leading jet
      //find corresponding genjet
      int idx_corresponding_genjet = -1;
      for(unsigned int i=0; i<event.genjets->size(); i++){
	if(debug) cout << idx_jet_matching_genjet[i] << endl;
	if(idx_jet_matching_genjet[i] == 0) idx_corresponding_genjet = i;
      }
      double response_jet1 = -1;
      if(idx_corresponding_genjet != -1) response_jet1 = event.jets->at(0).pt() / event.genjets->at(idx_corresponding_genjet).pt();



      event.set(tt_response_leadingjet,response_jet1);  
     


    } //isMC

    else{
      event.set(tt_flavorBarreljet,-1);   
      event.set(tt_flavorProbejet,-1);  
      event.set(tt_flavorLeadingjet,-1);  
      event.set(tt_flavorSubleadingjet,-1);  
      event.set(tt_response_leadingjet,-1.);  
    }
    
 


    return true;
  }



  // as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
  // make sure the ExampleModule is found by class name. This is ensured by this macro:
  UHH2_REGISTER_ANALYSIS_MODULE(TestModule)


