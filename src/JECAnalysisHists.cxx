#include "UHH2/BaconJets/include/JECAnalysisHists.h"
#include "UHH2/BaconJets/include/constants.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Jet.h"

#include "UHH2/BaconTrans/baconheaders/TJet.hh"
#include "UHH2/BaconTrans/baconheaders/TEventInfo.hh"
#include "UHH2/BaconTrans/baconheaders/BaconAnaDefs.hh"
#include "UHH2/BaconTrans/baconheaders/TVertex.hh"



#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <getopt.h>
using namespace std;
using namespace uhh2;
using namespace baconhep;
    uhh2::Event::Handle<TClonesArray> h_pv;
JECAnalysisHists::JECAnalysisHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
    // book all histograms here
    // jets
    TH1::SetDefaultSumw2();

    book<TH1F>("N_jets", "N_{jets}", 50, -0.5, 49.5);
    book<TH1F>("pt","p_{T} all jets; p_{T} (GeV)",100,0,1500);
    book<TH1F>("eta","#eta all jets; #eta",100,-5,5);
    book<TH1F>("phi","#phi all jets; #phi",50,-M_PI,M_PI);
    book<TH1F>("MET","MET all jets; MET",400,0,400);

    book<TH1F>("nPu","Number of PU events",60,0,60);
    book<TH1F>("N_PV","Number of PVtx",60,0,60);
    book<TH1F>("N_PV_sel","Number of PVtx after cuts on z, y, x",60,0,60);

    book<TH1F>("weight_histo","weight_histo ",20,0,2);

    book<TH1F>("pt_1","p_{T} jet 1",100,0,600);
    book<TH1F>("eta_1","#eta jet 1",100,-5,5);

    book<TH1F>("pt_2","p_{T} jet 2",100,0,600);
    book<TH1F>("eta_2","#eta jet 2",100,-5,5);

    book<TH1F>("pt_3","p_{T} jet 3",100,0,600);
    book<TH1F>("eta_3","#eta jet 3",100,-5,5);

    book<TH1F>("ptRaw_barrel","p^{Raw}_{T} barrel jet; p_{T}^{Raw,barrel} (GeV)",100,0,600);
    book<TH1F>("ptRaw_probe","p^{Raw}_{T} probe jet; p_{T}^{Raw,probe} (GeV)",100,0,600);
    book<TH1F>("pt_barrel","p_{T} barrel jet; p_{T}^{barrel} (GeV)",100,0,600);
    book<TH1F>("pt_probe","p_{T} probe jet; p_{T}^{probe} (GeV)",100,0,600);
    book<TH1F>("eta_barrel","#eta barrel jet; #eta^{barrel}",100,-5,5);
    book<TH1F>("eta_probe","#eta probe jet; #eta^{probe}",100,-5,5);
    book<TH1F>("eta_probe_pos","#eta probe jet >=0; #eta^{probe}",50,0,5);
    book<TH1F>("eta_probe_neg","#eta probe jet <0; |#eta^{probe}|",50,0,5);
    book<TH1F>("pt_ave","p_{T} ave jet; p_{T}^{ave} (GeV)",600,0,600);
    book<TH1F>("pt_ave_pthat","p_{T} ave jet; p_{T}^{ave} - p_{T}^{hat})/p_{T}^{hat}(GeV)",100,0,600);
    book<TH1F>("pt_ave_rebin","p_{T} ave jet; p_{T}^{ave} (GeV)",300,0,3000);

    book<TH1F>("pt_rel","p_{T}^{jet3} / p_{T}^{ave}; #alpha ", 50, 0, 1);
    book<TH1F>("generic_pt_rel","generic p_{T}^{jet3} / p_{T}^{ave}; #alpha ", 50, 0, 1);

    book<TH1F>("asym","asymmetrie jet 1 and jet 2; Asymmetry",150,-1.5,1.5);
    book<TH1F>("generic_asym","generic asymmetrie jet 1 and jet 2; Asymmetry",100,-1,1);
    // book<TH1F>("mpf","MPF response; MPF response",100,0.5,1.5);
    // book<TH1F>("generic_mpf","generic MPF response; MPF response",100,0.5,1.5);
    book<TH1F>("mpf","MPF response; MPF response",250,0.,2.5);
    book<TH1F>("generic_mpf","generic MPF response; MPF response",250,0.,2.5);
    // book<TH1F>("r_rel","R_{rel}; R_{rel}; Relative response",100,0.5,1.5);
    // book<TH1F>("generic_r_rel","generic R_{rel}; Relative response",100,0.5,1.5);
    book<TH1F>("r_rel","R_{rel}; R_{rel}; Relative response",250,0.,2.5);
    book<TH1F>("generic_r_rel","generic R_{rel}; Relative response",250,0.,2.5);

    book<TH1F>("DeltaPhi_Jet1_Jet2", "#Delta#Phi(first jet, second jet); #Delta #Phi", 100, 0, 7);
    book<TH2F>("ptrel_vs_deltaphi","delta phi vs pt_rel", 50, 0, 1 ,50, 0, 3.14);

    book<TH1F>("pt_ave_hltDiPFJetAve40","p_{T} ave40 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("pt_ave_hltDiPFJetAve60","p_{T} ave60 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("pt_ave_hltDiPFJetAve80","p_{T} ave80 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("pt_ave_hltDiPFJetAve140","p_{T} ave140 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("pt_ave_hltDiPFJetAve200","p_{T} ave200 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("pt_ave_hltDiPFJetAve260","p_{T} ave260 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("pt_ave_hltDiPFJetAve320","p_{T} ave320 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("pt_ave_hltDiPFJetAve400","p_{T} ave400 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("pt_ave_hltDiPFJetAve500","p_{T} ave500 jet; p_{T}^{ave} (GeV)",120,0,600);


    book<TH1F>("HLT_DiPFJetAve100_HFJEC","p_{T} ave100 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("HLT_DiPFJetAve160_HFJEC","p_{T} ave160 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("HLT_DiPFJetAve220_HFJEC","p_{T} ave220 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("HLT_DiPFJetAve300_HFJEC","p_{T} ave300 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("HLT_DiPFJetAve60_HFJEC","p_{T} ave60 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("HLT_DiPFJetAve80_HFJEC","p_{T} ave80 jet; p_{T}^{ave} (GeV)",120,0,600);

    book<TH1F>("probejet_eta_hltDiPFJetAve40","probejet_eta",100,-5,5);
    book<TH1F>("probejet_eta_hltDiPFJetAve60","probejet_eta",100,-5,5);
    book<TH1F>("probejet_eta_hltDiPFJetAve80","probejet_eta",100,-5,5);
    book<TH1F>("probejet_eta_hltDiPFJetAve140","probejet_eta",100,-5,5);
    book<TH1F>("probejet_eta_hltDiPFJetAve200","probejet_eta",100,-5,5);
    book<TH1F>("probejet_eta_hltDiPFJetAve260","probejet_eta",100,-5,5);
    book<TH1F>("probejet_eta_hltDiPFJetAve320","probejet_eta",100,-5,5);
    book<TH1F>("probejet_eta_hltDiPFJetAve400","probejet_eta",100,-5,5);
    book<TH1F>("probejet_eta_hltDiPFJetAve500","probejet_eta",100,-5,5);


    book<TH1F>("probejet_eta_HLT_DiPFJetAve100_HFJEC","probejet_eta",100,-5,5);
    book<TH1F>("probejet_eta_HLT_DiPFJetAve160_HFJEC","probejet_eta",100,-5,5);
    book<TH1F>("probejet_eta_HLT_DiPFJetAve220_HFJEC","probejet_eta",100,-5,5);
    book<TH1F>("probejet_eta_HLT_DiPFJetAve300_HFJEC","probejet_eta",100,-5,5);
    book<TH1F>("probejet_eta_HLT_DiPFJetAve60_HFJEC","probejet_eta",100,-5,5);
    book<TH1F>("probejet_eta_HLT_DiPFJetAve80_HFJEC","probejet_eta",100,-5,5);

    book<TH2D>("Rrel_vs_Npv","Rrel vs. Npv ", 50, 0, 50 ,100, 0.5,1.5);
    book<TH2D>("pt_ave_vs_weight","pt_ave_vs_weight", 100, 0,10 ,10, 0,1);
    book<TH2D>("ptjet3_vs_alpha","pt jet3 vs #alpha", 20, 0, 1., 30,0,300);
    book<TH2D>("pt_ave_vs_alpha","pt ave vs #alpha", 20, 0, 1., 200,0,1000);

    book<TH2D>("mpf_vs_etaProbe","MPF response vs. #eta probe jet; #eta probe; MPF response",100,-5,5,100,0.,2.);
    book<TH2D>("r_rel_vs_etaProbe","Relative response vs. #eta probe jet; #eta probe; R_{rel}",100,-5,5,100,0.,2.);
    book<TH2D>("pt_ave_vs_etaProbe","pt ave vs #eta probe jet; #eta probe; pT_{ave}, GeV",100,-5.2,5.2,200,0,1000);
    book<TH2D>("Rrel_vs_assym","Relative response vs. Asymmetry ; Asymmetry; R_{rel}",100,-1.5,1.5,200,0,10.);
    book<TH2D>("Rmpf_vs_assym","MPF response vs. Asymmetry ; Asymmetry; R_{MPF}",100,-1.5,1.5,200,0,10.);
    book<TH2D>("alpha_vs_assym","#alpha vs. Asymmetry ; Asymmetry; #alpha",100,-1.5,1.5,100,0,1.);
    book<TH2D>("alphaSum_vs_alpha3","#alpha_{sum} vs. #alpha_{3} ; #alpha_{3}; #alpha_{sum};",100,0.,2.,100,0,2.);
    book<TH2D>("Njets_vs_alpha","Number of jets vs. #alpha ; #alpha; Njets",100,0,1.,50, -0.5, 49.5);
    book<TH2D>("Pt_vs_JetN","pT of jets in event ; jet number; pT", 50, -0.5, 49.5, 200,0,500);
    book<TH2D>("eta_vs_JetN","#eta of jets in event ; jet number; #eta", 50, -0.5, 49.5, 100,-5,5);
    // book<TH2D>("diffPt1_vs_alpha3","(pT_{rec}-pT_{gen}) vs. #alpha_{3} ; #alpha_{3}; (pT_{rec}-pT_{gen})_{1st jet} ;",100,0.,2.,100,-200,200.);
    // book<TH2D>("diffPt1_vs_alphaSum","(pT_{rec}-pT_{gen}) vs. #alpha_{sum} ; #alpha_{sum}; (pT_{rec}-pT_{gen})_{1st jet} ;",100,0.,2.,100,-200,200.);
    // book<TH2D>("diffPt2_vs_alpha3","(pT_{rec}-pT_{gen}) vs. #alpha_{3} ; #alpha_{3}; (pT_{rec}-pT_{gen})_{2nd jet} ;",100,0.,2.,100,-200,200.);
    // book<TH2D>("diffPt2_vs_alphaSum","(pT_{rec}-pT_{gen}) vs. #alpha_{sum} ; #alpha_{sum}; (pT_{rec}-pT_{gen})_{2nd jet} ;",100,0.,2.,100,-200,200.);
    // // book<TH2D>("diffPtGen_vs_alpha3","((pT_{1}-pT_{2})_{gen})/pT^{gen}_{ave} vs. #alpha_{3} ; #alpha_{3}; (pT_{1}-pT_{2})_{gen}/pT^{gen}_{ave} ;",100,0.,2.,100,0,2.);
    // // book<TH2D>("diffPtGen_vs_alphaSum","(pT_{1}-pT_{2})_{gen}/pT^{gen}_{ave} vs. #alpha_{sum} ; #alpha_{sum}; (pT_{1}-pT_{2})_{gen}/pT^{gen}_{ave} ;",100,0.,2.,100,0,2.);
    book<TH2D>("diffPtGen_vs_alpha3","(pT_{1}-pT_{2})_{gen}/pT^{gen}_{ave} vs. #alpha_{3} ; #alpha_{3}; (pT_{1}-pT_{2})_{gen}/pT^{gen}_{ave} ;",200,0.,2.,400,-3.,3.);
    book<TH2D>("diffPtGen_vs_alphaSum","(pT_{1}-pT_{2})_{gen}/pT^{gen}_{ave} vs. #alpha_{sum} ; #alpha_{sum}; (pT_{1}-pT_{2})_{gen}/pT^{gen}_{ave} ;",200,0.,2.,400,-3.,3.);

    // book<TH2D>("AbsdiffPtGen_vs_alpha3","((pT_{1}-pT_{2})_{gen})/pT^{gen}_{ave} vs. #alpha_{3} ; #alpha_{3}; |(pT_{1}-pT_{2})_{gen}/pT^{gen}_{ave}|",100,0.,2.,300,0,3.);
    // book<TH2D>("AbsdiffPtGen_vs_alphaSum","(pT_{1}-pT_{2})_{gen}/pT^{gen}_{ave} vs. #alpha_{sum} ; #alpha_{sum}; |(pT_{1}-pT_{2})_{gen}/pT^{gen}_{ave}|",100,0.,2.,300,0,3.);
    // book_TH2F("ptjet3_vs_alpha","pt jet3 vs #alpha", 10, 0, 1., 100,0,1500);
    // book_TH2F("pt_ave_vs_alpha","pt ave vs #alpha", 10, 0, 1., 100,0,1500);
    // ptjet3_vs_alpha = new TH2F("ptjet3_vs_alpha","pt jet3 vs #alpha", 10, 0, 1., 100,0,1500);
    // pt_ave_vs_alpha = new TH2F("pt_ave_vs_alpha","pt ave vs #alpha", 10, 0, 1., 100,0,1500);


//     book<TProfile>("prof_Rrel_vs_Npv","Rrel vs. Npv ", 50, 0, 50 ,100, 0.5,1.5);
//     TProfile * prof_Rrel_vs_Npv  = new TProfile("prof_Rrel_vs_Npv","Rrel vs. Npv",100, 0, 50, 0.5, 1.5);


/*    uhh2::Event::Handle<TClonesArray> h_pv;*/
//    h_jets = ctx.get_handle<TClonesArray>("AK4PFCHS");
    //    h_jets = ctx.get_handle<TClonesArray>("AK4PFPUPPI");
    auto jetCollection = ctx.get("jetCollection");
    h_jets = ctx.declare_event_input<TClonesArray>(jetCollection);
    h_eventInfo = ctx.get_handle<baconhep::TEventInfo>("Info");
    h_pv = ctx.get_handle<TClonesArray>("PV");
    tt_gen_pthat  = ctx.get_handle<float>("gen_pthat");
    tt_gen_weight = ctx.get_handle<float>("gen_weight");
    tt_jet1_pt = ctx.get_handle<float>("jet1_pt");
    tt_jet2_pt = ctx.get_handle<float>("jet2_pt");
    tt_jet3_pt = ctx.get_handle<float>("jet3_pt");
    tt_jet1_ptRaw = ctx.get_handle<float>("jet1_ptRaw");
    tt_jet2_ptRaw = ctx.get_handle<float>("jet2_ptRaw");
    tt_jet3_ptRaw = ctx.get_handle<float>("jet3_ptRaw");
    tt_nvertices = ctx.get_handle<int>("nvertices");
    tt_probejet_eta = ctx.get_handle<float>("probejet_eta");
    tt_probejet_phi = ctx.get_handle<float>("probejet_phi");
    tt_probejet_pt = ctx.get_handle<float>("probejet_pt");
    tt_probejet_ptRaw = ctx.get_handle<float>("probejet_ptRaw");
    tt_barreljet_eta = ctx.get_handle<float>("barreljet_eta");
    tt_barreljet_phi = ctx.get_handle<float>("barreljet_phi");
    tt_barreljet_pt = ctx.get_handle<float>("barreljet_pt");
    tt_barreljet_ptRaw = ctx.get_handle<float>("barreljet_ptRaw");
    tt_pt_ave = ctx.get_handle<float>("pt_ave");
    tt_alpha = ctx.get_handle<float>("alpha");
    tt_alpha_sum = ctx.get_handle<float>("alpha_sum");
    tt_rel_r = ctx.get_handle<float>("rel_r");
    tt_mpf_r = ctx.get_handle<float>("mpf_r");
    tt_asymmetry = ctx.get_handle<float>("asymmetry");
    tt_nPU = ctx.get_handle<int>("nPU");

}

void JECAnalysisHists::fill(const uhh2::Event & ev){
    fill(ev, 0);
}
void JECAnalysisHists::fill(const uhh2::Event & ev, const int rand){
    // fill the histograms. Please note the comments in the header file:
    // 'hist' is used here a lot for simplicity, but it will be rather
    // slow when you have many histograms; therefore, better
    // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
    // Don't forget to always use the weight when filling.

    const TClonesArray & js = ev.get(h_jets);
    // const baconhep::TEventInfo & info = ev.get(h_eventInfo);
    // baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);
    const baconhep::TEventInfo& eventInfo = ev.get(h_eventInfo);

    const TClonesArray & pvs = ev.get(h_pv);

    double weight = ev.weight;
//    //cout << "weight_histo = " <<weight<<endl;
    Int_t njets = js.GetEntries();
    hist("N_jets")->Fill(njets, weight);

    for (int i=0; i<njets; i++){
        baconhep::TJet* jets = (baconhep::TJet*)js[i];
        hist("pt")->Fill(jets->pt, weight);
        hist("eta")->Fill(jets->eta, weight);
        hist("phi")->Fill(jets->phi, weight);
        // hist("MET")->Fill(eventInfo->pfMET, weight);
        // hist("nPu")->Fill(eventInfo->nPU, weight);
        hist("MET")->Fill(eventInfo.pfMET, weight);
        hist("nPu")->Fill(eventInfo.nPU, weight);
	((TH2D*)hist("Pt_vs_JetN"))->Fill(i+1,jets->pt, weight);
	((TH2D*)hist("eta_vs_JetN"))->Fill(i+1,jets->eta, weight);
        hist("weight_histo")->Fill(weight, 1);
    }
//         hist("nPu")->Fill(ev.nPU, weight);//for data only

    float nPrVer = 0;
    Int_t nvertices = pvs.GetEntries();
    //hist("N_PV")->Fill(ev.nvertices, weight);
    hist("N_PV")->Fill(ev.get(tt_nvertices), weight);

    // require in the event that there is at least one reconstructed vertex
    if(nvertices>0) {
        // pick the first (i.e. highest sum pt) verte
        for (int i=0; i<nvertices; i++){
            baconhep::TVertex* vertices = (baconhep::TVertex*)pvs[i];
            // require that the vertex meets certain criteria
            if((fabs(vertices->z) < s_n_Pv_z) && (fabs(vertices->y) < s_n_Pv_xy) && (fabs(vertices->x) < s_n_Pv_xy) ){
                nPrVer++;
            }
        }
    }
    hist("N_PV_sel")->Fill(nPrVer, weight);

    baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
    // hist("pt_1")->Fill(ev.jet1_pt, weight);
    hist("pt_1")->Fill(ev.get(tt_jet1_pt), weight);
    hist("eta_1")->Fill(jet1->eta, weight);

    baconhep::TJet* jet2 = (baconhep::TJet*)js[1];
    //    hist("pt_2")->Fill(ev.jet2_pt, weight);
    hist("pt_2")->Fill(ev.get(tt_jet2_pt), weight);
    hist("eta_2")->Fill(jet2->eta, weight);
//     cout << "ev.gen_pthat "<< ev.gen_pthat<< "ev.pt_ave"<< ev.pt_ave<<"(ev.pt_ave - ev.gen_pthat)/ev.gen_pthat "<<(ev.pt_ave - ev.gen_pthat)/ev.gen_pthat<< endl;
//    float ratio_pt = (ev.pt_ave - ev.gen_pthat)/ev.gen_pthat;
    float ratio_pt = (ev.get(tt_pt_ave) - ev.get(tt_gen_pthat))/ev.get(tt_gen_pthat);
    // hist("pt_ave_vs_weight") ->Fill( ratio_pt, ev.gen_weight);
    hist("pt_ave_vs_weight") ->Fill( ratio_pt, ev.get(tt_gen_weight));
    hist("pt_ave")          ->Fill(ev.get(tt_pt_ave), weight);
    hist("pt_ave_pthat")   ->Fill(ev.get(tt_pt_ave), weight);
    hist("pt_ave_rebin") ->Fill(ev.get(tt_pt_ave), weight);
    hist("ptRaw_barrel")    ->Fill(ev.get(tt_barreljet_ptRaw), weight);
    hist("ptRaw_probe")     ->Fill(ev.get(tt_probejet_ptRaw) , weight);
    hist("pt_barrel")   ->Fill(ev.get(tt_barreljet_pt), weight);
    hist("pt_probe")    ->Fill(ev.get(tt_probejet_pt) , weight);
    hist("eta_barrel")  ->Fill(ev.get(tt_barreljet_eta), weight);
    hist("eta_probe")   ->Fill(ev.get(tt_probejet_eta) , weight);
    if(ev.get(tt_probejet_eta)>=0) hist("eta_probe_pos")->Fill(ev.get(tt_probejet_eta) , weight);
    else hist("eta_probe_neg")->Fill(fabs(ev.get(tt_probejet_eta)) , weight);
    hist("mpf")         ->Fill(ev.get(tt_mpf_r), weight);
    hist("asym")        ->Fill(ev.get(tt_asymmetry), weight);
    hist("r_rel")       ->Fill(ev.get(tt_rel_r), weight);
    ((TH2D*)hist("Rrel_vs_assym"))->Fill(ev.get(tt_asymmetry),ev.get(tt_rel_r), weight);
    ((TH2D*)hist("Rmpf_vs_assym"))->Fill(ev.get(tt_asymmetry),ev.get(tt_mpf_r), weight);
    double deltaPhi = abs(jet1->phi - jet2->phi);
    hist("DeltaPhi_Jet1_Jet2")->Fill(deltaPhi, weight);

//     TString FileName[7] = {"pt_ave_hltDiPFJetAve40","pt_ave_hltDiPFJetAve80", "pt_ave_hltDiPFJetAve140", "pt_ave_hltDiPFJetAve200", "pt_ave_hltDiPFJetAve260", "pt_ave_hltDiPFJetAve320","pt_ave_hltDiPFJetAve400"};

//     HLT_DiPFJetAve140 = triggerBits[1]
//     HLT_DiPFJetAve200 = triggerBits[3]
//     HLT_DiPFJetAve260 = triggerBits[5]
//     HLT_DiPFJetAve320 = triggerBits[7]
//     HLT_DiPFJetAve40  = triggerBits[9]
//     HLT_DiPFJetAve400 = triggerBits[8]
//     HLT_DiPFJetAve500 = triggerBits[10]
//     HLT_DiPFJetAve60  = triggerBits[12]
//     HLT_DiPFJetAve80  = triggerBits[13]

    TString FileNameRun2[9] = {"pt_ave_hltDiPFJetAve140","pt_ave_hltDiPFJetAve200", "pt_ave_hltDiPFJetAve260", "pt_ave_hltDiPFJetAve320","pt_ave_hltDiPFJetAve40",  "pt_ave_hltDiPFJetAve400","pt_ave_hltDiPFJetAve500","pt_ave_hltDiPFJetAve60","pt_ave_hltDiPFJetAve80"};
  //   if(eventInfo->triggerBits[1]==1) hist(FileNameRun2[0])->Fill(ev.pt_ave, weight);
//     if(eventInfo->triggerBits[3]==1) hist(FileNameRun2[1])->Fill(ev.pt_ave, weight);
//     if(eventInfo->triggerBits[5]==1) hist(FileNameRun2[2])->Fill(ev.pt_ave, weight);
//     if(eventInfo->triggerBits[7]==1) hist(FileNameRun2[3])->Fill(ev.pt_ave, weight);
//     if(eventInfo->triggerBits[9]==1) hist(FileNameRun2[4])->Fill(ev.pt_ave, weight);
//     if(eventInfo->triggerBits[8]==1) hist(FileNameRun2[5])->Fill(ev.pt_ave, weight);
//     if(eventInfo->triggerBits[10]==1) hist(FileNameRun2[6])->Fill(ev.pt_ave, weight);
//     if(eventInfo->triggerBits[12]==1) hist(FileNameRun2[7])->Fill(ev.pt_ave, weight);
//     if(eventInfo->triggerBits[14]==1) hist(FileNameRun2[8])->Fill(ev.pt_ave, weight);
// //     if(eventInfo->triggerBits[12]==1 || eventInfo->triggerBits[11]==1) hist(FileNameRun2[7])->Fill(ev.pt_ave, weight);
// //     if(eventInfo->triggerBits[14]==1 || eventInfo->triggerBits[13]==1) hist(FileNameRun2[8])->Fill(ev.pt_ave, weight);

    if(eventInfo.triggerBits[1]==1) hist(FileNameRun2[0])->Fill(ev.get(tt_pt_ave), weight);
    if(eventInfo.triggerBits[3]==1) hist(FileNameRun2[1])->Fill(ev.get(tt_pt_ave), weight);
    if(eventInfo.triggerBits[5]==1) hist(FileNameRun2[2])->Fill(ev.get(tt_pt_ave), weight);
    if(eventInfo.triggerBits[7]==1) hist(FileNameRun2[3])->Fill(ev.get(tt_pt_ave), weight);
    if(eventInfo.triggerBits[9]==1) hist(FileNameRun2[4])->Fill(ev.get(tt_pt_ave), weight);
    if(eventInfo.triggerBits[8]==1) hist(FileNameRun2[5])->Fill(ev.get(tt_pt_ave), weight);
    if(eventInfo.triggerBits[10]==1) hist(FileNameRun2[6])->Fill(ev.get(tt_pt_ave), weight);
    if(eventInfo.triggerBits[12]==1) hist(FileNameRun2[7])->Fill(ev.get(tt_pt_ave), weight);
    if(eventInfo.triggerBits[14]==1) hist(FileNameRun2[8])->Fill(ev.get(tt_pt_ave), weight);


    if (fabs(ev.get(tt_probejet_eta)) < 2.6) {
        TString FileNameRun2HF[6] = {"HLT_DiPFJetAve100_HFJEC","HLT_DiPFJetAve160_HFJEC", "HLT_DiPFJetAve220_HFJEC", "HLT_DiPFJetAve300_HFJEC","HLT_DiPFJetAve60_HFJEC",  "HLT_DiPFJetAve80_HFJEC"};
        // if(eventInfo->triggerBits[0]==1) hist(FileNameRun2HF[0])->Fill(ev.pt_ave, weight);
        // if(eventInfo->triggerBits[2]==1) hist(FileNameRun2HF[1])->Fill(ev.pt_ave, weight);
        // if(eventInfo->triggerBits[4]==1) hist(FileNameRun2HF[2])->Fill(ev.pt_ave, weight);
        // if(eventInfo->triggerBits[6]==1) hist(FileNameRun2HF[3])->Fill(ev.pt_ave, weight);
        // if(eventInfo->triggerBits[11]==1) hist(FileNameRun2HF[4])->Fill(ev.pt_ave, weight);
        // if(eventInfo->triggerBits[13]==1) hist(FileNameRun2HF[5])->Fill(ev.pt_ave, weight);
	if(eventInfo.triggerBits[0]==1) hist(FileNameRun2HF[0])->Fill(ev.get(tt_pt_ave), weight);
        if(eventInfo.triggerBits[2]==1) hist(FileNameRun2HF[1])->Fill(ev.get(tt_pt_ave), weight);
        if(eventInfo.triggerBits[4]==1) hist(FileNameRun2HF[2])->Fill(ev.get(tt_pt_ave), weight);
        if(eventInfo.triggerBits[6]==1) hist(FileNameRun2HF[3])->Fill(ev.get(tt_pt_ave), weight);
        if(eventInfo.triggerBits[11]==1) hist(FileNameRun2HF[4])->Fill(ev.get(tt_pt_ave), weight);
        if(eventInfo.triggerBits[13]==1) hist(FileNameRun2HF[5])->Fill(ev.get(tt_pt_ave), weight);
    }
    TString FileNameRun2probejet_eta[9] = {"probejet_eta_hltDiPFJetAve140","probejet_eta_hltDiPFJetAve200", "probejet_eta_hltDiPFJetAve260", "probejet_eta_hltDiPFJetAve320","probejet_eta_hltDiPFJetAve40",  "probejet_eta_hltDiPFJetAve400","probejet_eta_hltDiPFJetAve500","probejet_eta_hltDiPFJetAve60","probejet_eta_hltDiPFJetAve80"};
   //  if(eventInfo->triggerBits[1]==1) hist(FileNameRun2probejet_eta[0])->Fill(ev.probejet_eta, weight);
//     if(eventInfo->triggerBits[3]==1) hist(FileNameRun2probejet_eta[1])->Fill(ev.probejet_eta, weight);
//     if(eventInfo->triggerBits[5]==1) hist(FileNameRun2probejet_eta[2])->Fill(ev.probejet_eta, weight);
//     if(eventInfo->triggerBits[7]==1) hist(FileNameRun2probejet_eta[3])->Fill(ev.probejet_eta, weight);
//     if(eventInfo->triggerBits[9]==1) hist(FileNameRun2probejet_eta[4])->Fill(ev.probejet_eta, weight);
//     if(eventInfo->triggerBits[8]==1) hist(FileNameRun2probejet_eta[5])->Fill(ev.probejet_eta, weight);
//     if(eventInfo->triggerBits[10]==1) hist(FileNameRun2probejet_eta[6])->Fill(ev.probejet_eta, weight);
//     if(eventInfo->triggerBits[12]==1) hist(FileNameRun2probejet_eta[7])->Fill(ev.probejet_eta, weight);
//     if(eventInfo->triggerBits[14]==1) hist(FileNameRun2probejet_eta[8])->Fill(ev.probejet_eta, weight);
// //     if(eventInfo->triggerBits[12]==1 || eventInfo->triggerBits[11]==1) hist(FileNameRun2probejet_eta[7])->Fill(ev.probejet_eta, weight);
// //     if(eventInfo->triggerBits[14]==1 || eventInfo->triggerBits[13]==1) hist(FileNameRun2probejet_eta[8])->Fill(ev.probejet_eta, weight);
    if(eventInfo.triggerBits[1]==1) hist(FileNameRun2probejet_eta[0])->Fill(ev.get(tt_probejet_eta), weight);
    if(eventInfo.triggerBits[3]==1) hist(FileNameRun2probejet_eta[1])->Fill(ev.get(tt_probejet_eta), weight);
    if(eventInfo.triggerBits[5]==1) hist(FileNameRun2probejet_eta[2])->Fill(ev.get(tt_probejet_eta), weight);
    if(eventInfo.triggerBits[7]==1) hist(FileNameRun2probejet_eta[3])->Fill(ev.get(tt_probejet_eta), weight);
    if(eventInfo.triggerBits[9]==1) hist(FileNameRun2probejet_eta[4])->Fill(ev.get(tt_probejet_eta), weight);
    if(eventInfo.triggerBits[8]==1) hist(FileNameRun2probejet_eta[5])->Fill(ev.get(tt_probejet_eta), weight);
    if(eventInfo.triggerBits[10]==1) hist(FileNameRun2probejet_eta[6])->Fill(ev.get(tt_probejet_eta), weight);
    if(eventInfo.triggerBits[12]==1) hist(FileNameRun2probejet_eta[7])->Fill(ev.get(tt_probejet_eta), weight);
    if(eventInfo.triggerBits[14]==1) hist(FileNameRun2probejet_eta[8])->Fill(ev.get(tt_probejet_eta), weight);


    TString FileNameRun2HFprobejet_eta[6] = {"probejet_eta_HLT_DiPFJetAve100_HFJEC","probejet_eta_HLT_DiPFJetAve160_HFJEC", "probejet_eta_HLT_DiPFJetAve220_HFJEC", "probejet_eta_HLT_DiPFJetAve300_HFJEC","probejet_eta_HLT_DiPFJetAve60_HFJEC",  "probejet_eta_HLT_DiPFJetAve80_HFJEC"};
    // if(eventInfo->triggerBits[0]==1) hist(FileNameRun2HFprobejet_eta[0])->Fill(ev.probejet_eta, weight);
    // if(eventInfo->triggerBits[2]==1) hist(FileNameRun2HFprobejet_eta[1])->Fill(ev.probejet_eta, weight);
    // if(eventInfo->triggerBits[4]==1) hist(FileNameRun2HFprobejet_eta[2])->Fill(ev.probejet_eta, weight);
    // if(eventInfo->triggerBits[6]==1) hist(FileNameRun2HFprobejet_eta[3])->Fill(ev.probejet_eta, weight);
    // if(eventInfo->triggerBits[11]==1) hist(FileNameRun2HFprobejet_eta[4])->Fill(ev.probejet_eta, weight);
    // if(eventInfo->triggerBits[13]==1) hist(FileNameRun2HFprobejet_eta[5])->Fill(ev.probejet_eta, weight);

    if(eventInfo.triggerBits[0]==1) hist(FileNameRun2HFprobejet_eta[0])->Fill(ev.get(tt_probejet_eta), weight);
    if(eventInfo.triggerBits[2]==1) hist(FileNameRun2HFprobejet_eta[1])->Fill(ev.get(tt_probejet_eta), weight);
    if(eventInfo.triggerBits[4]==1) hist(FileNameRun2HFprobejet_eta[2])->Fill(ev.get(tt_probejet_eta), weight);
    if(eventInfo.triggerBits[6]==1) hist(FileNameRun2HFprobejet_eta[3])->Fill(ev.get(tt_probejet_eta), weight);
    if(eventInfo.triggerBits[11]==1) hist(FileNameRun2HFprobejet_eta[4])->Fill(ev.get(tt_probejet_eta), weight);
    if(eventInfo.triggerBits[13]==1) hist(FileNameRun2HFprobejet_eta[5])->Fill(ev.get(tt_probejet_eta), weight);

    TVector2 pt, met;
    if(fabs(jet1->eta) < s_eta_barr){
      pt.Set(ev.get(tt_jet1_pt) * cos(jet1->phi),ev.get(tt_jet1_pt) * sin(jet1->phi));
      hist("generic_asym")->Fill((ev.get(tt_jet2_pt) - ev.get(tt_jet1_pt)) / (ev.get(tt_jet2_pt) + ev.get(tt_jet1_pt)), weight);
        //j(E_{jet}) = 1 + \frac{ {E^{\gamma}_{T}} \cdot { \slashed{E}_{T} } } { (E^{\gamma}_{T})^{2} } 
        hist("generic_mpf")->Fill(1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py()), weight);
        hist("generic_r_rel")->Fill(ev.get(tt_jet2_pt) / ev.get(tt_jet1_pt), weight);
        ((TH2D*)hist("Rrel_vs_Npv"))->Fill(nvertices,ev.get(tt_jet2_pt) / ev.get(tt_jet1_pt),weight);
//         hist("prof_Rrel_vs_Npv")->Fill(nvertices,ev.jet2_pt / ev.jet1_pt);

    }
    if(fabs(jet2->eta) < s_eta_barr){
      pt.Set(ev.get(tt_jet2_pt) * cos(jet2->phi),ev.get(tt_jet2_pt) * sin(jet2->phi));
      hist("generic_asym")->Fill((ev.get(tt_jet1_pt) - ev.get(tt_jet2_pt)) / (ev.get(tt_jet2_pt) + ev.get(tt_jet1_pt)), weight);
        //j(E_{jet}) = 1 + \frac{ {E^{\gamma}_{T}} \cdot { \slashed{E}_{T} } } { (E^{\gamma}_{T})^{2} }
        hist("generic_mpf")->Fill(1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py()), weight);
        hist("generic_r_rel")->Fill(ev.get(tt_jet1_pt) / ev.get(tt_jet2_pt), weight);
        ((TH2D*)hist("Rrel_vs_Npv"))->Fill(nvertices,ev.get(tt_jet1_pt) / ev.get(tt_jet2_pt),weight);
//       hist("prof_Rrel_vs_Npv")->Fill(nvertices,ev.jet1_pt / ev.jet2_pt);

    }

    ((TH2D*)hist("mpf_vs_etaProbe"))->Fill(ev.get(tt_probejet_eta),ev.get(tt_mpf_r),weight);
    ((TH2D*)hist("r_rel_vs_etaProbe"))->Fill(ev.get(tt_probejet_eta),ev.get(tt_rel_r),weight);
    float pt_ave = (0.5*(ev.get(tt_jet1_pt) + ev.get(tt_jet2_pt)));
    ((TH2D*)hist("pt_ave_vs_etaProbe"))->Fill(ev.get(tt_probejet_eta),pt_ave,weight);
    baconhep::TJet* jet3 = (baconhep::TJet*)js[2];
    if (njets > 2){
      hist("pt_3")->Fill(ev.get(tt_jet3_pt), weight);
        hist("eta_3")->Fill(jet3->eta, weight);
        hist("pt_rel")->Fill(ev.get(tt_jet3_pt)/(0.5*(ev.get(tt_barreljet_pt) + ev.get(tt_probejet_pt) )),weight);
	float pt_jet3 = ev.get(tt_jet3_pt);
	float alpha = ev.get(tt_jet3_pt)/(0.5*(ev.get(tt_jet1_pt) + ev.get(tt_jet2_pt)));
	hist("generic_pt_rel")->Fill(alpha,weight);
        hist("ptrel_vs_deltaphi")->Fill(alpha,deltaPhi);
	((TH2D*)hist("ptjet3_vs_alpha"))->Fill(alpha,pt_jet3,weight);
	((TH2D*)hist("pt_ave_vs_alpha"))->Fill(alpha,pt_ave,weight);
	((TH2D*)hist("alpha_vs_assym"))->Fill(ev.get(tt_asymmetry),alpha, weight);
	((TH2D*)hist("Njets_vs_alpha"))->Fill(alpha,njets,weight);
	((TH2D*)hist("alphaSum_vs_alpha3"))->Fill(alpha,ev.get(tt_alpha_sum),weight);
	// ((TH2D*)hist("diffPt1_vs_alpha3"))->Fill(alpha,((jet1->pt)-(jet1->genpt)),weight);
	// ((TH2D*)hist("diffPt2_vs_alpha3"))->Fill(alpha,((jet2->pt)-(jet2->genpt)),weight);
	// ((TH2D*)hist("diffPt1_vs_alphaSum"))->Fill(ev.get(tt_alpha_sum),((jet1->pt)-(jet1->genpt)),weight);
	// ((TH2D*)hist("diffPt2_vs_alphaSum"))->Fill(ev.get(tt_alpha_sum),((jet2->pt)-(jet2->genpt)),weight);
	// ((TH2D*)hist("diffPtGen_vs_alpha3"))->Fill(alpha,abs(2*((jet1->genpt)-(jet2->genpt))/((jet1->genpt)+(jet2->genpt))),weight);
	// ((TH2D*)hist("diffPtGen_vs_alphaSum"))->Fill(ev.get(tt_alpha_sum),abs(2*((jet1->genpt)-(jet2->genpt))/((jet1->genpt)+(jet2->genpt))),weight);
	double diffPt = (jet1->genpt)-(jet2->genpt);
	double avePt = 0.5*((jet1->genpt)+(jet2->genpt));
	if(jet2->genpt>0 && jet1->genpt>0){
	  //	if((diffPt/avePt)>1.99) 
	// cout<<"diffPt/avePt = "<<diffPt/avePt<<" diffPt = "<<diffPt<<" avePt = "<<avePt<<" jet1->genpt = "
	//     <<jet1->genpt<<" jet2->genpt = "<<jet2->genpt<<" jet2->pt = "<<jet2->pt<<endl;
	((TH2D*)hist("diffPtGen_vs_alpha3"))->Fill(alpha,(diffPt/avePt),weight);
	((TH2D*)hist("diffPtGen_vs_alphaSum"))->Fill(ev.get(tt_alpha_sum),(diffPt/avePt),weight);
	// ((TH2D*)hist("AbsdiffPtGen_vs_alpha3"))->Fill(alpha,fabs(diffPt/avePt),weight);
	// ((TH2D*)hist("AbsdiffPtGen_vs_alphaSum"))->Fill(ev.get(tt_alpha_sum),fabs(diffPt/avePt),weight);
	}
    }
}
JECAnalysisHists::~JECAnalysisHists(){}
