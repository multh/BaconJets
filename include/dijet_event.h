#pragma once
//#include "UHH2/core/include/Jet.h"
using namespace std;
using namespace uhh2;
#include <vector>
#include <UHH2/core/include/Jet.h>


/** \brief Dijet events to store in the tree for JEC L2 residuals */
class dijet_event{
 public:
  dijet_event(){
    gen_pthat=0;gen_weight=0;
    jet1_pt =0; jet2_pt = 0; jet3_pt = 0; 
    jet1_ptRaw = 0; jet2_ptRaw = 0; jet3_ptRaw = 0;
    nvertices = 0;
    probejet_eta = 0; probejet_phi = 0; probejet_pt = 0; probejet_ptRaw = 0;
    barreljet_eta = 0; barreljet_phi = 0; barreljet_pt = 0;  barreljet_ptRaw = 0;
    pt_ave = 0; alpha = 0; 
    rel_r =0; mpf_r = 0;
    asymmetry = 0; nPU = 0;
    ev_weight = 0; jets_pt = 0; nJets = 0;
  }
  void set_gen_pthat(const float x){gen_pthat = x;}
  void set_gen_weight(const float x){gen_weight = x;}
  void set_jet123_pt(const float x1, const float x2, const float x3){
    jet1_pt = x1; jet2_pt = x2; jet3_pt = x3;
  }
  void set_jet123_ptRaw(const float x1, const float x2, const float x3){
    jet1_ptRaw = x1; jet2_ptRaw = x2; jet3_ptRaw = x3;
  }
  void set_nvertices(const float x){nvertices = x;}
  void set_probejet_kinematics(const float eta, const float phi, const float pt, const float ptRaw){
    probejet_eta = eta; probejet_phi = phi;  probejet_pt = pt;  probejet_ptRaw = ptRaw;
  }
  void set_barreljet_kinematics(const float eta, const float phi, const float pt, const float ptRaw){
    barreljet_eta = eta;  barreljet_phi = phi; barreljet_pt = pt; barreljet_ptRaw = ptRaw;
  }
  void set_pt_ave(const float x){pt_ave = x;}
  void set_alpha(const float x){alpha = x;}
  void set_rel_r(const float x){rel_r = x;}
  void set_mpf_r(const float x){mpf_r = x;}
  void set_asymmetry(const float x){asymmetry = x;}
  void set_nPU(const float x){nPU = x;}
  void set_ev_weight(const float x){ev_weight = x;}
  void set_jets_pt(const float x){jets_pt = x;}
  void set_nJets(const int x){nJets = x;}
 private:
  float gen_pthat; //pt hat (from QCD simulation)
  float gen_weight;// weight from MC
  float jet1_pt, jet2_pt, jet3_pt; //leading, subleading and the 3rd jet pt (corrected)
  float jet1_ptRaw, jet2_ptRaw, jet3_ptRaw;//leading, subleading and the 3rd jet pt (not corrected)
  float nvertices;//number of vertices
  float probejet_eta, probejet_phi, probejet_pt, probejet_ptRaw;// probe jet parameters
  float barreljet_eta, barreljet_phi, barreljet_pt, barreljet_ptRaw;//reference jet parameters 
  float pt_ave;//pt average of leading and subleading jets
  float alpha;// pt of the 3rd jet divided by  pt_ave
  float rel_r, mpf_r; //responces from pT-balance and MPF method
  float asymmetry;//asymmetry=(p_{T}^{probe}-p_{T}^{barrel})/(p_{T}^{probe}+p_{T}^{barrel})
  float nPU;//number of pile-up vertices
  float ev_weight;//weight of the event
  float jets_pt;//sum of additional jets pT (does _not_ include leading and subleading jets)
  int nJets;//number of jets
};
