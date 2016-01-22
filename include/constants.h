#ifndef  CONSTANTS_H
#define  CONSTANTS_H

/** \brief Binning **/
// static std::vector<double>   eta_range  = {0, 0.261, 0.522, 0.763, 0.957, 1.131, 1.305, 1.479, 1.93, 2.322, 2.411, 2.5, 2.853, 2.964, 3.139, 3.489, 5.191};
// static std::vector<double>   pt_range   = {66, 107, 191, 240, 306, 379, 468, 900};
// static std::vector<double>   alpha_range= {0., 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25};


// static std::vector<double>   eta_range  = {0, 0.25, 0.5, 0.75, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.7, 3.0, 3.2, 5.0};//{0, 1.3, 2.5, 3.0, 5.0};
//static std::vector<double>   eta_range  = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.191};


// additional HF bin:
static std::vector<double>   eta_range  = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 5.191};
// static std::vector<double>   eta_range_mikko_barrel = {0, 0.8, 1.3};


static std::vector<double>   eta_range_mikko  = {0, 1.3, 1.9, 2.5, 3.0, 3.2, 5.0, 0, 0.8, 1.3};

//static std::vector<double>   pt_range   = {55, 76, 93, 172, 232, 300, 366, 452, 558};
static std::vector<double>   pt_range   = {56, 78, 100, 168, 232, 300, 366, 453, 562}; //2.11 nominal triggers



// static std::vector<double>   pt_range   = {43, 80, 88, 135, 223, 290, 365, 448, 561};
// static std::vector<double>   pt_range   = {57, 80, 100, 166, 231, 290, 365, 478, 556};

static std::vector<double>   alpha_range= {0., 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25};

/** \brief Dijet event selection **/
// barrel region (|eta| < 1.3)
static float s_eta_barr = 1.3;
// two back-to-back leading jets (delta_phi(j1,j2) = min(|phi1 - phi2|, 2PI - |phi2 - phi1|) > 2.7)
static float s_delta_phi = 2.7;
// cut on the asymmetry for events with two jets  |(j2->pt - j1->pt /(j2->pt + j1->pt)| < 0.70
static float s_asymm = 0.7;
// relative third jet fraction pt_rel = 2*j3_pt/(j1_pt + j2_pt) < 0.2
static float s_pt_rel = 0.4;

/** \brief good Primary Vertex reconstruction **/
// more than four tracks
static float s_n_PvTracks = 4;
// PV is located within 24cm in z vertex
static float s_n_Pv_z = 24.0;
// PV is located within 2cm in xy direction from the nominal interaction point
 static float s_n_Pv_xy = 2.0;


/** \brief The trigger thresholds of pt_ave **/
/// Used for the 2012 data
// static float s_Pt_Ave40_cut   = 66; //Hennings: 62
// static float s_Pt_Ave80_cut   = 107;//Hennings: 107
// static float s_Pt_Ave140_cut  = 191;//Hennings: 175
// static float s_Pt_Ave200_cut  = 240;//Hennings: 242
// static float s_Pt_Ave260_cut  = 306;//Hennings: 310
// static float s_Pt_Ave320_cut  = 379;//Hennings: 379
// static float s_Pt_Ave400_cut  = 468;//Hennings: 467
/// Used for the 2015 data
// static float s_Pt_Ave40_cut   = 43;
// static float s_Pt_Ave60_cut   = 80;
// static float s_Pt_Ave80_cut   = 88;
// static float s_Pt_Ave140_cut  = 135;
// static float s_Pt_Ave200_cut  = 223;
// static float s_Pt_Ave260_cut  = 290;
// static float s_Pt_Ave320_cut  = 365;
// static float s_Pt_Ave400_cut  = 448;
// static float s_Pt_Ave500_cut  = 561;

/// Used for the 2015  2.11/fb  25ns data
static float s_Pt_Ave40_cut   = 56;
static float s_Pt_Ave60_cut   = 78;
static float s_Pt_Ave80_cut   = 100;
static float s_Pt_Ave140_cut  = 168;
static float s_Pt_Ave200_cut  = 232;
static float s_Pt_Ave260_cut  = 300;
static float s_Pt_Ave320_cut  = 366;
static float s_Pt_Ave400_cut  = 453;
static float s_Pt_Ave500_cut  = 562;
static int n_pt_bins = 9;

/// HF thresholds Used for the 2015  2.11/fb  25ns data
static float s_Pt_Ave60HF_cut   = 77;
static float s_Pt_Ave80HF_cut   = 131;
static float s_Pt_Ave100HF_cut  = 154;
static float s_Pt_Ave160HF_cut  = 244;
static float s_Pt_Ave220HF_cut  = 321;
static float s_Pt_Ave300HF_cut  = 426;




/// Used for the 2015  40/pb  50ns data
/*
static float s_Pt_Ave40_cut   = 55;
static float s_Pt_Ave60_cut   = 76;
static float s_Pt_Ave80_cut   = 93;
static float s_Pt_Ave140_cut  = 172;
static float s_Pt_Ave200_cut  = 232;
static float s_Pt_Ave260_cut  = 300;
static float s_Pt_Ave320_cut  = 366;
static float s_Pt_Ave400_cut  = 453;
static float s_Pt_Ave500_cut  = 558;
static int n_pt_bins = 9;
*/
/// Used for the 2015  ~1/fb  25ns data
// static float s_Pt_Ave40_cut   = 55;
// static float s_Pt_Ave60_cut   = 77;
// static float s_Pt_Ave80_cut   = 99;
// static float s_Pt_Ave140_cut  = 165;
// static float s_Pt_Ave200_cut  = 231;
// static float s_Pt_Ave260_cut  = 298;
// static float s_Pt_Ave320_cut  = 365;
// static float s_Pt_Ave400_cut  = 451;
// static float s_Pt_Ave500_cut  = 561;

/// HF thresholds Used for the 2015  40/pb  50ns data
/*
static float s_Pt_Ave60HF_cut   = 76;
static float s_Pt_Ave80HF_cut   = 93;
static float s_Pt_Ave100HF_cut  = 120;
static float s_Pt_Ave160HF_cut  = 186;
static float s_Pt_Ave220HF_cut  = 253;
static float s_Pt_Ave300HF_cut  = 342;
*/
/// HF thresholds Used for the 2015  ~1/fb  25ns data
// static float s_Pt_Ave60HF_cut   = 77;
// static float s_Pt_Ave80HF_cut   = 99;
// static float s_Pt_Ave100HF_cut  = 121;
// static float s_Pt_Ave160HF_cut  = 187;
// static float s_Pt_Ave220HF_cut  = 253;
// static float s_Pt_Ave300HF_cut  = 342;

/** \brief Jet Resolution Smearering **/
// doing the matching from GEN to RECO
static float s_delta_R   = 0.3; 
//constant numberstaken from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
// from 8TeV JER measurement.
constexpr const size_t n = 7;
//static float eta_hi[n]    = {0.5, 1.1, 1.7, 2.3, 2.8, 3.2, 5.0};
static float eta_hi[n]    = {0.8, 1.3, 1.9, 2.5, 3.0, 3.2, 5.0}; ///RunII
//static float c_nominal[n] = {1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056};
// static float c_up[n]      = {1.105, 1.127, 1.150, 1.254, 1.316, 1.458, 1.247};
// static float c_down[n]    = {1.053, 1.071, 1.092, 1.162, 1.192, 1.332, 0.865};

static float c_nominal[n] = {1.061, 1.088, 1.106, 1.126, 1.343, 1.303, 1.320};///RunII official for 13 TeV, 25ns collisions 23/11/15

static float c_err_nominal[n] = {0.023, 0.029, 0.030, 0.094, 0.123, 0.111, 0.286};///RunII official for 13 TeV, 25ns collisions 23/11/15




/** \brief Dijet event weighting **/
// in each pt-region defined by the trigger thresholds # of events in MC and in data should equal

///Asympt Mc
static float scale_factor_noJERA[] = {3153, 11221.4, 10024.6, 9616.32, 465256, 667300, 4.138e+06, 1.43267e+07, 1.69729e+07};//run2 69000 V6 Asym
static float scale_factor_centrA[] = {3001.18, 10864.1, 9891.95, 9607.91, 464409, 667504, 4.1411e+06, 1.43601e+07, 1.69549e+07};//run2 69000 V6, cent.smear Asym
static float scale_factor_upA[] = {2896.24, 10658.4, 9750.38, 9567.18, 462895, 665397, 4.1283e+06, 1.43305e+07, 1.6918e+07};//run2 69000 V6, up.smear Asym
static float scale_factor_downA[] = {3097.2, 11091.8, 9998.62, 9635.96, 466336, 669405, 4.15346e+06, 1.43912e+07, 1.6992e+07};//run2 69000 V6, down.smear Asym
static float scale_factor_centrA_PU80[] = {3415.88, 11734.9, 10368.5, 10024.6, 483457, 693496, 4.28772e+06, 1.49127e+07, 1.76746e+07};//run2 80000 V6, cent.smear Asym

static float scale_factor_centrA_PU69_set1[] = {1119.55, 3035.73, 2778.49, 2689.78, 69284.6, 99612, 612942, 1.89493e+06, 1.83119e+06};//69mb V6, cent. Asym set1
static float scale_factor_centrA_PU69_set2[] = {424.986, 1632.72, 1493.05, 1433.51, 63826.9, 91486.2, 542606, 1.97444e+06, 1.9212e+06};//69mb V6, cent. Asym set2
static float scale_factor_centrA_PU69_set3[] = {383.799, 1531.58, 1362.94, 1349.3, 58697.6, 87290.4, 519621, 2.02935e+06, 1.9623e+06};//69mb V6, cent. Asym set3
static float scale_factor_centrA_PU69_set4[] = {302.178, 1239.91, 1109.87, 1082.25, 52135.2, 74471.4, 448186, 1.77749e+06, 1.72852e+06};//69mb V6, cent. Asym set4
static float scale_factor_centrA_PU69_set5[] = {190.486, 817.24, 742.009, 760.26, 45077.2, 65552.7, 405774, 1.50184e+06, 1.74933e+06};//69mb V6, cent. Asym set5
static float scale_factor_centrA_PU69_set6[] = {113.783, 537.28, 497.731, 477.569, 40197.2, 57013.3, 377918, 1.14542e+06, 1.84246e+06};//69mb V6, cent. Asym set6
static float scale_factor_centrA_PU69_set7[] = {260.369, 1026.12, 971.377, 927.87, 49901.1, 69802.7, 435178, 1.3947e+06, 1.6765e+06};//69mb V6, cent. Asym set7
static float scale_factor_centrA_PU69_set8[] = {127.443, 601.259, 542.876, 510.048, 46458.6, 66212.6, 435145, 1.46857e+06, 2.16278e+06};//69mb V6, cent. Asym set8
static float scale_factor_centrA_PU69_set9[] = {78.5873, 442.229, 393.594, 377.323, 38830.2, 56062.8, 363732, 1.17338e+06, 2.08058e+06};//69mb V6, cent. Asym set9


///Flat MC
static float scale_factor_noJERF[] = {2543.11, 9228.77, 8404.76, 7916.86, 385342, 550228, 3.42556e+06, 1.17561e+07, 1.41737e+07};//run2 69000 V6 Flat
static float scale_factor_centrF[] = {2413.77, 8916.56, 8200.09, 7844.02, 382407, 544738, 3.39795e+06, 1.17032e+07, 1.40901e+07};//run2 69000 V6, cent.smear Flat



//static float scale_factor1[] = {2.14563, 4.27422,7.57198,  10.7035,17.899, 15.6095, 4.64629};//run1
// static float scale_factor2[] = {242.617, 3514.77, 49549.8, 180540, 693177, 1.83449e+06, 1.63241e+07};//run1




#endif
