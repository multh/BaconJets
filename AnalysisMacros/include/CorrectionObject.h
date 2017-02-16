#pragma once

#include <cmath>
#include <iostream>
#include <TString.h>
#include <TFile.h>




using namespace std;

  class CorrectionObject {



  public:

    // Constructors, destructor
    CorrectionObject(const TString & runnr, const TString & generator, const TString & collection, const bool & split_JEC = true, const bool & closuretest = false);
    CorrectionObject(const CorrectionObject &) = default;
    CorrectionObject & operator = (const CorrectionObject &) = default;
    ~CorrectionObject() = default;
    inline void CloseFiles(){_MCFile->Close(); _DATAFile->Close();};
    
    // Setter and getter functions
    inline TString runnr(){ return _runnr;}
    inline TString collection(){ return _collection;}
    inline TString generator(){ return _generator;}
    inline TString jettag(){ return _jettag;}
    inline bool closuretest(){return _closuretest;}
    inline TString MCPath(){return _MCpath;}
    inline TString DATAPath(){return _DATApath;}
    inline TString OutPath(){return _outpath;}
    inline TString lumitag(){return _lumitag;}
    inline const TString runnr() const{return _runnr;}
    inline const TString collection() const{return _collection;}
    inline const TString generator() const{return _generator;}
    inline const TString jettag() const{return _jettag;}
    inline const bool closuretest() const {return _closuretest;}
    inline const TString MCPath() const {return _MCpath;}
    inline const TString DATAPath() const {return _DATApath;}
    inline const TString OutPath() const {return _outpath;}
    inline const TString lumitag() const {return _lumitag;}
    inline void set_runnr(TString x){_runnr = x;}
    inline void set_collection(TString x){_collection = x;}
    inline void set_generator(TString x){_generator = x;}
    inline void set_jettag(TString x){_jettag = x;}
    inline void set_closuretest(bool x){_closuretest = x;}
    inline void set_MCPath(TString x){_MCpath = x; _MCFile->Close(); _MCFile = new TFile(_MCpath,"READ");}
    inline void set_DATAPath(TString x){_DATApath = x; _DATAFile->Close(); _DATAFile = new TFile(_DATApath,"READ");}
    inline void set_outpath(TString x){_outpath = x;}
    inline void set_lumitag(TString x){_lumitag = x;}

    //Main functions for calculating L2 residuals, defined in CorrectionObject.cc
    void ControlPlots();
    void kFSR();
    void kFSR_CorrectFormulae();
    void Pt_Extrapolation(bool mpfMethod = true);
    void Pt_Extrapolation_Alternative(bool mpfMethod = true);
    void Pt_Extrapolation_Alternative_CorrectFormulae(bool mpfMethod = true);
    void L2ResOutput();
    void L2ResAllRuns();
    void L2ResOverlay();
    void InputForGlobalFit();
    void InputForGlobalFit_eta_0_13();
    void FinalControlPlots();
    void FinalControlPlots_CorrectFormulae();
    void CalculateMCWeights();
    void FullCycle_CorrectFormulae();

  private:
    TString _runnr;
    TString _collection;
    TString _generator, _generator_tag;
    TString _jettag;
    TString _lumitag;
    TString _MCpath, _MCpath_ForWeights, _DATApath, _DATApath_ForWeights;
    TString _outpath;
    TString _weightpath;
    TFile*  _MCFile;
    TFile*  _DATAFile;
    bool    _split_JEC;
    bool    _closuretest;
 
   

  }; // end of class CorrectionObject



// Equality operators

inline bool operator == (const CorrectionObject & a, const CorrectionObject & b) {
  return a.runnr() == b.runnr() && a.collection() == b.collection() && a.generator() == b.generator() && a.closuretest() == b.closuretest();
}


inline bool operator != (const CorrectionObject & a, const CorrectionObject & b) {
  return a.runnr() != b.runnr() || a.collection() != b.collection() || a.generator() != b.generator();
}

// I/O operations

inline ostream & operator << (ostream & os, const CorrectionObject & q) {
  return os << "(" << q.runnr() << "," << endl 
                  << q.generator() << ","  << endl
                  << q.collection() << "," << endl 
                  << q.closuretest() << "," << endl
                  << q.MCPath() << ","  << endl 
	          << q.DATAPath() << "," << endl
	          << q.OutPath() << ")";
}


