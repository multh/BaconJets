#include <TString.h>
#include <TFile.h>

#include "../include/CorrectionObject.h"


CorrectionObject::CorrectionObject(const TString & runnr, const TString & generator, const TString & collection, const bool & split_JEC, const bool & closuretest) :
  _runnr(runnr), _collection(collection), _generator(generator), _split_JEC(split_JEC), _closuretest(closuretest)
    {
      TString s_tmp = _collection;
      s_tmp.ReplaceAll("CHS", "PFchs");
      s_tmp.ReplaceAll("Puppi", "PFpuppi");
      _jettag = s_tmp;

      TString input_path;
      if(!_closuretest){  
	if(split_JEC){
	  input_path = "/nfs/dust/cms/user/reimersa/JEC/2016ReReco/Residuals/Summer16_23Sep2016_V3/" + _collection + "/MC_Reweighted_PtEta_Pt95/";
	  _weightpath =  "/nfs/dust/cms/user/reimersa/JEC/2016ReReco/Residuals/Summer16_23Sep2016_V3/" + _collection + "/";
	  _outpath =   input_path + "Run" + _runnr + "/";
	  //input_path =   "/nfs/dust/cms/user/reimersa/JEC/2016ReReco/Residuals/Spring16_23Sep2016_V1/" + _collection + "/";
	  //_outpath =     "/nfs/dust/cms/user/reimersa/JEC/2016ReReco/Residuals/Spring16_23Sep2016_V1/" + _collection + "/Run" + _runnr + "/";
	}
	else{
	  input_path = "/nfs/dust/cms/user/reimersa/JEC/2016ReReco/Residuals/Spring16_v8p2/GlobalJEC/" + _collection + "/";
	  _weightpath =  "/nfs/dust/cms/user/reimersa/JEC/2016ReReco/Residuals/Spring16_v8p2/GlobalJEC/";
	  _outpath = "/nfs/dust/cms/user/reimersa/JEC/2016ReReco/Residuals/Spring16_v8p2/GlobalJEC/" + _collection + "/Run" + _runnr + "/";
	}
      }
      else{
	input_path = "/nfs/dust/cms/user/reimersa/JEC/2016ReReco/ClosureTest/Summer16_23Sep2016_V3/" + _collection + "/";
	_weightpath =  "/nfs/dust/cms/user/reimersa/JEC/2016ReReco/ClosureTest/Summer16_23Sep2016_V3/";
	_outpath = "/nfs/dust/cms/user/reimersa/JEC/2016ReReco/ClosureTest/Summer16_23Sep2016_V3/" + _collection + "/Run" + _runnr + "/";
      }

      if(_generator == "pythia"){
	_MCpath = input_path + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_" + _collection  + "_Full_Run" + _runnr  + ".root";
	_MCpath_ForWeights = _weightpath + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_" + _collection  + "_Full.root";
	//_MCpath = input_path + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_" + _collection  + ".root";
	_generator_tag = "pythia8";
      }
      else if(_generator == "herwig"){
	_MCpath = input_path + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_herwigpp_"+ _collection  +".root";
	_MCpath_ForWeights = _weightpath + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_herwigpp_" + _collection  + ".root";
	_generator_tag = "herwigpp";
      }
      else if(_generator == "madgraph"){
	_MCpath = input_path + "uhh2.AnalysisModuleRunner.MC.QCDHtFULL_madgraph_"+ _collection  +".root";
	_MCpath_ForWeights = _weightpath + "uhh2.AnalysisModuleRunner.MC.QCDHtFULL_madgraph_" + _collection  + ".root";
	_generator_tag = "madgraphMLM";
      }
      _DATApath = input_path + "uhh2.AnalysisModuleRunner.DATA.DATA_Run" + _runnr + "_" + _collection + ".root";
      _DATApath_ForWeights = _weightpath + "uhh2.AnalysisModuleRunner.DATA.DATA_Run" + _runnr + "_" + _collection + ".root";
      //_DATApath = "/nfs/dust/cms/user/reimersa/JEC/2016ReReco/Residuals/Spring16_23Sep2016_V1/AK4CHS/uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_AK4CHS.root";
      cout << "Opening MC file:   " << _MCpath << endl;
      cout << "Opening DATA file: " << _DATApath << endl << endl;
      _MCFile = new TFile(_MCpath,"READ");
      _DATAFile = new TFile(_DATApath,"READ");

      //lumitags
      if(_runnr == "BCD")    _lumitag      = "RunBCD  12.9 fb^{-1}";
      else if(_runnr == "G") _lumitag      = "RunG  7.6 fb^{-1}";
      else if(_runnr == "B") _lumitag      = "RunB  5.8 fb^{-1}";
      else if(_runnr == "C") _lumitag      = "RunC  2.6 fb^{-1}";
      else if(_runnr == "D") _lumitag      = "RunD  4.3 fb^{-1}";
      else if(_runnr == "E") _lumitag      = "RunE  4.1 fb^{-1}";
      else if(_runnr == "F") _lumitag      = "RunF  3.2 fb^{-1}";
      else if(_runnr == "G") _lumitag      = "RunG  7.6 fb^{-1}";
      else if(_runnr == "Fearly") _lumitag = "RunFearly  2.7 fb^{-1}";
      else if(_runnr == "EFearly") _lumitag = "RunEFearly  6.8 fb^{-1}";
      else if(_runnr == "Flate") _lumitag = "RunFlate  0.4 fb^{-1}";
      else if(_runnr == "FlateG") _lumitag = "RunFlateG  8.0 fb^{-1}";
      else if(_runnr == "FlateGH") _lumitag = "RunFlateGH  X.X fb^{-1}";
      else if(_runnr == "H") _lumitag = "RunH  X.X fb^{-1}";
      else if(_runnr == "BCDEFearly") _lumitag = "RunBCDEFearly  19.7 fb^{-1}";
      else if(_runnr == "BCDEFGH") _lumitag = "RunBCDEFGH  36.8 fb^{-1}";
      else throw runtime_error("In constructor: Invalid RunNr. specified.");
    }

void CorrectionObject::FullCycle_CorrectFormulae(){
  CorrectionObject::ControlPlots();
  CorrectionObject::kFSR_CorrectFormulae();
  CorrectionObject::Pt_Extrapolation_Alternative_CorrectFormulae(true);
  CorrectionObject::Pt_Extrapolation_Alternative_CorrectFormulae(false);
  CorrectionObject::L2ResOutput();
  CorrectionObject::FinalControlPlots_CorrectFormulae();
}
