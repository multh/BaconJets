#include <TString.h>
#include <TFile.h>

#include "../include/CorrectionObject.h"



CorrectionObject::CorrectionObject(const TString & runnr, const TString & generator, const TString & collection, const TString & input_path, const TString & weight_path, const bool & closuretest) :
  _runnr(runnr), _collection(collection), _generator(generator),_input_path(input_path),_weight_path(weight_path), _closuretest(closuretest)
    {
      TString s_tmp = _collection;
      s_tmp.ReplaceAll("CHS", "PFchs");
      s_tmp.ReplaceAll("Puppi", "PFpuppi");
      _jettag = s_tmp;

      TString inputPath;
      if(!_closuretest){ 
	inputPath       = _input_path;
	_weightpath_FLAT = _weight_path+ "/CENTRAL/";
	_weightpath_FWD  = _weight_path+ "/FWD/";
	_outpath =   inputPath + "Run" + _runnr + "/";
      }
      else{
	inputPath  = _input_path;
	_weightpath_FLAT = _weight_path+"/CENTRAL/";
	_weightpath_FWD  = _weight_path+"/FWD/";
	_outpath    = inputPath + "Run" + _runnr + "/";
	cout<< "Enter Closure Test Process: "<<endl<<endl; 

      }

      if(_generator == "pythia"){
        _MCpath = inputPath + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_AK4CHS_Flat.root";

       	//_MCpath = inputPath + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_" + _collection  + "_Full_Run" + _runnr  + ".root";
	_MCpath_ForWeights_FLAT = _weightpath_FLAT + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_" + _collection  + "_Flat.root";
	_MCpath_ForWeights_FWD  = _weightpath_FWD + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_" + _collection  + "_Fwd.root";
	_generator_tag = "pythia8";
      }
      else if(_generator == "herwig"){
	_MCpath = inputPath + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_herwigpp_"+ _collection  +".root";
	_MCpath_ForWeights_FLAT = _weightpath_FLAT + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_herwigpp_" + _collection  + "_Flat.root";
	_MCpath_ForWeights_FWD = _weightpath_FWD + "uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_herwigpp_" + _collection  + "_Fwd.root";
	_generator_tag = "herwigpp";
      }
      else if(_generator == "madgraph"){
	_MCpath = inputPath + "uhh2.AnalysisModuleRunner.MC.QCDHtFULL_madgraph_"+ _collection  +".root";
	_MCpath_ForWeights_FLAT = _weightpath_FLAT + "uhh2.AnalysisModuleRunner.MC.QCDHtFULL_madgraph_" + _collection  + "_Flat.root";
	_MCpath_ForWeights_FWD = _weightpath_FWD + "uhh2.AnalysisModuleRunner.MC.QCDHtFULL_madgraph_" + _collection  + "_Fwd.root";
	_generator_tag = "madgraphMLM";
      }

      _DATApath = inputPath + "uhh2.AnalysisModuleRunner.DATA.DATA_Run" + _runnr + "_" + _collection + ".root";
      _DATApath_ForWeights_FLAT = _weightpath_FLAT + "uhh2.AnalysisModuleRunner.DATA.DATA_Run" + _runnr + "_" + _collection + ".root";
      _DATApath_ForWeights_FWD = _weightpath_FWD + "uhh2.AnalysisModuleRunner.DATA.DATA_Run" + _runnr + "_" + _collection + ".root";

      
      cout << "Opening MC file:   " << _MCpath << endl;
      cout << "Opening DATA file: " << _DATApath << endl << endl;


      _MCFile = new TFile(_MCpath,"READ");
      _DATAFile = new TFile(_DATApath,"READ");

      if(_MCFile->GetSize()==-1) throw runtime_error("In CorrectionObject.cc: File or Directory " + _MCpath+" does not exist!");
      if(_DATAFile->GetSize()==-1) throw runtime_error("In CorrectionObject.cc: File or Directory " + _DATApath+" does not exist!");

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
