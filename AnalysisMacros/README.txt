For 2016 DATA
1) hadd QCD Flat and QCD fwd to Full, e.g. uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_AK4CHS_Full_RunEFearly.root
2) hadd Data according to different periods
3) Run main with "CalculateMCWeights()" 
4) Run main with "FullCycle_CorrectFormulae()"

How to compile the program:
$> g++ -std=c++0x -Wall -o main main.C CorrectionObject.cc ControlPlots.cc kFSR_CorrectFormulae.cc Pt_Extrapolation_Alternative_CorrectFormulae.cc L2ResOutput.cc L2ResAllRuns.cc InputForGlobalFit.cc FinalControlPlots_CorrectFormulae.cc useful_functions.cc CalculateMCWeights.cc tdrstyle_mod15.C `root-config  --cflags --evelibs`

$> g++ -std=c++0x -Wall -o main main.C CorrectionObject.cc ControlPlots.cc kFSR_CorrectFormulae.cc Pt_Extrapolation_Alternative_CorrectFormulae.cc L2ResOutput.cc L2ResAllRuns.cc InputForGlobalFit.cc FinalControlPlots_CorrectFormulae.cc useful_functions.cc CalculateMCWeights_TriggerThresholds.cc tdrstyle_mod15.C `root-config  --cflags --evelibs`
How to execute the program:
./main
