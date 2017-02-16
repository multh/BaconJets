How to compile the program:
$> g++ -std=c++0x -Wall -o main main.C CorrectionObject.cc ControlPlots.cc kFSR.cc Pt_Extrapolation.cc L2ResOutput.cc L2ResAllRuns.cc InputForGlobalFit.cc FinalControlPlots.cc useful_functions.cc tdrstyle_mod15.C `root-config  --cflags --evelibs`

How to execute the program:
./main
