LIBRARY := SUHH2BaconJets
#DICT := include/ReconstructionHypothesis.h include/SUHH2common_LinkDef.h
#TEST := 1
USERLDFLAGS := -lSUHH2core -lSUHH2common -lSUHH2bacondataformats
USERCXXFLAGS := -Wno-reorder -Wno-unused-but-set-variable -Wno-unused-variable -Wno-unused-parameter
PAR := 1
include ../Makefile.common


