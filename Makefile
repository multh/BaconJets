LIBRARY := SUHH2BaconJets
#DICT := include/ReconstructionHypothesis.h include/SUHH2common_LinkDef.h
#TEST := 1
USERLDFLAGS := -lSUHH2core -lSUHH2bacondataformats
include ../Makefile.common

USERCXXFLAGS := -Wno-reorder -Wno-unused-but-set-variable -Wno-unused-variable -Wno-unused-parameter
