# BaconJets

This is a framework to derive jet energy corrections using dijet events. For the installation of the framework it's required to install SFrame and UHH2 (see step 1. and 2.). 

## 1. Install and compile SFrame
First of all one needs to install SFrame.
```
svn co https://svn.code.sf.net/p/sframe/code/SFrame/tags/SFrame-04-00-01/ SFramePhys14
```

Before compiling SFrame, remember to setup CMSSW first (to build SFrame against the right root version) and source `setup.sh` in the SFrame directory.


## 2. Install UHH2

```

cmsrel CMSSW_7_2_1_patch4
cd CMSSW_7_2_1_patch4/src
cmsenv
scram b -j 10

git clone -b https://github.com/UHH2/UHH2.git
```

Before compiling it is required to set up CMSSW ('cmsenv') and SFrame ('source setup.sh'). Afterwards, one can compile SFrame in the SFrame directory using 'make'.

In order to compile the UHH2 framework go to the UHH2 directory and type 'make'. This will compile both, CMSSW and SFrame. Please be aware that it's possible to compile only SFrame with 'make sframe' and only the CMSSW package with 'make scram'.

## 3. Get BaconJets
Once UHH2 is installed and compiled, get the BaconJets directory
```
git clone -b nataliias https://github.com/UHH2/BaconJets
```