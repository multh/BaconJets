# BaconJets

This is a framework to derive jet energy corrections using dijet events. For the installation of the framework it's required to install SFrame and UHH2 (see step 1. and 2.) which is the core of the framework. In order to execute jobs using the dijet event selection, one needs to install BaconJets (see step 3.).
Unfortunately, the ntuple production for the jet energy calibration process is not implemented within the UHH2 framework. Thus, it's required to add two further steps (step 4. and 5.). Step 4. shows a link to the JMEValidator which is used to produce ntuples. Step 5. shows how to convert the JME ntuples into the BaconJets ntuple format.
Finally, in order to calculate the L2 residuals execute the macros as shown in step 6.

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
Once UHH2 is installed and compiled, get the BaconJets directory and compile it
```
git clone -b nataliias https://github.com/UHH2/BaconJets
cd BaconJets
make
```
In order to start a job, use
```
sframe_main config/Example.xml
```


## 4. Ntuple production with the JMEValidator
In order to produce ntuples, follow the introduction here:
```
https://github.com/blinkseb/JMEValidator/tree/CMSSW_7_4_X
```

## 5. Conversion of the ntuples
Once one have produced nutples in the JME format, one needs to convert the ntuples. For this, a ROOT version with at least the CMSSW74X package is needed. Get the converter using
```
git clone https://github.com/stadie/BaconTrans
```
In the 'test.C' file one can loop over all JME ntuples which need to be converted. In order to convert the ntuples just use
```
root -l test.C
```

## 6. Calculation of the L2 residuals
If not already done before, get all the macros
```
git checkout -b marc
cd extrapolations
```
In order to get the correction factors the 'kFSR.C' and the 'PTextrapolations.C' macro needs to be executed. To execute both use
```
root -l run.C
```
and make sure that the extrapolation macros are commented in.