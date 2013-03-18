#!/bin/sh

# export VO_CMS_SW_DIR=/data/group/top/cmssw 
# source $VO_CMS_SW_DIR/cmsset_default.sh

###### TOP 2011  ########
CMS_version=CMSSW_4_2_8_patch7
export SCRAM_ARCH=slc5_amd64_gcc434

scram p -s CMSSW ${CMS_version}
cd ${CMS_version}/src
if [ $? != 0 ]; then
    exit 1
fi

echo export SCRAM_ARCH=${SCRAM_ARCH} > SCRAM_ARCH
chmod u+x SCRAM_ARCH

eval `scramv1 runtime -sh`

###### PAT #####

addpkg DataFormats/PatCandidates  V06-04-19
addpkg RecoJets/Configuration     V02-04-17
addpkg RecoTauTag/Configuration   V01-02-02
addpkg RecoTauTag/RecoTau         V01-02-06
addpkg RecoTauTag/TauTagTools     V01-02-00

###### Jet Energy Corrections #####

addpkg JetMETCorrections/Type1MET V04-04-04
echo "Fixing JetMETCorrections BuildFile"
mv JetMETCorrections/Type1MET/BuildFile.xml JetMETCorrections/Type1MET/BuildFile.xml.old
echo '<use   name="DataFormats/PatCandidates"/>' >  JetMETCorrections/Type1MET/BuildFile.xml
cat JetMETCorrections/Type1MET/BuildFile.xml.old >> JetMETCorrections/Type1MET/BuildFile.xml


###### Electron ID #####
addpkg RecoEgamma/ElectronIdentification V00-03-32
addpkg ElectroWeakAnalysis/WENu V00-04-00

###### Lumi Reweighting #####
addpkg PhysicsTools/Utilities V08-03-17
addpkg SimDataFormats/PileupSummaryInfo V00-02-04

###### Hamburg TOP package (plus dependencies) #####

cvs co -d TopAnalysis UserCode/Bromo/TopAnalysis
addpkg TopQuarkAnalysis/Configuration V06-01-04
addpkg TopQuarkAnalysis/TopEventSelection V06-07-13

### instead of deleting TopAnalysis/TopAnalyzer/src/FullHadQCDEstimation.cc
### in order to get the user code compiling, you could also do:
#addpkg PhysicsTools/KinFitter                V00-03-07
#addpkg TopQuarkAnalysis/TopKinFitter         V06-07-18
#addpkg TopQuarkAnalysis/TopObjectResolutions V06-07-09

checkdeps -a

scram b -j12

