#!/bin/zsh

ini () {
	eval "`/usr/share/NAF_profiles/ini.pl -b $*`"
}

ini cmssw

alias cmsenv='`scram runtime -sh`'

###### TOP 2011  ########
CMS_version=CMSSW_4_2_5
export SCRAM_ARCH=slc5_amd64_gcc434

scram p -s CMSSW ${CMS_version}
cd ${CMS_version}/src
if [ $? != 0 ]; then
    exit 1
fi

echo export SCRAM_ARCH=${SCRAM_ARCH} > SCRAM_ARCH
chmod u+x SCRAM_ARCH

cmsenv

###### PAT #####

addpkg DataFormats/PatCandidates  V06-04-19
addpkg PhysicsTools/PatAlgos      V08-06-38
addpkg PhysicsTools/PatExamples   V00-05-21
addpkg PhysicsTools/SelectorUtils V00-03-17
addpkg RecoJets/Configuration     V02-04-17
addpkg RecoTauTag/Configuration   V01-02-02
addpkg RecoTauTag/RecoTau         V01-02-06
addpkg RecoTauTag/TauTagTools     V01-02-00
addpkg FWCore/GuiBrowsers         V00-00-56

###### Jet Energy Corrections #####

addpkg JetMETCorrections/Type1MET V04-04-04
echo "Fixing JetMETCorrections BuildFile"
mv JetMETCorrections/Type1MET/BuildFile.xml JetMETCorrections/Type1MET/BuildFile.xml.old
echo '<use   name="DataFormats/PatCandidates"/>' >  JetMETCorrections/Type1MET/BuildFile.xml
cat JetMETCorrections/Type1MET/BuildFile.xml.old >> JetMETCorrections/Type1MET/BuildFile.xml

#cvs co -d JEC UserCode/KKousour/data/Jec11V2.db
#ln -s JEC/Jec11V2.db

###### Electron ID #####
addpkg RecoEgamma/ElectronIdentification V00-03-32
addpkg ElectroWeakAnalysis/WENu V00-04-00

###### Lumi Reweighting #####
addpkg PhysicsTools/Utilities V08-03-10

###### Hamburg TOP package (plus dependencies) #####

cvs co -d TopAnalysis UserCode/Bromo/TopAnalysis
addpkg TopQuarkAnalysis/Configuration V06-01-04
addpkg TopQuarkAnalysis/TopEventSelection V06-07-11

checkdeps -a

# to filter the madgraph/pythia momentum convservation bug
addpkg GeneratorInterface/GenFilters CMSSW_4_2_8_patch7

scram b -j8

