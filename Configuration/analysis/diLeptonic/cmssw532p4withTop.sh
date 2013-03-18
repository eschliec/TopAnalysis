#!/bin/sh

# export VO_CMS_SW_DIR=/data/group/top/cmssw 
# source $VO_CMS_SW_DIR/cmsset_default.sh

###### TOP 2012  ########
CMS_version=CMSSW_5_3_2_patch4
export SCRAM_ARCH=slc5_amd64_gcc462

scram p -s CMSSW ${CMS_version}
cd ${CMS_version}/src
if [ $? != 0 ]; then
    exit 1
fi

echo export SCRAM_ARCH=${SCRAM_ARCH} > SCRAM_ARCH
chmod u+x SCRAM_ARCH

eval `scramv1 runtime -sh`

###### PAT #####

addpkg DataFormats/PatCandidates  V06-05-00
addpkg RecoJets/Configuration     V02-07-02
addpkg RecoTauTag/Configuration   V01-03-04
addpkg RecoTauTag/RecoTau         V01-04-12
addpkg RecoTauTag/TauTagTools     V01-02-03
addpkg RecoParCcleFlow/PFProducer V15-01-11

###### Jet Energy Corrections #####

addpkg JetMETCorrections/Type1MET V04-06-03
echo "Fixing JetMETCorrections BuildFile"
mv JetMETCorrections/Type1MET/BuildFile.xml JetMETCorrections/Type1MET/BuildFile.xml.old
echo '<use   name="DataFormats/PatCandidates"/>' >  JetMETCorrections/Type1MET/BuildFile.xml
cat JetMETCorrections/Type1MET/BuildFile.xml.old >> JetMETCorrections/Type1MET/BuildFile.xml


###### Electron ID #####
addpkg RecoEgamma/ElectronIdentification V00-03-36
addpkg ElectroWeakAnalysis/WENu V00-04-02

###### Lumi Reweighting #####
addpkg PhysicsTools/Utilities V08-03-17
addpkg SimDataFormats/PileupSummaryInfo V00-02-04

###### ECAL Laser rejection in 2012A+B ####
addpkg RecoMET/METFilters V00-00-10

###### Hamburg TOP package (plus dependencies) #####

cvs co -d TopAnalysis UserCode/Bromo/TopAnalysis
cvs update -r Oct2012-Rewritten-Analysis TopAnalysis/TopAnalyzer/python
cvs update -r Oct2012-Rewritten-Analysis TopAnalysis/TopAnalyzer/plugins
cvs update -r B_5_2_X TopAnalysis/TopAnalyzer/python/FullLepGenAnalyzer_cfi.py
cvs update -r B_5_2_X TopAnalysis/TopFilter/python/filters/GeneratorTtDileptonFilter_cfi.py
cvs update -r B_5_2_X TopAnalysis/TopUtils/python/GenLevelBJetProducer_cfi.py
cvs update -r B_5_2_X TopAnalysis/TopUtils/python/genlevelbjetproducer_cfi.py
cvs update -r B_5_2_X TopAnalysis/TopUtils/python/HadronLevelBJetProducer_cfi.py
cvs update -r B_5_2_X TopAnalysis/TopUtils/python/usePatTupleWithParticleFlow_cff.py
cvs update -r Oct2012-Rewritten-Analysis TopAnalysis/Configuration/analysis/diLeptonic

#With the newest tag for this package, problem with elPFIsoValues (the other option is to change the "elPFIsoValue*03/04" in the file
addpkg RecoParticleFlow/PFProducer V15-01-03-02

addpkg TopQuarkAnalysis/Configuration V06-02-00
addpkg TopQuarkAnalysis/TopEventSelection V06-07-14

######## And things to change in the config file 
#In .py I had to comment out this part: Incompatible with Eike's new config file
#                                           'analyzersBeforeMuonIso': cms.Sequence(
#                                                process.unisolatedMuons *
#                                                process.analyzeMuons1),
#                                           'analyzersBeforeElecIso': cms.Sequence(
#                                                process.unisolatedElectrons *
#                                                process.analyzeMuons2*
#                                                process.analyzeElecs2
#                                                ),

#options.register('globalTag', 'START52_V9::All', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "which globalTag should be used") ===> write the globalTag, otherwise it complains about 
#"    from Configuration.PyReleaseValidation.autoCond import autoCond" ---> No module named autoCond
########

# setup the PDF reweighting, need Ewk/Utils to be recompiled to use "FULL MEMORY option" of LHAPDF
scram setup lhapdffull
addpkg ElectroWeakAnalysis/Utilities

#checkdeps -a
scram b -j12

