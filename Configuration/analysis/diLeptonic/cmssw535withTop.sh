#!/bin/sh

# export VO_CMS_SW_DIR=/data/group/top/cmssw 
# source $VO_CMS_SW_DIR/cmsset_default.sh

###### TOP 2012  ########
CMS_version=CMSSW_5_3_5
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

addpkg DataFormats/PatCandidates V06-05-06-03
addpkg PhysicsTools/PatAlgos     V08-09-43

###### Jet Energy Corrections #####

###### Electron ID #####

#electron mva id stuff (following top reference twiki)
cvs co -r V00-00-13 -d EGamma/EGammaAnalysisTools UserCode/EGamma/EGammaAnalysisTools
cd $CMSSW_BASE/src/EGamma/EGammaAnalysisTools/data
cat download.url | xargs wget
cd -


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

addpkg TopQuarkAnalysis


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

#checkdeps -a

# for full memory option of LHAPDF, we NEED to compile ElectroWeakAnalysis/Utilities after scram setup lhapdffull
scram setup lhapdffull
addpkg ElectroWeakAnalysis/Utilities


scram b -j12

