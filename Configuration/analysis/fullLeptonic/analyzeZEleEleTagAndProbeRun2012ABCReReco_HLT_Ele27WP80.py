jetType="PF"
dataSelector = 11
leptonTypeId = 11
triggerPathSelector = "HLT_Ele27_WP80_v*"

jsonFile = "../semiLeptonic/diffXSection/Cert_190456-198523_8TeV_201213Jul06Aug24AugReRecos_Collisions12_JSON.txt"

execfile("analyzeZMuMuTagAndProbe_cfg.py")

#process.load("TopAnalysis.Configuration.samples.Run2011A_SingleMu_PromptRecoV4_AOD_cff")
##  ? files
process.load("TopAnalysis/Configuration/samples/SingleElectron_ReRecoRun2012ABC_cff")

process.GlobalTag.globaltag = cms.string('GR_R_53_V17::All')
print "UPDATED global tag in data cfg file: ", process.GlobalTag.globaltag


# Lumi range
## ATTENTION!!! At the moment myLumis are filled in this data_cfg files again
## as otherwise overwritten by load("data_cff")
if(jsonFile==''):
    print "No JSON file used"
    #process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('160404:1-999999:1')
    #print "lumisToProcess = ", process.source.lumisToProcess
else:
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
    process.source.lumisToProcess.extend(myLumis)
    ## Use lumisToSkip if JSON file shall be constrained additionally to a certain run range
    #process.source.lumisToSkip = cms.untracked.VLuminosityBlockRange('165634:1-999999:1')
    #print "lumisToProcess = ", process.source.lumisToProcess
process.TFileService.fileName = 'analyzeZMuMuTagAndProbeRun2012ABCReReco_HLT_Ele27WP80.root'
process.maxEvents.input = 10000
process.MessageLogger.cerr.FwkReport.reportEvery = 10000