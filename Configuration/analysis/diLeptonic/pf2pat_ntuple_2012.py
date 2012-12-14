import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import sys
import os


####################################################################
# global job options

MAXEVENTS = -1
REPORTEVERY = 1000
WANTSUMMARY = True

####################################################################

process = cms.Process("pf2patDilepton")
#SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )

####################################################################
# setup command line options
options = VarParsing.VarParsing ('standard')
options.register('runOnMC', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "decide to run on MC or data")
options.register('runOnAOD', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "run on AOD")
options.register('globalTag', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "which globalTag should be used")
options.register('mode', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "which type of analysis to run")
options.register('samplename', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "which sample to run over")
options.register('inputScript', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "python file with input source")
options.register('outputFile', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "root output file")
options.register('systematicsName', 'Nominal', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "type of systematics")
options.register('json', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "limit to certain lumis")
options.register('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "skip N events")
options.register('triggerStudy',False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "do trigger efficiency study")
options.register('includePDFWeights', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "include the PDF weights *slow!!!*")


########################PF2Pat sequence

process.load("Configuration.EventContent.EventContent_cff")
process.out = cms.OutputModule("PoolOutputModule",
    process.FEVTEventContent,
    dataset = cms.untracked.PSet(dataTier = cms.untracked.string('RECO')),
     fileName = cms.untracked.string("eh.root"),
)

process.load("PhysicsTools.PatAlgos.patSequences_cff")

#pfpostfix = "PFlow"
pfpostfix = ""


from PhysicsTools.PatAlgos.tools.pfTools import *



##########Do primary vertex filtering###


from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )

# get and parse the command line arguments
if( hasattr(sys, "argv") ):
    for args in sys.argv :
        arg = args.split(',')
        for val in arg:
            val = val.split('=')
            if(len(val)==2):
                setattr(options,val[0], val[1])

print options.mode
if options.mode == '':
    print 'cannot run without specifying a mode'
    exit(8888)

if options.samplename == '':
    print 'cannot run without specifying a samplename'
    exit(8888)


if options.samplename == 'data':
    options.runOnMC = False

####################################################################
# define input

if options.inputScript != '':
    process.load(options.inputScript)
else:
    print 'need an input script'
    exit(8889)

####################################################################
# limit to json file (if passed as parameter)
if options.runOnMC:
    jetCorr =('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute'])
else:
    jetCorr = ('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute', 'L2L3Residual'])


usePF2PAT(process, runPF2PAT=True, jetAlgo='AK5', runOnMC=options.runOnMC, postfix=pfpostfix, jetCorrections=jetCorr, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'),typeIMetCorrections=True) 

if options.json != '':
    import FWCore.PythonUtilities.LumiList as LumiList
    import FWCore.ParameterSet.Types as CfgTypes
    myLumis = LumiList.LumiList(filename = options.json).getCMSSWString().split(',')
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    process.source.lumisToProcess.extend(myLumis)

if options.skipEvents > 0:
    process.source.skipEvents = cms.untracked.uint32(options.skipEvents)

####################################################################

## configure message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = REPORTEVERY

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(WANTSUMMARY)
)

print options.maxEvents
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

####################################################################
# Geometry and Detector Conditions

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

#there are different global tags for different epochs
#this won't be pretty but I'm not a python guy

if options.globalTag != '':
    process.GlobalTag.globaltag = cms.string( options.globalTag )
else:
    if options.runOnMC:
        process.GlobalTag.globaltag = cms.string('START53_V15::All')
    else:
        if options.outputFile == 'emu_run2012A.root' or options.outputFile == 'ee_run2012A.root' or options.outputFile == 'mumu_run2012A.root' or options.outputFile == 'emu_run2012B.root' or options.outputFile == 'ee_run2012B.root' or options.outputFile == 'mumu_run2012B.root':
            process.GlobalTag.globaltag = cms.string('FT_53_V6_AN3::All')
        elif options.outputFile == 'emu_run2012Arecover.root' or options.outputFile == 'ee_run2012Arecover.root' or options.outputFile == 'mumu_run2012Arecover.root':
            process.GlobalTag.globaltag = cms.string('FT_53_V6C_AN3::All')
        elif  options.outputFile == 'emu_run2012C_24Aug.root' or options.outputFile == 'ee_run2012C_24Aug.root' or options.outputFile == 'mumu_run2012C_24Aug.root':
            process.GlobalTag.globaltag = cms.string('FT53_V10A_AN3::All')
        else:
            process.GlobalTag.globaltag = cms.string('FT_P_V42C_AN3::All')
            
process.load("Configuration.StandardSequences.MagneticField_cff")

####################################################################

# trigger filtering
# get the central diLepton trigger lists
from TopAnalysis.TopFilter.sequences.diLeptonTriggers_cff import *

# setup filter
process.load("TopAnalysis.TopFilter.filters.TriggerFilter_cfi")
process.filterTrigger.TriggerResults = cms.InputTag('TriggerResults','','HLT')
process.filterTrigger.printTriggers = False
if options.mode == 'mumu':
    ttFilterChannelName = 'MuonMuon'
    process.filterTrigger.hltPaths  = mumuTriggers
elif options.mode == 'emu':
    ttFilterChannelName = 'ElectronMuon'
    process.filterTrigger.hltPaths  = emuTriggers
elif options.mode == 'ee':
    ttFilterChannelName = 'ElectronElectron'
    process.filterTrigger.hltPaths  = eeTriggers
else:
    print 'ERROR: unrecognised mode ' + options.mode +'\nuse ee, emu, or mumu'
    exit(8888)

if options.triggerStudy == True:
    process.filterTrigger.hltPaths  = METTriggers


####################################################################
# setup selections for PF2PAT & PAT objects

process.load('EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi')
process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 )
getattr(process,'patElectrons'+pfpostfix).electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
getattr(process, 'patPF2PATSequence'+pfpostfix).replace(getattr(process,'patElectrons'+pfpostfix),
                                              process.eidMVASequence *
                                              getattr(process,'patElectrons'+pfpostfix)
                                              )

    
process.pfIsolatedElectrons.isolationValueMapsCharged = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFId"))
process.pfIsolatedElectrons.deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFId")
process.pfIsolatedElectrons.isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFId"), cms.InputTag("elPFIsoValueGamma03PFId"))
process.patElectrons.isolationValues = cms.PSet( pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFId"), pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFId"), pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFId"), pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFId"), pfPhotons = cms.InputTag("elPFIsoValueGamma03PFId") )

process.selectedPatElectrons.cut = cms.string('electronID("mvaTrigV0") > 0.0')

process.pfIsolatedMuons.cut = cms.string('pt > 20'
                                         '&& abs(eta) < 2.4'
                                         '&& muonRef.isNonnull()'
                                         '&& (muonRef().isGlobalMuon() || muonRef.isTrackerMuon())')

process.selectedPatMuons.cut = cms.string('isPFMuon')

process.pfIsolatedMuons.isolationCut = cms.double(0.20)                 

# setup part running PAT objects

from PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi import *

print 'additional selection on top of PF selection (not used in top projection):'


# skip events (and jet calculation) if event is empty
skipIfNoElectrons = False
skipIfNoMuons     = False
if options.triggerStudy == True:
    print 'doing trigger study, will not reject events without identified leptons'
elif options.mode == 'mumu':
    skipIfNoMuons = True
elif options.mode == 'emu':
    skipIfNoElectrons = True
    skipIfNoMuons     = True
elif options.mode == 'ee':
    skipIfNoElectrons = True
else:
    print 'ERROR: unrecognised mode ' + options.mode +'\nuse ee, emu, or mumu'
    exit(8888)

####################################################################
# basic debugging analyzer

# process.load("TopAnalysis.TopAnalyzer.CheckDiLeptonAnalyzer_cfi")
# process.analyzeDiLepton.electrons = 'fullySelectedPatElectronsCiC'
# process.analyzeDiLepton.muons = 'fullySelectedPatMuons'

####################################################################
# create path

if options.outputFile == '':
    fn = options.mode + '_test.root'
else:
    fn = options.outputFile
print 'Using output file ' + fn

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(fn)
)

### OLD ANALYSIS STARTS HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

zfilter = False
topfilter = False
signal = False
alsoViaTau = False
useGenCutsInTopSignal = True

if options.samplename == 'ttbarsignal':
    topfilter = True
    signal = True
    viaTau = False
       
if options.samplename == 'ttbarsignalviatau':
    topfilter = True
    signal = True
    viaTau = True
       
if options.samplename == 'ttbarsignalplustau':
    topfilter = True
    signal = True
    viaTau = False
    alsoViaTau = True
       
if options.samplename == 'ttbarbg':
    topfilter = True
       
if options.samplename == 'dyee1050':
    zfilter = True
    zfilterValue = 11
    zrange = (10,50)

if options.samplename == 'dymumu1050':
    zfilter = True
    zfilterValue = 13
    zrange = (10,50)

if options.samplename == 'dytautau1050':
    zfilter = True
    zfilterValue = 15
    zrange = (10,50)

if options.samplename == 'dyee50inf':
    zfilter = True
    zfilterValue = 11
    zrange = (50,1e9)

if options.samplename == 'dymumu50inf':
    zfilter = True
    zfilterValue = 13
    zrange = (50,1e9)

if options.samplename == 'dytautau50inf':
    zfilter = True
    zfilterValue = 15
    zrange = (50,1e9)

if signal:
    print "Not skipping if no leptons -- need true level info\n"
    skipIfNoElectrons = False
    skipIfNoMuons = False

#-------------------------------------------------
# process configuration
#-------------------------------------------------

## define which collections and correction you want to be used
isolatedMuonCollection = "selectedPatMuons"
isolatedElecCollection = "selectedPatElectrons"

jetCollection = "hardJets"

metCollection = "scaledJetEnergy:patMETs"

PFElectronGSFMomentumCollection = "scaledJetEnergy:selectedPatElectrons"

#-------------------------------------------------
# modules
#-------------------------------------------------

## detector conditions and magnetic field

## select DY mass window and decay channel (**MC only**)
if zfilter:
        process.load("TopAnalysis.TopFilter.filters.GeneratorZFilter_cfi")
        process.generatorZFilter.zDecayModes = [zfilterValue]
        process.generatorZFilter.diLeptonMassIntervals = zrange

if topfilter:
        process.load("TopAnalysis.TopFilter.filters.GeneratorTopFilter_cfi")
        if signal:
                process.generatorTopFilter.invert_selection = False
                if viaTau:
                        process.generatorTopFilter.channels = [ttFilterChannelName + 'ViaTau']
                elif alsoViaTau:
                        process.generatorTopFilter.channels = [ttFilterChannelName, ttFilterChannelName + 'ViaTau']
                else:
                        process.generatorTopFilter.channels = [ttFilterChannelName]
        else:
                process.generatorTopFilter.channels = [ttFilterChannelName, ttFilterChannelName + 'ViaTau']
                process.generatorTopFilter.invert_selection = True


## produce pat trigger content
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")

###############################################################
# trigger matching for trigger studies
###############################################################

## build cut strings for the matchers
mumuMatchingCut = ''
eeMatchingCut = ''
emuMatchingCut = ''

for trig in mumuTriggers:
    if mumuMatchingCut == '':
        mumuMatchingCut += 'path("'+trig+'")'
    else:
        mumuMatchingCut += ' || path("'+trig+'")'

for trig in eeTriggers:
    if eeMatchingCut == '':
        eeMatchingCut += 'path("'+trig+'")'
    else:
        eeMatchingCut += ' || path("'+trig+'")'

for trig in emuTriggers:
    if emuMatchingCut == '':
        emuMatchingCut += 'path("'+trig+'")'
    else:
        emuMatchingCut += ' || path("'+trig+'")'


print "mumu trigger matching cut string:",mumuMatchingCut
print "ee trigger matching cut string:",eeMatchingCut
print "emu trigger matching cut string:",emuMatchingCut

## define input sources
triggerMatcherInputMuons     = isolatedMuonCollection  # for the moment use what the rest uses, but we have to think about this a bit more
triggerMatcherInputElectrons = isolatedElecCollection  # for the moment use what the rest uses, but we have to think about this a bit more

## generate the different matchers
process.muonTriggerMatchHLTmumu = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                  src     = cms.InputTag( triggerMatcherInputMuons ),
                                                  matched = cms.InputTag( "patTrigger" ),
                                                  matchedCuts = cms.string(mumuMatchingCut),
                                                  maxDPtRel = cms.double( 0.2 ),  #originally: 0.5
                                                  maxDeltaR = cms.double( 0.2 ),  #originally: 0.5
                                                  resolveAmbiguities    = cms.bool( True ),
                                                  resolveByMatchQuality = cms.bool( True )
                                                  )

process.electronTriggerMatchHLTee = process.muonTriggerMatchHLTmumu.clone()
process.electronTriggerMatchHLTee.src = cms.InputTag( triggerMatcherInputElectrons )
process.electronTriggerMatchHLTee.matchedCuts = cms.string(eeMatchingCut)

process.muonTriggerMatchHLTemu = process.muonTriggerMatchHLTmumu.clone()
process.muonTriggerMatchHLTemu.matchedCuts = cms.string(emuMatchingCut)

process.electronTriggerMatchHLTemu = process.electronTriggerMatchHLTee.clone()
process.electronTriggerMatchHLTemu.matchedCuts = cms.string(emuMatchingCut)

process.patTriggerEvent.patTriggerMatches = [ "muonTriggerMatchHLTmumu",
                                              "electronTriggerMatchHLTee",
                                              "muonTriggerMatchHLTemu", "electronTriggerMatchHLTemu"
                                              ]
# embed the results
process.triggerMatchedMuons = cms.EDProducer( "PATTriggerMatchMuonEmbedder",
                                              src = cms.InputTag( triggerMatcherInputMuons ),
                                              matches = cms.VInputTag( "muonTriggerMatchHLTmumu", "muonTriggerMatchHLTemu" )
                                              )

process.triggerMatchedElectrons = cms.EDProducer( "PATTriggerMatchElectronEmbedder",
                                                  src = cms.InputTag( triggerMatcherInputElectrons ),
                                                  matches = cms.VInputTag( "electronTriggerMatchHLTee", "electronTriggerMatchHLTemu" )
                                                  )

# sequence to run
process.triggerMatchers = cms.Sequence( process.muonTriggerMatchHLTmumu
                                        + process.electronTriggerMatchHLTee
                                        + process.muonTriggerMatchHLTemu + process.electronTriggerMatchHLTemu
                                        )

process.triggerEmbedders = cms.Sequence( process.triggerMatchedMuons
                                         + process.triggerMatchedElectrons
                                         )

process.triggerMatching = cms.Sequence( process.patTrigger
                                        * process.triggerMatchers
                                        * process.patTriggerEvent
                                        * process.triggerEmbedders )

###############################################################


## Build Jet Collections
###########################################################
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *

#-------------------------------------------------
# jet selection
#-------------------------------------------------
process.load("TopAnalysis.TopUtils.JetEnergyScale_cfi")

process.load("TopAnalysis.TopFilter.filters.JetIdFunctorFilter_cfi")
process.goodIdJets.jets    = cms.InputTag("scaledJetEnergy:selectedPatJets")
process.goodIdJets.jetType = cms.string('PF')
process.goodIdJets.version = cms.string('FIRSTDATA')
process.goodIdJets.quality = cms.string('LOOSE')

process.hardJets = selectedPatJets.clone(src = 'goodIdJets', cut = 'pt > 20 & abs(eta) < 2.4') 
process.buildJets = cms.Sequence(process.scaledJetEnergy * process.goodIdJets * process.hardJets)

## Lepton-Vertex matching
process.load("TopAnalysis.TopFilter.filters.LeptonVertexFilter_cfi")
process.filterLeptonVertexDistance.muons = isolatedMuonCollection
process.filterLeptonVertexDistance.elecs = PFElectronGSFMomentumCollection 


from TopAnalysis.TopFilter.filters.DiLeptonFilter_cfi import *
process.filterOppositeCharge = filterLeptonPair.clone(
    electrons    = PFElectronGSFMomentumCollection,
    muons        = isolatedMuonCollection,
    Cut          = (0.,0.),
    filterCharge = -1,
)

from PhysicsTools.PatAlgos.selectionLayer1.leptonCountFilter_cfi import *
process.filterChannel =  countPatLeptons.clone()
process.filterChannel.electronSource    = 'filterOppositeCharge'
process.filterChannel.muonSource        = 'filterOppositeCharge'
process.filterChannel.minNumber         = 2
process.filterChannel.countTaus         = False

leptons3 = 'filterDiLeptonMassQCDveto'
finalLeptons = 'filterDiLeptonMassQCDveto'
if options.mode == 'ee':
    process.filterChannel.countElectrons    = True
    process.filterChannel.countMuons        = False
elif options.mode == 'mumu':
    process.filterChannel.countElectrons    = False
    process.filterChannel.countMuons        = True
else:
    process.filterChannel.minNumber         = 1
    process.filterChannel1 = process.filterChannel.clone()
    process.filterChannel2 = process.filterChannel1.clone()
    process.filterChannel1.countElectrons    = True
    process.filterChannel1.countMuons        = False
    process.filterChannel2.countElectrons    = False
    process.filterChannel2.countMuons        = True
    process.filterChannel = cms.Sequence(process.filterChannel1 * process.filterChannel2)

process.filterDiLeptonMassQCDveto           = filterLeptonPair.clone()
process.filterDiLeptonMassQCDveto.muons     = 'filterOppositeCharge'
process.filterDiLeptonMassQCDveto.electrons = 'filterOppositeCharge'
process.filterDiLeptonMassQCDveto.Cut       = (0.,12.)

##Write Ntuple
from TopAnalysis.TopAnalyzer.NTupleWriter_cfi import writeNTuple

writeNTuple.sampleName = options.samplename
writeNTuple.channelName = options.mode
writeNTuple.systematicsName = options.systematicsName
writeNTuple.isMC = options.runOnMC
writeNTuple.isTtBarSample = signal
writeNTuple.includePDFWeights = options.includePDFWeights
writeNTuple.pdfWeights = "pdfWeights:cteq66"

process.writeNTuple = writeNTuple.clone(
    muons=isolatedMuonCollection,
    elecs=PFElectronGSFMomentumCollection,
    jets=jetCollection,
    met=metCollection,
    genMET="genMetTrue",
)

if options.triggerStudy:
    process.writeNTuple.muons = 'triggerMatchedMuons'
    process.writeNTuple.elecs = 'triggerMatchedElectrons'

if options.includePDFWeights:
    if not signal:
        print "PDF variations only supported for the signal"
        exit(5615)
    process.pdfWeights = cms.EDProducer("PdfWeightProducer",
                # Fix POWHEG if buggy (this PDF set will also appear on output,
                # so only two more PDF sets can be added in PdfSetNames if not "")
                #FixPOWHEG = cms.untracked.string("cteq66.LHgrid"),
                GenTag = cms.untracked.InputTag("genParticles"),
                PdfInfoTag = cms.untracked.InputTag("generator"),
                PdfSetNames = cms.untracked.vstring(
                        "cteq66.LHgrid"
                        #, "MRST2006nnlo.LHgrid"
                        #, "NNPDF10_100.LHgrid"
                        #"cteq6mE.LHgrid"
                        # ,"cteq6m.LHpdf"
                        #"cteq6m.LHpdf"
                ))
else:
    process.pdfWeights = cms.Sequence()


## std sequence to produce the ttFullLepEvent
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttFullLepEvtBuilder_cff")
from TopQuarkAnalysis.TopEventProducers.sequences.ttFullLepEvtBuilder_cff import *
if not signal:
    removeTtFullLepHypGenMatch(process)

setForAllTtFullLepHypotheses(process,"muons"    ,finalLeptons)
setForAllTtFullLepHypotheses(process,"electrons",finalLeptons)
setForAllTtFullLepHypotheses(process,"jets"     ,jetCollection)
setForAllTtFullLepHypotheses(process,"mets"     ,metCollection)
if options.runOnMC:
    setForAllTtFullLepHypotheses(process,"jetCorrectionLevel","L3Absolute")
    print "L3Absolute"
else:
    setForAllTtFullLepHypotheses(process,"jetCorrectionLevel","L2L3Residual")
    print "L2L3Residual"
setForAllTtFullLepHypotheses(process,"maxNJets",-1)

#use this?
#process.ttFullLepJetPartonMatch.jets = hardJets

process.kinSolutionTtFullLepEventHypothesis.maxNComb = -1
process.kinSolutionTtFullLepEventHypothesis.searchWrongCharge = True
process.kinSolutionTtFullLepEventHypothesis.tmassbegin = 100.0
process.kinSolutionTtFullLepEventHypothesis.tmassend   = 300.0
process.kinSolutionTtFullLepEventHypothesis.neutrino_parameters = (30.641, 57.941, 22.344, 57.533, 22.232)
#according to our MC 8 TeV values are:
#nu    mpv 40.567 sigma = 16.876
#nubar mpv 40.639 sigma = 17.021

process.kinSolutionTtFullLepEventHypothesis.mumuChannel = False
process.kinSolutionTtFullLepEventHypothesis.emuChannel  = False
process.kinSolutionTtFullLepEventHypothesis.eeChannel = False
if options.mode == 'mumu':
    process.kinSolutionTtFullLepEventHypothesis.mumuChannel = True
    process.ttFullLepEvent.decayChannel1 = cms.int32(2)
    process.ttFullLepEvent.decayChannel2 = cms.int32(2)
elif options.mode == 'emu':
    process.kinSolutionTtFullLepEventHypothesis.emuChannel = True
    process.ttFullLepEvent.decayChannel1 = cms.int32(1)
    process.ttFullLepEvent.decayChannel2 = cms.int32(2)
elif options.mode == 'ee':
    process.kinSolutionTtFullLepEventHypothesis.eeChannel = True
    process.ttFullLepEvent.decayChannel1 = cms.int32(1)
    process.ttFullLepEvent.decayChannel2 = cms.int32(1)


#-------------------------------------------------
# analysis path
#-------------------------------------------------

if zfilter:
        process.zsequence = cms.Sequence(process.generatorZFilter)
else:
        process.zsequence = cms.Sequence()


if topfilter:
	process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
	process.load("TopAnalysis.TopUtils.HadronLevelBJetProducer_cfi")

        process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi") # supplies PDG ID to real name resolution of MC particles, necessary for GenLevelBJetProducer
	process.load("TopAnalysis.TopUtils.GenLevelBJetProducer_cfi")
        process.produceGenLevelBJets.deltaR = 5.0
        process.produceGenLevelBJets.noBBbarResonances = True

	process.decaySubset.fillMode = "kME" # Status3, use kStable for Status2
        process.topsequence = cms.Sequence(
	        process.makeGenEvt *
		process.generatorTopFilter *
	        process.produceHadronLevelBJets *
                process.produceGenLevelBJets
	)

else:
        process.topsequence = cms.Sequence()



if signal:
    process.ntupleInRecoSeq = cms.Sequence()
else:
    if options.triggerStudy:
        process.ntupleInRecoSeq = cms.Sequence(process.triggerMatching * process.writeNTuple)
    else:
        process.ntupleInRecoSeq = cms.Sequence(process.writeNTuple)

process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
process.HBHENoiseFilter.minIsolatedNoiseSumE = cms.double(999999.)
process.HBHENoiseFilter.minNumIsolatedNoiseChannels = cms.int32(999999)
process.HBHENoiseFilter.minIsolatedNoiseSumEt = cms.double(999999.)

process.scrapingFilter = cms.EDFilter( "FilterOutScraping"
                                                 , applyfilter = cms.untracked.bool( True )
                                                 , debugOn     = cms.untracked.bool( False )
                                                 , numtrack    = cms.untracked.uint32( 10 )
                                                 , thresh      = cms.untracked.double( 0.25 )
                                                 )



process.p = cms.Path(
    process.goodOfflinePrimaryVertices *
    getattr(process,'patPF2PATSequence'+pfpostfix) *
    process.buildJets                     *
    process.filterOppositeCharge          *
    process.filterChannel                 *
    process.filterDiLeptonMassQCDveto     *
    process.makeTtFullLepEvent            *
    process.ntupleInRecoSeq
)

if signal:
    if options.triggerStudy:
        process.pNtuple = cms.Path(
            process.goodOfflinePrimaryVertices *
            getattr(process,'patPF2PATSequence'+pfpostfix) *
            process.buildJets *
            process.triggerMatching *
            process.writeNTuple
            )
    else:
        process.pNtuple = cms.Path(
            process.goodOfflinePrimaryVertices *
            getattr(process,'patPF2PATSequence'+pfpostfix) *
            process.buildJets *
            process.writeNTuple
            )


####################################################################
# prepend PF2PAT

from TopAnalysis.TopAnalyzer.CountEventAnalyzer_cfi import countEvents
process.EventsBeforeSelection = countEvents.clone()
process.EventsBeforeSelection.includePDFWeights = options.includePDFWeights
process.EventsBeforeSelection.pdfWeights = "pdfWeights:cteq66"
    

pathnames = process.paths_().keys()
print 'prepending trigger sequence to paths:', pathnames
for pathname in pathnames:
    getattr(process, pathname).insert(0, cms.Sequence(
        process.pdfWeights *
        process.EventsBeforeSelection *
        process.topsequence *
        process.zsequence *
        process.filterTrigger
        ))
if signal:
    process.pNtuple.remove(process.filterTrigger)

process.scaledJetEnergy.inputElectrons       = cms.InputTag("selectedPatElectrons")
process.scaledJetEnergy.inputJets            = cms.InputTag("selectedPatJets")
process.scaledJetEnergy.inputMETs            = cms.InputTag("patMETs")
process.scaledJetEnergy.JECUncSrcFile        = cms.FileInPath("TopAnalysis/TopUtils/data/Fall12_V6_DATA_UncertaintySources_AK5PFchs.txt")
process.scaledJetEnergy.scaleType = "abs"   #abs = 1, jes:up, jes:down

if options.runOnMC:
    process.scaledJetEnergy.resolutionEtaRanges  = cms.vdouble(0, 0.5, 0.5, 1.1, 1.1, 1.7, 1.7, 2.3, 2.3, -1)
    process.scaledJetEnergy.resolutionFactors    = cms.vdouble(1.052, 1.057, 1.096, 1.134, 1.288) # JER standard

    #please change this on the top where the defaults for the VarParsing are given
    if options.systematicsName == "JES_UP":
        process.scaledJetEnergy.scaleType = "jes:up"
    if options.systematicsName == "JES_DOWN":
        process.scaledJetEnergy.scaleType = "jes:down"
    if options.systematicsName == "JER_UP":
        process.scaledJetEnergy.resolutionFactors = cms.vdouble(1.114, 1.113, 1.159, 1.221, 1.443)
    if options.systematicsName == "JER_DOWN":
        process.scaledJetEnergy.resolutionFactors = cms.vdouble(0.991, 1.002, 1.034, 1.049, 1.135)
    

else:
    process.scaledJetEnergy.resolutionEtaRanges  = cms.vdouble(0, -1)
    process.scaledJetEnergy.resolutionFactors    = cms.vdouble(1.0) # JER standard
    process.load("RecoMET.METFilters.ecalLaserCorrFilter_cfi")
    for pathname in pathnames:
        getattr(process, pathname).replace(process.goodOfflinePrimaryVertices,
                                           process.HBHENoiseFilter * process.scrapingFilter * process.ecalLaserCorrFilter * process.goodOfflinePrimaryVertices)
        
process.load("TopAnalysis.TopUtils.SignalCatcher_cfi")

