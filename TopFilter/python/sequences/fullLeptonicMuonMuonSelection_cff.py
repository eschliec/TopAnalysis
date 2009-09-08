import FWCore.ParameterSet.Config as cms

## jet selector
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
## muon selector
from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi import *
## electron selector
from PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi import *

## jet count filter
from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *
## muon count filter
from PhysicsTools.PatAlgos.selectionLayer1.muonCountFilter_cfi import *
## electron count filter
from PhysicsTools.PatAlgos.selectionLayer1.electronCountFilter_cfi import *

## tight muon selection
tightMuons = selectedLayer1Muons.clone(src = 'selectedLayer1Muons', 
                                       cut = 'combinedMuon.isNull = 0'
					     '& pt > 20.'
					     '& abs(eta) < 2.1'    
                			     '& track.numberOfValidHits >= 11' 
					     '& abs(track.d0) < 0.2'
					     '& combinedMuon.normalizedChi2 < 10.0'
					     '& ecalIsoDeposit.candEnergy < 4'
					     '& hcalIsoDeposit.candEnergy < 6'
					     '& (trackIso+caloIso)/pt < 0.1' )

## tight electron selection for veto
looseElectrons = selectedLayer1Electrons.clone(src = 'selectedLayer1Electrons', 
                                               cut = 'et > 20.' 
					             '& abs(eta) < 2.4' 
					             '& (trackIso+caloIso)/et <  0.2')

  
## tight jet selection
bJetCandidates = selectedLayer1Jets.clone(src = 'selectedLayer1Jets', 
                                          cut = 'pt > 40.' 
					        '& abs(eta) < 2.3' 
						'& 0.05 < emEnergyFraction' 
						'& emEnergyFraction < 0.95')

## at least two tight muons
muonSelection  = countLayer1Muons.clone(src = 'tightMuons', minNumber = 2, maxNumber = 999)
## no additional tight electron
electronVeto   = countLayer1Electrons.clone(src = 'looseElectrons', maxNumber = 0)
## at least two jets
jetSelection = countLayer1Jets.clone(src = 'bJetCandidates', minNumber = 2)

fullLeptonicMuonMuonSelection = cms.Sequence((tightMuons      * muonSelection ) *
                                             (looseElectrons  * electronVeto  ) *
                                             (bJetCandidates  * jetSelection  )
                                            )
