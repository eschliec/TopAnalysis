#include "FWCore/Framework/interface/MakerMacros.h"
#include "TopAnalysis/TopAnalyzer/plugins/SingleObjectAnalyzer.h"
#include "TopAnalysis/TopAnalyzer/plugins/DoubleObjectAnalyzer.h"

/**
    Define modules needed for the analysis of the muon
**/

#include "TopAnalysis/TopAnalyzer/interface/MuonQuality.h"
typedef SingleObjectAnalyzer<edm::View<pat::Muon>, MuonQuality> MuonQualityAnalyzer;
DEFINE_FWK_MODULE(MuonQualityAnalyzer);

#include "TopAnalysis/TopAnalyzer/interface/MuonKinematics.h"
typedef SingleObjectAnalyzer<edm::View<pat::Muon>, MuonKinematics> MuonKinematicsAnalyzer;
DEFINE_FWK_MODULE(MuonKinematicsAnalyzer);

#include "TopAnalysis/TopAnalyzer/interface/MuonResolution.h"
typedef SingleObjectAnalyzer<edm::View<pat::Muon>, MuonResolution> MuonResolutionAnalyzer;
DEFINE_FWK_MODULE(MuonResolutionAnalyzer);

#include "TopAnalysis/TopAnalyzer/interface/MuonCrossSection.h"
typedef SingleObjectAnalyzer<edm::View<pat::Muon>, MuonCrossSection> MuonCrossSectionAnalyzer;
DEFINE_FWK_MODULE(MuonCrossSectionAnalyzer);

/**
    Define modules needed for the analysis of the jets
**/

#include "TopAnalysis/TopAnalyzer/interface/JetQuality.h"
typedef SingleObjectAnalyzer<edm::View<pat::Jet>, JetQuality> JetQualityAnalyzer;
DEFINE_FWK_MODULE(JetQualityAnalyzer);

#include "TopAnalysis/TopAnalyzer/interface/JetKinematics.h"
typedef SingleObjectAnalyzer<edm::View<pat::Jet>, JetKinematics> JetKinematicsAnalyzer;
DEFINE_FWK_MODULE(JetKinematicsAnalyzer);

#include "TopAnalysis/TopAnalyzer/interface/JetResolution.h"
typedef SingleObjectAnalyzer<edm::View<pat::Jet>, JetResolution> JetResolutionAnalyzer;
DEFINE_FWK_MODULE(JetResolutionAnalyzer);

#include "TopAnalysis/TopAnalyzer/plugins/JetEnergyCorrectionsAnalyzer.h"
DEFINE_FWK_MODULE(JetEnergyCorrectionsAnalyzer);

/**
    Define modules needed for the analysis of the electron
**/

#include "TopAnalysis/TopAnalyzer/interface/ElectronQuality.h"
typedef SingleObjectAnalyzer<edm::View<pat::Electron>, ElectronQuality> ElectronQualityAnalyzer;
DEFINE_FWK_MODULE(ElectronQualityAnalyzer);

#include "TopAnalysis/TopAnalyzer/interface/ElectronKinematics.h"
typedef SingleObjectAnalyzer<edm::View<pat::Electron>, ElectronKinematics> ElectronKinematicsAnalyzer;
DEFINE_FWK_MODULE(ElectronKinematicsAnalyzer);

/**
    Define modules needed for the analysis of event shapes
**/

#include "TopAnalysis/TopAnalyzer/interface/EventShapes.h"
typedef SingleObjectAnalyzer<edm::View<reco::Candidate>, EventShapes> EventShapeAnalyzer;
DEFINE_FWK_MODULE(EventShapeAnalyzer);

/**
    Define modules needed for the analysis of fully hadronic events
**/

#include "TopAnalysis/TopAnalyzer/interface/FullHadSpecial.h"
typedef SingleObjectAnalyzer<edm::View<pat::Jet>, FullHadSpecial> FullHadSpecialAnalyzer;
DEFINE_FWK_MODULE(FullHadSpecialAnalyzer);

#include "TopAnalysis/TopAnalyzer/interface/FullHadTopReco.h"
typedef DoubleObjectAnalyzer<TtFullHadronicEvent, edm::View<pat::Jet>, FullHadTopReco> FullHadTopRecoAnalyzer;
DEFINE_FWK_MODULE(FullHadTopRecoAnalyzer);


/**
    Define modules needed for the analysis of muon jet relations
**/

#include "TopAnalysis/TopAnalyzer/interface/MuonJetKinematics.h"
typedef DoubleObjectAnalyzer<edm::View<pat::Muon>, edm::View<pat::Jet>, MuonJetKinematics> MuonJetKinematicsAnalyzer;
DEFINE_FWK_MODULE(MuonJetKinematicsAnalyzer);


/**
    Define modules needed for the analysis of different event hypothesis
**/

#include "TopAnalysis/TopAnalyzer/interface/JetCombinatorics.h"
typedef SingleObjectAnalyzer<TtSemiLeptonicEvent, JetCombinatorics> JetCombinatoricsAnalyzer;
DEFINE_FWK_MODULE(JetCombinatoricsAnalyzer);

#include "TopAnalysis/TopAnalyzer/interface/HypothesisKinFit.h"
typedef SingleObjectAnalyzer<TtSemiLeptonicEvent, HypothesisKinFit> HypothesisKinFitAnalyzer;
DEFINE_FWK_MODULE(HypothesisKinFitAnalyzer);

#include "TopAnalysis/TopAnalyzer/interface/HypothesisKinFitJets.h"
typedef DoubleObjectAnalyzer<TtSemiLeptonicEvent, edm::View<pat::Jet>, HypothesisKinFitJets> HypothesisKinFitJetsAnalyzer;
DEFINE_FWK_MODULE(HypothesisKinFitJetsAnalyzer);

#include "TopAnalysis/TopAnalyzer/interface/HypothesisKinFitMuon.h"
typedef DoubleObjectAnalyzer<TtSemiLeptonicEvent, edm::View<pat::Muon>, HypothesisKinFitMuon> HypothesisKinFitMuonAnalyzer;
DEFINE_FWK_MODULE(HypothesisKinFitMuonAnalyzer);

#include "TopAnalysis/TopAnalyzer/interface/HypothesisKinFitMET.h"
typedef DoubleObjectAnalyzer<TtSemiLeptonicEvent, edm::View<pat::MET> , HypothesisKinFitMET> HypothesisKinFitMETAnalyzer;
DEFINE_FWK_MODULE(HypothesisKinFitMETAnalyzer);

#include "TopAnalysis/TopAnalyzer/plugins/FullHadHypothesisAnalyzer.h"
DEFINE_FWK_MODULE(FullHadHypothesisAnalyzer);

#include "TopAnalysis/TopAnalyzer/interface/KinFitQuality.h"
typedef DoubleObjectAnalyzer<TtFullHadronicEvent, edm::View<pat::Jet>, KinFitQuality> KinFitQualityAnalyzer;
DEFINE_FWK_MODULE(KinFitQualityAnalyzer);

#include "TopAnalysis/TopAnalyzer/interface/KinFitImprover.h"
typedef DoubleObjectAnalyzer<TtFullHadronicEvent, edm::View<pat::Jet>, KinFitImprover> KinFitImproveAnalyzer;
DEFINE_FWK_MODULE(KinFitImproveAnalyzer);

/**
    Define modules needed for the analysis of the top quarks
**/

#include "TopAnalysis/TopAnalyzer/interface/TopKinematics.h"
typedef SingleObjectAnalyzer<TtGenEvent, TopKinematics> TopKinematicsGenAnalyzer;
DEFINE_FWK_MODULE(TopKinematicsGenAnalyzer);

typedef SingleObjectAnalyzer<TtSemiLeptonicEvent, TopKinematics> TopKinematicsRecAnalyzer;
DEFINE_FWK_MODULE(TopKinematicsRecAnalyzer);

/**
    Define modules needed for the analysis of b-tags
**/

#include "TopAnalysis/TopAnalyzer/interface/BTags.h"
typedef SingleObjectAnalyzer<edm::View<pat::Jet>, BTags> BTagsAnalyzer;
DEFINE_FWK_MODULE(BTagsAnalyzer);

/**
    Define modules needed for the analysis of generator information
**/

#include "TopAnalysis/TopAnalyzer/interface/GenParticle.h"
typedef SingleObjectAnalyzer<edm::View<reco::GenParticle>, GenParticle> GenParticleAnalyzer;
DEFINE_FWK_MODULE(GenParticleAnalyzer);

/**
    Define modules needed for the analysis of trigger


#include "TopAnalysis/TopAnalyzer/interface/JetTrigger.h"
DEFINE_FWK_MODULE(JetTriggerAnalyzer);
**/
