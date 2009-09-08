#include "FWCore/Framework/interface/MakerMacros.h"
#include "TopAnalysis/TopAnalyzer/plugins/SingleObjectAnalyzer.h"
#include "TopAnalysis/TopAnalyzer/plugins/DoubleObjectAnalyzer.h"


#include "TopAnalysis/TopAnalyzer/interface/MuonQuality.h"
typedef SingleObjectAnalyzer<std::vector<pat::Muon>, MuonQuality> MuonQualityAnalyzer;
DEFINE_FWK_MODULE(MuonQualityAnalyzer);

#include "TopAnalysis/TopAnalyzer/interface/MuonKinematics.h"
typedef SingleObjectAnalyzer<std::vector<pat::Muon>, MuonKinematics> MuonKinematicsAnalyzer;
DEFINE_FWK_MODULE(MuonKinematicsAnalyzer);

#include "TopAnalysis/TopAnalyzer/interface/MuonResolution.h"
typedef SingleObjectAnalyzer<std::vector<pat::Muon>, MuonResolution> MuonResolutionAnalyzer;
DEFINE_FWK_MODULE(MuonResolutionAnalyzer);

#include "TopAnalysis/TopAnalyzer/interface/MuonJetKinematics.h"
typedef DoubleObjectAnalyzer<std::vector<pat::Muon>, std::vector<pat::Jet>, MuonJetKinematics> MuonJetKinematicsAnalyzer;
DEFINE_FWK_MODULE(MuonJetKinematicsAnalyzer);


#include "TopAnalysis/TopAnalyzer/interface/TopKinematics.h"
typedef SingleObjectAnalyzer<TtGenEvent, TopKinematics> TopKinematicsGenAnalyzer;
typedef SingleObjectAnalyzer<TtSemiLeptonicEvent, TopKinematics> TopKinematicsRecAnalyzer;
DEFINE_FWK_MODULE(TopKinematicsGenAnalyzer);
DEFINE_FWK_MODULE(TopKinematicsRecAnalyzer);

#include "TopAnalysis/TopAnalyzer/plugins/FullLepHypothesesAnalyzer.h"
DEFINE_FWK_MODULE(FullLepHypothesesAnalyzer);
