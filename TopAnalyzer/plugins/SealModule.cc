//
// ObjectAnalyzers
//
#include "TopAnalysis/TopAnalyzer/plugins/ObjectAnalyzer.h"

#include "TopAnalysis/TopAnalyzer/interface/JetId.h"
#include "TopAnalysis/TopAnalyzer/interface/JetKinematic.h"
#include "TopAnalysis/TopAnalyzer/interface/JetResolution.h"

typedef ObjectAnalyzer<std::vector<pat::Jet>, JetId, JetKinematic, JetResolution> JetAnalyzer;

#include "TopAnalysis/TopAnalyzer/interface/MuonId.h"
#include "TopAnalysis/TopAnalyzer/interface/MuonKinematic.h"
#include "TopAnalysis/TopAnalyzer/interface/MuonResolution.h"

typedef ObjectAnalyzer<std::vector<pat::Muon>, MuonId, MuonKinematic, MuonResolution> MuonAnalyzer;

#include "TopAnalysis/TopAnalyzer/interface/ElecId.h"
#include "TopAnalysis/TopAnalyzer/interface/ElecKinematic.h"
#include "TopAnalysis/TopAnalyzer/interface/ElecResolution.h"

typedef ObjectAnalyzer<std::vector<pat::Electron>, ElecId, ElecKinematic, ElecResolution> ElecAnalyzer;

//
// EventShapeAnayzers
//
#include "TopAnalysis/TopAnalyzer/plugins/EventShapeAnalyzer.h"

typedef EventShapeAnalyzer<std::vector<pat::Muon> > EventShapeAnalyzerMuon;
typedef EventShapeAnalyzer<std::vector<pat::Electron> > EventShapeAnalyzerElec;


//
// Analyzers for Jet Parton Association
//
#include "TopAnalysis/TopAnalyzer/plugins/PartonMatchAnalyzer.h"
#include "TopAnalysis/TopAnalyzer/plugins/SemiLepHypothesesAnalyzer.h"

//
// Analyzers for QCD BG Studies
//
#include "TopAnalysis/TopAnalyzer/plugins/LeptonNumberAnalyzer.h"
#include "TopAnalysis/TopAnalyzer/plugins/MatrixAnalyzer.h"
#include "TopAnalysis/TopAnalyzer/plugins/IsolationAnalyzer.h"

//
// plugins
//
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetAnalyzer);
DEFINE_FWK_MODULE(MuonAnalyzer);
DEFINE_FWK_MODULE(ElecAnalyzer);
DEFINE_FWK_MODULE(PartonMatchAnalyzer);
DEFINE_FWK_MODULE(SemiLepHypothesesAnalyzer);
DEFINE_FWK_MODULE(EventShapeAnalyzerMuon);
DEFINE_FWK_MODULE(EventShapeAnalyzerElec);
DEFINE_FWK_MODULE(LeptonNumberAnalyzer);
DEFINE_FWK_MODULE(IsolationAnalyzer);
DEFINE_FWK_MODULE(MatrixAnalyzer);
