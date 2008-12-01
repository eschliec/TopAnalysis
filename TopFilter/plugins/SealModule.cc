// simple event filter expansions
#include "TopAnalysis/TopFilter/interface/PtFilter.h"
#include "TopAnalysis/TopFilter/interface/EtaFilter.h"
#include "TopAnalysis/TopFilter/interface/DoubleFilter.h"
#include "TopAnalysis/TopFilter/plugins/SimpleEventFilter.h"
#include "TopAnalysis/TopFilter/plugins/SimpleFilter.h"

typedef SimpleEventFilter<PtFilter > PtEventFilter;
typedef SimpleEventFilter<EtaFilter> EtaEventFilter;
typedef SimpleFilter<DoubleFilter> DoubleEventFilter;

// specific event filter expansions
#include "TopAnalysis/TopFilter/plugins/SpecificEventFilter.h"
#include "TopAnalysis/TopFilter/interface/IsolationFilterFwd.h"
#include "TopAnalysis/TopFilter/interface/JetIsolationFilterFwd.h"
#include "TopAnalysis/TopFilter/plugins/TtSemiLepSignalSelectorMVATrainer.h"
#include "TopAnalysis/TopFilter/plugins/TtSemiLepSignalSelectorMVAComputer.h"

typedef SpecificEventFilter<std::vector<pat::Muon>, MuonIsolationFilter> MuonIsolationEventFilter;
typedef SpecificEventFilter<std::vector<pat::Electron>, ElecIsolationFilter> ElecIsolationEventFilter;
typedef SpecificEventFilter<std::vector<pat::Muon>, MuonJetIsolationFilter> MuonJetIsolationEventFilter;
typedef SpecificEventFilter<std::vector<pat::Electron>, ElecJetIsolationFilter> ElecJetIsolationEventFilter;

// n-1 filter expansions
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TopAnalysis/TopFilter/plugins/SemiLepEventFilter.h"

typedef SemiLepEventFilter<std::vector<pat::Muon> > SemiLepMuonEventFilter;

// plugin definition
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( PtEventFilter  );
DEFINE_FWK_MODULE( EtaEventFilter );
DEFINE_FWK_MODULE( SemiLepMuonEventFilter );
DEFINE_FWK_MODULE( MuonIsolationEventFilter );
DEFINE_FWK_MODULE( ElecIsolationEventFilter );
DEFINE_FWK_MODULE( MuonJetIsolationEventFilter );
DEFINE_FWK_MODULE( ElecJetIsolationEventFilter );
DEFINE_FWK_MODULE( TtSemiLepSignalSelectorMVATrainer );
DEFINE_FWK_MODULE( TtSemiLepSignalSelectorMVAComputer );
DEFINE_FWK_MODULE( DoubleEventFilter );
