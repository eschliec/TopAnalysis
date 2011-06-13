// simple event filter expansions
#include "TopAnalysis/TopFilter/plugins/SimpleEventFilter.h"

#include "TopAnalysis/TopFilter/interface/PtFilter.h"
#include "TopAnalysis/TopFilter/interface/EtaFilter.h"

typedef SimpleEventFilter<PtFilter > PtEventFilter;
typedef SimpleEventFilter<EtaFilter> EtaEventFilter;

// specific event filter expansions
#include "TopAnalysis/TopFilter/interface/IsolationFilterFwd.h"

#include "TopAnalysis/TopFilter/plugins/SpecificEventFilter.h"
#include "TopAnalysis/TopFilter/interface/DistanceFilterFwd.h"

typedef SpecificEventFilter<std::vector<pat::Muon>,     MuonIsolationFilter> MuonIsolationEventFilter;
typedef SpecificEventFilter<std::vector<pat::Muon>,     MuonDistanceFilter > MuonDistanceEventFilter;
typedef SpecificEventFilter<std::vector<pat::Electron>, ElecIsolationFilter> ElecIsolationEventFilter;
typedef SpecificEventFilter<std::vector<pat::Electron>, ElecDistanceFilter > ElecDistanceEventFilter;

// n-1 filter expansions
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TopAnalysis/TopFilter/plugins/SemiLepEventFilter.h"

typedef SemiLepEventFilter<std::vector<pat::Muon> > SemiLepMuonEventFilter;

// plugin definition
#include "FWCore/Framework/interface/MakerMacros.h"

// simple event filters
DEFINE_FWK_MODULE( PtEventFilter            );
DEFINE_FWK_MODULE( EtaEventFilter           );
DEFINE_FWK_MODULE( MuonDistanceEventFilter  );
DEFINE_FWK_MODULE( MuonIsolationEventFilter );
DEFINE_FWK_MODULE( ElecDistanceEventFilter  );
DEFINE_FWK_MODULE( ElecIsolationEventFilter );

// full event filters
DEFINE_FWK_MODULE( SemiLepMuonEventFilter   );

