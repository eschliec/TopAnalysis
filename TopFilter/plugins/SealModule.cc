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

#include "TopAnalysis/TopFilter/plugins/JetOverlapEventFilter.h"

// n-1 filter expansions
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TopAnalysis/TopFilter/plugins/SemiLepEventFilter.h"

typedef SemiLepEventFilter<std::vector<pat::Muon> > SemiLepMuonEventFilter;

// plugin definition
#include "FWCore/Framework/interface/MakerMacros.h"

// di-muon-mass filter
#include "TopAnalysis/TopFilter/interface/DiMuonMassFilter.h"

// filter for full leptonic event properties
#include "TopAnalysis/TopFilter/plugins/FullLepHypothesesFilter.h"

// filter to select evts with at least n objects of given type
#include "TopAnalysis/TopFilter/plugins/FullLepObjectCountFilter.h"

// simple event filters
DEFINE_FWK_MODULE( PtEventFilter            );
DEFINE_FWK_MODULE( EtaEventFilter           );
DEFINE_FWK_MODULE( MuonDistanceEventFilter  );
DEFINE_FWK_MODULE( MuonIsolationEventFilter );
DEFINE_FWK_MODULE( ElecDistanceEventFilter  );
DEFINE_FWK_MODULE( ElecIsolationEventFilter );
DEFINE_FWK_MODULE( JetOverlapEventFilter    );

DEFINE_FWK_MODULE( DiMuonMassFilter    );

// full event filters
DEFINE_FWK_MODULE( SemiLepMuonEventFilter   );
DEFINE_FWK_MODULE( FullLepHypothesesFilter  );
DEFINE_FWK_MODULE( FullLepObjectCountFilter  );
