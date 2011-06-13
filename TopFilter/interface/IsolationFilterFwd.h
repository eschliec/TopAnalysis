#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

// isolation filter expansions for muons and electrons
#include "TopAnalysis/TopFilter/interface/IsolationFilter.h"

typedef IsolationFilter   <std::vector<pat::Muon>     > MuonIsolationFilter;
typedef IsolationFilter   <std::vector<pat::Electron> > ElecIsolationFilter;
