#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

// isolation filter expansions for muons and electrons
#include "TopAnalysis/TopFilter/interface/DistanceFilter.h"

typedef DistanceFilter<std::vector<pat::Muon>     > MuonDistanceFilter;
typedef DistanceFilter<std::vector<pat::Electron> > ElecDistanceFilter;
