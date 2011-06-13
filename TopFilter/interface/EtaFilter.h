#ifndef EtaFilter_h
#define EtaFilter_h

#include "TopAnalysis/TopFilter/interface/ViewFilter.h"

class EtaFilter : public ViewFilter{

 public:
  explicit EtaFilter(const edm::ParameterSet&);
  ~EtaFilter(){};

  /// filter (worker class)
  bool filter(const std::vector<edm::View<reco::Candidate> >&);
};

#endif
  
