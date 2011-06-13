#ifndef PtFilter_h
#define PtFilter_h

#include "TopAnalysis/TopFilter/interface/ViewFilter.h"

class PtFilter : public ViewFilter{

 public:
  explicit PtFilter(const edm::ParameterSet&);
  ~PtFilter(){};

  /// filter (worker class)
  bool filter(const std::vector<edm::View<reco::Candidate> >&);
};

#endif
  
