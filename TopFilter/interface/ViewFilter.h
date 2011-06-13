#ifndef ViewFilter_h
#define ViewFilter_h

#include <memory>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


class ViewFilter {

 public:
  /// choose whether all cut (strong) or at least one cut (weak) should be passed
  enum modes {kStrong, kWeak}; 

  explicit ViewFilter(const edm::ParameterSet&);
  virtual ~ViewFilter(){};

  /// filter wrapper and counter
  bool operator()(edm::Event&, const std::vector<edm::View<reco::Candidate> >&, const double&);
  /// event summary of counter
  void summarize();

  // -----------------------------------------
  // implemet the following function
  // for a concrete event filter
  // -----------------------------------------

  /// filter (worker class)
  virtual bool filter(const std::vector<edm::View<reco::Candidate> >&) = 0;
  
 protected:

  unsigned beforeCut_, afterCut_;
  double beforeCutWeighted_, afterCutWeighted_;

  modes mode_;
  std::string name_;
  std::vector<double> min_;
  std::vector<double> max_;
};

#endif
  
