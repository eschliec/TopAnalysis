#ifndef DoubleFilter_h
#define DoubleFilter_h

#include <memory>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"


class DoubleFilter {

 public:

  explicit DoubleFilter(const edm::ParameterSet&);
  ~DoubleFilter(){};

 public:

  bool operator()(edm::Event&, const std::vector<double>&, const double&);
  bool filter(const std::vector<double>&);
  void summarize();

 private:

  std::string name_;
  std::vector<double> min_;
  std::vector<double> max_;

 private:

  unsigned int beforeCut_, afterCut_;
  double beforeCutWeighted_, afterCutWeighted_;
};

#endif

