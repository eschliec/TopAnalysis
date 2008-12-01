#ifndef SimpleFilter_h
#define SimpleFilter_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


template <typename Filter>
class SimpleFilter : public edm::EDFilter {

 public:

  explicit SimpleFilter(const edm::ParameterSet&);
  ~SimpleFilter(){};

 private:

  virtual void beginJob(const edm::EventSetup&);
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){filter_.summarize();};

 private:

  std::vector<edm::InputTag> src_;

  edm::InputTag wgt_;

 private:

  Filter filter_;
};

template <typename Filter>
SimpleFilter<Filter>::SimpleFilter(const edm::ParameterSet& cfg):
  src_( cfg.getParameter<std::vector<edm::InputTag> >( "input" ) ),
  wgt_( cfg.getParameter<edm::InputTag>( "weight" ) ),
  filter_( cfg.template getParameter<edm::ParameterSet> ("cuts" ) )
{
}

template <typename Filter>
bool SimpleFilter<Filter>::filter(edm::Event& evt, const edm::EventSetup& setup)
{
  // get weight
  edm::Handle<double> wgt;
  evt.getByLabel(wgt_, wgt);
  // get input objects
  std::vector<double> objs;
  for(std::vector<edm::InputTag>::const_iterator tag = src_.begin();
      tag!=src_.end(); ++tag){
    edm::Handle<double> src;
    evt.getByLabel(*tag, src);
    objs.push_back(*src);
  }
  // perform cuts
  return filter_(evt, objs, *wgt);
}

template <typename Filter>
void SimpleFilter<Filter>::beginJob(const edm::EventSetup& setup)
{
}

#endif
