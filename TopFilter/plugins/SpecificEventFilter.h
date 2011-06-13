#ifndef SpecificEventFilter_h
#define SpecificEventFilter_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


template <typename Collection, typename Filter> 
class SpecificEventFilter : public edm::EDFilter {

 public:

  explicit SpecificEventFilter(const edm::ParameterSet&);
  ~SpecificEventFilter(){};
  
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

template <typename Collection, typename Filter> 
SpecificEventFilter<Collection, Filter>::SpecificEventFilter(const edm::ParameterSet& cfg):
  src_( cfg.getParameter<std::vector<edm::InputTag> >( "input" ) ),
  wgt_( cfg.getParameter<edm::InputTag>( "weight" ) ),
  filter_( cfg.template getParameter<edm::ParameterSet> ("cuts" ) )
{
}

template <typename Collection, typename Filter> 
bool SpecificEventFilter<Collection, Filter>::filter(edm::Event& evt, const edm::EventSetup& setup)
{
  // get weight
  edm::Handle<double> wgt;
  evt.getByLabel(wgt_, wgt);
  // get input objects
  std::vector<Collection> objs;
  for(std::vector<edm::InputTag>::const_iterator tag = src_.begin(); 
      tag!=src_.end(); ++tag){
    edm::Handle<Collection> src;
    evt.getByLabel(*tag, src);
    objs.push_back(*src);
  }
  // perform cuts
  return filter_(evt, objs, *wgt);
}

template <typename Collection, typename Filter> 
void SpecificEventFilter<Collection, Filter>::beginJob(const edm::EventSetup& setup)
{
}

#endif
  
