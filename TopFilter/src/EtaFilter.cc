#include "TopAnalysis/TopFilter/interface/EtaFilter.h"


/// default constructor
EtaFilter::EtaFilter(const edm::ParameterSet& cfg): ViewFilter(cfg)
{
}

/// filter (worker class):
bool EtaFilter::filter(const std::vector<edm::View<reco::Candidate> >& objs)
{
  // recieves a vector of objects and filters according to a vector of cuts;
  // first cut acts on the leading object, second on the second leading aso; 
  // passedWeak: at least one cut was passed; passedStrong all cuts were 
  // passed

  bool passedWeak=false, passedStrong=true;
  for(unsigned int jdx=0; jdx<objs.size(); ++jdx){
    bool passedOnce=true;
    // skip if this collection has less members than required
    // by the length of the vectors of min or max 
    if( objs[jdx].size()<min_.size() || objs[jdx].size()<max_.size() )
      passedOnce=false;
    
    unsigned int idx=0;
    for(edm::View<reco::Candidate>::const_iterator obj=objs[jdx].begin();
	obj!=objs[jdx].end(); ++obj) {
      if( idx<min_.size() ) // check for min Eta as long as vector is long enough
	if( !(obj->eta()>min_[idx]) ) passedOnce=false;
      if( idx<max_.size() ) // check for max Eta as long as vector is long enough
	if( !(obj->eta()<max_[idx]) ) passedOnce=false;
      
      // break slope if both vector lengths are exceeded
      ++idx;
      if( idx>min_.size() && idx>max_.size())
	break;
    }
    if(  passedOnce ) passedWeak  =true;
    if( !passedOnce ) passedStrong=false;
  }

  // return filter result
  bool passed = false;
  switch(mode_){
  case kStrong: passed = passedWeak;
  case kWeak  : passed = passedWeak;
  }
  return passed;
}
