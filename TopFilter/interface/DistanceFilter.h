#ifndef DistanceFilter_h
#define DistanceFilter_h

#include <memory>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


template <typename Collection> 
class DistanceFilter {

 public:
  /// choose whether all cut (strong) or at least one cut (weak) should be passed
  enum modes {kStrong, kWeak}; 

  explicit DistanceFilter(const edm::ParameterSet&);
  ~DistanceFilter(){};

  /// filter wrapper and counter  
  bool operator()(edm::Event&, const std::vector<Collection>&, const double&);
  /// filter (worker classes)
  bool filter(const edm::View<reco::Candidate>& , const std::vector<Collection>&);  
  /// event summary of counter
  void summarize();
  
 private:
  
  unsigned int beforeCut_, afterCut_;
  double beforeCutWeighted_, afterCutWeighted_;

  modes mode_;
  std::string name_;
  std::vector<double> min_;
  std::vector<double> max_;

  /// extra information for distance filter
  edm::InputTag refs_;
  std::vector<double> thresh_;
};

/// default constructor
template <typename Collection> 
DistanceFilter<Collection>::DistanceFilter(const edm::ParameterSet& cfg):
  beforeCut_( 0 ), afterCut_( 0 ), beforeCutWeighted_( 0. ), afterCutWeighted_( 0. ),
  name_( cfg.getParameter<std::string>("name") ),
  min_ ( cfg.getParameter<std::vector<double> >("min") ),
  max_ ( cfg.getParameter<std::vector<double> >("max") )
{
  // read extra information for distance filter
  refs_  = cfg.getParameter<edm::InputTag>("refs");
  thresh_= cfg.getParameter<std::vector<double> >("thresh" );

  // choose operation mode (weak or strong filtering)
  switch(cfg.getParameter<unsigned>("mode")){
  case  0: mode_=kStrong; break;
  case  1: mode_=kWeak;   break;
  default: throw edm::Exception(edm::errors::Configuration, 
				"unknown choice for filter mode; choose between 0(strong) or 1(weak) \n");
  }

  // print cuts
  unsigned int maxSize=min_.size();
  if(max_.size()>maxSize) maxSize=max_.size();
  for(unsigned int idx=0; idx<maxSize; ++idx){
    std::cout << ::std::setw( 20 ) << ::std::left;
    if(idx==0) std::cout << name_; else std::cout << " ";
    std::cout << ": ";
    if(idx<min_.size()) std::cout << "min = " << ::std::setw( 8 ) << ::std::right 
				  <<  min_[idx] << "  "; else std::cout << "   "; 
    if(idx<max_.size()) std::cout << "max = " << ::std::setw( 8 ) << ::std::right
				  << max_[idx] 
				  << std::endl;
  }
}

/// counter and filter wrapper
template <typename Collection> 
bool DistanceFilter<Collection>::operator()(edm::Event& evt, const std::vector<Collection>& objs, const double& weight)
{
  edm::Handle<edm::View<reco::Candidate> > refs;
  evt.getByLabel(refs_, refs);
  ++beforeCut_;
  beforeCutWeighted_ += weight;
  if( filter(*refs, objs) ){
    ++afterCut_;
   afterCutWeighted_ += weight;
   return true;
  }
  return false;
}

/// summarize event counting
template <typename Collection> 
void DistanceFilter<Collection>::summarize()
{
  std::cout << ::std::setw( 20 ) << ::std::left  << name_ << " : "
	    << ::std::setw( 10 ) << ::std::right << afterCut_ << " (" 
	    << ::std::setw( 10 ) << ::std::right << afterCutWeighted_  << ") outof "
	    << ::std::setw( 10 ) << ::std::right << beforeCut_<< " (" 
	    << ::std::setw( 10 ) << ::std::right << beforeCutWeighted_ << ")" << std::endl;
}

/// filter (worker class):
template <typename Collection> 
bool DistanceFilter<Collection>::filter(const edm::View<reco::Candidate>& refs, const std::vector<Collection>& objs)
{
  // recieves a vector of objects and filters according to a vector of cuts;
  // first cut acts on the leading object, second on the second leading aso; 
  // passedWeak: at least one cut was passed; passedStrong all cuts were 
  // passed; cut criterion is the minimal distance of the obj to an object 
  // of refs above a certain threshold defined in thresh_; this is caluculated
  // on the fly

  bool passedWeak=false, passedStrong=true;
  for(unsigned int jdx=0; jdx<objs.size(); ++jdx){
    bool passedOnce=true;
    // skip if this collection has less members than required
    // by the length of the vector iso_
    if( objs[jdx].size()<min_.size() || objs[jdx].size()<max_.size() )
      passedOnce=false;
    
    unsigned int idx=0;
    for(typename Collection::const_iterator obj=objs[jdx].begin();
	obj!=objs[jdx].end(); ++obj) {
      // calculate minDR on the fly
      double minDR=-1.;
      for(edm::View<reco::Candidate>::const_iterator ref = refs.begin();
	    ref!=refs.end(); ++ref) {
	// skip jet if too soft as long as vector is long enough
	if( idx<thresh_.size() && ref->pt()<thresh_[idx] ) continue; 
	double dR=deltaR(obj->eta(), obj->phi(), ref->eta(), ref->phi());
	if( minDR<0 || dR<minDR ) minDR=dR;
      }
      if( minDR<0 ){
	edm::LogWarning ( "NoRefObjFound" ) << "no ref object nearest to object found";
	minDR=999.;
      }
      if( idx<min_.size() ) // check for min isolation as long as vector is long enough       
	if( !(minDR>min_[idx]) ) passedOnce=false;
      if( idx<max_.size() ) // check for min isolation as long as vector is long enough       
	if( !(minDR<max_[idx]) ) passedOnce=false;
   
      // break slope if both vector lengths are exceeded
      ++idx;
      if( idx>min_.size() && idx>max_.size())
	break;
    }
    if(  passedOnce ) passedWeak=true;
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

#endif
  
  
