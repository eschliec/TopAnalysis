#ifndef IsolationFilter_h
#define IsolationFilter_h

#include <memory>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


template <typename Collection> 
class IsolationFilter {

 public:
  /// choose whether all cut (strong) or at least one cut (weak) should be passed
  enum modes {kStrong, kWeak}; 
  /// choose from these filter types
  enum types {kAbsoluteTrack, kRelativeTrack, kAbsoluteCalo,  kRelativeCalo, kRelCombIso}; 

  explicit IsolationFilter(const edm::ParameterSet&);
  ~IsolationFilter(){};

  /// filter wrapper and counter  
  bool operator()(edm::Event&, const std::vector<Collection>&, const double&);
  /// filter (worker classes)
  bool filter(const std::vector<Collection>&); 
  /// event summary of counter
  void summarize();

 private:

  unsigned int beforeCut_, afterCut_;
  double beforeCutWeighted_, afterCutWeighted_;

  modes mode_;
  types type_;
  std::string name_;
  std::vector<double> min_;
  std::vector<double> max_;
};

/// default constructor
template <typename Collection> 
IsolationFilter<Collection>::IsolationFilter(const edm::ParameterSet& cfg):
  beforeCut_( 0 ), afterCut_( 0 ), beforeCutWeighted_( 0. ), afterCutWeighted_( 0. ),
  name_( cfg.getParameter<std::string>("name") ),
  min_ ( cfg.getParameter<std::vector<double> >("min") ),
  max_ ( cfg.getParameter<std::vector<double> >("max") )
{
  // choose fiter type (rel or abs isolation)
  switch(cfg.getParameter<unsigned>("type")){
  case  0: type_=kAbsoluteTrack; break;
  case  1: type_=kRelativeTrack; break;
  case  2: type_=kAbsoluteCalo;  break;
  case  3: type_=kRelativeCalo;  break;
  case  4: type_=kRelCombIso;    break;
  default: throw edm::Exception(edm::errors::Configuration, 
				"unknown choice for filter type; choose between 0(absolute track), 1(relative track), 2(absolute calo), 3(relative calo), 4(relative combined isolation) \n");
  }

  // choose operation mode (weak or strong filtering)
  switch(cfg.getParameter<unsigned>("mode")){
  case  0: mode_=kStrong; break;
  case  1: mode_=kWeak;   break;
  default: throw edm::Exception(edm::errors::Configuration, 
				"unknown choice for filter mode; choose between 0(strong) or 1(weak) \n");
  }

  // print cuts
  unsigned int maxSize=min_.size();
  edm::LogVerbatim log("topFilter");  
  if(max_.size()>maxSize) maxSize=max_.size();
  for(unsigned int idx=0; idx<maxSize; ++idx){
    log << ::std::setw( 20 ) << ::std::left;
    if(idx==0) log << name_; else log << " ";
    log << ": ";
    if(idx<min_.size()) log << "min = " << ::std::setw( 8 ) << ::std::right 
		            <<  min_[idx] << "  "; else std::cout << "   "; 
    if(idx<max_.size()) log << "max = " << ::std::setw( 8 ) << ::std::right
			    << max_[idx] 
			    << "\n";
  }
}

/// counter and filter wrapper
template <typename Collection> 
bool IsolationFilter<Collection>::operator()(edm::Event& evt, const std::vector<Collection>& objs, const double& weight)
{
  ++beforeCut_;
  beforeCutWeighted_ += weight;
  if( filter(objs) ) {
    ++afterCut_;
    afterCutWeighted_ += weight;
    return true;
  }
  return false; 
}

/// summarize event counting
template <typename Collection> 
void IsolationFilter<Collection>::summarize()
{
  // if no event weight or dummy weight is in 
  // use a simplified summary format is taken
  edm::LogVerbatim log("topFilter");
  if(beforeCut_ != beforeCutWeighted_){
  log << ::std::setw( 20 ) << ::std::left  << name_ << " : "
      << ::std::setw( 10 ) << ::std::right << afterCut_ << " (" 
      << ::std::setw( 10 ) << ::std::right << afterCutWeighted_  << ") out of "
      << ::std::setw( 10 ) << ::std::right << beforeCut_<< " (" 
      << ::std::setw( 10 ) << ::std::right << beforeCutWeighted_ << ")" << "\n";
  }
  else{
    log << ::std::setw( 20 ) << ::std::left  << name_ << " : "
	<< ::std::setw( 10 ) << ::std::right << afterCut_ << "  out of "
	<< ::std::setw( 10 ) << ::std::right << beforeCut_ << "\n";  
  }
}

/// filter (worker class):
template <typename Collection> 
bool IsolationFilter<Collection>::filter(const std::vector<Collection>& objs)
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
    for(typename Collection::const_iterator obj=objs[jdx].begin();
	obj!=objs[jdx].end(); ++obj) {
      switch(type_){
      case kAbsoluteTrack:
	if( idx<min_.size() ) // check for min iso as long as vector is long enough
	  if( !(obj->trackIso()>min_[idx]) ) 
	    passedOnce=false;
	if( idx<max_.size() ) // check for max iso as long as vector is long enough
	  if( !(obj->trackIso()<max_[idx]) ) 
	    passedOnce=false;
	break;
      case kRelativeTrack:
	if( idx<min_.size() ) // check for min iso as long as vector is long enough
	  if( !(obj->trackIso()/obj->pt()>min_[idx]) ) 
	    passedOnce=false;
	if( idx<max_.size() ) // check for max iso as long as vector is long enough
	  if( !(obj->trackIso()/obj->pt()<max_[idx]) ) 
	    passedOnce=false;
	break;
      case kAbsoluteCalo:
	if( idx<min_.size() ) // check for min iso as long as vector is long enough
	  if( !(obj->caloIso()>min_[idx]) ) 
	    passedOnce=false;
	if( idx<max_.size() ) // check for max iso as long as vector is long enough
	  if( !(obj->caloIso()<max_[idx]) ) 
	    passedOnce=false;
	break;
      case kRelativeCalo:
	if( idx<min_.size() ) // check for min iso as long as vector is long enough
	  if( !(obj->caloIso()/obj->pt()>min_[idx]) ) 
	    passedOnce=false;
	if( idx<max_.size() ) // check for max iso as long as vector is long enough
	  if( !(obj->caloIso()/obj->pt()<max_[idx]) ) 
	    passedOnce=false;
	break;
      case kRelCombIso:
	if( idx<min_.size() ) // check for min iso as long as vector is long enough
	  if( !( (obj->caloIso()+obj->trackIso())/(obj->pt() )>min_[idx]) ) 
	    passedOnce=false;
	if( idx<max_.size() ) // check for max iso as long as vector is long enough
	  if( !( (obj->caloIso()+obj->trackIso())/(obj->pt() )<max_[idx]) ) 
	    passedOnce=false;
	break;
      default:
	edm::LogWarning ( "UnknownIsoMode" ) << "Isolation mode: " << mode_ << "is not known";
	passedOnce=false;
	break;
      }
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
