#ifndef IsolationFilter_h
#define IsolationFilter_h

#include <memory>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"


template <typename Collection> 
class IsolationFilter {

 public:

  explicit IsolationFilter(const edm::ParameterSet&);
  ~IsolationFilter(){};

 public:
  
  bool operator()(edm::Event&, const std::vector<Collection>&, const double&);
  bool filter(const std::vector<Collection>&); 
  void summarize();

 private:

  std::string name_;         // cut name (for monitoring)
  std::vector<double> iso_;  // isolation threshold(s)
  unsigned int mode_;        // choose 0: for trk isolation
                             // choose 1: for cal isolation
			     // choose 2: for combined isolation

 private:

  unsigned int beforeCut_, afterCut_;
  double beforeCutWeighted_, afterCutWeighted_;
};

template <typename Collection> 
IsolationFilter<Collection>::IsolationFilter(const edm::ParameterSet& cfg):
  name_( cfg.getParameter<std::string> ("name") ),
  iso_ ( cfg.getParameter<std::vector<double> >( "isolation" ) ),
  mode_( cfg.getParameter<unsigned int>("mode") ),
  beforeCut_( 0 ), afterCut_( 0 ),
  beforeCutWeighted_( 0. ), afterCutWeighted_( 0. )
{
  // print applied cuts
  for(unsigned int idx=0; idx<iso_.size(); ++idx){
    std::cout << ::std::setw( 20 );
    if( idx  ==0 ) std::cout << name_;
    std::cout << ": ";
    if( mode_==0 ) std::cout << " Trk";
    if( mode_==1 ) std::cout << " Cal"; 
    if( mode_==2 ) std::cout << " Comb";
    if( mode_<2) std::cout << "Isolation < " << iso_[idx] << std::endl;
    else std::cout << "Isolation > " << iso_[idx] << std::endl;
  }
}

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

template <typename Collection> 
bool IsolationFilter<Collection>::filter(const std::vector<Collection>& objs)
{
  bool passedWeak=false, passedStrong=true;
  for(unsigned int jdx=0; jdx<objs.size(); ++jdx){
    bool passedOnce=true;
    // skip if this collection has less members than required
    // by the length of the vector iso_ 
    if( objs[jdx].size()<iso_.size() )
      passedOnce=false;
    
    unsigned int idx=0;
    for(typename Collection::const_iterator obj=objs[jdx].begin();
	obj!=objs[jdx].end(); ++obj) {
      if( idx<iso_.size() ){ // check for isolation as long as vector is long enough
	if( mode_==0 ) // trk isolation
	  if( !(obj->trackIso()<iso_[idx]) ) passedOnce=false;
	if( mode_==1 ) // cal isolation
	  if( !(obj->caloIso() <iso_[idx]) ) passedOnce=false;	  
	if( mode_==2 ) // cal isolation
	  if( !((obj->pt()/(obj->pt()+obj->trackIso()+obj->caloIso())) >iso_[idx]) ) passedOnce=false;	  
	if( !(mode_==0 || mode_==1 || mode_ ==2) ){
	  edm::LogWarning ( "UnknownIsoMode" ) << "Isolation mode: " << mode_ << "is not known";
	  passedOnce=false;
	}
      }
      
      // break slope if both vector lengths are exceeded
      ++idx;
      if( idx>iso_.size() )
	break;
    }
    if(  passedOnce ) passedWeak=true;
    if( !passedOnce ) passedStrong=false;
  }
  return passedWeak;
}
template <typename Collection> 
void IsolationFilter<Collection>::summarize()
{
  std::cout << ::std::setw( 20 ) << ::std::left  << name_ << " : "
	    << ::std::setw( 10 ) << ::std::right << afterCut_ << " (" 
	    << ::std::setw( 10 ) << ::std::right << afterCutWeighted_  << ") outof "
	    << ::std::setw( 10 ) << ::std::right << beforeCut_<< " (" 
	    << ::std::setw( 10 ) << ::std::right << beforeCutWeighted_ << ")" << std::endl;
}

#endif
  
