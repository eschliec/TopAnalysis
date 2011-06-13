#include "TopAnalysis/TopFilter/interface/ViewFilter.h"
#include "FWCore/Utilities/interface/EDMException.h"

/// default constructor
ViewFilter::ViewFilter(const edm::ParameterSet& cfg):
  beforeCut_( 0 ), afterCut_( 0 ), beforeCutWeighted_( 0. ), afterCutWeighted_( 0. ),
  name_( cfg.getParameter<std::string>("name") ),
  min_ ( cfg.getParameter<std::vector<double> >("min") ),
  max_ ( cfg.getParameter<std::vector<double> >("max") )
{
  // choose operation mode (weak or strong filtering)
  switch(cfg.getParameter<unsigned>("mode")){
  case  0: mode_=kStrong; break;
  case  1: mode_=kWeak;   break;
  default: throw edm::Exception(edm::errors::Configuration, 
				"unknown mode for filter; choose between 0(strong) or 1(weak) \n");
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
bool ViewFilter::operator()(edm::Event& evt, const std::vector<edm::View<reco::Candidate> >& objs, const double& weight)
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
void ViewFilter::summarize()
{
  std::cout << ::std::setw( 20 ) << ::std::left  << name_ << " : "
	    << ::std::setw( 10 ) << ::std::right << afterCut_ << " (" 
	    << ::std::setw( 10 ) << ::std::right << afterCutWeighted_  << ") outof "
	    << ::std::setw( 10 ) << ::std::right << beforeCut_<< " (" 
	    << ::std::setw( 10 ) << ::std::right << beforeCutWeighted_ << ")" << std::endl;
}
