#include "TopAnalysis/TopFilter/plugins/FullLepObjectCountFilter.h"

FullLepObjectCountFilter::FullLepObjectCountFilter(const edm::ParameterSet& cfg):
  wgt_        (cfg.getParameter<edm::InputTag>( "weight"  )),
  name_       (cfg.getParameter<std::string>  ( "name"    )),  
  objects_    (cfg.getParameter<edm::InputTag>( "objects" )),
  size_       (cfg.getParameter<unsigned int> ( "n"       )),      
  beforeCuts_( 0 ), afterCuts_( 0 ), beforeCutsWeighted_( 0. ), afterCutsWeighted_( 0. )
{
}

bool FullLepObjectCountFilter::filter(edm::Event& evt, const edm::EventSetup& setup)
{
  edm::Handle<double> wgt;
  evt.getByLabel(wgt_, wgt);
  double weight = *wgt;
  
  edm::Handle<edm::View<reco::Candidate> > objects; 
  evt.getByLabel(objects_, objects); 
    
  ++beforeCuts_;
  beforeCutsWeighted_ += weight;    
       
  if(objects->size()<size_) return false;
    
  ++afterCuts_;
  afterCutsWeighted_ += weight;

  return true;
}

void FullLepObjectCountFilter::beginJob(const edm::EventSetup& setup)
{   
  edm::LogVerbatim log("topFilter");  
  log << ::std::setw( 20 ) << ::std::left;  
  log << name_ << ": " << " n >= " << ::std::setw( 8 ) << ::std::right  <<  size_;   
}

void FullLepObjectCountFilter::endJob()
{
  edm::LogVerbatim log("topFilter");

  if(beforeCuts_ != beforeCutsWeighted_) {
    log << std::setw( 20 ) << std::left  << name_ << " : "
	<< std::setw( 10 ) << std::right << afterCuts_ << " (" 
	<< std::setw( 10 ) << std::right << afterCutsWeighted_  << ") out of "
	<< std::setw( 10 ) << std::right << beforeCuts_<< " (" 
	<< std::setw( 10 ) << std::right << beforeCutsWeighted_ << ")";	
  }
  else{
    log << std::setw( 20 ) << std::left  << name_ << " : "
	<< std::setw( 10 ) << std::right << afterCuts_ << "  out of "
	<< std::setw( 10 ) << std::right << beforeCuts_;    	
  }
}
