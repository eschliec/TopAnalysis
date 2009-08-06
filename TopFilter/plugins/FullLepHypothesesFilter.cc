#include "TopAnalysis/TopFilter/plugins/FullLepHypothesesFilter.h"
#include "TLorentzVector.h"

FullLepHypothesesFilter::FullLepHypothesesFilter(const edm::ParameterSet& cfg):
  hypoKey_    (cfg.getParameter<edm::InputTag>(       "hypoKey"        )),
  FullLepEvt_ (cfg.getParameter<edm::InputTag>(       "FullLepEvent"   )),
  wgt_        (cfg.getParameter<edm::InputTag>(       "weight"         )),
  jets_       (cfg.getParameter<edm::InputTag>(       "jets"           )),  
  qcdCut_     (cfg.getParameter<double>(              "qcdCut"         )),
  zCut_       (cfg.getParameter<std::vector<double> >("zCut"           )),        
  bAlgo_      (cfg.getParameter<std::string >(        "bAlgorithm"     )),  
  bDisc_      (cfg.getParameter<std::vector<double> >("bDiscriminator" )),     
  beforeCuts_( 0 ), validCounter_( 0 ), betweenCuts_( 0 ), afterCuts_( 0 ), 
  beforeCutsWeighted_( 0. ), validCounterWeighted_( 0. ), betweenCutsWeighted_( 0. ), afterCutsWeighted_( 0. )
{
}

bool FullLepHypothesesFilter::filter(edm::Event& evt, const edm::EventSetup& setup)
{
  edm::Handle<TtFullLeptonicEvent> FullLepEvt;
  evt.getByLabel(FullLepEvt_, FullLepEvt);

  edm::Handle<int> hypoKeyHandle;
  evt.getByLabel(hypoKey_, hypoKeyHandle);
  TtEvent::HypoClassKey& hypoKey = (TtEvent::HypoClassKey&) *hypoKeyHandle;

  edm::Handle<double> wgt;
  evt.getByLabel(wgt_, wgt);
  double weight = *wgt;
  
  edm::Handle<std::vector<pat::Jet> > jets; 
  evt.getByLabel(jets_, jets);   

  ++beforeCuts_;
  beforeCutsWeighted_ += weight;

  // -----------------------
  // check if hypothesis is valid in this event
  // and if it satisfies some quality criteria
  // -----------------------
         
  if( !FullLepEvt->isHypoAvailable(hypoKey) ){
    edm::LogInfo ( "NonValidHyp" ) << "Hypothesis not available for this event";
    return false;
  }
  if( !FullLepEvt->isHypoValid(hypoKey) ){
    edm::LogInfo ( "NonValidHyp" ) << "Hypothesis not valid for this event";
    return false; 
  }
 
  ++validCounter_;
  validCounterWeighted_ += weight; 
  
  // get objects needed for this filter
  const reco::Candidate* LepBar = FullLepEvt->leptonBar(hypoKey);    
  const reco::Candidate* Lep    = FullLepEvt->lepton(hypoKey);
  
  int bJetIdx    = FullLepEvt->jetLepComb(hypoKey)[0];
  int bBarJetIdx = FullLepEvt->jetLepComb(hypoKey)[1];  

  const pat::Jet B    = jets->at(bJetIdx);   
  const pat::Jet BBar = jets->at(bBarJetIdx);   

  // reconstruct invariant dilepton mass and make cuts on it
  TLorentzVector diLepLVector = TLorentzVector(LepBar->px()+Lep->px(), LepBar->py()+Lep->py(), 
                                               LepBar->pz()+Lep->pz(), LepBar->energy()+Lep->energy());   
  double dilepMass = (diLepLVector).M();  
  // cut to avoid qcd resonances like J/Psi or cascadic decays
  if(dilepMass < qcdCut_) return false;
  // cut to surpress Z
  if(dilepMass > zCut_[0] && dilepMass < zCut_[1]) return false;
 
  ++betweenCuts_;
  betweenCutsWeighted_ += weight; 
 
  // make cut(s) on b-tag discriminator
  if(bDisc_.size()==1){
    if(!(B.bDiscriminator(bAlgo_)>=bDisc_[0] || BBar.bDiscriminator(bAlgo_)>=bDisc_[0])) return false;
  }
  else if(bDisc_.size()==2){
    if((B.bDiscriminator(bAlgo_)<bDisc_[0] || BBar.bDiscriminator(bAlgo_)<bDisc_[1]) 
      && (B.bDiscriminator(bAlgo_)<bDisc_[1] || BBar.bDiscriminator(bAlgo_)<bDisc_[0])) return false;
  }   
  
  ++afterCuts_;
  afterCutsWeighted_ += weight;

  return true;
}

void FullLepHypothesesFilter::beginJob(const edm::EventSetup& setup)
{ 
  // sanity check for zCut
  edm::LogError err("topFilter");
  if(zCut_.size()!=2)   err << "zCut has wrong size. Size has to be 2!\n";
  if(zCut_[0]>zCut_[1]) err << "Lower zCut is higher than upper zCut. All events will be skipped!\n";
  
  edm::LogVerbatim log("topFilter");
  log << "ValidityFullLepHyp\n";
  
  log << ::std::setw( 20 ) << ::std::left;  
  log << "DiLepMassFullLepHyp" << ": "   
      << "min = " << ::std::setw( 8 ) << ::std::right  <<  zCut_[0] << "  "
      << "max = " << ::std::setw( 8 ) << ::std::right  <<  zCut_[1] << "\n";
  if(qcdCut_!=0.){ 
    log << ::std::setw( 20 ) << ::std::left << " " << ": ";
    log << "qcd = " << ::std::setw( 8 ) << ::std::right  <<  qcdCut_ << "\n";
  }
  
  for(unsigned int idx=0; idx<bDisc_.size(); ++idx){
    log << ::std::setw( 20 ) << ::std::left;
    if(idx==0) log << "bTagFullLepHyp"; else log << " ";
    log << ": " << "tag >= " <<  bDisc_[idx] << "\n";
  }    
}

void FullLepHypothesesFilter::endJob()
{
  edm::LogVerbatim log("topFilter");

  if(beforeCuts_ != beforeCutsWeighted_) {
    log << std::setw( 20 ) << std::left  << "ValidityFullLepHyp" << " : "
	<< std::setw( 10 ) << std::right << validCounter_ << " (" 
	<< std::setw( 10 ) << std::right << validCounterWeighted_  << ") out of "
	<< std::setw( 10 ) << std::right << beforeCuts_<< " (" 
	<< std::setw( 10 ) << std::right << beforeCutsWeighted_ << ")" << "\n";
	
    log << std::setw( 20 ) << std::left  << "DiLepMassFullLepHyp" << " : "
	<< std::setw( 10 ) << std::right << betweenCuts_ << " (" 
	<< std::setw( 10 ) << std::right << betweenCutsWeighted_  << ") out of "
	<< std::setw( 10 ) << std::right << validCounter_<< " (" 
	<< std::setw( 10 ) << std::right << validCounterWeighted_ << ")" << "\n";
	
    log << std::setw( 20 ) << std::left  << "bTagFullLepHyp" << " : "
	<< std::setw( 10 ) << std::right << afterCuts_ << " (" 
	<< std::setw( 10 ) << std::right << afterCutsWeighted_  << ") out of "
	<< std::setw( 10 ) << std::right << betweenCuts_<< " (" 
	<< std::setw( 10 ) << std::right << betweenCutsWeighted_ << ")" << "\n";		
  }
  else{
    log << std::setw( 20 ) << std::left  << "ValidityFullLepHyp" << " : "
	<< std::setw( 10 ) << std::right << validCounter_ << "  out of "
	<< std::setw( 10 ) << std::right << beforeCuts_ << "\n";  
	
    log << std::setw( 20 ) << std::left  << "DiLepMassFullLepHyp" << " : "
	<< std::setw( 10 ) << std::right << betweenCuts_ << "  out of "
	<< std::setw( 10 ) << std::right << validCounter_ << "\n";
	
    log << std::setw( 20 ) << std::left  << "bTagFullLepHyp" << " : "
	<< std::setw( 10 ) << std::right << afterCuts_ << "  out of "
	<< std::setw( 10 ) << std::right << betweenCuts_ << "\n"; 	  	
  }
}
