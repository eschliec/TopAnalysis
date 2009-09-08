#include "TLorentzVector.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TopAnalysis/TopAnalyzer/interface/DimuonAnalyzer.h"

DimuonAnalyzer::DimuonAnalyzer(const edm::ParameterSet& cfg):
  muons_      ( cfg.getParameter<edm::InputTag>( "muons"          ) ),
  ref_        ( cfg.getParameter<edm::InputTag>( "ref"            ) ),  
  matchMuons_ ( cfg.getParameter<bool>         ( "doMatch"        ) ),  
  useEvtWgt_  ( cfg.getParameter<bool>         ( "useEventWeight" ) ),  
  hist_       ( cfg.getParameter<std::string>  ( "hist"           ) )  
{
}

DimuonAnalyzer::~DimuonAnalyzer()
{
}

void
DimuonAnalyzer::beginJob(const edm::EventSetup&)
{
  if( hist_.empty() )
    return;

  edm::Service<TFileService> fs;
  if( !fs ){
    throw edm::Exception( edm::errors::Configuration,
                          "TFile Service is not registered in cfg file" );
  }

  ofstream hist(hist_.c_str(), std::ios::out);
  // define logarithmic bins for a histogram with 50 bins going from 10^0 to 10^2.4 =~ 250
  const int nbins = 50;
  double logmin =  0.0;
  double logmax =  2.4;
  double bins[nbins+1];
  for (int i = 0; i <= nbins; i++) {
    double log = logmin + (logmax-logmin)*i/nbins;
    bins[i] = std::pow(10.0, log);
  }

  NameScheme e("dimuon");
  dimassRC_    = fs->make<TH1F>(e.name(hist, "dimassRC"    ), e.name("dimassRC"    ), nbins, bins);  
  dimassWC_    = fs->make<TH1F>(e.name(hist, "dimassWC"    ), e.name("dimassWC"    ), nbins, bins);
         
  isoDimassCorrelation_ = fs->make<TH2F>(e.name(hist, "isoDimassCorrelation" ), e.name("isoDimassCorrelation" ), 10,   0.,   200., 10, 0., 1.);
  
  absCount_    = fs->make<TH1F>(e.name( hist, "absCount"    ), e.name("absCount"    ),100,  0.,  50.0);        
  relCount_    = fs->make<TH1F>(e.name( hist, "relCount"    ), e.name("relCount"    ),100,  0.,   5.0); 
  combCount_   = fs->make<TH1F>(e.name( hist, "combCount"   ), e.name("combCount"   ),100,  0.0,  1.0); 
  diCombCount_ = fs->make<TH1F>(e.name( hist, "diCombCount" ), e.name("diCombCount" ),100,  0.0,  2.0);
}

void
DimuonAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup&)
{
  edm::Handle<PatMuonCollection> muons; 
  evt.getByLabel(muons_, muons); 
                 
  double weight = 1.;
  if(useEvtWgt_){ 
    edm::Handle<double> weightHandle; 
    evt.getByLabel("eventWeight", weightHandle); 
    weight = *weightHandle;
  }

  if(muons->size()<2) return;
        
  PatMuonCollection::const_reference mu1 = muons->at(0);
  PatMuonCollection::const_reference mu2 = muons->at(1);   
  
  if(matchMuons_){       
    if(!mu1.genLepton() || !mu2.genLepton()){ 
      return;
    }
    edm::Handle<reco::GenParticleCollection> genParticles;
    evt.getByLabel(ref_, genParticles);
    
    const reco::GenParticle* genMuMinus = 0;
    const reco::GenParticle* genMuPlus = 0;
    
    for(reco::GenParticleCollection::const_iterator part = genParticles->begin(); part!=genParticles->end(); ++part ){
      if(part->pdgId()== 13) genMuMinus = &*part;
      if(part->pdgId()==-13) genMuPlus  = &*part;
    }
    
    if(!genMuMinus || !genMuPlus) return;
    
    const reco::GenParticle* genMu1 = mu1.genLepton();
    const reco::GenParticle* genMu2 = mu2.genLepton();    
    
    if( !(matches(genMu1, genMuPlus) && matches(genMu2, genMuMinus)) && !(matches(genMu2, genMuPlus) && matches(genMu1, genMuMinus)) ){ 
      return;
    }
  } 
      
  TLorentzVector diLepLVector = TLorentzVector(mu1.px()+mu2.px(), mu1.py()+mu2.py(), 
                                               mu1.pz()+mu2.pz(), mu1.energy()+mu2.energy());
   
  double dilepMass = (diLepLVector).M();
       
  
  if(mu1.charge()*mu2.charge()<0.)  
    dimassRC_  ->Fill( dilepMass, weight);    
  
  if(mu1.charge()*mu2.charge()>0.) 
    dimassWC_  ->Fill( dilepMass, weight);    
     
    
  double absTrackIso1 = mu1.trackIso();
  double absTrackIso2 = mu2.trackIso();
  
  double absCaloIso1 = mu1.caloIso();
  double absCaloIso2 = mu2.caloIso();  
  
  double relTrackIso1 = mu1.trackIso()/mu1.pt();
  double relTrackIso2 = mu2.trackIso()/mu2.pt();
  
  double relCaloIso1 = mu1.caloIso()/mu1.pt();
  double relCaloIso2 = mu2.caloIso()/mu2.pt();  
  
  double combIso1 = mu1.pt()/(mu1.pt()+mu1.trackIso()+mu1.caloIso());
  double combIso2 = mu2.pt()/(mu2.pt()+mu2.trackIso()+mu2.caloIso());     
  
  double diCombIso = sqrt(combIso1*combIso1+combIso2*combIso2);
    
  isoDimassCorrelation_->Fill( dilepMass, combIso1, weight);  
  isoDimassCorrelation_->Fill( dilepMass, combIso2, weight);
      
  for(int i=1;i<=100;++i){
    if(diCombIso>0.02*i){
      diCombCount_->SetBinContent(i,diCombCount_->GetBinContent(i)+1);
    }
    if(combIso1>0.01*i && combIso2>0.01*i){
      combCount_->SetBinContent(i,combCount_->GetBinContent(i)+1);
    }
    if(relTrackIso1<(0.05*i) && relCaloIso1<(0.05*i) && relTrackIso2<(0.05*i) && relCaloIso2<(0.05*i)){
      relCount_->SetBinContent(i,relCount_->GetBinContent(i)+1);
    }
    if(absTrackIso1<(0.5*i) && absCaloIso1<(0.5*i) && absTrackIso2<(0.5*i) && absCaloIso2<(0.5*i)){
      absCount_->SetBinContent(i,absCount_->GetBinContent(i)+1);
    }      
  }     
}

void
DimuonAnalyzer::endJob()
{
}

bool
DimuonAnalyzer::matches(const reco::GenParticle* p1, const reco::GenParticle* p2)
{
  if( deltaR(p1->phi(), p1->eta(), p2->phi(), p2->eta()) < 0.01 ) return true;
  return false;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( DimuonAnalyzer );
