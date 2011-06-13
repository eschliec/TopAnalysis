#include "TopAnalysis/TopAnalyzer/interface/MuonId.h"
#include "DataFormats/MuonReco/interface/MuonEnergy.h"


/// constructor for FWLite analyzer
MuonId::MuonId() : 
  fwLite_(true ) { }

/// constructor for full FW analyzer
MuonId::MuonId(const edm::ParameterSet& cfg): 
  fwLite_(false)    
{
}

/// fill interface for full FW analyzer
void
MuonId::fill(const edm::Event& evt, const std::vector<pat::Muon>& muons, const double& weight=1.)
{
  if(muons.begin()==muons.end()) return;
  fill(muons, weight);
}

/// fill interface for FWLite analyzer
void
MuonId::fill(const std::vector<pat::Muon>& muons, const double& weight=1.)
{
  std::vector<pat::Muon>::const_iterator muon=muons.begin();
  if(muon!=muons.end()){
    reco::MuonEnergy muEnergy = muon->calEnergy();   
    muComp_ ->Fill( muon->caloCompatibility(), weight );
    muEm_   ->Fill( muEnergy.em, weight );
    muEmS9_ ->Fill( muEnergy.emS9, weight );  
    muHad_  ->Fill( muEnergy.had, weight );
    muHadS9_->Fill( muEnergy.hadS9,weight );  
    muHo_   ->Fill( muEnergy.ho, weight );
    muHoS9_ ->Fill( muEnergy.hoS9, weight );  
  }
}

/// book for FWLite
void 
MuonId::book()
{
  NameScheme mu("id");
  muComp_ = new TH1F(mu.name("muComp" ), mu.name("muComp" ), 50,  0.,  1.);
  muEm_   = new TH1F(mu.name("muEm" ),   mu.name("muEm" ),   50,  0.,  5.); 
  muEmS9_ = new TH1F(mu.name("muEmS9" ), mu.name("muEmS9" ), 50,  0.,  5.); 
  muHad_  = new TH1F(mu.name("muHad" ),  mu.name("muHad" ),  50,  0.,  5.); 
  muHadS9_= new TH1F(mu.name("muHadS9" ),mu.name("muHadS9" ),50,  0.,  5.);
  muHo_   = new TH1F(mu.name("muHo" ),   mu.name("muHo" ),   50,  0.,  5.); 
  muHoS9_ = new TH1F(mu.name("muHoS9" ), mu.name("muHoS9" ), 50,  0.,  5.);   
}

/// book for full FW
void 
MuonId::book(edm::Service<TFileService>& fs)
{
  NameScheme mu("id");
  muComp_ = fs->make<TH1F>(mu.name("muComp" ), mu.name("muComp" ), 50,  0.,  1.);
  muEm_   = fs->make<TH1F>(mu.name("muEm" ),   mu.name("muEm" ),   50,  0.,  5.); 
  muEmS9_ = fs->make<TH1F>(mu.name("muEmS9" ), mu.name("muEmS9" ), 50,  0.,  5.); 
  muHad_  = fs->make<TH1F>(mu.name("muHad" ),  mu.name("muHad" ),  50,  0.,  5.); 
  muHadS9_= fs->make<TH1F>(mu.name("muHadS9" ),mu.name("muHadS9" ),50,  0.,  5.);
  muHo_   = fs->make<TH1F>(mu.name("muHo" ),   mu.name("muHo" ),   50,  0.,  5.); 
  muHoS9_ = fs->make<TH1F>(mu.name("muHoS9" ), mu.name("muHoS9" ), 50,  0.,  5.);   
}

/// book for full FW with output stream
void 
MuonId::book(edm::Service<TFileService>& fs, ofstream& file)
{
  NameScheme mu("id");
  muComp_ = fs->make<TH1F>(mu.name(file, "muComp" ), mu.name("muComp" ), 50,  0.,  1.);
  muEm_   = fs->make<TH1F>(mu.name(file, "muEm" ),   mu.name("muEm" ),   50,  0.,  5.); 
  muEmS9_ = fs->make<TH1F>(mu.name(file, "muEmS9" ), mu.name("muEmS9" ), 50,  0.,  5.); 
  muHad_  = fs->make<TH1F>(mu.name(file, "muHad" ),  mu.name("muHad" ),  50,  0.,  5.); 
  muHadS9_= fs->make<TH1F>(mu.name(file, "muHadS9" ),mu.name("muHadS9" ),50,  0.,  5.);
  muHo_   = fs->make<TH1F>(mu.name(file, "muHo" ),   mu.name("muHo" ),   50,  0.,  5.); 
  muHoS9_ = fs->make<TH1F>(mu.name(file, "muHoS9" ), mu.name("muHoS9" ), 50,  0.,  5.);   
}

/// write to file and free allocated space for FWLite
void 
MuonId::write(TFile& file, const char* directory)
{
  /// save histograms to file
  file.cd( directory );
  muComp_ ->Write( );
  muEm_   ->Write( );
  muEmS9_ ->Write( );
  muHad_  ->Write( );
  muHadS9_->Write( );
  muHo_   ->Write( );
  muHoS9_ ->Write( );
}
