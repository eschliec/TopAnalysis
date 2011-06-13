#include "TMath.h"
#include "TopAnalysis/TopAnalyzer/interface/ElecId.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"

/// constructor for FWLite analyzer
ElecId::ElecId()  
{
}

/// constructor for full FW analyzer
ElecId::ElecId(const edm::ParameterSet& cfg)
{
}

/// fill interface for full FW analyzer
void
ElecId::fill(const edm::Event& evt, const std::vector<pat::Electron>& elecs, const double& weight=1.)
{

 if(elecs.begin()==elecs.end()) return;
 fill( elecs, weight);
}

/// fill interface for FWLite analyzer
void
ElecId::fill(const std::vector<pat::Electron>& elecs, const double& weight=1.)
{
  std::vector<pat::Electron>::const_iterator elec=elecs.begin();
  if(elec!=elecs.end()){
    eops_->Fill( elec->eSuperClusterOverP(),             weight );
    eopb_->Fill( elec->eSeedClusterOverPout(),           weight );
    nocl_->Fill( elec->numberOfClusters(),               weight );
    hoem_->Fill( TMath::Log10( elec->hadronicOverEm() ), weight );
    deta_->Fill( elec->deltaEtaSeedClusterTrackAtCalo(), weight );
    dphi_->Fill( elec->deltaPhiSeedClusterTrackAtCalo(), weight );
    
    //reco::GsfTrackRef track = elec->gsfTrack();
    //reco::Particle::Vector iMom = track->innerMomentum();
    //reco::Particle::Vector oMom = track->outerMomentum();
    
    //tdpt_ ->Fill( (iMom.Rho()-oMom.Rho())/iMom.Rho(), weight );
    //tdeta_->Fill( (iMom.Eta()-oMom.Eta())/iMom.Eta(), weight );
    //tdphi_->Fill( (iMom.Phi()-oMom.Phi())/iMom.Phi(), weight );
    //nohit_->Fill( track->recHitsSize(),               weight );
    //nvhit_->Fill( track->found(),                     weight );
    //chi2_ ->Fill( track->chi2()/track->ndof(),        weight );

    reco::SuperClusterRef clus = elec->superCluster();
    math::XYZPoint tPos(elec->TrackPositionAtCalo().X(), 
			elec->TrackPositionAtCalo().Y(),
			elec->TrackPositionAtCalo().Z());
    drtk_->Fill( (clus->position()-tPos).R(), weight );  
    
  }
}

/// book for FWLite
void 
ElecId::book()
{
  NameScheme id("id");
  eops_ = new TH1F(id.name("eops" ), id.name("eops" ), 50,   0.,   5.);
  eopb_ = new TH1F(id.name("eopb" ), id.name("eopb" ), 50,   0.,   5.);
  nocl_ = new TH1F(id.name("nocl" ), id.name("nocl" ), 20,   0.,  20.);
  deta_ = new TH1F(id.name("deta" ), id.name("deta" ), 40,-0.02, 0.02);
  dphi_ = new TH1F(id.name("dphi" ), id.name("dphi" ), 40,-0.02, 0.02);
  hoem_ = new TH1F(id.name("hoem" ), id.name("hoem" ), 50,  -5.,   0.);
  tdpt_ = new TH1F(id.name("tdpt" ), id.name("dpt"  ), 30, -1.5,  1.5);
  tdeta_= new TH1F(id.name("tdeta"), id.name("deta" ), 30, -1.5,  1.5);
  tdphi_= new TH1F(id.name("tdphi"), id.name("dphi" ), 30, -1.5,  1.5);
  nohit_= new TH1F(id.name("nohit"), id.name("nohit"), 25,   0.,  25.);
  nvhit_= new TH1F(id.name("nvhit"), id.name("nvhit"), 25,   0.,  25.);
  chi2_ = new TH1F(id.name("chi2" ), id.name("chi2" ), 50,   0.,  10.);
  drtk_ = new TH1F(id.name("drtrk"), id.name("drtrk"), 31, -0.1,   3.);
  s1os9_= new TH1F(id.name("s1os9"), id.name("s1os9"), 25,   0.,   1.);
}

/// book for full FW
void 
ElecId::book(edm::Service<TFileService>& fs)
{
  NameScheme id("id");
  eops_ = fs->make<TH1F>(id.name("eops" ), id.name("eops" ), 50,   0.,   5.);
  eopb_ = fs->make<TH1F>(id.name("eopb" ), id.name("eopb" ), 50,   0.,   5.);
  nocl_ = fs->make<TH1F>(id.name("nocl" ), id.name("nocl" ), 20,   0.,  20.);
  deta_ = fs->make<TH1F>(id.name("deta" ), id.name("deta" ), 40,-0.02, 0.02);
  dphi_ = fs->make<TH1F>(id.name("dphi" ), id.name("dphi" ), 40,-0.02, 0.02);
  hoem_ = fs->make<TH1F>(id.name("hoem" ), id.name("hoem" ), 50,  -5.,   0.);
  tdpt_ = fs->make<TH1F>(id.name("tdpt" ), id.name("dpt"  ), 30, -1.5,  1.5);
  tdeta_= fs->make<TH1F>(id.name("tdeta"), id.name("deta" ), 30, -1.5,  1.5);
  tdphi_= fs->make<TH1F>(id.name("tdphi"), id.name("dphi" ), 30, -1.5,  1.5);
  nohit_= fs->make<TH1F>(id.name("nohit"), id.name("nohit"), 25,   0.,  25.);
  nvhit_= fs->make<TH1F>(id.name("nvhit"), id.name("nvhit"), 25,   0.,  25.);
  chi2_ = fs->make<TH1F>(id.name("chi2" ), id.name("chi2" ), 50,   0.,  10.);
  drtk_ = fs->make<TH1F>(id.name("drtrk"), id.name("drtrk"), 31, -0.1,   3.);
  s1os9_= fs->make<TH1F>(id.name("s1os9"), id.name("s1os9"), 25,   0.,   1.);
}

/// book for full FW with output stream
void 
ElecId::book(edm::Service<TFileService>& fs, ofstream& file)
{
  NameScheme id("id");
  eops_ = fs->make<TH1F>(id.name(file, "eops" ), id.name("eops" ), 50,   0.,   5.);
  eopb_ = fs->make<TH1F>(id.name(file, "eopb" ), id.name("eopb" ), 50,   0.,   5.);
  nocl_ = fs->make<TH1F>(id.name(file, "nocl" ), id.name("nocl" ), 20,   0.,  20.);
  deta_ = fs->make<TH1F>(id.name(file, "deta" ), id.name("deta" ), 40,-0.02, 0.02);
  dphi_ = fs->make<TH1F>(id.name(file, "dphi" ), id.name("dphi" ), 40,-0.02, 0.02);
  hoem_ = fs->make<TH1F>(id.name(file, "hoem" ), id.name("hoem" ), 50,  -5.,   0.);
  tdpt_ = fs->make<TH1F>(id.name(file, "tdpt" ), id.name("dpt"  ), 30, -1.5,  1.5);
  tdeta_= fs->make<TH1F>(id.name(file, "tdeta"), id.name("deta" ), 30, -1.5,  1.5);
  tdphi_= fs->make<TH1F>(id.name(file, "tdphi"), id.name("dphi" ), 30, -1.5,  1.5);
  nohit_= fs->make<TH1F>(id.name(file, "nohit"), id.name("nohit"), 25,   0.,  25.);
  nvhit_= fs->make<TH1F>(id.name(file, "nvhit"), id.name("nvhit"), 25,   0.,  25.);
  chi2_ = fs->make<TH1F>(id.name(file, "chi2" ), id.name("chi2" ), 50,   0.,  10.);
  drtk_ = fs->make<TH1F>(id.name(file, "drtrk"), id.name("drtrk"), 31, -0.1,   3.);
  s1os9_= fs->make<TH1F>(id.name(file, "s1os9"), id.name("s1os9"), 25,   0.,   1.);
}

/// write to file and free allocated space for FWLite
void 
ElecId::write(const char* filename, const char* directory)
{
  /// save histograms to file
  TFile outFile( filename, "recreate" );
  outFile.mkdir( directory );
  outFile.cd( directory );

  /// basic kinematic
  eops_ ->Write( );
  eopb_ ->Write( );
  nocl_ ->Write( );
  deta_ ->Write( );
  dphi_ ->Write( );
  hoem_ ->Write( );
  tdpt_ ->Write( );
  tdeta_->Write( );
  tdphi_->Write( );
  nohit_->Write( );
  nvhit_->Write( );
  chi2_ ->Write( );
  drtk_ ->Write( );
  s1os9_->Write( );

  outFile.Close();

  // free allocated space
  delete eops_; 
  delete eopb_; 
  delete nocl_; 
  delete deta_; 
  delete dphi_;
  delete hoem_; 
  delete tdpt_; 
  delete tdeta_; 
  delete tdphi_; 
  delete nohit_; 
  delete nvhit_; 
  delete chi2_;
  delete drtk_; 
  delete s1os9_; 
}
