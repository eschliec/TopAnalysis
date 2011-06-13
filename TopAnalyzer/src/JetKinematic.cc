#include "TopAnalysis/TopAnalyzer/interface/JetKinematic.h"


/// constructor for FWLite analyzer
JetKinematic::JetKinematic(int nJets):
  nJets_(nJets)
{
}

/// constructor for full FW analyzer
JetKinematic::JetKinematic(const edm::ParameterSet& cfg):
  nJets_( cfg.getParameter<int>( "nJets" ) )
{
}

/// fill interface for full FW analyzer
void
JetKinematic::fill(const edm::Event& evt, const std::vector<pat::Jet>& jets, const double& weight)
{
  fill(jets, weight);
}

/// fill interface for FWLite analyzer
void
JetKinematic::fill(const std::vector<pat::Jet>& jets, const double& weight)
{
  double sumEt=0;
  unsigned int idx=0;
  int n10=0, n20=0, n30=0, n40=0;
  for(std::vector<pat::Jet>::const_iterator jet = jets.begin(); 
      jet!=jets.end(); ++jet, ++idx) {
    allEn_ ->Fill( jet->energy(), weight );
    allEt_ ->Fill( jet->et(), weight );
    allEta_->Fill( jet->eta(), weight );
    allPhi_->Fill( jet->phi(), weight );
    
    if( idx< en_.size() ) en_ [idx]->Fill( jet->energy(), weight );
    if( idx< et_.size() ) et_ [idx]->Fill( jet->et(), weight );
    if( idx<eta_.size() ) eta_[idx]->Fill( jet->eta(), weight );
    if( idx<phi_.size() ) phi_[idx]->Fill( jet->phi(), weight );

    //---------------------------------------------
    // check scalar sum pt
    //---------------------------------------------    
    sumEt+=jet->et();
    if(idx == 3) scalSum4_->Fill( sumEt, weight );
    if(idx == 5) scalSum6_->Fill( sumEt, weight );

    //---------------------------------------------
    // check jet multiplicities
    //---------------------------------------------
    if( jet->et() > 10 ) ++n10;
    if( jet->et() > 20 ) ++n20;
    if( jet->et() > 30 ) ++n30;
    if( jet->et() > 40 ) ++n40;
  }
  mult_->Fill( jets.size(), weight );
  n10_ ->Fill( n10, weight );
  n20_ ->Fill( n20, weight );
  n30_ ->Fill( n30, weight );
  n40_ ->Fill( n40, weight );
}

/// book for FWLite
void 
JetKinematic::book()
{
  NameScheme kin("kin");
  allEn_ = new TH1F(kin.name("allE"  ), kin.name("allE"  ), 50, 0., 500.);
  allEt_ = new TH1F(kin.name("allEt" ), kin.name("allEt" ), 30, 0., 300.);
  allEta_= new TH1F(kin.name("allEta"), kin.name("allEta"), 34, -3.4,  3.4);
  allPhi_= new TH1F(kin.name("allPhi"), kin.name("allPhi"), 35, -3.5,  3.5);
  scalSum4_ = new TH1F(kin.name("scalSum4"), kin.name("scalSum4"), 40, 0., 800.);
  scalSum6_ = new TH1F(kin.name("scalSum6"), kin.name("scalSum6"), 40, 0., 800.);

  for(int idx=0; idx<nJets_; ++idx){
    en_.push_back ( new TH1F(kin.name("e",  idx), kin.name("en",  idx), 50,   0., 500.) );
    et_.push_back ( new TH1F(kin.name("et", idx), kin.name("et",  idx), 30,   0., 300.) );
    eta_.push_back( new TH1F(kin.name("eta",idx), kin.name("eta", idx), 34, -3.4,  3.4) );
    phi_.push_back( new TH1F(kin.name("phi",idx), kin.name("phi", idx), 35, -3.5,  3.5) );
  }

  mult_= new TH1F(kin.name("n"  ), kin.name("n"  ), 15, 0., 15.);
  n10_ = new TH1F(kin.name("n10"), kin.name("n10"), 15, 0., 15.);
  n20_ = new TH1F(kin.name("n20"), kin.name("n20"), 15, 0., 15.);
  n30_ = new TH1F(kin.name("n30"), kin.name("n30"), 15, 0., 15.);
  n40_ = new TH1F(kin.name("n40"), kin.name("n40"), 15, 0., 15.);
}

/// book for full FW
void 
JetKinematic::book(edm::Service<TFileService>& fs)
{
  NameScheme kin("kin");
  allEn_ = fs->make<TH1F>(kin.name("allE"  ), kin.name("allE"  ), 50, 0., 500.);
  allEt_ = fs->make<TH1F>(kin.name("allEt" ), kin.name("allEt" ), 30, 0., 300.);
  allEta_= fs->make<TH1F>(kin.name("allEta"), kin.name("allEta"), 34, -3.4,  3.4);
  allPhi_= fs->make<TH1F>(kin.name("allPhi"), kin.name("allPhi"), 35, -3.5,  3.5);
  scalSum4_ = fs->make<TH1F>(kin.name("scalSum4"), kin.name("scalSum4"), 40, 0., 800.);
  scalSum6_ = fs->make<TH1F>(kin.name("scalSum6"), kin.name("scalSum6"), 40, 0., 800.);

  for(int idx=0; idx<nJets_; ++idx){
    en_.push_back ( fs->make<TH1F>(kin.name("e",  idx), kin.name("en",  idx), 50,   0., 500.) );
    et_.push_back ( fs->make<TH1F>(kin.name("et", idx), kin.name("et",  idx), 30,   0., 300.) );
    eta_.push_back( fs->make<TH1F>(kin.name("eta",idx), kin.name("eta", idx), 34, -3.4,  3.4) );
    phi_.push_back( fs->make<TH1F>(kin.name("phi",idx), kin.name("phi", idx), 35, -3.5,  3.5) );
  }

  mult_= fs->make<TH1F>(kin.name("n"  ), "N_{jets}"                 , 15, 0., 15.);
  n10_ = fs->make<TH1F>(kin.name("n10"), "N_{jets} (p_{T} > 10 GeV)", 15, 0., 15.);
  n20_ = fs->make<TH1F>(kin.name("n20"), "N_{jets} (p_{T} > 20 GeV)", 15, 0., 15.);
  n30_ = fs->make<TH1F>(kin.name("n30"), "N_{jets} (p_{T} > 30 GeV)", 15, 0., 15.);
  n40_ = fs->make<TH1F>(kin.name("n40"), "N_{jets} (p_{T} > 40 GeV)", 15, 0., 15.);
}

/// book for full FW with output stream
void 
JetKinematic::book(edm::Service<TFileService>& fs, ofstream& file)
{
  NameScheme kin("kin");
  allEn_ = fs->make<TH1F>(kin.name(file, "allE"  ), kin.name("allE"  ), 50, 0., 500.);
  allEt_ = fs->make<TH1F>(kin.name(file, "allEt" ), kin.name("allEt" ), 30, 0., 300.);
  allEta_= fs->make<TH1F>(kin.name(file, "allEta"), kin.name("allEta"), 34, -3.4,  3.4);
  allPhi_= fs->make<TH1F>(kin.name(file, "allPhi"), kin.name("allPhi"), 35, -3.5,  3.5);
  scalSum4_ = fs->make<TH1F>(kin.name(file, "scalSum4"), kin.name("scalSum4"), 40, 0., 800.);
  scalSum6_ = fs->make<TH1F>(kin.name(file, "scalSum6"), kin.name("scalSum6"), 40, 0., 800.);

  for(int idx=0; idx<nJets_; ++idx){
    en_.push_back ( fs->make<TH1F>(kin.name(file, "e",  idx), kin.name("en",  idx), 50,   0., 500.) );
    et_.push_back ( fs->make<TH1F>(kin.name(file, "et", idx), kin.name("et",  idx), 30,   0., 300.) );
    eta_.push_back( fs->make<TH1F>(kin.name(file, "eta",idx), kin.name("eta", idx), 34, -3.4,  3.4) );
    phi_.push_back( fs->make<TH1F>(kin.name(file, "phi",idx), kin.name("phi", idx), 35, -3.5,  3.5) );
  }

  mult_= fs->make<TH1F>(kin.name(file, "n"  ), "N_{jets}"                 , 15, 0., 15.);
  n10_ = fs->make<TH1F>(kin.name(file, "n10"), "N_{jets} (p_{T} > 10 GeV)", 15, 0., 15.);
  n20_ = fs->make<TH1F>(kin.name(file, "n20"), "N_{jets} (p_{T} > 20 GeV)", 15, 0., 15.);
  n30_ = fs->make<TH1F>(kin.name(file, "n30"), "N_{jets} (p_{T} > 30 GeV)", 15, 0., 15.);
  n40_ = fs->make<TH1F>(kin.name(file, "n40"), "N_{jets} (p_{T} > 40 GeV)", 15, 0., 15.);
}

/// write to file and free allocated space for FWLite
void 
JetKinematic::write(const char* filename, const char* directory)
{
  /// save histograms to file
  TFile outFile( filename, "recreate" );
  outFile.mkdir( directory );
  outFile.cd( directory );

  /// basic kinematic
  mult_->Write( );
  n10_ ->Write( );
  n20_ ->Write( );
  n30_ ->Write( );
  n40_ ->Write( );

  outFile.Close();

  // free allocated space
  delete mult_;
  delete n10_;
  delete n20_;
  delete n30_;
  delete n40_;
}
