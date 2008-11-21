#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "TopAnalysis/TopUtils/interface/NameScheme.h"
#include "TopAnalysis/TopUtils/interface/RootHistograms.h"
#include "TopAnalysis/TopAnalyzer/interface/BJetKin.h"


BJetKin::BJetKin()
{
}

BJetKin::~BJetKin()
{
}

void
BJetKin::fill(const pat::Jet& jet, double weight)
{
  e_       ->Fill( jet.energy()          , weight );
  et_      ->Fill( jet.et()              , weight );  
  p_       ->Fill( jet.p()               , weight );
  pt_      ->Fill( jet.pt()              , weight );      
  eta_     ->Fill( jet.eta()             , weight );
  theta_   ->Fill( jet.theta()*180/3.1416, weight );  
  rapidity_->Fill( jet.rapidity()        , weight );  
  phi_     ->Fill( jet.phi()*180/3.1416  , weight );
  mass_    ->Fill( jet.mass()            , weight );
  mt_      ->Fill( jet.mt()              , weight );
  n60_     ->Fill( jet.n60()             , weight );
  n90_     ->Fill( jet.n90()             , weight );
  nConst_  ->Fill( jet.nConstituents()   , weight );
}

void 
BJetKin::book(ofstream& file, char* name)
{
  edm::Service<TFileService> fs;

  NameScheme e(name);
  e_       = fs->make<TH1F>(e.name( file, "e"        ), e.name("e"        ), 60,    0. , 300.  );
  et_      = fs->make<TH1F>(e.name( file, "et"       ), e.name("et"       ), 40,    0. , 200.  );
  p_       = fs->make<TH1F>(e.name( file, "p"        ), e.name("p"        ), 60,    0. , 300.  );
  pt_      = fs->make<TH1F>(e.name( file, "pt"       ), e.name("pt"       ), 40,    0. , 200.  );  
  eta_     = fs->make<TH1F>(e.name( file, "eta"      ), e.name("eta"      ), 50,   -2.5,   2.5 );
  theta_   = fs->make<TH1F>(e.name( file, "theta"    ), e.name("theta"    ), 36,    0. , 180.  );  
  rapidity_= fs->make<TH1F>(e.name( file, "rapidity" ), e.name("rapidity" ), 50,   -2.5,   2.5 );  
  phi_     = fs->make<TH1F>(e.name( file, "phi"      ), e.name("phi"      ), 36, -180. , 180.  );
  mass_    = fs->make<TH1F>(e.name( file, "mass"     ), e.name("mass"     ), 25,    0. ,  25.  );  
  mt_      = fs->make<TH1F>(e.name( file, "mt"       ), e.name("mt"       ), 50,    0. , 100.  );  
  n60_     = fs->make<TH1F>(e.name( file, "n60"      ), e.name("n60"      ), 10,    0.5,  10.5 ); 
  n90_     = fs->make<TH1F>(e.name( file, "n90"      ), e.name("n90"      ), 15,    0.5,  15.5 );     
  nConst_  = fs->make<TH1F>(e.name( file, "nConst"   ), e.name("nConst"   ), 25,    0.5,  25.5 );  
}
