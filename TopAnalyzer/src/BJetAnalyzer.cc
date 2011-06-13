#include "TopAnalysis/TopAnalyzer/interface/BJetAnalyzer.h"

BJetAnalyzer::BJetAnalyzer(const edm::ParameterSet& cfg):
  jets_       ( cfg.getParameter<edm::InputTag>( "jets"            ) ),
  hist_       ( cfg.getParameter<std::string>  ( "hist"            ) ), 
  onlyLeading_( cfg.getParameter<int>          ( "only_leading"    ) ),   
  doKin_      ( cfg.getParameter<bool>         ( "show_kinematics" ) ),
  doES_       ( cfg.getParameter<bool>         ( "show_energyscale") ),
  doBtag_     ( cfg.getParameter<bool>         ( "show_btagging"   ) )
{  
  if( doKin_ ){ 
    kin_uds_ = new BJetKin();
    kin_c_   = new BJetKin();
    kin_b_   = new BJetKin();  
    kin_g_   = new BJetKin();  
  }
  if( doBtag_){ 
    btg_uds_ = new BJetTag(cfg);
    btg_c_   = new BJetTag(cfg);
    btg_b_   = new BJetTag(cfg);
    btg_g_   = new BJetTag(cfg);
  }
}

BJetAnalyzer::~BJetAnalyzer()
{
  if( doKin_){  
    delete kin_uds_;
    delete kin_c_;  
    delete kin_b_;  
    delete kin_g_;  
  }
  if( doBtag_){
    delete btg_uds_;
    delete btg_c_;
    delete btg_b_;
    delete btg_g_;
  }
}

void
BJetAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{    
  edm::Handle<PatJetCollection> jets; 
  evt.getByLabel(jets_, jets); 
                    
  int i=1;    
  for(PatJetCollection::const_iterator jet = jets->begin(); jet!= jets->end(); ++jet) { 
    
    // only take the first muons from the collection
    if(i>onlyLeading_) break;
    ++i; 
    
    // distinguish between b jets and other flavors
    const reco::Particle* parton = jet->genParton();
    if(!parton){
      nFlav_->SetBinContent(1, nFlav_->GetBinContent(1)+1);
      continue;
    }  
      
    int flavor = abs(parton->pdgId()); 
    
    switch(flavor){
      case  1 : nFlav_->SetBinContent(3, nFlav_->GetBinContent(3)+1); break;
      case  2 : nFlav_->SetBinContent(3, nFlav_->GetBinContent(3)+1); break;
      case  3 : nFlav_->SetBinContent(3, nFlav_->GetBinContent(3)+1); break;
      case  4 : nFlav_->SetBinContent(4, nFlav_->GetBinContent(4)+1); break;
      case  5 : nFlav_->SetBinContent(5, nFlav_->GetBinContent(5)+1); break;
      case 21 : nFlav_->SetBinContent(2, nFlav_->GetBinContent(2)+1); break;
      default : cout << "Jet matched to parton " << flavor << endl;
    }
      	       
    if( doKin_ ){ 
      switch(flavor){
        case  1 : kin_uds_ ->fill( *jet ); break;
	case  2 : kin_uds_ ->fill( *jet ); break;
	case  3 : kin_uds_ ->fill( *jet ); break;
	case  4 : kin_c_   ->fill( *jet ); break;
        case  5 : kin_b_   ->fill( *jet ); break;
        case 21 : kin_g_   ->fill( *jet ); break;
	default : cout << "Jet matched to parton " << flavor << endl;
      }
    }
   
    if( doBtag_){ 
      switch(flavor){
        case  1 : btg_uds_ ->fill( *jet ); break;
	case  2 : btg_uds_ ->fill( *jet ); break;
	case  3 : btg_uds_ ->fill( *jet ); break;
	case  4 : btg_c_   ->fill( *jet ); break;
        case  5 : btg_b_   ->fill( *jet ); break;
        case 21 : btg_g_   ->fill( *jet ); break;
	default : cout << "Jet matched to parton " << flavor << endl;
      }
    }
  }                 
}

void 
BJetAnalyzer::beginJob(const edm::EventSetup&)
{
  if( hist_.empty() )
    return;

  edm::Service<TFileService> fs;
  if( !fs ){
    throw edm::Exception( edm::errors::Configuration,
                          "TFile Service is not registered in cfg file" );
  }

  ofstream hist(hist_.c_str(), std::ios::out);
  NameScheme n("N");

  if( doKin_ ){ 
    kin_uds_->book(hist, "kin_uds");
    kin_c_  ->book(hist, "kin_c"  );    
    kin_b_  ->book(hist, "kin_b"  );    
    kin_g_  ->book(hist, "kin_g"  );      
  }
  if( doBtag_){ 
    btg_uds_->book(hist, "btg_uds");
    btg_c_  ->book(hist, "btg_c"  );    
    btg_b_  ->book(hist, "btg_b"  );    
    btg_g_  ->book(hist, "btg_g"  );    
  }
  nFlav_ = fs->make<TH1F>(n.name( hist, "nFlav" ), n.name("nFlav" ), 5, 0, 5);

}


void 
BJetAnalyzer::endJob() 
{
  if( doBtag_){ 
    btg_uds_->evaluate();
    btg_c_  ->evaluate();    
    btg_b_  ->evaluate();    
    btg_g_  ->evaluate();    
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( BJetAnalyzer );
