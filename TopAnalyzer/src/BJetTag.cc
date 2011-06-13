#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "TopAnalysis/TopUtils/interface/NameScheme.h"
#include "TopAnalysis/TopUtils/interface/RootHistograms.h"
#include "TopAnalysis/TopAnalyzer/interface/BJetTag.h"

BJetTag::BJetTag(const edm::ParameterSet& cfg):
  trkCountHighEffCut_( cfg.getParameter<double>( "trk_count_high_eff_cut" ) ),
  trkCountHighPurCut_( cfg.getParameter<double>( "trk_count_high_pur_cut" ) ),
  jetProbCut_        ( cfg.getParameter<double>( "jet_prob_cut"          ) ),
  jetBProbCut_       ( cfg.getParameter<double>( "jet_b_prob_cut"        ) ),
  impactParMvaCut_   ( cfg.getParameter<double>( "impact_par_mva_cut"     ) ),
  simpleSvCut_       ( cfg.getParameter<double>( "simple_sv_cut"          ) ),
  combSvCut_         ( cfg.getParameter<double>( "comb_sv_cut"            ) ),
  combSvMvaCut_      ( cfg.getParameter<double>( "comb_sv_mva_cut"        ) ),
  softMuonCut_       ( cfg.getParameter<double>( "soft_mu_cut"            ) ),
  softMuonNoIPCut_   ( cfg.getParameter<double>( "soft_mu_no_ip_cut"      ) ),
  softElecCut_       ( cfg.getParameter<double>( "soft_elec_cut"          ) ),
  
  etaBins_( cfg.getParameter<vector<double> >( "eta_bins" ) ),
  ptBins_ ( cfg.getParameter<vector<double> >( "pt_bins"  ) )    
{ 
}

BJetTag::~BJetTag()
{
}

void
BJetTag::fill(const pat::Jet& jet)
{  
  trkCountHighEffBJetTags_ ->Fill(jet.bDiscriminator("trackCountingHighEffBJetTags"      ));
  trkCountHighPurBJetTags_ ->Fill(jet.bDiscriminator("trackCountingHighPurBJetTags"      ));  
  jetProbBJetTags_         ->Fill(jet.bDiscriminator("jetProbabilityBJetTags"            ));  
  jetBProbBJetTags_        ->Fill(jet.bDiscriminator("jetBProbabilityBJetTags"           ));  
  impactParMvaBJetTags_    ->Fill(jet.bDiscriminator("impactParameterMVABJetTags"        ));  
  simpleSvBJetTags_        ->Fill(jet.bDiscriminator("simpleSecondaryVertexBJetTags"     )); 
  combSvBJetTags_          ->Fill(jet.bDiscriminator("combinedSecondaryVertexBJetTags"   ));  
  combSvMvaBJetTags_       ->Fill(jet.bDiscriminator("combinedSecondaryVertexMVABJetTags"));  
  softMuonBJetTags_        ->Fill(jet.bDiscriminator("softMuonBJetTags"                  ));  
  softMuonNoIPBJetTags_    ->Fill(jet.bDiscriminator("softMuonNoIPBJetTags"              ));  
  softElecBJetTags_        ->Fill(jet.bDiscriminator("softElectronBJetTags"              ));  
  
  for(size_t i=0; i<etaBins_.size(); ++i){
    if(fabs(jet.eta())>fabs(etaBins_[i]) && fabs(jet.eta()) < fabs(etaBins_[i+1])){ 
      etaDiscrTrkCountEff_[i] ->Fill(jet.bDiscriminator("trackCountingHighEffBJetTags"      ));
      etaDiscrTrkCountPur_[i] ->Fill(jet.bDiscriminator("trackCountingHighPurBJetTags"      ));      
      etaDiscrJetProb_[i]     ->Fill(jet.bDiscriminator("jetProbabilityBJetTags"            ));      
      etaDiscrJetBProb_[i]    ->Fill(jet.bDiscriminator("jetBProbabilityBJetTags"           ));     
      etaDiscrImpactParMva_[i]->Fill(jet.bDiscriminator("impactParameterMVABJetTags"        ));
      etaDiscrSimpleSv_[i]    ->Fill(jet.bDiscriminator("simpleSecondaryVertexBJetTags"     ));      
      etaDiscrCombSv_[i]      ->Fill(jet.bDiscriminator("combinedSecondaryVertexBJetTags"   ));      
      etaDiscrCombSvMva_[i]   ->Fill(jet.bDiscriminator("combinedSecondaryVertexMVABJetTags"));        
      etaDiscrSoftMu_[i]      ->Fill(jet.bDiscriminator("softMuonBJetTags"                  ));      
      etaDiscrSoftMuNoIP_[i]  ->Fill(jet.bDiscriminator("softMuonNoIPBJetTags"              ));      
      etaDiscrSoftElec_[i]    ->Fill(jet.bDiscriminator("softElectronBJetTags"              ));        
      break;      
    }
  }
  
  for(size_t i=0; i<ptBins_.size(); ++i){
    if(fabs(jet.pt())>fabs(ptBins_[i]) && fabs(jet.pt()) < fabs(ptBins_[i+1])){ 
      ptDiscrTrkCountEff_[i] ->Fill(jet.bDiscriminator("trackCountingHighEffBJetTags"      ));
      ptDiscrTrkCountPur_[i] ->Fill(jet.bDiscriminator("trackCountingHighPurBJetTags"      ));      
      ptDiscrJetProb_[i]     ->Fill(jet.bDiscriminator("jetProbabilityBJetTags"            ));      
      ptDiscrJetBProb_[i]    ->Fill(jet.bDiscriminator("jetBProbabilityBJetTags"           ));     
      ptDiscrImpactParMva_[i]->Fill(jet.bDiscriminator("impactParameterMVABJetTags"        ));
      ptDiscrSimpleSv_[i]    ->Fill(jet.bDiscriminator("simpleSecondaryVertexBJetTags"     ));      
      ptDiscrCombSv_[i]      ->Fill(jet.bDiscriminator("combinedSecondaryVertexBJetTags"   ));      
      ptDiscrCombSvMva_[i]   ->Fill(jet.bDiscriminator("combinedSecondaryVertexMVABJetTags"));        
      ptDiscrSoftMu_[i]      ->Fill(jet.bDiscriminator("softMuonBJetTags"                  ));      
      ptDiscrSoftMuNoIP_[i]  ->Fill(jet.bDiscriminator("softMuonNoIPBJetTags"              ));      
      ptDiscrSoftElec_[i]    ->Fill(jet.bDiscriminator("softElectronBJetTags"              ));        
      break;          
    }
  }
}

void 
BJetTag::book(ofstream& file, char* name)
{
  nBins = 100;

  edm::Service<TFileService> fs;

  NameScheme e(name);
  trkCountHighEffBJetTags_ = fs->make<TH1F>(e.name( file, "trkCountHighEffBJetTags" ), e.name("trkCountHighEffBJetTags"), nBins, 0. , 5.0 ); 
  trkCountHighPurBJetTags_ = fs->make<TH1F>(e.name( file, "trkCountHighPurBJetTags" ), e.name("trkCountHighPurBJetTags"), nBins, 0. , 5.0 );  
  jetProbBJetTags_         = fs->make<TH1F>(e.name( file, "jetProbBJetTags"         ), e.name("jetProbBJetTags"        ), nBins, 0. , 1.0 );  
  jetBProbBJetTags_        = fs->make<TH1F>(e.name( file, "jetBProbBJetTags"        ), e.name("jetBProbBJetTags"       ), nBins, 0. , 1.0 );  
  impactParMvaBJetTags_    = fs->make<TH1F>(e.name( file, "impactParMvaBJetTags"    ), e.name("impactParMvaBJetTags"   ), nBins, 0. , 1.0 );  
  simpleSvBJetTags_        = fs->make<TH1F>(e.name( file, "simpleSvBJetTags"        ), e.name("simpleSvBJetTags"       ), nBins, 0. , 5.0 );  
  combSvBJetTags_          = fs->make<TH1F>(e.name( file, "combSvBJetTags"          ), e.name("combSvBJetTags"         ), nBins, 0. , 1.0 );  
  combSvMvaBJetTags_       = fs->make<TH1F>(e.name( file, "combSvMvaBJetTags"       ), e.name("combSvMvaBJetTags"      ), nBins, 0. , 1.0 );  
  softMuonBJetTags_        = fs->make<TH1F>(e.name( file, "softMuonBJetTags"        ), e.name("softMuonBJetTags"       ), nBins, 0. , 1.0 );  
  softMuonNoIPBJetTags_    = fs->make<TH1F>(e.name( file, "softMuonNoIPBJetTags"    ), e.name("softMuonNoIPBJetTags"   ), nBins, 0. , 1.0 );  
  softElecBJetTags_        = fs->make<TH1F>(e.name( file, "softElecBJetTags"        ), e.name("softElecBJetTags"       ), nBins, 0. , 1.0 );  
    
  for(size_t i=0; i<etaBins_.size()-1; ++i){    
    TH1F* hist;  
    hist = fs->make<TH1F>(e.name(file, "etaBinsTrkCountEff", i) , e.name("etaBinsTrkCountEff", i) , nBins, 0., 1.0);
    etaDiscrTrkCountEff_.push_back(hist); 
    hist = fs->make<TH1F>(e.name(file, "etaBinsTrkCountPur", i) , e.name("etaBinsTrkCountPur", i) , nBins, 0., 1.0);
    etaDiscrTrkCountPur_.push_back(hist);     
    hist = fs->make<TH1F>(e.name(file, "etaBinsJetProb", i)     , e.name("etaBinsJetProb", i)     , nBins, 0., 1.0);
    etaDiscrJetProb_.push_back(hist); 
    hist = fs->make<TH1F>(e.name(file, "etaBinsJetBProb", i)    , e.name("etaBinsJetBProb", i)    , nBins, 0., 1.0);
    etaDiscrJetBProb_.push_back(hist);     
    hist = fs->make<TH1F>(e.name(file, "etaBinsImpactParMva", i), e.name("etaBinsImpactParMva", i), nBins, 0., 1.0);
    etaDiscrImpactParMva_.push_back(hist); 
    hist = fs->make<TH1F>(e.name(file, "etaBinsSimpleSv", i)    , e.name("etaBinsSimpleSv", i)    , nBins, 0., 1.0);
    etaDiscrSimpleSv_.push_back(hist);       
    hist = fs->make<TH1F>(e.name(file, "etaBinsCombSv", i)      , e.name("etaBinsCombSv", i)      , nBins, 0., 1.0);
    etaDiscrCombSv_.push_back(hist); 
    hist = fs->make<TH1F>(e.name(file, "etaBinsCombSvMva", i)   , e.name("etaBinsCombSvMva", i)   , nBins, 0., 1.0);
    etaDiscrCombSvMva_.push_back(hist);       
    hist = fs->make<TH1F>(e.name(file, "etaBinsSoftMu", i)      , e.name("etaBinsSoftMu", i)      , nBins, 0., 1.0);
    etaDiscrSoftMu_.push_back(hist); 
    hist = fs->make<TH1F>(e.name(file, "etaBinsSoftMuNoIP", i)  , e.name("etaBinsSoftMuNoIP", i)  , nBins, 0., 1.0);
    etaDiscrSoftMuNoIP_.push_back(hist);     
    hist = fs->make<TH1F>(e.name(file, "etaBinsSoftElec", i)    , e.name("etaBinsSoftElec", i)    , nBins, 0., 1.0);
    etaDiscrSoftElec_.push_back(hist);                 
  }
  
  for(size_t i=0; i<ptBins_.size()-1; ++i){
    TH1F* hist;  
    hist = fs->make<TH1F>(e.name(file, "ptBinsTrkCountEff", i) , e.name("ptBinsTrkCountEff", i) , nBins, 0., 1.0);
    ptDiscrTrkCountEff_.push_back(hist); 
    hist = fs->make<TH1F>(e.name(file, "ptBinsTrkCountPur", i) , e.name("ptBinsTrkCountPur", i) , nBins, 0., 1.0);
    ptDiscrTrkCountPur_.push_back(hist);     
    hist = fs->make<TH1F>(e.name(file, "ptBinsJetProb", i)     , e.name("ptBinsJetProb", i)     , nBins, 0., 1.0);
    ptDiscrJetProb_.push_back(hist); 
    hist = fs->make<TH1F>(e.name(file, "ptBinsJetBProb", i)    , e.name("ptBinsJetBProb", i)    , nBins, 0., 1.0);
    ptDiscrJetBProb_.push_back(hist); 
    hist = fs->make<TH1F>(e.name(file, "ptBinsImpactParMva", i), e.name("ptBinsImpactParMva", i), nBins, 0., 1.0);
    ptDiscrImpactParMva_.push_back(hist);     
    hist = fs->make<TH1F>(e.name(file, "ptBinsSimpleSv", i)    , e.name("ptBinsSimpleSv", i)    , nBins, 0., 1.0);
    ptDiscrSimpleSv_.push_back(hist);      
    hist = fs->make<TH1F>(e.name(file, "ptBinsCombSv", i)      , e.name("ptBinsCombSv", i)      , nBins, 0., 1.0);
    ptDiscrCombSv_.push_back(hist); 
    hist = fs->make<TH1F>(e.name(file, "ptBinsCombSvMva", i)   , e.name("ptBinsCombSvMva", i)   , nBins, 0., 1.0);
    ptDiscrCombSvMva_.push_back(hist);     
    hist = fs->make<TH1F>(e.name(file, "ptBinsSoftMu", i)      , e.name("ptBinsSoftMu", i)      , nBins, 0., 1.0);
    ptDiscrSoftMu_.push_back(hist);         
    hist = fs->make<TH1F>(e.name(file, "ptBinsSoftMuNoIP", i)  , e.name("ptBinsSoftMuNoIP", i)  , nBins, 0., 1.0);
    ptDiscrSoftMuNoIP_.push_back(hist); 
    hist = fs->make<TH1F>(e.name(file, "ptBinsSoftElec", i)    , e.name("ptBinsSoftElec", i)    , nBins, 0., 1.0);
    ptDiscrSoftElec_.push_back(hist);        
  }     
  
  size_t nEta = etaBins_.size()-1;
  
  etaEffTrkCountEff_ = fs->make<TH1F>(e.name( file, "etaEffTrkCountEff" ), e.name("etaEffTrkCountEff" ), nEta, &etaBins_.front());
  etaEffTrkCountPur_ = fs->make<TH1F>(e.name( file, "etaEffTrkCountPur" ), e.name("etaEffTrkCountPur" ), nEta, &etaBins_.front());  
  etaEffJetProb_     = fs->make<TH1F>(e.name( file, "etaEffJetProb"     ), e.name("etaEffJetProb"     ), nEta, &etaBins_.front());  
  etaEffJetBProb_    = fs->make<TH1F>(e.name( file, "etaEffJetBProb"    ), e.name("etaEffJetBProb"    ), nEta, &etaBins_.front());  
  etaEffImpactParMva_= fs->make<TH1F>(e.name( file, "etaEffImpactParMva"), e.name("etaEffImpactParMva"), nEta, &etaBins_.front());  
  etaEffSimpleSv_    = fs->make<TH1F>(e.name( file, "etaEffSimpleSv"    ), e.name("etaEffSimpleSv"    ), nEta, &etaBins_.front()); 
  etaEffCombSv_      = fs->make<TH1F>(e.name( file, "etaEffCombSv"      ), e.name("etaEffCombSv"      ), nEta, &etaBins_.front());  
  etaEffCombSvMva_   = fs->make<TH1F>(e.name( file, "etaEffCombSvMva"   ), e.name("etaEffCombSvMva"   ), nEta, &etaBins_.front());  
  etaEffSoftMu_      = fs->make<TH1F>(e.name( file, "etaEffSoftMu"      ), e.name("etaEffSoftMu"      ), nEta, &etaBins_.front());  
  etaEffSoftMuNoIP_  = fs->make<TH1F>(e.name( file, "etaEffSoftMuNoIP"  ), e.name("etaEffSoftMuNoIP"  ), nEta, &etaBins_.front());  
  etaEffSoftElec_    = fs->make<TH1F>(e.name( file, "etaEffSoftElec"    ), e.name("etaEffSoftElec"    ), nEta, &etaBins_.front());   
  
  size_t nPt = ptBins_.size()-1;
   
  ptEffTrkCountEff_ = fs->make<TH1F>(e.name( file, "ptEffTrkCountEff" ), e.name("ptEffTrkCountEff" ), nPt, &ptBins_.front()); 
  ptEffTrkCountPur_ = fs->make<TH1F>(e.name( file, "ptEffTrkCountPur" ), e.name("ptEffTrkCountPur" ), nPt, &ptBins_.front());   
  ptEffJetProb_     = fs->make<TH1F>(e.name( file, "ptEffJetProb"     ), e.name("ptEffJetProb"     ), nPt, &ptBins_.front());   
  ptEffJetBProb_    = fs->make<TH1F>(e.name( file, "ptEffJetBProb"    ), e.name("ptEffJetBProb"    ), nPt, &ptBins_.front());   
  ptEffImpactParMva_= fs->make<TH1F>(e.name( file, "ptEffImpactParMva"), e.name("ptEffImpactParMva"), nPt, &ptBins_.front()); 
  ptEffSimpleSv_    = fs->make<TH1F>(e.name( file, "ptEffSimpleSv"    ), e.name("ptEffSimpleSv"    ), nPt, &ptBins_.front());   
  ptEffCombSv_      = fs->make<TH1F>(e.name( file, "ptEffCombSv"      ), e.name("ptEffCombSv"      ), nPt, &ptBins_.front());   
  ptEffCombSvMva_   = fs->make<TH1F>(e.name( file, "ptEffCombSvMva"   ), e.name("ptEffCombSvMva"   ), nPt, &ptBins_.front());     
  ptEffSoftMu_      = fs->make<TH1F>(e.name( file, "ptEffSoftMu"      ), e.name("ptEffSoftMu"      ), nPt, &ptBins_.front()); 
  ptEffSoftMuNoIP_  = fs->make<TH1F>(e.name( file, "ptEffSoftMuNoIP"  ), e.name("ptEffSoftMuNoIP"  ), nPt, &ptBins_.front());   
  ptEffSoftElec_    = fs->make<TH1F>(e.name( file, "ptEffSoftElec"    ), e.name("ptEffSoftElec"    ), nPt, &ptBins_.front());     
}

void
BJetTag::evaluate()
{  
  size_t nEta = etaBins_.size()-1;

  fillEfficencyHist(etaEffTrkCountEff_ , etaDiscrTrkCountEff_ , nEta, trkCountHighEffCut_);
  fillEfficencyHist(etaEffTrkCountPur_ , etaDiscrTrkCountPur_ , nEta, trkCountHighPurCut_);  
  fillEfficencyHist(etaEffJetProb_     , etaDiscrJetProb_     , nEta, jetProbCut_        );  
  fillEfficencyHist(etaEffJetBProb_    , etaDiscrJetBProb_    , nEta, jetBProbCut_       );  
  fillEfficencyHist(etaEffImpactParMva_, etaDiscrImpactParMva_, nEta, impactParMvaCut_   );  
  fillEfficencyHist(etaEffSimpleSv_    , etaDiscrSimpleSv_    , nEta, simpleSvCut_       );  
  fillEfficencyHist(etaEffCombSv_      , etaDiscrCombSv_      , nEta, combSvCut_         );  
  fillEfficencyHist(etaEffCombSvMva_   , etaDiscrCombSvMva_   , nEta, combSvMvaCut_      );  
  fillEfficencyHist(etaEffSoftMu_      , etaDiscrSoftMu_      , nEta, softMuonCut_       );
  fillEfficencyHist(etaEffSoftMuNoIP_  , etaDiscrSoftMuNoIP_  , nEta, softMuonNoIPCut_   );
  fillEfficencyHist(etaEffSoftElec_    , etaDiscrSoftElec_    , nEta, softElecCut_       ); 
  
  size_t nPt = ptBins_.size()-1;

  fillEfficencyHist(ptEffTrkCountEff_ , ptDiscrTrkCountEff_ , nPt, trkCountHighEffCut_);
  fillEfficencyHist(ptEffTrkCountPur_ , ptDiscrTrkCountPur_ , nPt, trkCountHighPurCut_);  
  fillEfficencyHist(ptEffJetProb_     , ptDiscrJetProb_     , nPt, jetProbCut_        );  
  fillEfficencyHist(ptEffJetBProb_    , ptDiscrJetBProb_    , nPt, jetBProbCut_       );  
  fillEfficencyHist(ptEffImpactParMva_, ptDiscrImpactParMva_, nPt, impactParMvaCut_   );  
  fillEfficencyHist(ptEffSimpleSv_    , ptDiscrSimpleSv_    , nPt, simpleSvCut_       );  
  fillEfficencyHist(ptEffCombSv_      , ptDiscrCombSv_      , nPt, combSvCut_         );  
  fillEfficencyHist(ptEffCombSvMva_   , ptDiscrCombSvMva_   , nPt, combSvMvaCut_      );  
  fillEfficencyHist(ptEffSoftMu_      , ptDiscrSoftMu_      , nPt, softMuonCut_       );
  fillEfficencyHist(ptEffSoftMuNoIP_  , ptDiscrSoftMuNoIP_  , nPt, softMuonNoIPCut_   );
  fillEfficencyHist(ptEffSoftElec_    , ptDiscrSoftElec_    , nPt, softElecCut_       );   
}

void
BJetTag::fillEfficencyHist(TH1F* hist, vector<TH1F*> src, size_t n, double cut)
{
  for(size_t i=0; i<n; ++i){ 
    int cutBin = cut*nBins+1;   
    double integral = src[i]->Integral(cutBin, nBins+1);      
    double norm = src[i]->GetEntries();                  
    if(norm==0) continue; 
    double efficency = integral/norm; 
    double error = binomerror(integral, norm); 
    hist->SetBinContent(i+1, efficency);
    hist->SetBinError(i+1, error); 
  }  
}

double
BJetTag::binomerror(double numerator, double denominator)
{
  if(numerator>denominator) return 999999.;
  
  double error;  
  error=sqrt((1-numerator/denominator)
   *numerator/denominator/denominator);
   
  return error;
}
