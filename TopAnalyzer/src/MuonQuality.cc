#include "TopAnalysis/TopAnalyzer/interface/MuonQuality.h"

/// default constructor for fw lite
MuonQuality::MuonQuality(const int index) : norm_(0.), index_(index)
{ 
}

/// constructor for full FW analyzer
MuonQuality::MuonQuality(const edm::ParameterSet& cfg) : norm_(0.), 
  index_( cfg.getParameter<int>( "index" ) )
{
}

/// histogramm booking for fwlite 
void 
MuonQuality::book()
{
  /** 
      Selection Variables
  **/
  // number of valid hits in silicon tracker
  hists_["nHit"     ] = new TH1F( "nHit"    ,  "nHit"     ,   35,    -0.5,  34.5 );
  // normalized chi2 of global muon track fit
  hists_["chi2"     ] = new TH1F( "chi2"    ,  "chi2"     ,   30,     0. ,  30.  );
  // d0 significance of track (still to nominal IP)
  hists_["dB"       ] = new TH1F( "dB"      ,  "dB"       ,  200,     0. ,  1.   );
  // dz significance of track in z-dimension 
  hists_["dz"       ] = new TH1F( "dz"      ,  "dz"       ,  150,   -15. ,  15.  );
  // energy in ecal attached to the candidate trajectory
  hists_["ecalEn"   ] = new TH1F( "ecalEn"  ,  "ecalEn"   ,   40,     0. ,  10.  );
  // energy in hcal attached to the candidate trajectory
  hists_["hcalEn"   ] = new TH1F( "hcalEn"  ,  "hcalEn"   ,   40,     0. ,  10.  );
  // relative isolation (tracker and calo combined)
  hists_["relIso"   ] = new TH1F( "relIso"  ,  "relIso"   ,   50,     0. ,  1.0  );
  // is GlobalMuonPromptTight?
  hists_["isGlobalMuonPromptTight"] = new TH1F("isGlobalMuonPromptTight","isGlobalMuonPromptTight"  ,   2,     0. ,  2.  ); 
  // is TrackerMuon?
  hists_["isTrackerMuon"          ] = new TH1F("isTrackerMuon"          ,  "isTrackerMuon"          ,   2,     0. ,  2.  ); 
  // number of matches
  hists_["matches"                ] = new TH1F("matches"                ,  "matches"                ,  10,     0. , 10.  );
  
  /** 
      Monitoring Variables
  **/
  // compatibility of the energy deposits in the calorimeter with the muon hypoothesis
  hists_["calCmp_" ] = new TH1F( "calCmp_"  ,  "calCmp_"  ,   50,   0.,   1. );
  // energy deposited in crossed ecal crystals (recHist based)
  hists_["ecalEn_" ] = new TH1F( "ecalEn_"  ,  "ecalEn_"  ,   40,   0.,  10. ); 
  // energy deposited in 3x3 crystal shape around cerntral crystal (recHits based)
  hists_["ecalS9_" ] = new TH1F( "ecalS9_"  ,  "ecalS9_"  ,   40,   0.,  10. ); 
  // energy deposited in crossed hcal towers (recHits based)
  hists_["hcalEn_" ] = new TH1F( "hcalEn_"  ,  "hcalEn_"  ,   40,   0.,  10. ); 
  // energy deposited in 3x3 tower shape around cerntral tower (recHits based)
  hists_["hcalS9_" ] = new TH1F( "hcalS9_"  ,  "hcalS9_"  ,   40,   0.,  10. );
  // energy deposited in crossed hcal outer towers (recHits based)
  hists_["houtEn_" ] = new TH1F( "houtEn_"  ,  "houtEn_"  ,   40,   0.,  10. ); 
  // energy deposited in 3x3 tower shape around cerntral tower (recHits based)
  hists_["houtS9_" ] = new TH1F( "houtS9_"  ,  "houtS9_"  ,   40,   0.,  10. );   

  /** 
      Isolation Monitoring Variables
  **/
  // number of tracks in isolation cone
  hists_["trkIsoN_"] = new TH1F( "trkIsoN_" ,  "isoTkN_"  ,   20,   0.,  20. );
  // summed track pt in isolation cone
  hists_["trkIso_" ] = new TH1F( "trkIso_"  ,  "trkIso_"  ,   40,   0.,  10. );
  // number of ecal crystals above noise threshold in isolation cone
  hists_["eclIsoN_"] = new TH1F( "eclIsoN_" ,  "eclIsoN_" ,   30,   0.,  30. );
  // number of hcal towers above noise threshold in isolation cone
  hists_["hclIsoN_"] = new TH1F( "hclIsoN_" ,  "hclIsoN_" ,   30,   0.,  30. );
  // summed energy in ecal crystals and hcal towers above noise threshold in isolation cone (recHit based)
  hists_["calIso_" ] = new TH1F( "calIso_"  ,  "calIso_"  ,   40,   0.,  10. );

  /** 
      Energy & Object Flow Variables
  **/
  // conserve the normalization information (for the use with hadd)
  hists_["norm_"   ] = new TH1F( "norm_"    ,  "norm_"    ,    1,   0.,   1. );
  // <number of tracks> as a function of deltaR 
  hists_["trkDRN_" ] = new TH1F( "trkDRN_"  ,  "trkDRN_"  ,   20,   0.,   1. );
  // <summerd pt of tracks> as a function of deltaR 
  hists_["trkDR_"  ] = new TH1F( "trkDR_"   ,  "trkDR_"   ,   20,   0.,   1. );
  // <number of ecal crystals above noise threshold> as a function of deltaR 
  hists_["eclDRN_" ] = new TH1F( "eclDRN_"  ,  "eclDRN_"  ,   20,   0.,   1. );
  // <number of hcal towers above noise threshold> as a function of deltaR 
  hists_["hclDRN_" ] = new TH1F( "hclDRN_"  ,  "hclDRN_"  ,   20,   0.,   1. );
  // <summerd crystal and tower energies above noise threshold> as a function of deltaR 
  hists_["calDR_"  ] = new TH1F( "calDR_"   ,  "calDR_"   ,   20,   0.,   1. );

  /** 
      uncorrelated variables for ABCD QCD estimation-method
  **/
  // impact parameter with respect to the beamspot vs. relative isolation
  hists2D_["relIsoVsDb_"] = new TH2F( "relIsoVsDb_"   , "relIsoVsdB_"   ,  100, 0., 1.0,  1000, 0., 1.);
}

/// histogramm booking for full fw
void 
MuonQuality::book(edm::Service<TFileService>& fs)
{
  /** 
      Selection Variables
  **/
  // number of valid hits in silicon tracker
  hists_["nHit"     ] = fs->make<TH1F>( "nHit"    ,  "nHit"     ,   35,  -0.5, 34.5  );
  // normalized chi2 of global muon track fit
  hists_["chi2"     ] = fs->make<TH1F>( "chi2"    ,  "chi2"     ,   30,    0.,  30.  );
  // d0 significance of track (still to nominal IP)
  hists_["dB"       ] = fs->make<TH1F>( "dB"      ,  "dB"       ,   200,   0.,  1.   );
  // dz significance of track in z-dimension
  hists_["dz"       ] = fs->make<TH1F>( "dz"      ,  "dz"       ,   150, -15., 15.   );
  // energy in ecal attached to the candidate trajectory
  hists_["ecalEn"   ] = fs->make<TH1F>( "ecalEn"  ,  "ecalEn"   ,   40,    0.,  10.  );
  // energy in hcal attached to the candidate trajectory
  hists_["hcalEn"   ] = fs->make<TH1F>( "hcalEn"  ,  "hcalEn"   ,   40,    0.,  10.  );
  // relative isolation (tracker and calo combined)
  hists_["relIso"   ] = fs->make<TH1F>( "relIso"  ,  "relIso"   ,   50,    0.,  1.0  );
  // is GlobalMuonPromptTight?
  hists_["isGlobalMuonPromptTight"] = fs->make<TH1F>("isGlobalMuonPromptTight" , "isGlobalMuonPromptTight" ,   2,     0. ,  2.  ); 
  // is TrackerMuon?
  hists_["isTrackerMuon"          ] = fs->make<TH1F>("isTrackerMuon"           , "isTrackerMuon"           ,   2,     0. ,  2.  ); 
  // number of matches
  hists_["matches"                ] = fs->make<TH1F>("matches"                 , "matches"                 ,  10,     0. ,  10. );
 

  /** 
      Monitoring Variables
  **/
  // compatibility of the energy deposits in the calorimeter with the muon hypothesis
  hists_["calCmp_" ] = fs->make<TH1F>( "calCmp_"   ,  "calCmp_"  ,  50,   0.,   1. );
  // energy deposited in crossed ecal crystals (recHist based)
  hists_["ecalEn_" ] = fs->make<TH1F>( "ecalEn_"  ,  "ecalEn_"  ,   40,   0.,  10. ); 
  // energy deposited in 3x3 crystal shape around cerntral crystal (recHits based)
  hists_["ecalS9_" ] = fs->make<TH1F>( "ecalS9_"  ,  "ecalS9_"  ,   40,   0.,  10. ); 
  // energy deposited in crossed hcal towers (recHits based)
  hists_["hcalEn_" ] = fs->make<TH1F>( "hcalEn_"  ,  "hcalEn_"  ,   40,   0.,  10. ); 
  // energy deposited in 3x3 tower shape around cerntral tower (recHits based)
  hists_["hcalS9_" ] = fs->make<TH1F>( "hcalS9_"  ,  "hcalS9_"  ,   40,   0.,  10. );
  // energy deposited in crossed hcal outer towers (recHits based)
  hists_["houtEn_" ] = fs->make<TH1F>( "houtEn_"  ,  "houtEn_"  ,   40,   0.,  10. ); 
  // energy deposited in 3x3 tower shape around cerntral tower (recHits based)
  hists_["houtS9_" ] = fs->make<TH1F>( "houtS9_"  ,  "houtS9_"  ,   40,   0.,  10. );   

  /** 
      Isolation Monitoring Variables
  **/
  // conserve the normalization information (for the use with hadd)
  hists_["norm_"   ] = fs->make<TH1F>( "norm_"    ,  "norm_"    ,    1,   0.,   1. );
  // number of tracks in isolation cone
  hists_["trkIsoN_"] = fs->make<TH1F>( "trkIsoN_" ,  "isoTkN_"  ,   20,   0.,  20. );
  // summed track pt in isolation cone
  hists_["trkIso_" ] = fs->make<TH1F>( "trkIso_"  ,  "trkIso_"  ,   40,   0.,  10. );
  // number of ecal crystals above noise threshold in isolation cone
  hists_["eclIsoN_"] = fs->make<TH1F>( "eclIsoN_" ,  "eclIsoN_" ,   30,   0.,  30. );
  // number of hcal towers above noise threshold in isolation cone
  hists_["hclIsoN_"] = fs->make<TH1F>( "hclIsoN_" ,  "hclIsoN_" ,   30,   0.,  30. );
  // summed energy in ecal crystals and hcal towers above noise threshold in isolation cone (recHit based)
  hists_["calIso_" ] = fs->make<TH1F>( "calIso_"  ,  "calIso_"  ,   40,   0.,  10. );

  /** 
      Energy & Object Flow Variables
  **/
  // <number of tracks> as a function of deltaR (differential)
  hists_["trkDRN_" ] = fs->make<TH1F>( "trkDRN_"  ,  "trkDRN_"  ,   20,   0.,   1. );
  // <summed pt of tracks> as a function of deltaR (differential)
  hists_["trkDR_"  ] = fs->make<TH1F>( "trkDR_"   ,  "trkDR_"   ,   20,   0.,   1. );
  // <number of ecal crystals above noise threshold> as a function of deltaR (differential)
  hists_["eclDRN_" ] = fs->make<TH1F>( "eclDRN_"  ,  "eclDRN_"  ,   20,   0.,   1. );
  // <number of hcal towers above noise threshold> as a function of deltaR (differential)
  hists_["hclDRN_" ] = fs->make<TH1F>( "hclDRN_"  ,  "hclDRN_"  ,   20,   0.,   1. );
  // <summed crystal and tower energies above noise threshold> as a function of deltaR (differential)
  hists_["calDR_"  ] = fs->make<TH1F>( "calDR_"   ,  "calDR_"   ,   20,   0.,   1. );

  /** 
      uncorrelated variables for ABCD QCD estimation-method
  **/
  // impact parameter with respect to the beamspot vs. relative isolation
  hists2D_["relIsoVsDb_"] = fs->make<TH2F>( "relIsoVsDb_", "relIsoVsdB_",  100, 0., 1.0,  1000, 0., 1.);
}

/// histogram filling for fwlite and for full fw from reco objects
void
MuonQuality::fill(const edm::View<pat::Muon>& muons, const double& weight)
{      
  int index=0;
  for(edm::View<pat::Muon>::const_iterator muon=muons.begin(); muon!=muons.end(); ++muon, ++index){
    // NOTE: against the common policy *not* to have any implicit selection cuts 
    // within analyzers these plots are restricted to global muons, as otherwise
    // some of the monitor histograms cannot be filled
    if( (index_<0 || index_==index) && muon->isGlobalMuon() ){

      /**
	 Fill Selection Variables
      **/
      // number of valid hits in silicon tracker
      hists_.find("nHit")->second->Fill( muon->track()->numberOfValidHits(), weight );  
      // normalized chi2 of global muon track fit
      hists_.find("chi2")->second->Fill( muon->combinedMuon()->normalizedChi2(), weight ); 
      // d0 significance of track (still to nominal IP)
      hists_.find("dB"  )->second->Fill( muon->dB(), weight );
      // dz significance of track in z-dimension
      hists_.find("dz"  )->second->Fill( muon->track()->dz(), weight );
      // energy in ecal attached to the candidate trajectory
      //hists_.find("ecalEn")->second->Fill( muon->ecalIsoDeposit()->candEnergy(), weight );
      // energy in hcal attached to the candidate trajectory
      //hists_.find("hcalEn")->second->Fill( muon->hcalIsoDeposit()->candEnergy(), weight );   
      // relative isolation (tracker and calo combined)
      //hists_.find("relIso")->second->Fill( (muon->trackIso()+muon->caloIso())/muon->pt() , weight );
      // PF relIso
      hists_.find("relIso")->second->Fill( (muon->chargedHadronIso() + muon->neutralHadronIso() + 
					    muon->photonIso()) / muon->pt() , weight );
      // is GlobalMuonPromptTight?
      hists_.find("isGlobalMuonPromptTight")->second->Fill(muon->isGood("GlobalMuonPromptTight") , weight );
      // is TrackerMuon?
      hists_.find("isTrackerMuon")->second->Fill(muon->isGood("AllTrackerMuons") , weight );
      // number of matches
      hists_.find("matches")->second->Fill(muon->numberOfMatches() , weight );
       
      /** 
	  Fill Monitoring Variables
      **/
      // compatibility of the energy deposits in the calorimeter with the muon hypothesis
      //hists_.find("calCmp_" )->second->Fill( muon->caloCompatibility() , weight );
      // energy deposited in crossed ecal crystals (recHist based)
      //hists_.find("ecalEn_" )->second->Fill( muon->calEnergy().em , weight );
      // energy deposited in 3x3 crystal shape around cerntral crystal (recHits based)    
      //hists_.find("ecalS9_" )->second->Fill( muon->calEnergy().emS9 , weight );  
      // energy deposited in crossed hcal towers (recHits based)
      //hists_.find("hcalEn_" )->second->Fill( muon->calEnergy().had , weight );
      // energy deposited in 3x3 tower shape around cerntral tower (recHits based)
      //hists_.find("hcalS9_" )->second->Fill( muon->calEnergy().hadS9 , weight );  
      // energy deposited in crossed hcal outer towers (recHits based)
      //hists_.find("houtEn_" )->second->Fill( muon->calEnergy().ho , weight );
      // energy deposited in 3x3 tower shape around cerntral tower (recHits based)
      //hists_.find("houtS9_" )->second->Fill( muon->calEnergy().hoS9 , weight );  

      /** 
	  Fill Isolation Monitoring Variables
      **/
      //static const double STANDARD_CONE_SIZE = 0.3;
      // number of tracks in isolation cone
      //hists_.find("trkIsoN_" )->second->Fill( muon->trackIsoDeposit()->countWithin( STANDARD_CONE_SIZE ) , weight );
      // summed track pt in isolation cone
      //hists_.find("trkIso_"  )->second->Fill( muon->trackIso() , weight );
      // number of ecal crystals above noise threshold in isolation cone
      //hists_.find("eclIsoN_" )->second->Fill( muon->ecalIsoDeposit()->countWithin( STANDARD_CONE_SIZE ) , weight );
      // number of hcal towers above noise threshold in isolation cone
      //hists_.find("hclIsoN_" )->second->Fill( muon->hcalIsoDeposit()->countWithin( STANDARD_CONE_SIZE ) , weight );
      // summed energy in ecal crystals and hcal towers above noise threshold in isolation cone (recHit based)
      //hists_.find("calIso_"  )->second->Fill( muon->caloIso () , weight );

      /** 
	  Fill Energy & Object Flow Variables
      **/
      // increment histogram normalization
      norm_+=weight;
      // <number of tracks> as a function of deltaR (differential) 
      //objectFlow( hists_.find("trkDRN_" )->second , muon->trackIsoDeposit() );
      // <summerd pt of tracks> as a function of deltaR (differential)
      //energyFlow( hists_.find("trkDR_"  )->second , muon->trackIsoDeposit() );
      // <number of ecal crystals above noise threshold> as a function of deltaR (differential)
      //objectFlow( hists_.find("eclDRN_" )->second , muon->ecalIsoDeposit() );
      // <number of hcal towers above noise threshold> as a function of deltaR (differential)
      //objectFlow( hists_.find("hclDRN_" )->second , muon->hcalIsoDeposit() );
      // <summerd crystal and tower energies above noise threshold> as a function of deltaR (differential)
      //energyFlow( hists_.find("calDR_"  )->second , muon->ecalIsoDeposit(), muon->hcalIsoDeposit() );

      /** 
	  fill relIso vs. dB hitogram for ABCD QCD-estimation method
      **/
      // transverse momentum of the muon
      //hists2D_.find("relIsoVsDb_")->second->Fill( (muon->trackIso()+muon->caloIso())/muon->pt(), muon->dB(), weight );
    }
  }
}

/// everything which needs to be done after the event loop
void 
MuonQuality::process()
{
  // fill normalization histogram with the number of muons, 
  // that went into the eflow histograms (weighted); we 
  // need this workaround to be able to normalize even 
  // after the use of hadd, when running in batch mode...
  hists_.find("norm_" )->second->SetBinContent( 1, norm_ );
}

/// get number of objects within a ring in deltaR corresponding to the bin width 
/// of the histogram 'hist' from 'deposit' and fill hist with it
void
MuonQuality::objectFlow(TH1* hist, const pat::IsoDeposit* deposit)
{
  if(hist && deposit){
    for(int bin=1; bin<=hist->GetNbinsX(); ++bin){
      double lowerEdge = hist->GetBinLowEdge(bin);
      double upperEdge = hist->GetBinLowEdge(bin)+hist->GetBinWidth(bin);
      hist->Fill(hist->GetBinCenter(bin), deposit->countWithin(upperEdge) - deposit->countWithin(lowerEdge));
    }
  }
}

/// get energy of objects within a ring in deltaR corresponding to the bin width 
/// of the histogram 'hist' from 'deposit' and fill hist with it
void
MuonQuality::energyFlow(TH1* hist, const pat::IsoDeposit* deposit)
{
  if(hist && deposit){
    for(int bin=1; bin<=hist->GetNbinsX(); ++bin){
      double lowerEdge = hist->GetBinLowEdge(bin);
      double upperEdge = hist->GetBinLowEdge(bin)+hist->GetBinWidth(bin);
      // restrict the differential energy deposit to less than 10GeV
      // in order not to have too high weight in the plot
      if(deposit->depositWithin(upperEdge) - deposit->depositWithin(lowerEdge) < 10)
	hist->Fill(hist->GetBinCenter(bin), deposit->depositWithin(upperEdge) - deposit->depositWithin(lowerEdge));
    }
  }
}

/// get energy of objects within a ring in deltaR corresponding to the bin width 
/// of the histogram 'hist' from 'deposit' and fill hist with it
void
MuonQuality::energyFlow(TH1* hist, const pat::IsoDeposit* ecalDeposit, const pat::IsoDeposit* hcalDeposit)
{
  if(hist && ecalDeposit && hcalDeposit){
    for(int bin=1; bin<=hist->GetNbinsX(); ++bin){
      double lowerEdge = hist->GetBinLowEdge(bin);
      double upperEdge = hist->GetBinLowEdge(bin)+hist->GetBinWidth(bin);
      hist->Fill(hist->GetBinCenter(bin), ecalDeposit->depositWithin(upperEdge)-ecalDeposit->depositWithin(lowerEdge) +  hcalDeposit->depositWithin(upperEdge)-hcalDeposit->depositWithin(lowerEdge) );
    }
  }
}
