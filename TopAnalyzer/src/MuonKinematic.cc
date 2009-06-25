#include "TopAnalysis/TopAnalyzer/interface/MuonKinematic.h"

/// constructor for FWLite analyzer
MuonKinematic::MuonKinematic():
  fwLite_(true ), jets_() { }

/// constructor for full FW analyzer
MuonKinematic::MuonKinematic(const edm::ParameterSet& cfg):
  fwLite_(false), jets_  ( cfg.getParameter<edm::InputTag>( "jets" ) )
{
}

/// fill interface for full FW analyzer
void
MuonKinematic::fill(const edm::Event& evt, const std::vector<pat::Muon>& muons, const double& weight=1.)
{
  if(muons.begin()==muons.end()) return;

  // jet isolation
  edm::Handle<std::vector<pat::Jet> > jets;
  evt.getByLabel(jets_, jets);

  fill(*jets, muons, weight);
}

/// fill interface for FWLite analyzer
void
MuonKinematic::fill(const std::vector<pat::Jet>& jets,
		    const std::vector<pat::Muon>& muons, const double& weight=1.)
{
  std::vector<pat::Muon>::const_iterator muon=muons.begin();
  if(muon!=muons.end()){
    //---------------------------------------------
    // basic kinematics
    //---------------------------------------------
    en_ ->Fill( muon->energy(), weight );
    pt_ ->Fill( muon->et(),     weight );
    eta_->Fill( muon->eta(),    weight );
    phi_->Fill( muon->phi(),    weight );

    //---------------------------------------------
    // jet isolation
    //---------------------------------------------
    double minDR   =-1.;
    double minDR_5 =-1.;
    double minDR_10=-1.;
    double minDR_15=-1.;
    double minDR_20=-1.;
    double minDR_25=-1.;
    double minDR_30=-1.;
    double minDR_35=-1.;
    double minDR_40=-1.;
    for(std::vector<pat::Jet>::const_iterator jet = jets.begin();
	jet!=jets.end(); ++jet) {
      double dR=deltaR(muon->eta(), muon->phi(), jet->eta(), jet->phi());
      if(  minDR   <0 || dR<minDR   )                   minDR   =dR;
      if( (minDR_5 <0 || dR<minDR_5 ) && jet->pt()>5  ) minDR_5 =dR;
      if( (minDR_10<0 || dR<minDR_10) && jet->pt()>10 ) minDR_10=dR;
      if( (minDR_15<0 || dR<minDR_15) && jet->pt()>15 ) minDR_15=dR;
      if( (minDR_20<0 || dR<minDR_20) && jet->pt()>20 ) minDR_20=dR;
      if( (minDR_25<0 || dR<minDR_25) && jet->pt()>25 ) minDR_25=dR;
      if( (minDR_30<0 || dR<minDR_30) && jet->pt()>30 ) minDR_30=dR;
      if( (minDR_35<0 || dR<minDR_35) && jet->pt()>35 ) minDR_35=dR;
      if( (minDR_40<0 || dR<minDR_40) && jet->pt()>40 ) minDR_40=dR;
    }
    if( minDR   >=0 ) isoJet_  ->Fill( minDR   , weight );
    if( minDR_5 >=0 ) isoJet5_ ->Fill( minDR_5 , weight );
    if( minDR_10>=0 ) isoJet10_->Fill( minDR_10, weight );
    if( minDR_15>=0 ) isoJet15_->Fill( minDR_15, weight );
    if( minDR_20>=0 ) isoJet20_->Fill( minDR_20, weight );
    if( minDR_25>=0 ) isoJet25_->Fill( minDR_25, weight );
    if( minDR_30>=0 ) isoJet30_->Fill( minDR_30, weight );
    if( minDR_35>=0 ) isoJet35_->Fill( minDR_35, weight );
    if( minDR_40>=0 ) isoJet40_->Fill( minDR_40, weight );

    //---------------------------------------------
    // track & calo isolation
    //---------------------------------------------
    const pat::IsoDeposit *trackDep=0, *ecalDep=0, *hcalDep=0;
    trackDep = muon->trackerIsoDeposit();
    ecalDep  = muon->ecalIsoDeposit();
    hcalDep  = muon->hcalIsoDeposit();

    // fill enegry flow and object flow histograms
    for(int bin = 1; bin <= dRTrkN_->GetNbinsX(); ++bin){
      dRTrkN_  ->Fill( dRTrkN_  ->GetBinCenter(bin),
		       trackDep ->countWithin  (dRTrkN_  ->GetBinLowEdge(bin)+dRTrkN_  ->GetBinWidth(bin)) -
		       trackDep ->countWithin  (dRTrkN_  ->GetBinLowEdge(bin))
		       );
      dREcalN_ ->Fill( dREcalN_ ->GetBinCenter(bin),
		       ecalDep  ->countWithin  (dREcalN_ ->GetBinLowEdge(bin)+dREcalN_ ->GetBinWidth(bin)) -
		       ecalDep  ->countWithin  (dREcalN_ ->GetBinLowEdge(bin))
		       );
      dRHcalN_ ->Fill( dRHcalN_ ->GetBinCenter(bin),
		       hcalDep  ->countWithin  (dRHcalN_ ->GetBinLowEdge(bin)+dRHcalN_ ->GetBinWidth(bin)) -
		       hcalDep  ->countWithin  (dRHcalN_ ->GetBinLowEdge(bin))
		       );
      dRTrkPt_ ->Fill( dRTrkPt_ ->GetBinCenter(bin),
		       trackDep ->depositWithin(dRTrkPt_ ->GetBinLowEdge(bin)+dRTrkPt_ ->GetBinWidth(bin)) -
		       trackDep ->depositWithin(dRTrkPt_ ->GetBinLowEdge(bin))
		       );
      dREcalPt_->Fill( dREcalPt_->GetBinCenter(bin),
		       ecalDep  ->depositWithin(dREcalPt_->GetBinLowEdge(bin)+dREcalPt_->GetBinWidth(bin)) -
		       ecalDep  ->depositWithin(dREcalPt_->GetBinLowEdge(bin))
		       );
      dRHcalPt_->Fill( dRHcalPt_->GetBinCenter(bin),
		       hcalDep  ->depositWithin(dRHcalPt_->GetBinLowEdge(bin)+dRHcalPt_->GetBinWidth(bin)) -
		       hcalDep  ->depositWithin(dRHcalPt_->GetBinLowEdge(bin))
		       );
    }
    // fill summed deposits
    double stdDR = 0.3;
    isoEcalN_->Fill( ecalDep ->countWithin(stdDR), weight );
    isoHcalN_->Fill( hcalDep ->countWithin(stdDR), weight );
    isoTrkN_ ->Fill( trackDep->countWithin(stdDR), weight );
    isoTrkPt_->Fill( muon->trackIso(), weight );
    isoCalPt_->Fill( muon->caloIso (), weight );

    double relIso = muon->pt() / ( muon->pt() + muon->trackIso() + muon->caloIso() );
    double relCombIso = (muon->trackIso() + muon->caloIso())/muon->pt();
    isoRelPt_->Fill( relIso, weight );
    isoRelComb_->Fill( relCombIso, weight );

    // fill correlation plots
    ptVsTrkIso_->Fill( muon->pt(), muon->trackIso() );
    ptVsCalIso_->Fill( muon->pt(), muon->caloIso () );
  }
}

/// book for FWLite
void
MuonKinematic::book()
{
  NameScheme kin("kin");
  en_ = new TH1F(kin.name("en" ), "en" , 60,   0., 300.);
  pt_ = new TH1F(kin.name("pt" ), "pt" , 30,   0., 150.);
  eta_= new TH1F(kin.name("eta"), "eta", 34, -3.4,  3.4);
  phi_= new TH1F(kin.name("phi"), "phi", 35, -3.5,  3.5);

  NameScheme iso("iso");
  isoJet_    = new TH1F(iso.name("isoJet"  ), "isoJet"  , 80,   0.,  4.);
  isoJet5_   = new TH1F(iso.name("isoJet5" ), "isoJet5" , 80,   0.,  4.);
  isoJet10_  = new TH1F(iso.name("isoJet10"), "isoJet10", 80,   0.,  4.);
  isoJet15_  = new TH1F(iso.name("isoJet15"), "isoJet15", 80,   0.,  4.);
  isoJet20_  = new TH1F(iso.name("isoJet20"), "isoJet20", 80,   0.,  4.);
  isoJet25_  = new TH1F(iso.name("isoJet25"), "isoJet25", 80,   0.,  4.);
  isoJet30_  = new TH1F(iso.name("isoJet30"), "isoJet30", 80,   0.,  4.);
  isoJet35_  = new TH1F(iso.name("isoJet35"), "isoJet35", 80,   0.,  4.);
  isoJet40_  = new TH1F(iso.name("isoJet40"), "isoJet40", 80,   0.,  4.);
  isoTrkPt_  = new TH1F(iso.name("isoTrkPt"), "isoTrkPt", 60,  -1.,  5.);
  isoCalPt_  = new TH1F(iso.name("isoCalPt"), "isoCalPt", 40, -10., 30.);
  isoRelPt_  = new TH1F(iso.name("isoRelPt"), "isoRelPt", 44,   0., 1.1);
  isoRelComb_ = new TH1F(iso.name("isoRelComb" ), "isoRelComb" , 44,   0.,  1.1);
  isoTrkN_   = new TH1F(iso.name("isoTrkN" ), "isoTrkN" , 20,   0., 20.);
  isoEcalN_  = new TH1F(iso.name("isoEcalN"), "isoEcalN", 30,   0., 30.);
  isoHcalN_  = new TH1F(iso.name("isoHcalN"), "isoHcalN", 30,   0., 30.);
  dRTrkPt_   = new TH1F(iso.name("dRTrkPt" ), "dRTrkPt" , 20,   0.,  1.);
  dRTrkN_    = new TH1F(iso.name("dRTrkN"  ), "dRTrkN"  , 20,   0.,  1.);
  dREcalPt_  = new TH1F(iso.name("dREcalPt"), "dREcalPt", 20,   0.,  1.);
  dREcalN_   = new TH1F(iso.name("dREcalN" ), "dREcalN" , 20,   0.,  1.);
  dRHcalPt_  = new TH1F(iso.name("dRHcalPt"), "dRHcalPt", 20,   0.,  1.);
  dRHcalN_   = new TH1F(iso.name("dRHcalN" ), "dRHcalN" , 20,   0.,  1.);

  NameScheme iso2d("iso2d");
  ptVsTrkIso_ = new TH2F(iso2d.name("ptVsTrkIso"), "ptVsTrkIso", 100, 0., 100., 50,   0., 25.);
  ptVsCalIso_ = new TH2F(iso2d.name("ptVsCalIso"), "ptVsCalIso", 100, 0., 100., 50, -10., 25.);
}

/// book for full FW
void
MuonKinematic::book(edm::Service<TFileService>& fs)
{
  NameScheme kin("kin");
  en_ = fs->make<TH1F>(kin.name("en" ), "E(muon) [GeV]"    , 60,   0., 300.);
  pt_ = fs->make<TH1F>(kin.name("pt" ), "p_{T}(muon) [GeV]", 24,   0., 120.);
  eta_= fs->make<TH1F>(kin.name("eta"), "#eta(muon)"       , 34, -3.4,  3.4);
  phi_= fs->make<TH1F>(kin.name("phi"), "#phi(muon)"       , 35, -3.5,  3.5);

  NameScheme iso("iso");
  isoJet_  = fs->make<TH1F>(iso.name("isoJet"  ), "smallest dR(muon,jet_{all})"   , 35,   0., 3.5);
  isoJet5_ = fs->make<TH1F>(iso.name("isoJet5" ), "smallest dR(muon,jet_{5 GeV})" , 35,   0., 3.5);
  isoJet10_= fs->make<TH1F>(iso.name("isoJet10"), "smallest dR(muon,jet_{10 GeV})", 35,   0., 3.5);
  isoJet15_= fs->make<TH1F>(iso.name("isoJet15"), "smallest dR(muon,jet_{15 GeV})", 35,   0., 3.5);
  isoJet20_= fs->make<TH1F>(iso.name("isoJet20"), "smallest dR(muon,jet_{20 GeV})", 35,   0., 3.5);
  isoJet25_= fs->make<TH1F>(iso.name("isoJet25"), "smallest dR(muon,jet_{25 GeV})", 35,   0., 3.5);
  isoJet30_= fs->make<TH1F>(iso.name("isoJet30"), "smallest dR(muon,jet_{30 GeV})", 35,   0., 3.5);
  isoJet35_= fs->make<TH1F>(iso.name("isoJet35"), "smallest dR(muon,jet_{35 GeV})", 35,   0., 3.5);
  isoJet40_= fs->make<TH1F>(iso.name("isoJet40"), "smallest dR(muon,jet_{40 GeV})", 35,   0., 3.5);
  isoTrkPt_= fs->make<TH1F>(iso.name("isoTrkPt"), "muon trk. isol. [GeV]"         , 20,   0., 10.);
  isoCalPt_= fs->make<TH1F>(iso.name("isoCalPt"), "muon cal. isol. [GeV]"         , 20,   0., 10.);
  isoRelPt_= fs->make<TH1F>(iso.name("isoRelPt"), "rel. comb. isol."              , 44,   0., 1.1);
  isoRelComb_ =fs->make<TH1F>(iso.name("isoRelComb" ), "rel. comb. isol." , 44,   0.,  1.1);
  isoTrkN_ = fs->make<TH1F>(iso.name("isoTrkN" ), "N_{track}^{dR<0.3}"            , 20,   0., 20.);
  isoEcalN_= fs->make<TH1F>(iso.name("isoEcalN"), "N_{tower,ECAL}^{dR<0.3}"       , 30,   0., 30.);
  isoHcalN_= fs->make<TH1F>(iso.name("isoHcalN"), "N_{tower,HCAL}^{dR<0.3}"       , 30,   0., 30.);
  dRTrkPt_ = fs->make<TH1F>(iso.name("dRTrkPt" ), "dR(muon,track)"                , 20,   0.,  1.);
  dRTrkN_  = fs->make<TH1F>(iso.name("dRTrkN"  ), "dR(muon,track)"                , 20,   0.,  1.);
  dREcalPt_= fs->make<TH1F>(iso.name("dREcalPt"), "dR(muon,ECAL-tower)"           , 20,   0.,  1.);
  dREcalN_ = fs->make<TH1F>(iso.name("dREcalN" ), "dR(muon,ECAL-tower)"           , 20,   0.,  1.);
  dRHcalPt_= fs->make<TH1F>(iso.name("dRHcalPt"), "dR(muon,HCAL-tower)"           , 20,   0.,  1.);
  dRHcalN_ = fs->make<TH1F>(iso.name("dRHcalN" ), "dR(muon,HCAL-tower)"           , 20,   0.,  1.);

  NameScheme iso2d("iso2d");
  ptVsTrkIso_ = fs->make<TH2F>(iso2d.name("ptVsTrkIso"), "ptVsTrkIso", 100, 0., 100., 50,   0., 25.);
  ptVsCalIso_ = fs->make<TH2F>(iso2d.name("ptVsCalIso"), "ptVsCalIso", 100, 0., 100., 50, -10., 25.);
}

/// book for full FW with output stream
void
MuonKinematic::book(edm::Service<TFileService>& fs, ofstream& file)
{
  NameScheme kin("kin");
  en_ = fs->make<TH1F>(kin.name(file, "en" ), "E(muon) [GeV]"    , 60,   0., 300.);
  pt_ = fs->make<TH1F>(kin.name(file, "pt" ), "p_{T}(muon) [GeV]", 24,   0., 120.);
  eta_= fs->make<TH1F>(kin.name(file, "eta"), "#eta(muon)"       , 34, -3.4,  3.4);
  phi_= fs->make<TH1F>(kin.name(file, "phi"), "#phi(muon)"       , 35, -3.5,  3.5);

  NameScheme iso("iso");
  isoJet_  = fs->make<TH1F>(iso.name(file, "isoJet"  ), "smallest dR(muon,jet_{all})"   , 35,   0., 3.5);
  isoJet5_ = fs->make<TH1F>(iso.name(file, "isoJet5" ), "smallest dR(muon,jet_{5 GeV})" , 35,   0., 3.5);
  isoJet10_= fs->make<TH1F>(iso.name(file, "isoJet10"), "smallest dR(muon,jet_{10 GeV})", 35,   0., 3.5);
  isoJet15_= fs->make<TH1F>(iso.name(file, "isoJet15"), "smallest dR(muon,jet_{15 GeV})", 35,   0., 3.5);
  isoJet20_= fs->make<TH1F>(iso.name(file, "isoJet20"), "smallest dR(muon,jet_{20 GeV})", 35,   0., 3.5);
  isoJet25_= fs->make<TH1F>(iso.name(file, "isoJet25"), "smallest dR(muon,jet_{25 GeV})", 35,   0., 3.5);
  isoJet30_= fs->make<TH1F>(iso.name(file, "isoJet30"), "smallest dR(muon,jet_{30 GeV})", 35,   0., 3.5);
  isoJet35_= fs->make<TH1F>(iso.name(file, "isoJet35"), "smallest dR(muon,jet_{35 GeV})", 35,   0., 3.5);
  isoJet40_= fs->make<TH1F>(iso.name(file, "isoJet40"), "smallest dR(muon,jet_{40 GeV})", 35,   0., 3.5);
  isoTrkPt_= fs->make<TH1F>(iso.name(file, "isoTrkPt"), "muon trk. isol. [GeV]"         , 20,   0., 10.);
  isoCalPt_= fs->make<TH1F>(iso.name(file, "isoCalPt"), "muon cal. isol. [GeV]"         , 20,   0., 10.);
  isoRelPt_= fs->make<TH1F>(iso.name(file, "isoRelPt"), "rel. comb. isol."              , 44,   0., 1.1);
  isoRelComb_ =fs->make<TH1F>(iso.name("isoRelComb" ), "rel. comb. isol."		 , 44,   0.,  1.1);
  isoTrkN_ = fs->make<TH1F>(iso.name(file, "isoTrkN" ), "N_{track}^{dR<0.3}"            , 20,   0., 20.);
  isoEcalN_= fs->make<TH1F>(iso.name(file, "isoEcalN"), "N_{tower,ECAL}^{dR<0.3}"       , 30,   0., 30.);
  isoHcalN_= fs->make<TH1F>(iso.name(file, "isoHcalN"), "N_{tower,HCAL}^{dR<0.3}"       , 30,   0., 30.);
  dRTrkPt_ = fs->make<TH1F>(iso.name(file, "dRTrkPt" ), "dR(muon,track)"                , 20,   0.,  1.);
  dRTrkN_  = fs->make<TH1F>(iso.name(file, "dRTrkN"  ), "dR(muon,track)"                , 20,   0.,  1.);
  dREcalPt_= fs->make<TH1F>(iso.name(file, "dREcalPt"), "dR(muon,ECAL-tower)"           , 20,   0.,  1.);
  dREcalN_ = fs->make<TH1F>(iso.name(file, "dREcalN" ), "dR(muon,ECAL-tower)"           , 20,   0.,  1.);
  dRHcalPt_= fs->make<TH1F>(iso.name(file, "dRHcalPt"), "dR(muon,HCAL-tower)"           , 20,   0.,  1.);
  dRHcalN_ = fs->make<TH1F>(iso.name(file, "dRHcalN" ), "dR(muon,HCAL-tower)"           , 20,   0.,  1.);

  NameScheme iso2d("iso2d");
  ptVsTrkIso_ = fs->make<TH2F>(iso2d.name(file, "ptVsTrkIso"), "ptVsTrkIso", 100, 0., 100., 50,   0., 25.);
  ptVsCalIso_ = fs->make<TH2F>(iso2d.name(file, "ptVsCalIso"), "ptVsCalIso", 100, 0., 100., 50, -10., 25.);
}

/// write to file and free allocated space for FWLite
void
MuonKinematic::write(TFile& file, const char* directory)
{
  file.cd( directory );
  en_        ->Write( );
  pt_        ->Write( );
  eta_       ->Write( );
  phi_       ->Write( );
  isoJet_    ->Write( );
  isoJet5_   ->Write( );
  isoJet10_  ->Write( );
  isoJet15_  ->Write( );
  isoJet20_  ->Write( );
  isoJet25_  ->Write( );
  isoJet30_  ->Write( );
  isoJet35_  ->Write( );
  isoJet40_  ->Write( );
  isoTrkPt_  ->Write( );
  isoCalPt_  ->Write( );
  isoRelPt_  ->Write( );
  isoRelComb_->Write();
  isoTrkN_   ->Write( );
  isoEcalN_  ->Write( );
  isoHcalN_  ->Write( );
  dRTrkPt_   ->Write( );
  dRTrkN_    ->Write( );
  dREcalPt_  ->Write( );
  dREcalN_   ->Write( );
  dRHcalPt_  ->Write( );
  dRHcalN_   ->Write( );
  ptVsTrkIso_->Write( );
  ptVsCalIso_->Write( );
}
