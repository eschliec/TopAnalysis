
#include "TopAnalysis/TopAnalyzer/plugins/TriggerAnalyzer.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace std;
using namespace edm;
using namespace trigger;

TriggerAnalyzer::TriggerAnalyzer( const edm::ParameterSet& ps ) {

  triggerResults_ = ps.getParameter<edm::InputTag>("TriggerResults");
  triggerEvent_   = ps.getParameter<edm::InputTag>("TriggerEvent");
  triggerFilter_  = ps.getParameter<edm::InputTag>("TriggerFilter");

  hltPaths_       = ps.getParameter<vector<string> >("hltPaths");
  hltPaths_sig_   = ps.getParameter<vector<string> >("hltPaths_sig");
  hltPaths_trig_  = ps.getParameter<vector<string> >("hltPaths_trig");

  vertex_       = ps.getParameter<edm::InputTag>("vertexCollection");
  vertex_X_cut_ = ps.getParameter<double>("vertex_X_cut");
  vertex_Y_cut_ = ps.getParameter<double>("vertex_Y_cut");
  vertex_Z_cut_ = ps.getParameter<double>("vertex_Z_cut");

  muons_        = ps.getParameter<edm::InputTag>("muonCollection");
  muon_pT_cut_  = ps.getParameter<double>("muon_pT_cut");
  muon_eta_cut_ = ps.getParameter<double>("muon_eta_cut");
  muon_iso_cut_ = ps.getParameter<double>("muon_iso_cut");

  MassWindow_up_   = ps.getParameter<double>("MassWindow_up");
  MassWindow_down_ = ps.getParameter<double>("MassWindow_down");

}

TriggerAnalyzer::~TriggerAnalyzer() {

}

void TriggerAnalyzer::beginJob() {

  Service<TFileService> fs;

  if( !fs )  throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );

  Trigs      = fs->make<TH1F>("Trigs",      "Fired muon/electron triggers", 15,  0., 15.);
  TriggerEfficiencies = fs->make<TH1F>("TriggerEfficiencies", "HL Trigger Efficiencies",      10,  0., 10.);

  TriggerEfficiencies->SetTitle("HL Trigger Efficiencies #epsilon_{signal} = #frac{[signal] && [control]}{[control]}");

  TriggerEfficiencies_sig  = fs->make<TH1F>("TriggerEfficiencies_sig",  "HL Trigger Signal && Control Counts", 10,  0., 10.);
  TriggerEfficiencies_trig = fs->make<TH1F>("TriggerEfficiencies_trig", "HL Trigger Control Counts",           10,  0., 10.);

  const int N_TriggerPaths = hltPaths_.size();
  const int N_SignalPaths  = hltPaths_sig_.size();

  for( int i = 0; i < N_TriggerPaths; i++ ) {

    const string &label = hltPaths_[i];
    Trigs->GetXaxis()->SetBinLabel( i+1,label.c_str() );

  }

  for( int j = 0; j < N_SignalPaths; ++j ) {

    const string &label_eff = "#frac{["+hltPaths_sig_[j]+"]}{vs. ["+hltPaths_trig_[j]+"]}";
    const string &label_sig = hltPaths_sig_[j]+" && "+hltPaths_trig_[j];
    TriggerEfficiencies->GetXaxis()->SetBinLabel(     j+1, label_eff.c_str() );
    TriggerEfficiencies_sig->GetXaxis()->SetBinLabel( j+1, label_sig.c_str() );

  }

  Correlations = fs->make<TH2D>("Correlations", "Trigger Correlations", 15, 0., 15., 15, 0., 15.);

  for( int j = 0; j < N_TriggerPaths; j++ ) {

    const string &label = hltPaths_[j];
    Correlations->GetXaxis()->SetBinLabel( j+1,label.c_str() );
    Correlations->GetYaxis()->SetBinLabel( j+1,label.c_str() );

  }

  const int nbins_Pt = 7;

  float bins_Pt[nbins_Pt+1] = { 0., 5., 10., 20., 40., 60., 80., 100. };

  MuonEfficiency_pT      = fs->make<TH1F>("MuonEfficiency_pT",      "Muon Efficiency P_{T}",           nbins_Pt, &bins_Pt[0]);
  MuonEfficiency_pT_sig  = fs->make<TH1F>("MuonEfficiency_Pt_sig",  "P^{#mu}_{T} (signal triggered)",  nbins_Pt, &bins_Pt[0]);
  MuonEfficiency_pT_trig = fs->make<TH1F>("MuonEfficiency_Pt_trig", "P^{#mu}_{T} (control triggered)", nbins_Pt, &bins_Pt[0]);

  const int nbins_eta = 8;

  float bins_eta[nbins_eta+1] = { -2.5, -2.1, -1.3, -0.8, 0., 0.8, 1.3, 2.1, 2.5 };

  MuonEfficiency_eta      = fs->make<TH1F>("MuonEfficiency_eta",      "Muon Efficiency  #eta",           nbins_eta, &bins_eta[0]);
  MuonEfficiency_eta_sig  = fs->make<TH1F>("MuonEfficiency_Eta_sig",  "#eta_{muon} (signal triggered)",  nbins_eta, &bins_eta[0]);
  MuonEfficiency_eta_trig = fs->make<TH1F>("MuonEfficiency_Eta_trig", "#eta_{muon} (control triggered)", nbins_eta, &bins_eta[0]);

  const int nbins_phi = 9;

  float bins_phi[nbins_phi+1] = { -3.5, -3.2, -2.6, -1.56, -0.52, 0.52, 1.56, 2.6, 3.2, 3.5 };

  MuonEfficiency_phi      = fs->make<TH1F>("MuonEfficiency_phi",      "Muon Efficiency  #phi",           nbins_phi, &bins_phi[0]);
  MuonEfficiency_phi_sig  = fs->make<TH1F>("MuonEfficiency_Phi_sig",  "#phi_{muon} (signal triggered)",  nbins_phi, &bins_phi[0]);
  MuonEfficiency_phi_trig = fs->make<TH1F>("MuonEfficiency_Phi_trig", "#phi_{muon} (control triggered)", nbins_phi, &bins_phi[0]);

  MuonEfficiency_eta_phi      = fs->make<TH2F>("MuonEfficiency_eta_phi", "Muon Efficiency  #eta vs. #phi",
					       nbins_eta, &bins_eta[0], nbins_phi, &bins_phi[0]);

  MuonEfficiency_eta_phi_sig  = fs->make<TH2F>("MuonEfficiency_Eta_Phi_sig", "#eta_{muon} vs. #phi_{muon} (signal triggered)",
					       nbins_eta, &bins_eta[0], nbins_phi, &bins_phi[0]);

  MuonEfficiency_eta_phi_trig = fs->make<TH2F>("MuonEfficiency_Eta_Phi_trig", "#eta_{muon} vs. #phi_{muon} (control triggered)",
					       nbins_eta, &bins_eta[0], nbins_phi, &bins_phi[0]);

  MuonEfficiency_pT_eta       = fs->make<TH2F>("MuonEfficiency_pT_eta","Muon Efficiency  P_{T} vs. #eta",
					       nbins_Pt, &bins_Pt[0], nbins_eta, &bins_eta[0]);

  MuonEfficiency_pT_eta_sig   = fs->make<TH2F>("MuonEfficiency_PT_Eta_sig","P_{T}^{muon} vs. #eta_{muon}  (signal triggered)",
					       nbins_Pt, &bins_Pt[0], nbins_eta, &bins_eta[0]);

  MuonEfficiency_pT_eta_trig  = fs->make<TH2F>("MuonEfficiency_PT_Eta_trig", "P_{T}^{muon} vs. #eta_{muon}  (control triggered)",
					       nbins_Pt, &bins_Pt[0], nbins_eta, &bins_eta[0]);

  VxVyVz_PV = fs->make<TH3F>("VxVyVz_PV", "Primary Vertex 3D-positon", 40, -20., 20., 40, -1., 1., 40, -1., 1.);

  NMuons        = fs->make<TH1F>("Nmuons",        "Number of muons",           10, -0.5,  9.5);
  NMuons_iso    = fs->make<TH1F>("Nmuons_Iso",    "Number of isolated muons",  10, -0.5,  9.5);
  NMuons_charge = fs->make<TH1F>("Nmuons_Charge", "N_{muons} * Q(#mu)",        19, -10.,  10.);
  NTracks       = fs->make<TH1F>("Ntracks",       "Number of tracks",          50, -0.5, 49.5);
  VxVy_muons    = fs->make<TH2F>("VxVy_muons",    "Vertex x-y-positon (muon)", 40,  -1.,   1., 40, -1., 1.);
  Vz_muons      = fs->make<TH1F>("Vz_muons",      "Vertex z-positon (muon)",   40, -20.,  20.);
  PtMuons       = fs->make<TH1F>("PtMuon",        "P^{#mu}_{T}",               50,   0., 100.);
  EtaMuons      = fs->make<TH1F>("EtaMuon",       "#eta_{muon}",               50,  -5.,   5.);
  PhiMuons      = fs->make<TH1F>("PhiMuon",       "#phi_{muon}",               40,  -4.,   4.);

  DeltaEtaMuonsRC   = fs->make<TH1F>("DeltaEtaMuonsRC",    "#Delta #eta of muon pair (RC)", 50, -5.,   5. );
  DeltaPhiMuonsRC   = fs->make<TH1F>("DeltaPhiMuonsRC",    "#Delta #phi of muon pair (RC)", 50, -5.,   5. );
  DeltaEtaMuonsWC   = fs->make<TH1F>("DeltaEtaMuonsWC",    "#Delta #eta of muon pair (WC)", 50, -5.,   5. );
  DeltaPhiMuonsWC   = fs->make<TH1F>("DeltaPhiMuonsWC",    "#Delta #phi of muon pair (WC)", 50, -5.,   5. );
  CombRelIso03      = fs->make<TH1F>("MuIso_CombRelIso03", "Muon CombRelIso dR=03",         50,  0.,   1. );
  PixelHits_muons   = fs->make<TH1F>("NPixelHits_muons",   "Number of pixel hits",          50, -0.5, 49.5);
  TrackerHits_muons = fs->make<TH1F>("NTrackerHits_muons", "Number of hits in the tracker", 50, -0.5, 49.5);

  DeltaR_Trig   = fs->make<TH1F>("DeltaRTrigger", "#Delta R of trigger muon pair",       50, 0., 5.);
  DeltaR_Reco   = fs->make<TH1F>("DeltaRReco",    "#Delta R of RECO muon pair",          50, 0., 5.);
  DeltaR_Match  = fs->make<TH1F>("DeltaRMatch",   "#Delta R of matched muon pair",       50, 0., 5.);
  Trigger_Match = fs->make<TH1F>("TriggerMatch",  "Number of Trigger-RECO assignements",  6, 0., 6.);

  const int nbins_Pt_Log = 20;

  double logmin = 0.;
  double logmax = 3.;  // 10^(3.)=1000

  float bins_Pt_Log[nbins_Pt_Log+1];

  for (int i = 0; i <= nbins_Pt_Log; i++) {
    double log = logmin + (logmax-logmin)*i/nbins_Pt_Log;
    bins_Pt_Log[i] = std::pow(10.0, log);
  }

  PtMuons_LOGX  = fs->make<TH1F>("Pt_muon_LOGX", "P^{#mu}_{T}", nbins_Pt_Log, &bins_Pt_Log[0]);

  MuonEfficiency_pT_LOGX      = fs->make<TH1F>("MuonEfficiency_pT_LOGX",      "Muon Efficiency P_{T}",           nbins_Pt_Log, &bins_Pt_Log[0]);
  MuonEfficiency_pT_LOGX_sig  = fs->make<TH1F>("MuonEfficiency_Pt_LOGX_sig",  "P^{#mu}_{T} (signal triggered)",  nbins_Pt_Log, &bins_Pt_Log[0]);
  MuonEfficiency_pT_LOGX_trig = fs->make<TH1F>("MuonEfficiency_Pt_LOGX_trig", "P^{#mu}_{T} (control triggered)", nbins_Pt_Log, &bins_Pt_Log[0]);

  const int nbins_mass = 50;

  float bins_mass[nbins_mass+1];

  for (int i = 0; i <= nbins_mass; i++) {
    double log = logmin + (logmax-logmin)*i/nbins_mass;
    bins_mass[i] = std::pow(10.0, log);
  }

  DiMuonMassRC       = fs->make<TH1F>("DiMuonMassRC",      "Invariant Dimuon Mass (Right Charge)", 50, 0., 200.);
  DiMuonMassWC       = fs->make<TH1F>("DiMuonMassWC",      "Invariant Dimuon Mass (Wrong Charge)", 50, 0., 200.);

  DiMuonMassRC_LOGX  = fs->make<TH1F>("DiMuonMassRC_LOGX", "Invariant Dimuon Mass (Right Charge)", nbins_mass, &bins_mass[0]);
  DiMuonMassWC_LOGX  = fs->make<TH1F>("DiMuonMassWC_LOGX", "Invariant Dimuon Mass (Wrong Charge)", nbins_mass, &bins_mass[0]);

}

void TriggerAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup&) {

  // ------------------------
  //  Global Event Variables
  // ------------------------

  vector<string> hltPaths = hltPaths_;

  vector<reco::Particle> Triggered_muons;
  reco::MuonCollection    Isolated_muons;
  reco::MuonCollection     Matched_muons;

  const int N_TriggerPaths = hltPaths.size();
  const int N_SignalPaths  = hltPaths_sig_.size();
  const int N_ControlPaths = hltPaths_trig_.size();

  bool Fired_Signal_Trigger[10]  = {false};
  bool Fired_Control_Trigger[10] = {false};

  double DilepMass = 0.;

  double deltaR_Trig  = 999.;
  double deltaR_Reco  = 999.;
  double deltaR_Match = 999.;

  int N_iso_mu = 0;

  double vertex_X = 100.;
  double vertex_Y = 100.;
  double vertex_Z = 100.;

  // -------------------------
  //  Analyze Trigger Results
  // -------------------------

  edm::Handle<TriggerResults> trigResults;
  evt.getByLabel(triggerResults_, trigResults);

  if( trigResults.failedToGet() ) {

    //    cout << endl << "-----------------------------" << endl;
    //    cout << "--- NO TRIGGER RESULTS !! ---" << endl;
    //    cout << "-----------------------------" << endl << endl;

  }

  if( !trigResults.failedToGet() ) {

    const edm::TriggerNames & trigName = evt.triggerNames(*trigResults);

    for( unsigned int i_Trig = 0; i_Trig < trigResults->size(); ++i_Trig ) {

      //      cout << "Trigger Name : " << trigName.triggerName(i_Trig) << endl;

      if(trigResults->accept(i_Trig)) {

	// Check for all trigger paths

	for( int i = 0; i < N_TriggerPaths; i++ ) {

	  if( trigName.triggerName(i_Trig) == hltPaths[i] )  Trigs->Fill(i);

	  // Check for trigger path correlations

	  for( unsigned int j_Trig = 0; j_Trig < trigResults->size(); ++j_Trig ) {

	    if(trigResults->accept(j_Trig)) {

	      for( int j = 0; j < N_TriggerPaths; j++ ) {

		if( trigName.triggerName(j_Trig) == hltPaths[j] )  Correlations->Fill(i,j);

	      }

	    }

	  }

	  // Check for signal & control trigger paths

	  for( int j = 0; j < N_SignalPaths; ++j ) {

	    if( trigName.triggerName(i_Trig) == hltPaths_sig_[j]  )  Fired_Signal_Trigger[j]  = true;

	  }

	  for( int k = 0; k < N_ControlPaths; ++k ) {

	    if( trigName.triggerName(i_Trig) == hltPaths_trig_[k] )  Fired_Control_Trigger[k] = true;

	  }

	}

      }

    }

  }

  // -----------------------
  //  Analyze Trigger Event
  // -----------------------

  edm::Handle<TriggerEvent> triggerEvent;
  evt.getByLabel(triggerEvent_, triggerEvent);

  if( triggerEvent.failedToGet() ) {

    cout << endl << "---------------------------" << endl;
    cout << "--- NO TRIGGER EVENT !! ---" << endl;
    cout << "---------------------------" << endl << endl;

  }

  if( !triggerEvent.failedToGet() ) {

    size_t filterIndex = triggerEvent->filterIndex( triggerFilter_ );
    TriggerObjectCollection triggerObjects = triggerEvent->getObjects();

    if( filterIndex >= triggerEvent->sizeFilters() ) {

      //      cout << endl << "------------------------------" << endl;
      //      cout << "--- NO FILTERED OBJECTS !! ---" << endl;
      //      cout << "------------------------------" << endl << endl;

      //      cout << "Trigger Filter Size  : " << triggerEvent->sizeFilters() << endl;
      //      cout << "Trigger Filter Index : " << filterIndex                 << endl;

    }

    if( filterIndex < triggerEvent->sizeFilters() ) {

      //      cout << "Trigger Filter Size  : " << triggerEvent->sizeFilters() << endl;
      //      cout << "Trigger Filter Index : " << filterIndex                 << endl;

      const Keys & keys = triggerEvent->filterKeys( filterIndex );

      int N_mu = 0;

      for( size_t j = 0; j < keys.size(); j++ ) {

	TriggerObject foundObject = triggerObjects[keys[j]];

	if(      foundObject.pt()   < muon_pT_cut_  )      continue;
	if( abs( foundObject.eta()) > muon_eta_cut_ )      continue;
	if( abs( foundObject.particle().pdgId() ) != 13 )  continue;

	++N_mu;
	Triggered_muons.push_back( foundObject.particle() );

      }

      if( Triggered_muons.size() == 2 ) {

	reco::Particle mu1 = Triggered_muons.at(0);
	reco::Particle mu2 = Triggered_muons.at(1);

	deltaR_Trig = deltaR( mu1.eta(), mu1.phi(), mu2.eta(), mu2.phi() );
	DeltaR_Trig->Fill(deltaR_Trig);

      }

    }

  }

  // ------------------------
  //  Analyze Primary Vertex
  // ------------------------

  edm::Handle<reco::VertexCollection> vertexs;
  evt.getByLabel(vertex_, vertexs);

  if( vertexs.failedToGet() ) {

    //      cout << endl << "----------------------------" << endl;
    //      cout << "--- NO PRIMARY VERTEX !! ---" << endl;
    //      cout << "----------------------------" << endl << endl;

  }

  if( !vertexs.failedToGet() ) {

    reco::Vertex primaryVertex = vertexs->front();

    int numberTracks = primaryVertex.tracksSize();
    //      double ndof      = primaryVertex.ndof();
    bool fake        = primaryVertex.isFake();

    NTracks->Fill(numberTracks);

    if( !fake && numberTracks > 3 ) {

      vertex_X = primaryVertex.x();
      vertex_Y = primaryVertex.y();
      vertex_Z = primaryVertex.z();

      VxVyVz_PV->Fill( vertex_Z, vertex_X, vertex_Y );

    }

  }

  // -------------------
  //  Analyze PAT Muons
  // -------------------

  //  edm::Handle<reco::MuonCollection> muons;
  edm::Handle< vector<pat::Muon> > muons;
  evt.getByLabel(muons_, muons);

  //  reco::MuonCollection::const_iterator muon;
  vector<pat::Muon>::const_iterator muon;

  if( muons.failedToGet() ) {

    cout << endl << "-----------------------" << endl;
    cout << "--- NO PAT MUONS !! ---" << endl;
    cout << "-----------------------" << endl << endl;

  }

  if( !muons.failedToGet() ) {

    NMuons->Fill( muons->size() );

    for(muon = muons->begin(); muon!= muons->end(); ++muon) {

      float N_muons = muons->size();
      float Q_muon  = muon->charge();

      NMuons_charge->Fill(N_muons*Q_muon);

      double track_X = 100.;
      double track_Y = 100.;
      double track_Z = 100.;

      double N_PixelHits   = 0.;
      double N_TrackerHits = 0.;

      if( muon->isGlobalMuon() && muon->isTrackerMuon() ) {

	reco::TrackRef track = muon->globalTrack();

	track_X       = track->vx();
	track_Y       = track->vy();
	track_Z       = track->vz();
	N_PixelHits   = track->hitPattern().numberOfValidPixelHits();
	N_TrackerHits = track->hitPattern().numberOfValidTrackerHits();

	VxVy_muons->Fill(track_X, track_Y);
	Vz_muons->Fill(track_Z);
	PixelHits_muons->Fill(N_PixelHits);
	TrackerHits_muons->Fill(N_TrackerHits);

      }

      // Vertex and kinematic cuts

      if(          track_X > vertex_X_cut_ )  continue;
      if(          track_Y > vertex_Y_cut_ )  continue;
      if(          track_Z > vertex_Z_cut_ )  continue;
      if(      N_PixelHits <  1.           )  continue;
      if(    N_TrackerHits < 11.           )  continue;
      if(     muon->pt()   < muon_pT_cut_  )  continue;
      if( abs(muon->eta()) > muon_eta_cut_ )  continue;

      reco::MuonIsolation muIso03 = muon->isolationR03();

      double muonCombRelIso = 1.;

      if ( muon->pt() != 0. )
	muonCombRelIso = ( muIso03.emEt + muIso03.hadEt + muIso03.hoEt + muIso03.sumPt ) / muon->pt();

      CombRelIso03->Fill( muonCombRelIso );

      if( muonCombRelIso < muon_iso_cut_ ) {

	++N_iso_mu;
	Isolated_muons.push_back(*muon);

      }

    }

    NMuons_iso->Fill(N_iso_mu);

    //    if( Isolated_muons.size() > 1 && Fired_Control_Trigger[0] ) {
    if( Isolated_muons.size() > 1 ) {

      // Vertex cut

      if( vertex_X < vertex_X_cut_ && vertex_Y < vertex_Y_cut_ && vertex_Z < vertex_Z_cut_ ) {

	for( int i = 0; i < (static_cast<int>(Isolated_muons.size()) - 1); ++i ) {

	  for( int j = i+1; j < static_cast<int>(Isolated_muons.size()); ++j ) {

	    reco::MuonCollection::const_reference mu1 = Isolated_muons.at(i);
	    reco::MuonCollection::const_reference mu2 = Isolated_muons.at(j);

	    DilepMass = sqrt( (mu1.energy()+mu2.energy())*(mu1.energy()+mu2.energy())
			      - (mu1.px()+mu2.px())*(mu1.px()+mu2.px())
			      - (mu1.py()+mu2.py())*(mu1.py()+mu2.py())
			      - (mu1.pz()+mu2.pz())*(mu1.pz()+mu2.pz())
			      );

	    if( DilepMass < 1. ) {

	      if( i > 0 ) {

		Isolated_muons.erase(Isolated_muons.begin()+i);
		--i;

	      }

	      continue;

	    }

	    // Opposite muon charges -> Right Charge (RC)

	    if( mu1.charge()*mu2.charge() < 0. ) {

	      DiMuonMassRC->Fill( DilepMass );
	      DiMuonMassRC_LOGX->Fill( DilepMass );

	      if( DilepMass > MassWindow_down_ && DilepMass < MassWindow_up_ ) {

		PtMuons->Fill(  mu1.pt()  );
		PtMuons->Fill(  mu2.pt()  );
		PtMuons_LOGX->Fill( mu1.pt() );
		PtMuons_LOGX->Fill( mu2.pt() );
		EtaMuons->Fill( mu1.eta() );
		EtaMuons->Fill( mu2.eta() );
		PhiMuons->Fill( mu1.phi() );
		PhiMuons->Fill( mu2.phi() );

		DeltaEtaMuonsRC->Fill(mu1.eta()-mu2.eta());
		DeltaPhiMuonsRC->Fill( deltaPhi(mu1.phi(),mu2.phi()) );

		// Determinating relative trigger efficiencies

		for( int k = 0; k < N_SignalPaths; ++k ) {

		  if( Fired_Signal_Trigger[k] && Fired_Control_Trigger[k] )  TriggerEfficiencies_sig->Fill(k);

		  if( Fired_Control_Trigger[k] )  TriggerEfficiencies_trig->Fill(k);

		}

		// Trigger object matching

		int N_Match = 0;
		double   DR = 0.1;

		if( Isolated_muons.size() == 2 && Triggered_muons.size() > 0 ) {

		  deltaR_Reco = deltaR( mu1.eta(), mu1.phi(), mu2.eta(), mu2.phi() );
		  DeltaR_Reco->Fill(deltaR_Reco);

		  if( deltaR_Reco > 2.*DR && deltaR_Trig > 2.*DR ) {

		    for( int i = 0; i < static_cast<int>(Isolated_muons.size()); ++i ) {

		      for( int j = 0; j < static_cast<int>(Triggered_muons.size()); ++j ) {

			reco::Particle & Trigger_mu = Triggered_muons.at(j);
			reco::Muon     &    Reco_mu =  Isolated_muons.at(i);

			deltaR_Match = deltaR( Trigger_mu.eta(), Trigger_mu.phi(), Reco_mu.eta(), Reco_mu.phi() );
			DeltaR_Match->Fill(deltaR_Match);

			if( deltaR_Match < DR) {

			  ++N_Match;
			  Matched_muons.push_back(Reco_mu);

			}

		      }

		    }

		  }

		  Trigger_Match->Fill( N_Match );

		}


		// Muon Tag & Probe Efficiency

		if( Matched_muons.size() == 1 ) {

		  reco::MuonCollection::const_reference matched_mu1 = Matched_muons.at(0);

		  MuonEfficiency_pT_trig->Fill( matched_mu1.pt() );
		  MuonEfficiency_pT_LOGX_trig->Fill( matched_mu1.pt() );
		  MuonEfficiency_eta_trig->Fill(matched_mu1.eta());
		  MuonEfficiency_phi_trig->Fill(matched_mu1.phi());

		  MuonEfficiency_eta_phi_trig->Fill( matched_mu1.eta(), matched_mu1.phi());
		  MuonEfficiency_pT_eta_trig->Fill(  matched_mu1.pt(),  matched_mu1.eta());

		}

		if( Matched_muons.size() == 2 ) {

		  reco::MuonCollection::const_reference matched_mu1 = Matched_muons.at(0);
		  reco::MuonCollection::const_reference matched_mu2 = Matched_muons.at(1);

		  MuonEfficiency_pT_trig->Fill( matched_mu1.pt() );
		  MuonEfficiency_pT_trig->Fill( matched_mu2.pt() );
		  MuonEfficiency_pT_LOGX_trig->Fill( matched_mu1.pt() );
		  MuonEfficiency_pT_LOGX_trig->Fill( matched_mu2.pt() );
		  MuonEfficiency_eta_trig->Fill(matched_mu1.eta());
		  MuonEfficiency_eta_trig->Fill(matched_mu2.eta());
		  MuonEfficiency_phi_trig->Fill(matched_mu1.phi());
		  MuonEfficiency_phi_trig->Fill(matched_mu2.phi());

		  MuonEfficiency_eta_phi_trig->Fill( matched_mu1.eta(), matched_mu1.phi());
		  MuonEfficiency_eta_phi_trig->Fill( matched_mu2.eta(), matched_mu2.phi());
		  MuonEfficiency_pT_eta_trig->Fill(  matched_mu1.pt(),  matched_mu1.eta());
		  MuonEfficiency_pT_eta_trig->Fill(  matched_mu2.pt(),  matched_mu2.eta());

		  MuonEfficiency_pT_sig->Fill( matched_mu1.pt() );
		  MuonEfficiency_pT_sig->Fill( matched_mu2.pt() );
		  MuonEfficiency_pT_LOGX_sig->Fill( matched_mu1.pt() );
		  MuonEfficiency_pT_LOGX_sig->Fill( matched_mu2.pt() );
		  MuonEfficiency_eta_sig->Fill(matched_mu1.eta());
		  MuonEfficiency_eta_sig->Fill(matched_mu2.eta());
		  MuonEfficiency_phi_sig->Fill(matched_mu1.phi());
		  MuonEfficiency_phi_sig->Fill(matched_mu2.phi());

		  MuonEfficiency_eta_phi_sig->Fill( matched_mu1.eta(), matched_mu1.phi());
		  MuonEfficiency_eta_phi_sig->Fill( matched_mu2.eta(), matched_mu2.phi());
		  MuonEfficiency_pT_eta_sig->Fill(  matched_mu1.pt(),  matched_mu1.eta());
		  MuonEfficiency_pT_eta_sig->Fill(  matched_mu2.pt(),  matched_mu2.eta());

		}

	      }

	    }

	    // Same muon charges -> Wrong Charge (WC)

	    if( mu1.charge()*mu2.charge() > 0. ) {

	      DiMuonMassWC->Fill( DilepMass );
	      DiMuonMassWC_LOGX->Fill( DilepMass );

	      if( DilepMass > MassWindow_down_ && DilepMass < MassWindow_up_ ) {

		DeltaEtaMuonsWC->Fill( mu1.eta()-mu2.eta() );
		DeltaPhiMuonsWC->Fill( deltaPhi(mu1.phi(),mu2.phi()) );

	      }

	    }

	  }

	}

      }

    }

  }

}

void TriggerAnalyzer::endJob() {

}

DEFINE_FWK_MODULE( TriggerAnalyzer );
