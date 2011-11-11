#include "TLorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TopAnalysis/TopAnalyzer/plugins/FullLepKinAnalyzer.h"
#include "TopAnalysis/TopUtils/interface/NameScheme.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "TopAnalysis/TopAnalyzer/interface/DileptonEventWeight.h"

/// default constructor
FullLepKinAnalyzer::FullLepKinAnalyzer(const edm::ParameterSet& cfg):
  isSignalMC_      (cfg.getParameter<bool>         ("isSignalMC"     )),
  useLeadingJets_  (cfg.getParameter<bool>         ("useLeadingJets" )),   
  useBtagging_     (cfg.getParameter<bool>         ("useBtagging"    )),  
  bAlgo_           (cfg.getParameter<std::string>  ("bAlgorithm"     )),      
  FullLepEvt_      (cfg.getParameter<edm::InputTag>("FullLepEvent"   )),
  hypoKey_         (cfg.getParameter<edm::InputTag>("hypoKey"        )),
  jets_            (cfg.getParameter<edm::InputTag>("jets"           )),
  puWeight_        (cfg.getParameter<edm::InputTag>("weightPU"       )),
  lepSfWeight_     (cfg.getParameter<edm::InputTag>("weightLepSF"    ))  
{
}

/// default destructor
FullLepKinAnalyzer::~FullLepKinAnalyzer()
{
}

/// everything that has to be done before the event loop: books histograms
void
FullLepKinAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;
  if( !fs ) throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );

  bookKinHistos      (fs);
  bookQualityHistos  (fs);
  bookRecCorrelHistos(fs);
  if(isSignalMC_){
    bookGenHistos (fs);
    bookPullHistos(fs);
    book2DHistos  (fs);
  }
  if(useLeadingJets_ && useBtagging_){
    edm::LogWarning ( "KinAnalyzerSetting" ) << "Both 'useLeadingJets' and 'useBtagging' are set which is probably not intended.";
  }
}

/// everything that has to be done during the event loop: filling plots
void
FullLepKinAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{
  double weight = getDileptonEventWeight(evt, puWeight_, lepSfWeight_);
  edm::Handle<TtFullLeptonicEvent> FullLepEvt;
  evt.getByLabel(FullLepEvt_, FullLepEvt);

  edm::Handle<int> hypoKeyHandle;
  evt.getByLabel(hypoKey_, hypoKeyHandle);
  TtEvent::HypoClassKey& hypoKey = (TtEvent::HypoClassKey&) *hypoKeyHandle;

  // -----------------------
  // check if hypothesis is valid in this event
  // and if it satisfies some quality criteria
  // -----------------------

  if( !FullLepEvt->isHypoAvailable(hypoKey) ){
    edm::LogInfo ( "NonValidHyp" ) << "Hypothesis not available for this event";
    return;
  }

  edm::Handle<edm::View< pat::Jet > > jets;
  evt.getByLabel(jets_, jets);

  int nTCHEL  = 0;
  int nTCHEM  = 0;
  int nSSVHEM = 0;

  for(edm::View<pat::Jet>::const_iterator jet = jets->begin();jet != jets->end(); ++jet) {
    if(jet->bDiscriminator("trackCountingHighEffBJetTags")>1.7          && nTCHEL<2 ) ++nTCHEL;
    if(jet->bDiscriminator("trackCountingHighEffBJetTags")>3.3          && nTCHEM<2 ) ++nTCHEM;
    if(jet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags")>1.74 && nSSVHEM<2) ++nSSVHEM;
  }

  int y = 0;
  if( FullLepEvt->isHypoValid(hypoKey) ) y = 1;

  kinTCHELcorrelation_ ->Fill(nTCHEL, y, weight);
  kinTCHEMcorrelation_ ->Fill(nTCHEM, y, weight);
  kinSSVHEMcorrelation_->Fill(nSSVHEM,y, weight);

  if( !FullLepEvt->isHypoValid(hypoKey) ){
    edm::LogInfo ( "NonValidHyp" ) << "Hypothesis not valid for this event";
    TopKin_[4]->SetBinContent(1,TopKin_[4]->GetBinContent(1)+1);
    return;  // return if any of the hypotheses is not valid
  }

  // if more than 1 solution is found they can be accessed via the index cmb
  // by default the first solution with the highest probability is used
  int cmb = 0;
  
  // use only the hypotheses which use the two leading jets
  if(useLeadingJets_){
    bool foundSolution = false;
    for(size_t i=0;i<FullLepEvt->numberOfAvailableHypos(hypoKey);++i){
      int idx1 = FullLepEvt->jetLeptonCombination(hypoKey,i)[0];
      int idx2 = FullLepEvt->jetLeptonCombination(hypoKey,i)[1];      
    
      if((idx1==0 && idx2==1) || (idx1==1 && idx2==0)){
        cmb = i;
	foundSolution = true;
        break;
      }
    }
    if(!foundSolution){
      TopKin_[4]->SetBinContent(1,TopKin_[4]->GetBinContent(1)+1);
      return;
    }
  }else if(useBtagging_){  // find best solution taking into account the b-tagging discriminators
  
    std::vector<int> bidcs;
    int idx=0;
    // find indices of b-tagged jets
    for(edm::View<pat::Jet>::const_iterator jet = jets->begin();jet != jets->end(); ++jet) {
      if(jet->bDiscriminator(bAlgo_) > bCut_) { 
        bidcs.push_back(idx);
      }
      idx++;
    }
    
    int btagsinhypo;
    bool foundOneTagSolution = false;
    for(size_t i=0;i<FullLepEvt->numberOfAvailableHypos(hypoKey);++i){
      btagsinhypo = 0;
      for(size_t j=0; j<bidcs.size(); ++j){
        if(FullLepEvt->jetLeptonCombination(hypoKey,i)[0]==bidcs[j]) btagsinhypo++;
        if(FullLepEvt->jetLeptonCombination(hypoKey,i)[1]==bidcs[j]) btagsinhypo++;		
      }
      if(btagsinhypo==2){ // stop if hypothesis has two b-jets
        cmb = i;
	break;
      }	else if(btagsinhypo==1 && !foundOneTagSolution){ // if one b-tag in hypothesis store index but go on and look for solution with 2 tags
        cmb = i;
        foundOneTagSolution = true;
      }
    }        
  }


  // -----------------------
  // fill histos related to quality of the TtFullLeptonicEvent
  // -----------------------
  bJetIdcs_   ->Fill(FullLepEvt->jetLeptonCombination(hypoKey)[0], weight );
  bBarJetIdcs_->Fill(FullLepEvt->jetLeptonCombination(hypoKey)[1], weight );   
  
  if( hypoKey==6 ) {
    deltaM_->Fill(FullLepEvt->top(hypoKey)->mass()-FullLepEvt->topBar(hypoKey)->mass(), weight); 
    kinSolWeight_->Fill( FullLepEvt->solWeight(), weight );
  }

  fillQualityHistos(*FullLepEvt, hypoKey, weight);


  // -----------------------
  // fill histos for basic kinematic variables
  // -----------------------
  const reco::Candidate* Top    = FullLepEvt->top(hypoKey, cmb);
  const reco::Candidate* Wplus  = FullLepEvt->wPlus(hypoKey, cmb);
  const reco::Candidate* B      = FullLepEvt->b(hypoKey, cmb);
  const reco::Candidate* LepBar = FullLepEvt->leptonBar(hypoKey, cmb);
  const reco::Candidate* Nu     = FullLepEvt->neutrino(hypoKey, cmb);

  const reco::Candidate* TopBar = FullLepEvt->topBar(hypoKey, cmb);
  const reco::Candidate* Wminus = FullLepEvt->wMinus(hypoKey, cmb);
  const reco::Candidate* BBar   = FullLepEvt->bBar(hypoKey, cmb);
  const reco::Candidate* Lep    = FullLepEvt->lepton(hypoKey, cmb);
  const reco::Candidate* NuBar  = FullLepEvt->neutrinoBar(hypoKey, cmb);

  AddFourMomenta addFourMomenta;

  // combined top pair object
  reco::CompositeCandidate* TtBar = new reco::CompositeCandidate();
  TtBar->addDaughter(*Top);
  TtBar->addDaughter(*TopBar);
  addFourMomenta.set(*TtBar);

  // combined lepton pair object
  reco::CompositeCandidate* LepPair = new reco::CompositeCandidate();
  LepPair->addDaughter(*Lep);
  LepPair->addDaughter(*LepBar);
  addFourMomenta.set(*LepPair);

  // dijet pair object
  reco::CompositeCandidate* JetPair = new reco::CompositeCandidate();
  JetPair->addDaughter(*B);
  JetPair->addDaughter(*BBar);
  addFourMomenta.set(*JetPair);

  fillKinHistos(TopKin_,    *Top,    weight);
  fillKinHistos(WplusKin_,  *Wplus,  weight);
  fillKinHistos(BKin_,      *B,      weight);
  fillKinHistos(LepBarKin_, *LepBar, weight);
  fillKinHistos(NuKin_,     *Nu,     weight);

  fillKinHistos(TopBarKin_, *TopBar, weight);
  fillKinHistos(WminusKin_, *Wminus, weight);
  fillKinHistos(BBarKin_,   *BBar,   weight);
  fillKinHistos(LepKin_,    *Lep,    weight);
  fillKinHistos(NuBarKin_,  *NuBar,  weight);

  fillKinHistos(TtBarKin_,   *TtBar,   weight);
  fillKinHistos(LepPairKin_, *LepPair, weight);
  fillKinHistos(JetPairKin_, *JetPair, weight);

  diLeptonMassVsLeptonsPt_    ->Fill(LepPair->mass(), Lep->pt(), weight);
  diLeptonMassVsLeptonsPt_    ->Fill(LepPair->mass(), LepBar->pt(), weight);
  diLeptonMassVsLeptonsEta_   ->Fill(LepPair->mass(), Lep->eta(), weight);
  diLeptonMassVsLeptonsEta_   ->Fill(LepPair->mass(), LepBar->eta(), weight);
  diLeptonMassVsDiLeptonPt_   ->Fill(LepPair->mass(), LepPair->pt(), weight);
  diLeptonMassVsJetsPt_       ->Fill(LepPair->mass(), B->pt(), weight);
  diLeptonMassVsJetsPt_       ->Fill(LepPair->mass(), BBar->pt(), weight);
  diLeptonMassVsJetsEta_      ->Fill(LepPair->mass(), B->eta(), weight);
  diLeptonMassVsJetsEta_      ->Fill(LepPair->mass(), BBar->eta(), weight);  
  diLeptonMassVsTopPt_	      ->Fill(LepPair->mass(), Top->pt(), weight);
  diLeptonMassVsTopPt_	      ->Fill(LepPair->mass(), TopBar->pt(), weight);
  diLeptonMassVsTopRapidity_  ->Fill(LepPair->mass(), Top->rapidity(), weight);
  diLeptonMassVsTopRapidity_  ->Fill(LepPair->mass(), TopBar->rapidity(), weight);
  diLeptonMassVsTtBarPt_      ->Fill(LepPair->mass(), TtBar->pt(), weight);
  diLeptonMassVsTtBarRapidity_->Fill(LepPair->mass(), TtBar->rapidity(), weight);
  diLeptonMassVsTtBarMass_    ->Fill(LepPair->mass(), TtBar->mass(), weight);

  // count again the number of b-jets this time only looking at the two jets used in the event reconstruction
  // note( cmb = 0 in this case -> choose hypotheses independent of b-tags)
  nTCHEL  = 0;
  nTCHEM  = 0;
  nSSVHEM = 0;

  std::vector<pat::Jet>::const_reference jet1 = jets->at(FullLepEvt->jetLeptonCombination(hypoKey)[0]);
  std::vector<pat::Jet>::const_reference jet2 = jets->at(FullLepEvt->jetLeptonCombination(hypoKey)[1]);	
  
  if(jet1.bDiscriminator("trackCountingHighEffBJetTags")>1.7)
    nTCHEL++;
  if(jet1.bDiscriminator("trackCountingHighEffBJetTags")>3.3)
    nTCHEM++;    
  if(jet1.bDiscriminator("simpleSecondaryVertexHighEffBJetTags")>1.74)
    nSSVHEM++;     
    
  if(jet2.bDiscriminator("trackCountingHighEffBJetTags")>1.7)
    nTCHEL++;
  if(jet2.bDiscriminator("trackCountingHighEffBJetTags")>3.3)
    nTCHEM++;    
  if(jet2.bDiscriminator("simpleSecondaryVertexHighEffBJetTags")>1.74)
    nSSVHEM++;

  nBtagsTCHELVsLeptonsPt_    ->Fill(nTCHEL, Lep->pt(), weight);
  nBtagsTCHELVsLeptonsPt_    ->Fill(nTCHEL, LepBar->pt(), weight);
  nBtagsTCHELVsLeptonsEta_   ->Fill(nTCHEL, Lep->eta(), weight);
  nBtagsTCHELVsLeptonsEta_   ->Fill(nTCHEL, LepBar->eta(), weight);
  nBtagsTCHELVsDiLeptonPt_   ->Fill(nTCHEL, LepPair->pt(), weight);
  nBtagsTCHELVsDiLeptonMass_ ->Fill(nTCHEL, LepPair->mass(), weight);
  nBtagsTCHELVsJetsPt_       ->Fill(nTCHEL, B->pt(), weight);
  nBtagsTCHELVsJetsPt_       ->Fill(nTCHEL, BBar->pt(), weight);
  nBtagsTCHELVsJetsEta_      ->Fill(nTCHEL, B->eta(), weight);
  nBtagsTCHELVsJetsEta_      ->Fill(nTCHEL, BBar->eta(), weight);  
  nBtagsTCHELVsTopPt_	     ->Fill(nTCHEL, Top->pt(), weight);
  nBtagsTCHELVsTopPt_	     ->Fill(nTCHEL, TopBar->pt(), weight);
  nBtagsTCHELVsTopRapidity_  ->Fill(nTCHEL, Top->rapidity(), weight);
  nBtagsTCHELVsTopRapidity_  ->Fill(nTCHEL, TopBar->rapidity(), weight);
  nBtagsTCHELVsTtBarPt_	     ->Fill(nTCHEL, TtBar->pt(), weight);
  nBtagsTCHELVsTtBarRapidity_->Fill(nTCHEL, TtBar->rapidity(), weight);
  nBtagsTCHELVsTtBarMass_    ->Fill(nTCHEL, TtBar->mass(), weight);

  nBtagsTCHEMVsLeptonsPt_    ->Fill(nTCHEM, Lep->pt(), weight);
  nBtagsTCHEMVsLeptonsPt_    ->Fill(nTCHEM, LepBar->pt(), weight);
  nBtagsTCHEMVsLeptonsEta_   ->Fill(nTCHEM, Lep->eta(), weight);
  nBtagsTCHEMVsLeptonsEta_   ->Fill(nTCHEM, LepBar->eta(), weight);
  nBtagsTCHEMVsDiLeptonPt_   ->Fill(nTCHEM, LepPair->pt(), weight);
  nBtagsTCHEMVsDiLeptonMass_ ->Fill(nTCHEM, LepPair->mass(), weight);
  nBtagsTCHEMVsJetsPt_       ->Fill(nTCHEM, B->pt(), weight);
  nBtagsTCHEMVsJetsPt_       ->Fill(nTCHEM, BBar->pt(), weight);
  nBtagsTCHEMVsJetsEta_      ->Fill(nTCHEM, B->eta(), weight);
  nBtagsTCHEMVsJetsEta_      ->Fill(nTCHEM, BBar->eta(), weight);    
  nBtagsTCHEMVsTopPt_	     ->Fill(nTCHEM, Top->pt(), weight);
  nBtagsTCHEMVsTopPt_	     ->Fill(nTCHEM, TopBar->pt(), weight);
  nBtagsTCHEMVsTopRapidity_  ->Fill(nTCHEM, Top->rapidity(), weight);
  nBtagsTCHEMVsTopRapidity_  ->Fill(nTCHEM, TopBar->rapidity(), weight);
  nBtagsTCHEMVsTtBarPt_	     ->Fill(nTCHEM, TtBar->pt(), weight);
  nBtagsTCHEMVsTtBarRapidity_->Fill(nTCHEM, TtBar->rapidity(), weight);
  nBtagsTCHEMVsTtBarMass_    ->Fill(nTCHEM, TtBar->mass(), weight);

  nBtagsSSVHEMVsLeptonsPt_    ->Fill(nSSVHEM, Lep->pt(), weight);
  nBtagsSSVHEMVsLeptonsPt_    ->Fill(nSSVHEM, LepBar->pt(), weight);
  nBtagsSSVHEMVsLeptonsEta_   ->Fill(nSSVHEM, Lep->eta(), weight);
  nBtagsSSVHEMVsLeptonsEta_   ->Fill(nSSVHEM, LepBar->eta(), weight);
  nBtagsSSVHEMVsDiLeptonPt_   ->Fill(nSSVHEM, LepPair->pt(), weight);
  nBtagsSSVHEMVsDiLeptonMass_ ->Fill(nSSVHEM, LepPair->mass(), weight);
  nBtagsSSVHEMVsJetsPt_       ->Fill(nSSVHEM, B->pt(), weight);
  nBtagsSSVHEMVsJetsPt_       ->Fill(nSSVHEM, BBar->pt(), weight);
  nBtagsSSVHEMVsJetsEta_      ->Fill(nSSVHEM, B->eta(), weight);
  nBtagsSSVHEMVsJetsEta_      ->Fill(nSSVHEM, BBar->eta(), weight);    
  nBtagsSSVHEMVsTopPt_	      ->Fill(nSSVHEM, Top->pt(), weight);
  nBtagsSSVHEMVsTopPt_	      ->Fill(nSSVHEM, TopBar->pt(), weight);
  nBtagsSSVHEMVsTopRapidity_  ->Fill(nSSVHEM, Top->rapidity(), weight);
  nBtagsSSVHEMVsTopRapidity_  ->Fill(nSSVHEM, TopBar->rapidity(), weight);
  nBtagsSSVHEMVsTtBarPt_      ->Fill(nSSVHEM, TtBar->pt(), weight);
  nBtagsSSVHEMVsTtBarRapidity_->Fill(nSSVHEM, TtBar->rapidity(), weight);
  nBtagsSSVHEMVsTtBarMass_    ->Fill(nSSVHEM, TtBar->mass(), weight);

  // combine masses of leptons and jets
  reco::CompositeCandidate* LepBarB = new reco::CompositeCandidate();
  LepBarB->addDaughter(*LepBar);
  LepBarB->addDaughter(*B);
  addFourMomenta.set(*LepBarB);
  LepBarBMass_->Fill(LepBarB->mass(), weight);

  reco::CompositeCandidate* LepB = new reco::CompositeCandidate();
  LepB->addDaughter(*Lep);
  LepB->addDaughter(*B);
  addFourMomenta.set(*LepB);
  LepBMass_->Fill(LepB->mass(), weight);

  reco::CompositeCandidate* LepBBar = new reco::CompositeCandidate();
  LepBBar->addDaughter(*Lep);
  LepBBar->addDaughter(*BBar);
  addFourMomenta.set(*LepBBar);
  LepBBarMass_->Fill(LepBBar->mass(), weight);

  reco::CompositeCandidate* LepBarBBar = new reco::CompositeCandidate();
  LepBarBBar->addDaughter(*LepBar);
  LepBarBBar->addDaughter(*BBar);
  addFourMomenta.set(*LepBarBBar);
  LepBarBBarMass_->Fill(LepBarBBar->mass(), weight);


  // -----------------------
  // fill generator and pull histos for kinematic variables
  // -----------------------
  if( !isSignalMC_ || !FullLepEvt->genEvent() || !FullLepEvt->genEvent()->isFullLeptonic()){
    delete TtBar;
    delete LepPair;
    delete JetPair;
    return;
  }

  const reco::Candidate* genTop    = FullLepEvt->genTop();
  const reco::Candidate* genWplus  = FullLepEvt->genWPlus();
  const reco::Candidate* genB      = FullLepEvt->genB();
  const reco::Candidate* genLepBar = FullLepEvt->genLeptonBar();
  const reco::Candidate* genNu     = FullLepEvt->genNeutrino();

  const reco::Candidate* genTopBar = FullLepEvt->genTopBar();
  const reco::Candidate* genWminus = FullLepEvt->genWMinus();
  const reco::Candidate* genBBar   = FullLepEvt->genBBar();
  const reco::Candidate* genLep    = FullLepEvt->genLepton();
  const reco::Candidate* genNuBar  = FullLepEvt->genNeutrinoBar();

  // combined top pair object
  reco::CompositeCandidate* genTtBar = new reco::CompositeCandidate();
  genTtBar->addDaughter(*genTop);
  genTtBar->addDaughter(*genTopBar);
  addFourMomenta.set(*genTtBar);

  // combined lepton pair object
  reco::CompositeCandidate* genLepPair = new reco::CompositeCandidate();
  genLepPair->addDaughter(*genLep);
  genLepPair->addDaughter(*genLepBar);
  addFourMomenta.set(*genLepPair);

  // combined lepton pair object
  reco::CompositeCandidate* genJetPair = new reco::CompositeCandidate();
  genJetPair->addDaughter(*genB);
  genJetPair->addDaughter(*genBBar);
  addFourMomenta.set(*genJetPair);

  // gen particle distributions
  fillKinHistos(TopGen_,    *genTop,    weight);
  fillKinHistos(WplusGen_,  *genWplus,  weight);
  fillKinHistos(BGen_,      *genB,      weight);
  fillKinHistos(LepBarGen_, *genLepBar, weight);
  fillKinHistos(NuGen_,     *genNu,     weight);

  fillKinHistos(TopBarGen_, *genTopBar, weight);
  fillKinHistos(WminusGen_, *genWminus, weight);
  fillKinHistos(BBarGen_,   *genBBar,   weight);
  fillKinHistos(LepGen_,    *genLep,    weight);
  fillKinHistos(NuBarGen_,  *genNuBar,  weight);

  fillKinHistos(TtBarGen_,   *genTtBar,   weight);
  fillKinHistos(LepPairGen_, *genLepPair, weight);
  fillKinHistos(JetPairGen_, *genJetPair, weight);

  // pull distributions
  fillPullHistos(TopPull_,    *TopPull2D_,    *Top,	*genTop,    weight);
  fillPullHistos(WplusPull_,  *WplusPull2D_,  *Wplus,	*genWplus,  weight);
  fillPullHistos(BPull_,      *BPull2D_,      *B,	*genB,      weight);
  fillPullHistos(LepBarPull_, *LepBarPull2D_, *LepBar,  *genLepBar, weight);
  fillPullHistos(NuPull_,     *NuPull2D_,     *Nu,	*genNu,     weight);

  fillPullHistos(TopBarPull_, *TopBarPull2D_, *TopBar,  *genTopBar, weight);
  fillPullHistos(WminusPull_, *WminusPull2D_, *Wminus,  *genWminus, weight);
  fillPullHistos(BBarPull_,   *BBarPull2D_,   *BBar,	*genBBar,   weight);
  fillPullHistos(LepPull_,    *LepPull2D_,    *Lep,	*genLep,    weight);
  fillPullHistos(NuBarPull_,  *NuBarPull2D_,  *NuBar,	*genNuBar,  weight);

  fillPullHistos(TtBarPull_,   *TtBarPull2D_,	*TtBar,   *genTtBar,   weight);
  fillPullHistos(LepPairPull_, *LepPairPull2D_, *LepPair, *genLepPair, weight);
  fillPullHistos(JetPairPull_, *JetPairPull2D_, *JetPair, *genJetPair, weight);

  // 2D distributions
  fill2DHistos(Top2D_,     *Top,     *genTop,     weight);
  fill2DHistos(Wplus2D_,   *Wplus,   *genWplus,   weight);
  fill2DHistos(B2D_,       *B,       *genB,       weight);
  fill2DHistos(LepBar2D_,  *LepBar,  *genLepBar,  weight);
  fill2DHistos(Nu2D_,      *Nu,      *genNu,      weight);

  fill2DHistos(TopBar2D_,  *TopBar,  *genTopBar,  weight);
  fill2DHistos(Wminus2D_,  *Wminus,  *genWminus,  weight);
  fill2DHistos(BBar2D_,    *BBar,    *genBBar,    weight);
  fill2DHistos(Lep2D_,     *Lep,     *genLep,     weight);
  fill2DHistos(NuBar2D_,   *NuBar,   *genNuBar,   weight);

  fill2DHistos(TtBar2D_,   *TtBar,   *genTtBar,   weight);
  fill2DHistos(LepPair2D_, *LepPair, *genLepPair, weight);
  fill2DHistos(JetPair2D_, *JetPair, *genJetPair, weight);

  // count jets matched to b flavour jets
  int nMatchedBjets = 0;
  if(isRealBJet(B,genB,genBBar)) 
    nMatchedBjets++;    
  if(isRealBJet(BBar,genB,genBBar)) 
    nMatchedBjets++;

  nMatchesVsLeptonsPt_    ->Fill(nMatchedBjets, Lep->pt());
  nMatchesVsLeptonsPt_    ->Fill(nMatchedBjets, LepBar->pt());
  nMatchesVsLeptonsEta_   ->Fill(nMatchedBjets, Lep->eta());
  nMatchesVsLeptonsEta_   ->Fill(nMatchedBjets, LepBar->eta());
  nMatchesVsDiLeptonPt_   ->Fill(nMatchedBjets, LepPair->pt());
  nMatchesVsDiLeptonMass_ ->Fill(nMatchedBjets, LepPair->mass());
  nMatchesVsJetsPt_       ->Fill(nMatchedBjets, B->pt());
  nMatchesVsJetsPt_       ->Fill(nMatchedBjets, BBar->pt());
  nMatchesVsJetsEta_      ->Fill(nMatchedBjets, B->eta());
  nMatchesVsJetsEta_      ->Fill(nMatchedBjets, BBar->eta());    
  nMatchesVsTopPt_	  ->Fill(nMatchedBjets, Top->pt());
  nMatchesVsTopPt_	  ->Fill(nMatchedBjets, TopBar->pt());
  nMatchesVsTopRapidity_  ->Fill(nMatchedBjets, Top->rapidity());
  nMatchesVsTopRapidity_  ->Fill(nMatchedBjets, TopBar->rapidity());
  nMatchesVsTtBarPt_	  ->Fill(nMatchedBjets, TtBar->pt());
  nMatchesVsTtBarRapidity_->Fill(nMatchedBjets, TtBar->rapidity());
  nMatchesVsTtBarMass_    ->Fill(nMatchedBjets, TtBar->mass());

  nMatchVsNtagsTCHEL_ ->Fill(nTCHEL, nMatchedBjets, weight); 
  nMatchVsNtagsTCHEM_ ->Fill(nTCHEM, nMatchedBjets, weight); 
  nMatchVsNtagsSSVHEM_->Fill(nSSVHEM,nMatchedBjets, weight);

  // reset counter to count now the number of jets matched to the correct jet
  nMatchedBjets = 0;
  if(isRealBJet(B,genB)) 
    nMatchedBjets++;    
  if(isRealBJet(BBar,genBBar)) 
    nMatchedBjets++;  
  
  nCorrectAssignments_->Fill(nMatchedBjets);

  delete TtBar;
  delete LepPair;
  delete JetPair;
  delete genTtBar;
  delete genLepPair;
  delete genJetPair;

  // indices Hist
  const TtEvent::HypoClassKey genKey = TtEvent::HypoClassKey(3);
  const TtEvent::HypoClassKey kinKey = TtEvent::HypoClassKey(6);

  if(FullLepEvt->isHypoValid(genKey) && FullLepEvt->isHypoValid(kinKey)){
    compare_->Fill(FullLepEvt->jetLeptonCombination(genKey)[0]   , FullLepEvt->jetLeptonCombination(kinKey, cmb)[0]   , weight);
    compare_->Fill(FullLepEvt->jetLeptonCombination(genKey)[1]+ 4, FullLepEvt->jetLeptonCombination(kinKey, cmb)[1]+ 4, weight);
    compare_->Fill(FullLepEvt->jetLeptonCombination(genKey)[2]+ 8, FullLepEvt->jetLeptonCombination(kinKey, cmb)[2]+ 8, weight);
    compare_->Fill(FullLepEvt->jetLeptonCombination(genKey)[3]+12, FullLepEvt->jetLeptonCombination(kinKey, cmb)[3]+12, weight);
    compare_->Fill(FullLepEvt->jetLeptonCombination(genKey)[4]+16, FullLepEvt->jetLeptonCombination(kinKey, cmb)[4]+16, weight);
    compare_->Fill(FullLepEvt->jetLeptonCombination(genKey)[5]+20, FullLepEvt->jetLeptonCombination(kinKey, cmb)[5]+20, weight);
  }
}

/// everything that has to be done after the event loop
void
FullLepKinAnalyzer::endJob()
{
}

/// book histograms for reconstructed particles properties: Pt, E, Eta, rapidity, Phi, m
void
FullLepKinAnalyzer::bookKinHistos(edm::Service<TFileService>& fs)
{

  NameScheme ns("kin");

  TopKin_.push_back( fs->make<TH1D>(ns.name("TopPt"      ), "p_{t} (t) [GeV]",100,  0. , 500. ) );
  TopKin_.push_back( fs->make<TH1D>(ns.name("TopEta"     ), "#eta (t)"       ,100, -5.0,   5.0) );
  TopKin_.push_back( fs->make<TH1D>(ns.name("TopRapidity"), "rapidity (t)"   ,100, -5.0,   5.0) );
  TopKin_.push_back( fs->make<TH1D>(ns.name("TopPhi"     ), "#phi (t)"       , 62, -3.1,   3.1) );
  TopKin_.push_back( fs->make<TH1D>(ns.name("TopMass"    ), "M (top) [GeV]"  ,100,-50. , 450. ) );

  WplusKin_.push_back( fs->make<TH1D>(ns.name("WplusPt"      ), "p_{t} (W^{+}) [GeV]",100,  0. , 500. ) );
  WplusKin_.push_back( fs->make<TH1D>(ns.name("WplusEta"     ), "#eta (W^{+})"       ,100, -5.0,   5.0) );
  WplusKin_.push_back( fs->make<TH1D>(ns.name("WplusRapidity"), "rapidity (W^{+})"   ,100, -5.0,   5.0) );
  WplusKin_.push_back( fs->make<TH1D>(ns.name("WplusPhi"     ), "#phi (W^{+})"       , 62, -3.1,   3.1) );
  WplusKin_.push_back( fs->make<TH1D>(ns.name("WplusMass"    ), "M (W^{+}) [GeV]"    ,100,  0. , 500. ) );

  BKin_.push_back( fs->make<TH1D>(ns.name("BPt"      ), "p_{t} (b) [GeV]",100,  0. , 500. ) );
  BKin_.push_back( fs->make<TH1D>(ns.name("BEta"     ), "#eta (b)"       ,100, -5.0,   5.0) );
  BKin_.push_back( fs->make<TH1D>(ns.name("BRapidity"), "rapidity (b)"   ,100, -5.0,   5.0) );
  BKin_.push_back( fs->make<TH1D>(ns.name("BPhi"     ), "#phi (b)"       , 62, -3.1,   3.1) );
  BKin_.push_back( fs->make<TH1D>(ns.name("BMass"    ), "M (b) [GeV]"    ,100,  0. , 500. ) );

  LepBarKin_.push_back( fs->make<TH1D>(ns.name("LepBarPt"      ), "p_{t} (l^{+}) [GeV]",100,  0. , 500. ) );
  LepBarKin_.push_back( fs->make<TH1D>(ns.name("LepBarEta"     ), "#eta (l^{+})"       ,100, -5.0,   5.0) );
  LepBarKin_.push_back( fs->make<TH1D>(ns.name("LepBarRapidity"), "rapidity (l^{+})"   ,100, -5.0,   5.0) );
  LepBarKin_.push_back( fs->make<TH1D>(ns.name("LepBarPhi"     ), "#phi (l^{+})"       , 62, -3.1,   3.1) );
  LepBarKin_.push_back( fs->make<TH1D>(ns.name("LepBarMass"    ), "M (l^{+}) [GeV]"    ,100,  0. , 500. ) );

  NuKin_.push_back( fs->make<TH1D>(ns.name("NuPt"      ), "p_{t} (#nu^{-}) [GeV]",100,  0. , 500. ) );
  NuKin_.push_back( fs->make<TH1D>(ns.name("NuEta"     ), "#eta (#nu^{-})"       ,100, -5.0,   5.0) );
  NuKin_.push_back( fs->make<TH1D>(ns.name("NuRapidity"), "rapidity (#nu^{-})"   ,100, -5.0,   5.0) );
  NuKin_.push_back( fs->make<TH1D>(ns.name("NuPhi"     ), "#phi (#nu^{-})"       , 62, -3.1,   3.1) );
  NuKin_.push_back( fs->make<TH1D>(ns.name("NuMass"    ), "M (#nu^{-}) [GeV]"    ,100,  0. , 500. ) );

  TopBarKin_.push_back( fs->make<TH1D>(ns.name("TopBarPt"      ), "p_{t} (#bar{t}) [GeV]",100,  0. , 500. ) );
  TopBarKin_.push_back( fs->make<TH1D>(ns.name("TopBarEta"     ), "#eta (#bar{t})"       ,100, -5.0,   5.0) );
  TopBarKin_.push_back( fs->make<TH1D>(ns.name("TopBarRapidity"), "rapidity (#bar{t})"   ,100, -5.0,   5.0) );
  TopBarKin_.push_back( fs->make<TH1D>(ns.name("TopBarPhi"     ), "#phi (#bar{t})"       , 62, -3.1,   3.1) );
  TopBarKin_.push_back( fs->make<TH1D>(ns.name("TopBarMass"    ), "M (#bar{t}) [GeV]"    ,100,  0. , 500. ) );

  WminusKin_.push_back( fs->make<TH1D>(ns.name("WminusPt"      ), "p_{t} (W^{-}) [GeV]",100,  0. , 500. ) );
  WminusKin_.push_back( fs->make<TH1D>(ns.name("WminusEta"     ), "#eta (W^{-})"       ,100, -5.0,   5.0) );
  WminusKin_.push_back( fs->make<TH1D>(ns.name("WminusRapidity"), "rapidity (W^{-})"   ,100, -5.0,   5.0) );
  WminusKin_.push_back( fs->make<TH1D>(ns.name("WminusPhi"     ), "#phi (W^{-})"       , 62, -3.1,   3.1) );
  WminusKin_.push_back( fs->make<TH1D>(ns.name("WminusMass"    ), "M (W^{-}) [GeV]"    ,100,  0. , 500. ) );

  BBarKin_.push_back( fs->make<TH1D>(ns.name("BBarPt"      ), "p_{t} (#bar{b}) [GeV]",100,  0. , 500. ) );
  BBarKin_.push_back( fs->make<TH1D>(ns.name("BBarEta"     ), "#eta (#bar{b})"       ,100, -5.0,   5.0) );
  BBarKin_.push_back( fs->make<TH1D>(ns.name("BBarRapidity"), "rapidity (#bar{b})"   ,100, -5.0,   5.0) );
  BBarKin_.push_back( fs->make<TH1D>(ns.name("BBarPhi"     ), "#phi (#bar{b})"       , 62, -3.1,   3.1) );
  BBarKin_.push_back( fs->make<TH1D>(ns.name("BBarMass"    ), "M (#bar{b}) [GeV]"    ,100,  0. , 500. ) );

  LepKin_.push_back( fs->make<TH1D>(ns.name("LepPt"      ), "p_{t} (l^{-}) [GeV]",100,  0. , 500. ) );
  LepKin_.push_back( fs->make<TH1D>(ns.name("LepEta"     ), "#eta (l^{-})"       ,100, -5.0,   5.0) );
  LepKin_.push_back( fs->make<TH1D>(ns.name("LepRapidity" ), "rapidity (l^{-})"   ,100, -5.0,   5.0) );
  LepKin_.push_back( fs->make<TH1D>(ns.name("LepPhi"     ), "#phi (l^{-})"       , 62, -3.1,   3.1) );
  LepKin_.push_back( fs->make<TH1D>(ns.name("LepMass"    ), "M (l^{-}) [GeV]"    ,100,  0. , 500. ) );

  NuBarKin_.push_back( fs->make<TH1D>(ns.name("NuBarPt"      ), "p_{t} (#bar{#nu}^{-}) [GeV]",100,  0. , 500. ) );
  NuBarKin_.push_back( fs->make<TH1D>(ns.name("NuBarEta"     ), "#eta (#bar{#nu}^{-})"       ,100, -5.0,   5.0) );
  NuBarKin_.push_back( fs->make<TH1D>(ns.name("NuBarRapidity"), "rapidity (#bar{#nu}^{-})"   ,100, -5.0,   5.0) );
  NuBarKin_.push_back( fs->make<TH1D>(ns.name("NuBarPhi"     ), "#phi (#bar{#nu}^{-})"       , 62, -3.1,   3.1) );
  NuBarKin_.push_back( fs->make<TH1D>(ns.name("NuBarMass"    ), "M (#bar{#nu}^{-}) [GeV]"    ,100,  0. , 500. ) );

  TtBarKin_.push_back( fs->make<TH1D>(ns.name("TtBarPt"      ), "p_{t} (t#bar{t}) [GeV]",100,  0. , 500. ) );
  TtBarKin_.push_back( fs->make<TH1D>(ns.name("TtBarEta"     ), "#eta (t#bar{t})"       ,100, -5.0,   5.0) );
  TtBarKin_.push_back( fs->make<TH1D>(ns.name("TtBarRapidity"), "rapidity (t#bar{t})"   ,100, -5.0,   5.0) );
  TtBarKin_.push_back( fs->make<TH1D>(ns.name("TtBarPhi"     ), "#phi (t#bar{t})"       , 62, -3.1,   3.1) );
  TtBarKin_.push_back( fs->make<TH1D>(ns.name("TtBarMass"    ), "M (t#bar{t}) [GeV]"    ,400,  0. ,2000. ) );

  LepPairKin_.push_back( fs->make<TH1D>(ns.name("LepPairPt"      ), "p_{t} (l^{+}l^{-}) [GeV]",100,  0. , 500. ) );
  LepPairKin_.push_back( fs->make<TH1D>(ns.name("LepPairEta"     ), "#eta (l^{+}l^{-})"       ,100, -5.0,   5.0) );
  LepPairKin_.push_back( fs->make<TH1D>(ns.name("LepPairRapidity"), "rapidity (l^{+}l^{-})"   ,100, -5.0,   5.0) );
  LepPairKin_.push_back( fs->make<TH1D>(ns.name("LepPairPhi"     ), "#phi (l^{+}l^{-})"       , 62, -3.1,   3.1) );
  LepPairKin_.push_back( fs->make<TH1D>(ns.name("LepPairMass"    ), "M (l^{+}l^{-}) [GeV]"    ,200,  0. , 400. ) );

  JetPairKin_.push_back( fs->make<TH1D>(ns.name("JetPairPt"      ), "p_{t} (jj) [GeV]",100,  0. , 500. ) );
  JetPairKin_.push_back( fs->make<TH1D>(ns.name("JetPairEta"     ), "#eta (jj)"       ,100, -5.0,   5.0) );
  JetPairKin_.push_back( fs->make<TH1D>(ns.name("JetPairRapidity"), "rapidity (jj)"   ,100, -5.0,   5.0) );
  JetPairKin_.push_back( fs->make<TH1D>(ns.name("JetPairPhi"     ), "#phi (jj)"       , 62, -3.1,   3.1) );
  JetPairKin_.push_back( fs->make<TH1D>(ns.name("JetPairMass"    ), "M (jj) [GeV]"    ,100,  0. , 500. ) );
}



/// book histograms for generated particles properties: Pt, E, Eta, Phi, m
void
FullLepKinAnalyzer::bookGenHistos(edm::Service<TFileService>& fs)
{

  NameScheme ns("gen");

  TopGen_.push_back( fs->make<TH1D>(ns.name("TopPt"      ), "p_{t} (t) [GeV]",100,  0. , 500. ) );
  TopGen_.push_back( fs->make<TH1D>(ns.name("TopEta"     ), "#eta (t)"       ,100, -5.0,   5.0) );
  TopGen_.push_back( fs->make<TH1D>(ns.name("TopRapidity"), "rapidity (t)"   ,100, -5.0,   5.0) );
  TopGen_.push_back( fs->make<TH1D>(ns.name("TopPhi"     ), "#phi (t)"       , 62, -3.1,   3.1) );
  TopGen_.push_back( fs->make<TH1D>(ns.name("TopMass"    ), "M (top) [GeV]"  ,100,  0. , 500. ) );

  WplusGen_.push_back( fs->make<TH1D>(ns.name("WplusPt"      ), "p_{t} (W^{+}) [GeV]",100,  0. , 500. ) );
  WplusGen_.push_back( fs->make<TH1D>(ns.name("WplusEta"     ), "#eta (W^{+})"       ,100, -5.0,   5.0) );
  WplusGen_.push_back( fs->make<TH1D>(ns.name("WplusRapidity"), "rapidity (W^{+})"   ,100, -5.0,   5.0) );
  WplusGen_.push_back( fs->make<TH1D>(ns.name("WplusPhi"     ), "#phi (W^{+})"       , 62, -3.1,   3.1) );
  WplusGen_.push_back( fs->make<TH1D>(ns.name("WplusMass"    ), "M (W^{+}) [GeV]"    ,100,  0. , 500. ) );

  BGen_.push_back( fs->make<TH1D>(ns.name("BPt"      ), "p_{t} (b) [GeV]",100,  0. , 500. ) );
  BGen_.push_back( fs->make<TH1D>(ns.name("BEta"     ), "#eta (b)"       ,100, -5.0,   5.0) );
  BGen_.push_back( fs->make<TH1D>(ns.name("BRapidity"), "rapidity (b)"   ,100, -5.0,   5.0) );
  BGen_.push_back( fs->make<TH1D>(ns.name("BPhi"     ), "#phi (b)"       , 62, -3.1,   3.1) );
  BGen_.push_back( fs->make<TH1D>(ns.name("BMass"    ), "M (b) [GeV]"    ,100,  0. , 500. ) );

  LepBarGen_.push_back( fs->make<TH1D>(ns.name("LepBarPt"      ), "p_{t} (l^{+}) [GeV]",100,  0. , 500. ) );
  LepBarGen_.push_back( fs->make<TH1D>(ns.name("LepBarEta"     ), "#eta (l^{+})"       ,100, -5.0,   5.0) );
  LepBarGen_.push_back( fs->make<TH1D>(ns.name("LepBarRapidity"), "rapidity (l^{+})"   ,100, -5.0,   5.0) );
  LepBarGen_.push_back( fs->make<TH1D>(ns.name("LepBarPhi"     ), "#phi (l^{+})"       , 62, -3.1,   3.1) );
  LepBarGen_.push_back( fs->make<TH1D>(ns.name("LepBarMass"    ), "M (l^{+}) [GeV]"    ,100,  0. , 500. ) );

  NuGen_.push_back( fs->make<TH1D>(ns.name("NuPt"      ), "p_{t} (#nu^{-}) [GeV]",100,  0. , 500. ) );
  NuGen_.push_back( fs->make<TH1D>(ns.name("NuEta"     ), "#eta (#nu^{-})"	 ,100, -5.0,   5.0) );
  NuGen_.push_back( fs->make<TH1D>(ns.name("NuRapidity"), "rapidity (#nu^{-})"   ,100, -5.0,   5.0) );
  NuGen_.push_back( fs->make<TH1D>(ns.name("NuPhi"     ), "#phi (#nu^{-})"	 , 62, -3.1,   3.1) );
  NuGen_.push_back( fs->make<TH1D>(ns.name("NuMass"    ), "M (#nu^{-}) [GeV]"	 ,100,  0. , 500. ) );

  TopBarGen_.push_back( fs->make<TH1D>(ns.name("TopBarPt"      ), "p_{t} (#bar{t}) [GeV]",100,  0. , 500. ) );
  TopBarGen_.push_back( fs->make<TH1D>(ns.name("TopBarEta"     ), "#eta (#bar{t})"	 ,100, -5.0,   5.0) );
  TopBarGen_.push_back( fs->make<TH1D>(ns.name("TopBarRapidity"), "rapidity (#bar{t})"   ,100, -5.0,   5.0) );
  TopBarGen_.push_back( fs->make<TH1D>(ns.name("TopBarPhi"     ), "#phi (#bar{t})"	 , 62, -3.1,   3.1) );
  TopBarGen_.push_back( fs->make<TH1D>(ns.name("TopBarMass"    ), "M (#bar{t}) [GeV]"	 ,100,  0. , 500. ) );

  WminusGen_.push_back( fs->make<TH1D>(ns.name("WminusPt"      ), "p_{t} (W^{-}) [GeV]",100,  0. , 500. ) );
  WminusGen_.push_back( fs->make<TH1D>(ns.name("WminusEta"     ), "#eta (W^{-})"       ,100, -5.0,   5.0) );
  WminusGen_.push_back( fs->make<TH1D>(ns.name("WminusRapidity"), "rapidity (W^{-})"   ,100, -5.0,   5.0) );
  WminusGen_.push_back( fs->make<TH1D>(ns.name("WminusPhi"     ), "#phi (W^{-})"       , 62, -3.1,   3.1) );
  WminusGen_.push_back( fs->make<TH1D>(ns.name("WminusMass"    ), "M (W^{-}) [GeV]"    ,100,  0. , 500. ) );

  BBarGen_.push_back( fs->make<TH1D>(ns.name("BBarPt"	   ), "p_{t} (#bar{b}) [GeV]",100,  0. , 500. ) );
  BBarGen_.push_back( fs->make<TH1D>(ns.name("BBarEta"     ), "#eta (#bar{b})"       ,100, -5.0,   5.0) );
  BBarGen_.push_back( fs->make<TH1D>(ns.name("BBarRapidity"), "rapidity (#bar{b})"   ,100, -5.0,   5.0) );
  BBarGen_.push_back( fs->make<TH1D>(ns.name("BBarPhi"     ), "#phi (#bar{b})"       , 62, -3.1,   3.1) );
  BBarGen_.push_back( fs->make<TH1D>(ns.name("BBarMass"    ), "M (#bar{b}) [GeV]"    ,100,  0. , 500. ) );

  LepGen_.push_back( fs->make<TH1D>(ns.name("LepPt"	 ), "p_{t} (l^{-}) [GeV]",100,  0. , 500. ) );
  LepGen_.push_back( fs->make<TH1D>(ns.name("LepEta"	 ), "#eta (l^{-})"	 ,100, -5.0,   5.0) );
  LepGen_.push_back( fs->make<TH1D>(ns.name("LepRapidity"), "rapidity (l^{-})"   ,100, -5.0,   5.0) );
  LepGen_.push_back( fs->make<TH1D>(ns.name("LepPhi"	 ), "#phi (l^{-})"	 , 62, -3.1,   3.1) );
  LepGen_.push_back( fs->make<TH1D>(ns.name("LepMass"	 ), "M (l^{-}) [GeV]"	 ,100,  0. , 500. ) );

  NuBarGen_.push_back( fs->make<TH1D>(ns.name("NuBarPt"      ), "p_{t} (#bar{#nu}^{-}) [GeV]",100,  0. , 500. ) );
  NuBarGen_.push_back( fs->make<TH1D>(ns.name("NuBarEta"     ), "#eta (#bar{#nu}^{-})"       ,100, -5.0,   5.0) );
  NuBarGen_.push_back( fs->make<TH1D>(ns.name("NuBarRapidity"), "rapidity (#bar{#nu}^{-})"   ,100, -5.0,   5.0) );
  NuBarGen_.push_back( fs->make<TH1D>(ns.name("NuBarPhi"     ), "#phi (#bar{#nu}^{-})"       , 62, -3.1,   3.1) );
  NuBarGen_.push_back( fs->make<TH1D>(ns.name("NuBarMass"    ), "M (#bar{#nu}^{-}) [GeV]"    ,100,  0. , 500. ) );

  TtBarGen_.push_back( fs->make<TH1D>(ns.name("TtBarPt"      ), "p_{t} (t#bar{t}) [GeV]",100,  0. , 500. ) );
  TtBarGen_.push_back( fs->make<TH1D>(ns.name("TtBarEta"     ), "#eta (t#bar{t})"	,100, -5.0,   5.0) );
  TtBarGen_.push_back( fs->make<TH1D>(ns.name("TtBarRapidity"), "rapidity (t#bar{t})"	,100, -5.0,   5.0) );
  TtBarGen_.push_back( fs->make<TH1D>(ns.name("TtBarPhi"     ), "#phi (t#bar{t})"	, 62, -3.1,   3.1) );
  TtBarGen_.push_back( fs->make<TH1D>(ns.name("TtBarMass"    ), "M (t#bar{t}) [GeV]"	,400,  0. ,2000. ) );

  LepPairGen_.push_back( fs->make<TH1D>(ns.name("LepPairPt"	 ), "p_{t} (l^{+}l^{-}) [GeV]",100,  0. , 500. ) );
  LepPairGen_.push_back( fs->make<TH1D>(ns.name("LepPairEta"	 ), "#eta (l^{+}l^{-})"       ,100, -5.0,   5.0) );
  LepPairGen_.push_back( fs->make<TH1D>(ns.name("LepPairRapidity"), "rapidity (l^{+}l^{-})"   ,100, -5.0,   5.0) );
  LepPairGen_.push_back( fs->make<TH1D>(ns.name("LepPairPhi"	 ), "#phi (l^{+}l^{-})"       , 62, -3.1,   3.1) );
  LepPairGen_.push_back( fs->make<TH1D>(ns.name("LepPairMass"	 ), "M (l^{+}l^{-}) [GeV]"    ,200,  0. , 400. ) );

  JetPairGen_.push_back( fs->make<TH1D>(ns.name("JetPairPt"	 ), "p_{t} (jj) [GeV]",100,  0. , 500. ) );
  JetPairGen_.push_back( fs->make<TH1D>(ns.name("JetPairEta"	 ), "#eta (jj)"       ,100, -5.0,   5.0) );
  JetPairGen_.push_back( fs->make<TH1D>(ns.name("JetPairRapidity"), "rapidity (jj)"   ,100, -5.0,   5.0) );
  JetPairGen_.push_back( fs->make<TH1D>(ns.name("JetPairPhi"	 ), "#phi (jj)"       , 62, -3.1,   3.1) );
  JetPairGen_.push_back( fs->make<TH1D>(ns.name("JetPairMass"	 ), "M (jj) [GeV]"    ,100,  0. , 500. ) );
}



/// book histograms for pull distributions of particle properties: Pt, E, Eta, Phi, m
void
FullLepKinAnalyzer::bookPullHistos(edm::Service<TFileService>& fs)
{

  NameScheme ns("pull");

  TopPull_.push_back( fs->make<TH1D>(ns.name("TopPt"	  ), "p_{t} (t) [GeV]", 100,  -1 , 1. ) );
  TopPull_.push_back( fs->make<TH1D>(ns.name("TopEta"	  ), "#eta (t)"       , 100,  -1 , 1. ) );
  TopPull_.push_back( fs->make<TH1D>(ns.name("TopRapidity"), "rapidity (t)"  ,  100,  -1 , 1. ) );
  TopPull_.push_back( fs->make<TH1D>(ns.name("TopPhi"	  ), "#phi (t)"       , 100,  -1 , 1. ) );
  TopPull_.push_back( fs->make<TH1D>(ns.name("TopMass"    ), "M (top) [GeV]"  , 100,  -1 , 1. ) );

  WplusPull_.push_back( fs->make<TH1D>(ns.name("WplusPt"      ), "p_{t} (W^{+}) [GeV]", 100,  -1 , 1. ) );
  WplusPull_.push_back( fs->make<TH1D>(ns.name("WplusEta"     ), "#eta (W^{+})"       , 100,  -1 , 1. ) );
  WplusPull_.push_back( fs->make<TH1D>(ns.name("WplusRapidity"), "rapidity (W^{+})"   , 100,  -1 , 1. ) );
  WplusPull_.push_back( fs->make<TH1D>(ns.name("WplusPhi"     ), "#phi (W^{+})"       , 100,  -1 , 1. ) );
  WplusPull_.push_back( fs->make<TH1D>(ns.name("WplusMass"    ), "M (W^{+}) [GeV]"    , 100,  -1 , 1. ) );

  BPull_.push_back( fs->make<TH1D>(ns.name("BPt"      ), "p_{t} (b) [GeV]", 100,  -1 , 1. ) );
  BPull_.push_back( fs->make<TH1D>(ns.name("BEta"     ), "#eta (b)"	  , 100,  -1 , 1. ) );
  BPull_.push_back( fs->make<TH1D>(ns.name("BRapidity"), "rapidity (b)"   , 100,  -1 , 1. ) );
  BPull_.push_back( fs->make<TH1D>(ns.name("BPhi"     ), "#phi (b)"	  , 100,  -1 , 1. ) );
  BPull_.push_back( fs->make<TH1D>(ns.name("BMass"    ), "M (b) [GeV]"    , 100,  -1 , 1. ) );

  LepBarPull_.push_back( fs->make<TH1D>(ns.name("LepBarPt"	), "p_{t} (l^{+}) [GeV]", 100,  -1 , 1. ) );
  LepBarPull_.push_back( fs->make<TH1D>(ns.name("LepBarEta"	), "#eta (l^{+})"	, 100,  -1 , 1. ) );
  LepBarPull_.push_back( fs->make<TH1D>(ns.name("LepBarRapidity"), "rapidity (l^{+})"	, 100,  -1 , 1. ) );
  LepBarPull_.push_back( fs->make<TH1D>(ns.name("LepBarPhi"	), "#phi (l^{+})"	, 100,  -1 , 1. ) );
  LepBarPull_.push_back( fs->make<TH1D>(ns.name("LepBarMass"	), "M (l^{+}) [GeV]"	, 100,  -1 , 1. ) );

  NuPull_.push_back( fs->make<TH1D>(ns.name("NuPt"	), "p_{t} (#nu^{-}) [GeV]", 100,  -1 , 1. ) );
  NuPull_.push_back( fs->make<TH1D>(ns.name("NuEta"	), "#eta (#nu^{-})"	  , 100,  -1 , 1. ) );
  NuPull_.push_back( fs->make<TH1D>(ns.name("NuRapidity"), "#eta (#nu^{-})"	  , 100,  -1 , 1. ) );
  NuPull_.push_back( fs->make<TH1D>(ns.name("NuPhi"	), "#phi (#nu^{-})"	  , 100,  -1 , 1. ) );
  NuPull_.push_back( fs->make<TH1D>(ns.name("NuMass"	), "M (#nu^{-}) [GeV]"    , 100,  -1 , 1. ) );

  TopBarPull_.push_back( fs->make<TH1D>(ns.name("TopBarPt"	), "p_{t} (#bar{t}) [GeV]", 100,  -1 , 1. ) );
  TopBarPull_.push_back( fs->make<TH1D>(ns.name("TopBarEta"	), "#eta (#bar{t})"	  , 100,  -1 , 1. ) );
  TopBarPull_.push_back( fs->make<TH1D>(ns.name("TopBarRapidity"), "rapidity (#bar{t})"   , 100,  -1 , 1. ) );
  TopBarPull_.push_back( fs->make<TH1D>(ns.name("TopBarPhi"	), "#phi (#bar{t})"	  , 100,  -1 , 1. ) );
  TopBarPull_.push_back( fs->make<TH1D>(ns.name("TopBarMass"	), "M (#bar{t}) [GeV]"    , 100,  -1 , 1. ) );

  WminusPull_.push_back( fs->make<TH1D>(ns.name("WminusPt"	), "p_{t} (W^{-}) [GeV]", 100,  -1 , 1. ) );
  WminusPull_.push_back( fs->make<TH1D>(ns.name("WminusEta"	), "#eta (W^{-})"	, 100,  -1 , 1. ) );
  WminusPull_.push_back( fs->make<TH1D>(ns.name("WminusRapidity"), "rapidity (W^{-})"	, 100,  -1 , 1. ) );
  WminusPull_.push_back( fs->make<TH1D>(ns.name("WminusPhi"	), "#phi (W^{-})"	, 100,  -1 , 1. ) );
  WminusPull_.push_back( fs->make<TH1D>(ns.name("WminusMass"	), "M (W^{-}) [GeV]"	, 100,  -1 , 1. ) );

  BBarPull_.push_back( fs->make<TH1D>(ns.name("BBarPt"      ), "p_{t} (#bar{b}) [GeV]", 100,  -1 , 1. ) );
  BBarPull_.push_back( fs->make<TH1D>(ns.name("BBarEta"     ), "#eta (#bar{b})"       , 100,  -1 , 1. ) );
  BBarPull_.push_back( fs->make<TH1D>(ns.name("BBarRapidity"), "rapidity (#bar{b})"   , 100,  -1 , 1. ) );
  BBarPull_.push_back( fs->make<TH1D>(ns.name("BBarPhi"     ), "#phi (#bar{b})"       , 100,  -1 , 1. ) );
  BBarPull_.push_back( fs->make<TH1D>(ns.name("BBarMass"    ), "M (#bar{b}) [GeV]"    , 100,  -1 , 1. ) );

  LepPull_.push_back( fs->make<TH1D>(ns.name("LepPt"	  ), "p_{t} (l^{-}) [GeV]", 100,  -1 , 1. ) );
  LepPull_.push_back( fs->make<TH1D>(ns.name("LepEta"	  ), "#eta (l^{-})"	  , 100,  -1 , 1. ) );
  LepPull_.push_back( fs->make<TH1D>(ns.name("LepRapidity"), "rapidity (l^{-})"   , 100,  -1 , 1. ) );
  LepPull_.push_back( fs->make<TH1D>(ns.name("LepPhi"	  ), "#phi (l^{-})"	  , 100,  -1 , 1. ) );
  LepPull_.push_back( fs->make<TH1D>(ns.name("LepMass"    ), "M (l^{-}) [GeV]"    , 100,  -1 , 1. ) );

  NuBarPull_.push_back( fs->make<TH1D>(ns.name("NuBarPt"      ), "p_{t} (#bar{#nu}^{-}) [GeV]", 100,  -1 , 1. ) );
  NuBarPull_.push_back( fs->make<TH1D>(ns.name("NuBarEta"     ), "#eta (#bar{#nu}^{-})"       , 100,  -1 , 1. ) );
  NuBarPull_.push_back( fs->make<TH1D>(ns.name("NuBarRapidity"), "rapidity (#bar{#nu}^{-})"   , 100,  -1 , 1. ) );
  NuBarPull_.push_back( fs->make<TH1D>(ns.name("NuBarPhi"     ), "#phi (#bar{#nu}^{-})"       , 100,  -1 , 1. ) );
  NuBarPull_.push_back( fs->make<TH1D>(ns.name("NuBarMass"    ), "M (#bar{#nu}^{-}) [GeV]"    , 100,  -1 , 1. ) );

  TtBarPull_.push_back( fs->make<TH1D>(ns.name("TtBarPt"      ), "p_{t} (t#bar{t}) [GeV]", 100,  -1 , 1. ) );
  TtBarPull_.push_back( fs->make<TH1D>(ns.name("TtBarEta"     ), "#eta (t#bar{t})"	 , 100,  -1 , 1. ) );
  TtBarPull_.push_back( fs->make<TH1D>(ns.name("TtBarRapidity"), "rapidity (t#bar{t})"   , 100,  -1 , 1. ) );
  TtBarPull_.push_back( fs->make<TH1D>(ns.name("TtBarPhi"     ), "#phi (t#bar{t})"	 , 100,  -1 , 1. ) );
  TtBarPull_.push_back( fs->make<TH1D>(ns.name("TtBarMass"    ), "M (t#bar{t}) [GeV]"	 , 100,  -1 , 1. ) );

  LepPairPull_.push_back( fs->make<TH1D>(ns.name("LepPairPt"	  ), "p_{t} (l^{+}l^{-}) [GeV]", 100,  -1 , 1. ) );
  LepPairPull_.push_back( fs->make<TH1D>(ns.name("LepPairEta"	  ), "#eta (l^{+}l^{-})"       , 100,  -1 , 1. ) );
  LepPairPull_.push_back( fs->make<TH1D>(ns.name("LepPairRapidity"), "rapidity (l^{+}l^{-})"   , 100,  -1 , 1. ) );
  LepPairPull_.push_back( fs->make<TH1D>(ns.name("LepPairPhi"	  ), "#phi (l^{+}l^{-})"       , 100,  -1 , 1. ) );
  LepPairPull_.push_back( fs->make<TH1D>(ns.name("LepPairMass"    ), "M (l^{+}l^{-}) [GeV]"    , 100,  -1 , 1. ) );

  JetPairPull_.push_back( fs->make<TH1D>(ns.name("JetPairPt"	  ), "p_{t} (jj) [GeV]", 100,  -1 , 1. ) );
  JetPairPull_.push_back( fs->make<TH1D>(ns.name("JetPairEta"	  ), "#eta (jj)"       , 100,  -1 , 1. ) );
  JetPairPull_.push_back( fs->make<TH1D>(ns.name("JetPairRapidity"), "rapidity (jj)"   , 100,  -1 , 1. ) );
  JetPairPull_.push_back( fs->make<TH1D>(ns.name("JetPairPhi"	  ), "#phi (jj)"       , 100,  -1 , 1. ) );
  JetPairPull_.push_back( fs->make<TH1D>(ns.name("JetPairMass"    ), "M (jj) [GeV]"    , 100,  -1 , 1. ) );

  TopPull2D_     = fs->make<TH2D>(ns.name("TopPull2D"),    "TopPull2D",     50, -1, 1, 50, -1, 1);
  WplusPull2D_   = fs->make<TH2D>(ns.name("WplusPull2D"),  "WplusPull2D",   50, -1, 1, 50, -1, 1);
  BPull2D_       = fs->make<TH2D>(ns.name("BPull2D"),      "BPull2D",       50, -1, 1, 50, -1, 1);
  LepBarPull2D_  = fs->make<TH2D>(ns.name("LepBarPull2D"), "LepBarPull2D",  50, -1, 1, 50, -1, 1);
  NuPull2D_      = fs->make<TH2D>(ns.name("NuPull2D"),     "NuPull2D",      50, -1, 1, 50, -1, 1);

  TopBarPull2D_  = fs->make<TH2D>(ns.name("TopBarPull2D"), "TopBarPull2D", 50, -1, 1, 50, -1, 1);
  WminusPull2D_  = fs->make<TH2D>(ns.name("WminusPull2D"), "WminusPull2D", 50, -1, 1, 50, -1, 1);
  BBarPull2D_    = fs->make<TH2D>(ns.name("BBarPull2D"),   "BBarPull2D",   50, -1, 1, 50, -1, 1);
  LepPull2D_     = fs->make<TH2D>(ns.name("LepPull2D"),    "LepPull2D",    50, -1, 1, 50, -1, 1);
  NuBarPull2D_   = fs->make<TH2D>(ns.name("NuBarPull2D"),  "NuBarPull2D",  50, -1, 1, 50, -1, 1);

  TtBarPull2D_     = fs->make<TH2D>(ns.name("TtBarPull2D"),     "TtBarPull2D", 50, -1, 1, 50, -1, 1);
  LepPairPull2D_ = fs->make<TH2D>(ns.name("LepPairPull2D"), "LepPairPull2D",   50, -1, 1, 50, -1, 1);
  JetPairPull2D_ = fs->make<TH2D>(ns.name("JetPairPull2D"), "JetPairPull2D",   50, -1, 1, 50, -1, 1);

}



/// book histograms for 2D distributions of particle properties: Pt, E, Eta, Phi, m
void
FullLepKinAnalyzer::book2DHistos(edm::Service<TFileService>& fs)
{

  NameScheme ns("2D");

  Top2D_.push_back( fs->make<TH2D>(ns.name("TopPt"	), "p_{t} (t) [GeV]",100,   0., 500. ,100,  0. , 500. ) );
  Top2D_.push_back( fs->make<TH2D>(ns.name("TopEta"	), "#eta (t)"	    ,100,  -5.,   5. ,100, -5. ,   5. ) );
  Top2D_.push_back( fs->make<TH2D>(ns.name("TopRapidity"), "rapidity (t)"   ,100,  -5.,   5. ,100, -5. ,   5. ) );
  Top2D_.push_back( fs->make<TH2D>(ns.name("TopPhi"	), "#phi (t)"	    , 62, -3.1,   3.1, 62, -3.1,   3.1) );
  Top2D_.push_back( fs->make<TH2D>(ns.name("TopMass"	), "M (top) [GeV]"  , 50,   0., 500.,  50,  0. , 500. ) );

  Wplus2D_.push_back( fs->make<TH2D>(ns.name("WplusPt"      ), "p_{t} (W^{+}) [GeV]",100,   0., 500. ,100,  0. , 500. ) );
  Wplus2D_.push_back( fs->make<TH2D>(ns.name("WplusEta"     ), "#eta (W^{+})"	    ,100,  -5.,   5. ,100, -5. ,   5. ) );
  Wplus2D_.push_back( fs->make<TH2D>(ns.name("WplusRapidity"), "rapidity (W^{+})"   ,100,  -5.,   5. ,100, -5. ,   5. ) );
  Wplus2D_.push_back( fs->make<TH2D>(ns.name("WplusPhi"     ), "#phi (W^{+})"	    , 62, -3.1,   3.1, 62, -3.1,   3.1) );

  B2D_.push_back( fs->make<TH2D>(ns.name("BPt"      ), "p_{t} (b) [GeV]",100,	0., 500. ,100,  0. , 500. ) );
  B2D_.push_back( fs->make<TH2D>(ns.name("BEta"     ), "#eta (b)"	,100,  -5.,   5. ,100, -5. ,   5. ) );
  B2D_.push_back( fs->make<TH2D>(ns.name("BRapidity"), "rapidity (b)"	,100,  -5.,   5. ,100, -5. ,   5. ) );
  B2D_.push_back( fs->make<TH2D>(ns.name("BPhi"     ), "#phi (b)"	, 62, -3.1,   3.1, 62, -3.1,   3.1) );

  LepBar2D_.push_back( fs->make<TH2D>(ns.name("LepBarPt"      ), "p_{t} (l^{+}) [GeV]",100,  0., 500. ,100,  0. , 500. ) );
  LepBar2D_.push_back( fs->make<TH2D>(ns.name("LepBarEta"     ), "#eta (l^{+})"       ,100, -5.,   5. ,100, -5. ,   5. ) );
  LepBar2D_.push_back( fs->make<TH2D>(ns.name("LepBarRapidity"), "rapidity (l^{+})"   ,100, -5.,   5. ,100, -5. ,   5. ) );
  LepBar2D_.push_back( fs->make<TH2D>(ns.name("LepBarPhi"     ), "#phi (l^{+})"       , 62, -3.1,  3.1, 62, -3.1,   3.1) );

  Nu2D_.push_back( fs->make<TH2D>(ns.name("NuPt"      ), "p_{t} (#nu^{-}) [GeV]",100,	0., 500. ,100,  0. , 500. ) );
  Nu2D_.push_back( fs->make<TH2D>(ns.name("NuEta"     ), "#eta (#nu^{-})"	,100,  -5.,   5. ,100, -5. ,   5. ) );
  Nu2D_.push_back( fs->make<TH2D>(ns.name("NuRapidity"), "rapidity (#nu^{-})"	,100,  -5.,   5. ,100, -5. ,   5. ) );
  Nu2D_.push_back( fs->make<TH2D>(ns.name("NuPhi"     ), "#phi (#nu^{-})"	, 62, -3.1,   3.1, 62, -3.1,   3.1) );

  TopBar2D_.push_back( fs->make<TH2D>(ns.name("TopBarPt"      ), "p_{t} (#bar{t}) [GeV]",100,	0., 500. ,100,  0. , 500. ) );
  TopBar2D_.push_back( fs->make<TH2D>(ns.name("TopBarEta"     ), "#eta (#bar{t})"	,100,  -5.,   5. ,100, -5. ,   5. ) );
  TopBar2D_.push_back( fs->make<TH2D>(ns.name("TopBarRapidity"), "rapidity (#bar{t})"	,100,  -5.,   5. ,100, -5. ,   5. ) );
  TopBar2D_.push_back( fs->make<TH2D>(ns.name("TopBarPhi"     ), "#phi (#bar{t})"	, 62, -3.1,   3.1, 62, -3.1,   3.1) );
  TopBar2D_.push_back( fs->make<TH2D>(ns.name("TopBarMass"    ), "M (#bar{t}) [GeV]"	, 50,	0., 500.,  50,  0. , 500. ) );

  Wminus2D_.push_back( fs->make<TH2D>(ns.name("WminusPt"      ), "p_{t} (W^{-}) [GeV]",100,  0., 500. ,100,  0. , 500. ) );
  Wminus2D_.push_back( fs->make<TH2D>(ns.name("WminusEta"     ), "#eta (W^{-})"       ,100, -5.,   5. ,100, -5. ,   5. ) );
  Wminus2D_.push_back( fs->make<TH2D>(ns.name("WminusRapidity"), "rapidity (W^{-})"   ,100, -5.,   5. ,100, -5. ,   5. ) );
  Wminus2D_.push_back( fs->make<TH2D>(ns.name("WminusPhi"     ), "#phi (W^{-})"       , 62, -3.1,  3.1, 62, -3.1,   3.1) );

  BBar2D_.push_back( fs->make<TH2D>(ns.name("BBarPt"	  ), "p_{t} (#bar{b}) [GeV]",100,  0., 500. ,100,  0. , 500. ) );
  BBar2D_.push_back( fs->make<TH2D>(ns.name("BBarEta"	  ), "#eta (#bar{b})"	    ,100, -5.,   5. ,100, -5. ,   5. ) );
  BBar2D_.push_back( fs->make<TH2D>(ns.name("BBarRapidity"), "rapidity (#bar{b})"   ,100, -5.,   5. ,100, -5. ,   5. ) );
  BBar2D_.push_back( fs->make<TH2D>(ns.name("BBarPhi"	  ), "#phi (#bar{b})"	    , 62, -3.1,  3.1, 62, -3.1,   3.1) );

  Lep2D_.push_back( fs->make<TH2D>(ns.name("LepPt"	), "p_{t} (l^{-}) [GeV]",100,  0., 500. ,100,  0. , 500. ) );
  Lep2D_.push_back( fs->make<TH2D>(ns.name("LepEta"	), "#eta (l^{-})"	,100, -5.,   5. ,100, -5. ,   5. ) );
  Lep2D_.push_back( fs->make<TH2D>(ns.name("LepRapidity"), "rapidity (l^{-})"	,100, -5.,   5. ,100, -5. ,   5. ) );
  Lep2D_.push_back( fs->make<TH2D>(ns.name("LepPhi"	), "#phi (l^{-})"	, 62, -3.1,  3.1, 62, -3.1,   3.1) );

  NuBar2D_.push_back( fs->make<TH2D>(ns.name("NuBarPt"      ), "p_{t} (#bar{#nu}^{-}) [GeV]",100,  0., 500. ,100,  0. , 500. ) );
  NuBar2D_.push_back( fs->make<TH2D>(ns.name("NuBarEta"     ), "#eta (#bar{#nu}^{-})"	    ,100, -5.,   5. ,100, -5. ,   5. ) );
  NuBar2D_.push_back( fs->make<TH2D>(ns.name("NuBarRapidity"), "rapidity (#bar{#nu}^{-})"   ,100, -5.,   5. ,100, -5. ,   5. ) );
  NuBar2D_.push_back( fs->make<TH2D>(ns.name("NuBarPhi"     ), "#phi (#bar{#nu}^{-})"	    , 62, -3.1,  3.1, 62, -3.1,   3.1) );

  TtBar2D_.push_back( fs->make<TH2D>(ns.name("TtBarPt"      ), "p_{t} (t#bar{t}) [GeV]",100,  0., 500. ,100,  0. , 500. ) );
  TtBar2D_.push_back( fs->make<TH2D>(ns.name("TtBarEta"     ), "#eta (t#bar{t})"       ,100, -5.,   5. ,100, -5. ,   5. ) );
  TtBar2D_.push_back( fs->make<TH2D>(ns.name("TtBarRapidity"), "rapidity (t#bar{t})"   ,100, -5.,   5. ,100, -5. ,   5. ) );
  TtBar2D_.push_back( fs->make<TH2D>(ns.name("TtBarPhi"     ), "#phi (t#bar{t})"       , 62, -3.1,  3.1, 62, -3.1,   3.1) );
  TtBar2D_.push_back( fs->make<TH2D>(ns.name("TtBarMass"    ), "M (t#bar{t}) [GeV]"    ,100,  0.,1000. ,100,  0. ,1000. ) );

  LepPair2D_.push_back( fs->make<TH2D>(ns.name("LepPairPt"	), "p_{t} (l^{+}l^{-}) [GeV]",100,  0., 500. ,100,  0. , 500. ) );
  LepPair2D_.push_back( fs->make<TH2D>(ns.name("LepPairEta"	), "#eta (l^{+}l^{-})"       ,100, -5.,   5. ,100, -5. ,   5. ) );
  LepPair2D_.push_back( fs->make<TH2D>(ns.name("LepPairRapidity"), "rapidity (l^{+}l^{-})"   ,100, -5.,   5. ,100, -5. ,   5. ) );
  LepPair2D_.push_back( fs->make<TH2D>(ns.name("LepPairPhi"	), "#phi (l^{+}l^{-})"       , 62, -3.1,  3.1, 62, -3.1,   3.1) );
  LepPair2D_.push_back( fs->make<TH2D>(ns.name("LepPairMass"	), "M (l^{+}l^{-}) [GeV]"    ,200,  0., 400. ,200,  0. , 400. ) );

  JetPair2D_.push_back( fs->make<TH2D>(ns.name("JetPairPt"	), "p_{t} (jj) [GeV]",100,  0., 500. ,100,  0. , 500. ) );
  JetPair2D_.push_back( fs->make<TH2D>(ns.name("JetPairEta"	), "#eta (jj)"       ,100, -5.,   5. ,100, -5. ,   5. ) );
  JetPair2D_.push_back( fs->make<TH2D>(ns.name("JetPairRapidity"), "rapidity (jj)"   ,100, -5.,   5. ,100, -5. ,   5. ) );
  JetPair2D_.push_back( fs->make<TH2D>(ns.name("JetPairPhi"	), "#phi (jj)"       , 62, -3.1,  3.1, 62, -3.1,   3.1) );
  JetPair2D_.push_back( fs->make<TH2D>(ns.name("JetPairMass"	), "M (jj) [GeV]"    , 50,  0., 500. , 50,  0. , 500. ) );
}


/// book histograms for 2D correlations between reconstructed quantities
void
FullLepKinAnalyzer::bookRecCorrelHistos(edm::Service<TFileService>& fs)
{
  NameScheme ns("correlation");
  diLeptonMassVsLeptonsPt_    = fs->make<TH2D>(ns.name("diLeptonMassVsLeptonsPt"),    "", 200,  0. , 400. , 250,  0. , 500. );
  diLeptonMassVsLeptonsEta_   = fs->make<TH2D>(ns.name("diLeptonMassVsLeptonsEta"),   "", 200,  0. , 400. , 100, -5.0,   5.0);
  diLeptonMassVsDiLeptonPt_   = fs->make<TH2D>(ns.name("diLeptonMassVsDiLeptonPt"),   "", 200,  0. , 400. , 250,  0. , 500. );
  diLeptonMassVsJetsPt_       = fs->make<TH2D>(ns.name("diLeptonMassVsJetsPt"),       "", 200,  0. , 400. , 250,  0. , 500. );
  diLeptonMassVsJetsEta_      = fs->make<TH2D>(ns.name("diLeptonMassVsJetsEta"),      "", 200,  0. , 400. , 100, -5.0,   5.0);  
  diLeptonMassVsTopPt_        = fs->make<TH2D>(ns.name("diLeptonMassVsTopPt"),        "", 200,  0. , 400. , 250,  0. , 500. );
  diLeptonMassVsTopRapidity_  = fs->make<TH2D>(ns.name("diLeptonMassVsTopRapidity"),  "", 200,  0. , 400. , 100, -5.0,   5.0);
  diLeptonMassVsTtBarPt_      = fs->make<TH2D>(ns.name("diLeptonMassVsTtBarPt"),      "", 200,  0. , 400. , 250,  0. , 500. );     
  diLeptonMassVsTtBarRapidity_= fs->make<TH2D>(ns.name("diLeptonMassVsTtBarRapidity"),"", 200,  0. , 400. , 100, -5.0,   5.0);
  diLeptonMassVsTtBarMass_    = fs->make<TH2D>(ns.name("diLeptonMassVsTtBarMass"),    "", 200,  0. , 400. , 500,  0. ,1000. );

  nBtagsTCHELVsLeptonsPt_    = fs->make<TH2D>(ns.name("nBtagsTCHELVsLeptonsPt"),      "", 3, -0.5, 2.5 , 250,  0. , 500. );
  nBtagsTCHELVsLeptonsEta_   = fs->make<TH2D>(ns.name("nBtagsTCHELVsLeptonsEta"),     "", 3, -0.5, 2.5 , 100, -5.0,   5.0);
  nBtagsTCHELVsDiLeptonPt_   = fs->make<TH2D>(ns.name("nBtagsTCHELVsDiLeptonPt"),     "", 3, -0.5, 2.5 , 250,  0. , 500. );  
  nBtagsTCHELVsDiLeptonMass_ = fs->make<TH2D>(ns.name("nBtagsTCHELVsDiLeptonMass"),   "", 3, -0.5, 2.5 , 250,  0. , 500. );
  nBtagsTCHELVsJetsPt_       = fs->make<TH2D>(ns.name("nBtagsTCHELVsJetsPt"),         "", 3, -0.5, 2.5 , 250,  0. , 500. );
  nBtagsTCHELVsJetsEta_      = fs->make<TH2D>(ns.name("nBtagsTCHELVsJetsEta"),        "", 3, -0.5, 2.5 , 100, -5.0,   5.0);  
  nBtagsTCHELVsTopPt_        = fs->make<TH2D>(ns.name("nBtagsTCHELVsTopPt"),          "", 3, -0.5, 2.5 , 250,  0. , 500. );
  nBtagsTCHELVsTopRapidity_  = fs->make<TH2D>(ns.name("nBtagsTCHELVsTopRapidity"),    "", 3, -0.5, 2.5 , 100, -5.0,   5.0);
  nBtagsTCHELVsTtBarPt_      = fs->make<TH2D>(ns.name("nBtagsTCHELVsTtBarPt"),        "", 3, -0.5, 2.5 , 250,  0. , 500. );   
  nBtagsTCHELVsTtBarRapidity_= fs->make<TH2D>(ns.name("nBtagsTCHELVsTtBarRapidity"),  "", 3, -0.5, 2.5 , 100, -5.0,   5.0);
  nBtagsTCHELVsTtBarMass_    = fs->make<TH2D>(ns.name("nBtagsTCHELVsTtBarMass"),      "", 3, -0.5, 2.5 , 500,  0. ,1000. );

  nBtagsTCHEMVsLeptonsPt_    = fs->make<TH2D>(ns.name("nBtagsTCHEMVsLeptonsPt"),      "", 3, -0.5, 2.5 , 250,  0. , 500. );
  nBtagsTCHEMVsLeptonsEta_   = fs->make<TH2D>(ns.name("nBtagsTCHEMVsLeptonsEta"),     "", 3, -0.5, 2.5 , 100, -5.0,   5.0);
  nBtagsTCHEMVsDiLeptonPt_   = fs->make<TH2D>(ns.name("nBtagsTCHEMVsDiLeptonPt"),     "", 3, -0.5, 2.5 , 250,  0. , 500. );
  nBtagsTCHEMVsDiLeptonMass_ = fs->make<TH2D>(ns.name("nBtagsTCHEMVsDiLeptonMass"),   "", 3, -0.5, 2.5 , 250,  0. , 500. );
  nBtagsTCHEMVsJetsPt_       = fs->make<TH2D>(ns.name("nBtagsTCHEMVsJetsPt"),         "", 3, -0.5, 2.5 , 250,  0. , 500. );
  nBtagsTCHEMVsJetsEta_      = fs->make<TH2D>(ns.name("nBtagsTCHEMVsJetsEta"),        "", 3, -0.5, 2.5 , 100, -5.0,   5.0);   
  nBtagsTCHEMVsTopPt_        = fs->make<TH2D>(ns.name("nBtagsTCHEMVsTopPt"),          "", 3, -0.5, 2.5 , 250,  0. , 500. );
  nBtagsTCHEMVsTopRapidity_  = fs->make<TH2D>(ns.name("nBtagsTCHEMVsTopRapidity"),    "", 3, -0.5, 2.5 , 100, -5.0,   5.0);
  nBtagsTCHEMVsTtBarPt_      = fs->make<TH2D>(ns.name("nBtagsTCHEMVsTtBarPt"),        "", 3, -0.5, 2.5 , 250,  0. , 500. );	
  nBtagsTCHEMVsTtBarRapidity_= fs->make<TH2D>(ns.name("nBtagsTCHEMVsTtBarRapidity"),  "", 3, -0.5, 2.5 , 100, -5.0,   5.0);
  nBtagsTCHEMVsTtBarMass_    = fs->make<TH2D>(ns.name("nBtagsTCHEMVsTtBarMass"),      "", 3, -0.5, 2.5 , 500,  0. ,1000. );

  nBtagsSSVHEMVsLeptonsPt_    = fs->make<TH2D>(ns.name("nBtagsSSVHEMVsLeptonsPt"),    "", 3, -0.5, 2.5 , 250,  0. , 500. );
  nBtagsSSVHEMVsLeptonsEta_   = fs->make<TH2D>(ns.name("nBtagsSSVHEMVsLeptonsEta"),   "", 3, -0.5, 2.5 , 100, -5.0,   5.0);
  nBtagsSSVHEMVsDiLeptonPt_   = fs->make<TH2D>(ns.name("nBtagsSSVHEMVsDiLeptonPt"),   "", 3, -0.5, 2.5 , 250,  0. , 500. );
  nBtagsSSVHEMVsDiLeptonMass_ = fs->make<TH2D>(ns.name("nBtagsSSVHEMVsDiLeptonMass"), "", 3, -0.5, 2.5 , 250,  0. , 500. );
  nBtagsSSVHEMVsJetsPt_       = fs->make<TH2D>(ns.name("nBtagsSSVHEMVsJetsPt"),        "", 3, -0.5, 2.5 , 250,  0. , 500. );
  nBtagsSSVHEMVsJetsEta_      = fs->make<TH2D>(ns.name("nBtagsSSVHEMVsJetsEta"),       "", 3, -0.5, 2.5 , 100, -5.0,   5.0);   
  nBtagsSSVHEMVsTopPt_        = fs->make<TH2D>(ns.name("nBtagsSSVHEMVsTopPt"),        "", 3, -0.5, 2.5 , 250,  0. , 500. );
  nBtagsSSVHEMVsTopRapidity_  = fs->make<TH2D>(ns.name("nBtagsSSVHEMVsTopRapidity"),  "", 3, -0.5, 2.5 , 100, -5.0,   5.0);
  nBtagsSSVHEMVsTtBarPt_      = fs->make<TH2D>(ns.name("nBtagsSSVHEMVsTtBarPt"),      "", 3, -0.5, 2.5 , 250,  0. , 500. );  
  nBtagsSSVHEMVsTtBarRapidity_= fs->make<TH2D>(ns.name("nBtagsSSVHEMVsTtBarRapidity"),"", 3, -0.5, 2.5 , 100, -5.0,   5.0);
  nBtagsSSVHEMVsTtBarMass_    = fs->make<TH2D>(ns.name("nBtagsSSVHEMVsTtBarMass"),    "", 3, -0.5, 2.5 , 500,  0. ,1000. );
  
  nMatchesVsLeptonsPt_	  = fs->make<TH2D>(ns.name("nMatchesVsLeptonsPt"),	"", 3, -0.5, 2.5 , 250,  0. , 500. );
  nMatchesVsLeptonsEta_   = fs->make<TH2D>(ns.name("nMatchesVsLeptonsEta"),	"", 3, -0.5, 2.5 , 100, -5.0,	5.0);
  nMatchesVsDiLeptonPt_   = fs->make<TH2D>(ns.name("nMatchesVsDiLeptonPt"),	"", 3, -0.5, 2.5 , 250,  0. , 500. );
  nMatchesVsDiLeptonMass_ = fs->make<TH2D>(ns.name("nMatchesVsDiLeptonMass"),	"", 3, -0.5, 2.5 , 250,  0. , 500. );
  nMatchesVsJetsPt_	  = fs->make<TH2D>(ns.name("nMatchesVsJetsPt"),	        "", 3, -0.5, 2.5 , 250,  0. , 500. );
  nMatchesVsJetsEta_      = fs->make<TH2D>(ns.name("nMatchesVsJetsEta"),	"", 3, -0.5, 2.5 , 100, -5.0,	5.0);  
  nMatchesVsTopPt_	  = fs->make<TH2D>(ns.name("nMatchesVsTopPt"),  	"", 3, -0.5, 2.5 , 250,  0. , 500. );
  nMatchesVsTopRapidity_  = fs->make<TH2D>(ns.name("nMatchesVsTopRapidity"),	"", 3, -0.5, 2.5 , 100, -5.0,	5.0);
  nMatchesVsTtBarPt_	  = fs->make<TH2D>(ns.name("nMatchesVsTtBarPt"),	"", 3, -0.5, 2.5 , 250,  0. , 500. );	
  nMatchesVsTtBarRapidity_= fs->make<TH2D>(ns.name("nMatchesVsTtBarRapidity"),  "", 3, -0.5, 2.5 , 100, -5.0,	5.0);
  nMatchesVsTtBarMass_    = fs->make<TH2D>(ns.name("nMatchesVsTtBarMass"),	"", 3, -0.5, 2.5 , 500,  0. ,1000. );


  nMatchVsNtagsTCHEL_ = fs->make<TH2D>(ns.name("nMatchVsNtagsTCHEL" ),"", 5, -0.5, 4.5, 3, -0.5, 2.5);  
  nMatchVsNtagsTCHEM_ = fs->make<TH2D>(ns.name("nMatchVsNtagsTCHEM" ),"", 5, -0.5, 4.5, 3, -0.5, 2.5);   
  nMatchVsNtagsSSVHEM_= fs->make<TH2D>(ns.name("nMatchVsNtagsSSVHEM"),"", 5, -0.5, 4.5, 3, -0.5, 2.5);    
}


/// book histograms for hypothesis specific histos and correlations between hypotheses
void
FullLepKinAnalyzer::bookQualityHistos(edm::Service<TFileService>& fs)
{
  NameScheme ns("qual");

  kinSolWeight_        = fs->make<TH1D>(ns.name("kinSolWeight"        ), "Weight of kin solution",           50,  0. , 1. ); 
  bJetIdcs_            = fs->make<TH1D>(ns.name("bJetIdcs"            ), "b jet indices used for hypo",       4, -0.5, 3.5);
  bBarJetIdcs_         = fs->make<TH1D>(ns.name("bBarJetIdcs"         ), "bbar jet indices used for hypo",    4, -0.5, 3.5);
  nCorrectAssignments_ = fs->make<TH1D>(ns.name("nCorrectAssignments" ), "N correctly assigned jets",         3, -0.5, 2.5);
  deltaM_              = fs->make<TH1D>(ns.name("deltaM"              ), "M_{top}",		             50, -25., 25.);
  LepBarBMass_         = fs->make<TH1D>(ns.name("LepBarBMass"	      ), "M_{l^{+},b}",		            250,   0.,250.);  
  LepBMass_            = fs->make<TH1D>(ns.name("LepBMass"	      ), "M_{l^{-},b}",		            250,   0.,250.);  
  LepBBarMass_         = fs->make<TH1D>(ns.name("LepBBarMass"	      ), "M_{l^{-},#bar{b}}",	            250,   0.,250.);  
  LepBarBBarMass_      = fs->make<TH1D>(ns.name("LepBarBBarMass"      ), "M_{l^{+},#bar{b}}",	            250,   0.,250.);  

  kinTCHELcorrelation_ = fs->make<TH2D>(ns.name("kinTCHELcorrelation" ), "mass reco vs. TCHEL",        3, -0.5, 2.5, 2, -0.5, 1.5);
  kinTCHEMcorrelation_ = fs->make<TH2D>(ns.name("kinTCHEMcorrelation" ), "mass reco vs. TCHEM",        3, -0.5, 2.5, 2, -0.5, 1.5);
  kinSSVHEMcorrelation_= fs->make<TH2D>(ns.name("kinSSVHEMcorrelation"), "mass reco vs. TCHEM",        3, -0.5, 2.5, 2, -0.5, 1.5);
  compare_             = fs->make<TH2D>(ns.name("compare"             ), "Indices",                   24, -1.5,22.5,24, -1.5,22.5);

  for(int i=1; i<=24;++i){
    char label[2];
    int j = i%4-2;

    std::ostringstream oss;
    oss << j;
    strcpy(label, oss.str().c_str());

    compare_->GetXaxis()->SetBinLabel(i,label);
    compare_->GetYaxis()->SetBinLabel(i,label);
  }
}


/// fill histograms for reconstructed or generated particle properties: Pt, E, Eta, Phi, m
void
FullLepKinAnalyzer::fillKinHistos(std::vector<TH1D*>& histos, const reco::Candidate& candidate, double weight)
{
  histos[0]->Fill( candidate.pt(), weight);
  histos[1]->Fill( candidate.eta(), weight);
  histos[2]->Fill( candidate.rapidity(), weight);
  histos[3]->Fill( candidate.phi(), weight);
  histos[4]->Fill( candidate.mass(), weight);
}


/// fill pull histograms: Pt, E, Eta, Phi, m
void
FullLepKinAnalyzer::fillPullHistos(std::vector<TH1D*>& histos, TH2D& hist2D, const reco::Candidate& candidate,  const reco::Candidate& gencandidate, double weight)
{
  histos[0]->Fill( ( candidate.pt()      - gencandidate.pt() )       / gencandidate.pt(),       weight);
  histos[1]->Fill( ( candidate.eta()     - gencandidate.eta() )      / gencandidate.eta(),      weight);
  histos[2]->Fill( ( candidate.rapidity()- gencandidate.rapidity() ) / gencandidate.rapidity(), weight);
  histos[3]->Fill( ( candidate.phi()     - gencandidate.phi() )      / gencandidate.phi(),      weight);
  histos[4]->Fill( ( candidate.mass()    - gencandidate.mass() )     / gencandidate.mass(),     weight);

  hist2D.Fill((candidate.eta()-gencandidate.eta())/gencandidate.eta(),(candidate.phi()-gencandidate.phi())/gencandidate.phi(), weight);
}


/// fill 2D histograms
void
FullLepKinAnalyzer::fill2DHistos(std::vector<TH2D*>& histos, const reco::Candidate& candidate,  const reco::Candidate& gencandidate, double weight)
{
    histos[0]->Fill( gencandidate.pt(),       candidate.pt() , weight      );
    histos[1]->Fill( gencandidate.eta(),      candidate.eta(), weight      );
    histos[2]->Fill( gencandidate.rapidity(), candidate.rapidity(), weight );
    histos[3]->Fill( gencandidate.phi(),      candidate.phi(), weight      );
    if (histos.size() > 4) {
        histos[4]->Fill( gencandidate.mass(), candidate.mass(), weight     );
    }
}


/// fill histograms for hypothesis specific histos and correlations between hypotheses
void
FullLepKinAnalyzer::fillQualityHistos(const TtFullLeptonicEvent& FullLepEvt, const TtEvent::HypoClassKey& hypoKey, double weight)
{
  bJetIdcs_   ->Fill(FullLepEvt.jetLeptonCombination(hypoKey)[0], weight );
  bBarJetIdcs_->Fill(FullLepEvt.jetLeptonCombination(hypoKey)[1], weight );

  if( hypoKey==6 ) {
    deltaM_->Fill(FullLepEvt.top(hypoKey)->mass()-FullLepEvt.topBar(hypoKey)->mass(), weight);
    kinSolWeight_->Fill( FullLepEvt.solWeight(), weight);
  }
}


bool 
FullLepKinAnalyzer::isRealBJet(const reco::Candidate* jet, const reco::Candidate* genB, const reco::Candidate* genBBar)
{
  double dist1 = deltaR(jet->eta(), jet->phi(), genB->eta(),    genB->phi()   );
  double dist2 = 999.;
  if(genBBar)
    dist2 = deltaR(jet->eta(), jet->phi(), genBBar->eta(), genBBar->phi());

  if(dist1 < .2 || dist2 < .2){  
    return true;
  }
  return false;
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( FullLepKinAnalyzer );
