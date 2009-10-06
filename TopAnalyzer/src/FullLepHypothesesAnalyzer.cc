#include "TopAnalysis/TopAnalyzer/interface/FullLepHypothesesAnalyzer.h"
#include "TopAnalysis/TopUtils/interface/NameScheme.h"
#include "AnalysisDataFormats/TopObjects/interface/TtEventPartons.h"

FullLepHypothesesAnalyzer::FullLepHypothesesAnalyzer(const edm::ParameterSet& cfg):
  FullLepEvt_      (cfg.getParameter<edm::InputTag>("FullLepEvent"    )),
  hypoKey_         (cfg.getParameter<edm::InputTag>("hypoKey"         )),
  wgt_             (cfg.getParameter<edm::InputTag>("weight"          )),
  useWrongCharge_  (cfg.getParameter<bool>         ("alsoWrongCharge" )),
  wantSummary_     (cfg.getParameter<bool>         ("wantSummary"     )),  
  hist_            (cfg.getParameter<std::string>  ("hist"            ))
{
  Nava = 0;
  Nval = 0;
}

void
FullLepHypothesesAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{

  edm::Handle<TtFullLeptonicEvent> FullLepEvt;
  evt.getByLabel(FullLepEvt_, FullLepEvt);

  edm::Handle<int> hypoKeyHandle;
  evt.getByLabel(hypoKey_, hypoKeyHandle);
  TtEvent::HypoClassKey& hypoKey = (TtEvent::HypoClassKey&) *hypoKeyHandle;

  edm::Handle<double> wgt;
  evt.getByLabel(wgt_, wgt);
  double weight = *wgt;

  // -----------------------
  // check if hypothesis is valid in this event
  // and if it satisfies some quality criteria
  // -----------------------
         
  if( !FullLepEvt->isHypoAvailable(hypoKey) ){
    edm::LogInfo ( "NonValidHyp" ) << "Hypothesis not available for this event";
    return;
  }
  else Nava++;
  
  if( !FullLepEvt->isHypoValid(hypoKey) )  {
    goodHypo_->Fill(0., weight); // not a good hypothesis
  }
  else goodHypo_->Fill(1., weight); // good hypothesis

   if( !FullLepEvt->isHypoValid(hypoKey) ){
     edm::LogInfo ( "NonValidHyp" ) << "Hypothesis not valid for this event";
     return;  // return if any of the hypotheses is not valid
   }
  else Nval++;  

  // -----------------------
  // fill histos related to quality of the TtFullLeptonicEvent
  // -----------------------
  fillQualityHistos(*FullLepEvt, hypoKey, weight);

  // -----------------------
  // fill histos for basic kinematic variables
  // -----------------------
  const reco::Candidate* Top    = FullLepEvt->top(hypoKey);
  const reco::Candidate* Wplus  = FullLepEvt->wPlus(hypoKey);
  const reco::Candidate* B      = FullLepEvt->b(hypoKey);
  const reco::Candidate* LepBar = FullLepEvt->leptonBar(hypoKey);
  const reco::Candidate* Nu     = FullLepEvt->neutrino(hypoKey);
  
  const reco::Candidate* TopBar = FullLepEvt->topBar(hypoKey);
  const reco::Candidate* Wminus = FullLepEvt->wMinus(hypoKey);
  const reco::Candidate* BBar   = FullLepEvt->bBar(hypoKey);
  const reco::Candidate* Lep    = FullLepEvt->lepton(hypoKey);
  const reco::Candidate* NuBar  = FullLepEvt->neutrinoBar(hypoKey); 
  
  if(!FullLepEvt->isWrongCharge()){
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
  }
  else{
    fillWrongChargeKinHistos(TopKin_,    *Top,    weight);
    fillWrongChargeKinHistos(WplusKin_,  *Wplus,  weight);
    fillWrongChargeKinHistos(BKin_,      *B,      weight);
    fillWrongChargeKinHistos(LepBarKin_, *LepBar, weight);
    fillWrongChargeKinHistos(NuKin_,     *Nu,     weight);

    fillWrongChargeKinHistos(TopBarKin_, *TopBar, weight);
    fillWrongChargeKinHistos(WminusKin_, *Wminus, weight);
    fillWrongChargeKinHistos(BBarKin_,   *BBar,   weight);
    fillWrongChargeKinHistos(LepKin_,    *Lep,    weight);
    fillWrongChargeKinHistos(NuBarKin_,  *NuBar,  weight);  
  }
  // -----------------------
  // fill resolution histos for kinematic variables
  // with respect to the generator particles
  // -----------------------
  if( FullLepEvt->genEvent()->isFullLeptonic() ) {   
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
    
    if(genTop)    fillKinGenHistos(TopKinGen_,    *genTop,    weight);
    if(genWplus)  fillKinGenHistos(WplusKinGen_,  *genWplus,  weight);
    if(genB)      fillKinGenHistos(BKinGen_,      *genB,      weight);
    if(genLepBar) fillKinGenHistos(LepBarKinGen_, *genLepBar, weight);
    if(genNu)     fillKinGenHistos(NuKinGen_,     *genNu,     weight);   
    
    if(genTopBar) fillKinGenHistos(TopBarKinGen_, *genTopBar, weight);
    if(genWminus) fillKinGenHistos(WminusKinGen_, *genWminus, weight);
    if(genBBar)   fillKinGenHistos(BBarKinGen_,   *genBBar,   weight);
    if(genLep)    fillKinGenHistos(LepKinGen_,    *genLep,    weight);
    if(genNuBar)  fillKinGenHistos(NuBarKinGen_,  *genNuBar,  weight);      

    if(!FullLepEvt->isWrongCharge()){
      if(genTop)    fillKinResHistos(TopKinRes_,    *Top,   *genTop,    weight);
      if(genWplus)  fillKinResHistos(WplusKinRes_,  *Wplus, *genWplus,  weight);
      if(genB)      fillKinResHistos(BKinRes_,      *B,     *genB,      weight);
      if(genLepBar) fillKinResHistos(LepBarKinRes_, *LepBar,*genLepBar, weight);
      if(genNu)     fillKinResHistos(NuKinRes_,     *Nu,    *genNu,     weight);

      if(genTopBar) fillKinResHistos(TopBarKinRes_, *TopBar,*genTopBar, weight);
      if(genWminus) fillKinResHistos(WminusKinRes_, *Wminus,*genWminus, weight);
      if(genBBar)   fillKinResHistos(BBarKinRes_,   *BBar,  *genBBar,   weight);
      if(genLep)    fillKinResHistos(LepKinRes_,    *Lep,   *genLep,    weight);
      if(genNuBar)  fillKinResHistos(NuBarKinRes_,  *NuBar, *genNuBar,  weight); 
    }
    else{ //should never be the case: right charge reconstructed as wrong charge
      if(genTop)    fillWrongChargeKinResHistos(TopKinRes_,    *Top,   *genTop,    weight);
      if(genWplus)  fillWrongChargeKinResHistos(WplusKinRes_,  *Wplus, *genWplus,  weight);
      if(genB)      fillWrongChargeKinResHistos(BKinRes_,      *B,     *genB,      weight);
      if(genLepBar) fillWrongChargeKinResHistos(LepBarKinRes_, *LepBar,*genLepBar, weight);
      if(genNu)     fillWrongChargeKinResHistos(NuKinRes_,     *Nu,    *genNu,     weight);

      if(genTopBar) fillWrongChargeKinResHistos(TopBarKinRes_, *TopBar,*genTopBar, weight);
      if(genWminus) fillWrongChargeKinResHistos(WminusKinRes_, *Wminus,*genWminus, weight);
      if(genBBar)   fillWrongChargeKinResHistos(BBarKinRes_,   *BBar,  *genBBar,   weight);
      if(genLep)    fillWrongChargeKinResHistos(LepKinRes_,    *Lep,   *genLep,    weight);
      if(genNuBar)  fillWrongChargeKinResHistos(NuBarKinRes_,  *NuBar, *genNuBar,  weight);   
    }
  }      
}

void 
FullLepHypothesesAnalyzer::beginJob(const edm::EventSetup&)
{
  if( hist_.empty() ) throw edm::Exception( edm::errors::Configuration, "Empty string for hist in cfi file" );
  ofstream hist(hist_.c_str(), std::ios::out);

  edm::Service<TFileService> fs;
  if( !fs ) throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );

  bookKinHistos      (fs, hist);
  bookKinGenHistos   (fs, hist);  
  bookKinResHistos   (fs, hist);
  bookQualityHistos  (fs, hist);
}

void
FullLepHypothesesAnalyzer::endJob() 
{
  if(!wantSummary_) return;
  
  std::cout << "++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "Summary for Hypothese " << hypoKey_ << std::endl;
  std::cout << "number of available events is " << Nava << std::endl;
  std::cout << "number of valid events is " << Nval << std::endl;
  std::cout << "++++++++++++++++++++++++++++++++" << std::endl;
}

void
FullLepHypothesesAnalyzer::bookKinHistos(edm::Service<TFileService>& fs, ofstream& hist)
{

  NameScheme ns("kin");

  TopKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopPt"    ), "p_{t} (t) [GeV]", 50,  0. , 500. ) );
  TopKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopEnergy"), "E (t) [GeV]"    , 50,  0. , 500. ) );  
  TopKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopEta"   ), "#eta (t)"       , 34, -3.4,   3.4) );
  TopKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopPhi"   ), "#phi (t)"       , 34, -3.4,   3.4) );
  TopKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopMass"  ), "M (top) [GeV]"  , 30, 50. , 350. ) );

  WplusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WplusPt"    ), "p_{t} (W^{+}) [GeV]", 50,  0. , 500. ) );
  WplusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WplusEnergy"), "E (W^{+}) [GeV]"    , 50,  0. , 500. ) );
  WplusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WplusEta"   ), "#eta (W^{+})"       , 34, -3.4,   3.4) );
  WplusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WplusPhi"   ), "#phi (W^{+})"       , 34, -3.4,   3.4) );
  WplusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WplusMass"  ), "M (W^{+}) [GeV]"    , 30,  0. , 600. ) );

  BKin_.push_back( fs->make<TH1F>(ns.name(hist, "BPt"    ), "p_{t} (b) [GeV]", 50,  0. , 500. ) );
  BKin_.push_back( fs->make<TH1F>(ns.name(hist, "BEnergy"), "E (b) [GeV]"    , 50,  0. , 500. ) );
  BKin_.push_back( fs->make<TH1F>(ns.name(hist, "BEta"   ), "#eta (b)"       , 34, -3.4,   3.4) );
  BKin_.push_back( fs->make<TH1F>(ns.name(hist, "BPhi"   ), "#phi (b)"       , 34, -3.4,   3.4) );
  BKin_.push_back( fs->make<TH1F>(ns.name(hist, "BMass"  ), "M (b) [GeV]"    , 30,  0. , 600. ) );
    
  LepBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarPt"    ), "p_{t} (l^{+}) [GeV]", 50,  0. , 500. ) );
  LepBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarEnergy"), "E (l^{+}) [GeV]"    , 50,  0. , 500. ) );
  LepBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarEta"   ), "#eta (l^{+})"       , 34, -3.4,   3.4) );
  LepBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarPhi"   ), "#phi (l^{+})"       , 34, -3.4,   3.4) );
  LepBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarMass"  ), "M (l^{+}) [GeV]"    , 30,  0. , 600. ) ); 

  NuKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuPt"    ), "p_{t} (#nu^{-}) [GeV]", 50,  0. , 500. ) );
  NuKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuEnergy"), "E (#nu^{-}) [GeV]"    , 50,  0. , 500. ) );
  NuKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuEta"   ), "#eta (#nu^{-})"       , 34, -3.4,   3.4) );
  NuKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuPhi"   ), "#phi (#nu^{-})"       , 34, -3.4,   3.4) );
  NuKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuMass"  ), "M (#nu^{-}) [GeV]"    , 30,  0. , 600. ) );

  TopBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarPt"    ), "p_{t} (#bar{t}) [GeV]", 50,  0. , 500. ) );
  TopBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarEnergy"), "E (#bar{t}) [GeV]"    , 50,  0. , 500. ) );
  TopBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarEta"   ), "#eta (#bar{t})"       , 34, -3.4,   3.4) );
  TopBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarPhi"   ), "#phi (#bar{t})"       , 34, -3.4,   3.4) );
  TopBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarMass"  ), "M (#bar{t}) [GeV]"    , 30, 50. , 350. ) );

  WminusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WminusPt"    ), "p_{t} (W^{-}) [GeV]", 50,  0. , 500. ) );
  WminusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WminusEnergy"), "E (W^{-}) [GeV]"    , 50,  0. , 500. ) );
  WminusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WminusEta"   ), "#eta (W^{-})"       , 34, -3.4,   3.4) );
  WminusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WminusPhi"   ), "#phi (W^{-})"       , 34, -3.4,   3.4) );
  WminusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WminusMass"  ), "M (W^{-}) [GeV]"    , 30,  0. , 600. ) );

  BBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "BBarPt"    ), "p_{t} (#bar{b}) [GeV]", 50,  0. , 500. ) );
  BBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "BBarEnergy"), "E (#bar{b}) [GeV]"    , 50,  0. , 500. ) );
  BBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "BBarEta"   ), "#eta (#bar{b})"       , 34, -3.4,   3.4) );
  BBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "BBarPhi"   ), "#phi (#bar{b})"       , 34, -3.4,   3.4) );
  BBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "BBarMass"  ), "M (#bar{b}) [GeV]"    , 30,  0. , 600. ) ); 
  
  LepKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepPt"    ), "p_{t} (l^{-}) [GeV]", 50,  0. , 500. ) );
  LepKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepEnergy"), "E (l^{-}) [GeV]"    , 50,  0. , 500. ) );
  LepKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepEta"   ), "#eta (l^{-})"       , 34, -3.4,   3.4) );
  LepKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepPhi"   ), "#phi (l^{-})"       , 34, -3.4,   3.4) );
  LepKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepMass"  ), "M (l^{-}) [GeV]"    , 30,  0. , 600. ) );
  
  NuBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarPt"    ), "p_{t} (#bar{#nu}^{-}) [GeV]", 50,  0. , 500. ) );
  NuBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarEnergy"), "E (#bar{#nu}^{-}) [GeV]"    , 50,  0. , 500. ) );
  NuBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarEta"   ), "#eta (#bar{#nu}^{-})"       , 34, -3.4,   3.4) );
  NuBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarPhi"   ), "#phi (#bar{#nu}^{-})"       , 34, -3.4,   3.4) );
  NuBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarMass"  ), "M (#bar{#nu}^{-}) [GeV]"    , 30,  0. , 600. ) );  
  
  // push back the same histograms for wrong charge event solutions if demanded
  if(useWrongCharge_){
    TopKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopPtWrong"    ), "p_{t} (t) [GeV]", 50,  0. , 500. ) );
    TopKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopEnergyWrong"), "E (t) [GeV]"    , 50,  0. , 500. ) );  
    TopKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopEtaWrong"   ), "#eta (t)"       , 34, -3.4,   3.4) );
    TopKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopPhiWrong"   ), "#phi (t)"       , 34, -3.4,   3.4) );
    TopKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopMassWrong"  ), "M (top) [GeV]"  , 30, 50. , 350. ) );

    WplusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WplusPtWrong"    ), "p_{t} (W^{+}) [GeV]", 50,  0. , 500. ) );
    WplusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WplusEnergyWrong"), "E (W^{+}) [GeV]"    , 50,  0. , 500. ) );
    WplusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WplusEtaWrong"   ), "#eta (W^{+})"       , 34, -3.4,   3.4) );
    WplusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WplusPhiWrong"   ), "#phi (W^{+})"       , 34, -3.4,   3.4) );
    WplusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WplusMassWrong"  ), "M (W^{+}) [GeV]"    , 30,  0. , 600. ) );

    BKin_.push_back( fs->make<TH1F>(ns.name(hist, "BPtWrong"    ), "p_{t} (b) [GeV]", 50,  0. , 500. ) );
    BKin_.push_back( fs->make<TH1F>(ns.name(hist, "BEnergyWrong"), "E (b) [GeV]"    , 50,  0. , 500. ) );
    BKin_.push_back( fs->make<TH1F>(ns.name(hist, "BEtaWrong"   ), "#eta (b)"       , 34, -3.4,   3.4) );
    BKin_.push_back( fs->make<TH1F>(ns.name(hist, "BPhiWrong"   ), "#phi (b)"       , 34, -3.4,   3.4) );
    BKin_.push_back( fs->make<TH1F>(ns.name(hist, "BMassWrong"  ), "M (b) [GeV]"    , 30,  0. , 600. ) );

    LepBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarPtWrong"    ), "p_{t} (l^{+}) [GeV]", 50,  0. , 500. ) );
    LepBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarEnergyWrong"), "E (l^{+}) [GeV]"    , 50,  0. , 500. ) );
    LepBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarEtaWrong"   ), "#eta (l^{+})"       , 34, -3.4,   3.4) );
    LepBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarPhiWrong"   ), "#phi (l^{+})"       , 34, -3.4,   3.4) );
    LepBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarMassWrong"  ), "M (l^{+}) [GeV]"    , 30,  0. , 600. ) ); 

    NuKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuPtWrong"    ), "p_{t} (#nu^{-}) [GeV]", 50,  0. , 500. ) );
    NuKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuEnergyWrong"), "E (#nu^{-}) [GeV]"    , 50,  0. , 500. ) );
    NuKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuEtaWrong"   ), "#eta (#nu^{-})"       , 34, -3.4,   3.4) );
    NuKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuPhiWrong"   ), "#phi (#nu^{-})"       , 34, -3.4,   3.4) );
    NuKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuMassWrong"  ), "M (#nu^{-}) [GeV]"    , 30,  0. , 600. ) );

    TopBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarPtWrong"    ), "p_{t} (#bar{t}) [GeV]", 50,  0. , 500. ) );
    TopBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarEnergyWrong"), "E (#bar{t}) [GeV]"    , 50,  0. , 500. ) );
    TopBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarEtaWrong"   ), "#eta (#bar{t})"       , 34, -3.4,   3.4) );
    TopBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarPhiWrong"   ), "#phi (#bar{t})"       , 34, -3.4,   3.4) );
    TopBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarMassWrong"  ), "M (#bar{t}) [GeV]"    , 30, 50. , 350. ) );

    WminusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WminusPtWrong"    ), "p_{t} (W^{-}) [GeV]", 50,  0. , 500. ) );
    WminusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WminusEnergyWrong"), "E (W^{-}) [GeV]"    , 50,  0. , 500. ) );
    WminusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WminusEtaWrong"   ), "#eta (W^{-})"       , 34, -3.4,   3.4) );
    WminusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WminusPhiWrong"   ), "#phi (W^{-})"       , 34, -3.4,   3.4) );
    WminusKin_.push_back( fs->make<TH1F>(ns.name(hist, "WminusMassWrong"  ), "M (W^{-}) [GeV]"    , 30,  0. , 600. ) );

    BBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "BBarPtWrong"    ), "p_{t} (#bar{b}) [GeV]", 50,  0. , 500. ) );
    BBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "BBarEnergyWrong"), "E (#bar{b}) [GeV]"    , 50,  0. , 500. ) );
    BBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "BBarEtaWrong"   ), "#eta (#bar{b})"       , 34, -3.4,   3.4) );
    BBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "BBarPhiWrong"   ), "#phi (#bar{b})"       , 34, -3.4,   3.4) );
    BBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "BBarMassWrong"  ), "M (#bar{b}) [GeV]"    , 30,  0. , 600. ) ); 

    LepKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepPtWrong"    ), "p_{t} (l^{-}) [GeV]", 50,  0. , 500. ) );
    LepKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepEnergyWrong"), "E (l^{-}) [GeV]"    , 50,  0. , 500. ) );
    LepKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepEtaWrong"   ), "#eta (l^{-})"       , 34, -3.4,   3.4) );
    LepKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepPhiWrong"   ), "#phi (l^{-})"       , 34, -3.4,   3.4) );
    LepKin_.push_back( fs->make<TH1F>(ns.name(hist, "LepMassWrong"  ), "M (l^{-}) [GeV]"    , 30,  0. , 600. ) );

    NuBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarPtWrong"    ), "p_{t} (#bar{#nu}^{-}) [GeV]", 50,  0. , 500. ) );
    NuBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarEnergyWrong"), "E (#bar{#nu}^{-}) [GeV]"    , 50,  0. , 500. ) );
    NuBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarEtaWrong"   ), "#eta (#bar{#nu}^{-})"       , 34, -3.4,   3.4) );
    NuBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarPhiWrong"   ), "#phi (#bar{#nu}^{-})"       , 34, -3.4,   3.4) );
    NuBarKin_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarMassWrong"  ), "M (#bar{#nu}^{-}) [GeV]"    , 30,  0. , 600. ) );    
  }
}

void
FullLepHypothesesAnalyzer::bookKinGenHistos(edm::Service<TFileService>& fs, ofstream& hist)
{

  NameScheme ns("gen");
  
  TopKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenTopPt"    ), "p_{t} (t) [GeV]", 50,  0. , 500. ) );
  TopKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenTopEnergy"), "E (t) [GeV]"    , 50,  0. , 500. ) );  
  TopKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenTopEta"   ), "#eta (t)"       , 34, -3.4,   3.4) );
  TopKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenTopPhi"   ), "#phi (t)"       , 34, -3.4,   3.4) );
  TopKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenTopMass"  ), "M (top) [GeV]"  , 30, 50. , 350. ) );

  WplusKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenWplusPt"    ), "p_{t} (W^{+}) [GeV]", 50,  0. , 500. ) );
  WplusKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenWplusEnergy"), "E (W^{+}) [GeV]"    , 50,  0. , 500. ) );
  WplusKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenWplusEta"   ), "#eta (W^{+})"       , 34, -3.4,   3.4) );
  WplusKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenWplusPhi"   ), "#phi (W^{+})"       , 34, -3.4,   3.4) );
  WplusKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenWplusMass"  ), "M (W^{+}) [GeV]"    , 30,  0. , 600. ) );

  BKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenBPt"    ), "p_{t} (b) [GeV]", 50,  0. , 500. ) );
  BKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenBEnergy"), "E (b) [GeV]"    , 50,  0. , 500. ) );
  BKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenBEta"   ), "#eta (b)"       , 34, -3.4,   3.4) );
  BKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenBPhi"   ), "#phi (b)"       , 34, -3.4,   3.4) );
  BKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenBMass"  ), "M (b) [GeV]"    , 30,  0. , 600. ) );
    
  LepBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenLepBarPt"    ), "p_{t} (l^{+}) [GeV]", 50,  0. , 500. ) );
  LepBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenLepBarEnergy"), "E (l^{+}) [GeV]"    , 50,  0. , 500. ) );
  LepBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenLepBarEta"   ), "#eta (l^{+})"       , 34, -3.4,   3.4) );
  LepBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenLepBarPhi"   ), "#phi (l^{+})"       , 34, -3.4,   3.4) );
  LepBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenLepBarMass"  ), "M (l^{+}) [GeV]"    , 30,  0. , 600. ) ); 

  NuKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenNuPt"    ), "p_{t} (#nu^{-}) [GeV]", 50,  0. , 500. ) );
  NuKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenNuEnergy"), "E (#nu^{-}) [GeV]"    , 50,  0. , 500. ) );
  NuKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenNuEta"   ), "#eta (#nu^{-})"       , 34, -3.4,   3.4) );
  NuKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenNuPhi"   ), "#phi (#nu^{-})"       , 34, -3.4,   3.4) );
  NuKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenNuMass"  ), "M (#nu^{-}) [GeV]"    , 30,  0. , 600. ) );

  TopBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenTopBarPt"    ), "p_{t} (#bar{t}) [GeV]", 50,  0. , 500. ) );
  TopBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenTopBarEnergy"), "E (#bar{t}) [GeV]"    , 50,  0. , 500. ) );
  TopBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenTopBarEta"   ), "#eta (#bar{t})"       , 34, -3.4,   3.4) );
  TopBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenTopBarPhi"   ), "#phi (#bar{t})"       , 34, -3.4,   3.4) );
  TopBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenTopBarMass"  ), "M (#bar{t}) [GeV]"    , 30, 50. , 350. ) );

  WminusKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenWminusPt"    ), "p_{t} (W^{-}) [GeV]", 50,  0. , 500. ) );
  WminusKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenWminusEnergy"), "E (W^{-}) [GeV]"    , 50,  0. , 500. ) );
  WminusKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenWminusEta"   ), "#eta (W^{-})"       , 34, -3.4,   3.4) );
  WminusKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenWminusPhi"   ), "#phi (W^{-})"       , 34, -3.4,   3.4) );
  WminusKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenWminusMass"  ), "M (W^{-}) [GeV]"    , 30,  0. , 600. ) );

  BBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenBBarPt"    ), "p_{t} (#bar{b}) [GeV]", 50,  0. , 500. ) );
  BBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenBBarEnergy"), "E (#bar{b}) [GeV]"    , 50,  0. , 500. ) );
  BBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenBBarEta"   ), "#eta (#bar{b})"       , 34, -3.4,   3.4) );
  BBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenBBarPhi"   ), "#phi (#bar{b})"       , 34, -3.4,   3.4) );
  BBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenBBarMass"  ), "M (#bar{b}) [GeV]"    , 30,  0. , 600. ) ); 
  
  LepKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenLepPt"    ), "p_{t} (l^{-}) [GeV]", 50,  0. , 500. ) );
  LepKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenLepEnergy"), "E (l^{-}) [GeV]"    , 50,  0. , 500. ) );
  LepKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenLepEta"   ), "#eta (l^{-})"       , 34, -3.4,   3.4) );
  LepKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenLepPhi"   ), "#phi (l^{-})"       , 34, -3.4,   3.4) );
  LepKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenLepMass"  ), "M (l^{-}) [GeV]"    , 30,  0. , 600. ) );
  
  NuBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenNuBarPt"    ), "p_{t} (#bar{#nu}^{-}) [GeV]", 50,  0. , 500. ) );
  NuBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenNuBarEnergy"), "E (#bar{#nu}^{-}) [GeV]"    , 50,  0. , 500. ) );
  NuBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenNuBarEta"   ), "#eta (#bar{#nu}^{-})"       , 34, -3.4,   3.4) );
  NuBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenNuBarPhi"   ), "#phi (#bar{#nu}^{-})"       , 34, -3.4,   3.4) );
  NuBarKinGen_.push_back( fs->make<TH1F>(ns.name(hist, "GenNuBarMass"  ), "M (#bar{#nu}^{-}) [GeV]"    , 30,  0. , 600. ) );  
}

void
FullLepHypothesesAnalyzer::bookKinResHistos(edm::Service<TFileService>& fs, ofstream& hist)
{

  NameScheme ns("res");

  TopKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopPtRes"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} (t)",    40, -1.0, 1.0 ) );
  TopKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopEnRes"  ), "E_{rec}-E_{gen}/E_{gen} (t)"      ,    40, -1.0, 1.0 ) );  
  TopKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopEtaRes" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} (t)", 40, -1.0, 1.0 ) );
  TopKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopPhiRes" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} (t)", 40, -1.0, 1.0 ) );
  TopKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopMassRes"), "M_{rec}-M_{gen}/M_{gen} (t)",          40, -1.0, 1.0 ) );
  
  WplusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WplusPtRes"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} (W^{+})",    40, -1.0, 1.0 ) );
  WplusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WplusEnRes"  ), "E_{rec}-E_{gen}/E_{gen} (W^{+})",          40, -1.0, 1.0 ) );  
  WplusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WplusEtaRes" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} (W^{+})", 40, -1.0, 1.0 ) );
  WplusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WplusPhiRes" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} (W^{+})", 40, -1.0, 1.0 ) );
  WplusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WplusMassRes"), "M_{rec}-M_{gen}/M_{gen} (W^{+})",          40, -1.0, 1.0 ) );  

  BKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BPtRes"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} (b)",    40, -1.0, 1.0 ) );
  BKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BEnRes"  ), "E_{rec}-E_{gen}/E_{gen} (b)",          40, -1.0, 1.0 ) );  
  BKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BEtaRes" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} (b)", 40, -1.0, 1.0 ) );
  BKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BPhiRes" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} (b)", 40, -1.0, 1.0 ) );
  BKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BMassRes"), "M_{rec}-M_{gen}/M_{gen} (b)",          40, -1.0, 1.0 ) );  

  LepBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarPtRes"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} (l^{+})",    40, -1.0, 1.0 ) );
  LepBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarEnRes"  ), "E_{rec}-E_{gen}/E_{gen} (l^{+})",          40, -1.0, 1.0 ) );  
  LepBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarEtaRes" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} (l^{+})", 40, -1.0, 1.0 ) );
  LepBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarPhiRes" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} (l^{+})", 40, -1.0, 1.0 ) );
  LepBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarMassRes"), "M_{rec}-M_{gen}/M_{gen} (l^{+})",          40, -1.0, 1.0 ) ); 

  NuKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuPtRes"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} #nu",    40, -1.0, 1.0 ) );
  NuKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuEnRes"  ), "E_{rec}-E_{gen}/E_{gen} #nu",          40, -1.0, 1.0 ) );  
  NuKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuEtaRes" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} #nu", 40, -1.0, 1.0 ) );
  NuKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuPhiRes" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} #nu", 40, -1.0, 1.0 ) );
  NuKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuMassRes"), "M_{rec}-M_{gen}/M_{gen} #nu",          40, -1.0, 1.0 ) ); 

  TopBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarPtRes"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} #bar{t}",    40, -1.0, 1.0 ) );
  TopBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarEnRes"  ), "E_{rec}-E_{gen}/E_{gen} #bar{t}",         40, -1.0, 1.0 ) );  
  TopBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarEtaRes" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} #bar{t}", 40, -1.0, 1.0 ) );
  TopBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarPhiRes" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} #bar{t}", 40, -1.0, 1.0 ) );
  TopBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarMassRes"), "M_{rec}-M_{gen}/M_{gen} #bar{t}",          40, -1.0, 1.0 ) ); 

  WminusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WminusPtRes"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} (W^{-})",    40, -1.0, 1.0 ) );
  WminusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WminusEnRes"  ), "E_{rec}-E_{gen}/E_{gen} (W^{-})",          40, -1.0, 1.0 ) );  
  WminusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WminusEtaRes" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} (W^{-})", 40, -1.0, 1.0 ) );
  WminusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WminusPhiRes" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} (W^{-})", 40, -1.0, 1.0 ) );
  WminusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WminusMassRes"), "M_{rec}-M_{gen}/M_{gen} (W^{-})",          40, -1.0, 1.0 ) ); 

  BBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BBarPtRes"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} (#bar{b})",    40, -1.0, 1.0 ) );
  BBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BBarEnRes"  ), "E_{rec}-E_{gen}/E_{gen} (#bar{b})",          40, -1.0, 1.0 ) );  
  BBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BBarEtaRes" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} (#bar{b})", 40, -1.0, 1.0 ) );
  BBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BBarPhiRes" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} (#bar{b})", 40, -1.0, 1.0 ) );
  BBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BBarMassRes"), "M_{rec}-M_{gen}/M_{gen} (#bar{b})",          40, -1.0, 1.0 ) ); 
  
  LepKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepPtRes"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} (l)",    40, -1.0, 1.0 ) );
  LepKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepEnRes"  ), "E_{rec}-E_{gen}/E_{gen} (l)",          40, -1.0, 1.0 ) );  
  LepKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepEtaRes" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} (l)", 40, -1.0, 1.0 ) );
  LepKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepPhiRes" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} (l)", 40, -1.0, 1.0 ) );
  LepKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepMassRes"), "M_{rec}-M_{gen}/M_{gen} (l)",          40, -1.0, 1.0 ) ); 
  
  NuBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarPtRes"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} (#bar{#nu})",    40, -1.0, 1.0 ) );
  NuBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarEnRes"  ), "E_{rec}-E_{gen}/E_{gen} (#bar{#nu})",          40, -1.0, 1.0 ) );  
  NuBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarEtaRes" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} (#bar{#nu})", 40, -1.0, 1.0 ) );
  NuBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarPhiRes" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} (#bar{#nu})", 40, -1.0, 1.0 ) );
  NuBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarMassRes"), "M_{rec}-M_{gen}/M_{gen} (#bar{#nu})",          40, -1.0, 1.0 ) );  
  
  // push back the same histograms for wrong charge event solutions if demanded
  if(useWrongCharge_){
    TopKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopPtResWrong"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} (t)",    40, -1.0, 1.0 ) );
    TopKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopEnResWrong"  ), "E_{rec}-E_{gen}/E_{gen} (t)"      ,    40, -1.0, 1.0 ) );  
    TopKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopEtaResWrong" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} (t)", 40, -1.0, 1.0 ) );
    TopKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopPhiResWrong" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} (t)", 40, -1.0, 1.0 ) );
    TopKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopMassResWrong"), "M_{rec}-M_{gen}/M_{gen} (t)",          40, -1.0, 1.0 ) );

    WplusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WplusPtResWrong"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} (W^{+})",    40, -1.0, 1.0 ) );
    WplusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WplusEnResWrong"  ), "E_{rec}-E_{gen}/E_{gen} (W^{+})",          40, -1.0, 1.0 ) );  
    WplusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WplusEtaResWrong" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} (W^{+})", 40, -1.0, 1.0 ) );
    WplusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WplusPhiResWrong" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} (W^{+})", 40, -1.0, 1.0 ) );
    WplusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WplusMassResWrong"), "M_{rec}-M_{gen}/M_{gen} (W^{+})",          40, -1.0, 1.0 ) );  

    BKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BPtResWrong"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} (b)",    40, -1.0, 1.0 ) );
    BKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BEnResWrong"  ), "E_{rec}-E_{gen}/E_{gen} (b)",          40, -1.0, 1.0 ) );  
    BKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BEtaResWrong" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} (b)", 40, -1.0, 1.0 ) );
    BKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BPhiResWrong" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} (b)", 40, -1.0, 1.0 ) );
    BKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BMassResWrong"), "M_{rec}-M_{gen}/M_{gen} (b)",          40, -1.0, 1.0 ) );  

    LepBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarPtResWrong"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} (l^{+})",    40, -1.0, 1.0 ) );
    LepBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarEnResWrong"  ), "E_{rec}-E_{gen}/E_{gen} (l^{+})",          40, -1.0, 1.0 ) );  
    LepBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarEtaResWrong" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} (l^{+})", 40, -1.0, 1.0 ) );
    LepBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarPhiResWrong" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} (l^{+})", 40, -1.0, 1.0 ) );
    LepBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepBarMassResWrong"), "M_{rec}-M_{gen}/M_{gen} (l^{+})",          40, -1.0, 1.0 ) ); 

    NuKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuPtResWrong"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} #nu",    40, -1.0, 1.0 ) );
    NuKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuEnResWrong"  ), "E_{rec}-E_{gen}/E_{gen} #nu",          40, -1.0, 1.0 ) );  
    NuKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuEtaResWrong" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} #nu", 40, -1.0, 1.0 ) );
    NuKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuPhiResWrong" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} #nu", 40, -1.0, 1.0 ) );
    NuKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuMassResWrong"), "M_{rec}-M_{gen}/M_{gen} #nu",          40, -1.0, 1.0 ) ); 

    TopBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarPtResWrong"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} #bar{t}",    40, -1.0, 1.0 ) );
    TopBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarEnResWrong"  ), "E_{rec}-E_{gen}/E_{gen} #bar{t}",         40, -1.0, 1.0 ) );  
    TopBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarEtaResWrong" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} #bar{t}", 40, -1.0, 1.0 ) );
    TopBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarPhiResWrong" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} #bar{t}", 40, -1.0, 1.0 ) );
    TopBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "TopBarMassResWrong"), "M_{rec}-M_{gen}/M_{gen} #bar{t}",          40, -1.0, 1.0 ) ); 

    WminusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WminusPtResWrong"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} (W^{-})",    40, -1.0, 1.0 ) );
    WminusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WminusEnResWrong"  ), "E_{rec}-E_{gen}/E_{gen} (W^{-})",          40, -1.0, 1.0 ) );  
    WminusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WminusEtaResWrong" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} (W^{-})", 40, -1.0, 1.0 ) );
    WminusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WminusPhiResWrong" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} (W^{-})", 40, -1.0, 1.0 ) );
    WminusKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "WminusMassResWrong"), "M_{rec}-M_{gen}/M_{gen} (W^{-})",          40, -1.0, 1.0 ) ); 

    BBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BBarPtResWrong"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} (#bar{b})",    40, -1.0, 1.0 ) );
    BBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BBarEnResWrong"  ), "E_{rec}-E_{gen}/E_{gen} (#bar{b})",          40, -1.0, 1.0 ) );  
    BBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BBarEtaResWrong" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} (#bar{b})", 40, -1.0, 1.0 ) );
    BBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BBarPhiResWrong" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} (#bar{b})", 40, -1.0, 1.0 ) );
    BBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "BBarMassResWrong"), "M_{rec}-M_{gen}/M_{gen} (#bar{b})",          40, -1.0, 1.0 ) ); 

    LepKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepPtResWrong"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} (l)",    40, -1.0, 1.0 ) );
    LepKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepEnResWrong"  ), "E_{rec}-E_{gen}/E_{gen} (l)",          40, -1.0, 1.0 ) );  
    LepKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepEtaResWrong" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} (l)", 40, -1.0, 1.0 ) );
    LepKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepPhiResWrong" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} (l)", 40, -1.0, 1.0 ) );
    LepKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "LepMassResWrong"), "M_{rec}-M_{gen}/M_{gen} (l)",          40, -1.0, 1.0 ) ); 

    NuBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarPtResWrong"  ), "p_{t,rec}-p_{t,gen}/p_{t,gen} (#bar{#nu})",    40, -1.0, 1.0 ) );
    NuBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarEnResWrong"  ), "E_{rec}-E_{gen}/E_{gen} (#bar{#nu})",          40, -1.0, 1.0 ) );  
    NuBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarEtaResWrong" ), "#eta_{rec}-#eta_{gen}/#eta_{gen} (#bar{#nu})", 40, -1.0, 1.0 ) );
    NuBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarPhiResWrong" ), "#phi_{rec}-#phi_{gen}/#phi_{gen} (#bar{#nu})", 40, -1.0, 1.0 ) );
    NuBarKinRes_.push_back( fs->make<TH1F>(ns.name(hist, "NuBarMassResWrong"), "M_{rec}-M_{gen}/M_{gen} (#bar{#nu})",          40, -1.0, 1.0 ) );   
  }   
}

void
FullLepHypothesesAnalyzer::bookQualityHistos(edm::Service<TFileService>& fs, ofstream& hist)
{

  NameScheme ns("qual");

  goodHypo_ = fs->make<TH1F>(ns.name(hist, "goodHypo"), "good hypothesis", 2, -0.5, 1.5);

  genMatchSumDR_ = fs->make<TH1F>(ns.name(hist, "genMatchSumDR"), "#Sigma #Delta R (genMatch)"          , 50, 0., 5.);
  genMatchSumPt_ = fs->make<TH1F>(ns.name(hist, "genMatchSumPt"), "#Sigma #Delta p_{T} (genMatch) [GeV]", 40, 0., 400.);
  
  kinSolWeight_      = fs->make<TH1F>(ns.name(hist, "kinSolWeight"), "Weight of kin solution", 50,  0. , 1. );
  kinSolWeightWrong_ = fs->make<TH1F>(ns.name(hist, "kinSolWeightWrong"), "Weight of kin solution", 50,  0. , 1. );  
  wrongCharge_       = fs->make<TH1F>(ns.name(hist, "wrongCharge"),  "is event wrong charge",   2, -0.5, 1.5);  
  
  bJetIdcs_   = fs->make<TH1F>(ns.name(hist, "bJetIdcs"   ),  "b jet indices used for hypo"     , 4, -0.5, 3.5);
  bBarJetIdcs_= fs->make<TH1F>(ns.name(hist, "bBarJetIdcs"),  "bbar jet indices used for hypo"  , 4, -0.5, 3.5);  
  elec1Idcs_  = fs->make<TH1F>(ns.name(hist, "elec1Idcs"   ),  "electron 1 indices used for hypo", 5, -1.5, 3.5);
  elec2Idcs_  = fs->make<TH1F>(ns.name(hist, "elec2Idcs"   ),  "electron 2 indices used for hypo", 5, -1.5, 3.5);    
  muon1Idcs_  = fs->make<TH1F>(ns.name(hist, "muon1Idcs"   ),  "muon 1 indices used for hypo"    , 5, -1.5, 3.5); 
  muon2Idcs_  = fs->make<TH1F>(ns.name(hist, "muon2Idcs"   ),  "muon 2 indices used for hypo"    , 5, -1.5, 3.5);  
  
  deltaM_  = fs->make<TH1F>(ns.name(hist, "deltaM"   ), "M_{top}-M{#bar{t}}", 50, -25., 25.);
  
  jetIdCorrectnes_  = fs->make<TH1F>(ns.name(hist, "jetIdCorrectnes" ), "Correctnes of Jet Id",     3, -0.5, 2.5);  
  jetFlCorrectnes_  = fs->make<TH1F>(ns.name(hist, "jetFlCorrectnes" ), "Correctnes of Jet Flavor", 3, -0.5, 2.5); 
    
  compare_  = fs->make<TH2F>(ns.name(hist, "compare"   ), "Indices", 15, -0.5, 14.5, 15, -0.5, 14.5);
  spectrum_ = fs->make<TH2F>(ns.name(hist, "spectrum"  ), "Nu spectrum", 50, 0., 250., 50, 0., 250.);       
}

void
FullLepHypothesesAnalyzer::fillKinHistos(std::vector<TH1F*>& histos, const reco::Candidate& candidate, const double& weight)
{
  histos[0]->Fill( candidate.pt()    , weight );
  histos[1]->Fill( candidate.energy(), weight );  
  histos[2]->Fill( candidate.eta()   , weight );
  histos[3]->Fill( candidate.phi()   , weight );
  histos[4]->Fill( candidate.mass()  , weight );
}

void
FullLepHypothesesAnalyzer::fillKinGenHistos(std::vector<TH1F*>& histos, const reco::Candidate& candidate, const double& weight)
{
  histos[0]->Fill( candidate.pt()    , weight );
  histos[1]->Fill( candidate.energy(), weight );  
  histos[2]->Fill( candidate.eta()   , weight );
  histos[3]->Fill( candidate.phi()   , weight );
  histos[4]->Fill( candidate.mass()  , weight );
}

void
FullLepHypothesesAnalyzer::fillWrongChargeKinHistos(std::vector<TH1F*>& histos, const reco::Candidate& candidate, const double& weight)
{
  histos[5]->Fill( candidate.pt()    , weight );
  histos[6]->Fill( candidate.energy(), weight );  
  histos[7]->Fill( candidate.eta()   , weight );
  histos[8]->Fill( candidate.phi()   , weight );
  histos[9]->Fill( candidate.mass()  , weight );
}

void
FullLepHypothesesAnalyzer::fillKinResHistos(std::vector<TH1F*>& histos, const reco::Candidate& candidate,
					    const reco::Candidate& genCandidate, const double& weight)
{
  histos[0]->Fill( (candidate.pt()    -genCandidate.pt()    )/genCandidate.pt()    , weight );
  histos[1]->Fill( (candidate.energy()-genCandidate.energy())/genCandidate.energy(), weight );  
  histos[2]->Fill( (candidate.eta()   -genCandidate.eta()   )/genCandidate.eta()   , weight );
  histos[3]->Fill( (candidate.phi()   -genCandidate.phi()   )/genCandidate.phi()   , weight );
  histos[4]->Fill( (candidate.mass()  -genCandidate.mass()  )/genCandidate.mass()  , weight );
}

void
FullLepHypothesesAnalyzer::fillWrongChargeKinResHistos(std::vector<TH1F*>& histos, const reco::Candidate& candidate,
					    const reco::Candidate& genCandidate, const double& weight)
{
  histos[5]->Fill( (candidate.pt()    -genCandidate.pt()    )/genCandidate.pt()    , weight );
  histos[6]->Fill( (candidate.energy()-genCandidate.energy())/genCandidate.energy(), weight );  
  histos[7]->Fill( (candidate.eta()   -genCandidate.eta()   )/genCandidate.eta()   , weight );
  histos[8]->Fill( (candidate.phi()   -genCandidate.phi()   )/genCandidate.phi()   , weight );
  histos[9]->Fill( (candidate.mass()  -genCandidate.mass()  )/genCandidate.mass()  , weight );
}

void
FullLepHypothesesAnalyzer::fillQualityHistos(const TtFullLeptonicEvent& FullLepEvt, 
                                             const TtEvent::HypoClassKey& hypoKey, 
					     const double& weight)
{
  bJetIdcs_   ->Fill(FullLepEvt.jetLeptonCombination(hypoKey)[0] );
  bBarJetIdcs_->Fill(FullLepEvt.jetLeptonCombination(hypoKey)[1] );  
  elec1Idcs_  ->Fill(FullLepEvt.jetLeptonCombination(hypoKey)[2] );  
  elec2Idcs_  ->Fill(FullLepEvt.jetLeptonCombination(hypoKey)[3] );  
  muon1Idcs_  ->Fill(FullLepEvt.jetLeptonCombination(hypoKey)[4] );  
  muon2Idcs_  ->Fill(FullLepEvt.jetLeptonCombination(hypoKey)[5] ); 
  
  spectrum_->Fill(FullLepEvt.neutrino(hypoKey)->energy(), FullLepEvt.neutrinoBar(hypoKey)->energy());
    
  // genMatch histos
  if( hypoKey==3 ) {
    genMatchSumDR_->Fill( FullLepEvt.genMatchSumDR() , weight );
    genMatchSumPt_->Fill( FullLepEvt.genMatchSumPt() , weight ); 
  } 
  
  if( hypoKey==6 ) {
    deltaM_->Fill(FullLepEvt.top(hypoKey)->mass()-FullLepEvt.topBar(hypoKey)->mass()); 
    if(!FullLepEvt.isWrongCharge()) kinSolWeight_->Fill( FullLepEvt.solWeight() , weight );
    else kinSolWeightWrong_->Fill( FullLepEvt.solWeight() , weight );
    wrongCharge_ ->Fill( FullLepEvt.isWrongCharge() , weight ); 
  } 
  
  // indices Hist (for muon muon only) 
  const TtEvent::HypoClassKey genKey = TtEvent::HypoClassKey(3);
  const TtEvent::HypoClassKey kinKey = TtEvent::HypoClassKey(6);
   
  if(FullLepEvt.isHypoValid(genKey) && FullLepEvt.isHypoValid(kinKey)){
    compare_->Fill(FullLepEvt.jetLeptonCombination(genKey)[0]  , FullLepEvt.jetLeptonCombination(kinKey)[0]  );
    compare_->Fill(FullLepEvt.jetLeptonCombination(genKey)[1]+4, FullLepEvt.jetLeptonCombination(kinKey)[1]+4);   
    compare_->Fill(FullLepEvt.jetLeptonCombination(genKey)[4]+8, FullLepEvt.jetLeptonCombination(kinKey)[4]+8);
    compare_->Fill(FullLepEvt.jetLeptonCombination(genKey)[5]+12, FullLepEvt.jetLeptonCombination(kinKey)[5]+12);  
    
    // are both b-jets identifierd right or only one or none of them?
    if(FullLepEvt.jetLeptonCombination(kinKey)[0] == FullLepEvt.jetLeptonCombination(genKey)[0] &&
     FullLepEvt.jetLeptonCombination(kinKey)[1] == FullLepEvt.jetLeptonCombination(genKey)[1]){
      jetIdCorrectnes_->SetBinContent(1,jetIdCorrectnes_->GetBinContent(1)+weight);
    }
    else if(FullLepEvt.jetLeptonCombination(kinKey)[0] == FullLepEvt.jetLeptonCombination(genKey)[0] ||
     FullLepEvt.jetLeptonCombination(kinKey)[1] == FullLepEvt.jetLeptonCombination(genKey)[1]){
      jetIdCorrectnes_->SetBinContent(2,jetIdCorrectnes_->GetBinContent(2)+weight);
    }
    else{
      jetIdCorrectnes_->SetBinContent(3,jetIdCorrectnes_->GetBinContent(3)+weight);
    } 
    
    //how many of the two jets have the right flavor?
    if((FullLepEvt.jetLeptonCombination(kinKey)[0] == FullLepEvt.jetLeptonCombination(genKey)[0] ||
     FullLepEvt.jetLeptonCombination(kinKey)[0] == FullLepEvt.jetLeptonCombination(genKey)[1]) &&
     (FullLepEvt.jetLeptonCombination(kinKey)[1] == FullLepEvt.jetLeptonCombination(genKey)[1] ||
      FullLepEvt.jetLeptonCombination(kinKey)[1] == FullLepEvt.jetLeptonCombination(genKey)[0])){
      jetFlCorrectnes_->SetBinContent(1,jetFlCorrectnes_->GetBinContent(1)+weight);
    }
    else if((FullLepEvt.jetLeptonCombination(kinKey)[0] == FullLepEvt.jetLeptonCombination(genKey)[0] ||
     FullLepEvt.jetLeptonCombination(kinKey)[0] == FullLepEvt.jetLeptonCombination(genKey)[1]) ||
     (FullLepEvt.jetLeptonCombination(kinKey)[1] == FullLepEvt.jetLeptonCombination(genKey)[1] ||
      FullLepEvt.jetLeptonCombination(kinKey)[1] == FullLepEvt.jetLeptonCombination(genKey)[0])){
      jetFlCorrectnes_->SetBinContent(2,jetFlCorrectnes_->GetBinContent(2)+weight);
    }
    else{
      jetFlCorrectnes_->SetBinContent(3,jetFlCorrectnes_->GetBinContent(3)+weight);
    }             
  }   
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( FullLepHypothesesAnalyzer );
