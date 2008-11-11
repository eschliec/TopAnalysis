#include "PhysicsTools/MVATrainer/interface/HelperMacros.h"

#include "TopAnalysis/TopFilter/plugins/TtSemiLepSignalSelectorMVATrainer.h"
#include "TopAnalysis/TopUtils/interface/TtSemiLepSignalSelectorEval.h"
#include "TopQuarkAnalysis/TopTools/interface/TtSemiEvtPartons.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/Math/interface/deltaR.h"
typedef std::vector<pat::Jet> TopJetCollection;

TtSemiLepSignalSelectorMVATrainer::TtSemiLepSignalSelectorMVATrainer(const edm::ParameterSet& cfg) :
	leptons_(cfg.getParameter<edm::InputTag> ("leptons")), jets_(cfg.getParameter<edm::InputTag> ("jets")), METs_(
			cfg.getParameter<edm::InputTag> ("METs")), nJetsMax_(cfg.getParameter<int> ("nJetsMax")), lepChannel_(
			cfg.getParameter<int> ("lepChannel")), minMuonPt_(cfg.getParameter<double> ("minMuonPt")), caloIso_(
			cfg.getParameter<double> ("caloIso")), trackIso_(cfg.getParameter<double> ("trackIso")), jetIso_(
			cfg.getParameter<double> ("jetIso")), maxEtaMuon_(cfg.getParameter<double> ("maxEtaMuon")), maxEtaJets_(
			cfg.getParameter<double> ("maxEtaJets")), minJetEt_(cfg.getParameter<double> ("minJetEt")),
			minLeadingJetEt_(cfg.getParameter<double> ("minLeadingJetEt")) {
}

TtSemiLepSignalSelectorMVATrainer::~TtSemiLepSignalSelectorMVATrainer() {
}

void TtSemiLepSignalSelectorMVATrainer::analyze(const edm::Event& evt, const edm::EventSetup& setup) {
	mvaComputer.update<TtSemiLepSignalSelectorMVARcd> ("trainer", setup, "ttSemiLepSignalSelectorMVA");

	// can occur in the last iteration when the
	// MVATrainer is about to save the result
	if (!mvaComputer)
		return;

	const pat::Muon *lepton;
	const pat::MET* MET;
	double weight = 0.;

	edm::Handle<TtGenEvent> genEvt;
	evt.getByLabel("genEvt", genEvt);

	edm::Handle<double> weightHandle;
	evt.getByLabel("eventWeight", weightHandle);
	weight = *weightHandle;

	edm::Handle<std::vector<pat::MET> > metH;
	evt.getByLabel(METs_, metH);
	MET = &(*metH)[0];

	edm::Handle<std::vector<pat::Muon> > leptons;
	evt.getByLabel(leptons_, leptons);

	//lepton cuts
	if (leptons->size() != 1)
		return;
	if (leptons->begin()->caloIso() >= caloIso_)
		return;
	if (leptons->begin()->trackIso() >= trackIso_)
		return;
	if (leptons->begin()->pt() <= minMuonPt_)
		return;
	if (fabs(leptons->begin()->eta()) >= maxEtaMuon_)
		return;

	lepton = &(*leptons)[0];

	edm::Handle< std::vector<pat::Jet> > jet_handle;
	evt.getByLabel(jets_, jet_handle);
	if (!jet_handle.isValid())
		return;
	const TopJetCollection jets = *jet_handle;
	if (jets.begin()->et() <= minLeadingJetEt_)
		return;

	double dRmin = 9999.;
	TopJetCollection seljets;
	for (std::vector<pat::Jet>::const_iterator it = jets.begin(); it != jets.end(); it++) {
		if (it->et() > minJetEt_ && fabs(it->eta()) < maxEtaJets_) {
			double tmpdR = deltaR(it->eta(), it->phi(), lepton->eta(), lepton->phi());
			if (tmpdR < dRmin)
				dRmin = tmpdR;
			seljets.push_back(*it);
		}
	}
	//jets more far from lepton than dRmin
	if (dRmin >= jetIso_)
		return;

	// skip events with less than 4 jets
	if (seljets.size() < 4)
		return;

	if (!(leptons->size() == 0)) {

		TtSemiLepSignalSelector selection(seljets, lepton, MET);

		if (genEvt->isSemiLeptonic() && genEvt->semiLeptonicChannel() == lepChannel_) {
			evaluateTtSemiLepSignalSelector(mvaComputer, selection, weight, true, true);
		} else {
			evaluateTtSemiLepSignalSelector(mvaComputer, selection, weight, true, false);
		}

	} else
		return;
}

// implement the plugins for the trainer
// -> defines TtSemiLepSignalSelMVAContainerSaveCondDB
// -> defines TtSemiLepSignalSelMVASaveFile
// -> defines TtSemiLepSignalSelMVATrainerLooper
MVA_TRAINER_IMPLEMENT(TtSemiLepSignalSelectorMVA);
