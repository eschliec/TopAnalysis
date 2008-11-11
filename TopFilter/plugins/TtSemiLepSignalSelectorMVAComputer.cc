#include "PhysicsTools/JetMCUtils/interface/combination.h"

#include "TopAnalysis/TopFilter/plugins/TtSemiLepSignalSelectorMVAComputer.h"
#include "TopQuarkAnalysis/TopTools/interface/TtSemiEvtPartons.h"
#include "TopAnalysis/TopUtils/interface/TtSemiLepSignalSelectorEval.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"

TtSemiLepSignalSelectorMVAComputer::TtSemiLepSignalSelectorMVAComputer(const edm::ParameterSet& cfg) :
	leptons_(cfg.getParameter<edm::InputTag> ("leptons")), jets_(cfg.getParameter<edm::InputTag> ("jets")), METs_(
			cfg.getParameter<edm::InputTag> ("METs")), nJetsMax_(cfg.getParameter<int> ("nJetsMax")), lepChannel_(
			cfg.getParameter<int> ("lepChannel")), minMuonPt_(cfg.getParameter<double> ("minMuonPt")), caloIso_(
			cfg.getParameter<double> ("caloIso")), trackIso_(cfg.getParameter<double> ("trackIso")), jetIso_(
			cfg.getParameter<double> ("jetIso")), maxEtaMuon_(cfg.getParameter<double> ("maxEtaMuon")), maxEtaJets_(
			cfg.getParameter<double> ("maxEtaJets")), minJetEt_(cfg.getParameter<double> ("minJetEt")),
			minLeadingJetEt_(cfg.getParameter<double> ("minLeadingJetEt")) {
	produces<double> ("DiscSel");
}

TtSemiLepSignalSelectorMVAComputer::~TtSemiLepSignalSelectorMVAComputer() {
}

void TtSemiLepSignalSelectorMVAComputer::produce(edm::Event& evt, const edm::EventSetup& setup) {

	std::auto_ptr< double> pOutDisc(new double);

	mvaComputer.update<TtSemiLepSignalSelectorMVARcd> (setup, "ttSemiLepSignalSelectorMVA");

	// read name of the last processor in the MVA calibration
	// (to be used as meta information)
	edm::ESHandle<PhysicsTools::Calibration::MVAComputerContainer> calibContainer;
	setup.get<TtSemiLepSignalSelectorMVARcd> ().get(calibContainer);
	std::vector<PhysicsTools::Calibration::VarProcessor*> processors = (calibContainer->find(
			"ttSemiLepSignalSelectorMVA")).getProcessors();

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


	double discrim;

	// skip events with no appropriate lepton candidate in
	if (leptons->size() != 1 && seljets.size() < 4 && leptons->begin()->caloIso() >= caloIso_
			&& leptons->begin()->trackIso() >= trackIso_ && jets.begin()->et() <= minLeadingJetEt_ && dRmin >= jetIso_
			|| leptons->begin()->pt() <= minMuonPt_ || fabs(leptons->begin()->eta()) >= maxEtaMuon_)
		discrim = -1.;
	else {
		TtSemiLepSignalSelector selection(seljets, lepton, MET);

		discrim = evaluateTtSemiLepSignalSelector(mvaComputer, selection, weight);
	}

	*pOutDisc = discrim;

	evt.put(pOutDisc, "DiscSel");
}

void TtSemiLepSignalSelectorMVAComputer::beginJob(const edm::EventSetup&) {
}

void TtSemiLepSignalSelectorMVAComputer::endJob() {
}

// implement the plugins for the computer container
// -> register TtSemiLepSignalSelMVARcd
// -> define TtSemiLepSignalSelMVAFileSource
MVA_COMPUTER_CONTAINER_IMPLEMENT(TtSemiLepSignalSelectorMVA);
