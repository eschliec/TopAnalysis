#include "TopAnalysis/TopAnalyzer/plugins/MatrixAnalyzer.h"

using std::cout;
using std::endl;
using reco::GenParticle;

MatrixAnalyzer::MatrixAnalyzer(const edm::ParameterSet& cfg) :
	hist_(cfg.getParameter<std::string> ("hist")),
	pmodulename_(cfg.getParameter<std::string> ("pmn")),
	before_(cfg.getParameter<bool> ("before")),
	muons_(cfg.getParameter<edm::InputTag> ("muons")),
	var_(cfg.getParameter<edm::InputTag> ("var")),
	varBins_(cfg.getParameter<std::vector<double> > ("varBins")) {
	debug_ = true;
	noBins_ = 5;
	beforeBin_ = 2;
	afterBin_ = 4;
	sampleweight_ = 0.;


	Counters_ = new LeptonCounter();
	countersBefore_ = new LeptonCounter();
	effCounter_ = new LeptonCounter();
	effErrors_ =  new LeptonCounter();
	for (unsigned int x = 0; x < varBins_.size() - 1; x++) {
		std::stringstream tmp;
		tmp << varBins_.at(x);
		Counters_->addCounter(tmp.str());
		countersBefore_->addCounter(tmp.str());
		effCounter_->addCounter(tmp.str());
		effErrors_->addCounter(tmp.str());
	}
	Counters_->addCounter("weighted");
	countersBefore_->addCounter("weighted");
	effCounter_->addCounter("weighted");
	effErrors_->addCounter("weighted");

	Counters_->addCounter("simple");
	countersBefore_->addCounter("simple");
	effCounter_->addCounter("simple");
	effErrors_->addCounter("simple");
}

MatrixAnalyzer::~MatrixAnalyzer() {
}

void MatrixAnalyzer::beginJob(const edm::EventSetup&) {

	if (hist_.empty())
		return;

	//TODO: Beim 2ten durchlauf sollen die histogramme geladen werden.
	edm::Service<TFileService> fs;
	if (!fs) {
		throw edm::Exception(edm::errors::Configuration,
				"TFile Service is not registered in cfg file");
	}
	f_ = &fs->file();
	ofstream hist(hist_.c_str(), std::ios::out);

	NameScheme nam("mbg");

	//make histograms
	background_ = fs->make<TH1F> (nam.name(hist, "N_background"), nam.name(
			"background"), 300, 0, 1.);
	lep_ = fs->make<TH1F> (nam.name(hist, "N_semilepton"), nam.name("lepton"),
			300, 0, 1.);
	llep_ = fs->make<TH1F> (nam.name(hist, "N_dilepton"), nam.name("ll"), 300,
			0, 1.);
	multilep_ = fs->make<TH1F> (nam.name(hist, "N_multilepton"),
			nam.name("ml"), 300, 0, 1.);
	overall_ = fs->make<TH1F> (nam.name(hist, "N_overall"),
			nam.name("overall"), noBins_, 0, noBins_);

	eff_ = fs->make<TH1F> (nam.name(hist, "efficency"), nam.name("eff"), 300,
			0, 1.);

	binnedBkg_ = fs->make<TH1F> (nam.name(hist, "binnedBkg"), nam.name(
			"binnedBkg"), varBins_.size() * 2 + 2, 0, 1.);
	binnedSemiLep_ = fs->make<TH1F> (nam.name(hist, "binnedSemiLep"), nam.name(
			"binnedSemiLep"), varBins_.size() * 2 + 2, 0, 1.);
	binnedDiLep_ = fs->make<TH1F> (nam.name(hist, "binnedDiLep"), nam.name(
			"binnedDiLep"), varBins_.size() * 2 + 2, 0, 1.);
	binnedMultiLep_ = fs->make<TH1F> (nam.name(hist, "binnedMultiLep"),
			nam.name("binnedMultiLep"), varBins_.size() * 2 + 2, 0, 1.);
	binnedEff_ = fs->make<TH1F> (nam.name(hist, "binnedEff"), nam.name(
			"binnedEff"), varBins_.size() + 2, 0, 1.);
	binnedOverall_ = fs->make<TH1F> (nam.name(hist, "binnedOverall"), nam.name(
			"binnedOverall"), varBins_.size()*2 + 2, 0, 1.);
	varPlot_ = fs->make<TH1F> (nam.name(hist, var_.label().c_str()), nam.name(var_.label().c_str()),
			500, 0, 500);
	//set number of bins (noBins_)
	background_->SetBins(noBins_, 0., 5.);
	lep_->SetBins(noBins_, 0., 5.);
	llep_->SetBins(noBins_, 0., 5.);
	multilep_->SetBins(noBins_, 0., 5.);
	eff_->SetBins(9, 0., 5.);
}

void MatrixAnalyzer::analyze(const edm::Event& evt,
		const edm::EventSetup& setup) {
	edm::Handle<reco::GenParticleCollection> genParticles;
	evt.getByLabel("genParticles", genParticles);

	edm::Handle<TopMuonCollection> muons;
	evt.getByLabel(muons_, muons);

	edm::Handle<edm::View<pat::Muon> > recMuons;
	evt.getByLabel(muons_, recMuons);

	edm::Handle<edm::View<reco::RecoCandidate> > recVars;
	evt.getByLabel(var_, recVars);

	edm::Handle<edm::Association<reco::GenParticleCollection> > genMatch;
	evt.getByLabel("muonMatch", genMatch);

	edm::Handle<double> weightHandle;
	evt.getByLabel("eventWeight", weightHandle);

	sampleweight_ = *weightHandle;

	Double_t pt = (*recVars)[0].pt();
	varPlot_->Fill(pt, sampleweight_);
	for (unsigned int i = 0; i < varBins_.size() - 1; i++) {
		if (pt >= varBins_[i] && pt < varBins_[i + 1]) {
			double t = varBins_.at(i);
			std::stringstream tmp;
			tmp << t;
			if (muons->size() < 1) {
				Counters_->addPureHadronic(tmp.str(), sampleweight_);
			}
			else if (muons->size() == 1) {
				Counters_->addSemiLeptonic(tmp.str(), sampleweight_);
			}
			else if (muons->size() == 2) {
				Counters_->addDiLeptonic(tmp.str(),sampleweight_);
			}
			else if (muons->size() > 2) {
				Counters_->addMultiLeptonic(tmp.str(), sampleweight_);
			}
		}
	}
	if (muons->size() < 1) {
		Counters_->addPureHadronic("simple");
		Counters_->addPureHadronic("weighted", sampleweight_);
	}

	if (muons->size() == 1) {
		Counters_->addSemiLeptonic("simple");
		Counters_->addSemiLeptonic("weighted", sampleweight_);
	}
	if (muons->size() == 2) {
		Counters_->addDiLeptonic("simple");
		Counters_->addDiLeptonic("weighted", sampleweight_);
	}
	if (muons->size() > 2) {
		Counters_->addMultiLeptonic("simple");
		Counters_->addMultiLeptonic("weighted", sampleweight_);
	}
}

void MatrixAnalyzer::endJob() {
	if (before_) {
		cout << "#----- before (" << "Module: " << pmodulename_ <<")-----#" << endl;
		setEnv(beforeBin_);
		cout << endl;
	} else {
		getBefore();
		cout << "#----- after (" << "Module: " << pmodulename_ <<")-----#" << endl;
		log(sampleweight_, "::endJob >> sampleweight", false);
		setEnv(afterBin_);
		cout.precision(10);
		if (countersBefore_->getPureHadronic("weighted") != 0) {
			double eff, err;
			eff = Counters_->getPureHadronic("weighted") / countersBefore_->getPureHadronic("weighted");
			err = sqrt(eff * (1 - eff) / countersBefore_->getPureHadronic("weighted"));
			cout << "hadEff " << eff << " +- " << err << endl;
			effCounter_->setPureHadronic("weighted", eff);
			effErrors_->setPureHadronic("weighted", err);
		}

		if (countersBefore_->getSemiLeptonic("weighted") != 0) {
			double eff, err;
			eff = Counters_->getSemiLeptonic("weighted") / countersBefore_->getSemiLeptonic("weighted");
			err = sqrt(eff * (1 - eff) / countersBefore_->getSemiLeptonic(
					"weighted"));
			cout << "semiEff " << eff << " +- " << err << endl;
			effCounter_->setSemiLeptonic("weighted", eff);
			effErrors_->setSemiLeptonic("weighted", err);
		}
		if (countersBefore_->getDiLeptonic("weighted") != 0) {
			double eff, err;
			eff = Counters_->getDiLeptonic("weighted") / countersBefore_->getDiLeptonic("weighted");
			err = sqrt(eff * (1 - eff) / countersBefore_->getDiLeptonic(
					"weighted"));
			cout << "diEff " << eff << " +- " << err << endl;
			effCounter_->setDiLeptonic("weighted", eff);
			effErrors_->setDiLeptonic("weighted", err);
		}
		if (countersBefore_->getMultiLeptonic("weighted") != 0) {
			double eff, err;
			eff = Counters_->getMultiLeptonic("weighted") / countersBefore_->getMultiLeptonic("weighted");
			err = sqrt(eff * (1 - eff) / countersBefore_->getMultiLeptonic(
					"weighted"));
			cout << "multiEff " << eff << " +- " << err << endl;
			effCounter_->setMultiLeptonic("weighted", eff);
			effErrors_->setMultiLeptonic("weighted", err);
		}

		for (unsigned int i = 0; i < varBins_.size() - 1; i++) {
			std::stringstream tmp;
			tmp << varBins_.at(i);
			double b, s, d, m, bb, sb, db, mb, o, ob;
			b = Counters_->getPureHadronic(tmp.str());
			s = Counters_->getSemiLeptonic(tmp.str());
			d = Counters_->getDiLeptonic(tmp.str());
			m = Counters_->getMultiLeptonic(tmp.str());
			o = Counters_->getAllLeptonic(tmp.str());

			bb = countersBefore_->getPureHadronic(tmp.str());
			sb = countersBefore_->getSemiLeptonic(tmp.str());
			db = countersBefore_->getDiLeptonic(tmp.str());
			mb = countersBefore_->getMultiLeptonic(tmp.str());
			ob = countersBefore_->getAllLeptonic(tmp.str());
			if (bb != 0.) {
				double eff,err;
				eff = b / bb;
				err = sqrt(eff * (1 - eff) / bb);
				effCounter_->setPureHadronic(tmp.str(), eff);
				effErrors_->setPureHadronic(tmp.str(), err);
				cout << tmp.str() << " (had) " << eff << " +- " << err << endl;
			}
			if (sb != 0.) {
				double eff,err;
				eff = s/sb;
				err = sqrt(eff * (1 - eff) / sb);
				effCounter_->setSemiLeptonic(tmp.str(), eff);
				effErrors_->setSemiLeptonic(tmp.str(), err);
				cout << tmp.str() << " (semi) " << eff << " +- " << err << endl;
			}
			if (db != 0.) {
				double eff, err;
				eff = d / db;
				err = sqrt(eff * (1 - eff) / db);
				effCounter_->setDiLeptonic(tmp.str(),eff);
				effErrors_->setDiLeptonic(tmp.str(), err);
				cout << tmp.str() << " (di) " << eff << " +- " << err << endl;
			}
			if (mb != 0.) {
				double eff, err;
				eff = m / mb;
				err = sqrt(eff * (1 - eff) / mb);
				effCounter_->setMultiLeptonic(tmp.str(), eff);
				effErrors_->setMultiLeptonic(tmp.str(), err);
				cout << tmp.str() << " (multi) " << eff << " +- " << err << endl;
			}

			if (ob != 0.) {
				double eff, err;
				eff = o / ob;
				err = sqrt(eff * (1 - eff) / ob);
//				effCounter_->setSemiLeptonic(tmp.str(), eff);
//				effErrors_->setSemiLeptonic(tmp.str(), err);
				cout << tmp.str() << " (overall) " << eff << " +- " << err << endl;
			}

			binnedEff_->SetBinContent(i+1, effCounter_->getSemiLeptonic(tmp.str()));
			binnedEff_->SetBinError(i+1, effErrors_->getSemiLeptonic(tmp.str()));
			binnedEff_->GetXaxis()->SetBinLabel(i + 1, tmp.str().c_str());
		}

		eff_->SetBinContent(2, effCounter_->getPureHadronic("weighted"));
		eff_->SetBinError(2, effErrors_->getPureHadronic("weighted"));

		eff_->SetBinContent(4, effCounter_->getSemiLeptonic("weighted"));
		eff_->SetBinError(4, effErrors_->getSemiLeptonic("weighted"));

		eff_->SetBinContent(6, effCounter_->getDiLeptonic("weighted"));
		eff_->SetBinError(6, effErrors_->getDiLeptonic("weighted"));

		eff_->SetBinContent(8, effCounter_->getMultiLeptonic("weighted"));
		eff_->SetBinError(8, effErrors_->getMultiLeptonic("weighted"));

		background_->SetBinContent(beforeBin_, countersBefore_->getPureHadronic("weighted"));
		lep_->SetBinContent(beforeBin_, countersBefore_->getSemiLeptonic("weighted"));
		llep_->SetBinContent(beforeBin_, countersBefore_->getDiLeptonic("weighted"));
		multilep_->SetBinContent(beforeBin_, countersBefore_->getMultiLeptonic("weighted"));
		overall_->SetBinContent(beforeBin_, countersBefore_->getAllLeptonic("weighted"));


	}
	//TODO: Fill in histogram.
//	for (map<int, TopMuonCollection>::iterator iter = mothermap_.begin(); iter
//			!= mothermap_.end(); iter++) {
//		cout << "total matched mu from " << (*iter).first << " : "
//				<< (*iter).second.size() << endl;
//	}

}

void MatrixAnalyzer::setEnv(int &bin) {
	if(!before_){
		cout << "type >> events with muons (after) <> weighted events with muons (before) <> gen. events" << endl;
		//cout << "hadr >> " << Counters_->getPureHadronic("weighted") << " <> " << countersBefore_->getPureHadronic("weighted") << " <> " << Counters_->getPureHadronic("simple") << ")" << endl;
		cout << "semi >> " << Counters_->getSemiLeptonic("weighted") << " <> " << countersBefore_->getSemiLeptonic("weighted") << " <> " << Counters_->getSemiLeptonic("simple") << ")" << endl;
		//cout << "llep >> " << Counters_->getDiLeptonic("weighted") << " <> " << countersBefore_->getPureHadronic("weighted") << " <> " << Counters_->getPureHadronic("simple") << ")" << endl;
		//cout << "mlep >> " << Counters_->getMultiLeptonic("weighted") << " <> " << countersBefore_->getPureHadronic("weighted") << " <> " << Counters_->getPureHadronic("simple") << ")" << endl;
		cout << "overall >> " << Counters_->getAllLeptonic("weighted") << " <> " << countersBefore_->getAllLeptonic("weighted") << " <> " << Counters_->getSemiLeptonic("simple") << ")" << endl;
	}


	for (unsigned int i = 0; i < varBins_.size() - 1; i++) {
		std::stringstream tmp;
		tmp << varBins_.at(i);
		double b, s, d, m, o, bb, sb, db, mb, ob;
		b = Counters_->getPureHadronic(tmp.str());
		s = Counters_->getSemiLeptonic(tmp.str());
		d = Counters_->getDiLeptonic(tmp.str());
		m = Counters_->getMultiLeptonic(tmp.str());
		o = Counters_->getAllLeptonic(tmp.str());

		bb = countersBefore_->getPureHadronic(tmp.str());
		sb = countersBefore_->getSemiLeptonic(tmp.str());
		db = countersBefore_->getDiLeptonic(tmp.str());
		mb = countersBefore_->getMultiLeptonic(tmp.str());
		ob = countersBefore_->getAllLeptonic(tmp.str());

		cout << varBins_[i] << " (semi):" << s << " <> " << sb << endl;
		cout << varBins_[i] << " (overall):" << o << " <> " << ob <<  endl;
		if(before_){
			binnedBkg_->SetBinContent(2 * i + 1, b);
			binnedSemiLep_->SetBinContent(2 * i + 1, s);
			binnedDiLep_->SetBinContent(2 * i + 1, d);
			binnedMultiLep_->SetBinContent(2 * i + 1, m);
			binnedOverall_->SetBinContent(2 * i + 1, o);
			binnedOverall_->GetXaxis()->SetBinLabel(2 * i + 1, "B");
		}
		else {
			binnedBkg_->SetBinContent(2 * i + 2, b);
			binnedSemiLep_->SetBinContent(2 * i + 2, s);
			binnedDiLep_->SetBinContent(2 * i + 2, d);
			binnedMultiLep_->SetBinContent(2 * i + 2, m);
			binnedOverall_->SetBinContent(2 * i + 2, o);

			binnedBkg_->SetBinContent(2 * i + 1, bb);
			binnedSemiLep_->SetBinContent(2 * i + 1, sb);
			binnedDiLep_->SetBinContent(2 * i + 1, db);
			binnedMultiLep_->SetBinContent(2 * i + 1, mb);
			binnedOverall_->SetBinContent(2 * i + 1, ob);

		}
		binnedBkg_->GetXaxis()->SetBinLabel(2 * i + 1, "B");
		binnedBkg_->GetXaxis()->SetBinLabel(2 * i + 2, "A");

		binnedSemiLep_->GetXaxis()->SetBinLabel(2 * i + 1, "B");
		binnedSemiLep_->GetXaxis()->SetBinLabel(2 * i + 2, "A");

		binnedDiLep_->GetXaxis()->SetBinLabel(2 * i + 1, "B");
		binnedDiLep_->GetXaxis()->SetBinLabel(2 * i + 2, "A");

		binnedMultiLep_->GetXaxis()->SetBinLabel(2 * i + 1, "B");
		binnedMultiLep_->GetXaxis()->SetBinLabel(2 * i + 2, "A");

		binnedOverall_->GetXaxis()->SetBinLabel(2 * i + 1, "B");
		binnedOverall_->GetXaxis()->SetBinLabel(2 * i + 2, "A");

	}

	background_->SetBinContent(bin, Counters_->getPureHadronic("weighted"));
	lep_->SetBinContent(bin, Counters_->getSemiLeptonic("weighted"));
	llep_->SetBinContent(bin, Counters_->getDiLeptonic("weighted"));
	multilep_->SetBinContent(bin, Counters_->getMultiLeptonic("weighted"));
	overall_->SetBinContent(bin, Counters_->getAllLeptonic("weighted"));
	if (before_) {
		background_->GetXaxis()->SetBinLabel(bin, "before");
		lep_->GetXaxis()->SetBinLabel(bin, "before");
		llep_->GetXaxis()->SetBinLabel(bin, "before");
		multilep_->GetXaxis()->SetBinLabel(bin, "before");
		overall_->GetXaxis()->SetBinLabel(bin, "before");
	}
	else {
		background_->GetXaxis()->SetBinLabel(bin, "after");
		lep_->GetXaxis()->SetBinLabel(bin, "after");
		llep_->GetXaxis()->SetBinLabel(bin, "after");
		multilep_->GetXaxis()->SetBinLabel(bin, "after");
		overall_->GetXaxis()->SetBinLabel(bin, "after");
	}

}

void MatrixAnalyzer::getBefore() {
	countersBefore_->setPureHadronic("weighted", getHist(pmodulename_, "mbg_N_background", beforeBin_));
	countersBefore_->setSemiLeptonic("weighted", getHist(pmodulename_, "mbg_N_semilepton", beforeBin_));
	countersBefore_->setDiLeptonic("weighted", getHist(pmodulename_, "mbg_N_dilepton", beforeBin_));
	countersBefore_->setMultiLeptonic("weighted",getHist(pmodulename_, "mbg_N_multilepton", beforeBin_));

	cout << "getBefore s" << countersBefore_->getSemiLeptonic("weighted") << endl;
	cout << "getBefore o" << countersBefore_->getAllLeptonic("weighted") << endl;
	for (unsigned int i = 0; i < varBins_.size() - 1; i++) {
		std::stringstream tmp;
		tmp << varBins_.at(i);
		double b, s, d, m;
		int bbin = 2 * i + 1;
		b = getHist(pmodulename_, "mbg_binnedBkg", bbin);
		s = getHist(pmodulename_, "mbg_binnedSemiLep", bbin);
		d = getHist(pmodulename_, "mbg_binnedDiLep", bbin);
		m = getHist(pmodulename_, "mbg_binnedMultiLep", bbin);

		countersBefore_->setPureHadronic(tmp.str(), b);
		countersBefore_->setSemiLeptonic(tmp.str(), s);
		countersBefore_->setDiLeptonic(tmp.str(), d);
		countersBefore_->setMultiLeptonic(tmp.str(), m);
	}
}

double MatrixAnalyzer::getHist(TString dir, TString hist, int & bin) {
	TH1F *dummy;
	TString directory(dir);
	directory += "/";
	TString name;
	name = hist;
	directory += name;

	dummy = (TH1F*) f_->Get(directory);

	if (!dummy) {
		cerr << "WARNING:" << " Didn't find indicated hist" << " ["
				<< directory << "]" << endl;
		return 2.0;
	}
	return dummy->GetBinContent(bin);
}

int MatrixAnalyzer::numberOfmatchedMuons(
		edm::Handle<edm::View<pat::Muon> > recMuons,
		edm::Handle<edm::Association<reco::GenParticleCollection> > genMatch) {
	int matchedmu = 0;
	int size = recMuons->size();
	cout << "<pdgId, status, mother.pdgId>" << endl;
	cout << "ID: ";
	for (int i = 0; i < size; ++i) {
		edm::RefToBase<reco::Candidate> muonRef =
				recMuons->refAt(i)->originalObjectRef();
		reco::GenParticleRef genMuon = (*genMatch)[muonRef];

		int mid = 0;
		int st = 0;
		if (genMuon->numberOfMothers() > 0) {
			mid = genMuon->mother(0)->pdgId();
			st = genMuon->mother(0)->status();
			cout << "mothers " << genMuon->numberOfMothers() << endl;
			if (genMuon->mother(0)->numberOfMothers() > 0) {
				cout << genMuon->mother(0)->pdgId() << endl;
				cout << genMuon->mother(0)->status() << endl;
				cout << genMuon->mother(0)->numberOfMothers() << endl;
			}
		}

		cout << "<" << genMuon->pdgId() << "," << genMuon->status() << ","
				<< mid << "," << st << "> ";
		if (fabs(mid) == 13)
			matchedmu++;
	}
	cout << "\n";
	return matchedmu;
}

