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
	notNeededHists_ = "QCDnotNeeded.hist";
	cout << "module " << pmodulename_ << endl;


	Counters_ = new LeptonCounter();
	countersBefore_ = new LeptonCounter();
	effCounter_ = new LeptonCounter();
	effErrors_ =  new LeptonCounter();
	for (unsigned int x = 0; x < varBins_.size() - 1; x++) {
		std::stringstream tmp;
		tmp << varBins_.at(x);
		Counters_->addCounter(tmp.str());
		Counters_->addCounter(tmp.str() + "simple");
		countersBefore_->addCounter(tmp.str());
		countersBefore_->addCounter(tmp.str() + "simple");
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

	double *t = &varBins_[0];
	for(unsigned int i = 0; i< varBins_.size();i++){
		cout << varBins_[i] << endl;
		cout << t[i] << endl;

	}
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
	ofstream non(notNeededHists_.c_str(), std::ios::out);

	NameScheme nam("mbg");

	//make histograms
	background_ = fs->make<TH1F> (nam.name(non, "N_background"), nam.name(
			"background"), 5, 0, 5.);
	lep_ = fs->make<TH1F> (nam.name(non, "N_semilepton"), nam.name("lepton"),
			5, 0, 5.);
	llep_ = fs->make<TH1F> (nam.name(non, "N_dilepton"), nam.name("ll"), 5, 0, 5.);
	multilep_ = fs->make<TH1F> (nam.name(non, "N_multilepton"),
			nam.name("ml"), 5, 0, 5.);
	overall_ = fs->make<TH1F> (nam.name(non, "N_overall"),
			nam.name("overall"), noBins_, 0, noBins_);

	if(!before_) eff_ = fs->make<TH1F> (nam.name(hist, "efficency"), nam.name("eff"), varBins_.size()-1, &varBins_[0]);

	binnedBkg_ = fs->make<TH1F> (nam.name(non, "binnedBkg"), nam.name(
			"binnedBkg"), varBins_.size() -1, &varBins_[0]);
	binnedSemiLep_ = fs->make<TH1F> (nam.name(non, "binnedSemiLep"), nam.name(
			"binnedSemiLep"), varBins_.size() -1, &varBins_[0]);
	binnedDiLep_ = fs->make<TH1F> (nam.name(non, "binnedDiLep"), nam.name(
			"binnedDiLep"), varBins_.size() -1, &varBins_[0]);
	binnedMultiLep_ = fs->make<TH1F> (nam.name(non, "binnedMultiLep"),
			nam.name("binnedMultiLep"), varBins_.size() -1, &varBins_[0]);
	if(!before_) binnedEff_ = fs->make<TH1F> (nam.name(hist, "binnedEff"), nam.name(
			"binnedEff"), varBins_.size() -1, &varBins_[0]);
	binnedOverall_ = fs->make<TH1F> (nam.name(non, "binnedOverall"), nam.name(
			"binnedOverall"), varBins_.size() -1, &varBins_[0]);

	binnedSimpleBkg_ = fs->make<TH1F> (nam.name(non, "binnedSimpleBkg"),
			nam.name("binnedSimpleBkg"), varBins_.size() -1, &varBins_[0]);
	binnedSimpleSemiLep_ = fs->make<TH1F> (
			nam.name(non, "binnedSimpleSemiLep"), nam.name(
					"binnedSimpleSemiLep"), varBins_.size() -1, &varBins_[0]);
	binnedSimpleDiLep_ = fs->make<TH1F> (nam.name(non, "binnedSimpleDiLep"),
			nam.name("binnedSimpleDiLep"), varBins_.size() -1, &varBins_[0]);
	binnedSimpleMultiLep_ = fs->make<TH1F> (nam.name(non,
			"binnedSimpleMultiLep"), nam.name("binnedSimpleMultiLep"),
			varBins_.size() -1, &varBins_[0]);
	binnedSimpleOverall_ = fs->make<TH1F> (
			nam.name(non, "binnedSimpleOverall"), nam.name(
					"binnedSimpleOverall"), varBins_.size() -1, &varBins_[0]);

	varPlot_ = fs->make<TH1F> (nam.name(hist, var_.label().c_str()), nam.name(var_.label().c_str()),
			250, 0, 500);
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
		if (pt >= varBins_[i] && (pt < varBins_[i + 1] || i == varBins_.size() - 2)) {
			double t = varBins_.at(i);
			std::stringstream tmp;
			tmp << t;
			if (muons->size() < 1) {
				Counters_->addPureHadronic(tmp.str() + "simple");
				Counters_->addPureHadronic(tmp.str(), sampleweight_);
			}
			else if (muons->size() == 1) {
				Counters_->addSemiLeptonic(tmp.str() + "simple");
				Counters_->addSemiLeptonic(tmp.str(), sampleweight_);
			}
			else if (muons->size() == 2) {
				Counters_->addDiLeptonic(tmp.str() + "simple");
				Counters_->addDiLeptonic(tmp.str(),sampleweight_);
			}
			else if (muons->size() > 2) {
				Counters_->addMultiLeptonic(tmp.str() + "simple");
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
		setEnv();
		cout << endl;
	} else {
		getBefore();
		cout << "#----- after (" << "Module: " << pmodulename_ <<")-----#" << endl;
		log(sampleweight_, "::endJob >> sampleweight", false);
		setEnv();
		cout.precision(10);
		if (countersBefore_->getPureHadronic("weighted") != 0 and countersBefore_->getPureHadronic("simple") !=0) {
			double eff, err;
			eff = Counters_->getPureHadronic("weighted") / countersBefore_->getPureHadronic("simple");
			err = sqrt(eff * (1 - eff) / countersBefore_->getPureHadronic("simple"));
			cout << "hadEff " << eff << " +- " << err << endl;
			effCounter_->setPureHadronic("weighted", eff);
			effErrors_->setPureHadronic("weighted", err);
		}

		if (countersBefore_->getSemiLeptonic("weighted") != 0 and countersBefore_->getSemiLeptonic("simple") !=0) {
			double eff, err;
			eff = Counters_->getSemiLeptonic("weighted") / countersBefore_->getSemiLeptonic("weighted");
			err = sqrt(eff * (1 - eff) / countersBefore_->getSemiLeptonic("simple"));
			cout << "semiEff " << eff << " +- " << err << endl;
			effCounter_->setSemiLeptonic("weighted", eff);
			effErrors_->setSemiLeptonic("weighted", err);
		}
		if (countersBefore_->getDiLeptonic("weighted") != 0 and countersBefore_->getDiLeptonic("simple") !=0) {
			double eff, err;
			eff = Counters_->getDiLeptonic("weighted") / countersBefore_->getDiLeptonic("weighted");
			err = sqrt(eff * (1 - eff) / countersBefore_->getDiLeptonic("simple"));
			cout << "diEff " << eff << " +- " << err << endl;
			effCounter_->setDiLeptonic("weighted", eff);
			effErrors_->setDiLeptonic("weighted", err);
		}
		if (countersBefore_->getMultiLeptonic("weighted") != 0 and countersBefore_->getMultiLeptonic(
				"simple") !=0) {
			double eff, err;
			eff = Counters_->getMultiLeptonic("weighted") / countersBefore_->getMultiLeptonic("weighted");
			err = sqrt(eff * (1 - eff) / countersBefore_->getMultiLeptonic("simple"));
			cout << "multiEff " << eff << " +- " << err << endl;
			effCounter_->setMultiLeptonic("weighted", eff);
			effErrors_->setMultiLeptonic("weighted", err);
		}

		for (unsigned int i = 0; i < varBins_.size() - 1; i++) {
			std::stringstream tmp;
			tmp << varBins_.at(i);
			double b, s, d, m, bb, sb, db, mb, o, ob;
			double sbb, ssb, sdb, smb, sob;
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

			sbb = countersBefore_->getPureHadronic(tmp.str() + "simple");
			ssb = countersBefore_->getSemiLeptonic(tmp.str() + "simple");
			sdb = countersBefore_->getDiLeptonic(tmp.str() + "simple");
			smb = countersBefore_->getMultiLeptonic(tmp.str() + "simple");
			sob = countersBefore_->getAllLeptonic(tmp.str() + "simple");

			if (bb != 0. and sbb != 0.) {
				double eff,err;
				eff = b / bb;
				err = sqrt(eff * (1 - eff) / sbb);
				effCounter_->setPureHadronic(tmp.str(), eff);
				effErrors_->setPureHadronic(tmp.str(), err);
//				cout << tmp.str() << " (had) " << eff << " +- " << err << endl;
			}
			if (sb != 0. and ssb != 0.) {
				double eff,err;
				eff = s/sb;
				err = sqrt(eff * (1 - eff) / ssb);
				effCounter_->setSemiLeptonic(tmp.str(), eff);
				effErrors_->setSemiLeptonic(tmp.str(), err);
				cout << tmp.str() << " (semi) " << eff << " +- " << err << endl;
			}
			if (db != 0. and sdb != 0.) {
				double eff, err;
				eff = d / db;
				err = sqrt(eff * (1 - eff) / sdb);
				effCounter_->setDiLeptonic(tmp.str(),eff);
				effErrors_->setDiLeptonic(tmp.str(), err);
//				cout << tmp.str() << " (di) " << eff << " +- " << err << endl;
			}
			if (mb != 0. and smb != 0.) {
				double eff, err;
				eff = m / mb;
				err = sqrt(eff * (1 - eff) / smb);
				effCounter_->setMultiLeptonic(tmp.str(), eff);
				effErrors_->setMultiLeptonic(tmp.str(), err);
//				cout << tmp.str() << " (multi) " << eff << " +- " << err << endl;
			}

			if (ob != 0. and sob != 0.) {
				double eff, err;
				eff = o / ob;
				err = sqrt(eff * (1 - eff) / sob);
				effCounter_->setSemiLeptonic(tmp.str(), eff);
				effErrors_->setSemiLeptonic(tmp.str(), err);
				cout << tmp.str() << " (overall) " << eff << " +- " << err
						<< endl;
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

//		background_->SetBinContent(beforeBin_, countersBefore_->getPureHadronic("weighted"));
//		lep_->SetBinContent(beforeBin_, countersBefore_->getSemiLeptonic("weighted"));
//		llep_->SetBinContent(beforeBin_, countersBefore_->getDiLeptonic("weighted"));
//		multilep_->SetBinContent(beforeBin_, countersBefore_->getMultiLeptonic("weighted"));
//		overall_->SetBinContent(beforeBin_, countersBefore_->getAllLeptonic("weighted"));


	}
	//TODO: Fill in histogram.
//	for (map<int, TopMuonCollection>::iterator iter = mothermap_.begin(); iter
//			!= mothermap_.end(); iter++) {
//		cout << "total matched mu from " << (*iter).first << " : "
//				<< (*iter).second.size() << endl;
//	}

}

void MatrixAnalyzer::setEnv() {
	if(!before_){
//		cout << "type >> events with muons (after) <> weighted events with muons (before) <> gen. events" << endl;
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
		double sbb, ssb, sdb, smb, sob;
		double _sb, _ss, _sd, _sm, _so;
		int bbin = i+1;

		b = Counters_->getPureHadronic(tmp.str());
		s = Counters_->getSemiLeptonic(tmp.str());
		d = Counters_->getDiLeptonic(tmp.str());
		m = Counters_->getMultiLeptonic(tmp.str());
		o = Counters_->getAllLeptonic(tmp.str());

		_sb = Counters_->getPureHadronic(tmp.str() + "simple");
		_ss = Counters_->getSemiLeptonic(tmp.str() + "simple");
		_sd = Counters_->getDiLeptonic(tmp.str() + "simple");
		_sm = Counters_->getMultiLeptonic(tmp.str() + "simple");
		_so = Counters_->getAllLeptonic(tmp.str() + "simple");

		bb = countersBefore_->getPureHadronic(tmp.str());
		sb = countersBefore_->getSemiLeptonic(tmp.str());
		db = countersBefore_->getDiLeptonic(tmp.str());
		mb = countersBefore_->getMultiLeptonic(tmp.str());
		ob = countersBefore_->getAllLeptonic(tmp.str());

		sbb = countersBefore_->getPureHadronic(tmp.str() + "simple");
		ssb = countersBefore_->getSemiLeptonic(tmp.str() + "simple");
		sdb = countersBefore_->getDiLeptonic(tmp.str() + "simple");
		smb = countersBefore_->getMultiLeptonic(tmp.str() + "simple");
		sob = countersBefore_->getAllLeptonic(tmp.str() + "simple");

		cout << varBins_[i] << " (semi):" << s << " <> " << sb << " <> " << _sb << endl;
		cout << varBins_[i] << " (overall):" << o << " <> " << ob << " <> " << _so << endl;
//		cout << varBins_[i] << " (semi):" << s << " <> " << sb << endl;
//		cout << varBins_[i] << " (overall):" << o << " <> " << ob <<  endl;
//		if(before_){
			binnedBkg_->SetBinContent(bbin, b);
			binnedSemiLep_->SetBinContent(bbin, s);
			binnedDiLep_->SetBinContent(bbin, d);
			binnedMultiLep_->SetBinContent(bbin, m);
			binnedOverall_->SetBinContent(bbin, o);

			binnedSimpleBkg_->SetBinContent(bbin, _sb);
			binnedSimpleSemiLep_->SetBinContent(bbin, _ss);
			binnedSimpleDiLep_->SetBinContent(bbin, _sd);
			binnedSimpleMultiLep_->SetBinContent(bbin, _sm);
			binnedSimpleOverall_->SetBinContent(bbin, _so);
//		}
//		else {
//			binnedBkg_->SetBinContent(2 * i + 2, b);
//			binnedSemiLep_->SetBinContent(2 * i + 2, s);
//			binnedDiLep_->SetBinContent(2 * i + 2, d);
//			binnedMultiLep_->SetBinContent(2 * i + 2, m);
//			binnedOverall_->SetBinContent(2 * i + 2, o);
//
//			binnedBkg_->SetBinContent(2 * i + 1, bb);
//			binnedSemiLep_->SetBinContent(2 * i + 1, sb);
//			binnedDiLep_->SetBinContent(2 * i + 1, db);
//			binnedMultiLep_->SetBinContent(2 * i + 1, mb);
//			binnedOverall_->SetBinContent(2 * i + 1, ob);
//
//			binnedSimpleBkg_->SetBinContent(2 * i + 2, _sb);
//			binnedSimpleSemiLep_->SetBinContent(2 * i + 2, _ss);
//			binnedSimpleDiLep_->SetBinContent(2 * i + 2, _sd);
//			binnedSimpleMultiLep_->SetBinContent(2 * i + 2, _sm);
//			binnedSimpleOverall_->SetBinContent(2 * i + 2, _so);
//
//			binnedSimpleBkg_->SetBinContent(2 * i + 1, sbb);
//			binnedSimpleSemiLep_->SetBinContent(2 * i + 1, ssb);
//			binnedSimpleDiLep_->SetBinContent(2 * i + 1, sdb);
//			binnedSimpleMultiLep_->SetBinContent(2 * i + 1, smb);
//			binnedSimpleOverall_->SetBinContent(2 * i + 1, sob);
//
//		}
//		binnedBkg_->GetXaxis()->SetBinLabel(bbin, tmp.str().c_str());
////		binnedBkg_->GetXaxis()->SetBinLabel(2 * i + 2, "A");
//
//		binnedSemiLep_->GetXaxis()->SetBinLabel(bbin, tmp.str().c_str());
////		binnedSemiLep_->GetXaxis()->SetBinLabel(2 * i + 2, "A");
//
//		binnedDiLep_->GetXaxis()->SetBinLabel(bbin, tmp.str().c_str());
////		binnedDiLep_->GetXaxis()->SetBinLabel(2 * i + 2, "A");
//
//		binnedMultiLep_->GetXaxis()->SetBinLabel(bbin, tmp.str().c_str());
////		binnedMultiLep_->GetXaxis()->SetBinLabel(2 * i + 2, "A");
//
//		binnedOverall_->GetXaxis()->SetBinLabel(bbin, tmp.str().c_str());
////		binnedOverall_->GetXaxis()->SetBinLabel(2 * i + 2, "A");

	}

	background_->SetBinContent(2, Counters_->getPureHadronic("weighted"));
	lep_->SetBinContent(2, Counters_->getSemiLeptonic("weighted"));
	llep_->SetBinContent(2, Counters_->getDiLeptonic("weighted"));
	multilep_->SetBinContent(2, Counters_->getMultiLeptonic("weighted"));
	overall_->SetBinContent(2, Counters_->getAllLeptonic("weighted"));
//	if (before_) {
//		background_->GetXaxis()->SetBinLabel(bin, "before");
//		lep_->GetXaxis()->SetBinLabel(bin, "before");
//		llep_->GetXaxis()->SetBinLabel(bin, "before");
//		multilep_->GetXaxis()->SetBinLabel(bin, "before");
//		overall_->GetXaxis()->SetBinLabel(bin, "before");
//	}
//	else {
//		background_->GetXaxis()->SetBinLabel(bin, "after");
//		lep_->GetXaxis()->SetBinLabel(bin, "after");
//		llep_->GetXaxis()->SetBinLabel(bin, "after");
//		multilep_->GetXaxis()->SetBinLabel(bin, "after");
//		overall_->GetXaxis()->SetBinLabel(bin, "after");
//	}

}

void MatrixAnalyzer::getBefore() {
	countersBefore_->setPureHadronic("weighted", getHist(pmodulename_, "mbg_N_background", 2));
	countersBefore_->setSemiLeptonic("weighted", getHist(pmodulename_, "mbg_N_semilepton", 2));
	countersBefore_->setDiLeptonic("weighted", getHist(pmodulename_, "mbg_N_dilepton", 2));
	countersBefore_->setMultiLeptonic("weighted",getHist(pmodulename_, "mbg_N_multilepton", 2));

	cout << "getBefore s" << countersBefore_->getSemiLeptonic("weighted") << endl;
	cout << "getBefore o" << countersBefore_->getAllLeptonic("weighted") << endl;
	for (unsigned int i = 0; i < varBins_.size() - 1; i++) {
		std::stringstream tmp;
		tmp << varBins_.at(i);
		double b, s, d, m;
		double sb, ss, sd, sm;
		int bbin = i + 1;
		b = getHist(pmodulename_, "mbg_binnedBkg", bbin);
		s = getHist(pmodulename_, "mbg_binnedSemiLep", bbin);
		d = getHist(pmodulename_, "mbg_binnedDiLep", bbin);
		m = getHist(pmodulename_, "mbg_binnedMultiLep", bbin);

		sb = getHist(pmodulename_, "mbg_binnedSimpleBkg", bbin);
		ss = getHist(pmodulename_, "mbg_binnedSimpleSemiLep", bbin);
		sd = getHist(pmodulename_, "mbg_binnedSimpleDiLep", bbin);
		sm = getHist(pmodulename_, "mbg_binnedSimpleMultiLep", bbin);

		countersBefore_->setPureHadronic(tmp.str(), b);
		countersBefore_->setSemiLeptonic(tmp.str(), s);
		countersBefore_->setDiLeptonic(tmp.str(), d);
		countersBefore_->setMultiLeptonic(tmp.str(), m);

		countersBefore_->setPureHadronic(tmp.str() + "simple", sb);
		countersBefore_->setSemiLeptonic(tmp.str() + "simple", ss);
		countersBefore_->setDiLeptonic(tmp.str() + "simple", sd);
		countersBefore_->setMultiLeptonic(tmp.str() + "simple", sm);
	}
}

double MatrixAnalyzer::getHist(TString dir, TString hist, int bin) {
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
	cout << "getHist ( " << bin << " ): " << directory << " = " <<  dummy->GetBinContent(bin) << endl;
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
		edm::Ptr<reco::Candidate> muonPtr =
				recMuons->refAt(i)->originalObjectRef();
		reco::GenParticleRef genMuon = (*genMatch)[muonPtr];

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

