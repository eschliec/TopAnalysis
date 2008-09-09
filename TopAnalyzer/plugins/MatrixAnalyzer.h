#ifndef MatrixAnalyzer_h
#define MatrixAnalyzer_h
/**
 *  class:   MatrixAnalyzer.h
 * @author: Lukas Kreczko, Uni Hamburg (lkreczko@mail.desy.de)
 * version $Id: MatrixAnalyzer.h,v 1.9 2008/09/08 08:44:42 kreczko Exp $

 ________________________________________________________________**/
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TopAnalysis/TopUtils/interface/NameScheme.h"
#include "TopAnalysis/TopUtils/interface/RootSystem.h"
#include "TopAnalysis/TopUtils/interface/RootHistograms.h"
#include "TopAnalysis/TopAnalyzer/plugins/LeptonCounter.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include <TGraphAsymmErrors.h>

class MatrixAnalyzer: public edm::EDAnalyzer {

public:
	explicit MatrixAnalyzer(const edm::ParameterSet&);
	~MatrixAnalyzer();

private:

	virtual void beginJob(const edm::EventSetup&);
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob();
	virtual void setEnv();
	virtual void getBefore();
	virtual TH1F* getHist(TString, TString);

	template<typename T1, typename T2> void log(T1 msg, T2 from, bool debug) {
		if (debug && debug_ || !debug)
			cout << from << ": " << msg << endl;
	}
	int numberOfmatchedMuons(edm::Handle<edm::View<pat::Muon> >,
			edm::Handle<edm::Association<reco::GenParticleCollection> >);
	std::string hist_, pmodulename_;
	bool before_, debug_;
	int noBins_, beforeBin_, afterBin_;
	double sampleweight_;
	edm::InputTag muons_;
	edm::InputTag var_;
	std::vector<double> varBins_;
	std::string notNeededHists_;
	typedef std::vector<pat::Muon> TopMuonCollection;

	std::map<int, TopMuonCollection> mothermap_;

	TH1F *background_, *binnedBkg_, *binnedSimpleBkg_;

	TH1F *lep_, *binnedSemiLep_, *binnedSimpleSemiLep_;

	TH1F *llep_, *binnedDiLep_, *binnedSimpleDiLep_;

	TH1F * multilep_, *binnedMultiLep_, *binnedSimpleMultiLep_;

	TH1F *overall_, *binnedOverall_, *binnedSimpleOverall_;

	TH1F *eff_, *binnedEff_;
	TGraphAsymmErrors *eff2_, *binnedEff2_;

	TFile *f_;

	TH1F *varPlot_;
	LeptonCounter *Counters_, *countersBefore_, *effCounter_, *effErrors_;
	//TODO: try to find a way to store double directly in root file
};
#endif
