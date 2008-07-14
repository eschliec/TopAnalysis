#ifndef LeptonCounter_h
#define LeptonCounter_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "TopAnalysis/TopUtils/interface/RootSystem.h"
#include "TopAnalysis/TopUtils/interface/RootHistograms.h"

#include "DataFormats/PatCandidates/interface/Muon.h"


class LeptonCounter : public edm::EDAnalyzer {
public:
	explicit LeptonCounter(const edm::ParameterSet&);
	~LeptonCounter();
	
private:
	virtual void beginJob(const edm::EventSetup&);
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob();
	
	edm::InputTag muons_;
	
	int eleCounter_, muCounter_, numberOfRatioBins_, numberOfMuonBins_,
			numberOfElecBins_ ;
	double minRatio_, maxRatio_, minNmuon_, maxNmuon_,minNelec_, maxNelec_;
	
	
	typedef std::vector<pat::Muon> TopMuonCollection;
	
	TH1F *muonElecRatio_, *numberOfMuons_, *numberOfElecs_;
};
#endif
