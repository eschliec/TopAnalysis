#ifndef SemiLepHypothesesAnalyzer_h
#define SemiLepHypothesesAnalyzer_h

#include "TH1F.h"
#include "TH2F.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "AnalysisDataFormats/TopObjects/interface/TtSemiLeptonicEvent.h"

class SemiLepHypothesesAnalyzer : public edm::EDAnalyzer {

 public:

  explicit SemiLepHypothesesAnalyzer(const edm::ParameterSet&);
  ~SemiLepHypothesesAnalyzer(){};
  
 private:

  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  edm::InputTag semiLepEvt_;
  edm::InputTag hypoKey_;
  edm::InputTag wgt_;
  unsigned int nJetsMax_;
  double maxSumDRGenMatch_;
  double minProbKinFit_;
  double minMVADisc_; 
  std::string hist_;

  void bookKinHistos      (edm::Service<TFileService>&, ofstream&);
  void bookKinResHistos   (edm::Service<TFileService>&, ofstream&);
  void bookJetCorrelHistos(edm::Service<TFileService>&, ofstream&);
  void bookQualityHistos  (edm::Service<TFileService>&, ofstream&);

  void fillKinHistos      (std::vector<TH1F*>&, const reco::Candidate&, const double&);
  void fillKinResHistos   (std::vector<TH1F*>&, const reco::Candidate&, const reco::Candidate&, const double&);
  void fillJetCorrelHistos(const std::vector<int>&, const std::vector<int>&, const double&);
  void fillQualityHistos  (const TtSemiLeptonicEvent&, const double&);

  std::vector<TH1F*> hadTopKin_;
  std::vector<TH1F*> hadWKin_;
  std::vector<TH1F*> hadBKin_;
  std::vector<TH1F*> hadQKin_;
  std::vector<TH1F*> hadPKin_;
  std::vector<TH1F*> lepTopKin_;
  std::vector<TH1F*> lepWKin_;
  std::vector<TH1F*> lepBKin_;
  std::vector<TH1F*> leptonKin_;
  std::vector<TH1F*> neutriKin_;

  TH1F* hadTopLepTopMassDiff;
  TH1F* hadWLepWMassDiff;

  std::vector<TH1F*> hadTopKinRes_;
  std::vector<TH1F*> hadWKinRes_;
  std::vector<TH1F*> hadBKinRes_;
  std::vector<TH1F*> hadQKinRes_;
  std::vector<TH1F*> hadPKinRes_;
  std::vector<TH1F*> lepTopKinRes_;
  std::vector<TH1F*> lepWKinRes_;
  std::vector<TH1F*> lepBKinRes_;
  std::vector<TH1F*> leptonKinRes_;
  std::vector<TH1F*> neutriKinRes_;

  TH2F* hadBJetCorrel_;
  TH2F* hadQJetCorrel_;
  TH2F* hadPJetCorrel_;
  TH2F* lepBJetCorrel_;

  TH1F* goodHypo_;

  TH1F* genMatchSumDR_;
  TH1F* genMatchSumPt_;
  TH1F* mvaDisc_;
  TH1F* fitChi2_;
  TH1F* fitProb_;

  TH2F* genMatchSumDRVsSumPt_;
  TH2F* genMatchSumDRVsHadWMass_;
  TH2F* genMatchSumDRVsHadTopMass_;
  TH2F* genMatchSumDRVsMVADisc_;
  TH2F* genMatchSumDRVsFitProb_;
  TH2F* mvaDiscVsHadWMass_;
  TH2F* mvaDiscVsHadTopMass_;
  TH2F* mvaDiscVsFitProb_;
  TH2F* fitProbVsHadWMass_;
  TH2F* fitProbVsHadTopMass_;

};

#endif
