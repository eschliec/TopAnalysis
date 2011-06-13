#include "TopAnalysis/TopAnalyzer/plugins/PartonMatchAnalyzer.h"

#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

PartonMatchAnalyzer::PartonMatchAnalyzer(const edm::ParameterSet& cfg):
  matches_(cfg.getParameter<std::vector<edm::InputTag> >("matches")),
  quarks_ (cfg.getParameter<std::vector<int > >("quarks")),
  jets_   (cfg.getParameter<std::vector<int > >("jets")),
  sumDR_  (cfg.getParameter<edm::InputTag>("sumDR" )),
  sumDPt_ (cfg.getParameter<edm::InputTag>("sumDPt")),
  hist_   (cfg.getParameter<std::string>("hist")),
  ref_ (cfg.getParameter<unsigned int>("ref")), 
  normQrks_(0), normJets_(0)
{
  if( ref_>matches_.size() )
    throw edm::Exception( edm::errors::Configuration, "ref match outof bounds of match vector" );
}

void
PartonMatchAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{
  // recieve all matches of interest
  std::vector<std::vector<int> > matches;
  for(std::vector<edm::InputTag>::const_iterator tag=matches_.begin();
      tag!=matches_.end(); ++tag){
    edm::Handle<std::vector<int> > match;
    evt.getByLabel(*tag, match);
    matches.push_back(*match);
  }

  // fill corr_ hist for each quark in the list of considered 
  // quarks; in case the list is empty fill it for each quark
  for(unsigned int q=0; q<matches[ref_].size(); ++q){
    if(quarks_.empty() || std::count(quarks_.begin(), quarks_.end(), (int)q)>0 ){
      for(unsigned int ref=0; ref<matches.size(); ++ref){
	for(unsigned int smp=0; smp<matches.size();++smp){
	  if( ( (matches[ref])[q]>=0 &&  (matches[smp])[q]>=0) || 
	      (!(matches[ref])[q]>=0 && !(matches[smp])[q]>=0) )
	    corr_->Fill( ref, smp, +1 );
	  else
	    corr_->Fill( ref, smp, -1 );
	}
      }
    }
  }
  
  // fill alg_ histograms
  for(unsigned int q=0; q<matches[ref_].size(); ++q){
    if(quarks_.empty() || std::count(quarks_.begin(), quarks_.end(), (int)q)>0 ){ 
      for(unsigned int alg=0; alg<matches.size(); ++alg){ 
	if( (matches[ref_])[q]>=0 && (matches[alg])[q]>=0 )
	  algQrks_->Fill(alg, +1);
      }
    }
    if(jets_.empty() || std::count(jets_.begin(), jets_.end(), (matches[ref_])[q])>0 ){ 
      for(unsigned int alg=0; alg<matches.size(); ++alg){ 
	if( (matches[ref_])[q]>=0 && (matches[alg])[q]>=0 )
	  algJets_->Fill(alg, +1);
      }
    }
  }


  // fill eff histograms
  for(unsigned int q=0; q<matches[ref_].size(); ++q){
    ++normQrks_;
    if( (matches[ref_])[q]>=0  )
      effQrks_->Fill(q, +1);
      purQrks_->Fill(q, +1);
  }
  for(unsigned int q=0; q<matches[ref_].size(); ++q){
    for(int j=0; j<6; ++j){ 
      if( (matches[ref_])[q]==j )
	effJets_->Fill(j, +1);
	purJets_->Fill(j, +1);
    }
  }  

  // fill mon histograms
  edm::Handle<double> sumDR;
  evt.getByLabel(sumDR_, sumDR);
  monDR_ ->Fill(*sumDR);

  edm::Handle<double> sumDPt;
  evt.getByLabel(sumDPt_, sumDPt);
  monDPt_->Fill(*sumDPt);
}

void 
PartonMatchAnalyzer::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  if( !fs )
    throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );

  NameScheme corr("corr"), mon("mon"), eff("eff"), pur("pur"), alg("alg");
  int bin = (int)matches_.size();
  double max = (double)matches_.size();
  if(hist_.empty()){
    corr_ = fs->make<TH2F>(corr.name("matches" ), corr.name("matches" ), bin, 0., max, bin, 0., max); 

    monDR_  = fs->make<TH1F>(mon.name(      "dr"  ), mon.name("dr"  ), 100, 0.,  1.); 
    monDPt_ = fs->make<TH1F>(mon.name(      "dpt" ), mon.name("dpt" ), 100, 0.,  1.); 
    algJets_= fs->make<TH1F>(alg.name(      "jets"), alg.name("jets"), bin, 0., max); 
    algQrks_= fs->make<TH1F>(alg.name(      "qrks"), alg.name("qrks"), bin, 0., max); 
    effJets_= fs->make<TH1F>(eff.name(      "jets"), eff.name("jets"),   6, 0.,  6.);
    purJets_= fs->make<TH1F>(pur.name(      "jets"), pur.name("jets"),   6, 0.,  6.);
    effQrks_= fs->make<TH1F>(eff.name(      "qrks"), eff.name("qrks"),   4, 0.,  4.);
    purQrks_= fs->make<TH1F>(pur.name(      "qrks"), pur.name("qrks"),   4, 0.,  4.);
  }
  else{
    ofstream file(hist_.c_str(), std::ios::out);
    corr_ = fs->make<TH2F>(corr.name(file, "matches" ), corr.name("matches" ), bin, 0., max, bin, 0., max); 

    monDR_  = fs->make<TH1F>(mon.name(file, "dr"  ), mon.name("dr"  ), 100, 0.,  1.); 
    monDPt_ = fs->make<TH1F>(mon.name(file, "dpt" ), mon.name("dpt" ), 100, 0.,  1.);
    algJets_= fs->make<TH1F>(alg.name(file, "jets"), alg.name("jets"), bin, 0., max); 
    algQrks_= fs->make<TH1F>(alg.name(file, "qrks"), alg.name("qrks"), bin, 0., max);  
    effJets_= fs->make<TH1F>(eff.name(file, "jets"), eff.name("jets"),   6, 0.,  6.);
    purJets_= fs->make<TH1F>(pur.name(file, "jets"), pur.name("jets"),   6, 0.,  6.);
    effQrks_= fs->make<TH1F>(eff.name(file, "qrks"), eff.name("qrks"),   4, 0.,  4.);
    purQrks_= fs->make<TH1F>(pur.name(      "qrks"), pur.name("qrks"),   4, 0.,  4.);
  }
}

void 
PartonMatchAnalyzer::endJob() 
{
  // get proper normalization of histograms
  corr_   ->Scale(1./corr_->GetBinContent(1, 1)   );
  algQrks_->Scale(1./algQrks_->GetBinContent(ref_));
  algJets_->Scale(1./algJets_->GetBinContent(ref_));  
  purQrks_->Scale(1./purQrks_->Integral());
  purJets_->Scale(1./purJets_->Integral());
  effQrks_->Scale(1./normQrks_);
  effJets_->Scale(1./normJets_);
}
