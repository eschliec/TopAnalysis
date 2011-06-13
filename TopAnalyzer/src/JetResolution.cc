#include "TopAnalysis/TopAnalyzer/interface/JetResolution.h"


/// constructor for FWLite analyzer
JetResolution::JetResolution(int nJets, double matchDR, std::vector<double> binsPt):
  nJets_(nJets), binsPt_(binsPt), matchDR_(matchDR)  
{
}

/// constructor for full FW analyzer
JetResolution::JetResolution(const edm::ParameterSet& cfg):
  nJets_( cfg.getParameter<int> ( "nJets" ) ),
  binsPt_( cfg.getParameter<std::vector<double> >( "binsPt" ) ),
  binsEta_( cfg.getParameter<std::vector<double> >( "binsEta") ),
  binsPhi_( cfg.getParameter<std::vector<double> >( "binsPhi") ),
  matchDR_( cfg.getParameter<double>( "matchDR" ) )
{
}

/// fill interface for full FW analyzer
void
JetResolution::fill(const edm::Event& evt, const std::vector<pat::Jet>& jets, const double& weight)
{
  fill(jets, weight);
}

/// fill interface for FWLite analyzer
void
JetResolution::fill(const std::vector<pat::Jet>& jets, const double& weight)
{
  unsigned int idx=0;
  for(std::vector<pat::Jet>::const_iterator jet = jets.begin(); 
      jet!=jets.end(); ++jet, ++idx) {
    if( jet->genJet() ){ 
      double dR=deltaR( jet->eta(), jet->phi(), jet->genJet()->eta(), jet->genJet()->phi() ); 
      if(dR<matchDR_){
	double dPt=(jet->pt()-jet->genJet()->pt())/(jet->genJet())->pt();
        double dPhi=deltaPhi( jet->phi(), jet->genJet()->phi() );
        double dEta=(jet->eta()-jet->genJet()->eta());
	for(unsigned int jdx=0; jdx<(binsPt_.size()-1); ++jdx) {
	  if( (binsPt_[jdx]<jet->genJet()->pt()) && (jet->genJet()->pt()<binsPt_[jdx+1]) ){
	    relPtAll_[jdx]->Fill(dPt, weight);
	    if( idx< relPtJet_.size() ){
	      (relPtJet_[idx])[jdx]->Fill(dPt, weight);
	    }
	  }
	}
	for(unsigned int jdx=0; jdx<(binsEta_.size()-1); ++jdx) {
	  if( (binsEta_[jdx]<jet->genJet()->eta()) && (jet->genJet()->eta()<binsEta_[jdx+1]) ){
	    EtaAll_[jdx]->Fill(dEta, weight);
	    if( idx< EtaJet_.size() ){
	      (EtaJet_[idx])[jdx]->Fill(dEta, weight);
	    }
	  }
	}
	for(unsigned int jdx=0; jdx<(binsPhi_.size()-1); ++jdx) {
	  if( (binsPhi_[jdx]<jet->genJet()->phi()) && (jet->genJet()->phi()<binsPhi_[jdx+1]) ){
	    PhiAll_[jdx]->Fill(dPhi, weight);
	    if( idx< PhiJet_.size() ){
	      (PhiJet_[idx])[jdx]->Fill(dPhi, weight);
	    }
	  }
	}
      }
    }
  }
}

/// book for FWLite
void
JetResolution::book()
{
  NameScheme fit("fit"), cal("cal"), res("res");
  for(int jdx=0; jdx<nJets_; ++jdx){
    std::vector<TH1F*> buffer;
    calPtJet_.push_back(new TH1F(cal.name("relPt",jdx), cal.name("relPt",jdx), ((int)binsPt_.size()-1), &binsPt_[0]));
    resPtJet_.push_back(new TH1F(res.name("relPt",jdx), res.name("relPt",jdx), ((int)binsPt_.size()-1), &binsPt_[0]));
    for(int idx=0; idx<((int)binsPt_.size()-1); ++idx)
      buffer.push_back( new TH1F(fit.name("relPt",jdx,idx), fit.name("relPt",jdx,idx), 100, -1., 1.) );
    relPtJet_.push_back(buffer); 
  }
  for(int idx = 0; idx < ((int)binsPt_.size()-1); ++idx) {
    relPtAll_.push_back( new TH1F(fit.name("relPt",idx), fit.name("relPt",idx), 100, -1., 1.) );
  }
  calPtAll_= new TH1F(cal.name("relPt"), cal.name("relPt"), ((int)binsPt_.size()-1), &binsPt_[0]);
  resPtAll_= new TH1F(res.name("relPt"), res.name("relPt"), ((int)binsPt_.size()-1), &binsPt_[0]);

  for(int jdx=0; jdx<nJets_; ++jdx){
    std::vector<TH1F*> buffer;
    calEtaJet_.push_back(new TH1F(cal.name("Eta",jdx), cal.name("Eta",jdx), ((int)binsEta_.size()-1), &binsEta_[0]));
    resEtaJet_.push_back(new TH1F(res.name("Eta",jdx), res.name("Eta",jdx), ((int)binsEta_.size()-1), &binsEta_[0]));
    for(int idx=0; idx<((int)binsEta_.size()-1); ++idx)
      buffer.push_back( new TH1F(fit.name("Eta",jdx,idx), fit.name("Eta",jdx,idx), 100, -0.5, 0.5) );
    EtaJet_.push_back(buffer); 
  }
  for(int idx = 0; idx < ((int)binsEta_.size()-1); ++idx) {
    EtaAll_.push_back( new TH1F(fit.name("Eta",idx), fit.name("Eta",idx), 100, -0.5, 0.5) );
  }
  calEtaAll_= new TH1F(cal.name("Eta"), cal.name("Eta"), ((int)binsEta_.size()-1), &binsEta_[0]);
  resEtaAll_= new TH1F(res.name("Eta"), res.name("Eta"), ((int)binsEta_.size()-1), &binsEta_[0]);

  for(int jdx=0; jdx<nJets_; ++jdx){
    std::vector<TH1F*> buffer;
    calPhiJet_.push_back(new TH1F(cal.name("Phi",jdx), cal.name("Phi",jdx), ((int)binsPhi_.size()-1), &binsPhi_[0]));
    resPhiJet_.push_back(new TH1F(res.name("Phi",jdx), res.name("Phi",jdx), ((int)binsPhi_.size()-1), &binsPhi_[0]));
    for(int idx=0; idx<((int)binsPhi_.size()-1); ++idx)
      buffer.push_back( new TH1F(fit.name("Phi",jdx,idx), fit.name("Phi",jdx,idx), 100, -0.5, 0.5) );
    PhiJet_.push_back(buffer); 
  }
  for(int idx = 0; idx < ((int)binsPhi_.size()-1); ++idx) {
    PhiAll_.push_back( new TH1F(fit.name("Phi",idx), fit.name("Phi",idx), 100, -0.5, 0.5) );
  }
  calPhiAll_= new TH1F(cal.name("Phi"), cal.name("Phi"), ((int)binsPhi_.size()-1), &binsPhi_[0]);
  resPhiAll_= new TH1F(res.name("Phi"), res.name("Phi"), ((int)binsPhi_.size()-1), &binsPhi_[0]);
}

/// book for full FW
void
JetResolution::book(edm::Service<TFileService>& fs)
{
  NameScheme fit("fit"), cal("cal"), res("res");
  for(int jdx=0; jdx<nJets_; ++jdx){
    std::vector<TH1F*> buffer;
    calPtJet_.push_back(fs->make<TH1F>(cal.name("relPt",jdx), cal.name("relPt",jdx), ((int)binsPt_.size()-1), &binsPt_[0]));
    resPtJet_.push_back(fs->make<TH1F>(res.name("relPt",jdx), res.name("relPt",jdx), ((int)binsPt_.size()-1), &binsPt_[0]));
    for(int idx=0; idx<((int)binsPt_.size()-1); ++idx)
      buffer.push_back( fs->make<TH1F>(fit.name("relPt",jdx,idx), fit.name("relPt",jdx,idx), 100, -0.5, 0.5) );
    relPtJet_.push_back(buffer); 
  }
  for(int idx = 0; idx < ((int)binsPt_.size()-1); ++idx) {
    relPtAll_.push_back( fs->make<TH1F>(fit.name("relPt",idx), fit.name("relPt",idx), 100, -0.5, 0.5) );
  }
  calPtAll_= fs->make<TH1F>(cal.name("relPt"), cal.name("relPt"), ((int)binsPt_.size()-1), &binsPt_[0]);
  resPtAll_= fs->make<TH1F>(res.name("relPt"), res.name("relPt"), ((int)binsPt_.size()-1), &binsPt_[0]);

  for(int jdx=0; jdx<nJets_; ++jdx){
    std::vector<TH1F*> buffer;
    calEtaJet_.push_back(fs->make<TH1F>(cal.name("Eta",jdx), cal.name("Eta",jdx), ((int)binsEta_.size()-1), &binsEta_[0]));
    resEtaJet_.push_back(fs->make<TH1F>(res.name("Eta",jdx), res.name("Eta",jdx), ((int)binsEta_.size()-1), &binsEta_[0]));
    for(int idx=0; idx<((int)binsEta_.size()-1); ++idx)
      buffer.push_back( fs->make<TH1F>(fit.name("Eta",jdx,idx), fit.name("Eta",jdx,idx), 100, -0.5, 0.5) );
    EtaJet_.push_back(buffer); 
  }
  for(int idx = 0; idx < ((int)binsEta_.size()-1); ++idx) {
    EtaAll_.push_back( fs->make<TH1F>(fit.name("Eta",idx), fit.name("Eta",idx), 100, -0.5, 0.5) );
  }
  calEtaAll_= fs->make<TH1F>(cal.name("Eta"), cal.name("Eta"), ((int)binsEta_.size()-1), &binsEta_[0]);
  resEtaAll_= fs->make<TH1F>(res.name("Eta"), res.name("Eta"), ((int)binsEta_.size()-1), &binsEta_[0]);

  for(int jdx=0; jdx<nJets_; ++jdx){
    std::vector<TH1F*> buffer;
    calPhiJet_.push_back(fs->make<TH1F>(cal.name("Phi",jdx), cal.name("Phi",jdx), ((int)binsPhi_.size()-1), &binsPhi_[0]));
    resPhiJet_.push_back(fs->make<TH1F>(res.name("Phi",jdx), res.name("Phi",jdx), ((int)binsPhi_.size()-1), &binsPhi_[0]));
    for(int idx=0; idx<((int)binsPhi_.size()-1); ++idx)
      buffer.push_back( fs->make<TH1F>(fit.name("Phi",jdx,idx), fit.name("Phi",jdx,idx), 100, -0.5, 0.5) );
    PhiJet_.push_back(buffer); 
  }
  for(int idx = 0; idx < ((int)binsPhi_.size()-1); ++idx) {
    PhiAll_.push_back( fs->make<TH1F>(fit.name("Phi",idx), fit.name("Phi",idx), 100, -0.5, 0.5) );
  }
  calPhiAll_= fs->make<TH1F>(cal.name("Phi"), cal.name("Phi"), ((int)binsPhi_.size()-1), &binsPhi_[0]);
  resPhiAll_= fs->make<TH1F>(res.name("Phi"), res.name("Phi"), ((int)binsPhi_.size()-1), &binsPhi_[0]);
}

/// book for full FW with output stream
void
JetResolution::book(edm::Service<TFileService>& fs, ofstream& file)
{
  NameScheme fit("fit"), cal("cal"), res("res");
  for(int jdx=0; jdx<nJets_; ++jdx){
    std::vector<TH1F*> buffer;
    calPtJet_.push_back(fs->make<TH1F>(cal.name(file,"relPt",jdx), cal.name("relPt",jdx), ((int)binsPt_.size()-1), &binsPt_[0]));
    resPtJet_.push_back(fs->make<TH1F>(res.name(file,"relPt",jdx), res.name("relPt",jdx), ((int)binsPt_.size()-1), &binsPt_[0]));
    for(int idx=0; idx<((int)binsPt_.size()-1); ++idx)
      buffer.push_back( fs->make<TH1F>(fit.name(file,"relPt",jdx,idx), fit.name("relPt",jdx,idx), 100, -0.5, 0.5) );
    relPtJet_.push_back(buffer); 
  }
  for(int idx = 0; idx < ((int)binsPt_.size()-1); ++idx) {
    relPtAll_.push_back( fs->make<TH1F>(fit.name(file,"relPt",idx), fit.name("relPt",idx), 100, -0.5, 0.5) );
  }
  calPtAll_= fs->make<TH1F>(cal.name(file,"relPt"), cal.name("relPt"), ((int)binsPt_.size()-1), &binsPt_[0]);
  resPtAll_= fs->make<TH1F>(res.name(file,"relPt"), res.name("relPt"), ((int)binsPt_.size()-1), &binsPt_[0]);

  for(int jdx=0; jdx<nJets_; ++jdx){
    std::vector<TH1F*> buffer;
    calEtaJet_.push_back(fs->make<TH1F>(cal.name(file,"Eta",jdx), cal.name("Eta",jdx), ((int)binsEta_.size()-1), &binsEta_[0]));
    resEtaJet_.push_back(fs->make<TH1F>(res.name(file,"Eta",jdx), res.name("Eta",jdx), ((int)binsEta_.size()-1), &binsEta_[0]));
    for(int idx=0; idx<((int)binsEta_.size()-1); ++idx)
      buffer.push_back( fs->make<TH1F>(fit.name(file,"Eta",jdx,idx), fit.name("Eta",jdx,idx), 100, -0.5, 0.5) );
    EtaJet_.push_back(buffer); 
  }
  for(int idx = 0; idx < ((int)binsEta_.size()-1); ++idx) {
    EtaAll_.push_back( fs->make<TH1F>(fit.name(file,"Eta",idx), fit.name("Eta",idx), 100, -0.5, 0.5) );
  }
  calEtaAll_= fs->make<TH1F>(cal.name(file,"Eta"), cal.name("Eta"), ((int)binsEta_.size()-1), &binsEta_[0]);
  resEtaAll_= fs->make<TH1F>(res.name(file,"Eta"), res.name("Eta"), ((int)binsEta_.size()-1), &binsEta_[0]);

  for(int jdx=0; jdx<nJets_; ++jdx){
    std::vector<TH1F*> buffer;
    calPhiJet_.push_back(fs->make<TH1F>(cal.name(file,"Phi",jdx), cal.name("Phi",jdx), ((int)binsPhi_.size()-1), &binsPhi_[0]));
    resPhiJet_.push_back(fs->make<TH1F>(res.name(file,"Phi",jdx), res.name("Phi",jdx), ((int)binsPhi_.size()-1), &binsPhi_[0]));
    for(int idx=0; idx<((int)binsPhi_.size()-1); ++idx)
      buffer.push_back( fs->make<TH1F>(fit.name(file,"Phi",jdx,idx), fit.name("Phi",jdx,idx), 100, -0.5, 0.5) );
    PhiJet_.push_back(buffer); 
  }
  for(int idx = 0; idx < ((int)binsPhi_.size()-1); ++idx) {
    PhiAll_.push_back( fs->make<TH1F>(fit.name(file,"Phi",idx), fit.name("Phi",idx), 100, -0.5, 0.5) );
  }
  calPhiAll_= fs->make<TH1F>(cal.name(file,"Phi"), cal.name("Phi"), ((int)binsPhi_.size()-1), &binsPhi_[0]);
  resPhiAll_= fs->make<TH1F>(res.name(file,"Phi"), res.name("Phi"), ((int)binsPhi_.size()-1), &binsPhi_[0]);
}

/// write to file and free allocated space for FWLite
void 
JetResolution::write(const char* filename, const char* directory)
{
  /// save histograms to file
  TFile outFile( filename, "recreate" );
  outFile.mkdir( directory );
  outFile.cd( directory );

  /// basic kinematic
  calPtAll_->Write( );
  resPtAll_->Write( );
  calEtaAll_->Write( );
  resEtaAll_->Write( );
  calPhiAll_->Write( );
  resPhiAll_->Write( );

  outFile.Close();

  // free allocated space
  delete calPtAll_;
  delete resPtAll_;
  delete calEtaAll_;
  delete resEtaAll_;
  delete calPhiAll_;
  delete resPhiAll_;
}
