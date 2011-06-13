#include "TopAnalysis/TopAnalyzer/interface/ElecResolution.h"


/// constructor for FWLite analyzer
ElecResolution::ElecResolution(std::vector<double> binsPt, double matchDR):
  binsPt_(binsPt), matchDR_(matchDR)  
{
}

/// constructor for full FW analyzer
ElecResolution::ElecResolution(const edm::ParameterSet& cfg):
  binsPt_ ( cfg.getParameter<std::vector<double> >( "binsPt" ) ),
  binsEta_( cfg.getParameter<std::vector<double> >( "binsEta") ),
  binsPhi_( cfg.getParameter<std::vector<double> >( "binsPhi") ),
  matchDR_( cfg.getParameter<double>( "matchDR" ) )
{
}

/// fill interface for full FW analyzer
void
ElecResolution::fill(const edm::Event& evt, const std::vector<pat::Electron>& elecs, const double& weight=1.)
{
  if(elecs.begin()==elecs.end()) return;
  fill(elecs, weight);
}

/// fill interface for FWLite analyzer
void
ElecResolution::fill(const std::vector<pat::Electron>& elecs, const double& weight=1.)
{
  std::vector<pat::Electron>::const_iterator elec=elecs.begin();
  if(elec!=elecs.end()){
    if( elec->genLepton()!=0 ){ 
      double dR=deltaR( elec->eta(), elec->phi(), elec->genLepton()->eta(), elec->genLepton()->phi() ); 
      if(dR<matchDR_){
	double dPt=(elec->pt()-elec->genLepton()->pt())/(elec->genLepton())->pt();
        double dPhi=deltaPhi( elec->phi(), elec->genLepton()->phi() );
        double dEta=(elec->eta()-elec->genLepton()->eta());
	for(unsigned int idx=0; idx<(binsPt_.size()-1); ++idx) {
	  if( (binsPt_[idx]<elec->genLepton()->pt()) && (elec->genLepton()->pt()<binsPt_[idx+1]) ){
	    relPt_[idx]->Fill(dPt, weight);
	  }
	}
	for(unsigned int idx=0; idx<(binsEta_.size()-1); ++idx) {
	  if( (binsEta_[idx]<elec->genLepton()->eta()) && (elec->genLepton()->eta()<binsEta_[idx+1]) ){
	    Eta_[idx]->Fill(dEta, weight);
	  }
	}
	for(unsigned int idx=0; idx<(binsPhi_.size()-1); ++idx) {
	  if( (binsPhi_[idx]<elec->genLepton()->phi()) && (elec->genLepton()->phi()<binsPhi_[idx+1]) ){
	    Phi_[idx]->Fill(dPhi, weight);
	  }
	}
      }
    }
  }
}

/// book for FWLite
void 
ElecResolution::book()
{
  NameScheme fit("fit"), cal("cal"), res("res");
  for(int idx = 0; idx < ((int)binsPt_.size()-1); ++idx) {
    relPt_.push_back( new TH1F(fit.name("relPt",idx), fit.name("relPt",idx), 100, -0.5, 0.5) );
  }
  calPt_= new TH1F(cal.name("relPt"), cal.name("relPt"), ((int)binsPt_.size()-1), &binsPt_[0]);
  resPt_= new TH1F(res.name("relPt"), res.name("relPt"), ((int)binsPt_.size()-1), &binsPt_[0]);
  for(int idx = 0; idx < ((int)binsEta_.size()-1); ++idx) {
    Eta_.push_back( new TH1F(fit.name("Eta",idx), fit.name("Eta",idx), 100, -0.005, 0.005) );
  }
  calEta_= new TH1F(cal.name("Eta"), cal.name("Eta"), ((int)binsEta_.size()-1), &binsEta_[0]);
  resEta_= new TH1F(res.name("Eta"), res.name("Eta"), ((int)binsEta_.size()-1), &binsEta_[0]);
  for(int idx = 0; idx < ((int)binsPhi_.size()-1); ++idx) {
    Phi_.push_back( new TH1F(fit.name("Phi",idx), fit.name("Phi",idx), 100, -0.005, 0.005) );
  }
  calPhi_= new TH1F(cal.name("Phi"), cal.name("Phi"), ((int)binsPhi_.size()-1), &binsPhi_[0]);
  resPhi_= new TH1F(res.name("Phi"), res.name("Phi"), ((int)binsPhi_.size()-1), &binsPhi_[0]);
}

/// book for full FW
void 
ElecResolution::book(edm::Service<TFileService>& fs)
{
  NameScheme fit("fit"), cal("cal"), res("res");
  for(int idx = 0; idx < ((int)binsPt_.size()-1); ++idx) {
    relPt_.push_back( fs->make<TH1F>(fit.name("relPt",idx), fit.name("relPt",idx), 100, -0.5, 0.5) );
  }
  calPt_= fs->make<TH1F>(cal.name("relPt"), cal.name("relPt"), ((int)binsPt_.size()-1), &binsPt_[0]);
  resPt_= fs->make<TH1F>(res.name("relPt"), res.name("relPt"), ((int)binsPt_.size()-1), &binsPt_[0]);
  for(int idx = 0; idx < ((int)binsEta_.size()-1); ++idx) {
    Eta_.push_back( fs->make<TH1F>(fit.name("Eta",idx), fit.name("Eta",idx), 100, -0.005, 0.005) );
  }
  calEta_= fs->make<TH1F>(cal.name("Eta"), cal.name("Eta"), ((int)binsEta_.size()-1), &binsEta_[0]);
  resEta_= fs->make<TH1F>(res.name("Eta"), res.name("Eta"), ((int)binsEta_.size()-1), &binsEta_[0]);
  for(int idx = 0; idx < ((int)binsPhi_.size()-1); ++idx) {
    Phi_.push_back( fs->make<TH1F>(fit.name("Phi",idx), fit.name("Phi",idx), 100, -0.005, 0.005) );
  }
  calPhi_= fs->make<TH1F>(cal.name("Phi"), cal.name("Phi"), ((int)binsPhi_.size()-1), &binsPhi_[0]);
  resPhi_= fs->make<TH1F>(res.name("Phi"), res.name("Phi"), ((int)binsPhi_.size()-1), &binsPhi_[0]);
}

/// book for full FW with output stream
void 
ElecResolution::book(edm::Service<TFileService>& fs, ofstream& file)
{
  NameScheme fit("fit"), cal("cal"), res("res");
  for(int idx = 0; idx < ((int)binsPt_.size()-1); ++idx) {
    relPt_.push_back( fs->make<TH1F>(fit.name(file,"relPt",idx), fit.name("relPt",idx), 100, -0.5, 0.5) );
  }
  calPt_= fs->make<TH1F>(cal.name(file,"relPt"), cal.name("relPt"), ((int)binsPt_.size()-1), &binsPt_[0]);
  resPt_= fs->make<TH1F>(res.name(file,"relPt"), res.name("relPt"), ((int)binsPt_.size()-1), &binsPt_[0]);
  for(int idx = 0; idx < ((int)binsEta_.size()-1); ++idx) {
    Eta_.push_back( fs->make<TH1F>(fit.name(file,"Eta",idx), fit.name("Eta",idx), 100, -0.005, 0.005) );
  }
  calEta_= fs->make<TH1F>(cal.name(file,"Eta"), cal.name("Eta"), ((int)binsEta_.size()-1), &binsEta_[0]);
  resEta_= fs->make<TH1F>(res.name(file,"Eta"), res.name("Eta"), ((int)binsEta_.size()-1), &binsEta_[0]);
  for(int idx = 0; idx < ((int)binsPhi_.size()-1); ++idx) {
    Phi_.push_back( fs->make<TH1F>(fit.name(file,"Phi",idx), fit.name("Phi",idx), 100, -0.005, 0.005) );
  }
  calPhi_= fs->make<TH1F>(cal.name(file,"Phi"), cal.name("Phi"), ((int)binsPhi_.size()-1), &binsPhi_[0]);
  resPhi_= fs->make<TH1F>(res.name(file,"Phi"), res.name("Phi"), ((int)binsPhi_.size()-1), &binsPhi_[0]);
}

/// write to file and free allocated space for FWLite
void 
ElecResolution::write(const char* filename, const char* directory)
{
  /// save histograms to file
  TFile outFile( filename, "recreate" );
  outFile.mkdir( directory );
  outFile.cd( directory );

  /// basic kinematic
  calPt_->Write( );
  resPt_->Write( );
  calEta_->Write( );
  resEta_->Write( );
  calPhi_->Write( );
  resPhi_->Write( );

  outFile.Close();

  // free allocated space
  delete calPt_;
  delete resPt_;
  delete calEta_;
  delete resEta_;
  delete calPhi_;
  delete resPhi_;
}
