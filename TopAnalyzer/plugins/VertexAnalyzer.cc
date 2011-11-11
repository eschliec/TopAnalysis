#include "TopAnalysis/TopAnalyzer/plugins/VertexAnalyzer.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "TopAnalysis/TopAnalyzer/interface/DileptonEventWeight.h"

VertexAnalyzer::VertexAnalyzer(const ParameterSet& cfg):
  vertices_   (cfg.getParameter<InputTag>     ( "vertices"    )),
  leptons_    (cfg.getParameter<InputTag>     ( "leptons"     )),
  puWeight_   (cfg.getParameter<edm::InputTag>( "weightPU"    )),
  lepSfWeight_(cfg.getParameter<edm::InputTag>( "weightLepSF" )),  
  ndof_       (cfg.getParameter<unsigned int> ( "ndof"        )),
  rho_        (cfg.getParameter<double>	      ( "rho"         )),
  z_          (cfg.getParameter<double>       ( "z"           ))  
{
}

VertexAnalyzer::~VertexAnalyzer()
{
}

void
VertexAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;
  if( !fs ){
    throw edm::Exception( edm::errors::Configuration,
                          "TFile Service is not registered in cfg file" );
  }

  multi_= fs->make<TH1D>( "multi", "N of Vertices", 30,-0.5,29.5);
  multi_->GetXaxis()->SetTitle("N_{vrtcs}");
  multi_->GetYaxis()->SetTitle("N");

  goodMulti_= fs->make<TH1D>( "goodMulti", "N of good Vertices", 30,-0.5,29.5);
  goodMulti_->GetXaxis()->SetTitle("N_{vrtcs}");
  goodMulti_->GetYaxis()->SetTitle("N");


  posX_= fs->make<TH1D>( "posX", "Vertex x-Position", 400,-2.,2.);
  posX_->GetXaxis()->SetTitle("x [cm]");
  posX_->GetYaxis()->SetTitle("N / 0.01cm");

  posY_= fs->make<TH1D>( "posY", "Vertex y-Position", 400,-2.,2.);
  posY_->GetXaxis()->SetTitle("y [cm]");
  posY_->GetYaxis()->SetTitle("N / 0.01cm");

  posZ_= fs->make<TH1D>( "posZ", "Vertex z-Position", 1000,-50.,50.);
  posZ_->GetXaxis()->SetTitle("z [cm]");
  posZ_->GetYaxis()->SetTitle("N / 0.1cm");

  posRho_= fs->make<TH1D>( "posRho", "Vertex #rho-Position", 500,-2.5,2.5);
  posRho_->GetXaxis()->SetTitle("#rho [cm]");
  posRho_->GetYaxis()->SetTitle("N / 0.01cm");


  posXerr_= fs->make<TH1D>( "posXerr", "Vertex x-Position Error", 100,0.,1.);
  posXerr_->GetXaxis()->SetTitle("#delta x [cm]");
  posXerr_->GetYaxis()->SetTitle("N / 0.01cm");

  posYerr_= fs->make<TH1D>( "posYerr", "Vertex y-Position Error", 100,0.,1.);
  posYerr_->GetXaxis()->SetTitle("#delta y [cm]");
  posYerr_->GetYaxis()->SetTitle("N / 0.01cm");

  posZerr_= fs->make<TH1D>( "posZerr", "Vertex z-Position Error", 100,0.,1.);
  posZerr_->GetXaxis()->SetTitle("#delta z [cm]");
  posZerr_->GetYaxis()->SetTitle("N / 0.01cm");


  nTracks_= fs->make<TH1D>( "nTracks", "Number of Tracks at Vertex", 200,0,200);
  nTracks_->GetXaxis()->SetTitle("N_{tracks}");
  nTracks_->GetYaxis()->SetTitle("N");

  nDof_= fs->make<TH1D>( "nDof", "N_{dof} at Vertex", 200,0,200);
  nDof_->GetXaxis()->SetTitle("N_{dof}");
  nDof_->GetYaxis()->SetTitle("N");

  chi2_= fs->make<TH1D>( "chi2", "#chi^2 of Vertices", 200,0,200);
  chi2_->GetXaxis()->SetTitle("#chi^2");
  chi2_->GetYaxis()->SetTitle("N");

  nchi2_= fs->make<TH1D>( "nchi2", "#chi^2/N_{dof} of Vertices", 200,0,20);
  nchi2_->GetXaxis()->SetTitle("#chi^2/N_{dof}");
  nchi2_->GetYaxis()->SetTitle("N");

  dzLep_= fs->make<TH1D>( "dzLep", "z-Distance Lepton,Vrtx", 200,0,20);
  dzLep_->GetXaxis()->SetTitle("#Delta z [cm]");
  dzLep_->GetYaxis()->SetTitle("N");

  position3D_= fs->make<TH3D>( "position3D", "Position of good Vertices", 1000,-50.,50.,100,-5.,5.,100,-5.,5.);
  position3D_->GetXaxis()->SetTitle("z [cm]"); // this is NO typo!
  position3D_->GetYaxis()->SetTitle("x [cm]");
  position3D_->GetZaxis()->SetTitle("y [cm]");
}

void
VertexAnalyzer::analyze(const Event& evt, const EventSetup&)
{
  double weight = getDileptonEventWeight(evt, puWeight_, lepSfWeight_);
  edm::Handle<std::vector<reco::Vertex> > vertices;
  evt.getByLabel(vertices_, vertices);

  multi_->Fill(vertices->size(), weight);

  int i=0;
  for(VertexCollection::const_iterator vrtx = vertices->begin(); vrtx!= vertices->end(); ++vrtx) {

    bool isFake = vrtx->isFake();
    unsigned int ndof    = vrtx->ndof();
    unsigned int ntracks = vrtx->tracksSize();
    double rho  = vrtx->position().rho();
    double z    = vrtx->position().z();
    double chi2 = vrtx->chi2();
    double nchi2= chi2/ndof;

    if(!isFake && ndof>ndof_ && abs(z)<z_ && rho<rho_){
      i++;

      posX_  ->Fill(vrtx->x(), weight);
      posY_  ->Fill(vrtx->y(), weight);
      posZ_  ->Fill(z, weight);
      posRho_->Fill(rho, weight);

      posXerr_->Fill(vrtx->xError(), weight);
      posYerr_->Fill(vrtx->yError(), weight);
      posZerr_->Fill(vrtx->zError(), weight);

      nTracks_->Fill(ntracks, weight);
      nDof_   ->Fill(ndof, weight);
      chi2_   ->Fill(chi2, weight);
      nchi2_  ->Fill(nchi2, weight);

      position3D_->Fill(vrtx->z(),vrtx->x(),vrtx->y(), weight); // this is NO typo!

      edm::Handle<edm::View<reco::Candidate> > leptons;
      evt.getByLabel(leptons_, leptons);

      for(edm::View<reco::Candidate>::const_iterator lepton = leptons->begin(); lepton!= leptons->end(); ++lepton) {
        dzLep_->Fill(abs(lepton->vz()-z), weight);
      }
    }
  }
  goodMulti_->Fill(i, weight);
}

void
VertexAnalyzer::endJob()
{
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( VertexAnalyzer );
