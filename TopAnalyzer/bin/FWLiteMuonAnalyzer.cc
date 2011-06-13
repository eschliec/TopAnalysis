#include <cstdlib>

#include "TopAnalysis/TopUtils/bin/NiceStyle.cc"
#include "TopAnalysis/TopAnalyzer/interface/MuonId.h"
#include "TopAnalysis/TopAnalyzer/interface/MuonKinematic.h"
#include "TopAnalysis/TopAnalyzer/interface/MuonResolution.h"
#include "TopAnalysis/TopAnalyzer/interface/FWLiteAnalyzer.h"

using namespace std;


class FWLiteMuonAnalyzer : public FWLiteAnalyzer{
  
 public:
  
  explicit FWLiteMuonAnalyzer();
  virtual ~FWLiteMuonAnalyzer();

  /// FWLiteAnalyzer begin
  virtual void begin(std::string&);
  /// FWLiteAnalyzer analyze 
  virtual void analyze(fwlite::Event&);
  /// FWLiteAnalyzer end
  virtual void end();
  
 protected:

  /// read Hist block from cfg file 
  void configBlockHist(ConfigFile&);

 private:
  
  bool doMuid_;
  bool doMukin_;
  bool doMures_;
  double matchDR_;
  std::vector<double> binsPt_;

  MuonId* muid_; 
  MuonKinematic* mukin_; 
  MuonResolution* mures_; 
};


FWLiteMuonAnalyzer::FWLiteMuonAnalyzer()
{
  if(doMuid_ ) muid_ = new MuonId();
  if(doMukin_) mukin_= new MuonKinematic();
  if(doMures_) mures_= new MuonResolution();
}

FWLiteMuonAnalyzer::~FWLiteMuonAnalyzer()
{
  if(doMuid_ ) {delete muid_ ; muid_ =0;};
  if(doMukin_) {delete mukin_; mukin_=0;};
  if(doMures_) {delete mures_; mures_=0;};
}

void
FWLiteMuonAnalyzer::begin(std::string& cfg)
{
  // read necessary parameters from cfg file
  ConfigFile file( cfg, "=", "#" );  
  configBlockIO  ( file );
  configBlockHist( file );

  // do histogram booking
  if(doMuid_ ){ 
    muid_ ->book(); 
  }
  if(doMures_){ 
    mukin_->book(); 
  }
  if(doMures_){
    mures_->setMatchDR(matchDR_); // these have to be called explicitely as
    mures_->setBinsPt (binsPt_ ); // passing them via the constructor fails
    mures_->book();
  }
}

void FWLiteMuonAnalyzer::analyze(fwlite::Event& event)
{
  // the following lines are commented
  // as weights are not in the event, yet
  //fwlite::Handle<double> weight;
  //weight.getByLabel(event, "csa07:weight");
  
  fwlite::Handle<std::vector<pat::Jet> > jets;
  jets.getByLabel(event, "selectedLayer1Jets");
  
  fwlite::Handle<std::vector<pat::Muon> > muons;
  muons.getByLabel(event, "selectedLayer1Muons");

  double weight=1.;  
  if(doMuid_ ){
    muid_ ->fill( *muons, weight); 
  }
  if(doMukin_){ 
    mukin_->fill( *jets, *muons, weight); 
  }
  if(doMures_){
    mures_->fill( *muons, weight); 
  }
}

void
FWLiteMuonAnalyzer::end()
{
  // normalize histograms if desired
  if(doMuid_ ) muid_ ->norm(); 
  if(doMukin_) mukin_->norm(); 
  if(doMures_) mures_->norm(); 
    
  // write output to file
  TFile file( outputFile_.c_str(), "recreate" );
  file.mkdir( outputDirectory_.c_str() );
  if(doMuid_ ) muid_ ->write( file, outputDirectory_.c_str() ); 
  if(doMukin_) mukin_->write( file, outputDirectory_.c_str() ); 
  if(doMures_) mures_->write( file, outputDirectory_.c_str() ); 
  file.Close();
}

void 
FWLiteMuonAnalyzer::configBlockHist(ConfigFile& cfg)
{
  //-----------------------------------------------
  // read configurables defined in FWAnalyzer from 
  // config file. Throw human readable exception 
  // when misspellings occure
  //-----------------------------------------------
  try {
    //-----------------------------------------------
    // input/output files
    //-----------------------------------------------
    doMuid_ = cfg.read<bool>( "doMuid"  );
    doMukin_= cfg.read<bool>( "doMukin" );
    doMures_= cfg.read<bool>( "doMures" );
    matchDR_= cfg.read<double>( "matchDR" );
    readVector( cfg.read<std::string>( "binsPt" ), binsPt_ ); 
  }
  catch(...) {
    std::cerr << "ERROR during reading of config file" << std::endl;
    std::cerr << "      misspelled variables in cfg ?" << std::endl;
    std::cerr << "      [--called in configBlockHist--]" << std::endl;
    std::exit(1);
  }
}



int main(int argc, char* argv[]) 
{
  if( argc<2 ){
    // ------------------------------------------------- 
    std::cerr << "ERROR:: " 
	      << "Wrong number of arguments! Please specify:" << std::endl
	      << "        * configfile" << std::endl;
    // -------------------------------------------------  
    return -1;
  }
  std::string cfgFile(argv[1]); 

  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  
  // set nice style for histograms
  setNiceStyle();

  // define worker class
  FWLiteMuonAnalyzer muana; 
  // configure
  muana.begin(cfgFile);
  // keep this! It is the evt loop
  // analyze is called internally
  muana.event();
  // clean up 
  muana.end();

  return 0;
}
