#ifndef FWAnalyzer_h
#define FWAnalyzer_h

#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include "DataFormats/FWLite/interface/Handle.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "TopAnalysis/TopUtils/interface/ConfigFile.h"
#include "TopAnalysis/TopUtils/interface/RootSystem.h"
#include "TopAnalysis/TopUtils/interface/RootHistograms.h"
#include "TopAnalysis/TopUtils/interface/RootPostScript.h"

#define DEBUG

class FWLiteAnalyzer{

 public:

  FWLiteAnalyzer(){};
  virtual ~FWLiteAnalyzer(){};

  /// event loop for FWLiteAnalyzer
  virtual void event();

  /// ------------------------------------------------------
  /// implement the following functions 
  /// in your derived class
  /// ------------------------------------------------------

  /// FWLiteAnalyzer begin
  virtual void begin(std::string&)=0;
  /// FWLiteAnalyzer analyze
  virtual void analyze(fwlite::Event&)=0;
  /// FWLiteAnalyzer end
  virtual void end()=0;


 protected:

  /// read IO block from cfg file
  void configBlockIO(ConfigFile&);
  /// return input files
  std::vector<std::string> inputFiles() { return inputFiles_; };
  /// break event loop after max events
  bool stop(int evt) { return (maxEvt_>0 ? evt>maxEvt_ : false); };

 protected:

  int maxEvt_;
  std::string outputFile_;
  std::string outputDirectory_;
  std::vector<std::string> inputFiles_;
};

#endif
