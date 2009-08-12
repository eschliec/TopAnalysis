#ifndef SingleObject_h
#define SingleObject_h

#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

//
// interface class for a common object analyzer, which should 
// be equally usable for full framework or framework lite 
//

template <typename Collection> 
class SingleObject{

 public:
  /// default constructor for fw lite
  explicit SingleObject(){};
  /// default constructor for full fw
  explicit SingleObject(const edm::ParameterSet& configFile){};
  /// default destructor
  virtual ~SingleObject(){};
  /// write histograms to file for fwlite
  void write(TFile& file, const char* directory);

  /// histogramm booking for fwlite 
  virtual void book() = 0;
  /// histogramm booking for full fw
  virtual void book(edm::Service<TFileService>& fileService) = 0;
  /// histogram filling for fwlite and for full fw
  virtual void fill(const Collection& inputCollection, const double& weight=1.) = 0;
  /// everything which needs to be done after the event loop
  virtual void process() = 0;

 protected:
  /// histogram container
  std::map<std::string, TH1*> hists_;
};


/// writing histograms to file in fwlite
template <typename Collection> 
void SingleObject<Collection>::write(TFile& file, const char* directory)
{
  file.cd( directory );
  for(std::map<std::string, TH1*>::const_iterator hist = hists_.begin(); hist !=hists_.end(); ++hist){
    hist->second->Write( );
  }
}

#endif
