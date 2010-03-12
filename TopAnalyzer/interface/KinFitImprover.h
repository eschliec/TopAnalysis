#ifndef KinFitImprover_h
#define KinFitImprover_h

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "TopAnalysis/TopAnalyzer/interface/DoubleObject.h"
#include "AnalysisDataFormats/TopObjects/interface/TtFullHadronicEvent.h"

/**
   \class   KinFitImprover KinFitImprover.h "TopAnalysis/TopAnalyzer/interface/KinFitImprover.h"

   \brief   Derived class to analyze and improve the kinematic fit hypothesis of fullhadronic ttbar events

   The structure keeps control histograms to analyze the relation between the fitted four vectors 
   of the kinematic fit hypothesis of fullhadronic ttbar events and the original input vectors 
   of the selected pat::Jets. These histograms are to be filled from the TtFullHadronicEvent and 
   from edm::View<pat::Jet>'s. The class is derived from the DoubleObject<CollectionA, CollectionB> 
   interface, which makes it usable in fwfull or fwlite.  This class is expected to be run with 
   the same jet collection from which the kinematic fit hypothesis was derived, otherwise it will 
   lead to wrong results!
*/

class KinFitImprover : public DoubleObject<TtFullHadronicEvent, const edm::View<pat::Jet> > {

 public:
  /// default constructor for fw lite
  explicit KinFitImprover();
  /// default constructor for fw full
  explicit KinFitImprover(const edm::ParameterSet& configFile);
  /// default destructor
  ~KinFitImprover(){};

  /**
     The following functions have to be implemented for any class
     derived from DoubleObject<CollectionA, CollectionB>
  **/
  /// histogramm booking for fwlite 
  void book();
  /// histogramm booking for fw
  void book(edm::Service<TFileService>& fileService);
  /// histogram filling interface for reconstruction level for access with fw or fwlite
  void fill(const TtFullHadronicEvent& tops, const edm::View<pat::Jet>& jets, const double& weight=1.);
  /// everything which needs to be done after the event loop
  void process(){};

 private:
  /// function to find types of jet-combinations in KinFits
  int comboType();
  /// vector to store the jet indices
  std::vector<int> jetIndex;
};

#endif
