#ifndef ControlPlots_h
#define ControlPlots_h

#include "THStack.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH1F.h"
#include <vector>
#include <iostream>
#include "TCanvas.h"
#include "TLegend.h"
#include <sstream>
#include "plotterclass_ControlPlots.h"   //modify if new controlplot ploterclass is defined
#include "plotterclass_ControlPlots.C"
#include "HHStyle.h"
#include "basicFunctions.h"

class ControlPlots {

 public:
  ControlPlots();
  ~ControlPlots();
  void MakePlots();

 private:

  //  Plotter h_generalPlot;
};

#endif
    
