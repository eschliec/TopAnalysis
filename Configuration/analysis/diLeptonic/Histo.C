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
#include "plotterclass.h"
#include "HHStyle.h"
#include "basicFunctions.h"

// you can run me like this:
//    root -b -q 'Histo.C+g("preunfold", "vertmulti")'
// or root -b -q 'Histo.C+g("unfold", "toppt")'
// or root -b -q 'Histo.C+g("preunfold")' //all preunfolded
// or root -b -q 'Histo.C+g' //everything

void Histo(TString type = "", TString oneHistoToProcess = "") {
  bool doPreunfold = 1;
  bool doUnfold = 1;
  if (type == "preunfold") doUnfold = 0;
  if (type == "unfold") doPreunfold = 0;

  const double lumi = 12100;

  gROOT->SetBatch(kTRUE);

  std::vector<double> binCenters;
  std::vector<double> Xbins;

  string histolist = "HistoList";
  ifstream HistStream;
  HistStream.open(histolist.c_str());

  while(!HistStream.eof() && doUnfold){
  	
  	// Read HistoList-File
    TString name, specialComment, YAxis, XAxis;
    bool logX, logY, DYScale;
    int bins, rebin;
    double ymin, ymax, xmin, xmax;
    HistStream>>name>>specialComment>>YAxis>>XAxis>>rebin>>DYScale>>logX>>logY>>ymin>>ymax>>xmin>>xmax>>bins;
    
    // Avoid running over empty lines in 'HistoList'-File
    if ( name.CompareTo("") == 0) continue;
    
    
    //read in the remaining bins
    Xbins.clear();
    binCenters.clear();    
    for(int i = 0; i < bins+1; i++){
      double temp;
      HistStream>>temp; 
      Xbins.push_back(temp);
    }
    for(int i = 0; i < bins; i++){//only until bincenter code is finalized
      double temp;
      HistStream>>temp;
      binCenters.push_back(temp);
    }
    
    cout << "checking " << name << endl;
    if (! name.Contains(oneHistoToProcess, TString::kIgnoreCase)) continue;
    
    // Create Plotter 
    Plotter h_generalPlot;
    h_generalPlot.setLumi(lumi);
     
    /////////////////////////////////////////////////////
    /////////   UNFOLDING OPTIONS     ///////////////////
    /////////////////////////////////////////////////////
    
    // Unfolding Options
    bool doSVD = true;
    TString outpath = "";
    h_generalPlot.UnfoldingOptions(doSVD);
    h_generalPlot.SetOutpath("");
    
    /////////////////////////////////////////////////////
    /////////////////////////////////////////////////////
    /////////////////////////////////////////////////////
    
    h_generalPlot.setOptions(name,specialComment,YAxis,XAxis, rebin, DYScale, logX, logY, ymin, ymax, xmin, xmax, bins, Xbins, binCenters);
    h_generalPlot.DYScaleFactor();
    h_generalPlot.preunfolding();
    h_generalPlot.unfolding();
    h_generalPlot.PlotDiffXSec("emu");
    h_generalPlot.PlotDiffXSec("mumu");
    h_generalPlot.PlotDiffXSec("ee");
    h_generalPlot.PlotDiffXSec("combined");
  }
  
  cout << "Done with the unfolding\n";
  
  string controlhistolist = "HistoList_control";
  ifstream controlHistStream;
  controlHistStream.open(controlhistolist.c_str());

  while(!controlHistStream.eof() && doPreunfold){
  	// Read HistoList-File
    TString name, specialComment, YAxis, XAxis;
    bool logX, logY, DYScale;
    int bins, rebin;
    double ymin, ymax, xmin, xmax;
    controlHistStream>>name>>specialComment>>YAxis>>XAxis>>rebin>>DYScale>>logX>>logY>>ymin>>ymax>>xmin>>xmax>>bins;

    // Avoid running over empty lines in 'HistoList'-File
    if (name.CompareTo("") == 0) continue;
    
    Xbins.clear();
    binCenters.clear();
    
    for(int i = 0; i < bins+1; i++){
      double temp;
      controlHistStream>>temp; 
      Xbins.push_back(temp);
    }
    for(int i = 0; i < bins; i++){//only until bincenter code is finalized
      double temp;
      controlHistStream>>temp;
      binCenters.push_back(temp);
    }
    
    cout << "checking " << name << endl;
    if (! name.Contains(oneHistoToProcess, TString::kIgnoreCase)) continue;

    // Create Plotter 
    Plotter h_generalPlot;
    h_generalPlot.setLumi(lumi);
    /////////////////////////////////////////////////////
    /////////   UNFOLDING OPTIONS     ///////////////////
    /////////////////////////////////////////////////////
    
    // Unfolding Options
    bool doSVD = false; 
    TString outpath = "";
    h_generalPlot.UnfoldingOptions(doSVD);
    h_generalPlot.SetOutpath("");
    
    /////////////////////////////////////////////////////
    /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////// 
    

    h_generalPlot.setOptions(name,specialComment,YAxis,XAxis, rebin, DYScale, logX, logY, ymin, ymax, xmin, xmax, bins, Xbins, binCenters);
    h_generalPlot.DYScaleFactor();
    h_generalPlot.preunfolding();
  }
}

