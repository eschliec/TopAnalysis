#include "Histo.h"

void Histo::MakePlots(){

  gROOT->SetBatch(kTRUE);

  std::vector<double> binCenters;
  std::vector<double> Xbins;

  string histolist = "HistoList";
  ifstream HistStream;
  HistStream.open(histolist.c_str());

  while(!HistStream.eof()){
  	
  	// Read HistoList-File
    TString name, specialComment, YAxis, XAxis;
    bool logX, logY, DYScale;
    int bins, rebin;
    double ymin, ymax, xmin, xmax;
    HistStream>>name>>specialComment>>YAxis>>XAxis>>rebin>>DYScale>>logX>>logY>>ymin>>ymax>>xmin>>xmax>>bins;
    
    // Avoid running over empty lines in 'HistoList'-File
    if ( name.CompareTo("") == 0 ) continue;
    
    // Create Plotter 
    Plotter h_generalPlot;
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
  
  string controlhistolist = "HistoList_control";
  ifstream controlHistStream;
  controlHistStream.open(controlhistolist.c_str());

  while(!controlHistStream.eof()){
  	// Read HistoList-File
    TString name, specialComment, YAxis, XAxis;
    bool logX, logY, DYScale;
    int bins, rebin;
    double ymin, ymax, xmin, xmax;
    controlHistStream>>name>>specialComment>>YAxis>>XAxis>>rebin>>DYScale>>logX>>logY>>ymin>>ymax>>xmin>>xmax>>bins;

    // Avoid running over empty lines in 'HistoList'-File
    if ( name.CompareTo("") == 0 ) continue;
    // Create Plotter 
    Plotter h_generalPlot;
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
  return;
}

Histo::Histo(){
  MakePlots();
}

Histo::~Histo(){
}

