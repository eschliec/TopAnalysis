#include "ControlPlots.h"


void ControlPlots::MakePlots(){

  gROOT->SetBatch(kTRUE);

  std::vector<double> Xbins;
  std::vector<double> binCenters;

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
    TString outpath = "";
    
    h_generalPlot.SetOutpath("");
    
    /////////////////////////////////////////////////////
    /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////// 
    
    h_generalPlot.SetDataLumi(5100);
    h_generalPlot.setOptions(name,specialComment,YAxis,XAxis, rebin, DYScale, logX, logY, ymin, ymax, xmin, xmax, bins, Xbins, binCenters);
//     h_generalPlot.DYScaleFactor();
//     h_generalPlot.setDataSet("mumu");
//     h_generalPlot.fillHisto();
//     h_generalPlot.write();
    h_generalPlot.DYScaleFactor();
    h_generalPlot.setDataSet("emu");
    h_generalPlot.fillHisto();
    h_generalPlot.write();
//     h_generalPlot.DYScaleFactor();
//     h_generalPlot.setDataSet("ee");
//     h_generalPlot.fillHisto();
//     h_generalPlot.write();
//     h_generalPlot.DYScaleFactor();
//     h_generalPlot.setDataSet("combined");
//     h_generalPlot.fillHisto();
//     h_generalPlot.write();
  }
  return;
}

ControlPlots::ControlPlots(){
  MakePlots();
}

ControlPlots::~ControlPlots(){
}

