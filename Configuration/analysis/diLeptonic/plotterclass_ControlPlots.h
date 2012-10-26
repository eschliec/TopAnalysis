#ifndef plotterclass_ControlPlots_h
#define plotterclass_ControlPlots_h


#include <vector>
#include <map>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

#include "TCanvas.h"
#include "THStack.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TExec.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include "basicFunctions.h"
#include "HHStyle.h"

class Plotter {

    public:
        Plotter();
        Plotter(TString, TString, TString, double, double);
        ~Plotter();
        void setOptions(TString, TString, TString, TString, int, bool, bool, bool, double, double, double, double, int, std::vector<double>, std::vector<double>);
        void setDataSet(std::vector<TString>, std::vector<double>, std::vector<TString>, std::vector<int>, TString);
        void setDataSet(TString);
        void fillSystHisto();
        void fillHisto();
        void setStyle(TH1D&, unsigned int);
        void write();
        void DYScaleFactor();
        // DAVID
        void SetOutpath(TString path); 
        // IVAN
        TLegend* ControlLegend(int HistsSize, TH1* drawhists[], std::vector<TString> legends, TLegend *leg);
        TLegend* ControlLegend(int HistsSize, TH1D* drawhists[], std::vector<TString> Legends, TLegend *leg);
        void DrawLabel(TString text, const double x1, const double y1, const double x2, const double y2, int centering, double textSize);
        void ApplyFlatWeights(TH1* varhists,   const double weight);
        void ApplyFlatWeights(TH1* varhists[], const double weight);
        double SampleXSection(TString filename);
        double CalcLumiWeight(TString WhichSample);
        void SetDataLumi(double Lumi);
        void DrawDecayChLabel(TString decaychannel = "", double textSize = 0.04);

    private:
        TString name;
        TString specialComment;
        TString DYEntry;
        TString YAxis;
        TString XAxis;
        TString channel;
        int bins, datafiles, rebin;
        int channelType; //0=ee 1==mumu 2==emu 3==combined  
        int signalHist;
        double rangemin, rangemax, ymin, ymax;
        std::vector<TFile> files;
        std::vector<TString> dataset;
        std::vector<TString> legends;
        std::vector<TString> channelLabel;
        std::vector<int> colors;
        std::vector<double> scales;
        std::vector<double> XAxisbins, XAxisbinCenters, DYScale;
        std::vector<TH1D> hists;
        bool initialized, logX, logY, doDYScale;

        
        // DAVID
        TString outpath;
        TString outpathPlots;
        TString subfolderChannel;
        TString subfolderSpecial;
        //IVAN
        double lumi;
};

#endif
