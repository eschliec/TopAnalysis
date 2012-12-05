#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>

#include <vector>
#include <iostream>
#include <set>

#include "plotterclass.h"
#include "HistoListReader.h"

// you can run me like this:
//    root -b -q 'Histo.C+g("preunfold", "vertmulti", "JERUP", "emu")'
//    root -b -q 'Histo.C+g("preunfold", "vertmulti", "JERUP")' //all channels
// or root -b -q 'Histo.C+g("preunfold", "vertmulti")' //all systematics and all channels
// or root -b -q 'Histo.C+g("unfold", "toppt")'
// or root -b -q 'Histo.C+g("preunfold")' //all preunfolded
// or root -b -q 'Histo.C+g' //everything

using namespace std;

std::set<TString> SetOfValidSystematics(){
    
    //Dummy function to create a 'set' containing the list of all systematics
    //Include also '' to be able to run in all the systematics
    static TString syst_array[]={"Nominal", "JERUP", "JERDOWN", "JESUP", "JESDOWN", "PU_UP", "PU_DOWN", "TRIG_UP", "TRIG_DOWN", "BTAG_UP", "BTAG_DOWN", "BTAG_PT_UP", "BTAG_PT_DOWN",
                          "BTAG_ETA_UP", "BTAG_ETA_DOWN", "MASSUP", "MASSDOWN", "MATCHUP", "MASSDOWN", "SCALEUP", "SCALEDOWN", "POWHEG", "MCATNLO", "SPINCORR", ""};
    std::set<TString> SetOfSystematics (syst_array, syst_array+sizeof(syst_array)/sizeof(syst_array[0]));
    return SetOfSystematics;
}


void Histo(TString type = "", TString oneHistoToProcess = "", TString systematic="", TString channel="") {
    const double lumi = 12100;

    //Take the list of systematica variations from 'SetOfValidSystematics()' and check if the systematic you want to run exists. If doesn't return
    set<TString> ListOfSysts = SetOfValidSystematics();

    //Check if the channel in which the code will run is valid: ee, emu, mumu, combined or '' (all the channels)
    if(channel != "ee" && channel != "emu" && channel != "mumu" && channel != "combined" && channel != ""){
        cout<<"\n\nWARNING (in Histo.C)!! You are using an unsupported channel type."<<endl;
        cout<<"Please use: ee, emu, mumu, combined or '' (all channels)"<<endl;
        cout<<"Returning"<<endl;
        return;
    }

    bool doPreunfold = 1;
    bool doUnfold = 1;
    if (type == "preunfold") doUnfold = 0;
    if (type == "unfold") doPreunfold = 0;

    if (doUnfold) {
        HistoListReader histoList("HistoList");
        if (histoList.IsZombie()) exit(11);
        for (map< TString, PlotProperties >::iterator it = histoList.begin(); it != histoList.end(); ++it) {
            const PlotProperties& p = it->second;
            cout << "checking " << p.name << endl;
            if (! p.name.Contains(oneHistoToProcess, TString::kIgnoreCase)) continue;

            // Create Plotter 
            Plotter h_generalPlot;
            h_generalPlot.setLumi(lumi);
            h_generalPlot.ListOfSystematics(ListOfSysts);

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

            h_generalPlot.setOptions(p.name,p.specialComment,p.ytitle,p.xtitle, 
                                     p.rebin, p.do_dyscale, p.logX, p.logY, 
                                     p.ymin, p.ymax, p.xmin, p.xmax, p.bins, p.xbinbounds, p.bincenters);
            h_generalPlot.DYScaleFactor();
            h_generalPlot.preunfolding(channel, systematic);
            h_generalPlot.unfolding();
            h_generalPlot.PlotDiffXSec("emu");
            h_generalPlot.PlotDiffXSec("mumu");
            h_generalPlot.PlotDiffXSec("ee");
            h_generalPlot.PlotDiffXSec("combined");
            
        }
    }
    cout << "Done with the unfolding\n";


    if (doPreunfold) {
        HistoListReader histoList("HistoList_control");
        if (histoList.IsZombie()) exit(12);
        for (map< TString, PlotProperties >::iterator it = histoList.begin(); it != histoList.end(); ++it) {
            const PlotProperties& p = it->second;
            cout << "checking " << p.name << endl;
            if (! p.name.Contains(oneHistoToProcess, TString::kIgnoreCase)) continue;

            // Create Plotter 
            Plotter h_generalPlot;
            h_generalPlot.setLumi(lumi);
            h_generalPlot.ListOfSystematics(ListOfSysts);
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
            
            h_generalPlot.setOptions(p.name,p.specialComment,p.ytitle,p.xtitle, 
                                     p.rebin, p.do_dyscale, p.logX, p.logY, 
                                     p.ymin, p.ymax, p.xmin, p.xmax, p.bins, p.xbinbounds, p.bincenters);
            h_generalPlot.DYScaleFactor();
            h_generalPlot.preunfolding(channel, systematic);
        }
    }
}

int main(int argc, const char * const argv[]) {
    TString type = argc > 1 ? argv[1] : "";
    TString oneHistoToProcess = argc > 2 ? argv[2] : "";
    TString systematic = argc > 3 ? argv[3] : "";
    TString channel = argc > 4 ? argv[4] : "";
    Histo(type, oneHistoToProcess, systematic, channel);
}
