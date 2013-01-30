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
#include "CommandLineParameters.hh"

using namespace std;

const std::vector<const char*> VectorOfValidSystematics 
    {"Nominal", 
    "JER_UP", "JER_DOWN", "JES_UP", "JES_DOWN", "PU_UP", "PU_DOWN", "TRIG_UP", "TRIG_DOWN", "LEPT_UP", "LEPT_DOWN",
    "BTAG_UP", "BTAG_DOWN", "BTAG_LJET_UP", "BTAG_LJET_DOWN",
    "BTAG_PT_UP", "BTAG_PT_DOWN", "BTAG_ETA_UP", "BTAG_ETA_DOWN",
    "BTAG_LJET_PT_UP", "BTAG_LJET_PT_DOWN", "BTAG_LJET_ETA_UP", "BTAG_LJET_ETA_DOWN",
    "BTAG_BEFF_UP", "BTAG_BEFF_DOWN", "BTAG_CEFF_UP", "BTAG_CEFF_DOWN", "BTAG_LEFF_UP", "BTAG_LEFF_DOWN",
    "MASS_UP", "MASS_DOWN", "MATCH_UP", "MASS_DOWN",
    "SCALE_UP", "SCALE_DOWN", 
    "POWHEG", "MCATNLO", "SPINCORR"};
    
void Histo(bool doControlPlots, bool doPreunfold, bool doUnfold, 
           std::vector<std::string> plots, 
           std::vector<std::string> systematics, 
           std::vector<std::string> channels) 
{
    //to stay compatible with old code
    std::set<TString> SetOfValidSystematics;
    for (auto s: VectorOfValidSystematics) SetOfValidSystematics.insert(s);

    const double lumi = 12210;

    HistoListReader histoList(doControlPlots ? "HistoList_control" : "HistoList");
    if (histoList.IsZombie()) exit(12);
    for (auto it = histoList.begin(); it != histoList.end(); ++it) {
        const PlotProperties& p = it->second;
        cout << "checking " << p.name << endl;
        bool found = false;
        for (auto plot : plots) {
            if (p.name.Contains(plot, TString::kIgnoreCase)) {
                found = true;
                break;
            }
        }
        if (!found) continue;

        // Create Plotter 
        Plotter h_generalPlot;
        h_generalPlot.setLumi(lumi);
        h_generalPlot.ListOfSystematics(SetOfValidSystematics);
        
        /////////////////////////////////////////////////////
        /////////   UNFOLDING OPTIONS     ///////////////////
        /////////////////////////////////////////////////////

        TString outpath = "";
        h_generalPlot.UnfoldingOptions(doUnfold);
        h_generalPlot.SetOutpath("");

        /////////////////////////////////////////////////////
        /////////////////////////////////////////////////////
        ///////////////////////////////////////////////////// 
        
        h_generalPlot.setOptions(p.name,p.specialComment,p.ytitle,p.xtitle, 
                                    p.rebin, p.do_dyscale, p.logX, p.logY, 
                                    p.ymin, p.ymax, p.xmin, p.xmax, p.bins, p.xbinbounds, p.bincenters);
        h_generalPlot.DYScaleFactor();
        for (auto channel : channels) {
            for (auto systematic : systematics) {
                if (doPreunfold || doControlPlots || doUnfold) {
                    h_generalPlot.preunfolding(channel, systematic);
                }
                if (doControlPlots) {
                    h_generalPlot.MakeTable();
                }
                if (doUnfold) {
                    h_generalPlot.unfolding();
                    h_generalPlot.PlotDiffXSec(channel);
                }
            }
        }
    }

    cout << "Done with the unfolding\n";
}

/**
 * Helper function to create a function which checks if a string found is in the
 * passed vector of string.
 * 
 * @param allowed a vector of allowed strings (char*s)
 * @return a function taking a std::string and returning a bool
 */
std::function<bool(const std::string &s)> makeStringChecker(const std::vector<const char *> allowed) {
    return [allowed](const std::string &test) {
        return std::find(begin(allowed), end(allowed), test) != end(allowed);
    };
}

int main(int argc, char** argv) {
    CLParameter<std::string> opt_type("t", "cp|preunfold|unfold|p+u - required, cp=contol plots, p+u = preunfold+unfold", true, 1, 1,
        makeStringChecker({"cp", "preunfold", "unfold", "p+u"}));    
    CLParameter<std::string> opt_plots("p", "Name (pattern) of plot; multiple patterns possible", false, 1, 100);
    CLParameter<std::string> opt_channel("c", "Specify channel(s), valid: emu, ee, mumu, combined. Default: all channels", false, 1, 4,
        makeStringChecker({"ee", "emu", "mumu", "combined"}));
    CLParameter<std::string> opt_sys("s", "Systematic variation - default is Nominal", false, 1, 100,
        makeStringChecker(VectorOfValidSystematics));
    CLAnalyser::interpretGlobal(argc, argv);
    
    std::vector<std::string> channels { "emu", "ee", "mumu", "combined" };
    if (opt_channel.isSet()) channels = opt_channel.getArguments();
    std::cout << "Processing channels: "; 
    for (auto ch: channels) cout << ch << " "; cout << "\n";
        
    std::vector<std::string> systematics { "Nominal" };
    if (opt_sys.isSet()) systematics = opt_sys.getArguments();
    std::cout << "Processing systematics: "; 
    for (auto sys: systematics) cout << sys << " "; cout << "\n";
    
    std::vector<std::string> plots { "" };
    if (opt_plots.isSet()) plots = opt_plots.getArguments();

    bool doControlPlots = opt_type[0] == "cp";
    bool doPreunfold = opt_type[0] == "preunfold" || opt_type[0] == "p+u";
    bool doUnfold = opt_type[0] == "unfold" || opt_type[0] == "p+u";
    Histo(doControlPlots, doPreunfold, doUnfold, plots, systematics, channels);
}
