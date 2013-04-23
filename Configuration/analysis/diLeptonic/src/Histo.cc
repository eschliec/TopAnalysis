#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TMath.h>

#include <vector>
#include <iostream>
#include <set>

#include <future>

#include "plotterclass.h"
#include "HistoListReader.h"
#include "CommandLineParameters.h"

//constexpr double lumi = 12210; 
constexpr double lumi = 19624.8;

using namespace std;

const std::vector<const char*> VectorOfValidSystematics 
    {"Nominal", 
    "JER_UP", "JER_DOWN", "JES_UP", "JES_DOWN",
    "PU_UP", "PU_DOWN", "TRIG_UP", "TRIG_DOWN", "LEPT_UP", "LEPT_DOWN",
    "DY_UP", "DY_DOWN", "BG_UP", "BG_DOWN", 
    "KIN_UP", "KIN_DOWN",
    "BTAG_UP", "BTAG_DOWN", "BTAG_LJET_UP", "BTAG_LJET_DOWN",
    "BTAG_PT_UP", "BTAG_PT_DOWN", "BTAG_ETA_UP", "BTAG_ETA_DOWN",
    "BTAG_LJET_PT_UP", "BTAG_LJET_PT_DOWN", "BTAG_LJET_ETA_UP", "BTAG_LJET_ETA_DOWN",
//     "BTAG_BEFF_UP", "BTAG_BEFF_DOWN", "BTAG_CEFF_UP", "BTAG_CEFF_DOWN", "BTAG_LEFF_UP", "BTAG_LEFF_DOWN",
    "MASS_UP", "MASS_DOWN", "MATCH_UP", "MATCH_DOWN", "SCALE_UP", "SCALE_DOWN", 
    "POWHEG", "MCATNLO",// "SPINCORR", 
    "all"};
    
void Histo(bool doControlPlots, bool doUnfold, bool doDiffXSPlotOnly,
           std::vector<std::string> plots, 
           std::vector<std::string> systematics, 
           std::vector<std::string> channels) 
{
    //to stay compatible with old code
    std::set<TString> SetOfValidSystematics;
    for (auto s: VectorOfValidSystematics) SetOfValidSystematics.insert(s);

    HistoListReader histoList(doControlPlots ? "HistoList_control" : "HistoList");
    if (histoList.IsZombie()) exit(12);
    for (auto it = histoList.begin(); it != histoList.end(); ++it) {
        const PlotProperties& p = it->second;
        cout << "checking " << p.name << endl;
        bool found = false;
        for (auto plot : plots) {
            if (plot.size() && plot[0] == '+') {
                if (p.name.CompareTo(&plot[1], TString::kIgnoreCase) == 0) {
                    found = true;
                    break;
                }
            } else if (p.name.Contains(plot, TString::kIgnoreCase)) {
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
        h_generalPlot.DYScaleFactor(p.specialComment);
        //need preunfolding for ALL channels before unfolding!!
        for (auto channel : channels) {
            for (auto systematic : systematics) {
                h_generalPlot.preunfolding(channel, systematic);
            }
        }
        if (doUnfold) {
            if (!doDiffXSPlotOnly) {
                //unfold all channels except combined in parallel
                std::vector<std::future<void>> unfoldJobs;
                for (auto channel:channels) {
                    if (channel != "combined") {
                        TString ch(channel.c_str());
                        unfoldJobs.push_back(std::async(std::launch::async, [ch](Plotter p) -> void { 
                            p.unfolding(ch);
                        }, h_generalPlot));
                        //FIXME:
                        //it seems unfolding is not thread safe! 
                        //fix that and remove the next line
                        unfoldJobs.at(unfoldJobs.size()-1).wait();
                    }
                }
                //wait for the 3 channels to finish
                for (auto &i : unfoldJobs) {
                    i.get();
                }
                //only now do combined (if requested)
                for (auto channel:channels) {
                    if (channel == "combined") h_generalPlot.unfolding(channel);
                }
            }
            for (auto channel:channels){
                h_generalPlot.PlotDiffXSec(channel);
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
    CLParameter<std::string> opt_type("t", "cp|unfold|plot - required, cp=contol plots, unfold, or only plot diffXS", true, 1, 1,
        makeStringChecker({"cp", "unfold", "plot"}));
    CLParameter<std::string> opt_plots("p", "Name (pattern) of plot; multiple patterns possible; use '+Name' to match name exactly", false, 1, 100);
    CLParameter<std::string> opt_channel("c", "Specify channel(s), valid: emu, ee, mumu, combined. Default: all channels", false, 1, 4,
        makeStringChecker({"ee", "emu", "mumu", "combined"}));
    CLParameter<std::string> opt_sys("s", "Systematic variation - default is Nominal, use 'all' for all", false, 1, 100,
        makeStringChecker(VectorOfValidSystematics));
    CLAnalyser::interpretGlobal(argc, argv);
    
    std::vector<std::string> channels { "emu", "ee", "mumu", "combined" };
    if (opt_channel.isSet()) channels = opt_channel.getArguments();
    std::cout << "Processing channels: "; 
    for (auto ch: channels) cout << ch << " "; cout << "\n";
        
    std::vector<std::string> systematics { "Nominal" };
    if (opt_sys.isSet()) {
        systematics = opt_sys.getArguments();
        if (systematics[0] == "all") {
            systematics.clear();
            for (string syst: VectorOfValidSystematics) {
                if (syst != "all") systematics.push_back(syst);
            }
        }
    }    
    std::cout << "Processing systematics (use >>-s all<< to process all knwon systematics): "; 
    for (auto sys: systematics) cout << sys << " "; cout << "\n";
    
    std::vector<std::string> plots { "" };
    if (opt_plots.isSet()) plots = opt_plots.getArguments();

    bool doControlPlots = opt_type[0] == "cp";
    bool doDiffXSPlotOnly = opt_type[0] == "plot";
    bool doUnfold = opt_type[0] == "unfold" || doDiffXSPlotOnly;
    Histo(doControlPlots, doUnfold, doDiffXSPlotOnly, plots, systematics, channels);
}
