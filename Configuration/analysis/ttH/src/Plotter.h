#ifndef plotterclass_h
#define plotterclass_h

#include <vector>
#include <map>

class TLegend;
class RootFileReader;
class TString;
class TH1;

#include "plotterHelpers.h"
#include "sampleHelpers.h"
#include "SamplesFwd.h"

class DyScaleFactors;




class Plotter {

public:
    
    /// Constructor
    Plotter(const Samples& samples,
            const double& luminosity,
            const DyScaleFactors& dyScaleFactors,
            const DrawMode::DrawMode& drawMode);
    
    /// Destructor
    ~Plotter(){};
    
    /// Set options specific for histogram under consideration
    void setOptions(const TString&, const TString&, const TString&, const TString&,
                    const int, const bool, const bool, const bool,
                    const double&, const double&, const double&, const double&,
                    const int, const std::vector<double>&, const std::vector<double>&);
    
    /// Produce the plots for histogram under consideration from all samples
    void producePlots();
    
    
    
private:
    
    /// Pair of a legend entry and the histogram for the corresponding sample
    typedef std::pair<TString, TH1*> LegendHistPair;
    
    /// Access histogram under consideration from all samples, do global scalings
    bool prepareDataset(const std::vector<Sample>&, const Systematic::Systematic& systematic);
    
    /// Do stacking, legending, and write in file
    void write(const Channel::Channel&, const Systematic::Systematic&);
    
    
    
    /// Prepare a separate legend for the plots, which can have e.g. a different order
    TLegend* controlLegend(const LegendHistPair& dataHist, const std::vector<LegendHistPair>& stackHists,
                           const std::vector<LegendHistPair>& higgsHists, TLegend* leg);
    
    /// Set the style of the plot
    void setStyle(SampleHistPair&, const bool isControlPlot=false);
    
    /// Draw label for decay channel in upper left corner of plot
    void drawDecayChannelLabel(const Channel::Channel&, const double& textSize=0.04);
    
    /// Draw official labels (CMS [Preliminary], luminosity and CM energy) above plot
    void drawCmsLabels(const int cmsprelim=1, const double& energy=8, const double& textSize=0.04);
    
    
    
    /// Samples to be analysed
    const Samples& samples_;
    
    /// Luminosity
    const double& luminosity_;
    
    /// Drell-Yan scale factors
    const DyScaleFactors& dyScaleFactors_;
    
    /// Draw mode for Higgs
    const DrawMode::DrawMode drawMode_;
    
    /// File reader for accessing specific histogram from given file
    RootFileReader* fileReader_;
    
    
    
    /// Pair of histogram to print and corresponding sample
    std::vector<SampleHistPair> v_sampleHistPair_;
    
    
    
    /// Name of histogram under consideration
    TString name_;
    
    /// Options for the histogram under consideration
    int bins_, rebin_;
    double rangemin_, rangemax_, ymin_, ymax_;
    std::vector<double> XAxisbins_, XAxisbinCenters_;
    TString YAxis_;
    TString XAxis_;
    bool logX_, logY_, doDYScale_; // The variable logX_ is not used at all...
};







#endif







