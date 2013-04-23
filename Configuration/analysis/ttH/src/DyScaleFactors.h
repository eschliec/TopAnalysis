#ifndef dyScaleFactors_h
#define dyScaleFactors_h

#include <vector>
#include <map>

class TString;
class TH1;

#include "sampleHelpers.h"

class RootFileReader;
class Sample;
class Samples;





class DyScaleFactors{

public:
    
    /// Typedef for the map containing the Drell-Yan scale factors
    typedef std::map<TString, std::map<Systematic::Systematic, std::map<Channel::Channel, double> > > DyScaleFactorMap;
    
    
    
    /// Constructor for producing Drell-Yan scale factors
    DyScaleFactors(const Samples& samples, const double& luminosity);
    
    /// Default destructor
    ~DyScaleFactors(){}
    
    
    
    /// Get the map containing the Drell-Yan scale factors
    DyScaleFactorMap getScaleFactors()const;
    
    /// Get Drell-Yan scale factor for given selection step, systematic and channel
    double dyScaleFactor(const TString& step, const Systematic::Systematic& systematic, const Channel::Channel& channel)const;
    
    /// Apply Drell-Yan scale factor to histogram automatically depending on its name
    void applyDyScaleFactor(TH1* histogram, const TString& histogramName,
                            const Sample& sample, const Systematic::Systematic& systematic,
                            const bool allowNonexistingStep)const;
    
    
    
private:
    
    /// Produce the Drell-Yan scale factors
    void produceScaleFactors(const Samples& samples);
    
    /// Produce the Drell-Yan scale factors for each selection step and each systematic
    void produceScaleFactors(const TString& step, const Systematic::Systematic& systematic, std::map<Channel::Channel, std::vector<Sample> >& channelSamples);
    
    /// Print full information about all ingoing numbers (deactivated by default)
    void printFullInformation(const double dyScaleFactor_ee, const double dyScaleFactor_mumu, 
                              const double k_ee, const double k_mumu,
                              const double rOutIn_ee, const double rOutIn_mumu,
                              const double nIn_ee_data_loose, const double nIn_mumu_data_loose,
                              const double nIn_ee_data, const double nIn_mumu_data, const double nIn_emu_data,
                              const double nIn_ee_mc, const double nIn_mumu_mc,
                              const double nIn_ee_dy, const double nIn_mumu_dy,
                              const double nOut_ee_mc, const double nOut_mumu_mc,
                              const double nOut_ee_dy, const double nOut_mumu_dy,
                              const TString& step)const;
    
    
    
    /// Luminosity
    const double& luminosity_;
    
    /// File reader for accessing specific histogram from given file
    RootFileReader* fileReader_;
    
    /// Map containing the Drell-Yan scale factors
    DyScaleFactorMap m_dyScaleFactors_;
};







#endif //dyScaleFactors_h





