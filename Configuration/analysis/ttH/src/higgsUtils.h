#ifndef tth_utils_h
#define tth_utils_h

#include <vector>
#include <utility>

class TString;





namespace tth{
    
    /// Get vector of pairs with < stepName, objectName from ROOT-file (e.g. histo or tree) >
    std::vector<std::pair<TString, TString> > nameStepPairs(const char* filename,
                                                            const char* objectNameFragment,
                                                            const std::vector<TString>& selectedSteps =std::vector<TString>());
    
    
    
    /// Assign a step name for a given short name of the step, and potentially for specific categories
    TString stepName(const TString& stepShort, const int& category =-1);
    
    /// Get from a TString the selection step of pattern "step*"
    TString extractSelectionStep(const TString& name);
    
    /// Get from a TString the jet category of pattern "cate*"
    TString extractJetCategory(const TString& name);
    
    /// Get from a TString the selection step of pattern "step*", combined with the category of pattern "bin_*",
    /// as it might be needed in case of binned control plots usage
    TString extractSelectionStepAndJetCategory(const TString& name);
    
    /// Helper functions only needed by functions defined in this file
    namespace helper{
        
        /// Helper function to get the fragment containing a searchPattern,
        /// fragments are separated by given token
        TString searchFragmentByToken(const TString& name, const TString& searchPattern, const TString& token);
    }
}





#endif




