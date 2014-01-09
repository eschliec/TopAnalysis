#ifndef BTAGBASE_H
#define BTAGBASE_H

#include <map>
#include "TH2D.h"
//#include "TString.h"
#include "../interface/miscUtils.h"
#include <iostream>
#include <string>


namespace ztop {
/**
 * workflow:
 *  -create object
 *  -set sample name              setSampleName(name)
 *  -fill efficiencies with jets  fillEff(p4,genPartonFlavour,PUweight)
 *  ---- get your event weight with your wrapper class!! ----
 *  ---functions to do that:
 *  - resetCounter()
 *  - countJet() (in jet loop)
 *  - getEventSF() (after jet loop - also resets the counter) returns 1 here
 *
 *
 *  -after event loop run makeEffs() once to create histograms
 *  - input / output should be organized by your wrapper class!
 *
 *  -load from file or use already filled object (implement in your wrapper!)
 *  -setMakeEff(false)
 * -(fillEff,makeEffs etc can remain in the loop/at the end of the loop, does nothing now)
 *  ---- get your event weight with your wrapper class!! ----
 *  ---functions to do that:
 *  - resetCounter()
 *  - countJet() (in jet loop)
 *  - getEventSF() (after jet loop - also resets the counter) now returns the proper SF!
 *
 *
 *
 *
 *  WHATEVER you add as functions, please don't use exit() in case an error occurs.
 *  replace it with either:
 *  - throw an exception (throw std::logic_error("sometext") or std::runtime_error("");)
 *  - return something (-1 or another int for dubugging)
 *
 */
class bTagBase {
public:

    bTagBase();

    ~bTagBase() {
        cleanptr();
    }

    /**
     * enum for working points
     * should be named <tagger><wp>_wp
     *
     */
    enum workingPoints {
        csvt_wp, csvm_wp, csvl_wp,
        /* some other WP/taggers ..*/
        //
        /*don't touch this: needs to be last entry*/
        length_wp
    };

    /**
     * enum for systematic variations
     */
    enum systematics {
        nominal, heavyup, heavydown, lightup, lightdown
    };

    void setWorkingPoint(workingPoints wp) { wp_ = wp; }
    const float& getWPDiscrValue()const{return wpvals_[wp_];}

    workingPoints getWorkingPoint() const { return wp_;}
    /**
     * switches on SF for 7 TeV data
     */
    void setIs2011(bool is) {is2011_ = is;}

    void setSystematic(systematics sys) { syst_ = sys;}
    systematics getSystematic() { return syst_;}


    int setSampleName(const std::string &); //checks if effs should be made, if sample exists,..

    /**
     * enables the filling of histograms for efficiencies and disables
     * scale factor output
     */
    void setMakeEff(bool makee) {makeeffs_ = makee; }
    bool getMakeEff() { return makeeffs_;}

    void fillEff(const float &, const float&, const int &, const float &,const float&);
    void makeEffs();

    /**
     * returns the number of nan sf
     */
    size_t getNanCount(){return nancount_;}

    bool showWarnings;
    bool debug;

protected:

    std::map<std::string, std::vector<TH2D> > histos_; //! bjets, btagged, cjets, ctagged, ljets, ltagged
    std::map<std::string, std::vector<TH2D> > effhistos_;   //! beff, ceff, leff

    // functions to be used in the wrapper classes inheriting from this class
    // example implementation:
    /*
     resetCounter();
     for(size_t i=0;i<jets.size();i++){
     ztop::NTJet *jet=jets.at(i);
     countJet(jet->p4(),jet->genPartonFlavour());
     }
     return getEventSF();
     */

    void resetCounter() {
        sumStuffEff_ = 0.99999999;
        sumStuffSfEff_ = 0.99999999;
    }
    void countJet(const float&, const float&, const int & genPartonFlavor);
    float getEventSF();

    bool jetIsTagged(const float&, const float&, const int & genPartonFlavor,
            const float & tagValue, const unsigned int & seed) const; // to be implmented

    void cleanptr() {
        histp_ = 0;
        effhistp_ = 0;
    }




private:
    enum jetTypes{bjet,cjet,lightjet,undefined};

    bool init_;

    workingPoints wp_;
    bool is2011_;
    systematics syst_;

    std::vector<float> wpvals_, minpt_, maxpt_;

    bool makeeffs_;

    std::string tempsamplename_;
    std::vector<TH2D> * histp_;
    std::vector<TH2D> * effhistp_;

    float sumStuffEff_;
    float sumStuffSfEff_;

    //counts the fraction of nans
    size_t nancount_;

    jetTypes jetType(const int & )const;

    float calcHeavySF(float *, float *, const size_t &, const float &,
            const float &, const float &) const;

    float jetSF(const float &pt, const float& abs_eta,const jetTypes & jettype) const;
    float jetEff(const float &pt, const float& abs_eta,const jetTypes & jettype) const;

    //SFs btv input
    float BJetSF(const float &pt, const float& abs_eta,
            float multiplier = 1) const;
    float CJetSF(const float &pt, const float &abs_eta,
            float multiplier = 1) const;
    float LJetSF(const float &pt, const float &abs_eta,
            float multiplier = 1) const;

    void initWorkingpoints();

};

///inlined functions for performance reasons:
inline bTagBase::jetTypes bTagBase::jetType(const int & partonflavor)const{
    if(fabs(fabs(partonflavor) - 5) < 0.1) return bjet;
    else if (fabs(fabs(partonflavor) - 4) < 0.1) return cjet;
    else if(fabs(partonflavor)>0) return lightjet;
    else return undefined;
}

inline float bTagBase::calcHeavySF(float* ptbins, float * SFb_error,
        const size_t & ptbinsSize, const float & pt, const float & SF,
        const float & multiplier) const {

    if (syst_ != heavyup && syst_ != heavydown)
        return SF;

    //this uses the standard histogram bin range definition from root
    size_t ptbin = std::lower_bound(ptbins, ptbins + ptbinsSize, pt)
    - &ptbins[0];
    if (ptbins[ptbin] == pt)
        --ptbin;

    if (syst_ == heavyup)
        return SF + (multiplier * SFb_error[ptbin]);
    if (syst_ == heavydown)
        return SF - (multiplier * SFb_error[ptbin]);

    return 0;    //never reaches
}

}    //namespace
////////////////////////////////////
#endif

