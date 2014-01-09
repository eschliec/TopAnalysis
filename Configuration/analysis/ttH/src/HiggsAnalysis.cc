#define HiggsAnalysis_cxx

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <algorithm>

#include <TTree.h>
#include <TString.h>
#include <TFile.h>
#include <TSystem.h>
#include <Math/VectorUtil.h>

#include "HiggsAnalysis.h"
#include "higgsUtils.h"
#include "DijetAnalyzer.h"
#include "analysisStructs.h"
#include "AnalysisHistograms.h"
#include "Playground.h"
#include "MvaTreeHandler.h"
#include "MvaTreeAnalyzer.h"
#include "MvaValidation.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/classes.h"
#include "../../common/include/ScaleFactors.h"
#include "../../common/include/BTagUtils.h"





/// Lepton eta selection (absolute value)
constexpr double LeptonEtaCUT = 2.4;

/// Lepton pt selection in GeV
constexpr double LeptonPtCut = 20.;

/// Jet eta selection (absolute value)
constexpr double JetEtaCUT = 2.4;

/// Jet pt selection in GeV
constexpr double JetPtCUT = 30.;


/// B-tag working point
/// Available options: 
///  csvl_wp (0.244) 
///  csvm_wp (0.679)
///  csvt_wp (0.898)
/// Used only if new btag scaling is used (useGenericBTagSF = true; in load_Analysis.cc)
constexpr BTagSFGeneric::workingPoints BtagWP = BTagSFGeneric::csvl_wp;

/// Used only if default old btag scaling is used (useGenericBTagSF = false; in load_Analysis.cc)
/// FIXME Should be removed when the old method is deprecated
constexpr double BtagWP_val = 0.244;


/// MET selection for same-flavour channels (ee, mumu)
constexpr double MetCUT = 40.;










HiggsAnalysis::HiggsAnalysis(TTree*):
isInclusiveHiggs_(false),
bbbarDecayFromInclusiveHiggs_(false),
runWithTtbb_(false),
retagBJets_(true),
mvaTreeHandler_(0)
{}



HiggsAnalysis::~HiggsAnalysis()
{}



void HiggsAnalysis::Begin(TTree*)
{
    // Defaults from AnalysisBase
    AnalysisBase::Begin(0);

    // Prepare things for analysis
    this->prepareJER_JES();

    // Set up selection steps of MVA tree handler
    if(mvaTreeHandler_) mvaTreeHandler_->book();
}




void HiggsAnalysis::Terminate()
{
    // Produce b-tag efficiencies
    // FIXME: runWithTtbb_ is dirty hack, since makeBtagEfficiencies() is in AnalysisBase
    // FIXME: Shouldn't we also clear b-tagging efficiency histograms if they are produced ?
    if(!runWithTtbb_ && this->makeBtagEfficiencies() && btagScaleFactors_) btagScaleFactors_->produceBtagEfficiencies(static_cast<std::string>(channel_));
    else if(!runWithTtbb_ && this->makeBtagEfficiencies() && bTagSFGeneric_) bTagSFGeneric_->produceBtagEfficiencies(static_cast<std::string>(channel_));

    // Do everything needed for MVA
    if(mvaTreeHandler_){
        // Produce and write tree
        mvaTreeHandler_->writeTrees(static_cast<std::string>(outputfilename_),
                                    Channel::convertChannel(static_cast<std::string>(channel_)),
                                    Systematic::convertSystematic(static_cast<std::string>(systematic_)));
        //mvaTreeHandler_->writeTrees(fOutput);

        // Create and store control plots in fOutput
        MvaTreeAnalyzer mvaTreeAnalyzer(mvaTreeHandler_->stepMvaVariablesMap());
        mvaTreeAnalyzer.plotVariables(fOutput);
        mvaTreeAnalyzer.clear();

        // Cleanup
        mvaTreeHandler_->clear();
    }


    // Defaults from AnalysisBase
    AnalysisBase::Terminate();
}



void HiggsAnalysis::SlaveBegin(TTree *)
{
    // Defaults from AnalysisBase
    AnalysisBase::SlaveBegin(0);

    // Histograms for b-tagging efficiencies
    if(!runWithTtbb_ && btagScaleFactors_ && this->makeBtagEfficiencies()) btagScaleFactors_->bookBtagHistograms(fOutput, static_cast<std::string>(channel_));
    else if(!runWithTtbb_ && bTagSFGeneric_) {
        bTagSFGeneric_->setWorkingPoint(BtagWP);
        bTagSFGeneric_->prepareBTags(fOutput, static_cast<std::string>(channel_));
    }

    // Book histograms of all analyzers
    this->bookAll();
}



void HiggsAnalysis::SlaveTerminate()
{
    this->clearAll();

    // Defaults from AnalysisBase
    AnalysisBase::SlaveTerminate();
}



Bool_t HiggsAnalysis::Process(Long64_t entry)
{
    // Defaults from AnalysisBase
    if(!AnalysisBase::Process(entry)) return kFALSE;


    // Use utilities without namespaces
    using namespace common;


    // Entry for object structs are not yet read, so reset
    this->resetObjectStructEntry();


    // Define the selection steps as strings
    std::string selectionStep("");


    //===CUT===
    // this is step0a, no cut application
    selectionStep = "0a";

    // ++++ Control Plots ++++



    //===CUT===
    // this is step0b, select events on generator level and access true level weights
    selectionStep = "0b";

    // Separate DY dilepton decays in lepton flavours
    if(this->failsDrellYanGeneratorSelection(entry)) return kTRUE;

    // Separate dileptonic ttbar decays via tau
    //if(this->failsTopGeneratorSelection(entry)) return kTRUE;

    // Separate inclusive ttH sample in decays H->bbbar and others
    if(this->failsHiggsGeneratorSelection(entry)) return kTRUE;

    // Separate tt+bb from tt+other
    if(this->failsAdditionalJetFlavourSelection(entry)) return kTRUE;

    // Correct for the MadGraph branching fraction being 1/9 for dileptons (PDG average is .108)
    const double weightMadgraphCorrection = this->madgraphWDecayCorrection(entry);

    // Get weight due to pileup reweighting
    const double weightPU = this->weightPileup(entry);

    // Get weight due to generator weights
    const double weightGenerator = this->weightGenerator(entry);

    // Get true level weights
    const double trueLevelWeightNoPileup = weightGenerator*weightMadgraphCorrection;
    const double trueLevelWeight = trueLevelWeightNoPileup*weightPU;

    const tth::GenLevelWeights genLevelWeights(weightMadgraphCorrection, weightPU, weightGenerator,
                                               trueLevelWeightNoPileup, trueLevelWeight);

    // ++++ Control Plots ++++



    //===CUT===
    selectionStep = "1";

    // Check if event was triggered with the same dilepton trigger as the specified analysis channel
    if(this->failsDileptonTrigger(entry)) return kTRUE;



    // === FULL RECO OBJECT SELECTION === (can thus be used at each selection step)

    // Access reco objects, and common generator objects
    const RecoObjects& recoObjects = this->getRecoObjects(entry);
    const CommonGenObjects& commonGenObjects = this->getCommonGenObjects(entry);

    // Create dummies for other objects, non-dummies are created only as soon as needed
    const KinRecoObjects kinRecoObjectsDummy;
    const TopGenObjects topGenObjectsDummy;
    const HiggsGenObjects higgsGenObjectsDummy;

    // Get allLepton indices, apply selection cuts and order them by pt (beginning with the highest value)
    const VLV& allLeptons = *recoObjects.allLeptons_;
    std::vector<int> allLeptonIndices = initialiseIndices(allLeptons);
    selectIndices(allLeptonIndices, allLeptons, LVeta, LeptonEtaCUT, false);
    selectIndices(allLeptonIndices, allLeptons, LVeta, -LeptonEtaCUT);
    selectIndices(allLeptonIndices, allLeptons, LVpt, LeptonPtCut);
    orderIndices(allLeptonIndices, allLeptons, LVpt);
    //const int numberOfAllLeptons = allLeptonIndices.size();

    // Get indices of leptons and antiLeptons separated by charge, and get the leading ones if they exist
    const std::vector<int>& lepPdgId = *recoObjects.lepPdgId_;
    std::vector<int> leptonIndices = allLeptonIndices;
    std::vector<int> antiLeptonIndices = allLeptonIndices;
    selectIndices(leptonIndices, lepPdgId, 0);
    selectIndices(antiLeptonIndices, lepPdgId, 0, false);
    const int numberOfLeptons = leptonIndices.size();
    const int numberOfAntiLeptons = antiLeptonIndices.size();
    const int leptonIndex = numberOfLeptons>0 ? leptonIndices.at(0) : -1;
    const int antiLeptonIndex = numberOfAntiLeptons>0 ? antiLeptonIndices.at(0) : -1;

    // In case of an existing opposite-charge dilepton system,
    // get their indices for leading and next-to-leading lepton
    int leadingLeptonIndex(-1);
    int nLeadingLeptonIndex(-1);
    if(numberOfLeptons>0 && numberOfAntiLeptons>0){
        leadingLeptonIndex = leptonIndex;
        nLeadingLeptonIndex = antiLeptonIndex;
        orderIndices(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, LVpt);
    }
    const bool hasLeptonPair = this->hasLeptonPair(leadingLeptonIndex, nLeadingLeptonIndex, lepPdgId);

    // Get two indices of the two leptons in the right order for trigger scale factor, if existing
    int leptonXIndex(leadingLeptonIndex);
    int leptonYIndex(nLeadingLeptonIndex);
    if(hasLeptonPair){
        //in ee and mumu channel leptonX must be the highest pt lepton, i.e. this is already correct
        // in emu channel leptonX must be electron
        if(std::abs(lepPdgId.at(leptonXIndex)) != std::abs(lepPdgId.at(leptonYIndex))){
            orderIndices(leptonYIndex, leptonXIndex, lepPdgId, true);
        }
    }

    // Get dilepton system, if existing
    const LV dummyLV(0.,0.,0.,0.);
    const LV dilepton(hasLeptonPair ? allLeptons.at(leadingLeptonIndex)+allLeptons.at(nLeadingLeptonIndex) : dummyLV);

    // Get jet indices, apply selection cuts and order them by pt (beginning with the highest value)
    const VLV& jets = *recoObjects.jets_;
    std::vector<int> jetIndices = initialiseIndices(jets);
    selectIndices(jetIndices, jets, LVeta, JetEtaCUT, false);
    selectIndices(jetIndices, jets, LVeta, -JetEtaCUT);
    selectIndices(jetIndices, jets, LVpt, JetPtCUT);
    orderIndices(jetIndices, jets, LVpt);
    const int numberOfJets = jetIndices.size();
    const bool has2Jets = numberOfJets > 1;

    // Fill a vector with all jet pair indices, while sorting each pair by the jet charge:
    // first entry is antiBIndex i.e. with higher jet charge, second entry is bIndex
    //const std::vector<double>& jetChargeGlobalPtWeighted = *recoObjects.jetChargeGlobalPtWeighted_;
    const std::vector<double>& jetChargeRelativePtWeighted = *recoObjects.jetChargeRelativePtWeighted_;
    const tth::IndexPairs& jetIndexPairs = this->chargeOrderedJetPairIndices(jetIndices, jetChargeRelativePtWeighted);

    // Get b-jet indices, apply selection cuts
    // and order b-jets by btag discriminator (beginning with the highest value)
    const std::vector<double>& jetBTagCSV = *recoObjects.jetBTagCSV_;
    const std::vector<int>& jetPartonFlavour = *commonGenObjects.jetPartonFlavour_;
    std::vector<int> bjetIndices = jetIndices;
    if(bTagSFGeneric_) selectIndices(bjetIndices, jetBTagCSV, (double)bTagSFGeneric_->getWPDiscrValue());
    else selectIndices(bjetIndices, jetBTagCSV, BtagWP_val);
    if(retagBJets_) {
        // FIXME Should be removed when the old method is deprecated
        if (isMC_ && btagScaleFactors_ && !(btagScaleFactors_->makeEfficiencies())){
            // Apply b-tag efficiency MC correction using random number based tag flipping
            btagScaleFactors_->indexOfBtags(bjetIndices, jetIndices,
                                            jets, jetPartonFlavour, jetBTagCSV,
                                            BtagWP_val, static_cast<std::string>(channel_));
        }
        else if (isMC_ && bTagSFGeneric_ && !(bTagSFGeneric_->makeEfficiencies())){
            // Apply b-tag efficiency MC correction using random number based tag flipping
            bTagSFGeneric_->indexOfBtags(bjetIndices, jetIndices,
                                         jets, jetPartonFlavour, jetBTagCSV);
        }
    }
    orderIndices(bjetIndices, jetBTagCSV);
    const int numberOfBjets = bjetIndices.size();
    const bool hasBtag = numberOfBjets > 0;

    // Get MET
    const LV& met = *recoObjects.met_;
    const bool hasMetOrEmu = channel_=="emu" || met.pt()>MetCUT;

    const tth::RecoObjectIndices recoObjectIndices(allLeptonIndices,
                                                   leptonIndices, antiLeptonIndices,
                                                   leptonIndex, antiLeptonIndex,
                                                   leadingLeptonIndex, nLeadingLeptonIndex,
                                                   leptonXIndex, leptonYIndex,
                                                   jetIndices, jetIndexPairs,
                                                   bjetIndices);

    const tth::GenObjectIndices genObjectIndicesDummy(-1, -1, -1, -1, -1, -1, -1, -1);


    // Determine all reco level weights
    const double weightLeptonSF = this->weightLeptonSF(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId);
    const double weightTriggerSF = this->weightTriggerSF(leptonXIndex, leptonYIndex, allLeptons);
    const double weightNoPileup = trueLevelWeightNoPileup*weightTriggerSF*weightLeptonSF;
    // We do not apply a b-tag scale factor
    //const double weightBtagSF = ReTagJet ? 1. : this->weightBtagSF(jetIndices, jets, jetPartonFlavour);
    double weightBtagSF = 1.0;
    if(!retagBJets_) {
        // FIXME Should be removed when the old method is deprecated
        if(btagScaleFactors_) weightBtagSF = this->weightBtagSF(jetIndices, jets, jetPartonFlavour);
        else if(bTagSFGeneric_ && isMC_) weightBtagSF = bTagSFGeneric_->calculateBtagSF(jetIndices, jets, jetPartonFlavour);
    }

    // The weight to be used for filling the histograms
    double weight = weightNoPileup*weightPU;


    tth::RecoLevelWeights recoLevelWeights(weightLeptonSF, weightTriggerSF, weightBtagSF,
                                           weightNoPileup, weight);


    // ++++ Control Plots ++++

    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjectsDummy, higgsGenObjectsDummy,
                  kinRecoObjectsDummy,
                  genObjectIndicesDummy, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  1.);



    //===CUT===
    selectionStep = "2";

    // we need an OS lepton pair matching the trigger selection...
    if (!hasLeptonPair) return kTRUE;

    // ++++ Control Plots ++++

    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjectsDummy, higgsGenObjectsDummy,
                  kinRecoObjectsDummy,
                  genObjectIndicesDummy, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);



    //===CUT===
    selectionStep = "3";

    // ...with at least 20 GeV invariant mass
    if(dilepton.M() < 20.) return kTRUE;

    const bool isZregion = dilepton.M() > 76 && dilepton.M() < 106;
    //const KinRecoObjects& kinRecoObjects = this->getKinRecoObjects(entry);
    const KinRecoObjects& kinRecoObjects = this->getKinRecoObjectsOnTheFly(leptonIndex, antiLeptonIndex, jetIndices,
                                                                           allLeptons, jets, jetBTagCSV, met);
    //const bool hasSolution = kinRecoObjects.valuesSet_;



    // ++++ Control Plots ++++

    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjectsDummy, higgsGenObjectsDummy,
                  kinRecoObjects,
                  genObjectIndicesDummy, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);



    // ****************************************
    //handle inverted Z cut
    // Z window plots need to be filled here, in order to rescale the contribution to data
    if(isZregion){
        selectionStep = "4zWindow";

        this->fillAll(selectionStep,
                      recoObjects, commonGenObjects,
                      topGenObjectsDummy, higgsGenObjectsDummy,
                      kinRecoObjects,
                      genObjectIndicesDummy, recoObjectIndices,
                      genLevelWeights, recoLevelWeights,
                      weight);

        if(has2Jets){
            selectionStep = "5zWindow";

            this->fillAll(selectionStep,
                          recoObjects, commonGenObjects,
                          topGenObjectsDummy, higgsGenObjectsDummy,
                          kinRecoObjects,
                          genObjectIndicesDummy, recoObjectIndices,
                          genLevelWeights, recoLevelWeights,
                          weight);

            if(hasMetOrEmu){
                selectionStep = "6zWindow";

                this->fillAll(selectionStep,
                              recoObjects, commonGenObjects,
                              topGenObjectsDummy, higgsGenObjectsDummy,
                              kinRecoObjects,
                              genObjectIndicesDummy, recoObjectIndices,
                              genLevelWeights, recoLevelWeights,
                              weight);

                if(hasBtag){
                    selectionStep = "7zWindow";

                    // FIXME: do not use b-tag scale factor
                    //weightBtagSF = isMC ? calculateBtagSF() : 1;
                    //fullWeights *= weightBtagSF;

                    this->fillAll(selectionStep,
                                  recoObjects, commonGenObjects,
                                  topGenObjectsDummy, higgsGenObjectsDummy,
                                  kinRecoObjects,
                                  genObjectIndicesDummy, recoObjectIndices,
                                  genLevelWeights, recoLevelWeights,
                                  weight);
                }
            }
        }
    }



    //=== CUT ===
    selectionStep = "4";

    //Exclude the Z window
    if(channel_!="emu" && isZregion) return kTRUE;

    // ++++ Control Plots ++++

    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjectsDummy, higgsGenObjectsDummy,
                  kinRecoObjects,
                  genObjectIndicesDummy, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);



    //=== CUT ===
    selectionStep = "5";

    //Require at least two jets
    if(!has2Jets) return kTRUE;

    // ++++ Control Plots ++++

    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjectsDummy, higgsGenObjectsDummy,
                  kinRecoObjects,
                  genObjectIndicesDummy, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);



    //=== CUT ===
    selectionStep = "6";

    //Require MET > 40 GeV in non-emu channels
    if(!hasMetOrEmu) return kTRUE;

    // ++++ Control Plots ++++

    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjectsDummy, higgsGenObjectsDummy,
                  kinRecoObjects,
                  genObjectIndicesDummy, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);

    // Fill the b-tagging efficiency plots
    if(!runWithTtbb_ && this->makeBtagEfficiencies() && btagScaleFactors_){
        btagScaleFactors_->fillBtagHistograms(jetIndices, bjetIndices,
                                              jets, jetPartonFlavour,
                                              weight, static_cast<std::string>(channel_));
    }
    else if(!runWithTtbb_ && this->makeBtagEfficiencies() && bTagSFGeneric_){
        bTagSFGeneric_->fillBtagHistograms(jetIndices, jetBTagCSV,
                                           jets, jetPartonFlavour,
                                           weight);
    }



    //=== CUT ===
    selectionStep = "7";

    //Require at least one b tagged jet
    if(!hasBtag) return kTRUE;

    weight *= weightBtagSF;



    // === FULL GEN OBJECT SELECTION ===

    // Access top generator object struct, and higgs generator object struct
    const TopGenObjects& topGenObjects = this->getTopGenObjects(entry);
    const HiggsGenObjects& higgsGenObjects = this->getHiggsGenObjects(entry);

    // Do jet matchings for ttbar system
    int genBjetFromTopIndex(-1);
    int genAntiBjetFromTopIndex(-1);
    int matchedBjetFromTopIndex(-1);
    int matchedAntiBjetFromTopIndex(-1);
    if(isTopSignal_){
        const VLV& allGenJets = *commonGenObjects.allGenJets_;
        // Find gen-level b jet and anti-b jet corresponding to (anti)b from (anti)top
        // FIXME: should one clean the genJetCollection to remove low-pt (or high-eta) jets?
        if(this->getGenBjetIndices(genBjetFromTopIndex, genAntiBjetFromTopIndex, topGenObjects, 6)){
            // Match recoJets to the two selected genJets from (anti)top
            this->matchRecoToGenJets(matchedBjetFromTopIndex, matchedAntiBjetFromTopIndex,
                                     jetIndices,
                                     jets,
                                     &allGenJets.at(genBjetFromTopIndex), &allGenJets.at(genAntiBjetFromTopIndex));
        }
    }

    // Do jet matchings for Higgs system
    int genBjetFromHiggsIndex(-1);
    int genAntiBjetFromHiggsIndex(-1);
    int matchedBjetFromHiggsIndex(-1);
    int matchedAntiBjetFromHiggsIndex(-1);
    // FIXME: should not use higgsDecayMode_, but specific method
    if(isHiggsSignal_ && higgsDecayMode_==5){
        const VLV& allGenJets = *commonGenObjects.allGenJets_;
        // Find gen-level b jet and anti-b jet corresponding to (anti)b from Higgs
        // FIXME: should one clean the genJetCollection to remove low-pt (or high-eta) jets?
        if(this->getGenBjetIndices(genBjetFromHiggsIndex, genAntiBjetFromHiggsIndex, topGenObjects, 25)){
            // Match recoJets to the two selected genJets from Higgs
            this->matchRecoToGenJets(matchedBjetFromHiggsIndex, matchedAntiBjetFromHiggsIndex,
                                     jetIndices,
                                     jets,
                                     &allGenJets.at(genBjetFromHiggsIndex), &allGenJets.at(genAntiBjetFromHiggsIndex));
        }
    }

    const tth::GenObjectIndices genObjectIndices(genBjetFromTopIndex, genAntiBjetFromTopIndex,
                                                 matchedBjetFromTopIndex, matchedAntiBjetFromTopIndex,
                                                 genBjetFromHiggsIndex, genAntiBjetFromHiggsIndex,
                                                 matchedBjetFromHiggsIndex, matchedAntiBjetFromHiggsIndex);



    // ++++ Control Plots ++++

    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjects, higgsGenObjects,
                  kinRecoObjects,
                  genObjectIndices, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);



    return kTRUE;
}



tth::IndexPairs HiggsAnalysis::chargeOrderedJetPairIndices(const std::vector<int>& jetIndices,
                                                           const std::vector<double>& jetCharges)
{
    tth::IndexPairs result;
    if(jetIndices.size() < 2) return result;

    // Loop over all jet combinations
    for(std::vector<int>::const_iterator i_jetIndex = jetIndices.begin(); i_jetIndex != --(jetIndices.end()); ++i_jetIndex){
        std::vector<int>::const_iterator incrementIterator(i_jetIndex);
        ++incrementIterator;
        for(std::vector<int>::const_iterator j_jetIndex = incrementIterator; j_jetIndex != jetIndices.end(); ++j_jetIndex){
            // Get the indices of b and anti-b jet defined by jet charge
            int bIndex = *i_jetIndex;
            int antiBIndex = *j_jetIndex;
            common::orderIndices(antiBIndex, bIndex, jetCharges);

            result.push_back(std::make_pair(antiBIndex, bIndex));
        }
    }

    return result;
}



bool HiggsAnalysis::getGenBjetIndices(int& genBjetIndex, int& genAntiBjetIndex,
                                      const TopGenObjects& topGenObjects, const int pdgId)
{
    if(!pdgId>0){
        std::cerr<<"ERROR! Method getGenBJetIndices needs a pdgId>0, but used is: "<<pdgId
                 <<"\n...break\n\n";
        exit(71);
    }

    const std::vector<int>& genBHadFlavour(*topGenObjects.genBHadFlavour_);
    const std::vector<int>& genBHadJetIndex(*topGenObjects.genBHadJetIndex_);

    for(size_t iBHadron=0; iBHadron<genBHadFlavour.size(); ++iBHadron){
        const int flavour = genBHadFlavour.at(iBHadron);
        if(std::abs(flavour) != std::abs(pdgId)) continue;     // Skipping hadrons with the wrong flavour
        // Assigning jet index of corresponding hadron. Set to -2 if >1 hadrons found for the same flavour
        if(flavour>0) genBjetIndex = (genBjetIndex==-1) ? genBHadJetIndex.at(iBHadron) : -2;
        else if(flavour<0) genAntiBjetIndex = (genAntiBjetIndex==-1) ? genBHadJetIndex.at(iBHadron) : -2;
    }

    // If no unique match of jets from (anti)b from (anti)top is found, return false
    if(genBjetIndex<0 || genAntiBjetIndex<0 || genBjetIndex==genAntiBjetIndex){
        return false;
    }
    return true;
}



bool HiggsAnalysis::matchRecoToGenJets(int& matchedBjetIndex, int& matchedAntiBjetIndex,
                                       const std::vector<int>& jetIndices,
                                       const VLV& jets,
                                       const LV* genBjet, const LV* genAntiBjet)
{
    using ROOT::Math::VectorUtil::DeltaR;

    // Find closest jet and its distance in deltaR
    double deltaRBjet(999.);
    double deltaRAntiBjet(999.);
    for(const auto& index : jetIndices){
        float deltaR = DeltaR(*genBjet, jets.at(index));
        if(deltaR < deltaRBjet){
            deltaRBjet = deltaR;
            matchedBjetIndex = index;
        }
        deltaR = DeltaR(*genAntiBjet, jets.at(index));
        if(deltaR < deltaRAntiBjet){
            deltaRAntiBjet = deltaR;
            matchedAntiBjetIndex = index;
        }
    }

    // Call a jet matched if it is close enough (should this be a configurable parameter?)
    if(deltaRBjet>0.5) matchedBjetIndex = -1;
    if(deltaRAntiBjet>0.5) matchedAntiBjetIndex = -1;

    // Check if both gen jets are successfully matched to different reco jets
    if(matchedBjetIndex==-1 || matchedAntiBjetIndex==-1 || matchedBjetIndex == matchedAntiBjetIndex) return false;

    return true;
}



void HiggsAnalysis::SetHiggsInclusiveSample(const bool isInclusiveHiggs)
{
    isInclusiveHiggs_ = isInclusiveHiggs;
}



void HiggsAnalysis::SetHiggsInclusiveSeparation(const bool bbbarDecayFromInclusiveHiggs)
{
    bbbarDecayFromInclusiveHiggs_ = bbbarDecayFromInclusiveHiggs;
}



void HiggsAnalysis::SetRunWithTtbb(const bool runWithTtbb)
{
    runWithTtbb_ = runWithTtbb;
}



void HiggsAnalysis::SetMvaInputProduction(MvaTreeHandler* mvaTreeHandler)
{
    mvaTreeHandler_ = mvaTreeHandler;
}



void HiggsAnalysis::SetAllAnalysisHistograms(std::vector<AnalysisHistogramsBase*> v_analysisHistograms)
{
    v_analysisHistograms_ = v_analysisHistograms;
}



bool HiggsAnalysis::failsHiggsGeneratorSelection(const Long64_t& entry)const
{
    if(!isHiggsSignal_) return false;
    GetHiggsDecayModeEntry(entry);
    // Separate ttH events from inclusve decay into H->bbbar and other decays
    // FIXME: do not use directly higgsDecayMode_, but a function
    if(isInclusiveHiggs_ && !bbbarDecayFromInclusiveHiggs_ && higgsDecayMode_==5) return true;
    if(isInclusiveHiggs_ && bbbarDecayFromInclusiveHiggs_ && higgsDecayMode_!=5) return true;
    return false;
}



bool HiggsAnalysis::failsAdditionalJetFlavourSelection(const Long64_t& entry)const
{
    if(!isTopSignal_) return false;
    if(isHiggsSignal_) return false;

    // FIXME: this is a workaround as long as there is no specific additional jet flavour info written to nTuple
    const TopGenObjects& topGenObjects = this->getTopGenObjects(entry);
    const int nGenBJets = topGenObjects.genBHadIndex_->size();
    if(runWithTtbb_ && nGenBJets<=2) return true;
    if(!runWithTtbb_ && nGenBJets>2) return true;

    return false;
}



void HiggsAnalysis::fillAll(const std::string& selectionStep,
                            const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                            const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                            const KinRecoObjects& kinRecoObjects,
                            const tth::GenObjectIndices& genObjectIndices, const tth::RecoObjectIndices& recoObjectIndices,
                            const tth::GenLevelWeights& genLevelWeights, const tth::RecoLevelWeights& recoLevelWeights,
                            const double& defaultWeight)const
{
    for(AnalysisHistogramsBase* analysisHistograms : v_analysisHistograms_){
        if(analysisHistograms) analysisHistograms->fill(recoObjects, commonGenObjects,
                                                   topGenObjects, higgsGenObjects,
                                                   kinRecoObjects,
                                                   recoObjectIndices, genObjectIndices,
                                                   genLevelWeights, recoLevelWeights,
                                                   defaultWeight, selectionStep);
    }

    if(mvaTreeHandler_) mvaTreeHandler_->fill(recoObjects, genObjectIndices, recoObjectIndices, defaultWeight, selectionStep);
}



void HiggsAnalysis::bookAll()
{
    for(AnalysisHistogramsBase* analysisHistograms : v_analysisHistograms_){
        if(analysisHistograms) analysisHistograms->book(fOutput);
    }
}



void HiggsAnalysis::clearAll()
{
    for(AnalysisHistogramsBase* analysisHistograms : v_analysisHistograms_){
        if(analysisHistograms) analysisHistograms->clear();
    }
}








