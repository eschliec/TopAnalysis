#define Analysis_cxx

#include "Analysis.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <TMath.h>
#include <TSystem.h>
#include <Math/VectorUtil.h>
#include <TLorentzVector.h>
#include <set>
#include <cmath>
#include <TString.h>

using namespace std;
using ROOT::Math::VectorUtil::DeltaPhi;
using ROOT::Math::VectorUtil::DeltaR;

double SampleXSection(TString sample){
    
    //MC cross sections taken from:
    //  https://twiki.cern.ch/twiki/bin/view/CMS/StandardModelCrossSectionsat8TeV
    //  AN-12/194    AN-12/228
    
    if(sample.Contains("data"))        {return 1;}
    if(sample.Contains("ttbar"))       {return 225.197;}
    if(sample.Contains("single"))      {return 22.2;}
    if(sample.Contains("ww"))          {return 54.838;}
    if(sample.Contains("wz"))          {return 33.21;}
    if(sample.Contains("zz"))          {return 17.654;}
    if(sample.Contains("1050"))        {return 860.5; } //5745.25;}
    if(sample.Contains("50inf"))       {return 3532.8; } //3503.71;}
    if(sample.Contains("wtolnu"))      {return 36257.2;}
    if(sample.Contains("wjets"))       {return 36257.2;}
    if(sample.Contains("qcdmu15"))     {return 3.640E8*3.7E-4;}
    if(sample.Contains("qcdmu2030"))   {return 2.870E8*6.500E-3;}
    if(sample.Contains("qcdmu3050"))   {return 6.609E7*12.20E-3;}
    if(sample.Contains("qcdmu5080"))   {return 8.802E6*21.80E-3;}
    if(sample.Contains("qcdmu80120"))  {return 1.024E6*39.50E-3;}
    if(sample.Contains("qcdmu120170")) {return 1.578E5*47.30E-3;}
    if(sample.Contains("qcdem2030"))   {return 2.886E8*10.10E-3;}
    if(sample.Contains("qcdem3080"))   {return 7.433E7*62.10E-3;}
    if(sample.Contains("qcdem80170"))  {return 1.191E6*153.9E-3;}
    if(sample.Contains("qcdbcem2030")) {return 2.886E8*5.800E-4;}
    if(sample.Contains("qcdbcem3080")) {return 7.424E7*2.250E-3;}
    if(sample.Contains("qcdbcem80170")){return 1.191E6*10.90E-3;}

    std::cerr << "No cross section found for sample: " << sample << std::endl;
    exit(2);
}

void Analysis::Begin ( TTree * )
{
    EventCounter = 0;
    bEff = 0;
    
    prepareTriggerSF();
    prepareLeptonIDSF();
    prepareBtagSF();
}

template<class T>
T* Analysis::store(T* obj)
{
    fOutput->Add(obj);
    return obj;
}

void Analysis::SlaveBegin ( TTree * )
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    h_step5 = store(new TH1D ( "step5", "event count at step 5", 10, 0, 10 ));
    h_step6 = store(new TH1D ( "step6", "event count at step 6", 10, 0, 10 ));
    h_step7 = store(new TH1D ( "step7", "event count at step 7", 10, 0, 10 ));
    h_step8 = store(new TH1D ( "step8", "event count at step 8", 10, 0, 10 ));
    h_step9 = store(new TH1D ( "step9", "event count at step 9", 10, 0, 10 ));

    //h_jetMultiAll = store(new TH1D ( "HypjetMultiAll", "Jet Multiplicity (AllJets)", 10, -0.5, 9.5 ));
    h_jetMultiXSec = store(new TH1D ( "HypjetMultiXSec", "Jet Multiplicity (for cross-section)", 10, -0.5, 9.5 ));
    h_jetMulti = store(new TH1D ( "HypjetMulti", "Jet Multiplicity", 10, -0.5, 9.5 ));
    h_jetMulti_diLep = store(new TH1D ( "HypjetMulti_diLep", "Jet Multiplicity (after dilepton)", 10, -0.5, 9.5 ));
    h_jetMultiNoPU = store(new TH1D ( "HypjetMultiNoPU", "Jet Multiplicity (No Pileup or lumi weight)", 10, -0.5, 9.5 ));
//     h_jetMultiVisTop = store(new TH1D ( "HypjetMultiVisTop", "Jet Multiplicity for Visible Top (No Pileup or lumi Weight)", 10, -0.5, 9.5 ));
    h_BjetMulti = store(new TH1D ( "HypBjetMulti", "B-Jet Multiplicity", 10, -0.5, 9.5 ));

    h_HypTTBarRapidity = store(new TH1D ( "HypTTBarRapidity", "Rapidity of TTbar System (HYP)", 100, -5, 5 ));
    h_HypTTBarpT = store(new TH1D ( "HypTTBarpT", "pT of TTbar System (HYP)", 500, 0, 500 ));
    h_HypTTBarMass = store(new TH1D ( "HypTTBarMass", "Mass of TTbar System (HYP)", 2000, 0, 2000 ));
    h_HypLLBarMass = store(new TH1D ( "HypLLBarMass", "Mass of LLbar System (HYP)", 500, 0, 1000 ));
    h_HypLLBarpT = store(new TH1D ( "HypLLBarpT", "pT of LLbar System (HYP)", 200, 0, 1000 ));

    h_VisGenTTBarMass = store(new TH1D ( "VisGenTTBarMass", "Mass of TTbar System(VisGEN)", 1200, 0, 1200 ));
    h_VisGenTTBarRapidity = store(new TH1D ( "VisGenTTBarRapidity", "Rapidity of TTbar System(VisGEN)", 100, -5, 5 ));
    h_VisGenTTBarpT = store(new TH1D ( "VisGenTTBarpT", "pT of TTbar System(VisGEN)", 1200, 0, 1200 ));
    h_VisGenTopRapidity = store(new TH1D ( "VisGenTopRapidity", "Rapidity of Top(VisGEN)", 100, -5, 5 ));
    h_VisGenAntiTopRapidity = store(new TH1D ( "VisGenAntiTopRapidity", "Rapidity of AntiTop(VisGEN)", 100, -5, 5 ));

    h_VisGenLLBarpT = store(new TH1D ( "VisGenLLBarpT", "pT of LLbar System(VisGEN)", 200, 0, 1000 ));
    h_VisGenLLBarMass = store(new TH1D ( "VisGenLLBarMass", "Mass of LLbar System(VisGEN)", 500, 0, 1000 ));

    h_RecoTTBarMass = store(new TH1D ( "RecoTTBarMass","Mass of TTbar System (HYP)",1200,0,1200 ));
    h_RecoTTBarRapidity = store(new TH1D ( "RecoTTBarRapidity","Rapidity of TTbar System (HYP)",100,-5,5 ));
    h_RecoTTBarpT = store(new TH1D ( "RecoTTBarpT","pT of TTbar System (HYP)",1200,0,1200 ));
    h_RecoToppT = store(new TH1D ( "RecoToppT","pT of Top (HYP)",1200,0,1200 ));
    h_RecoAntiToppT = store(new TH1D ( "RecoAntiToppT","pT of AntiTop (HYP)",1200,0,1200 ));
    h_RecoTopRapidity = store(new TH1D ( "RecoTopRapidity","Rapidity of Top (HYP)",100,-5,5 ));
    h_RecoAntiTopRapidity = store(new TH1D ( "RecoAntiTopRapidity","Rapidity of AntiTop (HYP)",100,-5,5 ));

    h_RecoBJetpT = store(new TH1D ( "RecoBJetpT","pT of BJet (HYP)",80,0,400 ));
    h_RecoAntiBJetpT = store(new TH1D ( "RecoAntiBJetpT","pT of AntiBJet (HYP)",80,0,400 ));
    h_RecoBJetRapidity = store(new TH1D ( "RecoBJetRapidity","Rapidity of BJet (HYP)",100,-5,5 ));
    h_RecoAntiBJetRapidity = store(new TH1D ( "RecoAntiBJetRapidity","Rapidity of AntiBJet (HYP)",100,-5,5 ));
    h_RecoBJetEta = store(new TH1D ( "RecoBJetEta","#eta of BJet (HYP)",100,-5,5 ));
    h_RecoAntiBJetEta = store(new TH1D ( "RecoAntiBJetEta","#eta of AntiBJet (HYP)",100,-5,5 ));

    h_RecoLLBarMass = store(new TH1D ( "RecoLLBarMass","Mass of LLbar System (HYP)",500,0,1000 ));
    h_RecoLLBarpT = store(new TH1D ( "RecoLLBarpT","pT of LLbar System (HYP)",200,0,1000 ));
    h_RecoLeptonpT = store(new TH1D ( "RecoLeptonpT","pT of Lepton (HYP)",240,0,1200 ));
    h_RecoAntiLeptonpT = store(new TH1D ( "RecoAntiLeptonpT","pT of AntiLepton (HYP)",240,0,1200 ));
    h_RecoLeptonEta = store(new TH1D ( "RecoLeptonEta","Eta of Lepton (HYP)",100,-5,5 ));
    h_RecoAntiLeptonEta = store(new TH1D ( "RecoAntiLeptonEta","Eta of AntiLepton (HYP)",100,-5,5 ));

    h_VisGenAll = store(new TH1D ( "VisGenAll", "All Visible Generated particles (IM)", 40, 0, 400 ));
    h_GenAll = store(new TH1D ( "GenAll", "AllGenerated particles (IM)", 40, 0, 400 ));
    Allh1 = store(new TH1D ( "Allh1", "DiLepton Mass", 40, 0, 400 ));
    h_diLepMassFull = store(new TH1D ( "DIMFull", "DiLepton Mass (Full Range)", 100, 0, 300 ));
    h_diLepMassFull_fullSel = store(new TH1D ( "DIMFull_fullSel", "DiLepton Mass (Full Range)", 100, 0, 300 ));
    Looseh1 = store(new TH1D ( "Looseh1", "DiLepton Mass", 40, 0, 400 ));
    Zh1 = store(new TH1D ( "Zh1", "DiLepton Mass in Z Window", 40, 0, 400 ));
    TTh1 = store(new TH1D ( "TTh1", "DiLepton Mass out of Z Window", 40, 0, 400 ));

    h_vertMulti = store(new TH1D ( "vertMulti", "Primary Vertex Multiplicity", 30, 0, 30 ));
    h_vertMulti_noPU = store(new TH1D ( "vertMulti_noPU", "Primary Vertex Multiplicity (no Pileup)", 30, 0, 30 ));
    h_MET = store(new TH1D ( "MET", "Missing Transverse Energy", 80, 0, 400 ));
    h_jetpT = store(new TH1D ( "jetpT", "jet pT", 80, 0, 400 ));
    h_jetHT = store(new TH1D ( "jetHT", "jet HT", 80, 0, 1000 ));

    h_MuonpT = store(new TH1D ( "MuonpT", "Muon pT (emu channel)", 80, 0, 400 ));
    h_MuonEta = store(new TH1D ( "MuonEta", "Muon Eta (emu channel)", 100, -5, 5 ));
    h_ElectronpT = store(new TH1D ( "ElectronpT", "Electron pT (emu channel)", 80, 0, 400 ));
    h_ElectronEta = store(new TH1D ( "ElectronEta", "Electron Eta (emu channel)", 100, -5, 5 ));

    h_LeptonpT = store(new TH1D ( "LeptonpT", "Lepton pT", 80, 0, 400 ));
    h_LeptonEta = store(new TH1D ( "LeptonEta", "Lepton Eta", 100, -5, 5 ));
    h_LeptonpT_diLep = store(new TH1D ( "LeptonpT_diLep", "Lepton pT (after dilepton cut)", 80, 0, 400 ));
    h_LeptonEta_diLep = store(new TH1D ( "LeptonEta_diLep", "Lepton Eta (after dilepton cut)", 100, -5, 5 ));

    h_AntiLeptonpT = store(new TH1D ( "AntiLeptonpT", "AntiLepton pT", 80, 0, 400 ));
    h_AntiLeptonEta = store(new TH1D ( "AntiLeptonEta", "AntiLepton Eta", 100, -5, 5 ));
    h_AntiLeptonpT_diLep = store(new TH1D ( "AntiLeptonpT_diLep", "Lepton pT (after dilepton cut)", 80, 0, 400 ));
    h_AntiLeptonEta_diLep = store(new TH1D ( "AntiLeptonEta_diLep", "Lepton Eta (after dilepton cut)", 100, -5, 5 ));

    h_HypToppT = store(new TH1D ( "HypToppT", "Top pT", 400, 0, 400 ));
    h_HypTopEta = store(new TH1D ( "HypTopEta", "Top pT", 100, -5, 5 ));
    h_HypTopMass = store(new TH1D ( "HypTopMass", "Top Mass", 80, 0, 400 ));
    h_HypTopRapidity = store(new TH1D ( "HypTopRapidity", "Top Rapidity", 100, -5, 5 ));

    h_HypAntiToppT = store(new TH1D ( "HypAntiToppT", "AntiTop pT", 400, 0, 400 ));
    h_HypAntiTopEta = store(new TH1D ( "HypAntiTopEta", "AntiTop pT", 100, -5, 5 ));
    h_HypAntiTopMass = store(new TH1D ( "HypAntiTopMass", "AntiTop Mass", 80, 0, 400 ));
    h_HypAntiTopRapidity = store(new TH1D ( "HypAntiTopRapidity", "Top Rapidity", 100, -5, 5 ));

    h_HypLeptonpT = store(new TH1D ( "HypLeptonpT", "Lepton Hypothesis pT", 80, 0, 400 ));
    h_HypLeptonEta = store(new TH1D ( "HypLeptonEta", "Lepton Eta", 100, -5, 5 ));

    h_HypAntiLeptonpT = store(new TH1D ( "HypAntiLeptonpT", "AntiLepton Hypothesis pT", 80, 0, 400 ));
    h_HypAntiLeptonEta = store(new TH1D ( "HypAntiLeptonEta", "AntiLepton Hypothesis Eta", 100, -5, 5 ));

    h_HypBJetpT = store(new TH1D ( "HypBJetpT", "B Hypothesis pT", 80, 0, 400 ));
    h_HypBJetEta = store(new TH1D ( "HypBJetEta", "B Hypothesis Eta", 100, -5, 5 ));
    h_HypBJetRapidity = store(new TH1D ( "HypBJetRapidity", "B Hypothesis Eta", 100, -5, 5 ));

    h_HypAntiBJetpT = store(new TH1D ( "HypAntiBJetpT", "AntiB Hypothesis pT", 80, 0, 400 ));
    h_HypAntiBJetEta = store(new TH1D ( "HypAntiBJetEta", "AntiB Hypothesis Eta", 100, -5, 5 ));
    h_HypAntiBJetRapidity = store(new TH1D ( "HypAntiBJetRapidity", "AntiB Hypothesis Eta", 100, -5, 5 ));

    h_VisGenToppT = store(new TH1D ( "VisGenToppT", "Top pT (VisGen)", 400, 0, 400 ));
    h_VisGenTopEta = store(new TH1D ( "VisGenTopEta", "Top Eta (VisGen)", 100, -5, 5 ));

    h_VisGenAntiToppT = store(new TH1D ( "VisGenAntiToppT", "AntiTop pT (VisGen)", 400, 0, 400 ));
    h_VisGenAntiTopEta = store(new TH1D ( "VisGenAntiTopEta", "AntiTop pT (VisGen)", 100, -5, 5 ));

    h_VisGenLeptonpT = store(new TH1D ( "VisGenLeptonpT", "Lepton VisGenothesis pT", 80, 0, 400 ));
    h_VisGenLeptonEta = store(new TH1D ( "VisGenLeptonEta", "Lepton Eta", 100, -5, 5 ));

    h_VisGenAntiLeptonpT = store(new TH1D ( "VisGenAntiLeptonpT", "AntiLepton VisGenothesis pT", 80, 0, 400 ));
    h_VisGenAntiLeptonEta = store(new TH1D ( "VisGenAntiLeptonEta", "AntiLepton VisGenothesis Eta", 100, -5, 5 ));

    h_VisGenBJetpT = store(new TH1D ( "VisGenBJetpT", "B VisGenothesis pT", 80, 0, 400 ));
    h_VisGenBJetEta = store(new TH1D ( "VisGenBJetEta", "B VisGenothesis Eta", 100, -5, 5 ));
    h_VisGenBJetRapidity = store(new TH1D ( "VisGenBJetRapidity", "B VisGenothesis Rapidity", 100, -5, 5 ));

    h_VisGenAntiBJetpT = store(new TH1D ( "VisGenAntiBJetpT", "AntiB VisGenothesis pT", 80, 0, 400 ));
    h_VisGenAntiBJetEta = store(new TH1D ( "VisGenAntiBJetEta", "AntiB VisGenothesis Eta", 100, -5, 5 ));
    h_VisGenAntiBJetRapidity = store(new TH1D ( "VisGenAntiBJetRapidity", "AntiB VisGenothesis Rapidity", 100, -5, 5 ));

    /*  h_VisGenBQuarkpT = store(new TH1D("VisGenBQuarkpT", "B Quark VisGenothesis pT", 80, 0, 400));
    h_VisGenBQuarkEta = store(new TH1D("VisGenBQuarkEta", "B Quark VisGenothesis Eta", 100, -5, 5));
    h_VisGenBQuarkRapidity = store(new TH1D("VisGenBQuarkRapidity", "B Quark VisGenothesis Rapidity", 100, -5, 5));

    h_VisGenAntiBQuarkpT = store(new TH1D("VisGenAntiBQuarkpT", "AntiB Quark VisGenothesis pT", 80, 0, 400));
    h_VisGenAntiBQuarkEta = store(new TH1D("VisGenAntiBQuarkEta", "AntiB Quark VisGenothesis Eta", 100, -5, 5));
    h_VisGenAntiBQuarkRapidity = store(new TH1D("VisGenAntiBQuarkRapidity", "AntiB Quark VisGenothesis Rapidity", 100, -5, 5));
    */
    /*h_GenToppT = store(new TH1D("GenToppT", "Top pT (Gen)", 80, 0, 400));
    h_GenTopEta = store(new TH1D("GenTopEta", "Top Eta (Gen)", 100, -5, 5));
    h_GenTopRapidity = store(new TH1D("GenTopRapidity", "Top Rapidity (Gen)", 100, -5, 5));

    h_GenAntiToppT = store(new TH1D("GenAntiToppT", "AntiTop pT (Gen)", 80, 0, 400));
    h_GenAntiTopEta = store(new TH1D("GenAntiTopEta", "AntiTop Eta (Gen)", 100, -5, 5));
    h_GenAntiTopRapidity = store(new TH1D("GenAntiTopRapidity", "AntiTop Rapidity (Gen)", 100, -5, 5));

    h_GenLeptonpT = store(new TH1D("GenLeptonpT", "Lepton Genothesis pT", 80, 0, 400));
    h_GenLeptonEta = store(new TH1D("GenLeptonEta", "Lepton Eta", 100, -5, 5));

    h_GenAntiLeptonpT = store(new TH1D("GenAntiLeptonpT", "AntiLepton Genothesis pT", 80, 0, 400));
    h_GenAntiLeptonEta = store(new TH1D("GenAntiLeptonEta", "AntiLepton Genothesis Eta", 100, -5, 5));

    h_GenBQuarkpT = store(new TH1D("GenBQuarkpT", "B Quark Genothesis pT", 80, 0, 400));
    h_GenBQuarkEta = store(new TH1D("GenBQuarkEta", "B Quark Genothesis Eta", 100, -5, 5));
    h_GenBQuarkRapidity = store(new TH1D("GenBQuarkRapidity", "B Quark Genothesis Rapidity", 100, -5, 5));

    h_GenAntiBQuarkpT = store(new TH1D("GenAntiBQuarkpT", "AntiB Quark Genothesis pT", 80, 0, 400));
    h_GenAntiBQuarkEta = store(new TH1D("GenAntiBQuarkEta", "AntiB Quark Genothesis Eta", 100, -5, 5));
    h_GenAntiBQuarkRapidity = store(new TH1D("GenAntiBQuarkRapidity", "AntiB Quark Genothesis Rapidity", 100, -5, 5));

    h_GenBJetpT = store(new TH1D("GenBJetpT", "B Genothesis pT", 80, 0, 400));
    h_GenBJetEta = store(new TH1D("GenBJetEta", "B Genothesis Eta", 100, -5, 5));
    h_GenBJetRapidity = store(new TH1D("GenBJetRapidity", "B Genothesis Rapidity", 100, -5, 5));

    h_GenAntiBJetpT = store(new TH1D("GenAntiBJetpT", "AntiB Genothesis pT", 80, 0, 400));
    h_GenAntiBJetEta = store(new TH1D("GenAntiBJetEta", "AntiB Genothesis Eta", 100, -5, 5));
    h_GenAntiBJetRapidity = store(new TH1D("GenAntiBJetRapidity", "Anti B Genothesis Rapidity", 100, -5, 5));
    */
    h_GenRecoBJetpT = store(new TH2D ( "GenRecoBJetpT", "Gen/Reco Matching", 80, 0, 400, 80, 0, 400 ));
    h_GenRecoBJetEta = store(new TH2D ( "GenRecoBJetEta", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 ));
    h_GenRecoBJetRapidity = store(new TH2D ( "GenRecoBJetRapidity", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 ));

    h_GenRecoAntiBJetpT = store(new TH2D ( "GenRecoAntiBJetpT", "Gen/Reco Matching", 80, 0, 400, 80, 0, 400 ));
    h_GenRecoAntiBJetEta = store(new TH2D ( "GenRecoAntiBJetEta", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 ));
    h_GenRecoAntiBJetRapidity = store(new TH2D ( "GenRecoAntiBJetRapidity", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 ));

    h_GenRecoLeptonEta = store(new TH2D ( "GenRecoLeptonEta", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 ));
    h_GenRecoAntiLeptonEta = store(new TH2D ( "GenRecoAntiLeptonEta", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 ));
    h_GenRecoLeptonpT = store(new TH2D ( "GenRecoLeptonpT", "Gen/Reco Matching", 80, 0, 400, 80, 0, 400 ));
    h_GenRecoAntiLeptonpT = store(new TH2D ( "GenRecoAntiLeptonpT", "Gen/Reco Matching", 80, 0, 400, 80, 0, 400 ));

    h_GenRecoTopRapidity = store(new TH2D ( "GenRecoTopRapidity", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 ));
    h_GenRecoAntiTopRapidity = store(new TH2D ( "GenRecoAntiTopRapidity", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 ));
    h_GenRecoToppT = store(new TH2D ( "GenRecoToppT", "Gen/Reco Matching", 400, 0, 400, 400, 0, 400 ));
    h_GenRecoAntiToppT = store(new TH2D ( "GenRecoAntiToppT", "Gen/Reco Matching", 400, 0, 400, 400, 0, 400 ));

    h_GenRecoTTBarRapidity = store(new TH2D ( "GenRecoTTBarRapidity", "Rapidity of TTbar System (HYP)", 100, -5, 5, 100, -5, 5 ));
    h_GenRecoTTBarpT = store(new TH2D ( "GenRecoTTBarpT", "pT of TTbar System (HYP)", 500, 0, 500, 500, 0, 500 ));
    h_GenRecoTTBarMass = store(new TH2D ( "GenRecoTTBarMass", "Mass of TTbar System (HYP)", 2000, 0, 2000, 2000, 0, 2000 ));
    h_GenRecoLLBarMass = store(new TH2D ( "GenRecoLLBarMass", "Mass of LLbar System (HYP)", 500, 0, 1000, 500, 0, 1000 ));
    h_GenRecoLLBarpT = store(new TH2D ( "GenRecoLLBarpT", "pT of LLbar System (HYP)", 200, 0, 1000, 200, 0, 1000 ));

    h_NJetMatching = store(new TH1D ( "NJetMatching", "NJet Gen/Reco Matching", 5, 0, 5 ));

    h_GenRecoLLBarDPhi = store(new TH2D ( "GenRecoLLBarDPhi", "Gen/Reco Matching", 100, 0., 3.2, 100, 0., 3.2 ));
    h_GenRecoLeptonantiBjetMass = store(new TH2D ( "GenRecoLeptonBjetMass", "Gen/Reco Matching", 500, 0, 1000, 500, 0, 1000 ));
    h_GenRecoAntiLeptonBjetMass = store(new TH2D ( "GenRecoAntiLeptonBjetMass", "Gen/Reco Matching", 500, 0, 1000, 500, 0, 1000 ));
    h_GenRecoJetMult = store(new TH2D ( "GenRecoJetMult", "Gen/REco Matching", 26, -0.5, 25.5, 26, -0.5, 25.5 ));

    h_HypLLBarDPhi = store(new TH1D ( "HypLLBarDPhi", "#Delta#phi(Lep, AntiLep) (HYP)",110, 0., 3.2 ));
    h_HypLeptonantiBjetMass = store(new TH1D ( "HypLeptonBjetMass", "Mass(Lep, AntiBJet) (HYP)", 500, 0, 1000 ));
    h_HypAntiLeptonBjetMass = store(new TH1D ( "HypAntiLeptonBjetMass", "Mass(AntiLep, BJet) (HYP)", 500, 0, 1000 ));
    h_HypJetMult = store(new TH1D ( "HypJetMult", "Jet Multiplicity (HYP)", 26, -0.5, 25.5 ));

    h_VisGenLLBarDPhi = store(new TH1D ( "VisGenLLBarDPhi", "#Delta #Phi (Lep, AntiLep) (VisGEN)", 100, 0., 3.2 ));
    h_VisGenLeptonantiBjetMass = store(new TH1D ( "VisGenLeptonBjetMass", "M(Lep, AntiBJet) (VisGEN)", 500, 0, 1000 ));
    h_VisGenAntiLeptonBjetMass = store(new TH1D ( "VisGenAntiLeptonBjetMass", "M(AntiLep, BJet) (VisGEN)", 500, 0, 1000 ));
    h_VisGenJetMult = store(new TH1D ( "VisGenJetMult", "Jet Multiplicty (VisGEN)", 26, -0.5, 25.5 ));

    h_RecoLLBarDPhi = store(new TH1D ( "RecoLLBarDPhi", "#Delta #Phi (Lep, AntiLep) (Reco)", 100, 0., 3.2 ));
    h_RecoLeptonantiBjetMass = store(new TH1D ( "RecoLeptonBjetMass", "M(Lep, AntiBJet) (Reco)", 500, 0, 1000 ));
    h_RecoAntiLeptonBjetMass = store(new TH1D ( "RecoAntiLeptonBjetMass", "M(AntiLep, BJet) (Reco)", 500, 0, 1000 ));
    h_RecoJetMult = store(new TH1D ( "RecoJetMult", "Jet Multiplicty (Reco)", 26, -0.5, 25.5 ));

    h_HypToppTLead = store(new TH1D ( "HypToppTLead","Leading pT Top pT",400,0,400 ));
    h_RecoToppTLead = store(new TH1D ( "RecoToppTLead","Leading pT Top pT",400,0,400 ));
    h_VisGenToppTLead = store(new TH1D ( "VisGenToppTLead","Leading pT Top pT",400,0,400 ));
    h_GenRecoToppTLead = store(new TH2D ( "GenRecoToppTLead", "Gen/Reco Matching", 400,0,400,400,0,400 ));

    h_HypToppTNLead = store(new TH1D ( "HypToppTNLead","NLeading pT Top pT",400,0,400 ));
    h_RecoToppTNLead = store(new TH1D ( "RecoToppTNLead","NLeading pT Top pT",400,0,400 ));
    h_VisGenToppTNLead = store(new TH1D ( "VisGenToppTNLead","NLeading pT Top pT",400,0,400 ));
    h_GenRecoToppTNLead = store(new TH2D ( "GenRecoToppTNLead", "Gen/Reco Matching", 400,0,400,400,0,400 ));

    h_HypTopRapidityLead = store(new TH1D ( "HypTopRapidityLead","Leading pT Top Rapidity",100,-5,5 ));
    h_RecoTopRapidityLead = store(new TH1D ( "RecoTopRapidityLead","Leading pT Top Rapidity",100,-5,5 ));
    h_VisGenTopRapidityLead = store(new TH1D ( "VisGenTopRapidityLead","Leading pT Top Rapidity",100,-5,5 ));
    h_GenRecoTopRapidityLead = store(new TH2D ( "GenRecoTopRapidityLead", "Gen/Reco Matching", 100,-5,5,100,-5,5 ));

    h_HypTopRapidityNLead = store(new TH1D ( "HypTopRapidityNLead","NLeading pT Top Rapidity",100,-5,5 ));
    h_RecoTopRapidityNLead = store(new TH1D ( "RecoTopRapidityNLead","NLeading pT Top Rapidity",100,-5,5 ));
    h_VisGenTopRapidityNLead = store(new TH1D ( "VisGenTopRapidityNLead","NLeading pT Top Rapidity",100,-5,5 ));
    h_GenRecoTopRapidityNLead = store(new TH2D ( "GenRecoTopRapidityNLead", "Gen/Reco Matching", 100,-5,5,100,-5,5 ));

    h_HypTopMassLead = store(new TH1D ( "HypTopMassLead","Leading pT Top Mass",80,0,400 ));
    h_RecoTopMassLead = store(new TH1D ( "RecoTopMassLead","Leading pT Top Mass",80,0,400 ));
    h_VisGenTopMassLead = store(new TH1D ( "VisGenTopMassLead","Leading pT Top Mass",80,0,400 ));
    h_GenRecoTopMassLead = store(new TH2D ( "GenRecoTopMassLead", "Gen/Reco Matching", 80,0,400,80,0,400 ));

    h_HypTopMassNLead = store(new TH1D ( "HypTopMassNLead","NLeading pT Top Mass",80,0,400 ));
    h_RecoTopMassNLead = store(new TH1D ( "RecoTopMassNLead","NLeading pT Top Mass",80,0,400 ));
    h_VisGenTopMassNLead = store(new TH1D ( "VisGenTopMassNLead","NLeading pT Top Mass",80,0,400 ));
    h_GenRecoTopMassNLead = store(new TH2D ( "GenRecoTopMassNLead", "Gen/Reco Matching", 80,0,400,80,0,400 ));


    h_HypLeptonpTLead = store(new TH1D ( "HypLeptonpTLead","Leading pT Lepton pT",400,0,400 ));
    h_RecoLeptonpTLead = store(new TH1D ( "RecoLeptonpTLead","Leading pT Lepton pT",400,0,400 ));
    h_VisGenLeptonpTLead = store(new TH1D ( "VisGenLeptonpTLead","Leading pT Lepton pT",400,0,400 ));
    h_GenRecoLeptonpTLead = store(new TH2D ( "GenRecoLeptonpTLead", "Gen/Reco Matching", 400,0,400,400,0,400 ));

    h_HypLeptonpTNLead = store(new TH1D ( "HypLeptonpTNLead","NLeading pT Lepton pT",400,0,400 ));
    h_RecoLeptonpTNLead = store(new TH1D ( "RecoLeptonpTNLead","NLeading pT Lepton pT",400,0,400 ));
    h_VisGenLeptonpTNLead = store(new TH1D ( "VisGenLeptonpTNLead","NLeading pT Lepton pT",400,0,400 ));
    h_GenRecoLeptonpTNLead = store(new TH2D ( "GenRecoLeptonpTNLead", "Gen/Reco Matching", 400,0,400,400,0,400 ));

    h_HypLeptonEtaLead = store(new TH1D ( "HypLeptonEtaLead","Leading pT Lepton Eta",100,-5,5 ));
    h_RecoLeptonEtaLead = store(new TH1D ( "RecoLeptonEtaLead","Leading pT Lepton Eta",100,-5,5 ));
    h_VisGenLeptonEtaLead = store(new TH1D ( "VisGenLeptonEtaLead","Leading pT Lepton Eta",100,-5,5 ));
    h_GenRecoLeptonEtaLead = store(new TH2D ( "GenRecoLeptonEtaLead", "Gen/Reco Matching", 100,-5,5,100,-5,5 ));

    h_HypLeptonEtaNLead = store(new TH1D ( "HypLeptonEtaNLead","NLeading pT Lepton Eta",100,-5,5 ));
    h_RecoLeptonEtaNLead = store(new TH1D ( "RecoLeptonEtaNLead","NLeading pT Lepton Eta",100,-5,5 ));
    h_VisGenLeptonEtaNLead = store(new TH1D ( "VisGenLeptonEtaNLead","NLeading pT Lepton Eta",100,-5,5 ));
    h_GenRecoLeptonEtaNLead = store(new TH2D ( "GenRecoLeptonEtaNLead", "Gen/Reco Matching", 100,-5,5,100,-5,5 ));

    h_HypBJetpTLead = store(new TH1D ( "HypBJetpTLead","Leading pT BJet pT",400,0,400 ));
    h_RecoBJetpTLead = store(new TH1D ( "RecoBJetpTLead","Leading pT BJet pT",400,0,400 ));
    h_VisGenBJetpTLead = store(new TH1D ( "VisGenBJetpTLead","Leading pT BJet pT",400,0,400 ));
    h_GenRecoBJetpTLead = store(new TH2D ( "GenRecoBJetpTLead", "Gen/Reco Matching", 400,0,400,400,0,400 ));

    h_HypBJetpTNLead = store(new TH1D ( "HypBJetpTNLead","NLeading pT BJet pT",400,0,400 ));
    h_RecoBJetpTNLead = store(new TH1D ( "RecoBJetpTNLead","NLeading pT BJet pT",400,0,400 ));
    h_VisGenBJetpTNLead = store(new TH1D ( "VisGenBJetpTNLead","NLeading pT BJet pT",400,0,400 ));
    h_GenRecoBJetpTNLead = store(new TH2D ( "GenRecoBJetpTNLead", "Gen/Reco Matching", 400,0,400,400,0,400 ));

    h_HypBJetEtaLead = store(new TH1D ( "HypBJetEtaLead","Leading pT BJet Eta",100,-5,5 ));
    h_RecoBJetEtaLead = store(new TH1D ( "RecoBJetEtaLead","Leading pT BJet Eta",100,-5,5 ));
    h_VisGenBJetEtaLead = store(new TH1D ( "VisGenBJetEtaLead","Leading pT BJet Eta",100,-5,5 ));
    h_GenRecoBJetEtaLead = store(new TH2D ( "GenRecoBJetEtaLead", "Gen/Reco Matching", 100,-5,5,100,-5,5 ));

    h_HypBJetEtaNLead = store(new TH1D ( "HypBJetEtaNLead","NLeading pT BJet Eta",100,-5,5 ));
    h_RecoBJetEtaNLead = store(new TH1D ( "RecoBJetEtaNLead","NLeading pT BJet Eta",100,-5,5 ));
    h_VisGenBJetEtaNLead = store(new TH1D ( "VisGenBJetEtaNLead","NLeading pT BJet Eta",100,-5,5 ));
    h_GenRecoBJetEtaNLead = store(new TH2D ( "GenRecoBJetEtaNLead", "Gen/Reco Matching", 100,-5,5,100,-5,5 ));

    /*
    //New plots from Carmen: Begin
    h_RecoLeadingJetpT   = store(new TH1D("RecoLeadingJetpT","pT of all Jets (HYP)",80,0,400));
    h_RecoLeadingJetEta  = store(new TH1D("RecoLeadingJetEta","#eta of all Jets (HYP)",100,-5,5));
    h_RecoNLeadingJetpT  = store(new TH1D("RecoNLeadingJetpT","pT of all Jets (HYP)",80,0,400));
    h_RecoNLeadingJetEta = store(new TH1D("RecoNLeadingJetEta","#eta of all Jets (HYP)",100,-5,5));

    h_GenRecoLeadingJetpT   = store(new TH2D("GenRecoLeadingJetpT","pT of all Jets (HYP)",80,0,400,80,0,400));
    h_GenRecoLeadingJetEta  = store(new TH2D("GenRecoLeadingJetEta","#eta of all Jets (HYP)",100,-5,5,100,-5,5));
    h_GenRecoNLeadingJetpT  = store(new TH2D("GenRecoNLeadingJetpT","pT of all Jets (HYP)",80,0,400,80,0,400));
    h_GenRecoNLeadingJetEta = store(new TH2D("GenRecoNLeadingJetEta","#eta of all Jets (HYP)",100,-5,5,100,-5,5));

    h_VisGenLeadingJetpT   = store(new TH1D("VisGenLeadingJetpT","pT of leading Jets (VisGEN)",80,0,400));
    h_VisGenLeadingJetEta  = store(new TH1D("VisGenLeadingJetEta","#eta of leading Jets (VisGEN)",100,-5,5));
    h_VisGenNLeadingJetpT  = store(new TH1D("VisGenNLeadingJetpT","pT of leading Jets (VisGEN)",80,0,400));
    h_VisGenNLeadingJetEta = store(new TH1D("VisGenNLeadingJetEta","#eta of leading Jets (VisGEN)",100,-5,5));

    h_HypLeadingJetpT   = store(new TH1D("HypLeadingJetpT","pT of leading Jets (HYP)",80,0,400));
    h_HypLeadingJetEta  = store(new TH1D("HypLeadingJetEta","#eta of leading Jets (HYP)",100,-5,5));
    h_HypNLeadingJetpT  = store(new TH1D("HypNLeadingJetpT","pT of leading Jets (HYP)",80,0,400));
    h_HypNLeadingJetEta = store(new TH1D("HypNLeadingJetEta","#eta of leading Jets (HYP)",100,-5,5));

    h_RecoExtraJetpT   = store(new TH1D("RecoExtraJetpT","pT of additional Jet (HYP)",80,0,400));
    h_RecoExtraJetEta  = store(new TH1D("RecoExtraJetEta","#eta of additional Jet (HYP)",100,-5,5));
    h_RecoExtraJetpT2  = store(new TH1D("RecoExtraJetpT2","pT of additional Jet (HYP)",80,0,400));
    h_RecoExtraJetEta2 = store(new TH1D("RecoExtraJetEta2","#eta of additional Jet (HYP)",100,-5,5));
    h_RecoExtraJetpT3  = store(new TH1D("RecoExtraJetpT3","pT of additional Jet (HYP)",80,0,400));
    h_RecoExtraJetEta3 = store(new TH1D("RecoExtraJetEta3","#eta of additional Jet (HYP)",100,-5,5));
    h_RecoExtraJetpT4  = store(new TH1D("RecoExtraJetpT4","pT of additional Jet (HYP)",80,0,400));
    h_RecoExtraJetEta4 = store(new TH1D("RecoExtraJetEta4","#eta of additional Jet (HYP)",100,-5,5));

    h_HypExtraJetpT   = store(new TH1D("HypExtraJetpT","pT of additional Jet",80,0,400));
    h_HypExtraJetEta  = store(new TH1D("HypExtraJetEta","#eta of additional Jet",100,-5,5));
    h_HypExtraJetpT2  = store(new TH1D("HypExtraJetpT2","pT of additional Jet",80,0,400));
    h_HypExtraJetEta2 = store(new TH1D("HypExtraJetEta2","#eta of additional Jet",100,-5,5));
    h_HypExtraJetpT3  = store(new TH1D("HypExtraJetpT3","pT of additional Jet",80,0,400));
    h_HypExtraJetEta3 = store(new TH1D("HypExtraJetEta3","#eta of additional Jet",100,-5,5));
    h_HypExtraJetpT4  = store(new TH1D("HypExtraJetpT4","pT of additional Jet",80,0,400));
    h_HypExtraJetEta4 = store(new TH1D("HypExtraJetEta4","#eta of additional Jet",100,-5,5));

    h_VisGenExtraJetpT   = store(new TH1D("VisGenExtraJetpT","pT of gen additional Jet",80,0,400));
    h_VisGenExtraJetEta  = store(new TH1D("VisGenExtraJetEta","#eta of gen additional Jet",100,-5,5));
    h_VisGenExtraJetpT2  = store(new TH1D("VisGenExtraJetpT2","pT of gen additional Jet",80,0,400));
    h_VisGenExtraJetEta2 = store(new TH1D("VisGenExtraJetEta2","#eta of gen additional Jet",100,-5,5));
    h_VisGenExtraJetpT3  = store(new TH1D("VisGenExtraJetpT3","pT of gen additional Jet",80,0,400));
    h_VisGenExtraJetEta3 = store(new TH1D("VisGenExtraJetEta3","#eta of gen additional Jet",100,-5,5));
    h_VisGenExtraJetpT4  = store(new TH1D("VisGenExtraJetpT4","pT of gen additional Jet",80,0,400));
    h_VisGenExtraJetEta4 = store(new TH1D("VisGenExtraJetEta4","#eta of gen additional Jet",100,-5,5));

    h_GenRecoExtraJetpT   = store(new TH2D("GenRecoExtraJetpT","Gen/Reco pT of additional Jet",80,0,400,80,0,400));
    h_GenRecoExtraJetEta  = store(new TH2D("GenRecoExtraJetEta","Gen/Reco #eta of additional Jet",100,-5,5,100,-5,5));
    h_GenRecoExtraJetpT2  = store(new TH2D("GenRecoExtraJetpT2","Gen/Reco pT of additional Jet",80,0,400,80,0,400));
    h_GenRecoExtraJetEta2 = store(new TH2D("GenRecoExtraJetEta2","Gen/Reco #eta of additional Jet",100,-5,5,100,-5,5));
    h_GenRecoExtraJetpT3  = store(new TH2D("GenRecoExtraJetpT3","Gen/Reco pT of additional Jet",80,0,400,80,0,400));
    h_GenRecoExtraJetEta3 = store(new TH2D("GenRecoExtraJetEta3","Gen/Reco #eta of additional Jet",100,-5,5,100,-5,5));
    h_GenRecoExtraJetpT4  = store(new TH2D("GenRecoExtraJetpT4","Gen/Reco pT of additional Jet",80,0,400,80,0,400));
    h_GenRecoExtraJetEta4 = store(new TH2D("GenRecoExtraJetEta4","Gen/Reco #eta of additional Jet",100,-5,5,100,-5,5));

    h_GenRecoJetMultpt40 = store(new TH2D("GenRecoJetMultpt40", "Gen/Reco Matching",10,-0.5,9.5,10,-0.5,9.5));
    h_RecoJetMultpt40    = store(new TH1D("RecoJetMultpt40", "Jet Multiplicity (HYP)",10,-0.5,9.5));
    h_HypJetMultpt40     = store(new TH1D("HypJetMultpt40", "Jet Multiplicity (HYP)",10,-0.5,9.5));
    h_GenJetMultpt40     = store(new TH1D("GenJetMultpt40", "Jet Multiplicty (GEN)",10,-0.5,9.5));
    h_VisGenJetMultpt40  = store(new TH1D("VisGenJetMultpt40", "Jet Multiplicty (VisGEN)",10,-0.5,9.5));

    h_GenRecoJetMultpt60 = store(new TH2D("GenRecoJetMultpt60", "Gen/Reco Matching",10,-0.5,9.5,10,-0.5,9.5));
    h_RecoJetMultpt60    = store(new TH1D("RecoJetMultpt60", "Jet Multiplicity (HYP)",10,-0.5,9.5));
    h_HypJetMultpt60     = store(new TH1D("HypJetMultpt60", "Jet Multiplicity (HYP)",10,-0.5,9.5));
    h_GenJetMultpt60     = store(new TH1D("GenJetMultpt60", "Jet Multiplicty (GEN)",10,-0.5,9.5));
    h_VisGenJetMultpt60  = store(new TH1D("VisGenJetMultpt60", "Jet Multiplicty (VisGEN)",10,-0.5,9.5));
    //New plots from Carmen: End
*/

    //btagSF
    const int PtMax = 17;
    const int EtaMax = 6;
    Double_t ptbins[PtMax+1] = {0.,30.,40.,50.,60.,70.,80.,100.,120.,160.,210.,260.,320.,400.,500.,670.,1000.,2000.};
    Double_t etabins[EtaMax+1] = {0.0,0.5,1.0,1.5,2.0,2.4,3.0};
    h_bjets = store(new TH2D("bjets2D", "unTagged Bjets", PtMax, ptbins, EtaMax, etabins));              h_bjets->Sumw2();
    h_btaggedjets = store(new TH2D("bjetsTagged2D", "Tagged Bjets", PtMax, ptbins, EtaMax, etabins));    h_btaggedjets->Sumw2();
    h_cjets = store(new TH2D("cjets2D", "unTagged Cjets", PtMax, ptbins, EtaMax, etabins));              h_cjets->Sumw2();
    h_ctaggedjets = store(new TH2D("cjetsTagged2D", "Tagged Cjets", PtMax, ptbins, EtaMax, etabins));    h_ctaggedjets->Sumw2();
    h_ljets = store(new TH2D("ljets2D", "unTagged Ljets", PtMax, ptbins, EtaMax, etabins));              h_ljets->Sumw2();
    h_ltaggedjets = store(new TH2D("ljetsTagged2D", "Tagged Ljets", PtMax, ptbins, EtaMax, etabins));    h_ltaggedjets->Sumw2();
    
    h_BTagSF = store(new TH1D ( "BTagSF", "BTagging SF per event", 100 , 0.95, 1.05 ));
    h_BTagSF->Sumw2();
}

double Analysis::Median(TH1 * h1)
{ 
   int n = h1->GetXaxis()->GetNbins();
   std::vector<double>  x(n);
   h1->GetXaxis()->GetCenter( &x[0] );
   TH1D* h1D = dynamic_cast<TH1D*>(h1);
   if (!h1D) { cerr << "Median needs a TH1D!\n"; exit(7); }
   const double * y = h1D->GetArray(); 
   // exclude underflow/overflows from bin content array y
   return TMath::Median(n, &x[0], &y[1]); 
}


Bool_t Analysis::Process ( Long64_t entry )
{
    static const double JETPTCUT = 30;
    
    if ( ++EventCounter % 100000 == 0 ) cout << "Event Counter: " << EventCounter << endl;
    GetRecoBranches(entry);
    
    if (isTtbarPlusTauSample) {
        bool isViaTau = decayMode > 40 || ( decayMode % 10 > 4 );
        if (runViaTau != isViaTau) return kTRUE;
    }
    
    //We must correct for the madGraph branching fraction being 1/9 for dileptons (PDG average is .108)
    if ( correctMadgraphBR ) {
        if ( decayMode == 11 ) { //all hadronic decay
            weightGenerator *= (0.676*1.5) * (0.676*1.5);
        } else if ( decayMode < 20 || ( decayMode % 10 == 1) ) { //semileptonic Decay
            weightGenerator *= (0.108*9) * (0.676*1.5);
        } else {//dileptonic decay (including taus!)
            weightGenerator *= (0.108*9) * (0.108*9);
        }
    }
    
    double weightPU = 1;
    if (isMC) { 
        //still have lumi weights for old plotterclass
//         weightGenerator *= lumiWeight;
        
        if (pureweighter) {
            weightPU = pureweighter->getPUweight(vertMultiTrue);
        }
    }
    
    int BHadronIndex=-1;
    int AntiBHadronIndex=-1;

    if ( isSignal ) {
        GetSignalBranches ( entry );

        std::vector<int> idx_leadbHadJet;
        std::vector<int> idx_nleadbHadJet;
        //To avoid recopying may code lines, we select HERE the BHadron JET Indices to cut on.

        //time to choose which genJet we actually want

        idx_leadbHadJet.insert ( idx_leadbHadJet.begin(), 4, -1 );
        idx_nleadbHadJet.insert ( idx_nleadbHadJet.begin(), 4, -1 );
        /*
          idx_bHadJet will have 4 jet indices
          [0] is the highest pT jet with a B-Hadron
          [1] is the highest pT jet with a B-Hadron also matched to a top quark
          [2] highest pT jet of those matched closest (in DeltaR) to the B-Hadron
          [3] highest pT jet of those matched closest (in DeltaR) to the B-Hadron also matched to a top quark
        */

        bool LeadBHadhighpTjet = false;
        bool LeadBHadhighpTjetfromtop = false;
        bool NLeadBHadhighpTjet = false;
        bool NLeadBHadhighpTjetfromtop = false;

        int hadron_index = -1;
        int antihadron_index = -1;
        int hadrontop_index = -1;
        int antihadrontop_index = -1;
        
        //Case 1: highest pT genJet matched to a BHadron
        //need to remove jets from the genJetCollection which are below the JETPTCUT
        //while (allGenJets->size() > 0 && allGenJets->back().Pt() < JETPTCUT) allGenJets->pop_back();
        //while (jet->size() > 0 && jet->back().Pt() < JETPTCUT) jet->pop_back();
        
        for ( size_t genJet = 0; 
              genJet < allGenJets->size() && allGenJets->at(genJet).pt() >= JETPTCUT; 
              ++genJet ) 
        {
            for ( size_t bHadron=0; bHadron < BHadrons->size(); bHadron++ ) {
                if ( (*BHadronVsJet)[genJet*BHadrons->size()+bHadron]==1 
                     && (!LeadBHadhighpTjet || !LeadBHadhighpTjetfromtop || !NLeadBHadhighpTjet || !NLeadBHadhighpTjetfromtop) )
                {
                    if ( LeadBHadhighpTjet==false ) {
                        idx_leadbHadJet[0] = genJet;
                        LeadBHadhighpTjet = true;
                        hadron_index = bHadron;
                        if ( ( *BHadronFromTopB ) [bHadron] == true ) {
                            idx_leadbHadJet[1] = genJet;
                            LeadBHadhighpTjetfromtop = true;
                            hadrontop_index = bHadron;
                        }
                    } else if ( LeadBHadhighpTjetfromtop == false ) {
                        if ( ( *BHadronFromTopB ) [bHadron] == true ) {
                            idx_leadbHadJet[1] = genJet;
                            LeadBHadhighpTjetfromtop = true;
                            hadrontop_index = bHadron;
                        }
                    } else if ( NLeadBHadhighpTjet==false && bHadron!=hadron_index && idx_leadbHadJet[0] != genJet ) {
                        idx_nleadbHadJet[0] = genJet;
                        NLeadBHadhighpTjet = true;
                        if ( ( *BHadronFromTopB ) [bHadron] == true && bHadron!=hadrontop_index && idx_leadbHadJet[1] != genJet ) {
                            idx_nleadbHadJet[1] = genJet;
                            NLeadBHadhighpTjetfromtop = true;
                        }
                    } else if ( NLeadBHadhighpTjetfromtop == false && bHadron!=hadrontop_index && idx_leadbHadJet[1] != genJet ) {
                        if ( ( *BHadronFromTopB ) [bHadron] == true ) {
                            idx_nleadbHadJet[1] = genJet;
                            LeadBHadhighpTjetfromtop = true;
                        }
                    }//series of if statements to find highest pT jet
                }
            }
            for ( size_t antibHadron=0; antibHadron < AntiBHadrons->size(); antibHadron++ ) {
                if ( (*AntiBHadronVsJet)[genJet*AntiBHadrons->size()+antibHadron]==1 
                    && ( LeadBHadhighpTjet ==false || LeadBHadhighpTjetfromtop == false || NLeadBHadhighpTjet ==false || NLeadBHadhighpTjetfromtop == false ) && idx_leadbHadJet[0] != genJet ) 
                {
                    if ( LeadBHadhighpTjet==false ) {
                        idx_leadbHadJet[0] = genJet;
                        LeadBHadhighpTjet = true;
                        antihadron_index = antibHadron;
                        if ( ( *AntiBHadronFromTopB ) [antibHadron] == true ) {
                            idx_leadbHadJet[1] = genJet;
                            LeadBHadhighpTjetfromtop = true;
                            antihadrontop_index = antibHadron;
                        }
                    } else if ( LeadBHadhighpTjetfromtop == false ) {
                        if ( ( *AntiBHadronFromTopB ) [antibHadron] == true ) {
                            idx_leadbHadJet[1] = genJet;
                            LeadBHadhighpTjetfromtop = true;
                            antihadrontop_index = antibHadron;
                        }
                    } else if ( NLeadBHadhighpTjet==false && antibHadron!=antihadron_index && idx_leadbHadJet[0] != genJet ) {
                        idx_nleadbHadJet[0] = genJet;
                        NLeadBHadhighpTjet = true;
                        if ( ( *AntiBHadronFromTopB ) [antibHadron] == true && antibHadron!=antihadrontop_index && idx_leadbHadJet[1] != genJet ) {
                            idx_nleadbHadJet[1] = genJet;
                            NLeadBHadhighpTjetfromtop = true;
                        }
                    } else if ( NLeadBHadhighpTjetfromtop == false && antibHadron!=antihadrontop_index && idx_leadbHadJet[1] != genJet ) {
                        if ( ( *AntiBHadronFromTopB ) [antibHadron] == true ) {
                            idx_nleadbHadJet[1] = genJet;
                            LeadBHadhighpTjetfromtop = true;
                        }
                    }
                }
            }
        }
 
        //Case 2: highest pT genJets matched closest to a BHadron
        //BHadJetIndex: vector containing the GetJet indices matched, in DeltaR, to a BHadron. Starting from the highest pT jet.
        if ( BHadJetIndex->size() != 0 ) idx_leadbHadJet[2] = ( *BHadJetIndex ) [0];
        for ( size_t i=0; i < BHadJetIndex->size(); ++i ) {
            //Only search for those jets matched in DeltaR with a BHadron
            for ( size_t j=0; j<BHadrons->size() ; ++j ) {
                if ( ( *BHadronVsJet ) [ ( ( *BHadJetIndex ) [i] ) * BHadrons->size()+j] == 1 && ( *BHadronFromTopB ) [j] == true ) {
                    idx_leadbHadJet[3] = ( *BHadJetIndex ) [i];
                }
            }
        }

        //AntiBHadJetIndex: vector containing the GetJet indices matched, in DeltaR, to a AntiBHadron. Starting from the highest pT jet.
        if ( AntiBHadJetIndex->size() != 0 ) idx_nleadbHadJet[2] = ( *AntiBHadJetIndex ) [0];
        for ( size_t i=0; i < AntiBHadJetIndex->size(); ++i ) {
            //Only search for those jets matched in DeltaR with a AntiBHadron
            for ( size_t j=0; j < AntiBHadrons->size() ; ++j ) {
                //if ((*AntiBHadronVsJet)[i*AntiBHadrons_+j] == 1 && (*AntiBHadronFromTopB)[j] == true) {idx_antibHadJet[3] = (*AntiBHadJetIndex)[i];}
                if ( ( *AntiBHadronVsJet ) [ ( ( *AntiBHadJetIndex ) [i] ) * AntiBHadrons->size()+j] == 1 && ( *AntiBHadronFromTopB ) [j] == true ) {
                    idx_nleadbHadJet[3] = ( *AntiBHadJetIndex ) [i];
                }
            }
        }


//     //To avoid recopying many code lines, we select HERE the BHadron JET Indices to cut on.
//     int BHadronIndex;
//     int AntiBHadronIndex;
        //Case 1A: highest pT genJet matched to a BHadron
        BHadronIndex = idx_leadbHadJet[0];
        AntiBHadronIndex = idx_nleadbHadJet[0];
//
        //   //Case 1B: highest pT genJet matched to a BHadron from Top
        //BHadronIndex = idx_bHadJet[1];
        //AntiBHadronIndex = idx_antibHadJet[1];

        //   //Case 2A: highest pT genJets matched closest to a BHadron
        //BHadronIndex = idx_bHadJet[2];
        //AntiBHadronIndex = idx_antibHadJet[2];
        //
        //   //Case 2B: highest pT genJets matched closest to a BHadron from Top
        //    BHadronIndex = idx_leadbHadJet[3];
        //AntiBHadronIndex = idx_nleadbHadJet[3];


    }

    
    double BtagWP = 0.244; //CSV Loose working point
    vector<int> BJetIndex;
    for ( vector<double>::iterator it = jetBTagCSV->begin(); it<jetBTagCSV->end(); it++ ) {
        if ( *it > BtagWP ) {
            //BJetIndex.push_back ( *it );
            BJetIndex.push_back((it-jetBTagCSV->begin())); //change asked by Tyler
        }
    }



    //Should we just look for two Bjets above 0.244 or the two highest bjets?:: Make this a function
    bool hasSolution = false;
    int solutionIndex = 0;
    for ( size_t i =0; i < HypTop->size(); ++i ) {
        if (jet->at((*HypJet0index)[i]).pt() < JETPTCUT || jet->at((*HypJet1index)[i]).pt() < JETPTCUT) continue;
        bool jet0tagged = jetBTagCSV->at( (*HypJet0index)[i] ) > BtagWP;
        bool jet1tagged = jetBTagCSV->at( (*HypJet1index)[i] ) > BtagWP;
        if (jet0tagged && jet1tagged) {   
            //solution with 2 tags found, so take it and stop
            solutionIndex = i;
            hasSolution = true;
            break;
        }
        if (!hasSolution && (jet0tagged || jet1tagged)) {
            //one btag found, save solution if it is the first one
            solutionIndex = i;
            hasSolution = true;
        }
    }
    
    
    LV LeadGenTop, NLeadGenTop;
    LV LeadGenLepton, NLeadGenLepton;
    LV LeadGenBJet, NLeadGenBJet;
    
    if ( isSignal ) {
    
        double trueLevelWeight = weightGenerator * weightPU;
        h_GenAll->Fill(GenTop->M(), trueLevelWeight);
        if ( GenLepton->pt() > 20 && GenAntiLepton->pt() > 20 
              && abs( GenLepton->eta() ) < 2.4 && abs ( GenAntiLepton->eta() ) < 2.4 ) {
            //if (LVGenBQuark.Pt()>JETPTCUT && LVGenAntiBQuark.Pt()>JETPTCUT && abs(LVGenBQuark.Eta())<2.4 && abs(LVGenAntiBQuark.Eta())<2.4){
            if ( BHadronIndex != -1 && AntiBHadronIndex != -1 ) {
                if ( allGenJets->at(BHadronIndex).pt() > JETPTCUT && 
                    abs ( allGenJets->at(BHadronIndex).eta() ) < 2.4 &&
                    allGenJets->at(AntiBHadronIndex).pt() > JETPTCUT &&
                    abs ( allGenJets->at(AntiBHadronIndex).Eta() ) < 2.4 )
                {

                    h_VisGenAll->Fill(GenTop->M(), trueLevelWeight);

                    h_VisGenLLBarpT->Fill(( *GenLepton + *GenAntiLepton ).Pt(), trueLevelWeight );
                    h_VisGenLLBarMass->Fill(( *GenLepton + *GenAntiLepton ).M(), trueLevelWeight );

                    h_VisGenLeptonpT->Fill(GenLepton->pt(), trueLevelWeight );
                    h_VisGenAntiLeptonpT->Fill(GenAntiLepton->Pt(), trueLevelWeight );

                    h_VisGenLeptonEta->Fill(GenLepton->Eta(), trueLevelWeight );
                    h_VisGenAntiLeptonEta->Fill(GenAntiLepton->Eta(), trueLevelWeight );

                    h_VisGenBJetEta->Fill(allGenJets->at(BHadronIndex).Eta(), trueLevelWeight );
                    h_VisGenAntiBJetEta->Fill(allGenJets->at(AntiBHadronIndex).Eta(), trueLevelWeight );
                    h_VisGenBJetRapidity->Fill(allGenJets->at(BHadronIndex).Rapidity(), trueLevelWeight );
                    h_VisGenAntiBJetRapidity->Fill(allGenJets->at(AntiBHadronIndex).Rapidity(), trueLevelWeight );
                    h_VisGenBJetpT->Fill(allGenJets->at(BHadronIndex).Pt(), trueLevelWeight );
                    h_VisGenAntiBJetpT->Fill(allGenJets->at(AntiBHadronIndex).Pt(), trueLevelWeight );

                    h_VisGenLLBarDPhi->Fill(abs( DeltaPhi(*GenLepton, *GenAntiLepton)), trueLevelWeight );
                    h_VisGenLeptonantiBjetMass->Fill(( *GenLepton + allGenJets->at(AntiBHadronIndex) ).M(), trueLevelWeight );
                    h_VisGenAntiLeptonBjetMass->Fill(( *GenAntiLepton + allGenJets->at(BHadronIndex) ).M(), trueLevelWeight );
                    h_VisGenJetMult->Fill(allGenJets->size(), trueLevelWeight );

                    //Begin: Select & Fill histograms with Leading pT and 2nd Leading pT: Lepton and BJet
                    if(GenLepton->Pt()>GenAntiLepton->Pt()){
                        LeadGenLepton  = *GenLepton;
                        NLeadGenLepton = *GenAntiLepton;
                    }
                    else{
                        LeadGenLepton  = *GenAntiLepton;
                        NLeadGenLepton = *GenLepton;

                    }
                    
                    if(allGenJets->at(BHadronIndex).Pt() > allGenJets->at(AntiBHadronIndex).Pt()){
                        LeadGenBJet  = (*allGenJets).at(BHadronIndex);
                        NLeadGenBJet = (*allGenJets).at(AntiBHadronIndex);
                    }
                    else{
                        LeadGenBJet  = (*allGenJets).at(AntiBHadronIndex);
                        NLeadGenBJet = (*allGenJets).at(BHadronIndex);
                    }
                    
                    h_VisGenLeptonpTLead->Fill(LeadGenLepton.Pt(), trueLevelWeight);
                    h_VisGenLeptonpTNLead->Fill(NLeadGenLepton.Pt(), trueLevelWeight);
                    h_VisGenLeptonEtaLead->Fill(LeadGenLepton.Eta(), trueLevelWeight);
                    h_VisGenLeptonEtaNLead->Fill(NLeadGenLepton.Eta(), trueLevelWeight);
                    
                    h_VisGenBJetpTLead->Fill(LeadGenBJet.Pt(), trueLevelWeight);
                    h_VisGenBJetpTNLead->Fill(NLeadGenBJet.Pt(), trueLevelWeight);
                    h_VisGenBJetEtaLead->Fill(LeadGenBJet.Eta(), trueLevelWeight);
                    h_VisGenBJetEtaNLead->Fill(NLeadGenBJet.Eta(), trueLevelWeight);
                    //End: Select & Fill histograms with Leading pT and 2nd Leading pT: Lepton and BJet
                    
//                     //New plots from Carmen: Begin
//                     bool firstJet = 0, secondJet = 0;
//                     for(int k=0; k<allGenJets->size(); k++){
//                         if(abs(allGenJets->at(k).Eta())>2.4 || allGenJets->at(k).Pt()< 30.0) {continue;}
//                         if(!firstJet) {
//                             h_VisGenLeadingJetpT->Fill(allGenJets->at(k).Pt(),trueLevelWeight);
//                             h_VisGenLeadingJetEta->Fill(allGenJets->at(k).Eta(),trueLevelWeight);
//                             firstJet=1;
//                             continue;
//                         }
//                         if(firstJet && !secondJet){
//                             h_VisGenNLeadingJetpT->Fill(allGenJets->at(k).Pt(),trueLevelWeight);
//                             h_VisGenNLeadingJetEta->Fill(allGenJets->at(k).Eta(),trueLevelWeight);
//                             secondJet=1;
//                             break;
//                         }
//                     }
//                     
//                     for(int genJet=0; genJet<allGenJets->size(); genJet++){
//                         if(abs(allGenJets->at(genJet).Eta() ) > 2.4 || TMath::Abs(DeltaR(*GenLepton, allGenJets->at(genJet))) < 0.4 
//                             || TMath::Abs(DeltaR(*GenAntiLepton, allGenJets->at(genJet))) < 0.4 ) {
//                             continue;
//                         }
//                         if(allGenJets->at(genJet).Pt()> 30) {
//                             GetJets_cut++; 
//                             if(allGenJets->at(BHadronIndex) != allGenJets->at(genJet) && allGenJets->at(AntiBHadronIndex) != allGenJets->at(genJet)) { 
//                                 jetHTGen+=allGenJets->at(genJet).Pt(); 
//                                 if(jetnum < 3) {
//                                     jetnum++;
//                                     extragenjet[jetnum] = genJet;
//                                 }
//                             }
//                             if(allGenJets->at(genJet).Pt()> 40) GetJets_cut40++;
//                             if(allGenJets->at(genJet).Pt()> 60) GetJets_cut60++;
//                         }
//                     }//for
// 
//                     h_VisGenJetMult->Fill(GetJets_cut,trueLevelWeight);
//                     h_VisGenJetMultpt40->Fill(GetJets_cut40,trueLevelWeight);
//                     h_VisGenJetMultpt60->Fill(GetJets_cut60,trueLevelWeight);
//                     if(jetnum>2){
//                         h_VisGenExtraJetpT4->Fill(allGenJets->at(extragenjet[3]).Pt(),trueLevelWeight);
//                         h_VisGenExtraJetEta4->Fill(allGenJets->at(extragenjet[3]).Eta(),trueLevelWeight);
//                     }
//                     else if(jetnum>1){
//                         h_VisGenExtraJetpT3->Fill(allGenJets->at(extragenjet[2]).Pt(),trueLevelWeight);
//                         h_VisGenExtraJetEta3->Fill(allGenJets->at(extragenjet[2]).Eta(),trueLevelWeight);
//                     }
//                     else if(jetnum>0){
//                         h_VisGenExtraJetpT2->Fill(allGenJets->at(extragenjet[1]).Pt(),trueLevelWeight);
//                         h_VisGenExtraJetEta2->Fill(allGenJets->at(extragenjet[1]).Eta(),trueLevelWeight);
//                     }
//                     else if(jetnum == 0){
//                         h_VisGenExtraJetpT->Fill(allGenJets->at(extragenjet[0]).Pt(),trueLevelWeight);
//                         h_VisGenExtraJetEta->Fill(allGenJets->at(extragenjet[0]).Eta(),trueLevelWeight);
//                     }
//                     //New plots from Carmen: End
                }
            }
        }
        
        LV genttbar(*GenTop + *GenAntiTop);
        h_VisGenTTBarMass->Fill(genttbar.M(), trueLevelWeight );
        h_VisGenTTBarRapidity->Fill(genttbar.Rapidity(), trueLevelWeight );
        h_VisGenTTBarpT->Fill(genttbar.Pt(), trueLevelWeight );

        h_VisGenToppT->Fill(GenTop->Pt(), trueLevelWeight );
        h_VisGenAntiToppT->Fill(GenAntiTop->Pt(), trueLevelWeight );
        h_VisGenTopRapidity->Fill(GenTop->Rapidity(), trueLevelWeight );
        h_VisGenAntiTopRapidity->Fill(GenAntiTop->Rapidity(), trueLevelWeight );
        h_VisGenTopEta->Fill(GenTop->Eta(), trueLevelWeight );
        h_VisGenAntiTopEta->Fill(GenAntiTop->Eta(), trueLevelWeight );
        
        //Begin: Fill histograms with Leading pT and 2nd Leading pT: Top
        if(GenTop->Pt()>GenAntiTop->Pt()){
            LeadGenTop  = *GenTop;
            NLeadGenTop = *GenAntiTop;
        }
        else{
            LeadGenTop  = *GenAntiTop;
            NLeadGenTop = *GenTop;
        }
        h_VisGenToppTLead->Fill(LeadGenTop.Pt(), trueLevelWeight);
        h_VisGenToppTNLead->Fill(NLeadGenTop.Pt(), trueLevelWeight);
        h_VisGenTopRapidityLead->Fill(LeadGenTop.Rapidity(), trueLevelWeight);
        h_VisGenTopRapidityNLead->Fill(NLeadGenTop.Rapidity(), trueLevelWeight);
        h_VisGenTopMassLead->Fill(LeadGenTop.M(), trueLevelWeight);
        h_VisGenTopMassNLead->Fill(NLeadGenTop.M(), trueLevelWeight);
        //End: Fill histograms with Leading pT and 2nd Leading pT: Top
        
    }//for visible top events

    //===CUT===
    // check if event was triggered (only needed for the signal sample which does not
    // contain trigger preselection cuts)
    if (isTtbarPlusTauSample) {
        if (!(((triggerBits & 0x0000FF) && channel == "mumu")    //mumu triggers in rightmost byte
           || ((triggerBits & 0x00FF00) && channel == "emu")     //emu in 2nd byte
           || ((triggerBits & 0xFF0000) && channel == "ee")))    //ee in 3rd byte
        {
            return kTRUE;
        }
    }

    size_t LeadLeptonNumber = 0;
    size_t NLeadLeptonNumber = 0;
    bool hasLeptonPair = getLeptonPair(LeadLeptonNumber, NLeadLeptonNumber);

    //===CUT===
    // we need an OS lepton pair
    if (! hasLeptonPair) return kTRUE;
    
    LV dilepton = lepton->at(LeadLeptonNumber) + lepton->at(NLeadLeptonNumber);
    
    //===CUT===
    //with at least 12 GeV invariant mass
    if (dilepton.M() < 12) return kTRUE;
    
    // find l+ and l-
    LV leptonPlus;
    LV leptonMinus;
    if (lepQ->at(LeadLeptonNumber) == +1) {
        leptonPlus = lepton->at(LeadLeptonNumber);
        leptonMinus = lepton->at(NLeadLeptonNumber);
    } else {
        leptonMinus = lepton->at(LeadLeptonNumber);
        leptonPlus = lepton->at(NLeadLeptonNumber);
    }

    //Now determine the lepton trigger and ID scale factors
    double weightLepSF = isMC ? getLeptonIDSF(leptonPlus, leptonMinus) : 1;
    double weightTrigSF = isMC ? getTriggerSF(leptonPlus, leptonMinus) : 1;
    
    //First control plots after dilepton selection (without Z cut)
    double weight = weightGenerator*weightTrigSF*weightLepSF;
    //weight even without PU reweighting
    h_vertMulti_noPU->Fill(vertMulti, weight);
    
    //apply PU reweighting - continue with control plots
    weight *= weightPU;
    h_vertMulti->Fill(vertMulti, weight);
    
    h_jetMulti_diLep->Fill(jet->size(), weight);
    h_diLepMassFull->Fill(dilepton.M(), weight);

    
    //****************************************
    //handle inverted Z cut
    // Fill loose dilepton mass histogram before any jet cuts
    bool isZregion = dilepton.M() > 76 && dilepton.M() < 106;
    bool hasJets = jet->size() > 1 && jet->at(1).Pt() > JETPTCUT;
    bool hasMetOrEmu = channel == "emu" || met->Et() > 30;
    bool hasBtag = BJetIndex.size() > 0;
    double weightKinFit = 1;
    double weightBtagSF = -1; //trick: initialize to -1 to avoid calculation of the btagSF twice
    
    if ( isZregion ) {
        Looseh1->Fill(dilepton.M(), weight);
        if ( hasJets && hasMetOrEmu && hasBtag && hasSolution) {
            weightBtagSF = isMC ? calculateBtagSF() : 1;
            double fullWeights = weightGenerator*weightPU*weightLepSF*weightBtagSF*weightTrigSF*weightKinFit;
            Zh1->Fill(dilepton.M(), fullWeights);
            Allh1->Fill(dilepton.M(), fullWeights);
        }
    }
    
    //=== CUT ===
    //Exclude the Z window
    if (channel != "emu" && isZregion) return kTRUE;
    
    h_step5->Fill(1, weight);
    h_LeptonpT_diLep->Fill(leptonMinus.Pt(), weight);
    h_AntiLeptonpT_diLep->Fill(leptonPlus.Pt(), weight);
    h_LeptonEta_diLep->Fill(leptonMinus.Eta(), weight);
    h_AntiLeptonEta_diLep->Fill(leptonPlus.Eta(), weight);
    
    //=== CUT ===
    //Require at least two jets > 30 GeV (check for > 30 needed because we might have 20 GeV jets in our NTuple)
    if (! hasJets) return kTRUE;
    h_step6->Fill(1, weight);
    
    //=== CUT ===
    //Require MET > 30 GeV in non-emu channels
    if (!hasMetOrEmu) return kTRUE;
    h_step7->Fill(1, weight);
 
    //=== CUT ===
    //Require at least one b tagged jet
    if (!hasBtag) return kTRUE;

    if (weightBtagSF == -1) weightBtagSF = isMC ? calculateBtagSF() : 1; //avoid calculation of the btagSF twice
    weight *= weightBtagSF;
    h_BTagSF->Fill(weightBtagSF );
    h_step8->Fill(1, weight );

    h_BjetMulti->Fill(BJetIndex.size(), weight);
    h_jetMulti->Fill(jet->size(), weight);
    
    //for HT, count only >= 30 GeV jets
    double jetHT = getJetHT(*jet, JETPTCUT);
    h_jetHT->Fill(jetHT, weight);
    for ( size_t i = 0; i < 2; ++i ) {
        h_jetpT->Fill(jet->at(i).Pt(), weight);
    }

    h_LeptonpT->Fill(leptonMinus.Pt(), weight);
    h_AntiLeptonpT->Fill(leptonPlus.Pt(), weight);
    h_LeptonEta->Fill(leptonMinus.Eta(), weight);
    h_AntiLeptonEta->Fill(leptonPlus.Eta(), weight);

    h_MET->Fill(met->Et(), weight);

    //loop over both leptons
    for (auto i : {LeadLeptonNumber, NLeadLeptonNumber}) {
        if ( lepType->at(i) == LEP_TYPE_ELECTRON ) {
            h_ElectronpT->Fill(lepton->at(i).Pt(), weight);
            h_ElectronEta->Fill(lepton->at(i).Eta(), weight);
        }
        if ( lepType->at(i) == LEP_TYPE_MUON ) {
            h_MuonpT->Fill(lepton->at(i).Pt(), weight);
            h_MuonEta->Fill(lepton->at(i).Eta(), weight);
        }
    }

    //=== CUT ===
    //Require at least one solution for the kinematic event reconstruction
    if (!hasSolution) return kTRUE;

    weight *= weightKinFit;
    h_step9->Fill(1, weight);
    h_jetMultiXSec->Fill(jet->size(), weight);
    h_jetMultiNoPU->Fill(jet->size(), weight / weightPU );
    h_diLepMassFull_fullSel->Fill(dilepton.M(), weight);
    
    //create helper variables
    
    //Begin: find 1st (and 2nd) leading pT particles: Top, Lepton, BJetIndex
    LV LeadHypTop, NLeadHypTop;
    LV LeadHypLepton, NLeadHypLepton;
    LV LeadHypBJet, NLeadHypBJet;
    
    if(HypTop->at(solutionIndex).Pt() > HypAntiTop->at(solutionIndex).Pt()){
        LeadHypTop  = HypTop->at(solutionIndex);
        NLeadHypTop = HypAntiTop->at(solutionIndex);
    }
    else{
        LeadHypTop  = HypAntiTop->at(solutionIndex);
        NLeadHypTop = HypTop->at(solutionIndex);
    }
    
    if(HypLepton->at(solutionIndex).Pt() > HypAntiLepton->at(solutionIndex).Pt()){
        LeadHypLepton  = HypLepton->at(solutionIndex);
        NLeadHypLepton = HypAntiLepton->at(solutionIndex);
    }
    else{
        LeadHypLepton  = HypAntiLepton->at(solutionIndex);
        NLeadHypLepton = HypLepton->at(solutionIndex);
    }
    
    if(HypBJet->at(solutionIndex).Pt() > HypAntiBJet->at(solutionIndex).Pt()){
        LeadHypBJet  = HypBJet->at(solutionIndex);
        NLeadHypBJet = HypAntiBJet->at(solutionIndex);
    }
    else{
        LeadHypBJet  = HypAntiBJet->at(solutionIndex);
        NLeadHypBJet = HypBJet->at(solutionIndex);
    }
    //End: find 1st (and 2nd) leading pT particles: Top, Lepton, BJetIndex
    
    //create ll and tt system
    LV hypllbar(HypLepton->at(solutionIndex) + HypAntiLepton->at(solutionIndex));
    LV hypttbar(HypTop->at(solutionIndex)+HypAntiTop->at(solutionIndex));
    
    //First fill the reco histograms (which have no scaling factors applied)
    double recoWeight = weightGenerator * weightPU;
    h_RecoTTBarMass->Fill(hypttbar.M(), recoWeight);
    h_RecoTTBarRapidity->Fill(hypttbar.Rapidity(), recoWeight);
    h_RecoTTBarpT->Fill(hypttbar.Pt(), recoWeight);
    h_RecoToppT->Fill(HypTop->at(solutionIndex).Pt(), recoWeight);
    h_RecoAntiToppT->Fill(HypAntiTop->at(solutionIndex).Pt(), recoWeight);
    h_RecoTopRapidity->Fill(HypTop->at(solutionIndex).Rapidity(), recoWeight);
    h_RecoAntiTopRapidity->Fill(HypAntiTop->at(solutionIndex).Rapidity(), recoWeight);

    h_RecoLLBarMass->Fill(hypllbar.M(), recoWeight);
    h_RecoLLBarpT->Fill(hypllbar.Pt(), recoWeight);
    h_RecoLeptonpT->Fill(HypLepton->at(solutionIndex).Pt(), recoWeight);
    h_RecoAntiLeptonpT->Fill(HypAntiLepton->at(solutionIndex).Pt(), recoWeight);
    h_RecoLeptonEta->Fill(HypLepton->at(solutionIndex).Eta(), recoWeight);
    h_RecoAntiLeptonEta->Fill(HypAntiLepton->at(solutionIndex).Eta(), recoWeight);

    h_RecoBJetpT->Fill(HypBJet->at(solutionIndex).Pt(), recoWeight);
    h_RecoAntiBJetpT->Fill(HypAntiBJet->at(solutionIndex).Pt(), recoWeight);
    h_RecoBJetRapidity->Fill(HypBJet->at(solutionIndex).Rapidity(), recoWeight);
    h_RecoAntiBJetRapidity->Fill(HypAntiBJet->at(solutionIndex).Rapidity(), recoWeight);
    h_RecoBJetEta->Fill(HypBJet->at(solutionIndex).Eta(), recoWeight);
    h_RecoAntiBJetEta->Fill(HypAntiBJet->at(solutionIndex).Eta(), recoWeight);

    h_RecoLLBarDPhi->Fill(abs ( DeltaPhi ( HypLepton->at(solutionIndex), HypAntiLepton->at(solutionIndex) ) ), recoWeight);
    h_RecoLeptonantiBjetMass->Fill(( HypLepton->at(solutionIndex)+HypAntiBJet->at(solutionIndex) ).M(), recoWeight);
    h_RecoAntiLeptonBjetMass->Fill(( HypAntiLepton->at(solutionIndex)+HypBJet->at(solutionIndex) ).M(), recoWeight);

//     //New plots from Carmen: Begin
//     h_RecoJetMult->Fill ( jet->size(),recoWeight);
//     h_RecoJetMultpt40->Fill(RecoJets_cut40,recoWeight);
//     h_RecoJetMultpt60->Fill(RecoJets_cut60,recoWeight);
//     h_HypJetMult->Fill ( jet->size(),weight);
//     h_HypJetMultpt40->Fill(RecoJets_cut40,weight);
//     h_HypJetMultpt60->Fill(RecoJets_cut60,weight);
// 
//     bool firstJet = 0, secondJet = 0;
//     for(int k=0; k<(int)jet->size(); k++){
//         if(abs(jet->at(k).Eta())>2.4 && jet->at(k).Pt()< 30.0) {continue;}
//         if(!firstJet) {
//             h_RecoLeadingJetpT->Fill(jet->at(k).Pt(),recoWeight);
//             h_RecoLeadingJetEta->Fill(jet->at(k).Eta(),recoWeight);
//             h_HypLeadingJetpT->Fill(jet->at(k).Pt(),weight);
//             h_HypLeadingJetEta->Fill(jet->at(k).Eta(),weight);
//             firstJet=1;
//             continue;
//         }
//         if(firstJet && !secondJet){
//             h_RecoNLeadingJetpT->Fill(jet->at(k).Pt(),recoWeight);
//             h_RecoNLeadingJetEta->Fill(jet->at(k).Eta(),recoWeight);
//             h_HypNLeadingJetpT->Fill(jet->at(k).Pt(),weight);
//             h_HypNLeadingJetEta->Fill(jet->at(k).Eta(),weight);
//             secondJet=1;
//             break;
//         }
//     }
//     
//     if(jetnumReco>2){
//         h_RecoExtraJetpT4->Fill(jet->at(extrarecojet[3]).Pt(),recoWeight);
//         h_RecoExtraJetEta4->Fill(jet->at(extrarecojet[3]).Eta(),recoWeight);
//         h_HypExtraJetpT4->Fill(jet->at(extrarecojet[3]).Pt(),weight);
//         h_HypExtraJetEta4->Fill(jet->at(extrarecojet[3]).Eta(),weight);
//     }
//     else if (jetnumReco>1){
//         h_RecoExtraJetpT3->Fill(jet->at(extrarecojet[2]).Pt(),recoWeight);
//         h_RecoExtraJetEta3->Fill(jet->at(extrarecojet[2]).Eta(),recoWeight);
//         h_HypExtraJetpT3->Fill(jet->at(extrarecojet[2]).Pt(),weight);
//         h_HypExtraJetEta3->Fill(jet->at(extrarecojet[2]).Eta(),weight);
//     }
//     else if (jetnumReco>0){
//         h_RecoExtraJetpT2->Fill(jet->at(extrarecojet[1]).Pt(),recoWeight);
//         h_RecoExtraJetEta2->Fill(jet->at(extrarecojet[1]).Eta(),recoWeight);
//         h_HypExtraJetpT2->Fill(jet->at(extrarecojet[1]).Pt(),weight);
//         h_HypExtraJetEta2->Fill(jet->at(extrarecojet[1]).Eta(),weight);
//     }
//     else if (jetnumReco >-1){
//         h_RecoExtraJetpT->Fill(jet->at(extrarecojet[0]).Pt(),recoWeight);
//         h_RecoExtraJetEta->Fill(jet->at(extrarecojet[0]).Eta(),recoWeight);
//         h_HypExtraJetpT->Fill(jet->at(extrarecojet[0]).Pt(),weight);
//         h_HypExtraJetEta->Fill(jet->at(extrarecojet[0]).Eta(),weight);
//     }
//     //New plots from Carmen: End
    
    h_RecoToppTLead->Fill(LeadHypTop.Pt(), recoWeight);
    h_RecoToppTNLead->Fill(NLeadHypTop.Pt(), recoWeight);
    h_RecoTopRapidityLead->Fill(LeadHypTop.Rapidity(), recoWeight);
    h_RecoTopRapidityNLead->Fill(NLeadHypTop.Rapidity(), recoWeight);
    h_RecoTopMassLead->Fill(LeadHypTop.M(), recoWeight);
    h_RecoTopMassNLead->Fill(NLeadHypTop.M(), recoWeight);
    
    h_RecoLeptonpTLead->Fill(LeadHypLepton.Pt(), recoWeight);
    h_RecoLeptonpTNLead->Fill(NLeadHypLepton.Pt(), recoWeight);
    h_RecoLeptonEtaLead->Fill(LeadHypLepton.Eta(), recoWeight);
    h_RecoLeptonEtaNLead->Fill(NLeadHypLepton.Eta(), recoWeight);
    
    h_RecoBJetpTLead->Fill(LeadHypBJet.Pt(), recoWeight);
    h_RecoBJetpTNLead->Fill(NLeadHypBJet.Pt(), recoWeight);
    h_RecoBJetEtaLead->Fill(LeadHypBJet.Eta(), recoWeight);
    h_RecoBJetEtaNLead->Fill(NLeadHypBJet.Eta(), recoWeight);
    
    //now go to the plots 
    h_HypTTBarMass->Fill(hypttbar.M(), weight);
    h_HypTTBarRapidity->Fill(hypttbar.Rapidity(), weight);
    h_HypTTBarpT->Fill(hypttbar.Pt(), weight);

    h_HypLLBarMass->Fill(hypllbar.M(), weight);
    h_HypLLBarpT->Fill(hypllbar.Pt(), weight);

    h_HypTopMass->Fill(HypTop->at(solutionIndex).M(), weight);
    h_HypAntiTopMass->Fill(HypAntiTop->at(solutionIndex).M(), weight);
    h_HypToppT->Fill(HypTop->at(solutionIndex).Pt(), weight);
    h_HypAntiToppT->Fill(HypAntiTop->at(solutionIndex).Pt(), weight);
    h_HypLeptonpT->Fill(HypLepton->at(solutionIndex).Pt(), weight);
    h_HypAntiLeptonpT->Fill(HypAntiLepton->at(solutionIndex).Pt(), weight);

    h_HypBJetpT->Fill(HypBJet->at(solutionIndex).Pt(), weight);
    h_HypAntiBJetpT->Fill(HypAntiBJet->at(solutionIndex).Pt(), weight);
    h_HypBJetRapidity->Fill(HypBJet->at(solutionIndex).Rapidity(), weight);
    h_HypAntiBJetRapidity->Fill(HypAntiBJet->at(solutionIndex).Rapidity(), weight);

    h_HypTopRapidity->Fill(HypTop->at(solutionIndex).Rapidity(), weight);
    h_HypAntiTopRapidity->Fill(HypAntiTop->at(solutionIndex).Rapidity(), weight);

    h_HypTopEta->Fill(HypTop->at(solutionIndex).Eta(), weight);
    h_HypAntiTopEta->Fill(HypAntiTop->at(solutionIndex).Eta(), weight);
    h_HypBJetEta->Fill(HypBJet->at(solutionIndex).Eta(), weight);
    h_HypAntiBJetEta->Fill(HypAntiBJet->at(solutionIndex).Eta(), weight);
    h_HypLeptonEta->Fill(HypLepton->at(solutionIndex).Eta(), weight);

    h_HypAntiLeptonEta->Fill(HypAntiLepton->at(solutionIndex).Eta(), weight);

    h_HypLLBarDPhi->Fill(abs ( DeltaPhi ( HypLepton->at(solutionIndex), HypAntiLepton->at(solutionIndex) ) ), weight);
    h_HypLeptonantiBjetMass->Fill(( HypLepton->at(solutionIndex) + HypAntiBJet->at(solutionIndex) ).M(), weight);
    h_HypAntiLeptonBjetMass->Fill(( HypAntiLepton->at(solutionIndex) + HypBJet->at(solutionIndex) ).M(), weight);

    h_HypToppTLead->Fill(LeadHypTop.Pt(), weight);
    h_HypToppTNLead->Fill(NLeadHypTop.Pt(), weight);
    h_HypTopRapidityLead->Fill(LeadHypTop.Rapidity(), weight);
    h_HypTopRapidityNLead->Fill(NLeadHypTop.Rapidity(), weight);
    h_HypTopMassLead->Fill(LeadHypTop.M(), weight);
    h_HypTopMassNLead->Fill(NLeadHypTop.M(), weight);

    h_HypLeptonpTLead->Fill(LeadHypLepton.Pt(), weight);
    h_HypLeptonpTNLead->Fill(NLeadHypLepton.Pt(), weight);
    h_HypLeptonEtaLead->Fill(LeadHypLepton.Eta(), weight);
    h_HypLeptonEtaNLead->Fill(NLeadHypLepton.Eta(), weight);
    
    h_HypBJetpTLead->Fill(LeadHypBJet.Pt(), weight);
    h_HypBJetpTNLead->Fill(NLeadHypBJet.Pt(), weight);
    h_HypBJetEtaLead->Fill(LeadHypBJet.Eta(), weight);
    h_HypBJetEtaNLead->Fill(NLeadHypBJet.Eta(), weight);
    
    if (!isZregion) { //also apply Z cut in emu!
        TTh1->Fill(dilepton.M(), weight);
        Allh1->Fill(dilepton.M(), weight);  //this is also filled in the Z region in the code above
    }

    //=== CUT ===
    //Following histograms only filled for the signal sample
    if (! isSignal) return kTRUE;

//     if ( GenLepton->Pt() > 20 && GenAntiLepton->Pt() > 20 
//          && abs(GenLepton->Eta()) < 2.4 && abs(GenAntiLepton->Eta()) < 2.4 ) 
//     {
//         //Comment the next 2 lines and uncomment the 3rd one for gen-level Vis PS cuts
//         //if (LVGenBQuark.Pt()>30 && LVGenAntiBQuark.Pt()>30 && abs(LVGenBQuark.Eta())<2.4 && abs(LVGenAntiBQuark.Eta())<2.4){
//         if ( BHadronIndex != -1 && allGenJets->at(BHadronIndex).Pt() > 30 
//             && abs( allGenJets->at(BHadronIndex).Eta() ) < 2.4 
//             && AntiBHadronIndex != -1 && allGenJets->at(AntiBHadronIndex).Pt() > 30 
//             && abs( allGenJets->at(AntiBHadronIndex).Eta() ) < 2.4 ) 
//         {
//             //            if(LVBHadronGenJet.Pt()>30 && abs(LVBHadronGenJet.Eta())<2.4 &&
//             //          LVAntiBHadronGenJet.Pt()>30 && abs(LVAntiBHadronGenJet.Eta())<2.4){
// 
//             h_jetMultiVisTop->Fill(jet->size(), weightLepSF*btagSF*trigEFF*weightKinFit );
//             //!! why these weights? , ah because label of the plot
//         }
//     }
    
    h_GenRecoLeptonEta->Fill(HypLepton->at(solutionIndex).Eta(), GenLepton->Eta(), weight );

    h_GenRecoAntiLeptonEta->Fill(HypAntiLepton->at(solutionIndex).Eta(), GenAntiLepton->Eta(), weight );
    h_GenRecoLeptonpT->Fill(HypLepton->at(solutionIndex).Pt(), GenLepton->Pt(), weight );
    h_GenRecoAntiLeptonpT->Fill(HypAntiLepton->at(solutionIndex).Pt(), GenAntiLepton->Pt(), weight );

    h_GenRecoTopRapidity->Fill(HypTop->at(solutionIndex).Rapidity(), GenTop->Rapidity(), weight );
    h_GenRecoAntiTopRapidity->Fill(HypAntiTop->at(solutionIndex).Rapidity(), GenAntiTop->Rapidity(), weight );
    h_GenRecoToppT->Fill(HypTop->at(solutionIndex).Pt(), GenTop->Pt(), weight );
    h_GenRecoAntiToppT->Fill(HypAntiTop->at(solutionIndex).Pt(), GenAntiTop->Pt(), weight );

    h_GenRecoLLBarDPhi->Fill(
        abs( DeltaPhi( HypLepton->at(solutionIndex), HypAntiLepton->at(solutionIndex) ) ), 
        abs( DeltaPhi( *GenLepton, *GenAntiLepton ) ), 
        weight );

    //Begin: find (and fill) the GenReco plots for Leading and NLeading Top, Lepton, BJet
    h_GenRecoToppTLead->Fill(LeadHypTop.Pt(), LeadGenTop.Pt(), weight);
    h_GenRecoToppTNLead->Fill(NLeadHypTop.Pt(), NLeadGenTop.Pt(), weight);
    h_GenRecoTopRapidityLead->Fill(LeadHypTop.Rapidity(), LeadGenTop.Rapidity(), weight);
    h_GenRecoTopRapidityNLead->Fill(NLeadHypTop.Rapidity(), NLeadGenTop.Rapidity(), weight);
    h_GenRecoTopMassLead->Fill(LeadHypTop.M(), LeadGenTop.M(), weight);
    h_GenRecoTopMassNLead->Fill(NLeadHypTop.M(), NLeadGenTop.M(), weight);

    h_GenRecoLeptonpTLead->Fill(LeadHypLepton.Pt(), LeadGenLepton.Pt(), weight);
    h_GenRecoLeptonpTNLead->Fill(NLeadHypLepton.Pt(), NLeadGenLepton.Pt(), weight);
    h_GenRecoLeptonEtaLead->Fill(LeadHypLepton.Eta(), LeadGenLepton.Eta(), weight);
    h_GenRecoLeptonEtaNLead->Fill(NLeadHypLepton.Eta(), NLeadGenLepton.Eta(), weight);
    
    if(BHadJetIndex >= 0 && AntiBHadJetIndex >= 0){
        h_GenRecoBJetpTLead->Fill(LeadHypBJet.Pt(), LeadGenBJet.Pt(), weight);
        h_GenRecoBJetpTNLead->Fill(NLeadHypBJet.Pt(), NLeadGenBJet.Pt(), weight);
        h_GenRecoBJetEtaLead->Fill(LeadHypBJet.Eta(), LeadGenBJet.Eta(), weight);
        h_GenRecoBJetEtaNLead->Fill(NLeadHypBJet.Eta(), NLeadGenBJet.Eta(), weight);
    }
    else{
        h_GenRecoBJetpTLead->Fill(LeadHypBJet.Pt(), -1000, weight);
        h_GenRecoBJetpTNLead->Fill(NLeadHypBJet.Pt(), -1000, weight);
        h_GenRecoBJetEtaLead->Fill(LeadHypBJet.Eta(), -1000, weight);
        h_GenRecoBJetEtaNLead->Fill(NLeadHypBJet.Eta(), -1000, weight);
    }
    //End: find (and fill) the GenReco plots for Leading and NLeading Leptons

    if ( BHadronIndex>=0 ) {
        h_GenRecoBJetpT->Fill(HypBJet->at(solutionIndex).Pt(), allGenJets->at(BHadronIndex).Pt(), weight );
        h_GenRecoBJetRapidity->Fill(HypBJet->at(solutionIndex).Rapidity(), allGenJets->at(BHadronIndex).Rapidity(), weight );
        h_GenRecoBJetEta->Fill(HypBJet->at(solutionIndex).Eta(), allGenJets->at(BHadronIndex).Eta(), weight );
        h_GenRecoAntiLeptonBjetMass->Fill(( HypAntiLepton->at(solutionIndex)+HypBJet->at(solutionIndex) ).M(), ( *GenAntiLepton+allGenJets->at(BHadronIndex) ).M(), weight );
    } else {
        h_GenRecoBJetpT->Fill(HypBJet->at(solutionIndex).Pt(), -1000., weight );
        h_GenRecoBJetRapidity->Fill(HypBJet->at(solutionIndex).Rapidity(), -1000., weight );
        h_GenRecoBJetEta->Fill(HypBJet->at(solutionIndex).Eta(), -1000., weight );
        h_GenRecoAntiLeptonBjetMass->Fill(( HypAntiLepton->at(solutionIndex) + HypBJet->at(solutionIndex) ).M(), -1000., weight );
    }
    if ( AntiBHadronIndex>=0 ) {
        h_GenRecoAntiBJetpT->Fill(HypAntiBJet->at(solutionIndex).Pt(), allGenJets->at(AntiBHadronIndex).Pt(), weight );
        h_GenRecoAntiBJetRapidity->Fill(HypAntiBJet->at(solutionIndex).Rapidity(), allGenJets->at(AntiBHadronIndex).Rapidity(), weight );
        h_GenRecoAntiBJetEta->Fill(HypAntiBJet->at(solutionIndex).Eta(), allGenJets->at(AntiBHadronIndex).Eta(), weight );
        h_GenRecoLeptonantiBjetMass->Fill(( HypLepton->at(solutionIndex) + HypAntiBJet->at(solutionIndex) ).M(), ( *GenLepton+allGenJets->at(AntiBHadronIndex) ).M(), weight );
    } else {
        h_GenRecoAntiBJetpT->Fill(HypAntiBJet->at(solutionIndex).Pt(), -1000., weight );
        h_GenRecoAntiBJetRapidity->Fill(HypAntiBJet->at(solutionIndex).Rapidity(), -1000., weight );
        h_GenRecoAntiBJetEta->Fill(HypAntiBJet->at(solutionIndex).Eta(), -1000., weight );
        h_GenRecoLeptonantiBjetMass->Fill(( HypLepton->at(solutionIndex) + HypAntiBJet->at(solutionIndex) ).M(), -1000., weight );
    }
    
//     int firstGenJet = -1, secondGenJet = -1;
//     int firstHypJet = -1, secondHypJet = -1;
//     for(int k=0; k<allGenJets->size(); k++){
//         if(abs(allGenJets->at(k).Eta())>2.4 || allGenJets->at(k).Pt()< 30.0) {continue;}
//         if(firstGenJet<0) {
//             firstGenJet=k;
//             continue;//jump to next iteration k -> k+1
//         }
//         if(firstGenJet>=0 && secondGenJet<0){
//             secondGenJet=k;
//             break; //exit 'for' loop, found interesting 1st and 2nd jets
//         }
//     }
//     for(int k=0; k<jet.size(); k++){
//         if(abs(jet->at(k).Eta())>2.4 || jet->at(k).Pt()< 30.0) {continue;}
//         if(firstHypJet<0) {
//             firstHypJet=k;
//             continue;//jump to next iteration k -> k+1
//         }
//         if(firstHypJet>=0 && secondHypJet<0){
//             secondHypJet=k;
//             break; //exit 'for' loop, found interesting 1st and 2nd jets
//         }
//     }
//     if(firstGenJet >=0 && secondGenJet>=0 && firstHypJet >=0 && secondHypJet >=0 ){
//         h_GenRecoLeadingJetpT->Fill(jet->at(firstHypJet).Pt(), allGenJets->at(firstGenJet).Pt(),weight);
//         h_GenRecoLeadingJetEta->Fill(jet->at(firstHypJet).Eta(), allGenJets->at(firstGenJet).Eta(),weight);
//         h_GenRecoNLeadingJetpT->Fill(jet->at(secondHypJet).Pt(), allGenJets->at(secondGenJet).Pt(),weight);
//         h_GenRecoNLeadingJetEta->Fill(jet->at(secondHypJet).Eta(), allGenJets->at(secondGenJet).Eta(),weight);
//     }
    
    if ( BHadronIndex>=0 && AntiBHadronIndex>=0 ) {
        h_GenRecoJetMult->Fill(jet->size(), allGenJets->size(), weight );
    } else {
        h_GenRecoJetMult->Fill(jet->size(), -1000., weight );
    }

    LV genllbar(*GenLepton + *GenAntiLepton);
    h_GenRecoLLBarMass->Fill(hypllbar.M(), genllbar.M(), weight );
    h_GenRecoLLBarpT->Fill(hypllbar.Pt(), genllbar.Pt(), weight );

    LV genttbar(*GenTop + *GenAntiTop);
    h_GenRecoTTBarMass->Fill(hypttbar.M(), genttbar.M(), weight );
    h_GenRecoTTBarpT->Fill(hypttbar.Pt(), genttbar.Pt(), weight );
    h_GenRecoTTBarRapidity->Fill(hypttbar.Rapidity(), genttbar.Rapidity(), weight );

    //finally do the btag SF calculation stuff
    for (int i=0; i < jet->size(); ++i) {
        if (jet->at(i).Pt() <= JETPTCUT) break;
        if (TMath::Abs(jet->at(i).Eta())<2.4) {
            int type = (*jetType)[i];
            if(type == 2){//b-quark
                h_bjets->Fill(jet->at(i).Pt(), TMath::Abs(jet->at(i).Eta()));
                if((*jetBTagCSV)[i]>BtagWP){
                    h_btaggedjets->Fill(jet->at(i).Pt(), TMath::Abs(jet->at(i).Eta()));
                }
            }
            else if (type == 1){//c-quark
                h_cjets->Fill(jet->at(i).Pt(), TMath::Abs(jet->at(i).Eta()));
                if((*jetBTagCSV)[i]>BtagWP){
                    h_ctaggedjets->Fill(jet->at(i).Pt(), TMath::Abs(jet->at(i).Eta()));
                }
            }
            else if (type == 0){//l-quark
                h_ljets->Fill(jet->at(i).Pt(), TMath::Abs(jet->at(i).Eta()));
                if((*jetBTagCSV)[i]>BtagWP){
                    h_ltaggedjets->Fill(jet->at(i).Pt(), TMath::Abs(jet->at(i).Eta()));
                }
            }
            else {
                cout<<"I found a jet in event "<<eventNumber<<" which is not b, c nor ligth"<<endl; 
                return kFALSE;
            }
        }
    }

    return kTRUE;
}

void Analysis::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.
    
}

void Analysis::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

    string f_savename = "selectionRoot/";
    gSystem->MakeDirectory( f_savename.c_str() );
    f_savename.append ( systematic );
    gSystem->MakeDirectory( f_savename.c_str() );
    f_savename.append ( "/" );
    f_savename.append ( channel );
    gSystem->MakeDirectory( f_savename.c_str() );
    f_savename.append ( "/" );
    f_savename.append ( outputfilename );
    //f_savename.append ( ".root" );

    std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!Finishing: "<<samplename<<"!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;

    //write stuff into file
    TFile f(f_savename.c_str(), "RECREATE");
    TIterator* it = fOutput->MakeIterator();
    while (TObject* obj = it->Next()) {
        obj->Write();
        //cout << obj->GetName() << "\n";
    }
    weightedEvents->Write();
    TObjString(channel).Write("channelName");
    TObjString(systematic).Write("systematicsName");
    TObjString(samplename).Write("sampleName");
    TObjString(isSignal ? "1" : "0").Write("isSignal");
    TObjString(isMC ? "1" : "0").Write("isMC");
    f.Close();
    
    cout<<"Created: "<<f_savename<<endl;
    
    if (isSignal) {
        cout << "Signal sample, writing out btag efficiencies\n";
        f_savename = "selectionRoot/BTagEff";
        gSystem->MakeDirectory(f_savename.c_str());
        f_savename.append("/");
        f_savename.append(systematic); 
        gSystem->MakeDirectory(f_savename.c_str());
        f_savename.append("/");
        f_savename.append(channel); 
        gSystem->MakeDirectory(f_savename.c_str());
        f_savename.append("/");
        f_savename.append(outputfilename);
        
        h_bjets = dynamic_cast<TH2*>( fOutput->FindObject("bjets2D") );
        h_btaggedjets = dynamic_cast<TH2*>( fOutput->FindObject("bjetsTagged2D") );
        h_cjets = dynamic_cast<TH2*>( fOutput->FindObject("cjets2D") );
        h_ctaggedjets = dynamic_cast<TH2*>( fOutput->FindObject("cjetsTagged2D") );
        h_ljets = dynamic_cast<TH2*>( fOutput->FindObject("ljets2D") );
        h_ltaggedjets = dynamic_cast<TH2*>( fOutput->FindObject("ljetsTagged2D") );
        if (!h_bjets || !h_btaggedjets || !h_cjets || !h_ctaggedjets || !h_ljets || !h_ltaggedjets) {
            cerr << "At least one of the btag histograms is missing\n";
            exit(4);
        }
        TFile fbtag(f_savename.c_str(),"RECREATE");
        h_bjets->Write();
        h_btaggedjets->Write();
        h_cjets->Write();
        h_ctaggedjets->Write();
        h_ljets->Write();
        h_ltaggedjets->Write();
        
        TH1 *btaggedPt = h_btaggedjets->ProjectionX(); TH1 *btaggedEta = h_btaggedjets->ProjectionY();
        TH1 *ctaggedPt = h_ctaggedjets->ProjectionX(); TH1 *ctaggedEta = h_ctaggedjets->ProjectionY();
        TH1 *ltaggedPt = h_ltaggedjets->ProjectionX(); TH1 *ltaggedEta = h_ltaggedjets->ProjectionY();
        
        TH1 *bUntaggedPt = h_bjets->ProjectionX(); TH1 *bEta = h_bjets->ProjectionY();
        TH1 *cUntaggedPt = h_cjets->ProjectionX(); TH1 *cEta = h_cjets->ProjectionY();
        TH1 *lUntaggedPt = h_ljets->ProjectionX(); TH1 *lEta = h_ljets->ProjectionY();
        
        //Calculate the medians and save them in a txt file
        double PtMedian = Median(btaggedPt);
        double EtaMedian = Median(btaggedEta);
        printf("Median: pT = %.0f, eta = %.2f\n", PtMedian, EtaMedian);
        TH1* medianHist = new TH1D("Medians", "medians", 2, -0.5, 1.5);
        medianHist->GetXaxis()->SetBinLabel(1, "pT");
        medianHist->GetXaxis()->SetBinLabel(2, "eta");
        medianHist->SetBinContent(1, PtMedian);
        medianHist->SetBinContent(2, EtaMedian);
        medianHist->Write();
        
        TH1 *beffPt =(TH1*) btaggedPt->Clone("beffPt");
        TH1 *ceffPt =(TH1*) ctaggedPt->Clone("ceffPt");
        TH1 *leffPt =(TH1*) ltaggedPt->Clone("leffPt");
        
        TH1 *beffEta =(TH1*) btaggedEta->Clone("beffEta");  
        TH1 *ceffEta =(TH1*) ctaggedEta->Clone("ceffEta");  
        TH1 *leffEta =(TH1*) ltaggedEta->Clone("leffEta");  
        
        //Calculate Efficiency: N_tageed/N_all
        //Calculate also the binomial error (option "B" does it)!!
        beffPt->Divide(btaggedPt, bUntaggedPt, 1, 1, "B"); 
        ceffPt->Divide(ctaggedPt, cUntaggedPt, 1, 1, "B"); 
        leffPt->Divide(ltaggedPt, lUntaggedPt, 1, 1, "B");
        beffEta->Divide(btaggedEta, bEta, 1, 1, "B"); 
        ceffEta->Divide(ctaggedEta, cEta, 1, 1, "B"); 
        leffEta->Divide(ltaggedEta, lEta, 1, 1, "B"); 
        h_btaggedjets->Divide(h_btaggedjets, h_bjets, 1, 1, "B"); 
        h_ctaggedjets->Divide(h_ctaggedjets, h_cjets, 1, 1, "B"); 
        h_ltaggedjets->Divide(h_ltaggedjets, h_ljets, 1, 1, "B"); 

        //Save histograms in ROOT file
        beffPt->Write("BEffPt"); 
        ceffPt->Write("CEffPt"); 
        leffPt->Write("LEffPt"); 
        beffEta->Write("BEffEta"); 
        ceffEta->Write("CEffEta"); 
        leffEta->Write("LEffEta"); 
        h_btaggedjets->Write("BEffPerJet");
        h_ctaggedjets->Write("CEffPerJet");
        h_ltaggedjets->Write("LEffPerJet");
        
        fbtag.Close();
    }
    fOutput->SetOwner();
    fOutput->Clear();
}

double Analysis::BJetSF( double pt, double eta )
{
    //CSVL b-jet SF
    //From BTV-11-004 and https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-mujet_payload.txt

    if ( abs(eta) > 2.4 ) {
        cout<<"Jet Eta out of the selected range. Check it"<<endl;
        exit(9);
    }
    
    if ( pt < 30 ) pt = 30;
    if ( pt > 670 ) pt = 670;

    return 1.02658*((1.+(0.0195388*pt))/(1.+(0.0209145*pt)));
}

double Analysis::CJetSF ( double pt, double eta )
{
    //CSVL c-jet SF
    //From BTV-11-004 and https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC
    return BJetSF( pt, eta );
}

double Analysis::LJetSF ( double pt, double eta )
{
    //CSVL ligth jet mistag SF.
    //From BTV-11-004 and https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlightFuncs.C

    double eta_abs = abs(eta);
    if (eta_abs > 2.4) {
        cout<<"There is a jet out of the selected ETA region. Check that!!!!!\n";
        return 1;
    }
    if ( pt > 670 ) {
        return (((0.956023+(0.000825106*pt))+(-3.18828e-06*(pt*pt)))+(2.81787e-09*(pt*(pt*pt)))) * (0.979396 + 0.000205898*pt + 2.49868e-07*pt*pt);
    } else {
        if ( eta_abs <= 0.5 ) {
            return (((0.994425+(-8.66392e-05*pt))+(-3.03813e-08*(pt*pt)))+(-3.52151e-10*(pt*(pt*pt)))) * (0.979396 + 0.000205898*pt + 2.49868e-07*pt*pt);
        } else if ( eta_abs <= 1.0 ) {
            return (((0.998088+(6.94916e-05*pt))+(-4.82731e-07*(pt*pt)))+(1.63506e-10*(pt*(pt*pt)))) * (0.979396 + 0.000205898*pt + 2.49868e-07*pt*pt);
        } else if ( eta_abs <= 1.5 ) {
            return (((1.00294+(0.000289844*pt))+(-7.9845e-07*(pt*pt)))+(5.38525e-10*(pt*(pt*pt)))) * (0.979396 + 0.000205898*pt + 2.49868e-07*pt*pt);
        } else {
            return (((0.979816+(0.000138797*pt))+(-3.14503e-07*(pt*pt)))+(2.38124e-10*(pt*(pt*pt)))) * (0.979396 + 0.000205898*pt + 2.49868e-07*pt*pt);
        }
    }
}

double Analysis::BJetSFAbsErr ( int ptbin )
{
    //c- and l-jets errors are not necessary for the calculation and are not implemented. If needed go to https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC

    //b-jet pt ranges {0, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 670};
    //this pt range matches the pTEff histogram!!

    double SFb_error[] = {0.1388743, 0.0188743, 0.0161816, 0.0139824, 0.0152644, 0.0161226, 0.0157396, 0.0161619, 0.0168747, 0.0257175, 0.026424, 0.0264928, 0.0315127, 0.030734, 0.0438259 };

    if ( ptbin > 14 ) {
        return 1.5 * 2 * SFb_error[14];
    } else {
        return 1.5 * SFb_error[ptbin];
    };
}

void Analysis::SetBTagFile(TString btagFile)
{
    this->btagFile = btagFile;
}

void Analysis::SetChannel(TString channel)
{
    this->channel = channel;
}

void Analysis::SetSignal(bool isSignal)
{
    this->isSignal = isSignal;
}

void Analysis::SetSystematic(TString systematic)
{
    this->systematic = systematic;
}

void Analysis::SetSamplename(TString samplename)
{
    this->samplename = samplename;
    isTtbarPlusTauSample = samplename.BeginsWith("ttbar") && !samplename.Contains("bg");
    correctMadgraphBR = samplename.BeginsWith("ttbar");
    //lumiWeight = 5100*SampleXSection(samplename)/weightedEvents->GetBinContent(1);
    lumiWeight = 12100*SampleXSection(samplename)/weightedEvents->GetBinContent(1);
}

void Analysis::SetMC(bool isMC)
{
    this->isMC = isMC;
}

void Analysis::SetOutputfilename(TString outputfilename)
{
    if (outputfilename.Contains('/')) {
        Ssiz_t last = outputfilename.Last('/');
        this->outputfilename = outputfilename.Data() + last + 1;
    } else {
        this->outputfilename = outputfilename;
    }
}

void Analysis::SetWeightedEvents(TH1* weightedEvents)
{
    this->weightedEvents = weightedEvents;
}

void Analysis::SetRunViaTau(bool runViaTau)
{
    this->runViaTau = runViaTau;
    if (runViaTau) isSignal = 0;
}

void Analysis::SetPUReweighter(PUReweighter* pu)
{
    pureweighter = pu;
}

void Analysis::Init ( TTree *tree )
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set object pointer
    lepton = 0;
    lepQ = 0;
    lepType = 0;
    lepPfIso = 0;
    lepCombIso = 0;
    jet = 0;
    jetBTagTCHE = 0;
    jetBTagCSV = 0;
    jetBTagSSVHE = 0;
    jetType = 0;
    met = 0;
    HypJet0index = 0;
    HypJet1index = 0;
    HypTop = 0;
    HypAntiTop = 0;
    HypLepton = 0;
    HypAntiLepton = 0;
    HypNeutrino = 0;
    HypAntiNeutrino = 0;
    HypBJet = 0;
    HypAntiBJet = 0;

    //for the signal
    genJet = 0;
    allGenJets = 0;
    BHadrons = 0;
    GenWPlus = 0;
    GenWMinus = 0;
    AntiBHadrons = 0;
    BHadJetIndex = 0;
    AntiBHadJetIndex = 0;
    BHadronFromTopB = 0;
    AntiBHadronFromTopB = 0;
    BHadronVsJet = 0;
    AntiBHadronVsJet = 0;
    GenNeutrino = 0;
    GenAntiNeutrino = 0;
    GenB = 0;
    GenAntiB = 0;
    GenLepton = 0;
    GenAntiLepton = 0;
    GenTop = 0;
    GenAntiTop = 0;

    // Set branch addresses and branch pointers
    if ( !tree ) return;
    fChain = tree;
    fChain->SetMakeClass ( 0 );
    fChain->SetBranchAddress("lepton", &lepton, &b_lepton );
    fChain->SetBranchAddress("lepQ", &lepQ, &b_lepQ );
    fChain->SetBranchAddress("lepType", &lepType, &b_lepType );
    fChain->SetBranchAddress("lepPfIso", &lepPfIso, &b_lepPfIso );
    fChain->SetBranchAddress("lepCombIso", &lepCombIso, &b_lepCombIso );
    fChain->SetBranchAddress("jet", &jet, &b_jet );
    fChain->SetBranchAddress("jetBTagTCHE", &jetBTagTCHE, &b_jetBTagTCHE );
    fChain->SetBranchAddress("jetBTagCSV", &jetBTagCSV, &b_jetBTagCSV );
    fChain->SetBranchAddress("jetBTagSSVHE", &jetBTagSSVHE, &b_jetBTagSSVHE );
    fChain->SetBranchAddress("jetType", &jetType, &b_jetType );
    fChain->SetBranchAddress("genJet", &genJet, &b_genJet );
    fChain->SetBranchAddress("met", &met, &b_met );
    fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber );
    fChain->SetBranchAddress("triggerBits", &triggerBits, &b_triggerBits );
    fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock );
    fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber );
    fChain->SetBranchAddress("weightGenerator", &weightGenerator, &b_weightGenerator );
    fChain->SetBranchAddress("vertMulti", &vertMulti, &b_vertMulti );
    fChain->SetBranchAddress("vertMultiTrue", &vertMultiTrue, &b_vertMultiTrue );


    fChain->SetBranchAddress("allGenJets", &allGenJets, &b_allGenJets );
    fChain->SetBranchAddress("HypTop", &HypTop, &b_HypTop );
    fChain->SetBranchAddress("HypAntiTop", &HypAntiTop, &b_HypAntiTop );
    fChain->SetBranchAddress("HypLepton", &HypLepton, &b_HypLepton );
    fChain->SetBranchAddress("HypAntiLepton", &HypAntiLepton, &b_HypAntiLepton );
    fChain->SetBranchAddress("HypNeutrino", &HypNeutrino, &b_HypNeutrino);
    fChain->SetBranchAddress("HypAntiNeutrino", &HypAntiNeutrino, &b_HypAntiNeutrino);
    fChain->SetBranchAddress("HypB", &HypBJet, &b_HypB );
    fChain->SetBranchAddress("HypAntiB", &HypAntiBJet, &b_HypAntiB );
    /*   fChain->SetBranchAddress("HypWPlus", &HypWPlus_, &b_HypWPlus_);
    fChain->SetBranchAddress("HypWMinus", &HypWMinus_, &b_HypWMinus_);
    */
    fChain->SetBranchAddress("HypJet0index", &HypJet0index, &b_HypJet0index );
    fChain->SetBranchAddress("HypJet1index", &HypJet1index, &b_HypJet1index );
    fChain->SetBranchAddress("decayMode", &decayMode, &b_decayMode );
    
    if (isSignal) {
        fChain->SetBranchAddress("GenTop", &GenTop, &b_GenTop );
        fChain->SetBranchAddress("GenAntiTop", &GenAntiTop, &b_GenAntiTop );
        fChain->SetBranchAddress("GenLepton", &GenLepton, &b_GenLepton );
        fChain->SetBranchAddress("GenAntiLepton", &GenAntiLepton, &b_GenAntiLepton );
        fChain->SetBranchAddress("GenNeutrino", &GenNeutrino, &b_GenNeutrino);
        fChain->SetBranchAddress("GenAntiNeutrino", &GenAntiNeutrino, &b_GenAntiNeutrino);
        fChain->SetBranchAddress("GenB", &GenB, &b_GenB );
        fChain->SetBranchAddress("GenAntiB", &GenAntiB, &b_GenAntiB );
        /*  
        fChain->SetBranchAddress("GenWPlus.fCoordinates.fX", &GenWPluspX, &b_GenWPluspX);
        fChain->SetBranchAddress("GenWMinus.fCoordinates.fX", &GenWMinuspX, &b_GenWMinuspX);
        */
        fChain->SetBranchAddress ( "BHadJetIndex", &BHadJetIndex, &b_BHadJetIndex );
        fChain->SetBranchAddress ( "AntiBHadJetIndex", &AntiBHadJetIndex, &b_AntiBHadJetIndex );
        fChain->SetBranchAddress ( "BHadrons", &BHadrons, &b_BHadrons );
        fChain->SetBranchAddress ( "AntiBHadrons", &AntiBHadrons, &b_AntiBHadrons);
        fChain->SetBranchAddress ( "BHadronFromTop", &BHadronFromTopB, &b_BHadronFromTopB );
        fChain->SetBranchAddress ( "AntiBHadronFromTopB", &AntiBHadronFromTopB, &b_AntiBHadronFromTopB );
        fChain->SetBranchAddress ( "BHadronVsJet", &BHadronVsJet, &b_BHadronVsJet );
        fChain->SetBranchAddress ( "AntiBHadronVsJet", &AntiBHadronVsJet, &b_AntiBHadronVsJet );

//         fChain->SetBranchAddress("GenJetHadronB.", &BHadronJet_, &b_BHadronJet_);
//         fChain->SetBranchAddress("GenJetHadronAntiB", &AntiBHadronJet_, &b_AntiBHadronJet_);
    }    
}

Bool_t Analysis::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}


void Analysis::GetRecoBranches ( Long64_t & entry )
{

    b_met->GetEntry(entry); //!
    b_eventNumber->GetEntry(entry); //!
    b_lepton->GetEntry(entry); //!
    b_jet->GetEntry(entry); //!
    b_lepQ->GetEntry(entry); //!
    b_lepType->GetEntry(entry); //!
    b_lepPfIso->GetEntry(entry); //!
    b_lepCombIso->GetEntry(entry); //!
    b_jetBTagTCHE->GetEntry(entry); //!
    b_jetBTagCSV->GetEntry(entry); //!
    b_jetBTagSSVHE->GetEntry(entry); //!
    b_jetType->GetEntry(entry); //!
    b_runNumber->GetEntry(entry); //!
    b_triggerBits->GetEntry(entry); //!
    b_lumiBlock->GetEntry(entry); //!
    b_weightGenerator->GetEntry(entry); //!
    b_vertMulti->GetEntry(entry); //!
    b_vertMultiTrue->GetEntry(entry); //!

    b_genJet->GetEntry(entry); //!
    b_allGenJets->GetEntry(entry); //!

    b_HypTop->GetEntry(entry); //!
    b_HypAntiTop->GetEntry(entry); //!
    b_HypLepton->GetEntry(entry); //!
    b_HypAntiLepton->GetEntry(entry); //!
    b_HypB->GetEntry(entry); //!
    b_HypAntiB->GetEntry(entry); //!
    b_HypNeutrino->GetEntry(entry);   //!
    b_HypAntiNeutrino->GetEntry(entry);   //!

    /* b_HypWPlus_->GetEntry(entry);   //!
    b_HypWPluspX->GetEntry(entry);   //!
    b_HypWPluspY->GetEntry(entry);   //!
    b_HypWPluspZ->GetEntry(entry);   //!
    b_HypWPlusE->GetEntry(entry);   //!

    b_HypWMinus_->GetEntry(entry);   //!
    b_HypWMinuspX->GetEntry(entry);   //!
    b_HypWMinuspY->GetEntry(entry);   //!
    b_HypWMinuspZ->GetEntry(entry);   //!
    b_HypWMinusE->GetEntry(entry);   //!
    */
    b_HypJet0index->GetEntry(entry);
    b_HypJet1index->GetEntry(entry);
    b_decayMode->GetEntry(entry);

}

void Analysis::GetSignalBranches ( Long64_t & entry )
{
    b_GenTop->GetEntry(entry, 1); //!
    b_GenAntiTop->GetEntry(entry); //!
    b_GenLepton->GetEntry(entry); //!
    b_GenAntiLepton->GetEntry(entry); //!
    b_GenB->GetEntry(entry); //!
    b_GenAntiB->GetEntry(entry); //!
    b_GenNeutrino->GetEntry(entry);   //!
    b_GenAntiNeutrino->GetEntry(entry);   //!

    /*
    b_GenWPluspX->GetEntry(entry);   //!
    b_GenWMinuspX->GetEntry(entry);   //!
    */
    b_BHadJetIndex->GetEntry(entry); //!
    b_AntiBHadJetIndex->GetEntry(entry); //!

    b_BHadrons->GetEntry(entry); //!
    b_AntiBHadrons->GetEntry(entry); //!

    b_BHadronFromTopB->GetEntry(entry); //!
    b_AntiBHadronFromTopB->GetEntry(entry); //!
    b_BHadronVsJet->GetEntry(entry); //!
    b_AntiBHadronVsJet->GetEntry(entry); //!

    /*  
    b_BHadronJet->GetEntry(entry);   //!
    b_AntiBHadronJet->GetEntry(entry);   //!
    */
}

bool Analysis::getLeptonPair(size_t &LeadLeptonNumber, size_t &NLeadLeptonNumber)
{
    if ( lepton->size() > 1 ) {
        if ( channel == "emu" ) { //quick and DIRTY!
            for ( size_t i = 1; i < lepton->size(); i++ ) {
                if ( ( ( *lepQ )[0] != ( *lepQ )[i] ) && ( (*lepType)[0] != (*lepType)[i] ) ) {
                    LeadLeptonNumber = 0;
                    NLeadLeptonNumber = i;
                    return true;
                }
            }
        }
        if ( channel == "ee" ) { //quick and DIRTY!
            for ( size_t i = 0; i < lepton->size(); i++ ) {
                if ( ( *lepType ) [i] == LEP_TYPE_ELECTRON ) {
                    LeadLeptonNumber=i;
                    break;
                }
            }
            for ( size_t i = LeadLeptonNumber+1; i < lepton->size(); i++ ) {
                if ( ( ( *lepQ ) [LeadLeptonNumber]!= ( *lepQ ) [i] ) && ( *lepType ) [i] == LEP_TYPE_ELECTRON ) {
                    NLeadLeptonNumber = i;
                    return true;
                }
            }
        }
        if ( channel == "mumu" ) { //quick and DIRTY!
            for ( size_t i = 0; i < lepton->size(); i++ ) {
                if ( ( *lepType )[i] == LEP_TYPE_MUON ) {
                    LeadLeptonNumber=i;
                    break;
                }
            }
            for ( size_t i = LeadLeptonNumber+1; i<lepton->size(); i++ ) {
                if ( ( ( *lepQ )[LeadLeptonNumber]!= ( *lepQ )[i] ) && ( *lepType )[i] == LEP_TYPE_MUON ) {
                    NLeadLeptonNumber = i;
                    return true;
                }
            }
        }
    }
    return false;
}

double Analysis::get2DSF(TH2* histo, const double x, const double y)
{
    int xbin, ybin, dummy;
    histo->GetBinXYZ(histo->FindBin(x, y), xbin, ybin, dummy);
    //overflow to last bin
    xbin = std::min(xbin, histo->GetNbinsX());
    ybin = std::min(ybin, histo->GetNbinsY());
    return histo->GetBinContent(xbin, ybin);
}


double Analysis::calculateBtagSF()
{
    if (!bEff) return 1; //no btag file given, so return 1
    
    double OneMinusEff=1;
    double OneMinusSEff=1;
    double SFPerJet=1, eff=1;
    for ( size_t i = 0; i < jet->size(); ++i ) {
        double pt = jet->at(i).Pt();
        double eta = abs(jet->at(i).Eta());
        if ( pt > 30 && eta < 2.4 ) {
            int ptbin, etabin, dummy;
            bEff->GetBinXYZ(bEff->FindBin(pt, eta), ptbin, etabin, dummy);
            //overflow to last bin
            ptbin = std::min(ptbin, bEff->GetNbinsX());
            etabin = std::min(etabin, bEff->GetNbinsY());
            //do the type-jet selection & Eff and SF obtention
            double SF_Error=0;
            if ( ( *jetType )[i] == 2 ) { //b-quark
                eff=bEff->GetBinContent ( ptbin, etabin );
                SFPerJet=BJetSF( pt, eta );
                SF_Error = BJetSFAbsErr ( ptbin );
            } else if ( ( *jetType )[i] == 1 ) { //c-quark
                SFPerJet=CJetSF( pt, eta );
                eff=cEff->GetBinContent ( ptbin, etabin );
            } else if ( ( *jetType )[i] == 0 ) { //l-quark
                SFPerJet=LJetSF( pt, eta );
                eff=lEff->GetBinContent ( ptbin, etabin );
            } else {
                cout<<"I found a jet in event "<<eventNumber<<" which is not b, c nor light"<<endl;
                return kFALSE;
            }
            if ( eff <= 0 ) eff = 1;
            //calculate both numerator and denominator for per-event SF calculation
            //consider also the UP and DOWN variation for systematics calculation. Same procedure as PU
            OneMinusEff = OneMinusEff* ( 1-eff );
            OneMinusSEff= OneMinusSEff* ( 1-SFPerJet*eff );
            double sf = SFPerJet;
            if ( systematic == "BTAG_UP" ) {
                sf = SFPerJet + SF_Error;
            }
            else if ( systematic == "BTAG_DOWN" ) {
                sf = SFPerJet - SF_Error;
            }
            else if ( systematic == "BTAG_PT_UP" ) {
                if ( pt>btag_ptmedian )  {
                    sf = SFPerJet - 0.5 * SF_Error;
                } else {
                    sf = SFPerJet + 0.5 * SF_Error;
                }
            }
            else if ( systematic == "BTAG_PT_DOWN" ) {
                if ( pt>btag_ptmedian )  {
                    sf = SFPerJet + 0.5 * SF_Error;
                } else {
                    sf = SFPerJet - 0.5 * SF_Error;
                }
            }
            else if ( systematic == "BTAG_ETA_UP" ) {
                if ( eta>btag_etamedian )  {
                    sf = SFPerJet - 0.5 * SF_Error;
                } else {
                    sf = SFPerJet + 0.5 * SF_Error;
                }
            }
            else if ( systematic == "BTAG_ETA_DOWN" ) {
                if ( eta>btag_etamedian )  {
                    sf = SFPerJet + 0.5 * SF_Error;
                } else {
                    sf = SFPerJet - 0.5 * SF_Error;
                }
            }
            
            OneMinusSEff *= 1 - eff * sf;
        }
    }
    //per-event SF calculation (also the UP and DOWN variations)
    return ( 1.-OneMinusSEff ) / ( 1.-OneMinusEff );
}

double Analysis::getJetHT(const VLV& jet, int pt_cut)
{
    double result = 0;
    for ( size_t i = 0; i < jet.size(); ++i ) {
        double pt = jet.at(i).Pt();
        if (pt < pt_cut) break;
        result += pt;
    }
    return result;
}

void Analysis::prepareTriggerSF()
{
    h_TrigSFeta = nullptr;
    
    TFile trigEfficiencies(TString("triggerSummary_").Append(channel).Append(".root"));
    if (trigEfficiencies.IsZombie()) {
        cout << "Trigger efficiencies not found. Assuming ScaleFactor = 1.\n";
        cout << "Currently triggerEfficieny files can be found in Jan's NAF public afs\n\n";
        return;
    }
    
    //Right now pT efficiency flat ==> Not used
    h_TrigSFeta = dynamic_cast<TH1*>(trigEfficiencies.Get("TH scalefactor eta incl corrErr"));
    if ( !h_TrigSFeta ) {
        cout<<"TH1 >>TH scalefactor eta<< is not in the file "<<trigEfficiencies.GetName()<<"\n";
        return;
    }
    
    if (systematic.BeginsWith("TRIG_")) {
        double factor = systematic.EndsWith("_UP") ? 1 : -1;
        for (int i = 1; i <= h_TrigSFeta->GetNbinsX(); ++i) {
            h_TrigSFeta->SetBinContent(i, 
                h_TrigSFeta->GetBinContent(i) + factor*h_TrigSFeta->GetBinError(i));
        }
    }
    
    h_TrigSFeta->SetDirectory(0);
    trigEfficiencies.Close();
}

double Analysis::getTriggerSF(const LV& lep1, const LV& lep2) {
    if (!h_TrigSFeta) return 1;
    return TMath::Sqrt(        
        h_TrigSFeta->GetBinContent(h_TrigSFeta->FindBin(lep1.eta()))
        * h_TrigSFeta->GetBinContent(h_TrigSFeta->FindBin(lep2.eta())));
}

double Analysis::getLeptonIDSF(const LV& lep1, const LV& lep2) {
    if (!h_LepIDSFpteta) return 1;
    return TMath::Sqrt(
        get2DSF(h_LepIDSFpteta, lep1.pt(), lep1.Eta())
        * get2DSF(h_LepIDSFpteta, lep2.pt(), lep2.Eta()));
}

void Analysis::prepareLeptonIDSF()
{
    h_LepIDSFpteta = nullptr;
    std::cout << "Please implement reading the TH2* for the lepton SF\n\n";
    //...
}


void Analysis::prepareBtagSF()
{
    //some defaults for the median, overwritten if btag files exist
    btag_ptmedian = 75;
    btag_etamedian = 0.75;

    //By now defined the per-jet SFs vary according to:
    //   BTag_Up   ==> pt>ptmedian vary DOWN, pt<ptmedian vary UP
    //   BTag_Down ==> pt>ptmedian vary UP, pt<ptmedian vary DOWN

    //load per-jet efficienciies file and Histograms
    TFile *bEfficiencies;
    if (btagFile!="") {
        bEfficiencies = TFile::Open(btagFile);
    } else {
        cout<<"WARNING!!! Provide b tag efficiencies before running"<<endl;
        return;
    }

    if (!bEfficiencies) {
        cout << "\n******************************************************\n"
             << "File " << btagFile << " does not exist. Running without btagsf!!!\n"
             << "To create the file, run:\n" 
             << "   ./load_Analysis -f ttbarsignal\n"
             << "and copy the selectionRoot/BTagEff directory to the cwd:\n"
             << "   cp -r selectionRoot/BTagEff .\n"
             << "This error is NOT fatal, using a btag SF = 1 everywhere\n"
             << "*******************************************************\n\n";
        return;
    }
    bEff = dynamic_cast<TH2*>(bEfficiencies->Get("BEffPerJet"));
    if (!bEff) {
        cout<<"Histogram bEff is not in the file "<<bEfficiencies->GetName();
        return;
    }
    cEff = dynamic_cast<TH2*>(bEfficiencies->Get("CEffPerJet"));
    if (!cEff) {
        cout<<"Histogram cEff is not in the file "<<bEfficiencies->GetName();
        return;
    }
    lEff = dynamic_cast<TH2*>(bEfficiencies->Get("LEffPerJet"));
    if (!lEff) {
        cout<<"Histogram lEff is not in the file "<<bEfficiencies->GetName();
        return;
    }
    
    TH1* medians = dynamic_cast<TH1*>(bEfficiencies->Get("Medians"));
    btag_ptmedian = medians->GetBinContent(1);
    btag_etamedian = medians->GetBinContent(2);
    printf("BTagSF: Using medians: pT = %.0f, eta = %.2f\n", btag_ptmedian, btag_etamedian);

    //load the histograms in memory, to avoid memory leaks
    bEff->SetDirectory(0);
    cEff->SetDirectory(0);
    lEff->SetDirectory(0);
    bEfficiencies->Close();
    bEfficiencies->Delete();
    // END: BTag SF calculation neccessary stuff

}
