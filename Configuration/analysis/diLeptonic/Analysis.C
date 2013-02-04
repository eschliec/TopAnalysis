#define Analysis_cxx

#include "Analysis.h"
#include "HistoListReader.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <TMath.h>
#include <TSystem.h>
#include <Math/VectorUtil.h>
#include <TLorentzVector.h>
#include <set>
#include <cmath>
#include <TString.h>
#include <limits>
#include <iomanip>
#include "utils.h"
#include "KinReco.h"

using namespace std;
using ROOT::Math::VectorUtil::DeltaPhi;
using ROOT::Math::VectorUtil::DeltaR;

double SampleXSection(TString sample){
    
    //MC cross sections taken from:
    //  https://twiki.cern.ch/twiki/bin/view/CMS/StandardModelCrossSectionsat8TeV
    //  AN-12/194    AN-12/228
    
    if(sample.Contains("data"))        {return 1;}
    if(sample.Contains("ttbar"))       {return 225.197;}
    if(sample.Contains("single"))      {return 11.1;}
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

    //lumiWeight = 5100*SampleXSection(samplename)/weightedEvents->GetBinContent(1);
    lumiWeight = 12100*SampleXSection(samplename)/weightedEvents->GetBinContent(1);
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

    binnedControlPlots = new std::map<std::string, std::pair<TH1*, std::vector<std::map<std::string, TH1*> > > >;
    
    h_step4 = store(new TH1D ( "step4", "event count at after 2lepton", 10, 0, 10 ));
    h_step5 = store(new TH1D ( "step5", "event count at after Zcut", 10, 0, 10 ));
    h_step6 = store(new TH1D ( "step6", "event count at after 2jets", 10, 0, 10 ));
    h_step7 = store(new TH1D ( "step7", "event count at after MET", 10, 0, 10 ));
    h_step8 = store(new TH1D ( "step8", "event count at after 1btag", 10, 0, 10 ));
    h_step9 = store(new TH1D ( "step9", "event count at step after KinReco", 10, 0, 10 ));

    //h_jetMultiAll = store(new TH1D ( "HypjetMultiAll", "Jet Multiplicity (AllJets)", 10, -0.5, 9.5 ));
    h_jetMultiXSec = store(new TH1D ( "HypjetMultiXSec", "Jet Multiplicity (for cross-section)", 10, -0.5, 9.5 ));
    h_jetMulti = store(new TH1D ( "HypjetMulti", "Jet Multiplicity", 10, -0.5, 9.5 ));
    h_jetMulti_diLep = store(new TH1D ( "HypjetMulti_diLep", "Jet Multiplicity (after dilepton)", 10, -0.5, 9.5 ));
    h_jetMulti_noBTag = store(new TH1D ( "HypjetMulti_noBTag", "Jet Multiplicity (after dilepton before btag)", 10, -0.5, 9.5 ));
    h_jetMultiNoPU = store(new TH1D ( "HypjetMultiNoPU", "Jet Multiplicity (No Pileup or lumi weight)", 10, -0.5, 9.5 ));
//     h_jetMultiVisTop = store(new TH1D ( "HypjetMultiVisTop", "Jet Multiplicity for Visible Top (No Pileup or lumi Weight)", 10, -0.5, 9.5 ));
    h_BjetMulti = store(new TH1D ( "HypBjetMulti", "B-Jet Multiplicity", 10, -0.5, 9.5 ));
    h_BjetMulti_noBTag = store(new TH1D ( "HypBjetMulti_noBTag", "B-Jet Multiplicity before B-Tag requirement", 10, -0.5, 9.5 ));

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
    Allh1_postKinReco = store(new TH1D ( "Allh1_postKinReco", "DiLepton Mass", 40, 0, 400 ));
    h_diLepMassFull = store(new TH1D ( "DIMFull", "DiLepton Mass (Full Range)", 100, 0, 300 ));
    h_diLepMassFull_fullSel = store(new TH1D ( "DIMFull_fullSel", "DiLepton Mass (Full Range)", 100, 0, 300 ));
    Looseh1 = store(new TH1D ( "Looseh1", "DiLepton Mass", 40, 0, 400 ));
    Zh1_postKinReco = store(new TH1D ( "Zh1_postKinReco", "DiLepton Mass in Z Window", 40, 0, 400 ));
    TTh1_postKinReco = store(new TH1D ( "TTh1_postKinReco", "DiLepton Mass out of Z Window", 40, 0, 400 ));

    Allh1_postZcut = store(new TH1D ( "Allh1_postZcut", "DiLepton Mass", 40, 0, 400 ));
    Zh1_postZcut = store(new TH1D ( "Zh1_postZcut", "DiLepton Mass in Z Window", 40, 0, 400 ));
    TTh1_postZcut = store(new TH1D ( "TTh1_postZcut", "DiLepton Mass out of Z Window", 40, 0, 400 ));

    Allh1_post2jets = store(new TH1D ( "Allh1_post2jets", "DiLepton Mass", 40, 0, 400 ));
    Zh1_post2jets = store(new TH1D ( "Zh1_post2jets", "DiLepton Mass in Z Window", 40, 0, 400 ));
    TTh1_post2jets = store(new TH1D ( "TTh1_post2jets", "DiLepton Mass out of Z Window", 40, 0, 400 ));

    Allh1_postMET = store(new TH1D ( "Allh1_postMET", "DiLepton Mass", 40, 0, 400 ));
    Zh1_postMET = store(new TH1D ( "Zh1_postMET", "DiLepton Mass in Z Window", 40, 0, 400 ));
    TTh1_postMET = store(new TH1D ( "TTh1_postMET", "DiLepton Mass out of Z Window", 40, 0, 400 ));

    Allh1_post1btag = store(new TH1D ( "Allh1_post1btag", "DiLepton Mass", 40, 0, 400 ));
    Zh1_post1btag = store(new TH1D ( "Zh1_post1btag", "DiLepton Mass in Z Window", 40, 0, 400 ));
    TTh1_post1btag = store(new TH1D ( "TTh1_post1btag", "DiLepton Mass out of Z Window", 40, 0, 400 ));

    
    h_vertMulti = store(new TH1D ( "vertMulti", "Primary Vertex Multiplicity", 30, 0, 30 ));
    h_vertMulti_noPU = store(new TH1D ( "vertMulti_noPU", "Primary Vertex Multiplicity (no Pileup)", 30, 0, 30 ));
    h_MET = store(new TH1D ( "MET", "Missing Transverse Energy", 80, 0, 400 ));
    h_jetpT = store(new TH1D ( "jetpT", "jet pT", 80, 0, 400 ));
    h_jetHT = store(new TH1D ( "jetHT", "jet HT", 80, 0, 1000 ));

    h_MuonpT = store(new TH1D ( "MuonpT", "Muon pT", 80, 0, 400 ));
    h_MuonEta = store(new TH1D ( "MuonEta", "Muon Eta", 100, -5, 5 ));
    h_ElectronpT = store(new TH1D ( "ElectronpT", "Electron pT", 80, 0, 400 ));
    h_ElectronEta = store(new TH1D ( "ElectronEta", "Electron Eta", 100, -5, 5 ));

    h_leptonPtBeforeKinReco = store(new TH1D ( "LeptonpTbkr", "Lepton pT (before kin reco)", 80, 0, 400 ));
    h_leptonPtAfterKinReco = store(new TH1D ( "LeptonpTakr", "Lepton pT (after kin reco)", 80, 0, 400 ));
    h_leptonEtaBeforeKinReco = store(new TH1D ( "LeptonEtabkr", "Lepton #eta (before kin reco)", 80, -2.5, 2.5 ));
    h_leptonEtaAfterKinReco = store(new TH1D ( "LeptonEtaakr", "Lepton #eta (after kin reco)", 80, -2.5, 2.5 ));
    
    h_LeptonpT_postMETcut = store(new TH1D ( "LeptonpT_postMETcut", "Lepton pT (post MET cut)", 80, 0, 400 ));
    h_LeptonEta_postMETcut = store(new TH1D ( "LeptonEta_postMETcut", "Lepton Eta (post MET cut)", 100, -5, 5 ));
    h_AntiLeptonpT_postMETcut = store(new TH1D ( "AntiLeptonpT_postMETcut", "AntiLepton pT (post MET cut)", 80, 0, 400 ));
    h_AntiLeptonEta_postMETcut = store(new TH1D ( "AntiLeptonEta_postMETcut", "AntiLepton Eta (post MET cut)", 100, -5, 5 ));
    h_MuonpT_postMETcut = store(new TH1D ( "MuonpT_postMETcut", "Muon pT (post MET cut)", 80, 0, 400 ));
    h_MuonEta_postMETcut = store(new TH1D ( "MuonEta_postMETcut", "Muon Eta (post MET cut)", 100, -5, 5 ));
    h_ElectronpT_postMETcut = store(new TH1D ( "ElectronpT_postMETcut", "Electron pT (post MET cut)", 80, 0, 400 ));
    h_ElectronEta_postMETcut = store(new TH1D ( "ElectronEta_postMETcut", "Electron Eta (post MET cut)", 100, -5, 5 ));
    
    h_LeptonpT = store(new TH1D ( "LeptonpT", "Lepton pT", 80, 0, 400 ));
    h_LeptonEta = store(new TH1D ( "LeptonEta", "Lepton Eta", 100, -5, 5 ));
    h_LeptonpT_diLep = store(new TH1D ( "LeptonpT_diLep", "Lepton pT (after dilepton cut)", 80, 0, 400 ));
    h_LeptonEta_diLep = store(new TH1D ( "LeptonEta_diLep", "Lepton Eta (after dilepton cut)", 100, -5, 5 ));

    h_AntiLeptonpT = store(new TH1D ( "AntiLeptonpT", "AntiLepton pT", 80, 0, 400 ));
    h_AntiLeptonEta = store(new TH1D ( "AntiLeptonEta", "AntiLepton Eta", 100, -5, 5 ));
    h_AntiLeptonpT_diLep = store(new TH1D ( "AntiLeptonpT_diLep", "Lepton pT (after dilepton cut)", 80, 0, 400 ));
    h_AntiLeptonEta_diLep = store(new TH1D ( "AntiLeptonEta_diLep", "Lepton Eta (after dilepton cut)", 100, -5, 5 ));

    h_HypToppT = store(new TH1D ( "HypToppT", "Top pT", 400, 0, 400 ));
    h_HypTopEta = store(new TH1D ( "HypTopEta", "Top #eta", 100, -5, 5 ));
    h_HypTopMass = store(new TH1D ( "HypTopMass", "Top Mass", 80, 0, 400 ));
    h_HypTopRapidity = store(new TH1D ( "HypTopRapidity", "Top Rapidity", 100, -5, 5 ));
    
    h_HypAntiToppT = store(new TH1D ( "HypAntiToppT", "AntiTop pT", 400, 0, 400 ));
    h_HypAntiTopEta = store(new TH1D ( "HypAntiTopEta", "AntiTop #eta", 100, -5, 5 ));
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

    h_GenRecoMet = store(new TH2D("GenRecoMet", "Missing ET in the event", 500, 0, 500, 500, 0, 500));
    h_VisGenMet = store(new TH1D("VisGenMet", "MET (VisGEN)", 500, 0, 500));
    h_RecoMet = store(new TH1D("RecoMet", "Reconstructed MET", 500, 0, 500));
    h_HypMet = store(new TH1D("HypMet","MET ", 500, 0, 500));
    
    h_GenRecoHT = store(new TH2D("GenRecoHT", "HT in the event", 800, 0, 800, 800, 0, 800));
    h_VisGenHT = store(new TH1D("VisGenHT", "HT (VisGEN)", 800, 0, 800));
    h_RecoHT = store(new TH1D("RecoHT", "Reconstructed HT", 800, 0, 800));
    h_HypHT = store(new TH1D("HypHT", "HT", 800, 0, 800));
    
    h_GenRecoNeutrinopT = store(new TH2D("GenRecoNeutrinopT", "Gen/Reco nu pt", 80, 0, 400, 80, 0, 400));
    h_VisGenNeutrinopT = store(new TH1D("VisGenNeutrinopT", "Nu pT (VisGEN)", 80, 0, 400));
    h_RecoNeutrinopT = store(new TH1D("RecoNeutrinopT", "reco nu pT", 80, 0, 400));
    h_HypNeutrinopT = store(new TH1D("HypNeutrinopT", "hyp nu pT", 80, 0, 400));
    
    h_GenRecoAntiNeutrinopT = store(new TH2D("GenRecoAntiNeutrinopT", "Gen/Reco nubar pt", 80, 0, 400, 80, 0, 400));
    h_VisGenAntiNeutrinopT = store(new TH1D("VisGenAntiNeutrinopT", "Nubar pT (VisGEN)", 80, 0, 400));
    h_RecoAntiNeutrinopT = store(new TH1D("RecoAntiNeutrinopT", "reco nubar pT", 80, 0, 400));
    h_HypAntiNeutrinopT = store(new TH1D("HypAntiNeutrinopT", "hyp nubar pT", 80, 0, 400));
    
    
    h_GenRecoTTBarRapidity = store(new TH2D ( "GenRecoTTBarRapidity", "Rapidity of TTbar System (HYP)", 100, -5, 5, 100, -5, 5 ));
    h_GenRecoTTBarpT = store(new TH2D ( "GenRecoTTBarpT", "pT of TTbar System (HYP)", 500, 0, 500, 500, 0, 500 ));
    h_GenRecoTTBarMass = store(new TH2D ( "GenRecoTTBarMass", "Mass of TTbar System (HYP)", 2000, 0, 2000, 2000, 0, 2000 ));
    h_GenRecoLLBarMass = store(new TH2D ( "GenRecoLLBarMass", "Mass of LLbar System (HYP)", 500, 0, 1000, 500, 0, 1000 ));
    h_GenRecoLLBarpT = store(new TH2D ( "GenRecoLLBarpT", "pT of LLbar System (HYP)", 200, 0, 1000, 200, 0, 1000 ));
    
    h_NJetMatching = store(new TH1D ( "NJetMatching", "NJet Gen/Reco Matching", 5, 0, 5 ));

    h_GenRecoLLBarDPhi = store(new TH2D ( "GenRecoLLBarDPhi", "Gen/Reco Matching", 100, 0., 3.2, 100, 0., 3.2 ));
    h_GenRecoLeptonantiBjetMass = store(new TH2D ( "GenRecoLeptonBjetMass", "Gen/Reco Matching", 500, 0, 1000, 500, 0, 1000 ));
    h_GenRecoAntiLeptonBjetMass = store(new TH2D ( "GenRecoAntiLeptonBjetMass", "Gen/Reco Matching", 500, 0, 1000, 500, 0, 1000 ));
    h_GenRecoJetMult = store(new TH2D ( "GenRecoJetMult", "Gen/Reco Matching", 26, -0.5, 25.5, 26, -0.5, 25.5 ));

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

    CreateBinnedControlPlots(h_HypToppT, h_LeptonpT);
    CreateBinnedControlPlots(h_HypToppT, h_LeptonEta);
    CreateBinnedControlPlots(h_HypToppT, h_MET);
    CreateBinnedControlPlots(h_HypToppT, h_diLepMassFull);
    
    CreateBinnedControlPlots(h_HypTopRapidity, h_LeptonpT);
    CreateBinnedControlPlots(h_HypTopRapidity, h_LeptonEta);
    CreateBinnedControlPlots(h_HypTopRapidity, h_MET);
    CreateBinnedControlPlots(h_HypTopRapidity, h_diLepMassFull);
    
    
    
    h_AllLeptonpT_step0 = store(new TH1D("AllLeptonpT_step0", "p_{T} of all leptons before trigger", 80, 0, 400));
    h_AllLeptonEta_step0 = store(new TH1D("AllLeptonEta_step0", "#eta of all leptons before trigger", 20, -2.5, 2.5));
    h_LeptonpT_step0 = store(new TH1D("LeptonpT_step0", "p_{T} of  leptons before trigger", 80, 0, 400));
    h_LeptonEta_step0 = store(new TH1D("LeptonEta_step0", "#eta of  leptons before trigger", 20, -2.5, 2.5));
    h_LeptonMult_step0 = store(new TH1D("LeptonMult_step0", "Number of all leptons before trigger", 20, -0.5, 20.5));
    h_AllJetspT_step0 = store(new TH1D("AllJetpT_step0", "p_{T} of all jets before trigger", 80, 0, 400));
    h_AllJetsEta_step0 = store(new TH1D("AllJetEta_step0", "#eta of all jets before trigger", 20, -2.5, 2.5));
    h_JetsMult_step0 = store(new TH1D("JetsMult_step0", "Number of the jets before trigger", 20, -0.5, 20.5));
    h_BJetsMult_step0 = store(new TH1D("BJetsMult_step0", "Number of the jets before trigger", 20, -0.5, 20.5));
    
    h_AllLeptonpT_step1 = store(new TH1D("AllLeptonpT_step1", "p_{T} of all leptons before 2 lepton", 80, 0, 400));
    h_AllLeptonEta_step1 = store(new TH1D("AllLeptonEta_step1", "#eta of all leptons before 2 lepton", 20, -2.5, 2.5));
    h_LeptonpT_step1 = store(new TH1D("LeptonpT_step1", "p_{T} of  leptons before 2 lepton", 80, 0, 400));
    h_LeptonEta_step1 = store(new TH1D("LeptonEta_step1", "#eta of  leptons before 2 lepton", 20, -2.5, 2.5));
    h_LeptonMult_step1 = store(new TH1D("LeptonMult_step1", "Number of all leptons before 2 lepton", 20, -0.5, 20.5));
    h_AllJetspT_step1 = store(new TH1D("AllJetpT_step1", "p_{T} of all jets before 2 lepton", 80, 0, 400));
    h_AllJetsEta_step1 = store(new TH1D("AllJetEta_step1", "#eta of all jets before 2 lepton", 20, -2.5, 2.5));
    h_JetspT_step1 = store(new TH1D("JetpT_step1", "p_{T} of  jets before 2 lepton", 80, 0, 400));
    h_JetsEta_step1 = store(new TH1D("JetEta_step1", "#eta of  jets before 2 lepton", 20, -2.5, 2.5));
    h_JetsMult_step1 = store(new TH1D("JetsMult_step1", "Number of the jets before 2 lepton", 20, -0.5, 20.5));
    h_BJetsMult_step1 = store(new TH1D("BJetsMult_step1", "Number of the jets before 2 lepton", 20, -0.5, 20.5));
    
    h_AllLeptonpT_step2 = store(new TH1D("AllLeptonpT_step2", "p_{T} of all leptons before m_{ll}>12 cut", 80, 0, 400));
    h_AllLeptonEta_step2 = store(new TH1D("AllLeptonEta_step2", "#eta of all leptons before m_{ll}>12 cut", 20, -2.5, 2.5));
    h_LeptonpT_step2 = store(new TH1D("LeptonpT_step2", "p_{T} of  leptons before m_{ll}>12 cut", 80, 0, 400));
    h_LeptonEta_step2 = store(new TH1D("LeptonEta_step2", "#eta of  leptons before m_{ll}>12 cut", 20, -2.5, 2.5));
    h_LeptonMult_step2 = store(new TH1D("LeptonMult_step2", "Number of all leptons before m_{ll}>12 cut", 20, -0.5, 20.5));
    h_AllJetspT_step2 = store(new TH1D("AllJetpT_step2", "p_{T} of all jets before m_{ll}>12 cut", 80, 0, 400));
    h_AllJetsEta_step2 = store(new TH1D("AllJetEta_step2", "#eta of all jets before m_{ll}>12 cut", 20, -2.5, 2.5));
    h_JetspT_step2 = store(new TH1D("JetpT_step2", "p_{T} of  jets before m_{ll}>12 cut", 80, 0, 400));
    h_JetsEta_step2 = store(new TH1D("JetEta_step2", "#eta of  jets before m_{ll}>12 cut", 20, -2.5, 2.5));
    h_JetsMult_step2 = store(new TH1D("JetsMult_step2", "Number of the jets before m_{ll}>12 cut", 20, -0.5, 20.5));
    h_BJetsMult_step2 = store(new TH1D("BJetsMult_step2", "Number of the jets before m_{ll}>12 cut", 20, -0.5, 20.5));
    
    h_AllLeptonpT_step3 = store(new TH1D("AllLeptonpT_step3", "p_{T} of all leptons before Zpeakcut", 80, 0, 400));
    h_AllLeptonEta_step3 = store(new TH1D("AllLeptonEta_step3", "#eta of all leptons before Zpeakcut", 20, -2.5, 2.5));
    h_LeptonpT_step3 = store(new TH1D("LeptonpT_step3", "p_{T} of  leptons before Zpeakcut", 80, 0, 400));
    h_LeptonEta_step3 = store(new TH1D("LeptonEta_step3", "#eta of  leptons before Zpeakcut", 20, -2.5, 2.5));
    h_LeptonMult_step3 = store(new TH1D("LeptonMult_step3", "Number of all leptons before Zpeakcut", 20, -0.5, 20.5));
    h_AllJetspT_step3 = store(new TH1D("AllJetpT_step3", "p_{T} of all jets before Zpeakcut", 80, 0, 400));
    h_AllJetsEta_step3 = store(new TH1D("AllJetEta_step3", "#eta of all jets before Zpeakcut", 20, -2.5, 2.5));
    h_JetspT_step3 = store(new TH1D("JetpT_step3", "p_{T} of  jets before Zpeakcut", 80, 0, 400));
    h_JetsEta_step3 = store(new TH1D("JetEta_step3", "#eta of  jets before Zpeakcut", 20, -2.5, 2.5));
    h_JetsMult_step3 = store(new TH1D("JetsMult_step3", "Number of the jets before Zpeakcut", 20, -0.5, 20.5));
    h_BJetsMult_step3 = store(new TH1D("BJetsMult_step3", "Number of the jets before Zpeakcut", 20, -0.5, 20.5));
    
    h_AllLeptonpT_step4 = store(new TH1D("AllLeptonpT_step4", "p_{T} of all leptons before 2jets", 80, 0, 400));
    h_AllLeptonEta_step4 = store(new TH1D("AllLeptonEta_step4", "#eta of all leptons before 2jets", 20, -2.5, 2.5));
    h_LeptonpT_step4 = store(new TH1D("LeptonpT_step4", "p_{T} of  leptons before 2jets", 80, 0, 400));
    h_LeptonEta_step4 = store(new TH1D("LeptonEta_step4", "#eta of  leptons before 2jets", 20, -2.5, 2.5));
    h_LeptonMult_step4 = store(new TH1D("LeptonMult_step4", "Number of all leptons before 2jets", 20, -0.5, 20.5));
    h_AllJetspT_step4 = store(new TH1D("AllJetpT_step4", "p_{T} of all jets before 2jets", 80, 0, 400));
    h_AllJetsEta_step4 = store(new TH1D("AllJetEta_step4", "#eta of all jets before 2jets", 20, -2.5, 2.5));
    h_JetspT_step4 = store(new TH1D("JetpT_step4", "p_{T} of  jets before 2jets", 80, 0, 400));
    h_JetsEta_step4 = store(new TH1D("JetEta_step4", "#eta of  jets before 2jets", 20, -2.5, 2.5));
    h_JetsMult_step4 = store(new TH1D("JetsMult_step4", "Number of the jets before 2jets", 20, -0.5, 20.5));
    h_BJetsMult_step4 = store(new TH1D("BJetsMult_step4", "Number of the jets before 2jets", 20, -0.5, 20.5));
    
    h_AllLeptonpT_step5 = store(new TH1D("AllLeptonpT_step5", "p_{T} of all leptons before MET", 80, 0, 400));
    h_AllLeptonEta_step5 = store(new TH1D("AllLeptonEta_step5", "#eta of all leptons before MET", 20, -2.5, 2.5));
    h_LeptonpT_step5 = store(new TH1D("LeptonpT_step5", "p_{T} of  leptons before MET", 80, 0, 400));
    h_LeptonEta_step5 = store(new TH1D("LeptonEta_step5", "#eta of  leptons before MET", 20, -2.5, 2.5));
    h_LeptonMult_step5 = store(new TH1D("LeptonMult_step5", "Number of all leptons before MET", 20, -0.5, 20.5));
    h_AllJetspT_step5 = store(new TH1D("AllJetpT_step5", "p_{T} of all jets before MET", 80, 0, 400));
    h_AllJetsEta_step5 = store(new TH1D("AllJetEta_step5", "#eta of all jets before MET", 20, -2.5, 2.5));
    h_JetspT_step5 = store(new TH1D("JetpT_step5", "p_{T} of  jets before MET", 80, 0, 400));
    h_JetsEta_step5 = store(new TH1D("JetEta_step5", "#eta of  jets before MET", 20, -2.5, 2.5));
    h_JetsMult_step5 = store(new TH1D("JetsMult_step5", "Number of the jets before MET", 20, -0.5, 20.5));
    h_BJetsMult_step5 = store(new TH1D("BJetsMult_step5", "Number of the jets before MET", 20, -0.5, 20.5));
    
    h_AllLeptonpT_step6 = store(new TH1D("AllLeptonpT_step6", "p_{T} of all leptons before 1b-tag", 80, 0, 400));
    h_AllLeptonEta_step6 = store(new TH1D("AllLeptonEta_step6", "#eta of all leptons before 1b-tag", 20, -2.5, 2.5));
    h_LeptonpT_step6 = store(new TH1D("LeptonpT_step6", "p_{T} of  leptons before 1b-tag", 80, 0, 400));
    h_LeptonEta_step6 = store(new TH1D("LeptonEta_step6", "#eta of  leptons before 1b-tag", 20, -2.5, 2.5));
    h_LeptonMult_step6 = store(new TH1D("LeptonMult_step6", "Number of all leptons before 1b-tag", 20, -0.5, 20.5));
    h_AllJetspT_step6 = store(new TH1D("AllJetpT_step6", "p_{T} of all jets before 1b-tag", 80, 0, 400));
    h_AllJetsEta_step6 = store(new TH1D("AllJetEta_step6", "#eta of all jets before 1b-tag", 20, -2.5, 2.5));
    h_JetspT_step6 = store(new TH1D("JetpT_step6", "p_{T} of  jets before 1b-tag", 80, 0, 400));
    h_JetsEta_step6 = store(new TH1D("JetEta_step6", "#eta of  jets before 1b-tag", 20, -2.5, 2.5));
    h_JetsMult_step6 = store(new TH1D("JetsMult_step6", "Number of the jets before 1b-tag", 20, -0.5, 20.5));
    h_BJetsMult_step6 = store(new TH1D("BJetsMult_step6", "Number of the jets before 1b-tag", 20, -0.5, 20.5));
    
    h_AllLeptonpT_step7 = store(new TH1D("AllLeptonpT_step7", "p_{T} of all leptons before kinReco", 80, 0, 400));
    h_AllLeptonEta_step7 = store(new TH1D("AllLeptonEta_step7", "#eta of all leptons before kinReco", 20, -2.5, 2.5));
    h_LeptonpT_step7 = store(new TH1D("LeptonpT_step7", "p_{T} of  leptons before kinReco", 80, 0, 400));
    h_LeptonEta_step7 = store(new TH1D("LeptonEta_step7", "#eta of  leptons before kinReco", 20, -2.5, 2.5));
    h_LeptonMult_step7 = store(new TH1D("LeptonMult_step7", "Number of all leptons before kinReco", 20, -0.5, 20.5));
    h_AllJetspT_step7 = store(new TH1D("AllJetpT_step7", "p_{T} of all jets before kinReco", 80, 0, 400));
    h_AllJetsEta_step7 = store(new TH1D("AllJetEta_step7", "#eta of all jets before kinReco", 20, -2.5, 2.5));
    h_JetspT_step7 = store(new TH1D("JetpT_step7", "p_{T} of  jets before kinReco", 80, 0, 400));
    h_JetsEta_step7 = store(new TH1D("JetEta_step7", "#eta of  jets before kinReco", 20, -2.5, 2.5));
    h_JetsMult_step7 = store(new TH1D("JetsMult_step7", "Number of the jets before kinReco", 20, -0.5, 20.5));
    h_BJetsMult_step7 = store(new TH1D("BJetsMult_step7", "Number of the jets before kinReco", 20, -0.5, 20.5));
    
    h_AllLeptonpT_step8 = store(new TH1D("AllLeptonpT_step8", "p_{T} of the leptons after kinReco", 80, 0, 400));
    h_AllLeptonEta_step8 = store(new TH1D("AllLeptonEta_step8", "#eta of the leptons after kinReco", 20, -2.5, 2.5));
    h_LeptonpT_step8 = store(new TH1D("LeptonpT_step8", "p_{T} of  leptons after kinReco", 80, 0, 400));
    h_LeptonEta_step8 = store(new TH1D("LeptonEta_step8", "#eta of  leptons after kinReco", 20, -2.5, 2.5));
    h_LeptonMult_step8 = store(new TH1D("LeptonMult_step8", "Number of the leptons after kinReco", 20, -0.5, 20.5));
    h_AllJetspT_step8 = store(new TH1D("AllJetpT_step8", "p_{T} of the jets after kinReco", 80, 0, 400));
    h_AllJetsEta_step8 = store(new TH1D("AllJetEta_step8", "#eta of the jets after kinReco", 20, -2.5, 2.5));
    h_JetspT_step8 = store(new TH1D("JetpT_step8", "p_{T} of jets after kinReco", 80, 0, 400));
    h_JetsEta_step8 = store(new TH1D("JetEta_step8", "#eta of jets after kinReco", 20, -2.5, 2.5));
    h_JetsMult_step8 = store(new TH1D("JetsMult_step8", "Number of the jets after kinReco", 20, -0.5, 20.5));
    h_BJetsMult_step8 = store(new TH1D("BJetsMult_step8", "Number of the jets after kinReco", 20, -0.5, 20.5));
    
    
    //btagSF
    const int PtMax = 11;
    const int EtaMax = 5;
    Double_t ptbins[PtMax+1] = {20.,30.,40.,50.,60.,70.,80.,100.,120.,160.,210.,800.};
    Double_t etabins[EtaMax+1] = {0.0,0.5,1.0,1.5,2.0,2.4};
    
    h_bjets = store(new TH2D("bjets2D", "unTagged Bjets", PtMax, ptbins, EtaMax, etabins));              h_bjets->Sumw2();
    h_btaggedjets = store(new TH2D("bjetsTagged2D", "Tagged Bjets", PtMax, ptbins, EtaMax, etabins));    h_btaggedjets->Sumw2();
    h_cjets = store(new TH2D("cjets2D", "unTagged Cjets", PtMax, ptbins, EtaMax, etabins));              h_cjets->Sumw2();
    h_ctaggedjets = store(new TH2D("cjetsTagged2D", "Tagged Cjets", PtMax, ptbins, EtaMax, etabins));    h_ctaggedjets->Sumw2();
    h_ljets = store(new TH2D("ljets2D", "unTagged Ljets", PtMax, ptbins, EtaMax, etabins));              h_ljets->Sumw2();
    h_ltaggedjets = store(new TH2D("ljetsTagged2D", "Tagged Ljets", PtMax, ptbins, EtaMax, etabins));    h_ltaggedjets->Sumw2();
    
    h_PUSF = store(new TH1D("PUSF", "PU SF per event", 200, 0.5, 1.5));
    h_TrigSF = store(new TH1D("TrigSF", "Trigger SF per event", 200, 0.5, 1.5));
    h_LepSF = store(new TH1D("LepSF", "Lep. Id and Isol. SF per event", 200, 0.75, 1.25));
    h_BTagSF = store(new TH1D("BTagSF", "BTagging SF per event", 200 , 0.95, 1.15 ));
    h_BTagSF->Sumw2();
    h_KinRecoSF = store(new TH1D("KinRecoSF", "Kinematic Reco. SF per event", 200, 0.5, 1.5));
    h_EventWeight = store(new TH1D("EventWeight", "Event SF", 600, 0, 3));

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

void Analysis::orderLVByPt(LV &leading, LV &Nleading, const LV &lv1, const LV &lv2) {
    if (lv1.Pt() > lv2.Pt()) {
        leading = lv1; Nleading = lv2;
    } else {
        leading = lv2; Nleading = lv1;
    }
}

int Analysis::NumberOfBJets(vector<double> *bjets){
    
    int nbjets = 0;
    for (int i=0; i<(int)bjets->size(); i++){
        if(bjets->at(i) > 0.244){nbjets++;}
    }
    return nbjets;
}


Bool_t Analysis::Process ( Long64_t entry )
{
    static const double JETPTCUT = 30;
    
    if ( ++EventCounter % 100000 == 0 ) cout << "Event Counter: " << EventCounter << endl;
    
    //do we have a DY true level cut?
    if (checkZDecayMode && !checkZDecayMode(entry)) return kTRUE;
    
    if (isTtbarPlusTauSample || correctMadgraphBR) b_TopDecayMode->GetEntry(entry);
    //decayMode contains the decay of the top (*10) + the decay of the antitop
    //1=hadron, 2=e, 3=mu, 4=tau->hadron, 5=tau->e, 6=tau->mu
    //i.e. 23 == top decays to e, tbar decays to mu
    if (isTtbarPlusTauSample) {
        bool isViaTau = topDecayMode > 40 || ( topDecayMode % 10 > 4 );
        bool isCorrectChannel = false;
        switch (channelPdgIdProduct) {
            case -11*13: isCorrectChannel = topDecayMode == 23 || topDecayMode == 32; break;
            case -11*11: isCorrectChannel = topDecayMode == 22; break;
            case -13*13: isCorrectChannel = topDecayMode == 33; break;
            default: cerr << "Invalid channel! Product = " << channelPdgIdProduct << "\n";
        };
        bool isBackgroundInSignalSample = !isCorrectChannel || isViaTau;
        if (runViaTau != isBackgroundInSignalSample) return kTRUE;
    }

    GetRecoBranches(entry);
    //We must correct for the madGraph branching fraction being 1/9 for dileptons (PDG average is .108)
    if ( correctMadgraphBR ) {
        if ( topDecayMode == 11 ) { //all hadronic decay
            weightGenerator *= (0.676*1.5) * (0.676*1.5);
        } else if ( topDecayMode< 20 || ( topDecayMode % 10 == 1) ) { //semileptonic Decay
            weightGenerator *= (0.108*9) * (0.676*1.5);
        } else {//dileptonic decay (including taus!)
            weightGenerator *= (0.108*9) * (0.108*9);
        }
    }
    
    if (pdf_no) {
        b_weightPDF->GetEntry(entry);
        weightGenerator *= weightPDF->at(pdf_no - 1); //vector is 0 based
    }
    
//     GetRecoBranches(entry);
    for (size_t i = 0; i < jets->size(); ++i) {
        if (jets->at(i).pt() < JETPTCUT) {
            jets->erase(jets->begin() + i, jets->end());
            jetBTagCSV->erase(jetBTagCSV->begin() + i, jetBTagCSV->end());
            break;
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
    
    if (isSignal) {
        GetSignalBranches(entry);

        std::vector<size_t> idx_leadbHadJet;
        std::vector<size_t> idx_nleadbHadJet;
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

        size_t hadron_index = std::numeric_limits<size_t>::max();
        size_t antihadron_index = std::numeric_limits<size_t>::max();
        size_t hadrontop_index = std::numeric_limits<size_t>::max();
        size_t antihadrontop_index = std::numeric_limits<size_t>::max();
        
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
        if ( *it > BtagWP) {
            //BJetIndex.push_back ( *it );
            BJetIndex.push_back((it-jetBTagCSV->begin())); //change asked by Tyler
        }
    }

    LV LeadGenTop, NLeadGenTop;
    LV LeadGenLepton, NLeadGenLepton;
    LV LeadGenBJet, NLeadGenBJet;
    double genHT = -1;
    
    if ( isSignal ) {
    
        double trueLevelWeight = weightGenerator * weightPU;
        h_GenAll->Fill(GenTop->M(), trueLevelWeight);
        
        //Begin: Select & Fill histograms with Leading pT and 2nd Leading pT: Lepton and BJet
        orderLVByPt(LeadGenLepton, NLeadGenLepton, *GenLepton, *GenAntiLepton);
        
        if (BHadronIndex != -1 && AntiBHadronIndex != -1) {
            orderLVByPt(LeadGenBJet, NLeadGenBJet, 
                        allGenJets->at(BHadronIndex), allGenJets->at(AntiBHadronIndex));
        }
        
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
                    h_VisGenMet->Fill(GenMet->Pt(), trueLevelWeight);
                    
                    //for HT, count only >= JETPTCUT
                    genHT = getJetHT(*allGenJets, JETPTCUT);
                    h_VisGenHT->Fill(genHT, trueLevelWeight);

                    h_VisGenLLBarDPhi->Fill(abs( DeltaPhi(*GenLepton, *GenAntiLepton)), trueLevelWeight );
                    h_VisGenLeptonantiBjetMass->Fill(( *GenLepton + allGenJets->at(AntiBHadronIndex) ).M(), trueLevelWeight );
                    h_VisGenAntiLeptonBjetMass->Fill(( *GenAntiLepton + allGenJets->at(BHadronIndex) ).M(), trueLevelWeight );
                    h_VisGenJetMult->Fill(allGenJets->size(), trueLevelWeight );

                    //Begin: Select & Fill histograms with Leading pT and 2nd Leading pT: Lepton and BJet
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
        
        h_VisGenNeutrinopT->Fill(GenNeutrino->Pt(), trueLevelWeight);
        h_VisGenAntiNeutrinopT->Fill(GenAntiNeutrino->Pt(), trueLevelWeight);
        
        //Begin: Fill histograms with Leading pT and 2nd Leading pT: Top
        orderLVByPt(LeadGenTop, NLeadGenTop, *GenTop, *GenAntiTop);
        h_VisGenToppTLead->Fill(LeadGenTop.Pt(), trueLevelWeight);
        h_VisGenToppTNLead->Fill(NLeadGenTop.Pt(), trueLevelWeight);
        h_VisGenTopRapidityLead->Fill(LeadGenTop.Rapidity(), trueLevelWeight);
        h_VisGenTopRapidityNLead->Fill(NLeadGenTop.Rapidity(), trueLevelWeight);
        h_VisGenTopMassLead->Fill(LeadGenTop.M(), trueLevelWeight);
        h_VisGenTopMassNLead->Fill(NLeadGenTop.M(), trueLevelWeight);
        //End: Fill histograms with Leading pT and 2nd Leading pT: Top
        
    }//for visible top events

    // ++++ Control Plots ++++
    for (int i=0; i<(int) leptons->size(); ++i){
        h_AllLeptonEta_step0->Fill(leptons->at(i).Eta(), 1);
        h_AllLeptonpT_step0->Fill(leptons->at(i).Pt(), 1);
    }
    
    for (int i=0; i<(int) jets->size(); ++i){
        h_AllJetsEta_step0->Fill(jets->at(i).Eta(), 1);
        h_AllJetspT_step0->Fill(jets->at(i).Pt(), 1);
    }
    h_LeptonMult_step0->Fill(leptons->size(), 1);
    h_JetsMult_step0->Fill(jets->size(), 1);
    int nbjets_step0 = NumberOfBJets(jetBTagCSV);
    h_BJetsMult_step0->Fill(nbjets_step0, 1);
    
    h_PUSF->Fill(weightPU, 1);

    
    //===CUT===
    // check if event was triggered
    //our triggers (bits: see the ntuplewriter!)    
    constexpr int mumuTriggers = 0x8 + 0x20; //17/8 + 17Tr8
    constexpr int emuTriggers = 0x2000 + 0x4000;
    constexpr int eeTriggers = 0x40000;
    
    if (!(((triggerBits & mumuTriggers) && channelPdgIdProduct == -13*13)    //mumu triggers in rightmost byte
          ||((triggerBits & emuTriggers) && channelPdgIdProduct == -11*13)     //emu in 2nd byte
          ||((triggerBits & eeTriggers) && channelPdgIdProduct == -11*11)))    //ee in 3rd byte
    {
        return kTRUE;
    }

    size_t LeadLeptonNumber = 0;
    size_t NLeadLeptonNumber = 0;
    bool hasLeptonPair = getLeptonPair(LeadLeptonNumber, NLeadLeptonNumber);

    // ++++ Control Plots ++++
    for (int i=0; i<(int) leptons->size(); ++i){
        h_AllLeptonEta_step1->Fill(leptons->at(i).Eta(), 1);
        h_AllLeptonpT_step1->Fill(leptons->at(i).Pt(), 1);
    }
    if(leptons->size()>1){
        h_LeptonEta_step1->Fill(leptons->at(0).Eta(), 1);
        h_LeptonpT_step1->Fill(leptons->at(0).Pt(), 1);
        h_LeptonEta_step1->Fill(leptons->at(1).Eta(), 1);
        h_LeptonpT_step1->Fill(leptons->at(1).Pt(), 1);
    }

    for (int i=0; i<(int) jets->size(); ++i){
        h_AllJetsEta_step1->Fill(jets->at(i).Eta(), 1);
        h_AllJetspT_step1->Fill(jets->at(i).Pt(), 1);
    }
    h_LeptonMult_step1->Fill(leptons->size(), 1);
    h_JetsMult_step1->Fill(jets->size(), 1);
    int nbjets_step1 = NumberOfBJets(jetBTagCSV);
    h_BJetsMult_step1->Fill(nbjets_step1, 1);
    
    
    //===CUT===
    // we need an OS lepton pair
    if (! hasLeptonPair) return kTRUE;
    // lepton pt cut > 20, sufficient to test the second leading lepton
    if (leptons->at(NLeadLeptonNumber).pt() <= 20) return kTRUE; 
    
    LV dilepton = leptons->at(LeadLeptonNumber) + leptons->at(NLeadLeptonNumber);
    
    // ++++ Control Plots ++++
    for (int i=0; i<(int) leptons->size(); ++i){
        h_AllLeptonEta_step2->Fill(leptons->at(i).Eta(), 1);
        h_AllLeptonpT_step2->Fill(leptons->at(i).Pt(), 1);
    }
    h_LeptonEta_step2->Fill(leptons->at(0).Eta(), 1);
    h_LeptonpT_step2->Fill(leptons->at(0).Pt(), 1);
    h_LeptonEta_step2->Fill(leptons->at(1).Eta(), 1);
    h_LeptonpT_step2->Fill(leptons->at(1).Pt(), 1);

    for (int i=0; i<(int) jets->size(); ++i){
        h_AllJetsEta_step2->Fill(jets->at(i).Eta(), 1);
        h_AllJetspT_step2->Fill(jets->at(i).Pt(), 1);
    }
    if(jets->size()>1){
        h_JetsEta_step2->Fill(jets->at(0).Eta(), 1);
        h_JetspT_step2->Fill(jets->at(0).Pt(), 1);
        h_JetsEta_step2->Fill(jets->at(1).Eta(), 1);
        h_JetspT_step2->Fill(jets->at(1).Pt(), 1);
    }
    h_LeptonMult_step2->Fill(leptons->size(), 1);
    h_JetsMult_step2->Fill(jets->size(), 1);
    int nbjets_step2 = NumberOfBJets(jetBTagCSV);
    h_BJetsMult_step2->Fill(nbjets_step2, 1);
    

    //===CUT===
    //with at least 12 GeV invariant mass
    if (dilepton.M() < 20) return kTRUE;
    
    // find l+ and l-
    LV leptonPlus;
    LV leptonMinus;
    if (lepPdgId->at(LeadLeptonNumber) < 0) {
        leptonPlus = leptons->at(LeadLeptonNumber);
        leptonMinus = leptons->at(NLeadLeptonNumber);
    } else {
        leptonMinus = leptons->at(LeadLeptonNumber);
        leptonPlus = leptons->at(NLeadLeptonNumber);
    }

    //This is necessary due to the ordering in the trigger 2D-plots
    LV leptonX, leptonY;
    if ( abs( lepPdgId->at(LeadLeptonNumber) ) == abs( lepPdgId->at(NLeadLeptonNumber) ) ){
        //in ee and mumu channel leptonX must be the highest pT lepton
        leptonX = leptons->at(LeadLeptonNumber);
        leptonY = leptons->at(NLeadLeptonNumber);
    } else {
        // in emu channel lepX should be electron
        if (abs(lepPdgId->at(LeadLeptonNumber)) == 11) {
            leptonX = leptons->at(LeadLeptonNumber);
            leptonY = leptons->at(NLeadLeptonNumber);
        } else {
            leptonX = leptons->at(NLeadLeptonNumber);
            leptonY = leptons->at(LeadLeptonNumber);
        }
    }
    
    int LleptonId = lepPdgId->at(LeadLeptonNumber);
    int NLleptonId = lepPdgId->at(NLeadLeptonNumber);
    //Now determine the lepton trigger and ID scale factors
    double weightLepSF = isMC ? getLeptonIDSF(leptons->at(LeadLeptonNumber), leptons->at(NLeadLeptonNumber), LleptonId, NLleptonId) : 1;  
    double weightTrigSF = isMC ? getTriggerSF(leptonX, leptonY) : 1;
    
    //First control plots after dilepton selection (without Z cut)
    double weight = weightGenerator*weightTrigSF*weightLepSF;

    //weight even without PU reweighting
    h_vertMulti_noPU->Fill(vertMulti, weight);
    
    //apply PU reweighting - continue with control plots
    weight *= weightPU;
    h_vertMulti->Fill(vertMulti, weight);
    
    h_step4->Fill(1, weight);
    h_TrigSF->Fill(weightTrigSF, 1);
    h_LepSF->Fill(weightLepSF, 1);
    
    h_jetMulti_diLep->Fill(jets->size(), weight);
    h_diLepMassFull->Fill(dilepton.M(), weight);

    
    //****************************************
    //handle inverted Z cut
    // Fill loose dilepton mass histogram before any jet cuts
    bool isZregion = dilepton.M() > 76 && dilepton.M() < 106;
    bool hasJets = jets->size() > 1 && jets->at(1).Pt() > JETPTCUT;
    bool hasMetOrEmu = channel == "emu" || met->Pt() > 40;
    bool hasBtag = BJetIndex.size() > 0;
    double weightKinFit = isMC ? 0.9866 : 1;
    double weightBtagSF = -1; //trick: initialize to -1 to avoid calculation of the btagSF twice
    
    bool hasSolution = HypTop->size() > 0;
    if (kinRecoOnTheFly || true) 
        hasSolution = calculateKinReco(leptonMinus, leptonPlus, JETPTCUT);
    
    if ( isZregion ) {
        double fullWeights = weightGenerator*weightPU*weightTrigSF*weightLepSF;
        Looseh1->Fill(dilepton.M(), fullWeights);
        Zh1_postZcut->Fill(dilepton.M(), fullWeights);
        Allh1_postZcut->Fill(dilepton.M(), fullWeights);
        
        if ( hasJets ) {
            Zh1_post2jets->Fill(dilepton.M(), fullWeights);
            Allh1_post2jets->Fill(dilepton.M(), fullWeights);
            
            if ( hasMetOrEmu ) {
                Zh1_postMET->Fill(dilepton.M(), fullWeights);
                Allh1_postMET->Fill(dilepton.M(), fullWeights);
                
                if ( hasBtag ) {
                    weightBtagSF = isMC ? calculateBtagSF() : 1;
                    fullWeights *= weightBtagSF;
                    Zh1_post1btag->Fill(dilepton.M(), fullWeights);
                    Allh1_post1btag->Fill(dilepton.M(), fullWeights);
                    
                    if ( hasSolution ) {
                        fullWeights *= weightKinFit;
                        Zh1_postKinReco->Fill(dilepton.M(), fullWeights);
                        Allh1_postKinReco->Fill(dilepton.M(), fullWeights);
                    }
                }
            }
        }
    }
    
    // ++++ Control Plots ++++
    for (int i=0; i<(int) leptons->size(); ++i){
        h_AllLeptonEta_step3->Fill(leptons->at(i).Eta(), weight);
        h_AllLeptonpT_step3->Fill(leptons->at(i).Pt(), weight);
    }
    h_LeptonEta_step3->Fill(leptons->at(0).Eta(), weight);
    h_LeptonpT_step3->Fill(leptons->at(0).Pt(), weight);
    h_LeptonEta_step3->Fill(leptons->at(1).Eta(), weight);
    h_LeptonpT_step3->Fill(leptons->at(1).Pt(), weight);

    for (int i=0; i<(int) jets->size(); ++i){
        h_AllJetsEta_step3->Fill(jets->at(i).Eta(), weight);
        h_AllJetspT_step3->Fill(jets->at(i).Pt(), weight);
    }
    if(jets->size()>1){
        h_JetsEta_step3->Fill(jets->at(0).Eta(), weight);
        h_JetspT_step3->Fill(jets->at(0).Pt(), weight);
        h_JetsEta_step3->Fill(jets->at(1).Eta(), weight);
        h_JetspT_step3->Fill(jets->at(1).Pt(), weight);
    }
    h_LeptonMult_step3->Fill(leptons->size(), weight);
    h_JetsMult_step3->Fill(jets->size(), weight);
    int nbjets_step3 = NumberOfBJets(jetBTagCSV);
    h_BJetsMult_step3->Fill(nbjets_step3, weight);
    
    //=== CUT ===
    //Exclude the Z window
    if (channel != "emu" && isZregion) return kTRUE;
    
    h_step5->Fill(1, weight);
    h_LeptonpT_diLep->Fill(leptonMinus.Pt(), weight);
    h_AntiLeptonpT_diLep->Fill(leptonPlus.Pt(), weight);
    h_LeptonEta_diLep->Fill(leptonMinus.Eta(), weight);
    h_AntiLeptonEta_diLep->Fill(leptonPlus.Eta(), weight);
    
    h_MET->Fill(met->Pt(), weight);
    //loop over both leptons
    for (auto i : {LeadLeptonNumber, NLeadLeptonNumber}) {
        if ( std::abs(lepPdgId->at(i)) == 11 ) {
            h_ElectronpT->Fill(leptons->at(i).Pt(), weight);
            h_ElectronEta->Fill(leptons->at(i).Eta(), weight);
        }
        if ( std::abs(lepPdgId->at(i)) == 13 ) {
            h_MuonpT->Fill(leptons->at(i).Pt(), weight);
            h_MuonEta->Fill(leptons->at(i).Eta(), weight);
        }
    }
    
    // ++++ Control Plots ++++
    for (int i=0; i<(int) leptons->size(); ++i){
        h_AllLeptonEta_step4->Fill(leptons->at(i).Eta(), weight);
        h_AllLeptonpT_step4->Fill(leptons->at(i).Pt(), weight);
    }
    h_LeptonEta_step4->Fill(leptons->at(0).Eta(), weight);
    h_LeptonpT_step4->Fill(leptons->at(0).Pt(), weight);
    h_LeptonEta_step4->Fill(leptons->at(1).Eta(), weight);
    h_LeptonpT_step4->Fill(leptons->at(1).Pt(), weight);

    for (int i=0; i<(int) jets->size(); ++i){
        h_AllJetsEta_step4->Fill(jets->at(i).Eta(), weight);
        h_AllJetspT_step4->Fill(jets->at(i).Pt(), weight);
    }
    if(jets->size()>1){
        h_JetsEta_step4->Fill(jets->at(0).Eta(), weight);
        h_JetspT_step4->Fill(jets->at(0).Pt(), weight);
        h_JetsEta_step4->Fill(jets->at(1).Eta(), weight);
        h_JetspT_step4->Fill(jets->at(1).Pt(), weight);
    }
    h_LeptonMult_step4->Fill(leptons->size(), weight);
    h_JetsMult_step4->Fill(jets->size(), weight);
    int nbjets_step4 = NumberOfBJets(jetBTagCSV);
    h_BJetsMult_step4->Fill(nbjets_step4, weight);
    
    if (!isZregion) { //also apply Z cut in emu!
        TTh1_postZcut->Fill(dilepton.M(), weight);
        Allh1_postZcut->Fill(dilepton.M(), weight);  //this is also filled in the Z region in the code above
    }
    
    //=== CUT ===
    //Require at least two jets > 30 GeV (check for > 30 needed because we might have 20 GeV jets in our NTuple)
    if (! hasJets) return kTRUE;
    h_step6->Fill(1, weight);
    
    // ++++ Control Plots ++++
    for (int i=0; i<(int) leptons->size(); ++i){
        h_AllLeptonEta_step5->Fill(leptons->at(i).Eta(), weight);
        h_AllLeptonpT_step5->Fill(leptons->at(i).Pt(), weight);
    }
    h_LeptonEta_step5->Fill(leptons->at(0).Eta(), weight);
    h_LeptonpT_step5->Fill(leptons->at(0).Pt(), weight);
    h_LeptonEta_step5->Fill(leptons->at(1).Eta(), weight);
    h_LeptonpT_step5->Fill(leptons->at(1).Pt(), weight);

    for (int i=0; i<(int) jets->size(); ++i){
        h_AllJetsEta_step5->Fill(jets->at(i).Eta(), weight);
        h_AllJetspT_step5->Fill(jets->at(i).Pt(), weight);
    }
    if(jets->size()>1){
        h_JetsEta_step5->Fill(jets->at(0).Eta(), weight);
        h_JetspT_step5->Fill(jets->at(0).Pt(), weight);
        h_JetsEta_step5->Fill(jets->at(1).Eta(), weight);
        h_JetspT_step5->Fill(jets->at(1).Pt(), weight);
    }
    h_LeptonMult_step5->Fill(leptons->size(), weight);
    h_JetsMult_step5->Fill(jets->size(), weight);
    int nbjets_step5 = NumberOfBJets(jetBTagCSV);
    h_BJetsMult_step5->Fill(nbjets_step5, weight);
    
    if (!isZregion) { //also apply Z cut in emu!
        TTh1_post2jets->Fill(dilepton.M(), weight);
        Allh1_post2jets->Fill(dilepton.M(), weight);  //this is also filled in the Z region in the code above
    }
    
    //=== CUT ===
    //Require MET > 30 GeV in non-emu channels
    if (!hasMetOrEmu) return kTRUE;
    h_step7->Fill(1, weight);
 
    h_LeptonpT_postMETcut->Fill(leptonMinus.Pt(), weight);
    h_AntiLeptonpT_postMETcut->Fill(leptonPlus.Pt(), weight);
    h_LeptonEta_postMETcut->Fill(leptonMinus.Eta(), weight);
    h_AntiLeptonEta_postMETcut->Fill(leptonPlus.Eta(), weight);

    //loop over both leptons
    for (auto i : {LeadLeptonNumber, NLeadLeptonNumber}) {
        if ( std::abs(lepPdgId->at(i)) == 11 ) {
            h_ElectronpT_postMETcut->Fill(leptons->at(i).Pt(), weight);
            h_ElectronEta_postMETcut->Fill(leptons->at(i).Eta(), weight);
        }
        if ( std::abs(lepPdgId->at(i)) == 13 ) {
            h_MuonpT_postMETcut->Fill(leptons->at(i).Pt(), weight);
            h_MuonEta_postMETcut->Fill(leptons->at(i).Eta(), weight);
        }
    }
    
    h_jetMulti_noBTag->Fill(jets->size(), weight);
    h_BjetMulti_noBTag->Fill(BJetIndex.size(), weight);

    double jetHT = getJetHT(*jets, JETPTCUT);
    h_jetHT->Fill(jetHT, weight);

    for ( size_t i = 0; i < 2; ++i ) {
        h_jetpT->Fill(jets->at(i).Pt(), weight);
    }

    // ++++ Control Plots ++++
    for (int i=0; i<(int) leptons->size(); ++i){
        h_AllLeptonEta_step6->Fill(leptons->at(i).Eta(), weight);
        h_AllLeptonpT_step6->Fill(leptons->at(i).Pt(), weight);
    }
    h_LeptonEta_step6->Fill(leptons->at(0).Eta(), weight);  
    h_LeptonpT_step6->Fill(leptons->at(0).Pt(), weight);    
    h_LeptonEta_step6->Fill(leptons->at(1).Eta(), weight);
    h_LeptonpT_step6->Fill(leptons->at(1).Pt(), weight);

    for (int i=0; i<(int) jets->size(); ++i){
        h_AllJetsEta_step6->Fill(jets->at(i).Eta(), weight);
        h_AllJetspT_step6->Fill(jets->at(i).Pt(), weight);
    }
    if(jets->size()>1){
        h_JetsEta_step6->Fill(jets->at(0).Eta(), weight);
        h_JetspT_step6->Fill(jets->at(0).Pt(), weight);
        h_JetsEta_step6->Fill(jets->at(1).Eta(), weight);
        h_JetspT_step6->Fill(jets->at(1).Pt(), weight);
    }
    h_LeptonMult_step6->Fill(leptons->size(), weight);
    h_JetsMult_step6->Fill(jets->size(), weight);
    int nbjets_step6 = NumberOfBJets(jetBTagCSV);
    h_BJetsMult_step6->Fill(nbjets_step6, weight);
    
    if (!isZregion) { //also apply Z cut in emu!
        TTh1_postMET->Fill(dilepton.M(), weight);
        Allh1_postMET->Fill(dilepton.M(), weight);  //this is also filled in the Z region in the code above
    }
    
    //=== CUT ===
    //Require at least one b tagged jet
    if (!hasBtag) return kTRUE;

    //finally do the btag efficiency calculation stuff
    if(isSignal){
        for (size_t i = 0; i < jets->size(); ++i) {
            if (jets->at(i).Pt() <= JETPTCUT) break;
            double absJetEta = abs(jets->at(i).Eta());
            if (absJetEta<2.4) {
                int partonFlavour = abs(jetPartonFlavour->at(i));
                if(partonFlavour == 5){//b-quark
                    h_bjets->Fill(jets->at(i).Pt(), absJetEta);
                    if((*jetBTagCSV)[i]>BtagWP){
                        h_btaggedjets->Fill(jets->at(i).Pt(), absJetEta);
                    }
                }
                else if (partonFlavour == 4){//c-quark
                    h_cjets->Fill(jets->at(i).Pt(), absJetEta);
                    if((*jetBTagCSV)[i]>BtagWP){
                        h_ctaggedjets->Fill(jets->at(i).Pt(), absJetEta);
                    }
                }
                else if (partonFlavour != 0){//l-quark
                    h_ljets->Fill(jets->at(i).Pt(), absJetEta);
                    if((*jetBTagCSV)[i]>BtagWP){
                        h_ltaggedjets->Fill(jets->at(i).Pt(), absJetEta);
                    }
                }
            }
        }
    }
    
    if (weightBtagSF == -1) weightBtagSF = isMC ? calculateBtagSF() : 1; //avoid calculation of the btagSF twice
    weight *= weightBtagSF;
    h_BTagSF->Fill(weightBtagSF );
    h_step8->Fill(1, weight );

    h_BjetMulti->Fill(BJetIndex.size(), weight);
    h_jetMulti->Fill(jets->size(), weight);
    
    h_leptonPtBeforeKinReco->Fill(leptonMinus.Pt(), weight);
    h_leptonPtBeforeKinReco->Fill(leptonPlus.Pt(), weight);
    h_leptonEtaBeforeKinReco->Fill(leptonMinus.Eta(), weight);
    h_leptonEtaBeforeKinReco->Fill(leptonPlus.Eta(), weight);

    // ++++ Control Plots ++++
    for (int i=0; i<(int) leptons->size(); ++i){
        h_AllLeptonEta_step7->Fill(leptons->at(i).Eta(), weight);
        h_AllLeptonpT_step7->Fill(leptons->at(i).Pt(), weight);
    }
    h_LeptonEta_step7->Fill(leptons->at(0).Eta(), weight);
    h_LeptonpT_step7->Fill(leptons->at(0).Pt(), weight);
    h_LeptonEta_step7->Fill(leptons->at(1).Eta(), weight);
    h_LeptonpT_step7->Fill(leptons->at(1).Pt(), weight);

    for (int i=0; i<(int) jets->size(); ++i){
        h_AllJetsEta_step7->Fill(jets->at(i).Eta(), weight);
        h_AllJetspT_step7->Fill(jets->at(i).Pt(), weight);
    }
    if(jets->size()>1){
        h_JetsEta_step7->Fill(jets->at(0).Eta(), weight);
        h_JetspT_step7->Fill(jets->at(0).Pt(), weight);
        h_JetsEta_step7->Fill(jets->at(1).Eta(), weight);
        h_JetspT_step7->Fill(jets->at(1).Pt(), weight);
    }
    h_LeptonMult_step7->Fill(leptons->size(), weight);
    h_JetsMult_step7->Fill(jets->size(), weight);
    int nbjets_step7 = NumberOfBJets(jetBTagCSV);
    h_BJetsMult_step7->Fill(nbjets_step7, weight);
    
    if (!isZregion) { //also apply Z cut in emu!
        TTh1_post1btag->Fill(dilepton.M(), weight);
        Allh1_post1btag->Fill(dilepton.M(), weight);  //this is also filled in the Z region in the code above
    }
    
    //=== CUT ===
    //Require at least one solution for the kinematic event reconstruction
    if (!hasSolution) return kTRUE;
    weight *= weightKinFit;
    
    h_KinRecoSF->Fill(weightKinFit, 1);
    h_EventWeight->Fill(weight, 1);
    
     // ++++ Control Plots ++++
    for (int i=0; i<(int) leptons->size(); ++i){
        h_AllLeptonEta_step8->Fill(leptons->at(i).Eta(), weight);
        h_AllLeptonpT_step8->Fill(leptons->at(i).Pt(), weight);
    }
    h_LeptonEta_step8->Fill(leptons->at(0).Eta(), weight);
    h_LeptonpT_step8->Fill(leptons->at(0).Pt(), weight);
    h_LeptonEta_step8->Fill(leptons->at(1).Eta(), weight);
    h_LeptonpT_step8->Fill(leptons->at(1).Pt(), weight);

    for (int i=0; i<(int) jets->size(); ++i){
        h_AllJetsEta_step8->Fill(jets->at(i).Eta(), weight);
        h_AllJetspT_step8->Fill(jets->at(i).Pt(), weight);
    }
    if(jets->size()>1){
        h_JetsEta_step8->Fill(jets->at(0).Eta(), weight);
        h_JetspT_step8->Fill(jets->at(0).Pt(), weight);
        h_JetsEta_step8->Fill(jets->at(1).Eta(), weight);
        h_JetspT_step8->Fill(jets->at(1).Pt(), weight);
    }
    h_LeptonMult_step8->Fill(leptons->size(), weight);
    h_JetsMult_step8->Fill(jets->size(), weight);
    int nbjets_step8 = NumberOfBJets(jetBTagCSV);
    h_BJetsMult_step8->Fill(nbjets_step8, weight);
    
    h_step9->Fill(1, weight);
    h_jetMultiXSec->Fill(jets->size(), weight);
    h_jetMultiNoPU->Fill(jets->size(), weight / weightPU );
    h_diLepMassFull_fullSel->Fill(dilepton.M(), weight);
    
    h_leptonPtAfterKinReco->Fill(leptonMinus.Pt(), weight);
    h_leptonPtAfterKinReco->Fill(leptonPlus.Pt(), weight);
    h_leptonEtaAfterKinReco->Fill(leptonMinus.Eta(), weight);
    h_leptonEtaAfterKinReco->Fill(leptonPlus.Eta(), weight);
    
    //create helper variables
    
    //Begin: find 1st (and 2nd) leading pT particles: Top, Lepton, BJetIndex
    LV LeadHypTop, NLeadHypTop;
    LV LeadHypLepton, NLeadHypLepton;
    LV LeadHypBJet, NLeadHypBJet;
    
    size_t solutionIndex = 0; //always zero!
    orderLVByPt(LeadHypTop, NLeadHypTop, HypTop->at(solutionIndex), HypAntiTop->at(solutionIndex));
    orderLVByPt(LeadHypLepton, NLeadHypLepton, HypLepton->at(solutionIndex), HypAntiLepton->at(solutionIndex));
    orderLVByPt(LeadHypBJet, NLeadHypBJet, HypBJet->at(solutionIndex), HypAntiBJet->at(solutionIndex));
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
    
    h_RecoMet->Fill(met->Pt(), recoWeight);
    h_RecoHT->Fill(jetHT, recoWeight);
    
    h_RecoNeutrinopT->Fill(HypNeutrino->at(solutionIndex).Pt(), recoWeight);
    h_RecoAntiNeutrinopT->Fill(HypAntiNeutrino->at(solutionIndex).Pt(), recoWeight);

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
    
    h_HypMet->Fill(met->Pt(), weight);
    h_HypHT->Fill(jetHT, weight);

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
    
    h_HypNeutrinopT->Fill(HypNeutrino->at(solutionIndex).Pt(), weight);
    h_HypAntiNeutrinopT->Fill(HypAntiNeutrino->at(solutionIndex).Pt(), weight);
    
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

    //make sure you have called CreateBinnedControlPlots in the SlaveBegin first
    for (const auto& i : { HypTop->at(solutionIndex), HypAntiTop->at(solutionIndex) } ) {
        FillBinnedControlPlot(h_HypToppT, i.Pt(), h_LeptonpT, leptonMinus.Pt(), weight);
        FillBinnedControlPlot(h_HypToppT, i.Pt(), h_LeptonpT, leptonPlus.Pt(), weight);
        FillBinnedControlPlot(h_HypToppT, i.Pt(), h_diLepMassFull, dilepton.M(), weight);
        FillBinnedControlPlot(h_HypToppT, i.Pt(), h_LeptonEta, leptonMinus.Eta(), weight);
        FillBinnedControlPlot(h_HypToppT, i.Pt(), h_LeptonEta, leptonPlus.Eta(), weight);
        FillBinnedControlPlot(h_HypToppT, i.Pt(), h_MET, met->Pt(), weight);
        
        FillBinnedControlPlot(h_HypTopRapidity, i.Rapidity(), h_LeptonpT, leptonMinus.Pt(), weight);
        FillBinnedControlPlot(h_HypTopRapidity, i.Rapidity(), h_LeptonpT, leptonPlus.Pt(), weight);
        FillBinnedControlPlot(h_HypTopRapidity, i.Rapidity(), h_diLepMassFull, dilepton.M(), weight);
        FillBinnedControlPlot(h_HypTopRapidity, i.Rapidity(), h_LeptonEta, leptonMinus.Eta(), weight);
        FillBinnedControlPlot(h_HypTopRapidity, i.Rapidity(), h_LeptonEta, leptonPlus.Eta(), weight);
        FillBinnedControlPlot(h_HypTopRapidity, i.Rapidity(), h_MET, met->Pt(), weight);
    }

    if (!isZregion) { //also apply Z cut in emu!
        TTh1_postKinReco->Fill(dilepton.M(), weight);
        Allh1_postKinReco->Fill(dilepton.M(), weight);  //this is also filled in the Z region in the code above
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
    h_GenRecoMet->Fill(met->Pt(), GenMet->Pt(), weight);
    h_GenRecoHT->Fill(jetHT, genHT, weight);
    h_GenRecoNeutrinopT->Fill(HypNeutrino->at(solutionIndex).Pt(), GenNeutrino->Pt(), weight);
    h_GenRecoAntiNeutrinopT->Fill(HypAntiNeutrino->at(solutionIndex).Pt(), GenAntiNeutrino->Pt(), weight);

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
    
    if(BHadronIndex >= 0 && AntiBHadronIndex >= 0){
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
        h_GenRecoJetMult->Fill(jets->size(), allGenJets->size(), weight );
    } else {
        h_GenRecoJetMult->Fill(jets->size(), -1000., weight );
    }

    LV genllbar(*GenLepton + *GenAntiLepton);
    h_GenRecoLLBarMass->Fill(hypllbar.M(), genllbar.M(), weight );
    h_GenRecoLLBarpT->Fill(hypllbar.Pt(), genllbar.Pt(), weight );

    LV genttbar(*GenTop + *GenAntiTop);
    h_GenRecoTTBarMass->Fill(hypttbar.M(), genttbar.M(), weight );
    h_GenRecoTTBarpT->Fill(hypttbar.Pt(), genttbar.Pt(), weight );
    h_GenRecoTTBarRapidity->Fill(hypttbar.Rapidity(), genttbar.Rapidity(), weight );

    return kTRUE;
}

void Analysis::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.
    for (auto it = binnedControlPlots->begin(); it != binnedControlPlots->end(); ++it) {
        delete (*it).second.first;
    }
    delete binnedControlPlots;
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
    cout<<"Created: \033[1;1m"<<f_savename<<"\033[1;m\n\n";
    
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
        cout << "Done.\n\n\n";
    }
    fOutput->SetOwner();
    fOutput->Clear();
}

double Analysis::BJetSF( double pt, double eta )
{
    //CSVL b-jet SF
    //From BTV-11-004 and https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-mujet_payload.tptt (ICHEP 2012 prescription)
    //From: https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-pt_payload_Moriond13.tptt  (Moriond 2013 prescription)

    if ( abs(eta) > 2.4 ) {
        cout<<"Jet Eta="<<eta<<" Out of the selected range (|eta|<2.4). Check it"<<endl;
        exit(9);
    }
    if ( pt < 30 ){
        cout<<"Jet pT="<<pt<<" Out of the selected range (pt>30GeV). Check it"<<endl;
        exit(9);
    }

    if ( pt < 20 ) pt = 20;
    if ( pt > 800 ) pt = 800;

    return 0.981149*((1.+(-0.000713295*pt))/(1.+(-0.000703264*pt)));
}

double Analysis::CJetSF ( double pt, double eta )
{
    //CSVL c-jet SF
    //From BTV-11-004 and https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC (ICHEP 2012 prescription)
    //From https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_Moriond13_presc  (Moriond 2013 prescription)

    return 2 * BJetSF( pt, eta );
}

double Analysis::LJetSF ( double pt, double eta, TString typevar )
{
    //CSVL ligth jet mistag SF. Includes also the SF for variations up and down
    //From BTV-11-004 and https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlightFuncs.C (ICHEP 2012 prescription)
    //From https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlightFuncs_Moriond2013.C  (Moriond 2013 prescription)

    double eta_abs = abs(eta);
    if ( eta_abs > 2.4 ) {
        cout<<"Jet Eta="<<eta<<" Out of the selected range (|eta|<2.4). Check it"<<endl;
        exit(91);
    }
    if ( pt < 30 ){
        cout<<"Jet pT="<<pt<<" out of the selected range. Check it"<<endl;
        exit(91);
    }

    if ( pt > 800 && eta_abs <= 1.5 ) {
        pt = 800;
    }    else if (pt > 700 && eta_abs > 1.5 ) {
        pt = 700;
    }

    if ( typevar == "central" ){
        if ( eta_abs <= 0.5 ) {
            return ((1.04901+(0.00152181*pt))+(-3.43568e-06*(pt*pt)))+(2.17219e-09*(pt*(pt*pt)));
        } else if ( eta_abs <= 1.0 ) {
            return ((0.991915+(0.00172552*pt))+(-3.92652e-06*(pt*pt)))+(2.56816e-09*(pt*(pt*pt)));
        } else if ( eta_abs <= 1.5 ) {
            return ((0.962127+(0.00192796*pt))+(-4.53385e-06*(pt*pt)))+(3.0605e-09*(pt*(pt*pt)));
        } else {
            return ((1.06121+(0.000332747*pt))+(-8.81201e-07*(pt*pt)))+(7.43896e-10*(pt*(pt*pt)));
        }
    }  else if ( typevar == "up" ){
        if ( eta_abs <= 0.5 ) {
            return ((1.12424+(0.00201136*pt))+(-4.64021e-06*(pt*pt)))+(2.97219e-09*(pt*(pt*pt)));
        } else if ( eta_abs <= 1.0 ) {
            return ((1.06231+(0.00215815*pt))+(-4.9844e-06*(pt*pt)))+(3.27623e-09*(pt*(pt*pt)));
        } else if ( eta_abs <= 1.5 ) {
            return ((1.02883+(0.00231985*pt))+(-5.57924e-06*(pt*pt)))+(3.81235e-09*(pt*(pt*pt)));
        } else {
            return ((1.1388+(0.000468418*pt))+(-1.36341e-06*(pt*pt)))+(1.19256e-09*(pt*(pt*pt)));
        }
    } else if ( typevar == "down" ){
        if ( eta_abs <= 0.5 ) {
            return ((0.973773+(0.00103049*pt))+(-2.2277e-06*(pt*pt)))+(1.37208e-09*(pt*(pt*pt)));
        } else if ( eta_abs <= 1.0 ) {
            return ((0.921518+(0.00129098*pt))+(-2.86488e-06*(pt*pt)))+(1.86022e-09*(pt*(pt*pt)));
        } else if ( eta_abs <= 1.5 ) {
            return ((0.895419+(0.00153387*pt))+(-3.48409e-06*(pt*pt)))+(2.30899e-09*(pt*(pt*pt)));
        } else {
            return ((0.983607+(0.000196747*pt))+(-3.98327e-07*(pt*pt)))+(2.95764e-10*(pt*(pt*pt)));
        }
    } else {
        cout<<"Type of variation not valid. Check it"<<endl;
        exit(92);
    }

}

double Analysis::BJetSFAbsErr ( double pt )
{
    //If needed go to https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC
    //this pt range MUST match the binning of the error array provided by the BTV in the above link

    double ptarray[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
    double SFb_error[] = { 0.0484285,  0.0126178,  0.0120027,  0.0141137, 0.0145441, 0.0131145, 0.0168479, 0.0160836, 0.0126209, 0.0136017, 0.019182, 0.0198805, 0.0386531, 0.0392831, 0.0481008, 0.0474291 };

    int ptbin = -1;
    int ptarray_Size = sizeof(ptarray) / sizeof(ptarray[0]);
    int SFb_error_array_Size = sizeof(SFb_error) / sizeof(SFb_error[0]);
    
    if( ptarray_Size != SFb_error_array_Size+1 ) {
        cout<<"You wrote 2 arrays with sizes that don't agree together!! Fix this."<<endl;
        exit(100);
    }
    for ( int i=0; i<(ptarray_Size-1); i++){
        if ( pt > ptarray[i] && pt < ptarray[i+1]) {
            ptbin = i;
            break;
        }
    }

    if ( ptbin > ptarray_Size ) {
        return 2 * SFb_error[ptarray_Size];
    } else if ( ptbin < 0){
        return 2 * SFb_error[0];
    } else {
        return SFb_error[ptbin];
    };
}

double Analysis::CJetSFAbsErr ( double pt)
{
    return 2 * BJetSFAbsErr ( pt );
}

void Analysis::SetBTagFile(TString btagFile)
{
    this->btagFile = btagFile;
}

void Analysis::SetChannel(TString channel)
{
    this->channel = channel;
    this->channelPdgIdProduct = 
        channel == "ee" ? -11*11
        : channel == "emu" ? -11*13
        : -13*13;
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
}

void Analysis::SetMC(bool isMC)
{
    this->isMC = isMC;
}

void Analysis::SetTrueLevelDYChannel(int dy)
{
    trueDYchannelCut = dy;
    if (dy) {
        std::cout << "Include true-level filter for Z decay to pdgid " << dy << "\n";
        
        //create function to check the DY decay channel
        checkZDecayMode = [&, dy](Long64_t entry) -> bool {
            b_ZDecayMode->GetEntry(entry);
            bool pass = false;
            for (const auto decayMode : *ZDecayMode) {
                if ((dy == 15 && decayMode > 15110000) ||
                    (dy == 13 && decayMode == 1313) ||
                    (dy == 11 && decayMode == 1111))
                {
                    pass = true;
                    break;
                }
            }
            return pass;
        };
        
    } else {
        checkZDecayMode = nullptr;
    }
}

    


void Analysis::SetPDF(int pdf_no)
{
    this->pdf_no = pdf_no;
}


void Analysis::SetOutputfilename(TString outputfilename)
{
    this->outputfilename = outputfilename;
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
    leptons = 0;
    lepPdgId = 0;
    lepPfIso = 0;
    lepCombIso = 0;
    lepDxyVertex0 = 0;
    jets = 0;
    jetBTagTCHE = 0;
    jetBTagCSV = 0;
    jetBTagSSVHE = 0;
    jetPartonFlavour = 0;
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
//     genJets = 0;
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
    GenMet = 0;

    // Set branch addresses and branch pointers
    if ( !tree ) return;
    fChain = tree;
    fChain->SetMakeClass ( 0 );
    fChain->SetBranchAddress("leptons", &leptons, &b_lepton );
    fChain->SetBranchAddress("lepPdgId", &lepPdgId, &b_lepPdgId );
    fChain->SetBranchAddress("lepPfIso", &lepPfIso, &b_lepPfIso );
    fChain->SetBranchAddress("lepCombIso", &lepCombIso, &b_lepCombIso );
    fChain->SetBranchAddress("lepDxyVertex0", &lepDxyVertex0, &b_lepDxyVertex0);
    fChain->SetBranchAddress("jets", &jets, &b_jet );
    fChain->SetBranchAddress("jetBTagTCHE", &jetBTagTCHE, &b_jetBTagTCHE );
    fChain->SetBranchAddress("jetBTagCSV", &jetBTagCSV, &b_jetBTagCSV );
    fChain->SetBranchAddress("jetBTagSSVHE", &jetBTagSSVHE, &b_jetBTagSSVHE );
    fChain->SetBranchAddress("jetPartonFlavour", &jetPartonFlavour, &b_jetPartonFlavour );
    //fChain->SetBranchAddress("genJet", &genJets, &b_genJet );
    fChain->SetBranchAddress("met", &met, &b_met );
    fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber );
    fChain->SetBranchAddress("triggerBits", &triggerBits, &b_triggerBits );
    fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock );
    fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber );
    fChain->SetBranchAddress("weightGenerator", &weightGenerator, &b_weightGenerator );
    if (pdf_no) fChain->SetBranchAddress("pdfWeights", &weightPDF, &b_weightPDF);
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
    fChain->SetBranchAddress("TopDecayMode", &topDecayMode, &b_TopDecayMode );
    fChain->SetBranchAddress("ZDecayMode", &ZDecayMode, &b_ZDecayMode);
    
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
        fChain->SetBranchAddress("BHadJetIndex", &BHadJetIndex, &b_BHadJetIndex );
        fChain->SetBranchAddress("AntiBHadJetIndex", &AntiBHadJetIndex, &b_AntiBHadJetIndex );
        fChain->SetBranchAddress("BHadrons", &BHadrons, &b_BHadrons );
        fChain->SetBranchAddress("AntiBHadrons", &AntiBHadrons, &b_AntiBHadrons);
        fChain->SetBranchAddress("BHadronFromTop", &BHadronFromTopB, &b_BHadronFromTopB );
        fChain->SetBranchAddress("AntiBHadronFromTopB", &AntiBHadronFromTopB, &b_AntiBHadronFromTopB );
        fChain->SetBranchAddress("BHadronVsJet", &BHadronVsJet, &b_BHadronVsJet );
        fChain->SetBranchAddress("AntiBHadronVsJet", &AntiBHadronVsJet, &b_AntiBHadronVsJet );

        fChain->SetBranchAddress("GenMET", &GenMet, &b_GenMet);
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
    b_lepton->GetEntry(entry); //!
    b_lepPdgId->GetEntry(entry); //!
    b_jet->GetEntry(entry); //!
    b_met->GetEntry(entry); //!
    
    b_eventNumber->GetEntry(entry); //!
    b_runNumber->GetEntry(entry); //!
    b_lumiBlock->GetEntry(entry); //!
    
    //special variables, not used currently
//     b_lepPfIso->GetEntry(entry); //!
//     b_lepCombIso->GetEntry(entry); //!
//     b_lepDxyVertex0->GetEntry(entry);
    
    
    //We only need CSV
//     b_jetBTagTCHE->GetEntry(entry); //!
    b_jetBTagCSV->GetEntry(entry); //!
//     b_jetBTagSSVHE->GetEntry(entry); //!
    
    b_jetPartonFlavour->GetEntry(entry); //!
    b_triggerBits->GetEntry(entry); //!
    b_weightGenerator->GetEntry(entry); //!
    b_vertMulti->GetEntry(entry); //!
    b_vertMultiTrue->GetEntry(entry); //!

    //b_genJet->GetEntry(entry); //!
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

}

void Analysis::GetSignalBranches ( Long64_t & entry )
{
    b_GenTop->GetEntry(entry); //!
    b_GenAntiTop->GetEntry(entry); //!
    b_GenLepton->GetEntry(entry); //!
    b_GenAntiLepton->GetEntry(entry); //!
    b_GenB->GetEntry(entry); //!
    b_GenAntiB->GetEntry(entry); //!
    b_GenNeutrino->GetEntry(entry);   //!
    b_GenAntiNeutrino->GetEntry(entry);   //!

    b_GenMet->GetEntry(entry);
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
    //find opposite-charge leading two leptons
    //the first lepton is always at index 0 because we only have two different charges
    
    // ---> this would be true if the electrons and muons had all cuts applied
    auto leptonPassesCut = [](const LV &lep){ return abs(lep.eta()) < 2.4 && lep.pt() > 20;};
    
    //find first lepton:
    for (LeadLeptonNumber = 0; LeadLeptonNumber < leptons->size(); ++LeadLeptonNumber) {
        if (leptonPassesCut(leptons->at(LeadLeptonNumber))) break;
    }
    
    //find second lepton
    for (size_t i = LeadLeptonNumber + 1; i < leptons->size(); ++i) {
        if (!leptonPassesCut(leptons->at(i))) continue;
        int product = lepPdgId->at(LeadLeptonNumber) * lepPdgId->at(i);
        if (product < 0) {
            NLeadLeptonNumber = i;
            return product == channelPdgIdProduct;
        }            
    }
    
    return false;
}

double Analysis::get2DSF(TH2* histo, double x, double y)
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
    for ( size_t i = 0; i < jets->size(); ++i ) {
        double pt = jets->at(i).Pt();
        double eta = abs(jets->at(i).Eta());
        if ( pt > 30 && eta < 2.4 ) {
            int partonFlavour = abs(jetPartonFlavour->at(i)); //store absolute value
            if (partonFlavour == 0) continue;
            int ptbin, etabin, dummy;
            bEff->GetBinXYZ(bEff->FindBin(pt, eta), ptbin, etabin, dummy);
            //overflow to last bin
            ptbin = std::min(ptbin, bEff->GetNbinsX());
            etabin = std::min(etabin, bEff->GetNbinsY());
            //do the type-jet selection & Eff and SF obtention
            double SF_Error=0;
            if ( partonFlavour == 5 ) { //b-quark
                eff=bEff->GetBinContent ( ptbin, etabin );
                if(systematic.Contains("BEFF_UP")){eff = eff + bEff->GetBinError ( ptbin, etabin );}
                else if(systematic.Contains("BEFF_DOWN")){eff = eff - bEff->GetBinError ( ptbin, etabin );} 
                SFPerJet = BJetSF( pt, eta );
                SF_Error = BJetSFAbsErr ( pt );
            } else if ( partonFlavour == 4 ) { //c-quark
                eff=cEff->GetBinContent ( ptbin, etabin );
                if(systematic.Contains("CEFF_UP")){eff = eff + cEff->GetBinError ( ptbin, etabin );}
                else if(systematic.Contains("CEFF_DOWN")){eff = eff - cEff->GetBinError ( ptbin, etabin );} 
                SFPerJet = CJetSF( pt, eta );
                SF_Error = CJetSFAbsErr ( pt );
            } else if ( partonFlavour != 0 ) { //l-quark
                eff=lEff->GetBinContent ( ptbin, etabin );
                if(systematic.Contains("LEFF_UP")){eff = eff + lEff->GetBinError ( ptbin, etabin );}
                else if(systematic.Contains("LEFF_DOWN")){eff = eff - lEff->GetBinError ( ptbin, etabin );} 
                SFPerJet = LJetSF( pt, eta, "central");
                if ( systematic.BeginsWith("BTAG_LJET_UP") ) {   //systematic variation of l-jets for inclusive XSection measurement
                    SFPerJet = LJetSF ( pt, eta, "up");
                } else if ( systematic.BeginsWith("BTAG_LJET_DOWN") ) { //systematic variation of l-jets for inclusive XSection measurement
                    SFPerJet = LJetSF ( pt, eta, "down");
                } else if ( systematic.BeginsWith("BTAG_LJET_") ){ //systematic variations for differential XSection
                    if ( (systematic.Contains("PT_UP") && pt>btag_ptmedian) || (systematic.Contains("PT_DOWN") && pt<btag_ptmedian)
                        || (systematic.Contains("ETA_UP") && eta>btag_etamedian) || (systematic.Contains("ETA_DOWN") && eta<btag_etamedian) ){
                        SFPerJet = LJetSF ( pt, eta, "up");
                    }
                    else {
                        SFPerJet = LJetSF ( pt, eta, "down");
                    }
                }
            } else {
                cout<<"I found a jet in event "<<eventNumber<<" which is not b, c nor light: "<<partonFlavour<<endl;
                return kFALSE;
            }
            if ( eff <= 0 ) eff = 1;
            //calculate both numerator and denominator for per-event SF calculation
            //consider also the UP and DOWN variation for systematics calculation. Same procedure as PU
            OneMinusEff = OneMinusEff* ( 1-eff );
            double sf = SFPerJet;
            if ( systematic.Contains("LJET") ) { SF_Error = 0; } //For l-jet systematics set to 0 the b- and c-jet errors
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

    if( abs(1.-OneMinusEff) < 1e-8 || abs(1.-OneMinusSEff) < 1e-8 ) return 1;

    //per-event SF calculation (also the UP and DOWN variations)
    double scale_factor = ( 1.-OneMinusSEff ) / ( 1.-OneMinusEff );
    
    if ( abs(scale_factor - 1.)>0.05 ){scale_factor = 1;}
    
    return scale_factor;
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
        cout << "Currently triggerEfficieny files can be found in the HEAD version of diLeptonic folder\n\n";
        return;
    }
    
    //Right now pT efficiency flat ==> Not used
    h_TrigSFeta = dynamic_cast<TH2*>(trigEfficiencies.Get("scalefactor eta2d with syst"));
    if ( !h_TrigSFeta ) {
        cout<<"TH2 >>TH scalefactor eta<< is not in the file "<<trigEfficiencies.GetName()<<"\n";
        return;
    }
    
    if (systematic.BeginsWith("TRIG_")) {
        double factor = systematic.EndsWith("_UP") ? 1 : -1;
        for (int i = 1; i <= h_TrigSFeta->GetNbinsX(); ++i) {
            for (int j = 1; j <= h_TrigSFeta->GetNbinsY(); ++j) {
                h_TrigSFeta->SetBinContent(i, j,
                    h_TrigSFeta->GetBinContent(i,j) + factor*h_TrigSFeta->GetBinError(i,j));
            }
        }
    }
    
    h_TrigSFeta->SetDirectory(0);
    trigEfficiencies.Close();
}

double Analysis::getTriggerSF(const LV& lep1, const LV& lep2) {
    
    //For 'ee' and 'mumu' channels Xaxis of the 2D plots is the highest pT lepton
    // for the 'emu' channel Xaxis is the electron and Y axis muon
    
    if (!h_TrigSFeta) return 1;
    return get2DSF(h_TrigSFeta, abs(lep1.eta()), abs(lep2.eta()));
}

double Analysis::getLeptonIDSF(const LV& lep1, const LV& lep2, int lep1pdgId, int lep2pdgId) {
    if (!h_MuonIDSFpteta || !h_ElectronIDSFpteta) return 1;
    if (abs(lep1pdgId)==11 && abs(lep2pdgId)==11) return get2DSF(h_ElectronIDSFpteta, lep1.Eta(), lep1.pt()) * get2DSF(h_ElectronIDSFpteta, lep2.Eta(), lep2.pt());
    else if (abs(lep1pdgId)==13 && abs(lep2pdgId)==13) return get2DSF(h_MuonIDSFpteta, lep1.Eta(), lep1.pt()) * get2DSF(h_MuonIDSFpteta, lep2.Eta(), lep2.pt());
    else if (abs(lep1pdgId)==13 && abs(lep2pdgId)== 11) return get2DSF(h_MuonIDSFpteta, lep1.Eta(), lep1.pt()) * get2DSF(h_ElectronIDSFpteta, lep2.Eta(), lep2.pt());
    else return get2DSF(h_ElectronIDSFpteta, lep1.Eta(), lep1.pt())
        * get2DSF(h_MuonIDSFpteta, lep2.Eta(), lep2.pt());
}

void Analysis::prepareLeptonIDSF() {
    h_MuonIDSFpteta = nullptr; h_ElectronIDSFpteta = nullptr;
    
    TFile MuonEfficiencies(TString("MuonSFtop12028.root"));
    if (MuonEfficiencies.IsZombie()) {
        cout << "Muon Id/Iso efficiencies not found. Assuming ScaleFactor = 1.\n";
///        cout << "Currently triggerEfficieny files can be found in the HEAD version of diLeptonic folder\n\n";
        return;
    }

    TFile ElecEfficiencies(TString("ElectronSFtop12028.root"));
    if (ElecEfficiencies.IsZombie()) {
        cout << "Electron Id/Iso efficiencies not found. Assuming ScaleFactor = 1.\n";
///        cout << "Currently triggerEfficieny files can be found in the HEAD version of diLeptonic folder\n\n";
        return;
    }
    
    //Right now pT efficiency flat ==> Not used
    h_MuonIDSFpteta = dynamic_cast<TH2*>(MuonEfficiencies.Get("MuonSFtop12028"));
    if ( !h_MuonIDSFpteta ) {
        cout<<"TH2 >>TH scalefactor << is not in the file "<<MuonEfficiencies.GetName()<<"\n";
        return;
    }

    h_ElectronIDSFpteta = dynamic_cast<TH2*>(ElecEfficiencies.Get("GlobalSF"));
    if ( !h_ElectronIDSFpteta ) {
        cout<<"TH2 >>TH scalefactor << is not in the file "<<ElecEfficiencies.GetName()<<"\n";
        return;
    }

    if (systematic.BeginsWith("LEPT_")) {
        double factor = systematic.EndsWith("_UP") ? 1 : -1;
        for (int i = 1; i <= h_MuonIDSFpteta->GetNbinsX(); ++i) {
            for (int j = 1; j <= h_MuonIDSFpteta->GetNbinsY(); ++j) {
                h_MuonIDSFpteta->SetBinContent(i, j,
                    h_MuonIDSFpteta->GetBinContent(i,j) + factor*h_MuonIDSFpteta->GetBinError(i,j));
            }
        }

        for (int i = 1; i <= h_ElectronIDSFpteta->GetNbinsX(); ++i) {
            for (int j = 1; j <= h_ElectronIDSFpteta->GetNbinsY(); ++j) {
                h_ElectronIDSFpteta->SetBinContent(i, j,h_ElectronIDSFpteta->GetBinContent(i,j) + factor*h_ElectronIDSFpteta->GetBinError(i,j));
            }
        }
    }

    h_MuonIDSFpteta->SetDirectory(0);
    MuonEfficiencies.Close();
    h_ElectronIDSFpteta->SetDirectory(0);
    ElecEfficiencies.Close();


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
             << "   ./load_Analysis -f ttbarsignalplustau.root -c emu\n"
             << "   ./load_Analysis -f ttbarsignalplustau.root -c ee\n"
             << "   ./load_Analysis -f ttbarsignalplustau.root -c mumu\n"
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

void Analysis::FillBinnedControlPlot(TH1* h_differential, double binvalue, 
                                     TH1* h_control, double value, double weight)
{
    auto pair = (*binnedControlPlots)[h_differential->GetName()];
    auto bin = pair.first->FindBin(binvalue);
    auto m = pair.second.at(bin);
    TH1* h = m[h_control->GetName()];
    if (!h) { 
        std::cerr << "Error: please CreateBinnedControlPlots for " << h_differential->GetName() 
            << " and " << h_control->GetName() << std::endl;
        exit(911);
    }
    h->Fill(value, weight);
}

void Analysis::CreateBinnedControlPlots(TH1* h_differential, TH1* h_control)
{
    HistoListReader histoList("HistoList");
    if (histoList.IsZombie()) { cout << "Need a HistoList to create binned control plots!\n"; exit(273); }
    auto &pair = (*binnedControlPlots)[h_differential->GetName()];
    pair.first = histoList.getPlotProperties(h_differential->GetName()).getClonedHistogram();
    std::string name = "bcp_";
    name.append(h_differential->GetName()).append("_bin_");
    //create maps if we are called for the first time with a certain h_differential
    if (pair.second.size() == 0) {
        for (int i = 0; i <= pair.first->GetNbinsX() + 1; ++i)
            pair.second.push_back(std::map<std::string, TH1*>());
    }
    //now really create the histograms
    for (int i = 0; i <= pair.first->GetNbinsX() + 1; ++i) {
        std::string binning = 
            i == 0 ? "underflow" :
            i == pair.first->GetNbinsX() + 1 ? "overflow" :
            d2s(pair.first->GetBinLowEdge(i)) + " to " + d2s(pair.first->GetBinLowEdge(i+1));
        binning = std::string(" (") + h_differential->GetName() + " " + binning + ")";
        std::string n = name + std::to_string(i) + "_" + h_control->GetName();
        pair.second[i][h_control->GetName()] = store(
            new TH1D(n.c_str(), (std::string(h_control->GetName())+binning).c_str(), 
                     h_control->GetNbinsX(), h_control->GetBinLowEdge(1), 
                     h_control->GetBinLowEdge(h_control->GetNbinsX()+1)));
    }
}

bool Analysis::calculateKinReco(const LV& leptonMinus, const LV& leptonPlus, double JETPTCUT)
{
    auto sols = GetKinSolutions(leptonMinus, leptonPlus, jets, jetBTagCSV, met, JETPTCUT);
    if (sols.size()) {
        const TtDilepEvtSolution &sol = sols[0];
        HypTop->clear();        HypTop->push_back(TLVtoLV(sol.top));
        HypAntiTop->clear();    HypAntiTop->push_back(TLVtoLV(sol.topBar));
        HypLepton->clear();     HypLepton->push_back(TLVtoLV(sol.lm));
        HypAntiLepton->clear(); HypAntiLepton->push_back(TLVtoLV(sol.lp));
        HypBJet->clear();       HypBJet->push_back(TLVtoLV(sol.jetB));
        HypAntiBJet->clear();   HypAntiBJet->push_back(TLVtoLV(sol.jetBbar));
        HypNeutrino->clear();   HypNeutrino->push_back(TLVtoLV(sol.neutrino));
        HypAntiNeutrino->clear(); HypAntiNeutrino->push_back(TLVtoLV(sol.neutrinoBar));
        HypJet0index->clear();  HypJet0index->push_back(sol.jetB_index);
        HypJet1index->clear();  HypJet1index->push_back(sol.jetBbar_index);
    }
    //check for strange events
    if (false && HypTop->size()) {
        double Ecm = (HypTop->at(0) + HypAntiTop->at(0) 
                        + HypLepton->at(0) + HypAntiLepton->at(0)
                        + HypNeutrino->at(0) + HypAntiNeutrino->at(0)).E();
        //does event's CM energy exceed LHC's energy?
        if (Ecm > 8000) {
            static int seCounter = 0;
            cout << "\n" << ++seCounter << " - Strange event: " << runNumber<<":"<<lumiBlock<<":"<<eventNumber
            << "\ntop:  " << HypTop->at(0) << " tbar: " << HypAntiTop->at(0)
            << "\nNeutrino:  " << HypNeutrino->at(0) << " NeutrinoBar: " << HypAntiNeutrino->at(0)
            <<"\n";
        }
    }
    return sols.size() > 0;
}
