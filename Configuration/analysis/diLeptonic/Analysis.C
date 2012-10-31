#define Analysis_cxx

#include "Analysis.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <fstream>
#include <iostream>
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

    return -1;
}

void Analysis::Begin ( TTree * )
{
    EventCounter = 0;
    bEff = 0;

    //By now defined the per-jet SFs vary according to:
    //   BTag_Up   ==> pt>ptmedian vary DOWN, pt<ptmedian vary UP
    //   BTag_Down ==> pt>ptmedian vary UP, pt<ptmedian vary DOWN

    //load per-jet efficienciies file and Histograms
    TFile *bEfficiencies;
    if ( btagFile!="" ) {
        bEfficiencies = TFile::Open( btagFile );
    } else {
        cout<<"WARNING!!! Provide b tag efficiencies before running"<<endl;
        return;
    }

    if ( bEfficiencies ) {
        cout<<"File "<< btagFile << " does not exist. Running without btagsf!!!"<<endl;
        return;
    }
    bEff = dynamic_cast<TH2*>(bEfficiencies->Get( "BEffPerJet" ));
    if (!bEff) {
        cout<<"Histogram bEff is not in the file "<<bEfficiencies->GetName();
        return;
    }
    cEff = dynamic_cast<TH2*>(bEfficiencies->Get ( "CEffPerJet" ));
    if (!cEff) {
        cout<<"Histogram cEff is not in the file "<<bEfficiencies->GetName();
        return;
    }
    lEff = dynamic_cast<TH2*>(bEfficiencies->Get ( "LEffPerJet" ));
    if (!lEff) {
        cout<<"Histogram lEff is not in the file "<<bEfficiencies->GetName();
        return;
    }

    //load the histograms in memory, to avoid memory leaks
    bEff->SetDirectory ( 0 );
    cEff->SetDirectory ( 0 );
    lEff->SetDirectory ( 0 );
    bEfficiencies->Close();
    bEfficiencies->Delete();
    // END: BTag SF calculation neccessary stuff
}

void Analysis::SlaveBegin ( TTree * )
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    h_step5 = new TH1D ( "step5", "event count at step 5", 10, 0, 10 );
    h_step6 = new TH1D ( "step6", "event count at step 6", 10, 0, 10 );
    h_step7 = new TH1D ( "step7", "event count at step 7", 10, 0, 10 );
    h_step8 = new TH1D ( "step8", "event count at step 8", 10, 0, 10 );
    h_step9 = new TH1D ( "step9", "event count at step 9", 10, 0, 10 );

    //h_jetMultiAll = new TH1D ( "HypjetMultiAll", "Jet Multiplicity (AllJets)", 10, -0.5, 9.5 );
    h_jetMultiXSec = new TH1D ( "HypjetMultiXSec", "Jet Multiplicity (for cross-section)", 10, -0.5, 9.5 );
    h_jetMulti = new TH1D ( "HypjetMulti", "Jet Multiplicity", 10, -0.5, 9.5 );
    h_jetMulti_diLep = new TH1D ( "HypjetMulti_diLep", "Jet Multiplicity (after dilepton)", 10, -0.5, 9.5 );
    h_jetMultiNoPU = new TH1D ( "HypjetMultiNoPU", "Jet Multiplicity (No Pileup or lumi weight)", 10, -0.5, 9.5 );
//     h_jetMultiVisTop = new TH1D ( "HypjetMultiVisTop", "Jet Multiplicity for Visible Top (No Pileup or lumi Weight)", 10, -0.5, 9.5 );
    h_BjetMulti = new TH1D ( "HypBjetMulti", "B-Jet Multiplicity", 10, -0.5, 9.5 );

    h_HypTTBarRapidity = new TH1D ( "HypTTBarRapidity", "Rapidity of TTbar System (HYP)", 100, -5, 5 );
    h_HypTTBarpT = new TH1D ( "HypTTBarpT", "pT of TTbar System (HYP)", 500, 0, 500 );
    h_HypTTBarMass = new TH1D ( "HypTTBarMass", "Mass of TTbar System (HYP)", 2000, 0, 2000 );
    h_HypLLBarMass = new TH1D ( "HypLLBarMass", "Mass of LLbar System (HYP)", 500, 0, 1000 );
    h_HypLLBarpT = new TH1D ( "HypLLBarpT", "pT of LLbar System (HYP)", 200, 0, 1000 );

    h_GenTTBarMass = new TH1D ( "GenTTBarMass", "Mass of TTbar System(GEN)", 1200, 0, 1200 );
    h_GenTTBarRapidity = new TH1D ( "GenTTBarRapidity", "Rapidity of TTbar System(GEN)", 100, -5, 5 );
    h_GenTTBarpT = new TH1D ( "GenTTBarpT", "pT of TTbar System(GEN)", 1200, 0, 1200 );
    h_GenLLBarpT = new TH1D ( "GenLLBarpT", "pT of LLbar System(GEN)", 200, 0, 1000 );
    h_GenLLBarMass = new TH1D ( "GenLLBarMass", "Mass of LLbar System(GEN)", 500, 0, 1000 );

    h_VisGenTTBarMass = new TH1D ( "VisGenTTBarMass", "Mass of TTbar System(VisGEN)", 1200, 0, 1200 );
    h_VisGenTTBarRapidity = new TH1D ( "VisGenTTBarRapidity", "Rapidity of TTbar System(VisGEN)", 100, -5, 5 );
    h_VisGenTTBarpT = new TH1D ( "VisGenTTBarpT", "pT of TTbar System(VisGEN)", 1200, 0, 1200 );
    h_VisGenTopRapidity = new TH1D ( "VisGenTopRapidity", "Rapidity of Top(VisGEN)", 100, -5, 5 );
    h_VisGenAntiTopRapidity = new TH1D ( "VisGenAntiTopRapidity", "Rapidity of AntiTop(VisGEN)", 100, -5, 5 );

    h_VisGenLLBarpT = new TH1D ( "VisGenLLBarpT", "pT of LLbar System(VisGEN)", 200, 0, 1000 );
    h_VisGenLLBarMass = new TH1D ( "VisGenLLBarMass", "Mass of LLbar System(VisGEN)", 500, 0, 1000 );

    h_RecoTTBarMass = new TH1D ( "RecoTTBarMass","Mass of TTbar System (HYP)",1200,0,1200 );
    h_RecoTTBarRapidity = new TH1D ( "RecoTTBarRapidity","Rapidity of TTbar System (HYP)",100,-5,5 );
    h_RecoTTBarpT = new TH1D ( "RecoTTBarpT","pT of TTbar System (HYP)",1200,0,1200 );
    h_RecoToppT = new TH1D ( "RecoToppT","pT of Top (HYP)",1200,0,1200 );
    h_RecoAntiToppT = new TH1D ( "RecoAntiToppT","pT of AntiTop (HYP)",1200,0,1200 );
    h_RecoTopRapidity = new TH1D ( "RecoTopRapidity","Rapidity of Top (HYP)",100,-5,5 );
    h_RecoAntiTopRapidity = new TH1D ( "RecoAntiTopRapidity","Rapidity of AntiTop (HYP)",100,-5,5 );

    h_RecoBJetpT = new TH1D ( "RecoBJetpT","pT of BJet (HYP)",80,0,400 );
    h_RecoAntiBJetpT = new TH1D ( "RecoAntiBJetpT","pT of AntiBJet (HYP)",80,0,400 );
    h_RecoBJetRapidity = new TH1D ( "RecoBJetRapidity","Rapidity of BJet (HYP)",100,-5,5 );
    h_RecoAntiBJetRapidity = new TH1D ( "RecoAntiBJetRapidity","Rapidity of AntiBJet (HYP)",100,-5,5 );
    h_RecoBJetEta = new TH1D ( "RecoBJetEta","#eta of BJet (HYP)",100,-5,5 );
    h_RecoAntiBJetEta = new TH1D ( "RecoAntiBJetEta","#eta of AntiBJet (HYP)",100,-5,5 );

    h_RecoLLBarMass = new TH1D ( "RecoLLBarMass","Mass of LLbar System (HYP)",500,0,1000 );
    h_RecoLLBarpT = new TH1D ( "RecoLLBarpT","pT of LLbar System (HYP)",200,0,1000 );
    h_RecoLeptonpT = new TH1D ( "RecoLeptonpT","pT of Lepton (HYP)",240,0,1200 );
    h_RecoAntiLeptonpT = new TH1D ( "RecoAntiLeptonpT","pT of AntiLepton (HYP)",240,0,1200 );
    h_RecoLeptonEta = new TH1D ( "RecoLeptonEta","Eta of Lepton (HYP)",100,-5,5 );
    h_RecoAntiLeptonEta = new TH1D ( "RecoAntiLeptonEta","Eta of AntiLepton (HYP)",100,-5,5 );

    h_RecoLLBarDPhi = new TH1D ( "RecoLLBarDPhi", "#Delta #Phi (Lep, AntiLep) (Reco)", 112, -0.1, 3.25 );
    h_RecoLeptonantiBjetMass = new TH1D ( "RecoLeptonantiBjetMass", "M(Lep, AntiBJet) (Reco)", 500, 0, 1000 );
    h_RecoAntiLeptonBjetMass = new TH1D ( "RecoAntiLeptonBjetMass", "M(AntiLep, BJet) (Reco)", 500, 0, 1000 );
    h_RecoJetMult = new TH1D ( "RecoJetMult", "Jet Multiplicty (Reco)", 26, -0.5, 25.5 );
    
    h_VisGenAll = new TH1D ( "VisGenAll", "All Visible Generated particles (IM)", 40, 0, 400 );
    h_GenAll = new TH1D ( "GenAll", "AllGenerated particles (IM)", 40, 0, 400 );
    Allh1 = new TH1D ( "Allh1", "DiLepton Mass", 40, 0, 400 );
    h_diLepMassFull = new TH1D ( "DIMFull", "DiLepton Mass (Full Range)", 100, 0, 300 );
    h_diLepMassFull_fullSel = new TH1D ( "DIMFull_fullSel", "DiLepton Mass (Full Range)", 100, 0, 300 );
    Looseh1 = new TH1D ( "Looseh1", "DiLepton Mass", 40, 0, 400 );
    Zh1 = new TH1D ( "Zh1", "DiLepton Mass in Z Window", 40, 0, 400 );
    TTh1 = new TH1D ( "TTh1", "DiLepton Mass out of Z Window", 40, 0, 400 );

    h_vertMulti = new TH1D ( "vertMulti", "Primary Vertex Multiplicity", 30, 0, 30 );
    h_vertMulti_noPU = new TH1D ( "vertMulti_noPU", "Primary Vertex Multiplicity (no Pileup)", 30, 0, 30 );
    h_MET = new TH1D ( "MET", "Missing Transverse Energy", 80, 0, 400 );
    h_jetpT = new TH1D ( "jetpT", "jet pT", 80, 0, 400 );
    h_jetHT = new TH1D ( "jetHT", "jet HT", 80, 0, 1000 );

    h_MuonpT = new TH1D ( "MuonpT", "Muon pT (emu channel)", 80, 0, 400 );
    h_MuonEta = new TH1D ( "MuonEta", "Muon Eta (emu channel)", 100, -5, 5 );
    h_ElectronpT = new TH1D ( "ElectronpT", "Electron pT (emu channel)", 80, 0, 400 );
    h_ElectronEta = new TH1D ( "ElectronEta", "Electron Eta (emu channel)", 100, -5, 5 );

    h_LeptonpT = new TH1D ( "LeptonpT", "Lepton pT", 80, 0, 400 );
    h_LeptonEta = new TH1D ( "LeptonEta", "Lepton Eta", 100, -5, 5 );
    h_LeptonpT_diLep = new TH1D ( "LeptonpT_diLep", "Lepton pT (after dilepton cut)", 80, 0, 400 );
    h_LeptonEta_diLep = new TH1D ( "LeptonEta_diLep", "Lepton Eta (after dilepton cut)", 100, -5, 5 );

    h_AntiLeptonpT = new TH1D ( "AntiLeptonpT", "AntiLepton pT", 80, 0, 400 );
    h_AntiLeptonEta = new TH1D ( "AntiLeptonEta", "AntiLepton Eta", 100, -5, 5 );
    h_AntiLeptonpT_diLep = new TH1D ( "AntiLeptonpT_diLep", "Lepton pT (after dilepton cut)", 80, 0, 400 );
    h_AntiLeptonEta_diLep = new TH1D ( "AntiLeptonEta_diLep", "Lepton Eta (after dilepton cut)", 100, -5, 5 );

    h_HypToppT = new TH1D ( "HypToppT", "Top pT", 400, 0, 400 );
    h_HypTopEta = new TH1D ( "HypTopEta", "Top pT", 100, -5, 5 );
    h_HypTopMass = new TH1D ( "HypTopMass", "Top Mass", 80, 0, 400 );
    h_HypTopRapidity = new TH1D ( "HypTopRapidity", "Top Rapidity", 100, -5, 5 );

    h_HypAntiToppT = new TH1D ( "HypAntiToppT", "AntiTop pT", 400, 0, 400 );
    h_HypAntiTopEta = new TH1D ( "HypAntiTopEta", "AntiTop pT", 100, -5, 5 );
    h_HypAntiTopMass = new TH1D ( "HypAntiTopMass", "AntiTop Mass", 80, 0, 400 );
    h_HypAntiTopRapidity = new TH1D ( "HypAntiTopRapidity", "Top Rapidity", 100, -5, 5 );

    h_HypLeptonpT = new TH1D ( "HypLeptonpT", "Lepton Hypothesis pT", 80, 0, 400 );
    h_HypLeptonEta = new TH1D ( "HypLeptonEta", "Lepton Eta", 100, -5, 5 );

    h_HypAntiLeptonpT = new TH1D ( "HypAntiLeptonpT", "AntiLepton Hypothesis pT", 80, 0, 400 );
    h_HypAntiLeptonEta = new TH1D ( "HypAntiLeptonEta", "AntiLepton Hypothesis Eta", 100, -5, 5 );

    h_HypBJetpT = new TH1D ( "HypBJetpT", "B Hypothesis pT", 80, 0, 400 );
    h_HypBJetEta = new TH1D ( "HypBJetEta", "B Hypothesis Eta", 100, -5, 5 );
    h_HypBJetRapidity = new TH1D ( "HypBJetRapidity", "B Hypothesis Eta", 100, -5, 5 );

    h_HypAntiBJetpT = new TH1D ( "HypAntiBJetpT", "AntiB Hypothesis pT", 80, 0, 400 );
    h_HypAntiBJetEta = new TH1D ( "HypAntiBJetEta", "AntiB Hypothesis Eta", 100, -5, 5 );
    h_HypAntiBJetRapidity = new TH1D ( "HypAntiBJetRapidity", "AntiB Hypothesis Eta", 100, -5, 5 );

    h_VisGenToppT = new TH1D ( "VisGenToppT", "Top pT (VisGen)", 400, 0, 400 );
    h_VisGenTopEta = new TH1D ( "VisGenTopEta", "Top Eta (VisGen)", 100, -5, 5 );

    h_VisGenAntiToppT = new TH1D ( "VisGenAntiToppT", "AntiTop pT (VisGen)", 400, 0, 400 );
    h_VisGenAntiTopEta = new TH1D ( "VisGenAntiTopEta", "AntiTop pT (VisGen)", 100, -5, 5 );

    h_VisGenLeptonpT = new TH1D ( "VisGenLeptonpT", "Lepton VisGenothesis pT", 80, 0, 400 );
    h_VisGenLeptonEta = new TH1D ( "VisGenLeptonEta", "Lepton Eta", 100, -5, 5 );

    h_VisGenAntiLeptonpT = new TH1D ( "VisGenAntiLeptonpT", "AntiLepton VisGenothesis pT", 80, 0, 400 );
    h_VisGenAntiLeptonEta = new TH1D ( "VisGenAntiLeptonEta", "AntiLepton VisGenothesis Eta", 100, -5, 5 );

    h_VisGenBJetpT = new TH1D ( "VisGenBJetpT", "B VisGenothesis pT", 80, 0, 400 );
    h_VisGenBJetEta = new TH1D ( "VisGenBJetEta", "B VisGenothesis Eta", 100, -5, 5 );
    h_VisGenBJetRapidity = new TH1D ( "VisGenBJetRapidity", "B VisGenothesis Rapidity", 100, -5, 5 );

    h_VisGenAntiBJetpT = new TH1D ( "VisGenAntiBJetpT", "AntiB VisGenothesis pT", 80, 0, 400 );
    h_VisGenAntiBJetEta = new TH1D ( "VisGenAntiBJetEta", "AntiB VisGenothesis Eta", 100, -5, 5 );
    h_VisGenAntiBJetRapidity = new TH1D ( "VisGenAntiBJetRapidity", "AntiB VisGenothesis Rapidity", 100, -5, 5 );

    /*  h_VisGenBQuarkpT = new TH1D("VisGenBQuarkpT", "B Quark VisGenothesis pT", 80, 0, 400);
    h_VisGenBQuarkEta = new TH1D("VisGenBQuarkEta", "B Quark VisGenothesis Eta", 100, -5, 5);
    h_VisGenBQuarkRapidity = new TH1D("VisGenBQuarkRapidity", "B Quark VisGenothesis Rapidity", 100, -5, 5);

    h_VisGenAntiBQuarkpT = new TH1D("VisGenAntiBQuarkpT", "AntiB Quark VisGenothesis pT", 80, 0, 400);
    h_VisGenAntiBQuarkEta = new TH1D("VisGenAntiBQuarkEta", "AntiB Quark VisGenothesis Eta", 100, -5, 5);
    h_VisGenAntiBQuarkRapidity = new TH1D("VisGenAntiBQuarkRapidity", "AntiB Quark VisGenothesis Rapidity", 100, -5, 5);
    */
    /*h_GenToppT = new TH1D("GenToppT", "Top pT (Gen)", 80, 0, 400);
    h_GenTopEta = new TH1D("GenTopEta", "Top Eta (Gen)", 100, -5, 5);
    h_GenTopRapidity = new TH1D("GenTopRapidity", "Top Rapidity (Gen)", 100, -5, 5);

    h_GenAntiToppT = new TH1D("GenAntiToppT", "AntiTop pT (Gen)", 80, 0, 400);
    h_GenAntiTopEta = new TH1D("GenAntiTopEta", "AntiTop Eta (Gen)", 100, -5, 5);
    h_GenAntiTopRapidity = new TH1D("GenAntiTopRapidity", "AntiTop Rapidity (Gen)", 100, -5, 5);

    h_GenLeptonpT = new TH1D("GenLeptonpT", "Lepton Genothesis pT", 80, 0, 400);
    h_GenLeptonEta = new TH1D("GenLeptonEta", "Lepton Eta", 100, -5, 5);

    h_GenAntiLeptonpT = new TH1D("GenAntiLeptonpT", "AntiLepton Genothesis pT", 80, 0, 400);
    h_GenAntiLeptonEta = new TH1D("GenAntiLeptonEta", "AntiLepton Genothesis Eta", 100, -5, 5);

    h_GenBQuarkpT = new TH1D("GenBQuarkpT", "B Quark Genothesis pT", 80, 0, 400);
    h_GenBQuarkEta = new TH1D("GenBQuarkEta", "B Quark Genothesis Eta", 100, -5, 5);
    h_GenBQuarkRapidity = new TH1D("GenBQuarkRapidity", "B Quark Genothesis Rapidity", 100, -5, 5);

    h_GenAntiBQuarkpT = new TH1D("GenAntiBQuarkpT", "AntiB Quark Genothesis pT", 80, 0, 400);
    h_GenAntiBQuarkEta = new TH1D("GenAntiBQuarkEta", "AntiB Quark Genothesis Eta", 100, -5, 5);
    h_GenAntiBQuarkRapidity = new TH1D("GenAntiBQuarkRapidity", "AntiB Quark Genothesis Rapidity", 100, -5, 5);

    h_GenBJetpT = new TH1D("GenBJetpT", "B Genothesis pT", 80, 0, 400);
    h_GenBJetEta = new TH1D("GenBJetEta", "B Genothesis Eta", 100, -5, 5);
    h_GenBJetRapidity = new TH1D("GenBJetRapidity", "B Genothesis Rapidity", 100, -5, 5);

    h_GenAntiBJetpT = new TH1D("GenAntiBJetpT", "AntiB Genothesis pT", 80, 0, 400);
    h_GenAntiBJetEta = new TH1D("GenAntiBJetEta", "AntiB Genothesis Eta", 100, -5, 5);
    h_GenAntiBJetRapidity = new TH1D("GenAntiBJetRapidity", "Anti B Genothesis Rapidity", 100, -5, 5);
    */
    h_GenRecoBJetpT = new TH2D ( "GenRecoBJetpT", "Gen/Reco Matching", 80, 0, 400, 80, 0, 400 );
    h_GenRecoBJetEta = new TH2D ( "GenRecoBJetEta", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 );
    h_GenRecoBJetRapidity = new TH2D ( "GenRecoBJetRapidity", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 );

    h_GenRecoAntiBJetpT = new TH2D ( "GenRecoAntiBJetpT", "Gen/Reco Matching", 80, 0, 400, 80, 0, 400 );
    h_GenRecoAntiBJetEta = new TH2D ( "GenRecoAntiBJetEta", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 );
    h_GenRecoAntiBJetRapidity = new TH2D ( "GenRecoAntiBJetRapidity", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 );

    h_GenRecoLeptonEta = new TH2D ( "GenRecoLeptonEta", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 );
    h_GenRecoAntiLeptonEta = new TH2D ( "GenRecoAntiLeptonEta", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 );
    h_GenRecoLeptonpT = new TH2D ( "GenRecoLeptonpT", "Gen/Reco Matching", 80, 0, 400, 80, 0, 400 );
    h_GenRecoAntiLeptonpT = new TH2D ( "GenRecoAntiLeptonpT", "Gen/Reco Matching", 80, 0, 400, 80, 0, 400 );

    h_GenRecoTopRapidity = new TH2D ( "GenRecoTopRapidity", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 );
    h_GenRecoAntiTopRapidity = new TH2D ( "GenRecoAntiTopRapidity", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 );
    h_GenRecoToppT = new TH2D ( "GenRecoToppT", "Gen/Reco Matching", 400, 0, 400, 400, 0, 400 );
    h_GenRecoAntiToppT = new TH2D ( "GenRecoAntiToppT", "Gen/Reco Matching", 400, 0, 400, 400, 0, 400 );

    h_GenRecoTTBarRapidity = new TH2D ( "GenRecoTTBarRapidity", "Rapidity of TTbar System (HYP)", 100, -5, 5, 100, -5, 5 );
    h_GenRecoTTBarpT = new TH2D ( "GenRecoTTBarpT", "pT of TTbar System (HYP)", 500, 0, 500, 500, 0, 500 );
    h_GenRecoTTBarMass = new TH2D ( "GenRecoTTBarMass", "Mass of TTbar System (HYP)", 2000, 0, 2000, 2000, 0, 2000 );
    h_GenRecoLLBarMass = new TH2D ( "GenRecoLLBarMass", "Mass of LLbar System (HYP)", 500, 0, 1000, 500, 0, 1000 );
    h_GenRecoLLBarpT = new TH2D ( "GenRecoLLBarpT", "pT of LLbar System (HYP)", 200, 0, 1000, 200, 0, 1000 );

    h_NJetMatching = new TH1D ( "NJetMatching", "NJet Gen/Reco Matching", 5, 0, 5 );

    h_GenRecoLLBarDPhi = new TH2D ( "GenRecoLLBarDPhi", "Gen/Reco Matching", 112, -0.1, 3.25, 112, -0.1, 3.25 );
    h_GenRecoLeptonBjetMass = new TH2D ( "GenRecoLeptonBjetMass", "Gen/Reco Matching", 500, 0, 1000, 500, 0, 1000 );
    h_GenRecoAntiLeptonBjetMass = new TH2D ( "GenRecoAntiLeptonBjetMass", "Gen/Reco Matching", 500, 0, 1000, 500, 0, 1000 );
    h_GenRecoJetMult = new TH2D ( "GenRecoJetMult", "Gen/REco Matching", 26, -0.5, 25.5, 26, -0.5, 25.5 );

    h_HypLLBarDPhi = new TH1D ( "HypLLBarDPhi", "#Delta#phi(Lep, AntiLep) (HYP)", 111, -0.1, 3.25 );
    h_HypLeptonBjetMass = new TH1D ( "HypLeptonBjetMass", "Mass(Lep, AntiBJet) (HYP)", 500, 0, 1000 );
    h_HypAntiLeptonBjetMass = new TH1D ( "HypAntiLeptonBjetMass", "Mass(AntiLep, BJet) (HYP)", 500, 0, 1000 );
    h_HypJetMult = new TH1D ( "HypJetMult", "Jet Multiplicity (HYP)", 26, -0.5, 25.5 );

    h_VisGenLLBarDPhi = new TH1D ( "VisGenLLBarDPhi", "#Delta #Phi (Lep, AntiLep) (VisGEN)", 112, -0.1, 3.25 );
    h_VisGenLeptonBjetMass = new TH1D ( "VisGenLeptonBjetMass", "M(Lep, AntiBJet) (VisGEN)", 500, 0, 1000 );
    h_VisGenAntiLeptonBjetMass = new TH1D ( "VisGenAntiLeptonBjetMass", "M(AntiLep, BJet) (VisGEN)", 500, 0, 1000 );
    h_VisGenJetMult = new TH1D ( "VisGenJetMult", "Jet Multiplicty (VisGEN)", 26, -0.5, 25.5 );

    h_HypLLBarpTDPhi = new TH2D ( "HypLLBarpTDPhi", "DiLep: pT(ll) vs DPhi(ll)", 500, 0, 1000, 112, -0.1, 3.25 );

    h_BTagSF = new TH1D ( "BTagSF", "BTagging SF per event", 100 , 0.95, 1.05 );
    h_BTagSF->Sumw2();
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
    
    if (isMC) { //still have lumi weights for old plotterclass
        weightGenerator *= lumiWeight;
        
        if (pureweighter) {
            weightPU = pureweighter->getPUweight(vertMultiTrue);
        } else {
            if ( systematic == "PU_UP" ) {
                weightPU = weightPU_Up;   //only for PU systematic run
            } else if ( systematic == "PU_DOWN" ) {
                weightPU = weightPU_Down;   //only for PU systematic run
            }
        }
    }

    double weightLepSF = isMC ? leptonSF : 1;
    double trigEFF = 1.0; //set trigger efficiency scale factor to 1.
    
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
                    h_VisGenLeptonBjetMass->Fill(( *GenLepton + allGenJets->at(AntiBHadronIndex) ).M(), trueLevelWeight );
                    h_VisGenAntiLeptonBjetMass->Fill(( *GenAntiLepton + allGenJets->at(BHadronIndex) ).M(), trueLevelWeight );
                    h_VisGenJetMult->Fill(allGenJets->size(), trueLevelWeight );

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
    }//for visible top events

    //===CUT===
    // check if event was triggered (only needed for the signal sample which does not
    // contain trigger preselection cuts)
    if (false && isTtbarPlusTauSample) {
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
    
    //First control plots after dilepton selection (without Z cut)
    double weight = weightGenerator*trigEFF*weightLepSF;
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
    bool hasMetOrEmu = channel == "emu" || *(metEt->begin()) > 30;
    bool hasBtag = BJetIndex.size() > 0;
    double weightKinFit = 1;
    double btagSF = -1; //trick: initialize to -1 to avoid calculation of the btagSF twice
    
    if ( isZregion ) {
        Looseh1->Fill(dilepton.M(), weight);
        if ( hasJets && hasMetOrEmu && hasBtag && hasSolution) {
            btagSF = isMC ? calculateBtagSF() : 1;
            double fullWeights = weightGenerator*weightPU*weightLepSF*btagSF*trigEFF*weightKinFit;
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

    if (btagSF == -1) btagSF = isMC ? calculateBtagSF() : 1; //avoid calculation of the btagSF twice
    weight *= btagSF;
    h_BTagSF->Fill(btagSF );                    
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

    h_MET->Fill(*(metEt->begin()), weight);

    //loop over both leptons
    FillLeptonHisto(LeadLeptonNumber, weight);
    FillLeptonHisto(NLeadLeptonNumber, weight);
                
    //=== CUT ===
    //Require at least one solution for the kinematic event reconstruction
    if (!hasSolution) return kTRUE;

    weight *= weightKinFit;
    h_step9->Fill(1, weight);
    h_jetMultiXSec->Fill(jet->size(), weight);
    h_jetMultiNoPU->Fill(jet->size(), weight / weightPU );
    h_diLepMassFull_fullSel->Fill(dilepton.M(), weight);

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
    h_RecoJetMult->Fill(jet->size(), recoWeight);    

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
    h_HypLeptonBjetMass->Fill(( HypLepton->at(solutionIndex) + HypAntiBJet->at(solutionIndex) ).M(), weight);
    h_HypAntiLeptonBjetMass->Fill(( HypAntiLepton->at(solutionIndex) + HypBJet->at(solutionIndex) ).M(), weight);
    h_HypJetMult->Fill(jet->size(), weight);

    h_HypLLBarpTDPhi->Fill(hypllbar.Pt(), 
                           abs(DeltaPhi(HypLepton->at(solutionIndex), HypAntiLepton->at(solutionIndex))),
                           weight);
    
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
        h_GenRecoLeptonBjetMass->Fill(( HypLepton->at(solutionIndex) + HypAntiBJet->at(solutionIndex) ).M(), ( *GenLepton+allGenJets->at(AntiBHadronIndex) ).M(), weight );
    } else {
        h_GenRecoAntiBJetpT->Fill(HypAntiBJet->at(solutionIndex).Pt(), -1000., weight );
        h_GenRecoAntiBJetRapidity->Fill(HypAntiBJet->at(solutionIndex).Rapidity(), -1000., weight );
        h_GenRecoAntiBJetEta->Fill(HypAntiBJet->at(solutionIndex).Eta(), -1000., weight );
        h_GenRecoLeptonBjetMass->Fill(( HypLepton->at(solutionIndex) + HypAntiBJet->at(solutionIndex) ).M(), -1000., weight );
    }

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

    return kTRUE;
}

void Analysis::FillLeptonHisto(size_t i, double weight) {
    if ( lepType->at(i) == LEP_TYPE_ELECTRON ) {
        h_ElectronpT->Fill(lepton->at(i).Pt(), weight);
        h_ElectronEta->Fill(lepton->at(i).Eta(), weight);
    }
    if ( lepType->at(i) == LEP_TYPE_MUON ) {
        h_MuonpT->Fill(lepton->at(i).Pt(), weight);
        h_MuonEta->Fill(lepton->at(i).Eta(), weight);
    }
}


void Analysis::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.
    
    fOutput->Add(h_GenRecoLeptonEta);
    fOutput->Add(h_GenRecoAntiLeptonEta);
    fOutput->Add(h_GenRecoLeptonpT);
    fOutput->Add(h_GenRecoAntiLeptonpT);

    fOutput->Add(h_GenRecoBJetpT);
    fOutput->Add(h_GenRecoAntiBJetpT);
    fOutput->Add(h_GenRecoBJetRapidity);
    fOutput->Add(h_GenRecoAntiBJetRapidity);
    fOutput->Add(h_GenRecoBJetEta);
    fOutput->Add(h_GenRecoAntiBJetEta);

    fOutput->Add(h_GenRecoTopRapidity);
    fOutput->Add(h_GenRecoAntiTopRapidity);
    fOutput->Add(h_GenRecoToppT);
    fOutput->Add(h_GenRecoAntiToppT);

    fOutput->Add(h_GenRecoLLBarpT);
    fOutput->Add(h_GenRecoLLBarMass);
    fOutput->Add(h_GenRecoTTBarpT);
    fOutput->Add(h_GenRecoTTBarMass);
    fOutput->Add(h_GenRecoTTBarRapidity);

    fOutput->Add(h_GenRecoLLBarDPhi);
    fOutput->Add(h_GenRecoLeptonBjetMass);
    fOutput->Add(h_GenRecoAntiLeptonBjetMass);
    fOutput->Add(h_GenRecoJetMult);

    fOutput->Add(h_NJetMatching);
    fOutput->Add(h_diLepMassFull);
    fOutput->Add(h_diLepMassFull_fullSel);
    fOutput->Add(Allh1);
    fOutput->Add(Looseh1);
    fOutput->Add(h_GenAll);
    fOutput->Add(h_VisGenAll);
    fOutput->Add(h_vertMulti);
    fOutput->Add(h_vertMulti_noPU);
    fOutput->Add(h_jetMulti);
    fOutput->Add(h_jetMulti_diLep);
    fOutput->Add(h_BjetMulti);
    fOutput->Add(h_jetMultiXSec);
    fOutput->Add(h_jetMultiNoPU);
    fOutput->Add(Zh1);
    fOutput->Add(TTh1);

    fOutput->Add(h_HypTTBarpT);
    fOutput->Add(h_HypTTBarRapidity);
    fOutput->Add(h_HypTTBarMass);
    fOutput->Add(h_HypLLBarpT);
    fOutput->Add(h_HypLLBarMass);

    fOutput->Add(h_GenTTBarMass);
    fOutput->Add(h_GenTTBarRapidity);
    fOutput->Add(h_GenTTBarpT);
    fOutput->Add(h_GenLLBarpT);
    fOutput->Add(h_GenLLBarMass);

    fOutput->Add(h_VisGenTTBarMass);
    fOutput->Add(h_VisGenTTBarRapidity);
    fOutput->Add(h_VisGenTTBarpT);
    fOutput->Add(h_VisGenLLBarpT);
    fOutput->Add(h_VisGenLLBarMass);

    fOutput->Add(h_RecoTTBarMass);
    fOutput->Add(h_RecoTTBarRapidity);
    fOutput->Add(h_RecoTTBarpT);
    fOutput->Add(h_RecoToppT);
    fOutput->Add(h_RecoAntiToppT);
    fOutput->Add(h_RecoTopRapidity);
    fOutput->Add(h_RecoAntiTopRapidity);

    fOutput->Add(h_RecoLLBarMass);
    fOutput->Add(h_RecoLLBarpT);
    fOutput->Add(h_RecoLeptonpT);
    fOutput->Add(h_RecoAntiLeptonpT);
    fOutput->Add(h_RecoLeptonEta);
    fOutput->Add(h_RecoAntiLeptonEta);

    fOutput->Add(h_RecoBJetpT);
    fOutput->Add(h_RecoAntiBJetpT);
    fOutput->Add(h_RecoBJetRapidity);
    fOutput->Add(h_RecoAntiBJetRapidity);
    fOutput->Add(h_RecoBJetEta);
    fOutput->Add(h_RecoAntiBJetEta);

    fOutput->Add(h_RecoLLBarDPhi);
    fOutput->Add(h_RecoLeptonantiBjetMass);
    fOutput->Add(h_RecoAntiLeptonBjetMass);
    fOutput->Add(h_RecoJetMult);
    
    fOutput->Add(h_jetpT);
    fOutput->Add(h_jetHT);
    fOutput->Add(h_MET);
    fOutput->Add(h_LeptonpT);
    fOutput->Add(h_LeptonpT_diLep);
    fOutput->Add(h_LeptonEta_diLep);
    fOutput->Add(h_AntiLeptonpT_diLep);
    fOutput->Add(h_AntiLeptonEta_diLep);
    fOutput->Add(h_LeptonEta);
    fOutput->Add(h_AntiLeptonpT);
    fOutput->Add(h_AntiLeptonEta);
    fOutput->Add(h_ElectronpT);
    fOutput->Add(h_ElectronEta);
    fOutput->Add(h_MuonpT);
    fOutput->Add(h_MuonEta);

    fOutput->Add(h_VisGenBJetpT);
    fOutput->Add(h_VisGenAntiBJetpT);
    fOutput->Add(h_VisGenBJetRapidity);
    fOutput->Add(h_VisGenAntiBJetRapidity);
    fOutput->Add(h_VisGenBJetEta);
    fOutput->Add(h_VisGenAntiBJetEta);

    /*  fOutput->Add(h_VisGenBQuarkpT);
    fOutput->Add(h_VisGenAntiBQuarkpT);
    fOutput->Add(h_VisGenBQuarkRapidity);
    fOutput->Add(h_VisGenAntiBQuarkRapidity);
    fOutput->Add(h_VisGenBQuarkEta);
    fOutput->Add(h_VisGenAntiBQuarkEta);
    */
    fOutput->Add(h_VisGenLeptonpT);
    fOutput->Add(h_VisGenAntiLeptonpT);
    fOutput->Add(h_VisGenLeptonEta);
    fOutput->Add(h_VisGenAntiLeptonEta);

    fOutput->Add(h_VisGenToppT);
    fOutput->Add(h_VisGenAntiToppT);
    fOutput->Add(h_VisGenTopEta);
    fOutput->Add(h_VisGenAntiTopEta);
    fOutput->Add(h_VisGenTopRapidity);
    fOutput->Add(h_VisGenAntiTopRapidity);

    fOutput->Add(h_VisGenLLBarDPhi);
    fOutput->Add(h_VisGenLeptonBjetMass);
    fOutput->Add(h_VisGenAntiLeptonBjetMass);
    fOutput->Add(h_VisGenJetMult);

    /*  fOutput->Add(h_GenBJetpT);
    fOutput->Add(h_GenAntiBJetpT);
    fOutput->Add(h_GenBJetRapidity);
    fOutput->Add(h_GenAntiBJetRapidity);
    fOutput->Add(h_GenBJetEta);
    fOutput->Add(h_GenAntiBJetEta);

    fOutput->Add(h_GenLeptonpT);
    fOutput->Add(h_GenAntiLeptonpT);
    fOutput->Add(h_GenLeptonEta);
    fOutput->Add(h_GenAntiLeptonEta);

    fOutput->Add(h_GenToppT);
    fOutput->Add(h_GenAntiToppT);
    fOutput->Add(h_GenTopEta);
    fOutput->Add(h_GenAntiTopEta);
    fOutput->Add(h_GenTopRapidity);
    fOutput->Add(h_GenAntiTopRapidity);

    fOutput->Add(h_GenLLBarDPhi);
    fOutput->Add(h_GenLeptonantiBjetMass);
    fOutput->Add(h_GenAntiLeptonBjetMass);
    fOutput->Add(h_GenJetMult);

    fOutput->Add(h_GenBQuarkpT);
    fOutput->Add(h_GenAntiBQuarkpT);
    fOutput->Add(h_GenBQuarkRapidity);
    fOutput->Add(h_GenAntiBQuarkRapidity);
    fOutput->Add(h_GenBQuarkEta);
    fOutput->Add(h_GenAntiBQuarkEta);
    */
    fOutput->Add(h_HypBJetpT);
    fOutput->Add(h_HypAntiBJetpT);
    fOutput->Add(h_HypBJetRapidity);
    fOutput->Add(h_HypAntiBJetRapidity);
    fOutput->Add(h_HypBJetEta);
    fOutput->Add(h_HypAntiBJetEta);

    fOutput->Add(h_HypLeptonpT);
    fOutput->Add(h_HypAntiLeptonpT);
    fOutput->Add(h_HypLeptonEta);

    fOutput->Add(h_HypAntiLeptonEta);

    fOutput->Add(h_HypTopRapidity);
    fOutput->Add(h_HypAntiTopRapidity);

    fOutput->Add(h_HypTopMass);
    fOutput->Add(h_HypAntiTopMass);
    fOutput->Add(h_HypToppT);
    fOutput->Add(h_HypAntiToppT);
    fOutput->Add(h_HypTopEta);
    fOutput->Add(h_HypAntiTopEta);

    fOutput->Add(h_HypLLBarDPhi);
    fOutput->Add(h_HypLeptonBjetMass);
    fOutput->Add(h_HypAntiLeptonBjetMass);
    fOutput->Add(h_HypJetMult);

    fOutput->Add(h_step5);
    fOutput->Add(h_step6);
    fOutput->Add(h_step7);
    fOutput->Add(h_step8);
    fOutput->Add(h_step9);


    fOutput->Add(h_BTagSF);

    fOutput->Add(h_HypLLBarpTDPhi);
    
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
    }
    weightedEvents->Write();
    TObjString(channel).Write("channelName");
    TObjString(systematic).Write("systematicsName");
    TObjString(samplename).Write("sampleName");
    TObjString(isSignal ? "1" : "0").Write("isSignal");
    TObjString(isMC ? "1" : "0").Write("isMC");
    f.Close();
    
    fOutput->SetOwner();
    fOutput->Clear();
    cout<<"Created: "<<f_savename<<endl;
}

double Analysis::BJetSF ( double pt, double eta )
{
    //CSVL b-jet SF
    //From BTV-11-004 and https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-mujet_payload.txt

    if ( abs(eta) > 2.4 ) {
        cout<<"Jet Eta out of the selected range. Check it"<<endl;
        return 0.0;
    }
    
    if ( pt < 30 ) pt = 30;
    if ( pt > 670 ) pt = 670;

    return 1.02658*((1.+(0.0195388*pt))/(1.+(0.0209145*pt)));
}

double Analysis::CJetSF ( double pt, double eta )
{
    //CSVL c-jet SF
    //From BTV-11-004 and https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2011_Data_and_MC
    return BJetSF( pt, eta );
}

double Analysis::LJetSF ( double pt, double eta )
{
    //CSVL ligth jet mistag SF.
    //From BTV-11-004 and https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlightFuncs.C

    double eta_abs = abs(eta);
    if (eta_abs > 2.4) {
        cout<<"There is a jet out of the selected ETA region. Check that!!!!!";
        return 1;
    }
    if ( pt > 670 ) {
        return ((0.956023+(0.000825106*pt))+(-3.18828e-06*(pt*pt)))+(2.81787e-09*(pt*(pt*pt)));
    } else {
        if ( eta_abs <= 0.5 ) {
            return ((0.994425+(-8.66392e-05*pt))+(-3.03813e-08*(pt*pt)))+(-3.52151e-10*(pt*(pt*pt)));
        } else if ( eta_abs <= 1.0 ) {
            return ((0.998088+(6.94916e-05*pt))+(-4.82731e-07*(pt*pt)))+(1.63506e-10*(pt*(pt*pt)));
        } else if ( eta_abs <= 1.5 ) {
            return ((1.00294+(0.000289844*pt))+(-7.9845e-07*(pt*pt)))+(5.38525e-10*(pt*(pt*pt)));
        } else {
            return ((0.979816+(0.000138797*pt))+(-3.14503e-07*(pt*pt)))+(2.38124e-10*(pt*(pt*pt)));
        }
    }
}

double Analysis::BJetSFAbsErr ( int ptbin )
{
    //c- and l-jets errors are not necessary for the calculation and are not implemented. If needed go to https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2011_Data_and_MC

    //b-jet pt ranges {0, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 670};
    //this pt range matches the pTEff histogram!!

    double SFb_error[] = {0.1388743, 0.0188743, 0.0161816, 0.0139824, 0.0152644, 0.0161226, 0.0157396, 0.0161619, 0.0168747, 0.0257175, 0.026424, 0.0264928, 0.0315127, 0.030734, 0.0438259 };

    if ( ptbin > 14 ) {
        return 2 * SFb_error[14];
    } else {
        return SFb_error[ptbin];
    };
}

void Analysis::SetBTagFile(TString btagFile)
{
    this->btagFile = btagFile;
}

void Analysis::SetChannel(TString channel)
{
    this->channel = channel;
    if (channel == "emu") leptonSF = 0.957;
    else if (channel == "ee") leptonSF = 0.935;
    else leptonSF = 0.987;

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
    lumiWeight = 5100*SampleXSection(samplename)/weightedEvents->GetBinContent(1);
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
    metEt = 0;
    metPhi = 0;
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
    fChain->SetBranchAddress("metEt", &metEt, &b_metEt );
    fChain->SetBranchAddress("metPhi", &metPhi, &b_metPhi );
    fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber );
    fChain->SetBranchAddress("triggerBits", &triggerBits, &b_triggerBits );
    fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock );
    fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber );
    fChain->SetBranchAddress("weightGenerator", &weightGenerator, &b_weightGenerator );
    fChain->SetBranchAddress("weightPU", &weightPU, &b_weightPU );
    fChain->SetBranchAddress("weightPU_Up", &weightPU_Up, &b_weightPU_Up );
    fChain->SetBranchAddress("weightPU_Down", &weightPU_Down, &b_weightPU_Down );
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

    b_metEt->GetEntry(entry); //!
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
    b_metPhi->GetEntry(entry); //!
    b_runNumber->GetEntry(entry); //!
    b_triggerBits->GetEntry(entry); //!
    b_lumiBlock->GetEntry(entry); //!
    b_weightGenerator->GetEntry(entry); //!
    b_weightPU->GetEntry(entry); //!
    b_weightPU_Up->GetEntry(entry); //!
    b_weightPU_Down->GetEntry(entry); //!
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

double Analysis::calculateBtagSF()
{
    return 1;
    
    if (!bEff) return 1; //no btag file given, so return 1
    
    //pt efficiency median value can be obtained running and reading the output of: root -l -b -q CalcMedian.C
    const double ptmedian = 65;
    const double etamedian = 0.75;

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
            if ( ( *jetType ) [i] == 2 ) { //b-quark
                eff=bEff->GetBinContent ( ptbin, etabin );
                SFPerJet=BJetSF ( pt, eta );
                SF_Error = BJetSFAbsErr ( ptbin );
            } else if ( ( *jetType ) [i] == 1 ) { //c-quark
                SFPerJet=CJetSF ( pt, eta );
                eff=cEff->GetBinContent ( ptbin, etabin );
            } else if ( ( *jetType ) [i] == 0 ) { //l-quark
                SFPerJet=LJetSF ( pt, eta );
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
                if ( pt>ptmedian )  {
                    sf = SFPerJet - 0.5 * SF_Error;
                } else {
                    sf = SFPerJet + 0.5 * SF_Error;
                }
            }
            else if ( systematic == "BTAG_PT_DOWN" ) {
                if ( pt>ptmedian )  {
                    sf = SFPerJet + 0.5 * SF_Error;
                } else {
                    sf = SFPerJet - 0.5 * SF_Error;
                }
            }
            else if ( systematic == "BTAG_ETA_UP" ) {
                if ( eta>etamedian )  {
                    sf = SFPerJet - 0.5 * SF_Error;
                } else {
                    sf = SFPerJet + 0.5 * SF_Error;
                }
            }
            else if ( systematic == "BTAG_ETA_DOWN" ) {
                if ( eta>etamedian )  {
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


