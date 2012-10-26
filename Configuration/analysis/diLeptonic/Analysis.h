#ifndef Analysis_h
#define Analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TExec.h"
#include "TStyle.h"
#include "TMath.h"
#include "TROOT.h"
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "classes.h"

using namespace std;

class Analysis : public TSelector
{
    static const int LEP_TYPE_ELECTRON = -1;
    static const int LEP_TYPE_MUON = 1;
    
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           EventCounter;
    VLV             *lepton;
    vector<int>     *lepQ;
    vector<int>     *lepType;
    vector<double>  *lepPfIso;
    vector<double>  *lepCombIso;
    VLV             *jet;
    vector<double>  *jetBTagTCHE;
    vector<double>  *jetBTagCSV;
    vector<double>  *jetBTagSSVHE;
    vector<int>     *jetType;
    VLV             *genJet;
    vector<double>  *metEt;
    vector<double>  *metPhi;
    UInt_t          runNumber;
    UInt_t          lumiBlock;
    UInt_t          eventNumber;
    UInt_t          triggerBits;
    Double_t        weightGenerator;
    Double_t        weightPU;
    Double_t        weightPU_Up;
    Double_t        weightPU_Down;
    Int_t           vertMulti;

    LV              *GenWPlus;
    LV              *GenWMinus;
    LV              *GenNeutrino;
    LV              *GenAntiNeutrino;
    LV              *GenB;
    LV              *GenAntiB;
    LV              *GenLepton;
    LV              *GenAntiLepton;
    LV              *GenTop;
    LV              *GenAntiTop;
    VLV             *allGenJets;
    vector<int>     *BHadJetIndex;
    vector<int>     *AntiBHadJetIndex;

    VLV             *BHadrons;
    VLV             *AntiBHadrons;
    vector<bool>    *BHadronFromTopB;
    vector<bool>    *AntiBHadronFromTopB;
    vector<int>     *BHadronVsJet;
    vector<int>     *AntiBHadronVsJet;

    /*   VLV           *BHadronJet;
    VLV           *AntiBHadronJet;
    */
    VLV             *HypTop;
    VLV             *HypAntiTop;
    VLV             *HypLepton;
    VLV             *HypAntiLepton;
    VLV             *HypNeutrino;
    VLV             *HypAntiNeutrino;
    VLV             *HypBJet;
    VLV             *HypAntiBJet;
//     VLV           *HypWPlus;
//     VLV           *HypWMinus;
    vector<int>     *HypJet0index;
    vector<int>     *HypJet1index;
    Int_t           decayMode;

    // List of branches
    TBranch        *b_lepton;   //!
    TBranch        *b_lepQ;   //!
    TBranch        *b_lepType;   //!
    TBranch        *b_lepPfIso;   //!
    TBranch        *b_lepCombIso;   //!
    TBranch        *b_jet;   //!
    TBranch        *b_jetBTagTCHE;   //!
    TBranch        *b_jetBTagCSV;   //!
    TBranch        *b_jetBTagSSVHE;   //!
    TBranch        *b_jetType;   //!
    TBranch        *b_genJet;   //!
    TBranch        *b_metEt;   //!
    TBranch        *b_metPhi;   //!
    TBranch        *b_runNumber;   //!
    TBranch        *b_lumiBlock;   //!
    TBranch        *b_eventNumber;   //!
    TBranch        *b_triggerBits;   //!
    TBranch        *b_weightGenerator;   //!
    TBranch        *b_weightPU;   //!
    TBranch        *b_weightPU_Up;   //!
    TBranch        *b_weightPU_Down;   //!
    TBranch        *b_vertMulti;   //!

    TBranch        *b_GenTop;   //!
    TBranch        *b_GenAntiTop;   //!
    TBranch        *b_GenLepton;   //!
    TBranch        *b_GenAntiLepton;   //!
    TBranch        *b_GenNeutrino;   //!
    TBranch        *b_GenAntiNeutrino;   //!
    TBranch        *b_GenB;   //!
    TBranch        *b_GenAntiB;   //!
    TBranch        *b_GenWPlus;   //!
    TBranch        *b_GenWMinus;   //!
    TBranch        *b_allGenJets;   //!
    TBranch        *b_BHadJetIndex;   //!
    TBranch        *b_AntiBHadJetIndex;   //!

    TBranch        *b_BHadrons;   //!
    TBranch        *b_AntiBHadrons;   //!
    TBranch        *b_BHadronFromTopB;   //!
    TBranch        *b_AntiBHadronFromTopB;   //!
    TBranch        *b_BHadronVsJet;   //!
    TBranch        *b_AntiBHadronVsJet;   //!

    /*   TBranch	  *b_BHadronJet_;   //!
    TBranch	  *b_AntiBHadronJet_;   //!
    */

    TBranch        *b_HypTop;   //!
    TBranch        *b_HypAntiTop;   //!
    TBranch        *b_HypLepton;   //!
    TBranch        *b_HypAntiLepton;   //!
    TBranch        *b_HypNeutrino;   //!
    TBranch        *b_HypAntiNeutrino;   //!
    TBranch        *b_HypB;   //!
    TBranch        *b_HypAntiB;   //!
    TBranch        *b_HypWPlus;   //!
    TBranch        *b_HypWMinus;   //!
    TBranch        *b_HypJet0index;   //!
    TBranch        *b_HypJet1index;   //!
    TBranch        *b_decayMode;   //!

    virtual ~Analysis() { }
    virtual Int_t   Version() const {
        return 3;
    }
    virtual void    Begin ( TTree* );
    virtual void    SlaveBegin ( TTree* );
    virtual void    Init ( TTree *tree );
    virtual Bool_t  Notify();
    virtual Bool_t  Process ( Long64_t entry );
    virtual Int_t   GetEntry ( Long64_t entry, Int_t getall = 0 ) {
        return fChain ? fChain->GetTree()->GetEntry ( entry, getall ) : 0;
    }
    virtual void    SetOption ( const char *option ) {
        fOption = option;
    }
    virtual void    SetObject ( TObject *obj ) {
        fObject = obj;
    }
    virtual void    SetInputList ( TList *input ) {
        fInput = input;
    }
    virtual TList  *GetOutputList() const {
        return fOutput;
    }
    virtual void    SlaveTerminate();
    virtual void    Terminate();

    void GetRecoBranches ( Long64_t & );
    void GetSignalBranches ( Long64_t & );

    TH2 *h_GenRecoLeptonpT,*h_GenRecoAntiLeptonpT,*h_GenRecoLeptonEta,*h_GenRecoAntiLeptonEta, *h_GenRecoLLBarMass, *h_GenRecoLLBarpT;
    TH2 *h_GenRecoBJetpT,*h_GenRecoAntiBJetpT, *h_GenRecoBJetEta,*h_GenRecoAntiBJetEta, *h_GenRecoBJetRapidity, *h_GenRecoAntiBJetRapidity;//, *h_GenRecoBJetE, *h_GenRecoAntiBJetE;;
    TH2 *h_GenRecoToppT,*h_GenRecoAntiToppT,*h_GenRecoTopRapidity,*h_GenRecoAntiTopRapidity, *h_GenRecoTTBarMass, *h_GenRecoTTBarpT, *h_GenRecoTTBarRapidity;

    TH1 *h_NJetMatching;

    TH1 *Looseh1, *Allh1, *Zh1, *TTh1, *h_GenAll, *h_jetMulti,
         *h_jetMulti_diLep, *h_BjetMulti,*h_jetMultiXSec,*h_jetMultiAll, 
         *h_jetMultiNoPU, *h_jetMultiVisTop, *h_VisGenAll, *h_diLepMassFull, 
         *h_diLepMassFull_fullSel;

    TH1 *h_HypTTBarMass, *h_HypTTBarRapidity, *h_HypTTBarpT;
    TH1 *h_HypLLBarMass, *h_HypLLBarpT;
    TH1 *h_GenTTBarMass, *h_GenTTBarRapidity, *h_GenTTBarpT;
    TH1 *h_GenLLBarMass, *h_GenLLBarpT;

    TH1 *h_VisGenTTBarMass,*h_VisGenTTBarRapidity,*h_VisGenTTBarpT;
    TH1 *h_VisGenTopRapidity,*h_VisGenAntiTopRapidity;
    TH1 *h_VisGenLLBarMass,*h_VisGenLLBarpT;

    TH1 *h_RecoTTBarMass, *h_RecoTTBarRapidity,*h_RecoTTBarpT;
    TH1 *h_RecoToppT,*h_RecoAntiToppT,*h_RecoTopRapidity,*h_RecoAntiTopRapidity;
    TH1 *h_RecoLLBarMass, *h_RecoLLBarpT;
    TH1 *h_RecoLeptonpT,*h_RecoAntiLeptonpT,*h_RecoLeptonEta,*h_RecoAntiLeptonEta;
    TH1 *h_RecoBJetpT,*h_RecoAntiBJetpT, *h_RecoBJetRapidity,*h_RecoAntiBJetRapidity,*h_RecoBJetEta,*h_RecoAntiBJetEta;
    TH1 *h_RecoLLBarDPhi, *h_RecoLeptonantiBjetMass, *h_RecoAntiLeptonBjetMass, *h_RecoJetMult;

    TH1 *h_vertMulti, *h_vertMulti_noPU, *h_MET;

    TH1 *h_jetpT,*h_jetHT;
    TH1 *h_MuonpT, *h_MuonEta;
    TH1 *h_ElectronpT, *h_ElectronEta;
    TH1 *h_LeptonpT, *h_LeptonEta;
    TH1 *h_AntiLeptonpT, *h_AntiLeptonEta;
    TH1 *h_LeptonpT_diLep, *h_LeptonEta_diLep;
    TH1 *h_AntiLeptonpT_diLep, *h_AntiLeptonEta_diLep;

    TH1 *h_HypAntiToppT, *h_HypAntiTopEta, *h_HypAntiTopMass,*h_HypAntiTopRapidity;
    TH1 *h_HypToppT, *h_HypTopEta,*h_HypTopMass, *h_HypTopRapidity ;

    TH1 *h_HypAntiBJetpT, *h_HypAntiBJetEta, *h_HypAntiBJetRapidity;
    TH1 *h_HypBJetpT, *h_HypBJetEta, *h_HypBJetRapidity;

    TH1 *h_HypAntiLeptonpT, *h_HypAntiLeptonEta;
    TH1 *h_HypLeptonpT, *h_HypLeptonEta;

    TH1 *h_step5,*h_step6,*h_step7,*h_step8,*h_step9;

    TH1 *h_VisGenAntiToppT, *h_VisGenAntiTopEta;
    TH1 *h_VisGenToppT, *h_VisGenTopEta;

    TH1 *h_VisGenAntiBJetpT, *h_VisGenAntiBJetEta, *h_VisGenAntiBJetRapidity;
    TH1 *h_VisGenBJetpT, *h_VisGenBJetEta, *h_VisGenBJetRapidity;

    TH1 *h_VisGenAntiBQuarkpT, *h_VisGenAntiBQuarkEta, *h_VisGenAntiBQuarkRapidity;
    TH1 *h_VisGenBQuarkpT, *h_VisGenBQuarkEta, *h_VisGenBQuarkRapidity;

    TH1 *h_VisGenAntiLeptonpT, *h_VisGenAntiLeptonEta;
    TH1 *h_VisGenLeptonpT, *h_VisGenLeptonEta;

    TH1 *h_GenAntiToppT, *h_GenAntiTopEta, *h_GenAntiTopRapidity;
    TH1 *h_GenToppT, *h_GenTopEta, *h_GenTopRapidity;

    TH1 *h_GenAntiBJetpT, *h_GenAntiBJetEta, *h_GenAntiBJetRapidity;
    TH1 *h_GenBJetpT, *h_GenBJetEta, *h_GenBJetRapidity;

    TH1 *h_GenAntiBQuarkpT, *h_GenAntiBQuarkEta, *h_GenAntiBQuarkRapidity;
    TH1 *h_GenBQuarkpT, *h_GenBQuarkEta, *h_GenBQuarkRapidity;

    TH1 *h_GenAntiLeptonpT, *h_GenAntiLeptonEta;
    TH1 *h_GenLeptonpT, *h_GenLeptonEta;

    TH2 *h_GenRecoLLBarDPhi, *h_GenRecoLeptonBjetMass, *h_GenRecoAntiLeptonBjetMass, *h_GenRecoJetMult;
    TH1 *h_VisGenLLBarDPhi,  *h_VisGenLeptonBjetMass,  *h_VisGenAntiLeptonBjetMass,  *h_VisGenJetMult;
    TH1 *h_HypLLBarDPhi,     *h_HypLeptonBjetMass,     *h_HypAntiLeptonBjetMass,     *h_HypJetMult;

    TH2 *h_HypLLBarpTDPhi;
    
    
    // BEGIN of btag SF stuff
    vector<double>  ptbinning, etabinning;

    TH2 *bEff, *cEff, *lEff;
    TH1 *h_BTagSF;

    double BJetSF ( double, double );
    double CJetSF ( double, double );
    double LJetSF ( double, double );
    double BJetSFAbsErr ( int );
    // END of btag SF stuff
    
    //other SF
    double leptonSF;
    double lumiWeight; //needed while using old plotterclass

    // Variables added from the outside
    TString btagFile;
    TString channel;
    TString systematic;
    TString samplename;
    bool isTtbarPlusTauSample;
    bool correctMadgraphBR;
    TString outputfilename;
    bool isSignal;
    bool isMC;
    bool getLeptonPair(size_t& LeadLeptonNumber, size_t& NLeadLeptonNumber);
    bool runViaTau;
    TH1* weightedEvents;
    
    double calculateBtagSF();
    double getJetHT(const VLV& jet, int pt_cut);
    
public:
    Analysis ( TTree * = 0 ) { runViaTau = 0; }
    void SetBTagFile(TString btagFile);
    void SetChannel(TString channel);
    void SetSignal(bool isSignal);
    void SetSystematic(TString systematic);
    void SetSamplename(TString samplename);
    void SetOutputfilename(TString outputfilename);
    void SetMC(bool isMC);
    void SetWeightedEvents(TH1* weightedEvents);
    void SetRunViaTau(bool runViaTau);
    ClassDef ( Analysis,0 );
};

#endif

