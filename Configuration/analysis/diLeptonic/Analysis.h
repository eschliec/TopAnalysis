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
#include <map>
#include <utility>
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
#include "PUReweighter.h"

class HistoListReader;
using namespace std;

class Analysis : public TSelector
{
    static const int LEP_TYPE_ELECTRON = -1;
    static const int LEP_TYPE_MUON = 1;
    
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           EventCounter;
    VLV             *leptons;
    vector<int>     *lepQ;
    vector<int>     *lepType;
    vector<double>  *lepPfIso;
    vector<double>  *lepCombIso;
    VLV             *jets;
    vector<double>  *jetBTagTCHE;
    vector<double>  *jetBTagCSV;
    vector<double>  *jetBTagSSVHE;
    vector<int>     *jetType;
    //VLV             *genJets;
    LV              *met;
    UInt_t          runNumber;
    UInt_t          lumiBlock;
    UInt_t          eventNumber;
    UInt_t          triggerBits;
    Double_t        weightGenerator;
    vector<double>  *weightPDF;
    Int_t           vertMulti;
    Int_t           vertMultiTrue;

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
    LV              *GenMet;
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
    //TBranch        *b_genJet;   //!
    TBranch        *b_met;   //!
    TBranch        *b_runNumber;   //!
    TBranch        *b_lumiBlock;   //!
    TBranch        *b_eventNumber;   //!
    TBranch        *b_triggerBits;   //!
    TBranch        *b_weightGenerator;   //!
    TBranch        *b_weightPDF;
    TBranch        *b_vertMulti;   //!
    TBranch        *b_vertMultiTrue;

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
    TBranch        *b_GenMet;

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
    TH2 *h_GenRecoMet;
    
    TH1 *h_NJetMatching;

    TH1 *Looseh1, *Allh1, *Zh1, *TTh1, *h_GenAll, *h_jetMulti, *h_jetMulti_noBTag,
      *h_jetMulti_diLep, *h_BjetMulti, *h_BjetMulti_noBTag,*h_jetMultiXSec,*h_jetMultiAll, 
         *h_jetMultiNoPU, *h_jetMultiVisTop, *h_VisGenAll, *h_diLepMassFull, 
         *h_diLepMassFull_fullSel;

    TH1 *h_HypTTBarMass, *h_HypTTBarRapidity, *h_HypTTBarpT;
    TH1 *h_HypLLBarMass, *h_HypLLBarpT;
    TH1 *h_HypMet;

    TH1 *h_VisGenTTBarMass,*h_VisGenTTBarRapidity,*h_VisGenTTBarpT;
    TH1 *h_VisGenTopRapidity,*h_VisGenAntiTopRapidity;
    TH1 *h_VisGenLLBarMass,*h_VisGenLLBarpT;
    TH1 *h_VisGenMet;

    TH1 *h_RecoTTBarMass, *h_RecoTTBarRapidity,*h_RecoTTBarpT;
    TH1 *h_RecoToppT,*h_RecoAntiToppT,*h_RecoTopRapidity,*h_RecoAntiTopRapidity;
    TH1 *h_RecoLLBarMass, *h_RecoLLBarpT;
    TH1 *h_RecoLeptonpT,*h_RecoAntiLeptonpT,*h_RecoLeptonEta,*h_RecoAntiLeptonEta;
    TH1 *h_RecoBJetpT,*h_RecoAntiBJetpT, *h_RecoBJetRapidity,*h_RecoAntiBJetRapidity,*h_RecoBJetEta,*h_RecoAntiBJetEta;
    TH1 *h_RecoMet;
    
    TH2 *h_GenRecoHT;
    TH1 *h_VisGenHT, *h_HypHT, *h_RecoHT;

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
    
    TH1* h_HypTopptSonnenschein;

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

    TH2 *h_GenRecoLLBarDPhi, *h_GenRecoLeptonantiBjetMass, *h_GenRecoAntiLeptonBjetMass, *h_GenRecoJetMult;
    TH1 *h_VisGenLLBarDPhi,  *h_VisGenLeptonantiBjetMass,  *h_VisGenAntiLeptonBjetMass,  *h_VisGenJetMult;
    TH1 *h_HypLLBarDPhi,     *h_HypLeptonantiBjetMass,     *h_HypAntiLeptonBjetMass,     *h_HypJetMult;
    TH1 *h_RecoLLBarDPhi,    *h_RecoLeptonantiBjetMass,    *h_RecoAntiLeptonBjetMass,    *h_RecoJetMult;

    TH1 *h_HypToppTLead,    *h_HypToppTNLead,    *h_HypTopRapidityLead, *h_HypTopRapidityNLead, *h_HypTopMassLead, *h_HypTopMassNLead;
    TH1 *h_HypLeptonpTLead, *h_HypLeptonpTNLead, *h_HypLeptonEtaLead,   *h_HypLeptonEtaNLead;
    TH1 *h_HypBJetpTLead,   *h_HypBJetpTNLead,   *h_HypBJetEtaLead,     *h_HypBJetEtaNLead;

    TH1 *h_RecoToppTLead,    *h_RecoToppTNLead,    *h_RecoTopRapidityLead, *h_RecoTopRapidityNLead, *h_RecoTopMassLead, *h_RecoTopMassNLead;
    TH1 *h_RecoLeptonpTLead, *h_RecoLeptonpTNLead, *h_RecoLeptonEtaLead,   *h_RecoLeptonEtaNLead;
    TH1 *h_RecoBJetpTLead,   *h_RecoBJetpTNLead,   *h_RecoBJetEtaLead,     *h_RecoBJetEtaNLead;

    TH1 *h_VisGenToppTLead,    *h_VisGenToppTNLead,    *h_VisGenTopRapidityLead, *h_VisGenTopRapidityNLead, *h_VisGenTopMassLead, *h_VisGenTopMassNLead;
    TH1 *h_VisGenLeptonpTLead, *h_VisGenLeptonpTNLead, *h_VisGenLeptonEtaLead,   *h_VisGenLeptonEtaNLead;
    TH1 *h_VisGenBJetpTLead,   *h_VisGenBJetpTNLead,   *h_VisGenBJetEtaLead,     *h_VisGenBJetEtaNLead;

    TH2 *h_GenRecoToppTLead,    *h_GenRecoToppTNLead,    *h_GenRecoTopRapidityLead, *h_GenRecoTopRapidityNLead, *h_GenRecoTopMassLead, *h_GenRecoTopMassNLead;
    TH2 *h_GenRecoLeptonpTLead, *h_GenRecoLeptonpTNLead, *h_GenRecoLeptonEtaLead,   *h_GenRecoLeptonEtaNLead;
    TH2 *h_GenRecoBJetpTLead,   *h_GenRecoBJetpTNLead,   *h_GenRecoBJetEtaLead,     *h_GenRecoBJetEtaNLead;

    
//     //Begin: Plots for Carmen
//     TH1 *h_RecoLeadingJetpT,    *h_RecoNLeadingJetpT,    *h_RecoLeadingJetEta,    *h_RecoNLeadingJetEta;
//     TH1 *h_HypLeadingJetpT,     *h_HypNLeadingJetpT,     *h_HypLeadingJetEta,     *h_HypNLeadingJetEta;
//     TH2 *h_GenRecoLeadingJetpT, *h_GenRecoLeadingJetEta, *h_GenRecoNLeadingJetpT, *h_GenRecoNLeadingJetEta;
//     TH1 *h_VisGenLeadingJetpT,  *h_VisGenLeadingJetEta,  *h_VisGenNLeadingJetpT,  *h_VisGenNLeadingJetEta;
// 
//     //Begin: Plots for Carmen
//     TH1 *h_RecoExtraJetpT,  *h_HypExtraJetpT, *h_VisGenExtraJetpT, *h_RecoExtraJetEta, *h_HypExtraJetEta, *h_VisGenExtraJetEta;
//     TH1 *h_RecoExtraJetpT2, *h_HypExtraJetpT2, *h_VisGenExtraJetpT2, *h_RecoExtraJetEta2, *h_HypExtraJetEta2, *h_VisGenExtraJetEta2;
//     TH1 *h_RecoExtraJetpT3, *h_HypExtraJetpT3, *h_VisGenExtraJetpT3, *h_RecoExtraJetEta3, *h_HypExtraJetEta3, *h_VisGenExtraJetEta3;
//     TH1 *h_RecoExtraJetpT4, *h_HypExtraJetpT4, *h_VisGenExtraJetpT4, *h_RecoExtraJetEta4, *h_HypExtraJetEta4, *h_VisGenExtraJetEta4;
//     TH2 *h_GenRecoExtraJetpT, *h_GenRecoExtraJetEta, *h_GenRecoExtraJetpT2, *h_GenRecoExtraJetEta2, *h_GenRecoExtraJetpT3, *h_GenRecoExtraJetEta3, *h_GenRecoExtraJetpT4, *h_GenRecoExtraJetEta4;
// 
//     TH1 *h_RecoJetMultpt40, *h_HypJetMultpt40, *h_VisGenJetMultpt40, *h_RecoJetMultpt60, *h_HypJetMultpt60, *h_VisGenJetMultpt60;
//     TH2 *h_GenRecoJetMultpt40, *h_GenRecoJetMultpt60, *h_GenRecoJetMultQ0, *h_GenRecoJetMultTotal;
//     //End: Plots for Carmen
    
    // BEGIN of btag SF stuff
    TH2 *h_bjets, *h_btaggedjets;
    TH2 *h_cjets, *h_ctaggedjets;
    TH2 *h_ljets, *h_ltaggedjets;
    
    //btag calculation
    TH2 *bEff, *cEff, *lEff;
    TH1 *h_BTagSF;

    double Median(TH1 *); 
    double BJetSF ( double, double );
    double CJetSF ( double, double );
    double LJetSF ( double, double );
    double BJetSFAbsErr ( int );
    double btag_ptmedian, btag_etamedian;
    // END of btag SF stuff
    
    //other SF
    double lumiWeight; //needed while using old plotterclass
    
    TH1* h_TrigSFeta;
    TH2* h_LepIDSFpteta;
    
    void prepareTriggerSF();
    void prepareBtagSF();
    void prepareLeptonIDSF();
    double getTriggerSF(const LV& lep1, const LV& lep2);
    double getLeptonIDSF(const LV& lep1, const LV& lep2);
    double get2DSF(TH2* histo, const double x, const double y);
    
    // store the object in the output list and return it
    template<class T> T* store(T* obj);
    
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
    int pdf_no;
    TH1* weightedEvents;
    PUReweighter *pureweighter;
    
    double calculateBtagSF();
    double getJetHT(const VLV& jet, int pt_cut);
    
    //order two LorentzVectors by transverse momentum
    //first two parameters are output, 3rd and 4th input
    void orderLVByPt(LV &leading, LV &Nleading, const LV &lv1, const LV &lv2);
    
    //binnedControlPlots contains:
    //map of name of differential distribution
    // -> pair( histogram with the binning of the differential distribution,
    //          vector(bin) -> map( control plot name -> TH1*))
    std::map<std::string, std::pair<TH1*, std::vector<std::map<std::string, TH1*> > > > *binnedControlPlots;
    
    // Create Nbins control plots for the differential distribution h_differential
    // Use h_control for the control plot name and binning
    void CreateBinnedControlPlots(TH1* h_differential, TH1* h_control);
    
    // h: differential distribution histogram
    // binvalue: the value of the quantity in the differential distribution histogram
    // the control plot histogram
    // the value for the control plot
    // weight: event weight
    void FillBinnedControlPlot(TH1* h_differential, double binvalue, 
                               TH1 *h_control, double value, double weight);
        
public:
    Analysis ( TTree * = 0 ) : runViaTau {0}, pureweighter {nullptr} {};
    void SetBTagFile(TString btagFile);
    void SetChannel(TString channel);
    void SetSignal(bool isSignal);
    void SetSystematic(TString systematic);
    void SetSamplename(TString samplename);
    void SetOutputfilename(TString outputfilename);
    void SetMC(bool isMC);
    void SetWeightedEvents(TH1* weightedEvents);
    void SetRunViaTau(bool runViaTau);
    void SetPUReweighter(PUReweighter *pu);
    void SetPDF(int pdf_no);
    ClassDef ( Analysis,0 );
};

#endif
