#include "basicFunctions.h"

void analyzeTopDiffXSecMonitoring(double luminosity = 1143, bool save = false, int verbose=0, TString inputFolderName="TOP2011/110819_AnalysisRun", 
				  TString dataFile= "/afs/naf.desy.de/group/cms/scratch/tophh/TOP2011/110819_AnalysisRun/analyzeDiffXData2011A_Muon_160404_167913_1fb.root"
				  , const std::string decayChannel = "muon", bool withRatioPlot = true)
{
  // set root style
	
  TStyle myStyle("HHStyle","HHStyle");
  setHHStyle(myStyle);
  TGaxis::SetMaxDigits(2);
  myStyle.cd();
  gROOT->SetStyle("HHStyle");
  gROOT->ForceStyle();
		
  //  ---
  //     name conventions
  //  ---
  // a) enumerator sample convention (as defined in basicFunctions.h)
  /*0:*/  /*1:*/  /*2:*/    /*3:*/    /*4:*/   /*5:*/    /*6:*/  /*7:*/  /*8,  9,  10*/ /* 11   ,  12     ,   13:  */
  // kSig  , kBkg  , kZjets  , kWjets  , kQCD   , kSTop   , kDiBos, kData , kWW, kWZ, kZZ, kSTops  , kSTopt  , kSToptW 
  // b) file name convention (implemented in basicFunctions.h)
  // "muonDiffXSec"+sampleName+GeneratorName+GeneratorTune+MCProductionCycle+systematicVariation+"PF.root"
  // sampleName = "Sig", "Bkg", Wjets", "Zjets", "WW", "WZ", "ZZ", "VV", "SingleTopSchannel", 
  //              "SingleTopTchannel", "SingleTopTWchannel", "QCD"
  // GeneratorName= "Mad", "Pythia6"
  // GeneratorTune= "Z2", "D6T"
  // MCProductionCycle= "Fall10"
	
  //  ---
  //     options
  //  ---
  // a) options directly entered when calling function
  // luminosity: [/pb]
  TString lumi = getTStringFromInt(roundToInt((luminosity), false));
  double lumiError = 0.045;  // relative luminosity error
  // save: save plots?
  // verbose: set detail level of output 
  // 0: no output, 1: std output 2: output for debugging
  // b) options to be configured only once
  // get the .root files from the following folder:
  //  TString inputFolder = "./diffXSecFromSignal/analysisRootFilesWithKinFit";
  TString inputFolder = "/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName;
  // see if its 2010 or 2011 data from luminosity
  TString dataSample="";
  if(luminosity<36) dataSample="2010";
  if(luminosity>36) dataSample="2011";
  // save all plots into the following folder
  TString outputFolder = "./diffXSecFromSignal/plots/"+decayChannel+"/";
  if(dataSample!="") outputFolder+=dataSample+"/";
  outputFolder += (withRatioPlot ? "monitoring/" : "monitoring/withoutRatioPlots/") ;
  // save all plots within a root file named:
  TString outputFileName="diffXSecTopSemi";
  if(decayChannel=="muon"    ) outputFileName+="Mu";
  if(decayChannel=="electron") outputFileName+="Elec";
  outputFileName+=dataSample+".root";
  // choose name of the output .pdf file
  TString pdfName="differentialXSecMonitoring"+lumi+"pb";
  /* systematicVariation: which systematic shift do you want to make? from basicFunctions.h:
     0:sysNo              1:sysLumiUp          2:sysLumiDown          3:sysJESUp      
     4:sysJESDown         5:sysJERUp           6:sysJERDown           7:sysTopScaleUp 
     8:sysTopScaleDown    9:sysVBosonScaleUp  10:sysVBosonScaleDown  11:sysTopMatchUp 
     12:sysTopMatchDown  13:sysVBosonMatchUp  14:sysVBosonMatchDown  15:sysMuEffSFup  
     16:sysMuEffSFdown   17:sysISRFSRup       18:sysISRFSRdown       19:sysPileUp    
     20:sysQCDup         21:sysQCDdown        22:sysSTopUp           23:sysSTopDown  
     24:sysBtagUp        25:sysBtagDown       26:sysShapeUp          27:sysShapeUp 
     28:sysPUup          29:sysPUdown         30:sysflatTrigSF       31:sysTrigEffSFNormUp
     32:sysTrigEffSFNormDown     33:sysTriggerEffSFShapeUpEta   34:sysTriggerEffSFShapeDownEta
     35:sysTriggerEffSFShapeUpPt 36:sysTriggerEffSFShapeDownPt  37:sysMisTagSFup     38:sysMisTagSFdown     
     39:sysDiBosUp       40:sysDiBosDown
  */
  int systematicVariation=sysNo;

  //  ---
  //     choose plots
  //  ---
  // a) list plots you would like to see ("folder/plotName") - same as in .root files (for 1D and 2D)

  TString plots1D[ ] = { 
    // (I) preselection
    // (ii) jet monitoring
    //"tightJetKinematicsPreSel/n"  ,
    //"tightJetKinematicsPreSel/en" ,
    //"tightJetKinematicsPreSel/pt" ,
    //"tightJetKinematicsPreSel/eta",
    //"tightJetKinematicsPreSel/phi",
    // (II) before btagging
    // (ii) jet monitoring
    "tightJetKinematics/n"  ,
    "tightJetKinematicsNjets1/n",
    "tightJetKinematics/en" ,
    "tightJetKinematics/pt" ,
    "tightJetKinematics/eta",
    "tightJetKinematics/phi",
    "tightJetKinematics/ht",
    "tightJetQuality/n_"    ,
    "tightJetQuality/charge",
    "tightJetQuality/nhf"   ,
    "tightJetQuality/nef"   ,
    "tightJetQuality/chf"   ,
    "tightJetQuality/cef"   ,
    "tightJetQuality/ncp"   ,
    "tightLead_0_JetKinematics/en" ,
    "tightLead_0_JetKinematics/pt" ,
    "tightLead_0_JetKinematics/eta",
    "tightLead_0_JetKinematics/phi",
    "tightLead_1_JetKinematics/en" ,
    "tightLead_1_JetKinematics/pt" ,
    "tightLead_1_JetKinematics/eta",
    "tightLead_1_JetKinematics/phi",
    "tightLead_2_JetKinematics/en" ,
    "tightLead_2_JetKinematics/pt" ,
    "tightLead_2_JetKinematics/eta",
    "tightLead_2_JetKinematics/phi",
    "tightLead_3_JetKinematics/en" ,
    "tightLead_3_JetKinematics/pt" ,
    "tightLead_3_JetKinematics/eta",
    "tightLead_3_JetKinematics/phi",
    // (iii) btag monitoring
    "tightJetQuality/btagTrkCntHighPurity",
    "tightJetQuality/btagTrkCntHighEff_"  ,
    "tightJetQuality/btagSimpleSecVtxHighEff_",
    "tightJetQuality/btagSimpleSecVtxHighPur_",
    "tightJetQuality/btagCombSecVtx_",
    "tightJetQuality/btagCombSecVtxMVA_",
    "tightJetQuality/btagJetBProbability_",
    "tightJetQuality/btagJetProbability_",
    "tightJetQuality/btagSoftMuon_",
    "tightJetQuality/btagSoftMuonByPt_",
    "tightJetQuality/btagSoftMuonByIP3d_",
    "bottomJetKinematics/n",
    // (iv) MET monitoring
    "analyzeMETMuon/metEt"   ,
    "analyzeMETMuon/metSumEt",
    // (v) Vertices and pileup
    "PUControlDistributionsBeforeBtagging/pileup",
    "PUControlDistributionsBeforeBtagging/pileup_reweighted",
    "PUControlDistributionsBeforeBtagging/pileup_reweighted3BX",
    "PUControlDistributionsBeforeBtagging/npvertex",
    "PUControlDistributionsBeforeBtagging/npvertex_reweighted",
    "PUControlDistributionsBeforeBtagging/npvertex_reweighted3BX",
    // (III) after btagging 
    // (ii) jet monitoring
    "tightJetKinematicsTagged/n"  ,
    "tightJetKinematicsTagged/pt" ,
    "tightJetKinematicsTagged/eta",
    "tightJetKinematicsTagged/phi",
    "tightJetKinematicsTagged/ht",
    "tightLead_0_JetKinematicsTagged/pt" ,
    "tightLead_0_JetKinematicsTagged/eta",
    "tightLead_1_JetKinematicsTagged/pt" ,
    "tightLead_1_JetKinematicsTagged/eta",
    "tightLead_2_JetKinematicsTagged/pt" ,
    "tightLead_2_JetKinematicsTagged/eta",
    "tightLead_3_JetKinematicsTagged/pt" ,
    "tightLead_3_JetKinematicsTagged/eta",
    "bottomLead_0_JetKinematicsTagged/pt",
    "bottomLead_1_JetKinematicsTagged/pt",
    "bottomLead_0_JetKinematicsTagged/eta",
    "bottomLead_1_JetKinematicsTagged/eta",
    // (iv) MET monitoring
    "analyzeMETMuonTagged/metEt"   ,
    "analyzeMETMuonTagged/metSumEt",
    // (v) Vertices and pileup
    "PUControlDistributionsAfterBtagging/pileup",
    "PUControlDistributionsAfterBtagging/pileup_reweighted",
    "PUControlDistributionsAfterBtagging/pileup_reweighted3BX",
    "PUControlDistributionsAfterBtagging/npvertex",
    "PUControlDistributionsAfterBtagging/npvertex_reweighted",
    "PUControlDistributionsAfterBtagging/npvertex_reweighted3BX", 
    // (III) after kinematic fit 
    "analyzeTopRecoKinematicsKinFit/topPt",
    "analyzeTopRecoKinematicsKinFit/topY",
    "analyzeTopRecoKinematicsKinFit/ttbarPt",
    "analyzeTopRecoKinematicsKinFit/ttbarY",
    "analyzeTopRecoKinematicsKinFit/ttbarMass",
    "analyzeTopRecoKinematicsKinFit/lepPt",
    "analyzeTopRecoKinematicsKinFit/lepEta"
    
  };

  TString plots1Dmu[ ] = { 
    // (I) preselection
    // (i) muon monitoring
    //"kinematicMuonQualityPreSel/nHit"   ,
    //"kinematicMuonQualityPreSel/chi2"   ,
    //"kinematicMuonQualityPreSel/dB"     ,
    //"kinematicMuonQualityPreSel/dz"     ,
    //"kinematicMuonQualityPreSel/matches",  
    //"trackMuontightJetsKinematicsPreSel/dist30_",
    //"goldenMuonQualityPreSel/relIso"    , 
    //"tightMuonKinematicsPreSel/n"       ,
    // (II) before btagging
    // (i) muon monitoring
    "tightMuonKinematics/n" ,
    "tightMuonKinematics/en" ,
    "tightMuonKinematics/pt" ,
    "tightMuonKinematics/eta",
    "tightMuonKinematics/y"  ,
    "tightMuonKinematics/phi",
    "tightMuonQuality/nHit"   ,
    "tightMuonQuality/chi2"   ,
    "tightMuonQuality/dB"     ,
    "tightMuonQuality/dz"     ,
    "tightMuonQuality/ecalEn" ,
    "tightMuonQuality/hcalEn" ,
    "tightMuonQuality/relIso" ,
    "tightMuonQuality/matches",  			 
    // (III) after btagging 
    // (i) muon monitoring
    "tightMuonQualityTagged/relIso" ,
    "tightMuonKinematicsTagged/pt" ,
    "tightMuonKinematicsTagged/eta",
    "tightMuonKinematicsTagged/phi"
  };
	
  TString plots1De[ ] = { 
    // (ib) electron monitoring
    "tightElectronKinematics/n" ,
    "tightElectronKinematics/en" ,
    "tightElectronKinematics/et" ,
    "tightElectronKinematics/eta",
    "tightElectronKinematics/phi",
    "tightElectronQuality/etaSC"  ,
    "tightElectronQuality/dB"     ,
    "tightElectronQuality/simpleEleId70cIso", 
    "tightElectronQuality/nHitsInner",
    "tightElectronQuality/convDcot" ,
    "tightElectronQuality/convDist" ,
    "tightElectronQuality/relIso"   ,
    // (ib) electron monitoring
    "tightElectronQualityTagged/relIso",
    "tightElectronKinematicsTagged/et" ,
    "tightElectronKinematicsTagged/eta",
    "tightElectronKinematicsTagged/phi"
  };
	
  TString plots2D[ ] = { 
  };	
	
  // b) list plot axes style
  // 1D: "x-axis title"/"y-axis title"/log/rebin-factor
  // log = 0 or 1 for linear or logarithmic axis 
	
  TString axisLabel1D[ ] = { 
    // (I) preselection
    // (ii) jet monitoring
    //"N_{jets}/events/0/1",
    //"E(jets)/jets/1/1",
    //"p_{t}(jets)/jets/1/1",
    //"#eta(jets)/jets/0/5",
    //"#phi(jets)/jets/0/10",
    // (II) before btagging
    // (ii) jet monitoring
    "N_{jets}/events/1/1",
    "N_{jets}/events/1/1",
    "E(jets) [GeV]/jets/1/1",
    "p_{t}(jets) #left[#frac{GeV}{c}#right]/jets/1/1",
    "#eta(jets)/jets/0/5",
    "#phi(jets)/jets/0/10",
    "H_{T} #left[#frac{GeV}{c}#right]/events/0/50",
    "N(jet constituents)/jets/0/10",
    "jet charge/jets/0/10"         ,
    "neutral hadron fraction (jets)/jets/1/1"         ,
    "neutral electromagnetic fraction (jets)/jets/0/2",
    "charged hadron fraction (jets)/jets/0/1"         ,
    "charged electromagnetic fraction (jets)/jets/1/1",
    "N_{charged particles} (jets)/jets/0/2"           ,
    "E(lead 1^{st} jet) [GeV]/events/1/2",
    "p_{t}(lead 1^{st} jet) #left[#frac{GeV}{c}#right]/events/0/5",
    "#eta(lead 1^{st} jet)/events/0/5",
    "#phi(lead 1^{st} jet)/events/0/10",
    "E(lead 2^{nd} jet) [GeV]/events/1/2",
    "p_{t}(lead 2^{nd} jet) #left[#frac{GeV}{c}#right]/events/0/5",
    "#eta(lead 2^{nd} jet)/events/0/5",
    "#phi(lead 2^{nd} jet)/events/0/10",
    "E(lead 3^{rd} jet) [GeV]/events/1/2",
    "p_{t}(lead 3^{rd} jet) #left[#frac{GeV}{c}#right]/events/0/5",
    "#eta(lead 3^{rd} jet)/events/0/5",
    "#phi(lead 3^{rd} jet)/events/0/10",
    "E(lead 4^{th} jet) [GeV]/events/1/2",
    "p_{t}(lead 4^{th} jet) #left[#frac{GeV}{c}#right]/events/0/5",
    "#eta(lead 4^{th} jet)/events/0/5",
    "#phi(lead 4^{th} jet)/events/0/10",
    // (iii) btag monitoring
    "b-discr.(TCHP)/jets/0/2"        ,
    "b-discr.(TCHE)/jets/0/2"	,
    "b-discr.(SSV HEff)/jets/1/2"	,
    "b-discr.(SSV HPur)/jets/0/2"	,
    "b-discr.(CSV)/jets/0/8"	,
    "b-discr.(CSVMVA)/jets/0/8"	,
    "b-discr.(JetBProb)/jets/0/10"	,
    "b-discr.(JetProb)/jets/0/10"	,
    "b-discr.(soft#mu)/jets/0/10"	,
    "b-discr.(soft#muPt)/jets/0/10"  ,                  
    "b-discr.(soft#muIP3d)/jets/0/10",
    "N_{b-jets}(SSVHE)/events/1/1"      ,
    // (iv) MET monitoring 
    "#slash{E}_{T} #left[#frac{GeV}{c}#right]/events/0/10",
    "#SigmaE_{T} [GeV]/events/0/50"  ,
    // (v) Vertices and pileup
    "Number of PU Events/Frequency/1/1",
    "Number of PU Events (Reweighted 1BX)/Frequency/1/1",
    "Number of PU Events (Reweighted 3BX)/Frequency/1/1",
    "Number of Vertices/Frequency/1/1",
    "Number of Vertices (Reweighted 1BX)/Frequency/1/1",
    "Number of Vertices (Reweighted 3BX)/Frequency/1/1",
    // (III) after btagging 
    // (ii) jet monitoring
    "N_{jets}/events/1/1",
    "p_{t}(jets) #left[#frac{GeV}{c}#right]/jets/1/2",
    "#eta(jets)/jets/0/5" ,
    "#phi(jets)/jets/0/10",
    "H_{T} [#frac{GeV}{c}]/events/0/100",
    "p_{t}(lead 1^{st} jet) #left[#frac{GeV}{c}#right]/events/1/5",
    "#eta(lead 1^{st} jet)/events/0/5" ,
    "p_{t}(lead 2^{nd} jet) #left[#frac{GeV}{c}#right]/events/1/5",
    "#eta(lead 2^{nd} jet)/events/0/5" ,
    "p_{t}(lead 3^{rd} jet) #left[#frac{GeV}{c}#right]/events/1/5",
    "#eta(lead 3^{rd} jet)/events/0/5" ,
    "p_{t}(lead 4^{th} jet) #left[#frac{GeV}{c}#right]/events/1/5",
    "#eta(lead 4^{th} jet)/events/0/5" ,
    "p_{t}(lead 1^{st} b-tagged jet) #left[#frac{GeV}{c}#right]/events/1/5",
    "p_{t}(lead 2^{nd} b-tagged jet) #left[#frac{GeV}{c}#right]/events/1/5",
    "#eta(lead 1^{st} b-tagged jet)/events/0/5" ,
    "#eta(lead 2^{nd} b-tagged jet)/events/0/5" ,
    // (iv) MET monitoring 
    "#slash{E}_{T} #left[#frac{GeV}{c}#right]/events/0/20",
    "#SigmaE_{T} [GeV]/events/0/30",
    // (v) Vertices and pileup
    "Number of PU Events/Frequency/0/1",
    "Number of PU Events (Reweighted 1BX)/Frequency/0/1",
    "Number of PU Events (Reweighted 3BX)/Frequency/0/1",
    "Number of Vertices/Frequency/0/1",
    "Number of Vertices (Reweighted 1BX)/Frequency/0/1",
    "Number of Vertices (Reweighted 3BX)/Frequency/0/1", 
    // (III) after kinematic fit 
    "p_{t}^{t and #bar{t}} #left[#frac{GeV}{c}#right]/Frequency/0/20",
    "y^{t and #bar{t}}/Frequency/0/1",
    "p_{t}^{t#bar{t}} #left[#frac{GeV}{c}#right]/Frequency/0/20",
    "y^{t#bar{t}}/Frequency/0/1",
    "m_{t#bar{t}} #left[#frac{GeV}{c^{2}}#right]/Frequency/0/20",
    "p_{t}^{l} #left[#frac{GeV}{c}#right]/Frequency/0/20",    
    "#eta^{l}/Frequency/0/1"
  };

  TString axisLabel1De[ ] = {
    // (iv) electron monitoring
    "N_{e}/events/0/1" ,
    "E(e) [GeV]/events/0/2",
    "E_{t}(e) [GeV]/events/0/1" ,
    "#eta(e)/events/0/5",
    "#phi(e)/events/0/5",
    "#eta(S.C.)/events/0/1"  ,
    "d_{xy} (e wrt. beamspot) [cm]/events/0/1",
    "simpleEleId70cIso/events/0/1", 
    "nHitsInner(conv)/events/0/1",
    "convCot/events/0/5",
    "convDist/events/0/5" ,
    "PF relIso(e)/events/0/1" ,
    // (iv) electron monitoring
    "PF relIso(e)/events/0/1" ,
    "E_{t}(e) [GeV]/events/0/2",
    "#eta(e)/events/0/1",
    "#phi(e)/events/0/1"
  };
	
  TString axisLabel1Dmu[ ] = {
    // (I) preselection
    // (i) muon monitoring
    //"N_{hits}(inner tracker #mu)/events/0/1"          ,
    //"#chi^{2} (global trackfit #mu(pt,#eta))/events/1/1",
    //"d_{xy} (#mu(pt,#eta) wrt. beamspot)/events/0/1" ,
    //"d_{z} (#mu(pt,#eta))/events/0/10",
    //"N_{matched #mu segments}(#mu(pt,#eta))/events/0/1",
    //"#DeltaR(jet(pt,#eta,ID), #mu(pt, #eta, track criteria))/events/0/1",
    //"relIso(#mu(no isolation))/events/0/1",
    //"N_{#mu}/events/0/1",
    // (II) before btagging
    // (i) muon monitoring
    "N_{#mu}/events/0/1"   ,
    "E(#mu)/events/0/2"    ,
    "p_{t}(#mu)/events/0/1",
    "#eta(#mu)/events/0/5" ,
    "y(#mu)/events/0/5"    ,	
    "#phi(#mu)/events/0/5" ,
    "N_{hits}(inner tracker #mu)/events/0/1"          ,
    "#chi^{2} (global trackfit #mu)/events/1/1",
    "d_{xy} (#mu wrt. beamspot) [cm]/events/0/1" ,
    "d_{z} (#mu) [cm]/events/0/10"               ,
    "E_{Ecal} (#mu)/events/1/1",
    "E_{Hcal} (#mu)/events/1/1",
    "PF relIso(#mu)/events/0/1",
    "N_{matched #mu segments}(#mu)/events/0/1",
    // (III) after btagging 
    // (i) muon monitoring
    "PF relIso(#mu)/events/0/1",
    "p_{t}(#mu)/events/0/2",
    "#eta(#mu)/events/0/10",
    "#phi(#mu)/events/0/10",
  };
	
  // 2D: "x-axis title"/"y-axis title"
  TString axisLabel2D[ ] = {
  };
  // count # plots
  unsigned int N1Dplots = sizeof(plots1D)/sizeof(TString);
  if(decayChannel=="electron") N1Dplots+=(sizeof(plots1De )/sizeof(TString));
  if(decayChannel=="muon"    ) N1Dplots+=(sizeof(plots1Dmu)/sizeof(TString));
  unsigned int N2Dplots = sizeof(plots2D)/sizeof(TString);
  // check if all axis labels exist
  unsigned int Naxislabels=sizeof(axisLabel1D)/sizeof(TString);
  if(decayChannel=="electron") Naxislabels+=(sizeof(axisLabel1De )/sizeof(TString));
  if(decayChannel=="muon"    ) Naxislabels+=(sizeof(axisLabel1Dmu)/sizeof(TString));
  if(N1Dplots != Naxislabels) std::cout << "ERROR: some 1D plots or axis label are missing" << std::endl;
  if(N2Dplots != sizeof(axisLabel2D)/sizeof(TString)) std::cout << "ERROR: some 2D plots or axis label are missing" << std::endl;
  if((N1Dplots != Naxislabels)||(N2Dplots != sizeof(axisLabel2D)/sizeof(TString))) exit (1);
  // run automatically in batch mode if there are many canvas
  if((N1Dplots+N2Dplots)>15) gROOT->SetBatch();
	
  // ---
  //    open our standard analysis files
  // ---
  std::map<unsigned int, TFile*> files_           = getStdTopAnalysisFiles(inputFolder, systematicVariation, dataFile, decayChannel);
  std::map<unsigned int, TFile*> filesUncJESUp_   = getStdTopAnalysisFiles(inputFolder, sysJESUp,   dataFile, decayChannel);
  std::map<unsigned int, TFile*> filesUncJESDown_ = getStdTopAnalysisFiles(inputFolder, sysJESDown, dataFile, decayChannel);
	
  // ---
  //    loading histos
  // ---
  // collect all plot names in vector (first 1D, then 2D)
  std::vector<TString> plotList_;
  plotList_.insert( plotList_.begin(), plots1D, plots1D + sizeof(plots1D)/sizeof(TString) );
  if(decayChannel=="electron") plotList_.insert( plotList_.end(), plots1De , plots1De  + sizeof(plots1De )/sizeof(TString) );
  if(decayChannel=="muon"    ) plotList_.insert( plotList_.end(), plots1Dmu, plots1Dmu + sizeof(plots1Dmu)/sizeof(TString) );
  plotList_.insert( plotList_.end(), plots2D, plots2D + sizeof(plots2D)/sizeof(TString) );
	
	
  // container for all histos (1D&2D)
  // example for acess: histo_["plotName"][sampleNr]
  std::map< TString, std::map <unsigned int, TH1F*> > histo_;
  std::map< TString, std::map <unsigned int, TH1F*> > histoUncJESUp_;
  std::map< TString, std::map <unsigned int, TH1F*> > histoUncJESDown_;
  std::map< TString, TH1F* > histoErrorBand_;
  std::map< TString, std::map <unsigned int, TH2F*> > histo2_;
  std::map< TString, std::map <unsigned int, TH2F*> > histo2UncJESUp_;
  std::map< TString, std::map <unsigned int, TH2F*> > histo2UncJESDown_;
  std::map< TString, TH2F* > histo2ErrorBand_;
  // total # plots 
  int Nplots=0;
  // save all histos from plotList_ that exist in files_ into 
  // histo_ and histo2_ and count total # of plots as Nplots
  if(verbose>0) std::cout << std::endl;
  std::vector<TString> vecRedundantPartOfNameInData;
  vecRedundantPartOfNameInData.push_back("_reweighted3BX");
  vecRedundantPartOfNameInData.push_back("_reweighted");
  getAllPlots( files_,           plotList_, histo_,           histo2_,           N1Dplots, Nplots, verbose, decayChannel, &vecRedundantPartOfNameInData);
  getAllPlots( filesUncJESUp_,   plotList_, histoUncJESUp_,   histo2UncJESUp_,   N1Dplots, Nplots, verbose, decayChannel, &vecRedundantPartOfNameInData);
  getAllPlots( filesUncJESDown_, plotList_, histoUncJESDown_, histo2UncJESDown_, N1Dplots, Nplots, verbose, decayChannel, &vecRedundantPartOfNameInData);
	
  // ---
  //    lumiweighting for choosen luminosity
  // ---
  // scale every histo in histo_ and histo2_ to the corresponding luminosity
  // Additionally the mu eff SF is applied
  // NOTE: luminosity [/pb]
  scaleByLuminosity(plotList_, histo_,           histo2_,           N1Dplots, luminosity, verbose, systematicVariation, decayChannel);
  scaleByLuminosity(plotList_, histoUncJESUp_,   histo2UncJESUp_,   N1Dplots, luminosity, verbose, sysJESUp,            decayChannel);
  scaleByLuminosity(plotList_, histoUncJESDown_, histo2UncJESDown_, N1Dplots, luminosity, verbose, sysJESDown,          decayChannel);
	
  // ---
  //    add single top channels and DiBoson contributions
  // ---
  // for SingleTop and DiBoson samples:
  // every plot plotList_ existing in sTop/ diBoson .root file
  // will be combined and saved in the histo_ and histo2_ map
  // reCreate: reCreate combined plots if they are already existing
  bool reCreate=false;
  AddSingleTopAndDiBoson(plotList_, histo_,           histo2_,           N1Dplots, verbose, reCreate, decayChannel);
  AddSingleTopAndDiBoson(plotList_, histoUncJESUp_,   histo2UncJESUp_,   N1Dplots, verbose, reCreate, decayChannel);
  AddSingleTopAndDiBoson(plotList_, histoUncJESDown_, histo2UncJESDown_, N1Dplots, verbose, reCreate, decayChannel);

  // ---
  //    configure histograms
  // ---
  // needs: plotList_, histo_, histo2_, N1Dplots, axisLabel_, axisLabel1D, axisLabel2D
  std::vector<TString> axisLabel_;
  axisLabel_.insert( axisLabel_.begin(), axisLabel1D, axisLabel1D + sizeof(axisLabel1D)/sizeof(TString) );
  if(decayChannel=="electron") axisLabel_.insert( axisLabel_.end(), axisLabel1De , axisLabel1De +sizeof(axisLabel1De )/sizeof(TString) );
  if(decayChannel=="muon"    ) axisLabel_.insert( axisLabel_.end(), axisLabel1Dmu, axisLabel1Dmu+sizeof(axisLabel1Dmu)/sizeof(TString) );
  axisLabel_.insert( axisLabel_.end()  , axisLabel2D, axisLabel2D + sizeof(axisLabel2D)/sizeof(TString) );
  if(verbose>1){
    std::cout << "(plot, x Axis label , y Axis label , log scale?, rebinning factor):" << std::endl;
    // loop plots
    for(unsigned int plot=0; plot<plotList_.size(); ++plot){
      std::cout << plotList_[plot] << ": " << getStringEntry(axisLabel_[plot],1);
      std::cout << " , " << getStringEntry(axisLabel_[plot], 2);
      if(plot<N1Dplots){
	std::cout<< " , " << getStringEntry(axisLabel_[plot], 3);
	std::cout<< " , " << getStringEntry(axisLabel_[plot], 4);
      }
      std::cout << std::endl;
    }
  }

  // loop samples
  for(unsigned int sample=kSig; sample<=kData; ++sample){
    // loop plots
    for(unsigned int plot=0; plot<plotList_.size(); ++plot){
      // a) 1D
      if((plot<N1Dplots)&&(histo_.count(plotList_[plot])>0)&&(histo_[plotList_[plot]].count(sample)>0)){ 
	// default
	histogramStyle( *histo_[plotList_[plot]][sample],           sample, true);
	histogramStyle( *histoUncJESUp_[plotList_[plot]][sample],   sample, true);
	histogramStyle( *histoUncJESDown_[plotList_[plot]][sample], sample, true);
	// special configurations
	if(getStringEntry(plotList_[plot], 2)=="PartonJetDRall"){
	  histo_[plotList_[plot]][sample]           -> SetNdivisions(816);
	  histoUncJESUp_[plotList_[plot]][sample]   -> SetNdivisions(816);
	  histoUncJESDown_[plotList_[plot]][sample] -> SetNdivisions(816);
	}
      }
      // b) 2D
      if((plot>=N1Dplots)&&(histo2_.count(plotList_[plot])>0)&&(histo2_[plotList_[plot]].count(sample)>0)){
	histStyle2D( *histo2_[plotList_[plot]][sample], sampleLabel(sample,decayChannel), getStringEntry(axisLabel_[plot],1), getStringEntry(axisLabel_[plot],2));
 	histStyle2D( *histo2UncJESUp_[plotList_[plot]][sample], sampleLabel(sample,decayChannel), getStringEntry(axisLabel_[plot],1), getStringEntry(axisLabel_[plot],2));
      	histStyle2D( *histo2UncJESDown_[plotList_[plot]][sample], sampleLabel(sample,decayChannel), getStringEntry(axisLabel_[plot],1), getStringEntry(axisLabel_[plot],2));
      }
    }
  }
  
  // ---
  //    event composition
  // ---
  std::map< TString, std::map <unsigned int, double> > events_;
  // a) get event numbers
  std::vector<TString> selection_;
  selection_.push_back("tightJetKinematics/n"      );
  selection_.push_back("tightJetKinematicsTagged/n");
  selection_.push_back("analyzeTopRecoKinematicsKinFit/ttbarMass");  
  unsigned int MCBG=42;
  events_[selection_[0]][MCBG]=0;
  events_[selection_[1]][MCBG]=0; 
  events_[selection_[2]][MCBG]=0;
  // loop pretagged/tagged
  for(unsigned int step=0; step<selection_.size(); ++step){
    // loop samples
    for(unsigned int sample=kSig; sample<=kData; ++sample){
      // save number
      events_[selection_[step]][sample]=histo_[selection_[step]][sample]->Integral(0,histo_[selection_[step]][sample]->GetNbinsX()+1);
      // add non ttbar MC
      if(sample>kSig&&sample<kData) events_[selection_[step]][MCBG]+=events_[selection_[step]][sample];
    }
  }
  // b) print composition
  if(verbose>=0){
    // loop pretagged/tagged
    for(unsigned int step=0; step<selection_.size(); ++step){    
      // print label
      switch (step){
          case 0 : std::cout << std::endl << "Event composition ---- pre-tagged, derived from: "+selection_[step]          << std::endl; break;
          case 1 : std::cout << std::endl << "Event composition ---- after b-tagging, derived from: "+selection_[step]     << std::endl; break;
          case 2 : std::cout << std::endl << "Event composition ---- after kinematic fit, derived from: "+selection_[step] << std::endl; break;
          default: break;
      }
      // all MC events
      double NAllMC=events_[selection_[step]][kSig]+events_[selection_[step]][MCBG];
      // print yield and composition
      std::cout << "events observed in data: " << events_[selection_[step]][kData] << std::endl;
      std::cout << "events expected from MC: " << NAllMC << std::endl;
      std::cout << "expected event composition:"   << std::endl; 
      std::cout << "ttbar SG:   " << std::setprecision(4) << std::fixed << events_[selection_[step]][kSig  ] / NAllMC << std::endl;
      std::cout << "ttbar BG:   " << std::setprecision(4) << std::fixed << events_[selection_[step]][kBkg  ] / NAllMC << std::endl;
      std::cout << "W+jets:     " << std::setprecision(4) << std::fixed << events_[selection_[step]][kWjets] / NAllMC << std::endl; 
      std::cout << "Z+jets:     " << std::setprecision(4) << std::fixed << events_[selection_[step]][kZjets] / NAllMC << std::endl;
      std::cout << "QCD:        " << std::setprecision(4) << std::fixed << events_[selection_[step]][kQCD  ] / NAllMC << std::endl;
      std::cout << "single top: " << std::setprecision(4) << std::fixed << events_[selection_[step]][kSTop ] / NAllMC << std::endl;
      std::cout << "diboson:    " << std::setprecision(4) << std::fixed << events_[selection_[step]][kDiBos] / NAllMC << std::endl;
    }
  }
	
  // ---
  //    rebinning 1D histograms
  // ---
  // loop samples
  for(unsigned int sample=kSig; sample<=kData; ++sample){
    // loop plots
    for(unsigned int plot=0; plot<plotList_.size(); ++plot){
      TString plotName=plotList_[plot];
      // check if plot exists and is 1D
      if((plot<N1Dplots)&&(plotExists(histo_, plotName, sample))){
	// equidistant binning
	double reBinFactor = atof(((string)getStringEntry(axisLabel_[plot],4)).c_str());
	if(reBinFactor>1){
	  equalReBinTH1(reBinFactor, histo_,           plotName, sample);
	  equalReBinTH1(reBinFactor, histoUncJESUp_,   plotName, sample);
	  equalReBinTH1(reBinFactor, histoUncJESDown_, plotName, sample);
	}
      }
    }
  }
  // ========================================================
  //  Errors for uncertainty bands from JES and luminosity
  // ========================================================
  
  std::cout << std::endl << " Start calculating error bands .... " << std::endl;

  for(unsigned int plot=0; plot<plotList_.size(); ++plot)
  {
    TString plotName     = plotList_[plot];

    // Initialize and reset histograms

    TH1F* histoSumRef  = (TH1F*)histo_[plotName][kSig]->Clone();
    TH1F* histoSumUp   = (TH1F*)histo_[plotName][kSig]->Clone();
    TH1F* histoSumDown = (TH1F*)histo_[plotName][kSig]->Clone();
    TH1F* histoError   = (TH1F*)histo_[plotName][kSig]->Clone();

    histoSumRef    -> Reset("ICESM");
    histoSumUp     -> Reset("ICESM");
    histoSumDown   -> Reset("ICESM");
    histoError     -> Reset("ICESM");

    // Integral over all samples before accessing the differences

    for(unsigned int sample=kSig; sample<kData; ++sample)
    {
      if((plot<N1Dplots)&&(plotExists(histo_, plotName, sample)))
      {
	histoSumRef  -> Add(histo_[plotName][sample]);
	histoSumUp   -> Add(histoUncJESUp_[plotName][sample]);
	histoSumDown -> Add(histoUncJESDown_[plotName][sample]);
      }
    }

    // Compare summed histograms, symmetrize deviations and store relative deviation to error histogram after adding constant contributions
        
    for (int bin = 0; bin < histoSumRef->GetNbinsX(); bin++)
    {
      double diffUp   = fabs(histoSumUp->GetBinContent(bin+1)   - histoSumRef->GetBinContent(bin+1));
      double diffDown = fabs(histoSumDown->GetBinContent(bin+1) - histoSumRef->GetBinContent(bin+1));
      double maxDiff  = (diffUp > diffDown) ? diffUp : diffDown; // symmetrize error

      maxDiff = (histoSumRef->GetBinContent(bin+1) > 0) ? (maxDiff/histoSumRef->GetBinContent(bin+1)) : 0; 
      maxDiff = sqrt(maxDiff*maxDiff+lumiError*lumiError); // add constant error for luminosity (defined at beginning of macro)
      histoError -> SetBinContent(bin+1,maxDiff);
      
      histoSumRef->SetBinError(bin+1,maxDiff*histoSumRef->GetBinContent(bin+1));
    }

    histoErrorBand_[plotName] = (TH1F*)histoSumRef->Clone();
  }

  std::cout << " .... Finished calculation of error bands." << std::endl;  

  // ---
  //    create legends
  // ---
  unsigned int Nlegends=0; 
  TLegend *leg  = new TLegend(); 
  TLegend *leg0 = new TLegend(0.05, 0.15, 1.05, 0.9);
  TLegend *leg1 = new TLegend(0.05, 0.15, 1.05, 0.9);
  Nlegends+=3;
  leg ->SetFillStyle(0);
  leg ->SetBorderSize(0);
  leg ->SetTextSize(0.03);
  leg ->SetTextAlign(12);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  leg0->SetHeader("After Selection, Before b-Tagging");
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetHeader("After Selection & b-Tagging");
  // fill in contributing sample
  // data is to be first entry
  bool TwoThousandElevenData=false;
  TString lumilabel = Form("%3.2f fb^{-1}",luminosity/1000);
  if(luminosity>36.0) TwoThousandElevenData=true;
  else lumilabel=Form("%2.0f pb^{-1}",luminosity);    
  leg ->AddEntry(histo_[plotList_[plotList_.size()-1]][kData], sampleLabel(kData,decayChannel,TwoThousandElevenData),"P");
  leg0->AddEntry(histo_[plotList_[plotList_.size()-1]][kData], sampleLabel(kData,decayChannel,TwoThousandElevenData)+", "+lumilabel,"PL");
  leg1->AddEntry(histo_[plotList_[plotList_.size()-1]][kData], sampleLabel(kData,decayChannel,TwoThousandElevenData)+", "+lumilabel,"PL");
  // now loop over MC samples
  for(unsigned int sample=kSig; sample<kData; ++sample){
    // check if sample exists in at least one plot
    bool exit=false;
    // loop plots
    for(unsigned int plot=0; plot<plotList_.size()-1; ++plot){  // <plotList_.size()-1, because last entry is for data (see above)
      // if found: add entry to legend
      if((histo_.count(plotList_[plot])>0)&&(histo_[plotList_[plot]].count(sample)>0)&&(!exit)){
	leg ->AddEntry(histo_[plotList_[plot]][sample], sampleLabel(sample,decayChannel), "F");
	leg0->AddEntry(histo_[plotList_[plot]][sample], sampleLabel(sample,decayChannel), "F");
	leg1->AddEntry(histo_[plotList_[plot]][sample], sampleLabel(sample,decayChannel), "F");
	exit=true;
      }
    }
  }
  // add entry for uncertainty to legends
  if (histoErrorBand_.size() > 0 && plotList_.size() > 0){
      leg ->AddEntry(histoErrorBand_[plotList_[0]],"Uncertainty","F");
      leg0->AddEntry(histoErrorBand_[plotList_[0]],"Uncertainty","F");
      leg1->AddEntry(histoErrorBand_[plotList_[0]],"Uncertainty","F");
  }
  // ---
  //    create canvas
  // ---
  std::vector<TCanvas*> plotCanvas_;
  //  int NCanvas = Nplots+Nlegends;
  // a) create canvas for all plots + legends
  //  for(int sample=0; sample<NCanvas; sample++){
  for(unsigned int sample=0; sample<N1Dplots+N2Dplots+Nlegends; sample++){
    char canvname[10];
    sprintf(canvname,"canv%i",sample);    
    plotCanvas_.push_back( new TCanvas( canvname, canvname, 600, 600) );
    //canvasStyle(*plotCanvas_[sample]);
  }
	
  // ---
  //    change 1D plots into stack plots
  // ---
  // loop plots -> all 1D plots will become stacked plots
  if(verbose>1) std::cout << std::endl;
  for(unsigned int plot=0; plot<plotList_.size(); ++plot){
    createStackPlot(plotList_, histo_, plot, N1Dplots, verbose, decayChannel);
  }
  if(verbose>1) std::cout << std::endl;
	
  // ---
  //    do the printing
  // ---
  // a) for plots
  int canvasNumber=0;
  // loop plots
  for(unsigned int plot=0; plot<plotList_.size(); ++plot){
    bool first=true;
    // open canvas and set titel corresponding to plotname in .root file
    plotCanvas_[canvasNumber]->cd(0);
    plotCanvas_[canvasNumber]->SetTitle(getStringEntry(plotList_[plot], 2)+getStringEntry(plotList_[plot], 1));	    
    // loop samples
    for(unsigned int sample=kSig; sample<=kData; ++sample){
      // a1) for 1D event yields, efficiency and cross section plots (existing)
      if((plot<N1Dplots)||(plot>=N1Dplots+N2Dplots)){
	// check if plot is existing
	if((histo_.count(plotList_[plot])>0)&&(histo_[plotList_[plot]].count(sample)>0)){
	  histo_[plotList_[plot]][sample]->GetXaxis()->SetNoExponent(true);
	  if(verbose>0){
	    std::cout << "plotting " << plotList_[plot];
	    std::cout << " from sample " << sampleLabel(sample,decayChannel);
	    std::cout << " to canvas " << canvasNumber << " ( ";
	    std::cout << plotCanvas_[canvasNumber]->GetTitle() << " )" << std::endl;
	  }
	  // first plot
	  if(first){
	    // plot to take maximum from
	    TString maxFrom=plotList_[plot];
	    if(plotList_[plot].Contains("_JetKinematicsTagged/pt")) maxFrom="tightLead_3_JetKinematicsTagged/pt";
	    if(plotList_[plot].Contains("_JetKinematics/pt")) maxFrom="tightLead_3_JetKinematics/pt";
	    if(!histo_[maxFrom].count(sample)>0||!histo_[maxFrom].count(kData)>0) maxFrom=plotList_[plot];
	    if(verbose>1) std:: cout << "take max for histo " << plotList_[plot] << " from histo " << maxFrom << std::endl;
	    // min / max
	    double max = 1.3*histo_[maxFrom][sample]->GetMaximum();
	    // if data file exists
	    if(histo_[maxFrom].count(kData)>0){
	      // and has larger maximum
	      if(max < 1.3*histo_[maxFrom][kData]->GetMaximum()){
		// take this maximum
		max = 1.3*histo_[maxFrom][kData]->GetMaximum();
		if(verbose>1) std:: cout << "take max from data! " << std::endl;
	      }
	    }
	    double min = 0;
	    // log plots
	    if(getStringEntry(axisLabel_[plot],3)=="1"){
	      plotCanvas_[canvasNumber]->SetLogy(1);
	      min=1;
	      max=exp(1.3*(std::log(max)-std::log(min))+std::log(min));
	      if(plotList_[plot]=="analyzeTopRecoKinematicsKinFit/prob") min=0.1; 
	      if(plotList_[plot].Contains("_JetKinematicsTagged/pt")) min=0.1;
	    }
	    // get nicer int values if maximum is large enough
	    if(max>3) max = (double)roundToInt(max);
	    if(plotList_[plot].Contains("btagSimpleSecVtx"))max*=0.8;
	    // restrict x axis for special plots
	    if(getStringEntry(plotList_[plot], 2)=="nHit"  ) histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(10,30 );
	    if(getStringEntry(plotList_[plot], 2)=="chi2"  ) histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(0,10  );
	    if(getStringEntry(plotList_[plot], 2)=="dB"    ) histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(0,0.02);
	    if(plotList_[plot].Contains("relIso")) histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(0,0.15);
	    if(plotList_[plot].Contains("tightJetKinematics/n")||plotList_[plot].Contains("tightJetKinematicsTagged/n")) histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(4,10);
	    if(plotList_[plot].Contains("tightJetKinematicsNjets1/n")) histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(1,10);
	    if(plotList_[plot].Contains("_JetKinematics/en")) histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(0,500);
	    if(plotList_[plot].Contains("_JetKinematics/pt")) histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(0,300);
	    if(plotList_[plot].Contains("analyzeMETMuon/metSumEt")) histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(200,1400);
	    if(plotList_[plot].Contains("analyzeMETMuon/metEt")) histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(0,300);
	    if(plotList_[plot].Contains("bottomJetKinematics/n")) histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(0,5);
	    if(plotList_[plot].Contains("JetKinematicsTagged")&&plotList_[plot].Contains("pt")) histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(0,350);
	    if(plotList_[plot].Contains("btagSimpleSecVtx")) histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(-1,7);			
	    // axis style
	    axesStyle(*histo_[plotList_[plot]][sample], getStringEntry(axisLabel_[plot],1), getStringEntry(axisLabel_[plot],2), min, max);
	    histo_[plotList_[plot]][sample]->GetXaxis()->SetNoExponent(true);
	    if(max<100) histo_[plotList_[plot]][sample]->GetYaxis()->SetNoExponent(true);
	    else histo_[plotList_[plot]][sample]->GetYaxis()->SetNoExponent(false);
	    // draw histos (as stack)
	    histo_[plotList_[plot]][sample]->Draw("hist X0");
	    histo_[plotList_[plot]][42] = (TH1F*)(histo_[plotList_[plot]][sample]->Clone());
	    histo_[plotList_[plot]][42]->GetXaxis()->SetNoExponent(true);
	  }
	  // draw other plots into same canvas
	  else{ 
	    // draw data as points
	    if(sample==kData) histo_[plotList_[plot]][sample]->Draw("p e X0 same");
	    // draw others as histo (stack)
	    else histo_[plotList_[plot]][sample]->Draw("hist X0 same");
	  }
	  first=false;
	  // at the end:
	  if((histo_.count(plotList_[plot])>0)&&(sample==kData)){
	    // configure style of and draw uncertainty bands
	    histoErrorBand_[plotList_[plot]]->SetMarkerStyle(0);
	    histoErrorBand_[plotList_[plot]]->SetFillColor(1);
	    histoErrorBand_[plotList_[plot]]->SetFillStyle(3004);
	    gStyle->SetErrorX(0.5);  
	    histoErrorBand_[plotList_[plot]]->Draw("E2 SAME");	 
	    // redraw axis 
	    histo_[plotList_[plot]][42]->Draw("axis X0 same");
	    if((unsigned int)canvasNumber<plotCanvas_.size()-Nlegends){
	      // draw label indicating event selection, common labels and legend
	      TString label = "pre-Tagged";
	      if(plotList_[plot].Contains("Tagged")) label = "Tagged";
	      if(plotList_[plot].Contains("PreSel")) label = "Pre-Selected";
	      if(plotList_[plot].Contains("Njets1")) label = "#geq 1 Jet";
	      DrawLabel(label, 1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.2, 1.0 - gStyle->GetPadTopMargin() - gStyle->GetTickLength() - 0.05,
		 	       1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength(),       1.0 - gStyle->GetPadTopMargin() - gStyle->GetTickLength(), 12    );
	      leg->SetX1NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.20);
	      leg->SetY1NDC(1.0 - gStyle->GetPadTopMargin()   - gStyle->GetTickLength() - 0.32);
	      leg->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength());
	      leg->SetY2NDC(1.0 - gStyle->GetPadTopMargin()   - gStyle->GetTickLength() - 0.05);
	      leg->Draw("SAME");
	      // add labels for decay channel, luminosity, energy and CMS preliminary (if applicable)
	      if (decayChannel=="muon") DrawDecayChLabel("#mu + Jets");
	      else if (decayChannel=="electron") DrawDecayChLabel("e + Jets");
	      else DrawDecayChLabel("e/#mu + Jets Combined");
	      // set first parameter to false once "CMS Preliminary" is not required anymore
	      DrawCMSLabels(true,luminosity); 
	      //draw data/MC ratio
	      if((histo_[plotList_[plot]].count(kSig)>0) && withRatioPlot){
		drawRatio(histo_[plotList_[plot]][kData], histo_[plotList_[plot]][kSig], 0.1, 1.9, myStyle, verbose);	       
	      }
	    }
	  }
	}
      }
      // a3) for 2D plots
      if((plot>=N1Dplots)&&(histo2_.count(plotList_[plot])>0)&&(histo2_[plotList_[plot]].count(sample)>0)){
	// new Canvas for every plot
	plotCanvas_[canvasNumber]->cd(0);
	plotCanvas_[canvasNumber]->SetRightMargin ( 0.15 );
	plotCanvas_[canvasNumber]->SetTitle(getStringEntry(plotList_[plot],2)+getStringEntry(plotList_[plot],1)+getTStringFromInt(sample)); 
	if(verbose>1){
	  std::cout << "plotting " << plotList_[plot];
	  std::cout << " from sample " << sampleLabel(sample,decayChannel);
	  std::cout << " to canvas " << canvasNumber  << " ( ";
	  std::cout << plotCanvas_[canvasNumber]->GetTitle() << " )"  << std::endl;
	}
	++canvasNumber;
	// draw histo
	histo2_[plotList_[plot]][sample]->Draw("colz");
	// print correlation factor
	double d = histo2_[plotList_[plot]][sample]->GetCorrelationFactor();
	char correlation[20];
	sprintf(correlation, "%f", d);
	TString corr = (TString)correlation;
	DrawLabel("correlation: "+corr, 0.35, 0.92, 0.75, 0.99, 0.7);
      }
    }
    // for 1D hists: next canvas
    if((plot<N1Dplots)&&(histo_.count(plotList_[plot])>0)) ++canvasNumber;
  }
  // b) for legends
  // b1) before btagging
  plotCanvas_[canvasNumber]->cd(0);
  plotCanvas_[canvasNumber]->SetTitle("legendMonitoringBeforeBtagging");
  leg0->Draw("");
  ++canvasNumber;
  // b2) after btagging
  plotCanvas_[canvasNumber]->cd(0);
  plotCanvas_[canvasNumber]->SetTitle("legendMonitoringAfterBtagging");
  leg1->Draw("");
  ++canvasNumber;
	
  // ---
  // saving
  // ---
  if(save){
    // pdf, eps and png
    if(verbose==0) gErrorIgnoreLevel=kWarning;
    saveCanvas(plotCanvas_, outputFolder, pdfName, true, false);
    for(unsigned int idx=0; idx<plotCanvas_.size(); idx++){
      TString title=(TString)(plotCanvas_[idx])->GetTitle();
      if(!title.Contains("canv")){
	plotCanvas_[idx]->Print(outputFolder+title+".eps"); 
	plotCanvas_[idx]->Print(outputFolder+title+".png");
      }
    }
    // root file
    if(verbose>0) std::cout << "will save into outputfile named " << outputFileName << std::endl;
    for(unsigned int idx=0; idx<plotCanvas_.size(); ++idx){
      saveToRootFile(outputFileName, plotCanvas_[idx], true, verbose, "monitoring");
    }
  }
}
