#include "basicFunctions.h"

void analyzeHypothesisKinFit(double luminosity = 1143.22, bool save = true, int systematicVariation=sysNo, unsigned int verbose=1, TString inputFolderName="TOP2011/110819_AnalysisRun",
			     //TString dataFile= "/afs/naf.desy.de/group/cms/scratch/tophh/TOP2011/110819_AnalysisRun/analyzeDiffXData2011A_Elec_160404_167913_1fb.root",
			     TString dataFile= "/afs/naf.desy.de/group/cms/scratch/tophh/TOP2011/110819_AnalysisRun/analyzeDiffXData2011A_Muon_160404_167913_1fb.root",
			     std::string decayChannel = "muon" )
{
  // c) set root style
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
  // systematicVariation= "JESup", "JESdown", "JERup", "JERdown", "PileUp", "ScaleUp", 
  //                      "ScaleDown", "MatchUp", "MatchDown"
  
  //  ---
  //     options
  //  ---
  // a) options directly entered when calling function
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
  if(luminosity<40.&&systematicVariation==sysLumiUp  ) luminosity*=1.04;
  else if(luminosity<40.&&systematicVariation==sysLumiDown) luminosity*=0.96;
  if(luminosity>40.&&systematicVariation==sysLumiUp  ) luminosity*=1.045;
  else if(luminosity>40.&&systematicVariation==sysLumiDown) luminosity*=0.955;
  // verbose: set detail level of output 
  // 0: no output, 1: std output 2: more output 3: output for debugging
  // data file: relative path of .root file
  // save: save plots?
  // luminosity: [/pb]
  TString lumi = getTStringFromInt(roundToInt((luminosity), false));
  // b) options to be configured only once
  // get the .root files from the following folder:
  TString inputFolder = "/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName;
  // see if its 2010 or 2011 data from luminosity
  TString dataSample="";
  if(luminosity<50) dataSample="2010";
  if(luminosity>50) dataSample="2011";
  // save all plots into the following folder
  TString outputFolder = "./diffXSecFromSignal/plots/"+decayChannel+"/";
  if(dataSample!="") outputFolder+=dataSample+"/";
  // save all plots within a root file named:
  TString outputFileName="diffXSecTopSemi";
  if(decayChannel=="muon"    ) outputFileName+="Mu";
  if(decayChannel=="electron") outputFileName+="Elec";
  outputFileName+=dataSample+".root";
  // choose name of the output .pdf file
  TString pdfName="kinFitpbHypothesis"+lumi+"pb";
  // choose if differential xSecs are extrapolated to whole phase space
  bool extrapolate=false;
  TString PS="";
  if(!extrapolate)PS="PhaseSpace";
  // choose plot input folder corresponding to systematic Variation  
  TString sysInputFolderExtension="";
  TString sysInputGenFolderExtension="";
  if(systematicVariation==sysPUup  ){ 
    sysInputFolderExtension="PUup";
    sysInputGenFolderExtension=sysInputFolderExtension;
  }
  else if(systematicVariation==sysPUdown){
    sysInputFolderExtension="PUdown";
    sysInputGenFolderExtension=sysInputFolderExtension;
  }
  else if(systematicVariation==sysflatTrigSF       ) sysInputFolderExtension="FlatTriggerSF";
  else if(systematicVariation==sysTrigEffSFNormUp  ) sysInputFolderExtension="TriggerEffSFNormUp";
  else if(systematicVariation==sysTrigEffSFNormDown) sysInputFolderExtension="TriggerEffSFNormDown";
  else if(systematicVariation==sysTriggerEffSFShapeUpEta  ) sysInputFolderExtension="TriggerEffSFShapeUpEta";
  else if(systematicVariation==sysTriggerEffSFShapeDownEta) sysInputFolderExtension="TriggerEffSFShapeDownEta";
  else if(systematicVariation==sysTriggerEffSFShapeUpPt   ) sysInputFolderExtension="TriggerEffSFShapeUpPt";
  else if(systematicVariation==sysTriggerEffSFShapeDownPt ) sysInputFolderExtension="TriggerEffSFShapeDownPt";
  else if(systematicVariation==sysMuEffSFup   ) sysInputFolderExtension="SelectionEffSFNormUp";
  else if(systematicVariation==sysMuEffSFdown ) sysInputFolderExtension="SelectionEffSFNormDown";
  else if(systematicVariation==sysBtagSFUp    ) sysInputFolderExtension="BtagSFup";
  else if(systematicVariation==sysBtagSFDown  ) sysInputFolderExtension="BtagSFdown";
  else if(systematicVariation==sysMisTagSFup  ) sysInputFolderExtension="MisTagSFup";
  else if(systematicVariation==sysMisTagSFdown) sysInputFolderExtension="MisTagSFdown";
  // additional (control plot folders): NoWeight, OnlyPUWeight, NoBtagSFWeight;

  //  ---
  //     choose plots
  //  ---
  // a) list plots you would like to see ("folder/plotName") - same as in .root files (for 1D and 2D)
  TString plots1D[ ] = { // general fit performance
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/prob"       , 
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/chi2"       , 
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/delChi2"    ,
    // combinatorics and Kinfit Hypothesis Quality(ttbar signal only)
    "analyzeHypoKinFit/hadBQuark"               , 
    "analyzeHypoKinFit/lepBQuark"               , 
    "analyzeHypoKinFit/lightQuark"              , 
    "analyzeHypoKinFit/wrongAssign"             ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/qAssignment",
    "analyzeHypoKinFit/PartonJetDRall"          ,
    // pull distributions
    "analyzeHypoKinFit/hadBQuarkPt"             ,
    "analyzeHypoKinFit/hadBQuarkEta"            ,
    "analyzeHypoKinFit/hadBQuarkPhi"            ,
    "analyzeHypoKinFit/lepBQuarkPt"             ,
    "analyzeHypoKinFit/lepBQuarkEta"            ,
    "analyzeHypoKinFit/lepBQuarkPhi"            ,
    "analyzeHypoKinFit/lightQuarkPt"            ,
    "analyzeHypoKinFit/lightQuarkEta"           ,
    "analyzeHypoKinFit/lightQuarkPhi"           ,
    "analyzeHypoKinFit/leptonPt"                ,
    "analyzeHypoKinFit/leptonEta"               ,
    "analyzeHypoKinFit/leptonPhi"               ,
    "analyzeHypoKinFit/neutrinoPt"              ,
    "analyzeHypoKinFit/neutrinoEta"             ,
    "analyzeHypoKinFit/neutrinoPhi"             ,
    "analyzeHypoKinFit/lepWPt"                  ,
    "analyzeHypoKinFit/lepWEta"                 ,
    "analyzeHypoKinFit/lepWPhi"                 ,
    "analyzeHypoKinFit/hadWPt"                  ,
    "analyzeHypoKinFit/hadWEta"                 ,
    "analyzeHypoKinFit/hadWPhi"                 ,
    // reconstructed top quantities
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topMass"    ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topPt"      ,                         
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topPhi"     ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topY"       ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topPtHad"   ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topPhiHad"  ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topYHad"    ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topPtLep"   ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topPhiLep"  ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topYLep"    ,
    // reconstructed angular distributions
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/bbbarAngle"  ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/bbbarAngleTtRF",
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/WWAngle"     ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topWAngleLep",
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topWAngleHad",
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topBAngleLep",
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topBAngleHad",
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/bWAngleLep"  ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/bWAngleHad"  ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/qqbarAngle"  ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/qBlepAngle"  ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/qBhadAngle"  ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/lepBlepAngle",
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/lepBlepAngleTtRF",
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/lepBhadAngle",
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/lepQAngle"   ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/MuonNeutrinoAngle",
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/lepBNeutrinoAngle",
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/hadBNeutrinoAngle",
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/qNeutrinoAngle"   ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/lepWDir"     ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/nuWDir"      ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/qWDir"       ,
    // reconstructed event shape variables
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/aplanarity" ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/sphericity" ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/C"          ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/D"          ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/circularity",
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/isotropy"   ,
    // generated top quantities
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/topMass"      , 
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/topPt"        ,    
    "analyzeTopPartonLevelKinematicsPhaseSpace"+sysInputGenFolderExtension+"/topPt", 
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/topPhi"       ,
    "analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/topY",
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/topPtHad"     ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/topPhiHad"    ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/topYHad"      ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/topPtLep"     ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/topPhiLep"    ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/topYLep"      ,
    // generated angular distributions
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/bbbarAngle"   ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/bbbarAngleTtRF",
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/WWAngle"     ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/topWAngleLep",
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/topWAngleHad",
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/topBAngleLep",
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/topBAngleHad",
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/bWAngleLep"  ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/bWAngleHad"  ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/qqbarAngle"  ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/qBlepAngle"  ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/qBhadAngle"  ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/lepBlepAngle",
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/lepBlepAngleTtRF",
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/lepBhadAngle",
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/lepQAngle"   ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/MuonNeutrinoAngle",
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/lepBNeutrinoAngle",
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/hadBNeutrinoAngle",
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/qNeutrinoAngle"   ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/lepWDir"     ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/qWDir"       ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/nuWDir"      ,
    // reconstructed ttbar quantities
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/ttbarMass"  ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/ttbarPt"    ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/ttbarY"     ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/ttbarHT"    ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/ttbarSumY"  ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/ttbarDelPhi",
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/ttbarDelY"  ,
    // generated ttbar quantities
    "analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/ttbarMass",
    "analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/ttbarPt"  ,
    "analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/ttbarY"   ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/ttbarHT"    ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/ttbarSumY"  ,
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/ttbarDelPhi",
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/ttbarDelY"  ,
    // reconstructed lepton quantities
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/lepPt" ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/lepEta",
    // generated lepton quantities
    "analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/lepPt" ,
    "analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/lepEta"
  };
  TString plots2D[ ] = { // reco - gen Match correlation plots (ttbar signal only)
    // a) combinatorics and Kinfit Hypothesis Quality(ttbar signal only)
    "analyzeHypoKinFit/mapKinFit_"               ,
    // b) reconstructed Top quantities
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topPt_"      ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topPhi_"     ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topY_"       ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/bbbarAngle_" ,
    // c) reconstructed ttbar quantities
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/ttbarMass_"  ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/ttbarPt_"    ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/ttbarY_"     ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/ttbarHT_"    ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/ttbarSumY_"  ,
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/ttbarDelPhi_",
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/ttbarDelY_"  
  };

  // b) list plot axes style
  // 1D: "x-axis title"/"y-axis title"/log/rebin-factor
  // log = 0 or 1 for linear or logarithmic axis 

  TString axisLabel1D[ ] = { // general fit performance
    "probability (best fit hypothesis)/events/1/25"                   , 
    "#chi^{2} (best fit hypothesis)/events/0/10"                      ,
    "#Delta#chi^{2} (1^{st} - 2^{nd} best fit hypothesis)/events/0/10",
    // combinatorics and Kinfit Hypothesis Quality(ttbar signal only)
    "#Deltai_{lead jet}(genMatch - kinFit), hadronic b-quark/events/0/1", 
    "#Deltai_{lead jet}(genMatch - kinFit), leptonic b-quark/events/0/1", 
    "#Deltai_{lead jet}(genMatch - kinFit), light quarks/events/0/1"    ,
    "N(wrong assigned jets)/events/0/1"                                 , 
    "permutation/events/0/1"                                            ,
    "#DeltaR(parton, reco jet assigned from Kinfit)/partons/1/10"       ,
    // pull distributions
    "(p_{t gen}-p_{t fit}) (#sigmap_{t})^{-1} (hadronic b-quark)/events/0/1"  ,
    "(#eta_{gen}-#eta_{fit}) (#sigma#eta)^{-1} (hadronic b-quark)/events/0/1" ,
    "(#phi_{gen}-#phi_{fit}) (#sigma#phi)^{-1} (hadronic b-quark)/events/0/1" ,
    "(p_{t gen}-p_{t fit}) (#sigmap_{t})^{-1} (leptonic b-quark)/events/0/1"  ,
    "(#eta_{gen}-#eta_{fit}) (#sigma#eta)^{-1} (leptonic b-quark)/events/0/1" ,
    "(#phi_{gen}-#phi_{fit}) (#sigma#phi)^{-1} (leptonic b-quark)/events/0/1" ,
    "(p_{t gen}-p_{t fit}) (#sigmap_{t})^{-1} (light quark)/events/0/1"       ,
    "(#eta_{gen}-#eta_{fit}) (#sigma#eta)^{-1} (light quark)/events/0/1"      ,
    "(#phi_{gen}-#phi_{fit}) (#sigma#phi)^{-1} (light quark)/events/0/1"      ,
    "(p_{t gen}-p_{t fit}) (#sigmap_{t})^{-1} (lepton)/events/0/1"            ,
    "(#eta_{gen}-#eta_{fit}) (#sigma#eta)^{-1} (lepton)/events/0/1"           ,
    "(#phi_{gen}-#phi_{fit}) (#sigma#phi)^{-1} (lepton)/events/0/1"           ,
    "(p_{t gen}-p_{t fit}) (#sigmap_{t})^{-1} (neutrino)/events/0/1"          ,
    "(#eta_{gen}-#eta_{fit}) (#sigma#eta)^{-1} (neutrino)/events/0/1"         ,
    "(#phi_{gen}-#phi_{fit}) (#sigma#phi)^{-1} (neutrino)/events/0/1"         ,
    "(p_{t fit}-p_{t gen}) p_{t gen}^{-1} (leptonic W)/events/0/1"          ,
    "(#eta_{fit}-#eta_{gen}) #eta_{gen}^{-1} (leptonic W)/events/0/1"       ,
    "(#phi_{fit}-#phi_{gen}) #phi_{gen}^{-1} (leptonic W)/events/0/1"       ,
    "(p_{t fit}-p_{t gen}) p_{t gen}^{-1} (hadronic W)/events/0/1"          ,
    "(#eta_{fit}-#eta_{gen}) #eta_{gen}^{-1} (hadronic W)/events/0/1"       ,
    "(#phi_{fit}-#phi_{gen}) #phi_{gen}^{-1} (hadronic W)/events/0/1"       ,
    // reconstructed top quantities
    "m^{t and #bar{t}} #left[#frac{GeV}{c^{2}}#right]/top quarks #left[(#frac{GeV}{c^{2}})^{-1}#right]/0/10",
    "p_{T}^{t and #bar{t}} #left[#frac{GeV}{c}#right]/top quarks #left[(#frac{GeV}{c})^{-1}#right]/0/1", //20"
    "#phi^{t and #bar{t}}/top quarks/0/4"         ,
    "y^{t and #bar{t}}/top quarks/0/1"             ,//5"
    "p_{T}(hadronic t) #left[#frac{GeV}{c}#right]/events #left[(#frac{GeV}{c})^{-1}#right]/0/20",                         
    "#phi(hadronic t)/events/0/4",
    "y(hadronic t)/events/0/5"    ,
    "p_{T}(leptonic t) #left[#frac{GeV}{c}#right]/events #left[(#frac{GeV}{c})^{-1}#right]/0/20",                         
    "#phi(leptonic t)/events/0/4",
    "y(leptonic t)/events/0/5"   ,
    //  reconstructed angular distributions
    "#angle(b,#bar{b}) (detector rest frame)/events/0/21",
    "#angle(b,#bar{b}) (t#bar{t} rest frame)/events/0/21",
    "#angle(W,W) (t#bar{t} RF)/events/0/21"              ,
    "#angle(lep.t,W) (t:t#bar{t} RF, W:t RF)/events/0/21",
    "#angle(had.t,W) (t:t#bar{t} RF, W:t RF)/events/0/21",
    "#angle(lep.t,b) (t#bar{t} rest frame)/events/0/21"  ,
    "#angle(had.t,b) (t#bar{t} rest frame)/events/0/21"  ,
    "#angle(lep.b,W) (t#bar{t} rest frame)/events/0/21"  ,
    "#angle(had.b,W) (t#bar{t} rest frame)/events/0/21"  ,
    "#angle(q,#bar{q}) (top rest frame)/events/0/21"     ,
    "#angle(q,lep.b) (t#bar{t} rest frame)/events/0/21"  ,
    "#angle(q,had.b) (top rest frame)/events/0/21"       ,
    "#angle(#mu,lep.b) (top rest frame)/events/0/21"     , 
    "#angle(#mu,lep.b) (t#bar{t} rest frame)/events/0/21",
    "#angle(#mu,lep.b) (top rest frame)/events/0/21"     ,
    "#angle(#mu,q) (t#bar{t} rest frame)/events/0/21"    ,
    "#angle(#mu,#nu) (t#bar{t} rest frame)/events/0/21"  ,
    "#angle(#nu,lep.b) (top rest frame)/events/0/21"     ,
    "#angle(#nu,had.b) (t#bar{t} rest frame)/events/0/21",
    "#angle(#nu,q) (t#bar{t} rest frame)/events/0/21"    ,
    "#angle(#mu,lep.W) (#mu:W RF, W:Det RF)/events/0/21" ,
    "#angle(#nu,lep.W) (#nu:W RF, W:Det RF)/events/0/21" ,
    "#angle(q,had.W) (q:W RF, W:Det RF)/events/0/21"     ,
    // reconstructed event shape variables
    "aplanarity/events/0/2"  ,
    "sphericity/events/0/10" ,
    "C/events/0/10"          ,
    "D/events/0/4"           ,
    "circularity/events/0/10",
    "isotropy/events/0/10"   ,
    // generated top quantities
    "m^{t and #bar{t}} parton truth #left[#frac{GeV}{c^{2}}#right]/events/0/10"      ,
    "p_{T}^{t and #bar{t}} #left[#frac{GeV}{c}#right] parton truth/events/0/1"          ,//20"
    "p_{T}^{t and #bar{t}} #left[#frac{GeV}{c}#right] parton truth Phase Space/events/0/1",//20"
    "#phi(t) parton truth/events/0/4"           ,
    "y^{t and #bar{t}} parton truth/events/0/1"              ,//5"
    "p_{T}(hadronic t) #left[#frac{GeV}{c}#right] parton truth/events/0/20",                         
    "#phi(hadronic t) parton truth/events/0/4"  ,
    "y(hadronic t) parton truth/events/0/5"     ,
    "p_{T}(leptonic t) #left[#frac{GeV}{c}#right] parton truth/events/0/20",                         
    "#phi(leptonic t) parton truth/events/0/4"  ,
    "y(leptonic t) parton truth/events/0/5"     ,
    // generated angular distributions
    "#angle(b,#bar{b}) parton truth (detector rest frame)/events/0/21",
    "#angle(b,#bar{b}) parton truth (t#bar{t} rest frame)/events/0/21",
    "#angle(W,W) parton truth (t#bar{t} RF)/events/0/21"              ,
    "#angle(lep.t,W) parton truth (t:t#bar{t} RF, W:t RF)/events/0/21",
    "#angle(had.t,W) parton truth (t:t#bar{t} RF, W:t RF)/events/0/21",
    "#angle(lep.t,b) parton truth (t#bar{t} rest frame)/events/0/21"  ,
    "#angle(had.t,b) parton truth (t#bar{t} rest frame)/events/0/21"  ,
    "#angle(lep.b,W) parton truth (t#bar{t} rest frame)/events/0/21"  ,
    "#angle(had.b,W) parton truth (t#bar{t} rest frame)/events/0/21"  ,
    "#angle(q,#bar{q}) parton truth (top rest frame)/events/0/21"     ,
    "#angle(q,lep.b) parton truth (t#bar{t} rest frame)/events/0/21"  ,
    "#angle(q,had.b) parton truth (top rest frame)/events/0/21"       ,
    "#angle(#mu,lep.b) parton truth (top rest frame)/events/0/21"     , 
    "#angle(#mu,lep.b) parton truth (t#bar{t} rest frame)/events/0/21",
    "#angle(#mu,lep.b) parton truth (top rest frame)/events/0/21"     ,
    "#angle(#mu,q) parton truth (t#bar{t} rest frame)/events/0/21"    ,
    "#angle(#mu,#nu) parton truth (t#bar{t} rest frame)/events/0/21"  ,
    "#angle(#nu,lep.b) parton truth (top rest frame)/events/0/21"     ,
    "#angle(#nu,had.b) parton truth (t#bar{t} rest frame)/events/0/21",
    "#angle(#nu,q) parton truth (t#bar{t} rest frame)/events/0/21"    ,
    "#angle(#mu,W) parton truth (#mu:W RF, W:Det RF)/events/0/21"     ,
    "#angle(q,W) parton truth (q:W RF, W:Det RF)/events/0/21"         ,
    "#angle(#nu,W) parton truth (#nu:W RF, W:Det RF)/events/0/21"     ,
    // 			     // generated event shape variables
    // 			     "aplanarity parton truth/events/0/4"  ,
    // 			     "sphericity parton truth/events/0/10" ,
    // 			     "C parton truth/events/0/10"          ,
    // 			     "D parton truth/events/0/4"           ,
    // 			     "circularity parton truth/events/0/10",
    // 			     "isotropy parton truth/events/0/10"   ,
    // reconstructed ttbar quantities	                            
    "m_{t#bar{t}} #left[#frac{GeV}{c^{2}}#right] /events/1/1"                         ,//60"
    "p_{T}^{t#bar{t}} #left[#frac{GeV}{c}#right]/events/0/1"                     ,//10"
    "y^{t#bar{t}}/events/0/1"                         ,//2
    "H_{T}^{t#bar{t}}=#Sigma(E_{T}(jets)) #left[#frac{GeV}{c}#right]/events/0/20",
    "y^{t and #bar{t}}/events/0/10"                   ,
    "#phi(leptonic t)-#phi(hadronic t)/events/0/4"    ,                
    "y(leptonic t)-y(hadronic t)/events/0/4"          ,  
    // generated ttbar quantities	                            
    "m_{t#bar{t}} #left[#frac{GeV}{c^{2}}#right] parton truth/events/1/1"                         ,//60"
    "p_{T}^{t#bar{t}} #left[#frac{GeV}{c}#right] parton truth/events/0/1"                     ,//10"
    "y^{t#bar{t}} parton truth/events/0/1"                         ,//2
    "H_{T}^{t#bar{t}}=#Sigma(E_{T}(jets)) #left[#frac{GeV}{c}#right] parton truth/events/0/20",
    "y^{t and #bar{t}} parton truth/events/0/10"                   ,
    "#phi(leptonic t)-#phi(hadronic t) parton truth/events/0/4"    ,                
    "y(leptonic t)-y(hadronic t) parton truth/events/0/4"          ,
    // reconstructed lepton quantities
    "p_{T}^{#mu} #left[#frac{GeV}{c}#right]/events #left[(#frac{GeV}{c})^{-1}#right]/0/1",
    "#eta^{#mu}/events/0/1" ,
    // generated lepton quantities
    "p_{T}^{#mu} #left[#frac{GeV}{c}#right] parton truth/events/0/1",
    "#eta^{#mu} parton truth/events/0/1"
  };
  // 2D: "x-axis title"/"y-axis title"
  TString axisLabel2D[ ] = {// reco - gen Match correlation plots (ttbar signal only)
    // a) combinatorics and KinFit Hypothesis Quality(ttbar signal only)
    "i_{lead jet} parton truth/i_{lead jet} hypothesis fit",
    // b) reconstructed Top quantities
    "p_{T}^{t and #bar{t}} #left[#frac{GeV}{c}#right] gen/p_{T}^{t and #bar{t}} #left[#frac{GeV}{c}#right] reco"                    ,
    "#phi^{t and #bar{t}} gen/#phi^{t and #bar{t}} reco"                      ,
    "y^{t and #bar{t}} gen/y^{t and #bar{t}} reco"                            ,
    "angle(b,#bar{b}) gen (t#bar{t} rest frame)/angle(b,#bar{b}) reco (t#bar{t} rest frame)",
    // c) reconstructed ttbar quantities
    "m(t#bar{t}) #left[#frac{GeV}{c^{2}}#right] gen/m(t#bar{t}) #left[#frac{GeV}{c^{2}}#right] reco"              ,
    "p_{T}(t#bar{t}) #left[#frac{GeV}{c}#right] gen/p_{T}(t#bar{t}) #left[#frac{GeV}{c}#right] reco"      ,
    "y(t#bar{t}) gen/y(t#bar{t}) reco"              ,
    "H_{T}(t#bar{t}) #left[#frac{GeV}{c}#right] gen/H_{T}(t#bar{t}) #left[#frac{GeV}{c}#right] reco"      ,
    "#Sigmay(t#bar{t}) gen/#Sigmay(t#bar{t}) reco"  ,
    "#phi(leptonic t)-#phi(hadronic t) gen/#phi(leptonic t)-#phi(hadronic t) Kinfit",
    "y(leptonic t)-y(hadronic t) gen/y(leptonic t)-y(hadronic t) Kinfit"            
  };
  // count # plots
  unsigned int N1Dplots = sizeof(plots1D)/sizeof(TString);
  unsigned int N2Dplots = sizeof(plots2D)/sizeof(TString);
  // check if all axis labels exist
  if(N1Dplots != sizeof(axisLabel1D)/sizeof(TString)) std::cout << "ERROR: some 1D plots or axis label are missing" << std::endl;
  if(N2Dplots != sizeof(axisLabel2D)/sizeof(TString)) std::cout << "ERROR: some 2D plots or axis label are missing" << std::endl;
  if((N1Dplots != sizeof(axisLabel1D)/sizeof(TString))||(N2Dplots != sizeof(axisLabel2D)/sizeof(TString))) exit (1);
  // run automatically in batch mode if there are many canvas
  if((N1Dplots+N2Dplots)>15) gROOT->SetBatch();

  // ---
  //    open our standard analysis files
  // ---
  std::map<unsigned int, TFile*> files_ = getStdTopAnalysisFiles(inputFolder, systematicVariation, dataFile, decayChannel);

  // ---
  //    loading histos
  // ---
  // collect all plot names in vector (first 1D, then 2D)
  std::vector<TString> plotList_;
  plotList_.insert( plotList_.begin(), plots1D, plots1D + sizeof(plots1D)/sizeof(TString) );
  plotList_.insert( plotList_.end()  , plots2D, plots2D + sizeof(plots2D)/sizeof(TString) );
  // container for all histos (1D&2D)
  // example for acess: histo_["plotName"][sampleNr]
  std::map< TString, std::map <unsigned int, TH1F*> > histo_;
  std::map< TString, std::map <unsigned int, TH2F*> > histo2_;
  // total # plots 
  int Nplots=0;
  // save all histos from plotList_ that exist in files_ into 
  // histo_ and histo2_ and count total # of plots as Nplots
  if(verbose>1) std::cout << std::endl;
  std::vector<TString> vecRedundantPartOfNameInData;
  vecRedundantPartOfNameInData.push_back(sysInputFolderExtension);
  getAllPlots(files_, plotList_, histo_, histo2_, N1Dplots, Nplots, verbose-1, decayChannel, &vecRedundantPartOfNameInData);
  // ---
  //    lumiweighting for choosen luminosity
  // ---
  // scale every histo in histo_ and histo2_ to the corresponding luminosity
  // Additionally the mu eff SF is applied
  // NOTE: luminosity [/pb]
  scaleByLuminosity(plotList_, histo_, histo2_, N1Dplots, luminosity, verbose-1, systematicVariation, decayChannel);

  // ---
  //    add single top channels and DiBoson contributions
  // ---
  // for SingleTop and DiBoson samples:
  // every plot plotList_ existing in sTop/ diBoson .root file
  // will be combined and saved in the histo_ and histo2_ map
  // reCreate: recreate combined plots if they are already existing
  bool reCreate=false;
  AddSingleTopAndDiBoson(plotList_, histo_, histo2_, N1Dplots, verbose-1, reCreate, decayChannel);

  // ---
  //    copy event yields for total xSec calculation
  // ---
  TH1F* dataYield =(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topPt"][kData ]->Clone());
  TH1F* SigYield  =(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topPt"][kSig  ]->Clone());
  TH1F* BkgYield  =(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topPt"][kBkg  ]->Clone());
  TH1F* STopYield =(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topPt"][kSTop ]->Clone());
  TH1F* WjetYield =(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topPt"][kWjets]->Clone());
  TH1F* ZjetYield =(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topPt"][kZjets]->Clone());
  TH1F* DiBosYield=(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topPt"][kDiBos]->Clone());
  TH1F* QCDYield  =(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topPt"][kQCD  ]->Clone());
  TH1F* GenInclusive =(TH1F*)(histo_["analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/topPt"          ][kSig]->Clone());
  TH1F* GenPhaseSpace=(TH1F*)(histo_["analyzeTopPartonLevelKinematicsPhaseSpace"+sysInputGenFolderExtension+"/topPt"][kSig]->Clone());

  // ---
  //    configure histogras
  // ---
  // needs: plotList_, histo_, histo2_, N1Dplots, axisLabel_, axisLabel1D, axisLabel2D
  std::vector<TString> axisLabel_;
  axisLabel_.insert( axisLabel_.begin(), axisLabel1D, axisLabel1D + sizeof(axisLabel1D)/sizeof(TString) );
  axisLabel_.insert( axisLabel_.end()  , axisLabel2D, axisLabel2D + sizeof(axisLabel2D)/sizeof(TString) );
  // loop plots
  for(unsigned int plot=0; plot<plotList_.size(); ++plot){
    if(decayChannel=="electron"&&axisLabel_[plot].Contains("#mu")) axisLabel_[plot].ReplaceAll("#mu", "e");
    if(verbose>1){
      std::cout << "(plot, x Axis label , y Axis label , log scale?, rebinning factor):" << std::endl;
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
	histogramStyle( *histo_[plotList_[plot]][sample], sample, true);
	// special configurations
	if(getStringEntry(plotList_[plot], 2)=="PartonJetDRall")histo_[plotList_[plot]][sample]->SetNdivisions(816);
      }
      // b) 2D
      if((plot>=N1Dplots)&&(histo2_.count(plotList_[plot])>0)&&(histo2_[plotList_[plot]].count(sample)>0)) histStyle2D( *histo2_[plotList_[plot]][sample], sampleLabel(sample,decayChannel), getStringEntry(axisLabel_[plot],1), getStringEntry(axisLabel_[plot],2));
    }
  }
  // ---
  //    rebinning 1D histograms
  // ---
  // create variable bin edges
  std::map<TString, std::vector<double> > binning_ = makeVariableBinning();
  // loop samples
  for(unsigned int sample=kSig; sample<=kData; ++sample){
    // loop plots
    for(unsigned int plot=0; plot<plotList_.size(); ++plot){
      TString plotName=plotList_[plot];
      // check if plot exists and is 1D
      if((plot<N1Dplots)&&(plotExists(histo_, plotName, sample))){
	// equidistant binning
	  double reBinFactor = atof(((string)getStringEntry(axisLabel_[plot],4)).c_str());
	if(reBinFactor>1&&binning_.count(plotName)==0){
	  equalReBinTH1(reBinFactor, histo_, plotName, sample);
	}
	// variable binning
	else{
	  // check existence of e.g. binning_["topPt"]
	  if(binning_.count(getStringEntry(plotName, 2))>0){
	    reBinTH1F(*histo_[plotName][sample], binning_[getStringEntry(plotName, 2)], verbose-1);
	    if(verbose>1){
	      for(int i=1; i<= histo_[plotName][sample]->GetNbinsX()+1; i++){
		std::cout << "bin " << i << ": " << histo_[plotName][sample]->GetBinContent(i) << " / ";
		std::cout << histo_[plotName][sample]->GetBinWidth(i) << std::endl;
	      }
	    }
	    // divide by binwidth
	    histo_[plotName][sample] = divideByBinwidth(histo_[plotName][sample], verbose-1);
	    if(verbose>1){
	      for(int i=1; i<= histo_[plotName][sample]->GetNbinsX()+1; i++){
		std::cout << "bin " << i << "= " << histo_[plotName][sample]->GetBinContent(i) << std::endl;
	      }
	    }
	    // change y axis label
	    if(!axisLabel_[plot].Contains("binwidth")){
	      //axisLabel_[plot].ReplaceAll(getStringEntry(axisLabel_[plot],2), getStringEntry(axisLabel_[plot],2)+" per binwidth");
	      // check if unit is GeV
	      if(getStringEntry(axisLabel_[plot],1).Contains("GeV")){
		//axisLabel_[plot].ReplaceAll(getStringEntry(axisLabel_[plot],2), getStringEntry(axisLabel_[plot],2)+" [#frac{c}{GeV}]");
	      }
	    }
	  }
	}
      }
    }
  }
  
  // ---
  //    MC efficiency determination
  // ---
  // count # of efficiency histos
  unsigned int NMCeff=0;
  // create list of variables you would like to create the efficiency / cross section for
  std::vector<TString> xSecVariables_, xSecLabel_;
  TString xSecVariables[] ={"topPt", "topY", "ttbarPt", "ttbarMass", "ttbarY", "lepPt" ,"lepEta"};
  TString xSecLabel    [] ={"p_{T}^{t and #bar{t}}/[#frac{GeV}{c}]" , "y^{t and #bar{t}}/ ", "p_{T}^{t#bar{t}}/[#frac{GeV}{c}]", "m_{t#bar{t}}/[#frac{GeV}{c^{2}}]", "y^{t#bar{t}}/ ", "p_{T}^{#mu}/[#frac{GeV}{c}]" , "#eta^{#mu}/ "};
  xSecVariables_ .insert( xSecVariables_.begin(), xSecVariables, xSecVariables + sizeof(xSecVariables)/sizeof(TString) );
  xSecLabel_     .insert( xSecLabel_    .begin(), xSecLabel    , xSecLabel     + sizeof(xSecLabel    )/sizeof(TString) );
  // loop all variables
  for(unsigned int number=0; number<xSecVariables_.size(); ++number){
    TString variable=xSecVariables_[number];
    TString efficiency="efficiency/"+variable;
    // check if gen and reco plots are available
    if(plotExists(histo_, "analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/"+variable, kSig)&&plotExists(histo_, "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable, kSig)){
      // std::cout << "found gen and reco" << std::endl;
      // get reco plot
      histo_[efficiency][kSig]=(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig]->Clone());
      // divide by gen plot
      histo_[efficiency][kSig]->Divide((TH1F*)(histo_["analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/"+variable][kSig]->Clone()));
//       if(sys==sysHadronisationUp){
// 	double SF=1;
// 	for(int bin=1; bin<=histo_[efficiency][kSig]->GetNbinsX(); ++bin){
// 	  histo_[efficiency][kSig]->SetBinContent(bin histo_[efficiency][kSig]->GetBinContent(bin)*SF);
// 	}
//       }
      // std::cout << "gen(bin1): "  <<  histo_["analyzeTopPartonLevelKinematics/"+variable  ][kSig]->GetBinContent(1) << std::endl;
      // std::cout << "reco(bin1): " <<  histo_["analyzeTopRecoKinematicsKinFit/"+variable][kSig]->GetBinContent(1) << std::endl;
      // std::cout << "eff(bin1): "  <<  histo_[efficiency][kSig]->GetBinContent(1) << std::endl;

      // add plot to list of plots
      plotList_.push_back(efficiency);
      // add axis configuration
      unsigned int positionOfRecoAxisLabel = positionInVector(plotList_, "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable);
      TString recoAxisLabel =axisLabel_[positionOfRecoAxisLabel];
      axisLabel_.push_back(""+getStringEntry(recoAxisLabel,1)+"/"+"#epsilon #times A (MC t#bar{t}#rightarrow#mu)/"+getStringEntry(recoAxisLabel,3)+"/"+getStringEntry(recoAxisLabel,4));
      if(decayChannel=="electron") axisLabel_[axisLabel_.size()-1].ReplaceAll("#mu", "e");
      // set binomial errors
      if(verbose>1) std::cout << "       eff, events, width, sqrt(eff*(1.-eff)/events" << std::endl;
      for(int bin=1; bin<=histo_[efficiency][kSig]->GetNbinsX(); ++bin){
	double eff=histo_[efficiency][kSig]->GetBinContent(bin); 
	double N=histo_["analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/"+variable][kSig]->GetBinContent(bin);
	double width=histo_["analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/"+variable][kSig]->GetBinWidth(bin);
	N*=width*effSFAB(systematicVariation,decayChannel);
	N/=(lumiweight(kSig, luminosity, systematicVariation, decayChannel));
	if(verbose>1){
	  std::cout << "bin " << bin << ": " << eff << ", " << N << ", " << width;
	  std::cout << ", " <<  sqrt(eff*(1.-eff)/N) << std::endl;
	}
	histo_[efficiency][kSig]->SetBinError(bin, sqrt(eff*(1.-eff)/N));
	// set error to 0 for bins with width 0 or 0 entries
	if(N==0) histo_[efficiency][kSig]->SetBinError(bin, 0);
      }
      // change histo style
      histogramStyle(*histo_[efficiency][kSig], kSig, false, 0.6);
      ++NMCeff;
    }
  }

  // ---
  //    differential cross section (extrapolated or PS) determination
  // ---
  // count # of cross section histos
  unsigned int NXSec=0;
  // loop all variables
  for(unsigned int number=0; number<xSecVariables_.size(); ++number){
    TString variable=xSecVariables_[number];
    TString label =getStringEntry(xSecLabel_[number],1);
    TString label2=getStringEntry(xSecLabel_[number],2);
    if(label2=="[#frac{GeV}{c}]")          label2="#left[pb #left(#frac{GeV}{c}#right)^{-1}#right]";
    else if(label2=="[#frac{GeV}{c^{2}}]") label2="#left[pb #left(#frac{GeV}{c^{2}}#right)^{-1}#right]";
    else if(label2==" ") label2="[pb]";
    TString xSec ="xSec/"+variable;
    bool calculateXSec=true;
    // loop samples 
    for(unsigned int sample=kSig; sample<=kData; ++sample){
      // check if reco MC plot and data event yield plots are available
      if(!plotExists(histo_, "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable, sample)){
	calculateXSec=false;
      }
    }
    // check if efficiency plot is available
    if(!plotExists(histo_, "efficiency/"+variable, kSig)){
      calculateXSec=false;
    }
    if(calculateXSec){
      // a) differential XSec from data
      // get data plot
      TString name=TString(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->GetName())+"kData";
      histo_[xSec][kData]=(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->Clone(name));
      // subtract BG(MC)
      //histo_[xSec][kData]->Add((TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kBkg  ]->Clone()), -1); // other ttbar BG is accessed by a signal fraction (see below) 
      histo_[xSec][kData]->Add((TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kZjets]->Clone()), -1);
      histo_[xSec][kData]->Add((TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kWjets]->Clone()), -1);
      histo_[xSec][kData]->Add((TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSTop ]->Clone()), -1);
      histo_[xSec][kData]->Add((TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kDiBos]->Clone()), -1);
      histo_[xSec][kData]->Add((TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kQCD  ]->Clone()), -1);
      // correct for fraction of ttbar events that are signal
      histo_[xSec+"AllTtBar"][kSig]=     (TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig]->Clone());
      histo_[xSec+"AllTtBar"][kSig]->Add((TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kBkg]->Clone()));
      histo_[xSec+"SignalFraction"][kSig]=(TH1F*)histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig]->Clone();
      histo_[xSec+"SignalFraction"][kSig]->Divide(histo_[xSec+"AllTtBar"][kSig]);
      histo_[xSec][kData]->Multiply(histo_[xSec+"SignalFraction"][kSig]);
      // apply efficiency correction
      DivideYieldByEfficiencyAndLumi(histo_[xSec][kData], ((TH1F*)(histo_["efficiency/"+variable][kSig])), luminosity, 0);
      // b) differential XSec from Signal(MC prediction)
      histo_[xSec][kSig]=(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig]->Clone());
      DivideYieldByEfficiencyAndLumi(histo_[xSec][kSig], (TH1F*)(histo_["efficiency/"+variable][kSig]), luminosity, 0);
      // add plot to list of plots
      plotList_.push_back(xSec);
      // add axis configuration
      unsigned int positionOfRecoAxisLabel = positionInVector(plotList_, "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable);
      TString recoAxisLabel =axisLabel_[positionOfRecoAxisLabel];
      recoAxisLabel.ReplaceAll("KinFit ","");
      axisLabel_.push_back(""+getStringEntry(recoAxisLabel,1)+"/"+"#frac{d#sigma}{d"+label+"} "+label2+" (t#bar{t}#rightarrow #mu prompt)/"+getStringEntry(recoAxisLabel,3)+"/"+getStringEntry(recoAxisLabel,4));
      if(decayChannel=="electron") axisLabel_[axisLabel_.size()-1].ReplaceAll("#mu", "e");
      // configure xSec plot histo style
      histogramStyle(*histo_[xSec][kData], kData, false);
      histogramStyle(*histo_[xSec][kSig ], kSig , false );
      // restrict axis
      if(variable=="topPt"){ 
	histo_[xSec][kData]->GetXaxis()->SetRange(1,5);
	histo_[xSec][kSig ]->GetXaxis()->SetRange(1,5);
      }
      if(variable=="topY"){ 
	histo_[xSec][kData]->GetXaxis()->SetRange(2,9);
	histo_[xSec][kSig ]->GetXaxis()->SetRange(2,9);
      }
      if(variable=="lepEta"){ 
	histo_[xSec][kData]->GetXaxis()->SetRange(2, histo_[xSec][kData]->GetNbinsX()-2);
	histo_[xSec][kSig ]->GetXaxis()->SetRange(2, histo_[xSec][kSig ]->GetNbinsX()-2);
      }
      ++NXSec;
    }
  }

  // ---
  //    total inclusive cross section (ttbar->X) determination
  // ---
  double xSecResult=0;
  double effA=0;
  double BR=0;
  ++NXSec;
  // get event yields
  // careful: plot contains leptonic and hadronic top -> *0.5
  double Ndata  = 0.5 * dataYield ->Integral(0, dataYield ->GetNbinsX()+1);
  double NSig   = 0.5 * SigYield  ->Integral(0, SigYield  ->GetNbinsX()+1);
  double NBGtop = 0.5 * BkgYield  ->Integral(0, BkgYield  ->GetNbinsX()+1); 
  double NBGZ   = 0.5 * ZjetYield ->Integral(0, ZjetYield ->GetNbinsX()+1);
  double NBGW   = 0.5 * WjetYield ->Integral(0, WjetYield ->GetNbinsX()+1);
  double NBGsTop= 0.5 * STopYield ->Integral(0, STopYield ->GetNbinsX()+1);
  double NBGVV  = 0.5 * DiBosYield->Integral(0, DiBosYield->GetNbinsX()+1);
  double NBGQCD = 0.5 * QCDYield  ->Integral(0, QCDYield  ->GetNbinsX()+1);
  double NBG= NBGZ+NBGW+NBGsTop+NBGVV+NBGQCD; //other ttbar BG is handled by a signal fraction (see below sigFrac)
  double NAllMC=NSig+NBGtop+NBG;
  double sigFrac=NSig/(NSig+NBGtop);
  // print event composition
  if(verbose>0&&systematicVariation==sysNo){ 
    std::cout << std::endl;
    std::cout << "events expected from MC: " << NAllMC << std::endl;
    std::cout << "expected event composition:"   << std::endl;
    std::cout << "ttbar SG:  " << std::setprecision(4) << std::fixed << NSig   / NAllMC;
    std::cout << " (" << lumiweight(kSig, luminosity, systematicVariation, decayChannel) << "*";
    std::cout << NSig/lumiweight(kSig, luminosity, systematicVariation, decayChannel) << "=" << NSig    << ") " << std::endl;
    std::cout << " (" << "mean event SF: " << NSig/(lumiweight(kSig, luminosity, systematicVariation, decayChannel)*(0.5*SigYield->GetEntries())) << ")" << std::endl;
    std::cout << "ttbar BG:  " << std::setprecision(4) << std::fixed << NBGtop / NAllMC; 
    std::cout << " (" << lumiweight(kBkg, luminosity, systematicVariation, decayChannel) << "*";
    std::cout << NBGtop/lumiweight(kBkg, luminosity, systematicVariation, decayChannel) << "=" << NBGtop << ") " << std::endl;
    std::cout << "single top:" << std::setprecision(4) << std::fixed << NBGsTop/ NAllMC; 
    std::cout << " (" << NBGsTop << ") " << std::endl;
    std::cout << "W+jets    :" << std::setprecision(4) << std::fixed << NBGW   / NAllMC; 
    std::cout << " (" << lumiweight(kWjets, luminosity, systematicVariation, decayChannel) << "*";
    std::cout << NBGW/lumiweight(kWjets, luminosity, systematicVariation, decayChannel) << "=" << NBGW << ") " << std::endl;
    std::cout << "QCD       :" << std::setprecision(4) << std::fixed << NBGQCD / NAllMC; 
    std::cout << " (" << lumiweight(kQCD, luminosity, systematicVariation, decayChannel) << "*";
    std::cout << NBGQCD/lumiweight(kQCD, luminosity, systematicVariation, decayChannel) << "=" << NBGQCD << ") " << std::endl;
    std::cout << "Z+jets    :" << std::setprecision(4) << std::fixed << NBGZ   / NAllMC;
    std::cout << " (" << lumiweight(kZjets, luminosity, systematicVariation, decayChannel) << "*";
    std::cout << NBGZ/lumiweight(kZjets, luminosity, systematicVariation, decayChannel) << "=" << NBGZ << ") " << std::endl;
    std::cout << "Diboson   :" << std::setprecision(4) << std::fixed << NBGVV  / NAllMC; 
    std::cout << " (" << NBGVV << ") " << std::endl;
  }
  // efficiency calculation
  double NGenPhaseSpace=0.5 * GenPhaseSpace->Integral(0, GenPhaseSpace->GetNbinsX()+1);
  double eff=NSig/NGenPhaseSpace;
  // acceptance
  double NGen          =0.5 * GenInclusive ->Integral(0, GenInclusive ->GetNbinsX()+1);
  effA=NSig/NGen;
  double A= effA/eff;
  // branching ratio
  // careful: the xSec here needs to be 
  // consistent with the one in lumiweight()
  BR=NGen*1.0/(ttbarCrossSection*luminosity);
  //  double BR2=12.0/81.0;
  // calculate xSec
  xSecResult= ( Ndata-NBG ) * sigFrac / ( eff*A*luminosity*BR );
  double sigmaxSec = sqrt( Ndata ) * sigFrac / ( eff*A*luminosity*BR );
  // text output
  if(verbose>0){ 
    std::cout << std::endl;
    std::cout << "systematic variation:" << sysLabel(systematicVariation) << std::endl;
    std::cout << "lumi[pb]:" << luminosity << std::endl;
    std::cout << "N(Data):" << Ndata << std::endl;
    std::cout << "N(BG): "  << NBG   << std::endl;
    std::cout << "eff: "    << eff   << std::endl;
    std::cout << "A: "      << A     << std::endl;
    std::cout << "BR MC: "  << BR    << std::endl;
    std::cout << "ttbar sigfrac: " << sigFrac << std::endl;
    std::cout << "inclusive cross section";
    if(extrapolate) std::cout << " (extrapolated)";
    std::cout << " [pb]: ";
    std::cout << xSecResult << " +/- " << sigmaxSec << "(stat.)" << std::endl;
    std::cout << std::endl;
  }
  // create histo
  TString inclName = "xSec/inclusive";
  histo_[inclName][kData] = new TH1F( "inclusivekData", "inclusivekData", 1, 0., 1.0);
  histo_[inclName][kData]->SetBinContent(1, xSecResult     );
  histo_[inclName][kData]->SetBinError  (1, sigmaxSec);
  // add plot to list of all plots
  plotList_.push_back(inclName);
  // configure xSec plot histo style
  histogramStyle(*histo_[inclName][kData], kData, false);
  axisLabel_.push_back(" /#sigma(t#bar{t}#rightarrowX)[pb]/0/1");
  // theory xSec
  histo_[inclName][kSig]=(TH1F*)histo_[inclName][kData]->Clone("xSec/inclusiveTheory");
  histo_[inclName][kSig]->SetBinContent(1, 157.5);
  histo_[inclName][kSig]->SetBinError  (1, 0);
  histogramStyle(*histo_[inclName][kSig ], kSig , false);

  // ---
  //    differential normalized cross section (phase space) determination
  // ---
  // loop all variables
  for(unsigned int number=0; number<xSecVariables_.size(); ++number){
    TString variable=xSecVariables_[number];
    TString label =getStringEntry(xSecLabel_[number],1);
    TString label2=getStringEntry(xSecLabel_[number],2);
    if(label2=="[#frac{GeV}{c}]") label2="#left[#left(#frac{GeV}{c}#right)^{-1}#right]";
    else if(label2=="[#frac{GeV}{c^{2}}]") label2="#left[#left(#frac{GeV}{c^{2}}#right)^{-1}#right]";
    else if(label2==" ") label2="";
    TString xSec ="xSecNorm/"+variable;
    bool calculateXSec=true;
    // loop samples 
    for(unsigned int sample=kSig; sample<=kData; ++sample){
      // check if reco MC plot and data event yield plots are available
      if(!plotExists(histo_, "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable, sample)){
	calculateXSec=false;
      }
    }
    // check if efficiency plot is available
    if(!plotExists(histo_, "efficiency/"+variable, kSig)){
      calculateXSec=false;
    }
    if(calculateXSec){
      // a) differential norm XSec from data
      // get data plot
      TString name=TString(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->GetName())+"kData";
      histo_[xSec][kData]=(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->Clone(name));
      // total number of data events for this variable
      double NdataVariable=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData],verbose-1);
      // total number of BG events for this variable
      double NBGVariable=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kWjets],verbose-1);
      // NBGVariable+=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kBkg],verbose-1);
      // BG from other ttbar events is treated in the signal fraction
      NBGVariable+=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSTop ],verbose-1);
      NBGVariable+=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kDiBos],verbose-1);
      NBGVariable+=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kQCD  ],verbose-1);
      NBGVariable+=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kZjets],verbose-1);
      // ttbar signal fraction
      double NttbarAll=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig],verbose-1);
      NttbarAll+=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kBkg],verbose-1);
      double sigFracVariable=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig],verbose-1);
      sigFracVariable/=NttbarAll;
      // inclusive efficiency (* A) for this variable
      double effIncl= getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable ][kSig],verbose-1);
      effIncl/= getInclusiveXSec(histo_["analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/"+variable][kSig],verbose-1);
      // loop bins
      for(int bin =0; bin<=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->GetNbinsX()+1; ++bin){
	// for inclusive values of xSecResult, effA, BR, Ndata and NBG 
	// which are used in the calculation are taken from section
	// "total inclusive cross section (ttbar->X) determination" above
	// remark: event numbers are taken from yield plots-> *binwidth
	// get data entries in bin

	// get binwidth
	double binwidth=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->GetBinWidth(bin);
	double NdataBin=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->GetBinContent(bin);
	NdataBin*=binwidth;
	// get BG in bin
	double NBGBin=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kZjets]->GetBinContent(bin);
	//NBGBin+=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kBkg]->GetBinContent(bin);
	// BG from other ttbar events is treated in the signal fraction
	NBGBin+=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kWjets]->GetBinContent(bin);
	NBGBin+=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSTop ]->GetBinContent(bin);
	NBGBin+=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kDiBos]->GetBinContent(bin);
	NBGBin+=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kQCD  ]->GetBinContent(bin);
	NBGBin*=binwidth;
	// ttbar signal fraction
	double NttbarAllBin=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig]->GetBinContent(bin);
	NttbarAllBin      +=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kBkg]->GetBinContent(bin);
	double sigFracBin=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig]->GetBinContent(bin);
	sigFracBin/=NttbarAllBin;
	if(!isfinite(sigFracBin)||isnan(sigFracBin)) sigFracBin=1;
	// calculate Ni = number of ttbar only events in bin from data
	double Ni= NdataBin-NBGBin;
	// calculate Ntot = sum_i(Ni) number of signal only events from data for this variable
	double Ntot =NdataVariable-NBGVariable;
	// get eff*A in bin
	double effABin=(histo_["efficiency/"+variable][kSig])->GetBinContent(bin);
	// calculate cross section for this bin
	double xSecBin = (Ni/Ntot) * (sigFracBin/sigFracVariable) * (effIncl/effABin) * (1/binwidth);
	// calculate complete error
 	double Nierror=sqrt(NdataBin); // the error of NBG is treated as systematic variation
	double sumErrorNjSquare = 0; // sum(j,i!=j) [ errorNj*errorNj ]
	for(int bin2=1; bin2<=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->GetNbinsX()+1; ++bin2){
	  if(bin2!=bin){    
	    double errorNj = sqrt( histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->GetBinContent(bin2) 
				   * (histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->GetBinWidth(bin2)) );   
	    // remark: Nj=NdataBinj-NBGBinj, as above, the error of NBGBinj is treated 
	    // as systematic variation and therefore not incorporated here
	    sumErrorNjSquare += errorNj*errorNj;
	  }
	}
	// calculate complete error
	double xSecBinError = sqrt( ((Ntot - Ni)*(Ntot - Ni)*Nierror*Nierror) + (Ni*Ni*sumErrorNjSquare) );
	xSecBinError/= (Ntot*Ntot);
	xSecBinError*= (effIncl/effABin) * (sigFracBin/sigFracVariable) * (1/binwidth);
	// check for ambiguities
	if(!isfinite(xSecBin)||isnan(xSecBin)){
	  if(verbose>1){ 
	    std::cout << "bin " << bin << " in ";
	    std::cout << "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable;
	    std::cout << "gives no valid cross section!" << std::endl;
	    std::cout << "!will set this bin to 0!" << std::endl;
	  }
	  xSecBin=0;	 
	  xSecBinError=0;
	}
	// set underflow/overflow to 0
	// because of infinite binwidth
	if(bin==0||bin==histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->GetNbinsX()+1){
	  xSecBinError=0;
	  xSecBin=0;
	}
	// check calculation for one example bin
	TString testVar = "ttbarPt";
	TString testVar2 = "ttbarMass";
	int testBin=0;
	if(verbose>0&&xSecBin!=0&&(variable==testVar||variable==testVar2)&&bin==testBin){
	  std::cout << std::endl << "differential normalized xSec ";
	  std::cout << variable+" bin " << bin << ":" << std::endl;
	  std::cout << "binwidth: " << std::setprecision(3) << std::fixed << binwidth << std::endl;
	  std::cout << "Ndata(bin): " << std::setprecision(3) << std::fixed << NdataBin << std::endl;
	  std::cout << "NBG(bin): " << std::setprecision(3) << std::fixed << NBGBin << std::endl;
	  std::cout << "eff*A(bin): " << std::setprecision(3) << std::fixed << effABin << std::endl;
	  std::cout << "ttbar sig frac(bin): " << std::setprecision(3) << std::fixed << sigFracBin << std::endl;
	  std::cout << "Ndata(incl): " << std::setprecision(3) << std::fixed << NdataVariable << std::endl;
	  std::cout << "NBG(incl): " << std::setprecision(3) << std::fixed << NBGVariable << std::endl;
	  std::cout << "eff*A(incl): " << std::setprecision(3) << std::fixed << effIncl << std::endl;
	  std::cout << "ttbar sig frac(incl): " << std::setprecision(3) << std::fixed << sigFracVariable << std::endl;
	  std::cout << "result: " << std::setprecision(3) << std::fixed << xSecBin << std::endl;
	  if(extrapolate){
	    double xSecBinFromDiffAndIncl= histo_["xSec/"+variable][kData]->GetBinContent(bin)/(xSecResult*BR);
	    std::cout << "from diff and incl. xSec: ";
	    std::cout << std::setprecision(3) << std::fixed << xSecBinFromDiffAndIncl << std::endl;
	    std::cout << "ratio: " << xSecBin/xSecBinFromDiffAndIncl << std::endl;
	  }
	}
	// apply cross section and error to plot
	histo_[xSec][kData]->SetBinContent(bin, xSecBin     );
	histo_[xSec][kData]->SetBinError  (bin, xSecBinError);
	//saveToRootFile(outputFileName, histo_[xSec][kData], true, verbose-1,"xSec/NormXSec");
      }
    }
    if(verbose>1){
      std::cout << "integral diff.norm.data distribution (" << variable << "): ";
      std::cout << getInclusiveXSec(histo_[xSec][kData],verbose-1) << std::endl;
    }
    // b) differential normalized XSec from Signal(MC prediction)
    histo_[xSec][kSig]=(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig]->Clone());
    DivideYieldByEfficiencyAndLumi(histo_[xSec][kSig], (TH1F*)(histo_["efficiency/"+variable][kSig]), luminosity, 0);
    double XSecInclTheoPS= getInclusiveXSec(histo_[xSec][kSig], verbose-1);
    // exclude underflow and overflow bins because they are negligible and treated wrong
    XSecInclTheoPS-=histo_[xSec][kSig]->GetBinContent(0);
    XSecInclTheoPS-=histo_[xSec][kSig]->GetBinContent(histo_[xSec][kSig]->GetNbinsX()+1);
    histo_[xSec][kSig]->Scale(1/(XSecInclTheoPS));
    if(verbose>1){
      std::cout << "integral diff.norm.MC distribution (" << variable << "): ";
      std::cout << getInclusiveXSec(histo_[xSec][kSig],verbose-1) << std::endl;
    }
    // add plot to list of plots
    plotList_.push_back(xSec);
    // add axis configuration
    unsigned int positionOfRecoAxisLabel = positionInVector(plotList_, "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable);
    TString recoAxisLabel =axisLabel_[positionOfRecoAxisLabel];
    recoAxisLabel.ReplaceAll("KinFit ","");
    axisLabel_.push_back(""+getStringEntry(recoAxisLabel,1)+"/"+"#frac{1}{#sigma}"+" #frac{d#sigma}{d"+label+"} "+label2+" (t#bar{t}#rightarrow #mu prompt"+")/"+getStringEntry(recoAxisLabel,3)+"/"+getStringEntry(recoAxisLabel,4));
    if(decayChannel=="electron") axisLabel_[axisLabel_.size()-1].ReplaceAll("#mu", "e");
    // configure xSec plot histo style
    histogramStyle(*histo_[xSec][kData], kData, false);
    histogramStyle(*histo_[xSec][kSig ], kSig , false);
    // restrict axis
    if(variable=="topPt"){ 
      histo_[xSec][kData]->GetXaxis()->SetRange(1,5);
      histo_[xSec][kSig ]->GetXaxis()->SetRange(1,5);
    }
    if(variable=="topY"){ 
      histo_[xSec][kData]->GetXaxis()->SetRange(2,9);
      histo_[xSec][kSig ]->GetXaxis()->SetRange(2,9);
    }
    if(variable=="lepEta"){ 
      histo_[xSec][kData]->GetXaxis()->SetRange(2, histo_[xSec][kData]->GetNbinsX()-2);
      histo_[xSec][kSig ]->GetXaxis()->SetRange(2, histo_[xSec][kSig ]->GetNbinsX()-2);
    }
    ++NXSec;
  }
  //   std::cout << std::endl;
  //   for(int bin=1; bin<=histo_["xSecNorm/ttbarPt"][kData]->GetNbinsX()+1; ++bin){
  //     double value =histo_["xSecNorm/ttbarPt"][kData]->GetBinContent(bin);
  //     std::cout << "bin " << bin << ": " << value << std::endl;
  //     if(!isfinite(value)||isnan(value)){
  //       histo_["xSecNorm/ttbarPt"][kData]->SetBinContent(bin,0.);
  //       std::cout << "bin " << bin << ": " << histo_["xSecNorm/ttbarPt"][kData]->GetBinContent(bin)<< std::endl;

  //     }	
  //   }

  // ---
  //    create one legend for all 1D histos
  // ---
  unsigned int Nlegends=0;  
  TLegend *leg  = new TLegend();  
  TLegend *leg0 = new TLegend(0.05, 0.15, 1.05, 0.9);
  Nlegends+=2;
  leg ->SetFillStyle(0);
  leg ->SetBorderSize(0);
  leg ->SetTextSize(0.03);
  leg ->SetTextAlign(12);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  leg0->SetHeader("Kin. Fit (after selection)");
  // fill in contributing sample
  // data is to be first entry
  bool TwoThousandElevenData=false;
  TString lumilabel = Form("%3.2f fb^{-1}",luminosity/1000);
  if(luminosity>50.0) TwoThousandElevenData=true;
  else lumilabel=Form("%2.0f pb^{-1}",luminosity);    
  leg ->AddEntry(histo_[plotList_[plotList_.size()-1]][kData], sampleLabel(kData,decayChannel,TwoThousandElevenData),"P");
  leg0->AddEntry(histo_[plotList_[plotList_.size()-1]][kData], sampleLabel(kData,decayChannel,TwoThousandElevenData)+", "+lumilabel,"P");
  // now loop over MC samples
  for(unsigned int sample=kSig; sample<kData; ++sample){
    // check if sampe exists in at least one plot
    bool exit=false;
    // loop plots
    for(unsigned int plot=0; plot<plotList_.size()-1; ++plot){  // <plotList_.size()-1, because last entry is for data (see above)
      // if found: add entry to legend
      if((histo_.count(plotList_[plot])>0)&&(histo_[plotList_[plot]].count(sample)>0)&&(!exit)){
	leg ->AddEntry(histo_[plotList_[plot]][sample], sampleLabel(sample,decayChannel), "F");
	leg0->AddEntry(histo_[plotList_[plot]][sample], sampleLabel(sample,decayChannel), "F");
	exit=true;
      }
    }
  }

  // ---
  //     print information about plots, canvas and legends
  // ---
  if(verbose>1){
    std::cout << std::endl << "# different plots requested: 1D + 2D = ";
    std::cout << N1Dplots + N2Dplots << " ( " << N1Dplots << " + " << N2Dplots << " )" << std::endl;
    std::cout << "1D plots for different samples will be drawn ";
    std::cout << "into same canvas as stacked plot" << std::endl;
    std::cout << "every 2D plot will be drawn into a seperate canvas" << std::endl;
    std::cout << "#Canvas for plots: " << Nplots << std::endl;
    std::cout << "#legends: " << Nlegends << std::endl;
    std::cout << "every legend will be drawn into a seperate canvas" << std::endl;
    std::cout << "#Canvas = #plots + #legends: " << Nplots+Nlegends << std::endl<< std::endl;
  }

  // ---
  //    create canvas
  // ---
  std::vector<TCanvas*> plotCanvas_;
  //  int NCanvas = Nplots+Nlegends;
  // a) create canvas for all plots + legends
  //  for(int sample=0; sample<NCanvas; sample++){
  for(unsigned int sample=0; sample<Nplots+NMCeff+NXSec+Nlegends; sample++){
    char canvname[10];
    sprintf(canvname,"canv%i",sample);    
    plotCanvas_.push_back( new TCanvas( canvname, canvname, 600, 600) );
    //canvasStyle(*plotCanvas_[sample]);
  }

  // ---
  //    change 1D plots into stack plots
  // ---
  // loop plots -> all 1D plots will become stacked plots
  if(verbose>2) std::cout << std::endl;
  for(unsigned int plot=0; plot<plotList_.size(); ++plot){
    createStackPlot(plotList_, histo_, plot, N1Dplots, verbose-1, decayChannel);
  }
  if(verbose>2) std::cout << std::endl;

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
	  if(verbose>1){
	    std::cout << "plotting " << plotList_[plot];
	    std::cout << " from sample " << sampleLabel(sample,decayChannel);
	    std::cout << " to canvas " << canvasNumber << " ( ";
	    std::cout << plotCanvas_[canvasNumber]->GetTitle() << " )" << std::endl;
	  }
	  // for efficiency plots: draw grid
	  if(getStringEntry(plotList_[plot], 1)=="efficiency") plotCanvas_[canvasNumber]->SetGrid(1,1);
	  // for histos with variable binning:
	  if(binning_.count("analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+getStringEntry(plotList_[plot], 2))>0){
	    // get variable binning
	    TString plotName=getStringEntry(plotList_[plot], 2);
	    std::vector<double> binEdges_=binning_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+plotName];
	    // set maximum of histo to last bin considered 
	    // in variable binning (overflow excluded)
	    double firstBin=0;
	    double lastBin=binEdges_.size()-2;
	    // -2 for subtracting overflow bin
	    if(verbose>2){
	      std::cout << "1st bin, last bin: " << firstBin << " , " << lastBin << std::endl;
	    }
	    histo_[plotList_[plot]][sample]->GetXaxis()->SetRange(firstBin, lastBin);
	  }
	  // first plot
	  if(first){
	    // min / max
	    double max = 1.3*histo_[plotList_[plot]][sample]->GetMaximum();
	    // if data file exists
	    if(histo_[plotList_[plot]].count(kData)>0){
	      // and has larger maximum
	      if(max < 1.3*histo_[plotList_[plot]][kData]->GetMaximum()){
		// take this maximum
		max = 1.3*histo_[plotList_[plot]][kData]->GetMaximum();
		if(verbose>1) std:: cout << "take max from data! " << std::endl;
	      }
	    }
	    double min = 0;
	    // log plots
	    if(getStringEntry(axisLabel_[plot],3)=="1"){
	      plotCanvas_[canvasNumber]->SetLogy(1);
	      min=1;
	      max=exp(1.3*(std::log(max)-std::log(min))+std::log(min));
	      if(plotList_[plot]=="analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/prob") min=0.1; 
	      if(plotList_[plot].Contains("ttbarMass")&&plotList_[plot].Contains("xSec")){
		min=0.0001;
		max=1.2*exp(1.3*(std::log(max)-std::log(min))+std::log(min));
		if(plotList_[plot].Contains("Norm")){
		  min=0.00001;
		  max=0.06;
		}
	      }
	    }
	    // get nicer int values if maximum is large enough
	    if(max>3) max = (double)roundToInt(max);
	    // axis style
	    axesStyle(*histo_[plotList_[plot]][sample], getStringEntry(axisLabel_[plot],1), getStringEntry(axisLabel_[plot],2), min, max); 
	    histo_[plotList_[plot]][sample]->GetXaxis()->SetNoExponent(true);
	    if(max>1&&max<100) histo_[plotList_[plot]][sample]->GetYaxis()->SetNoExponent(true);
	    else histo_[plotList_[plot]][sample]->GetYaxis()->SetNoExponent(false);
	    if(getStringEntry(plotList_[plot], 1).Contains("xSec")) histo_[plotList_[plot]][sample]->GetYaxis()->SetTitleOffset(1.6);
	    // restrict x axis for different plots
	    if(getStringEntry(plotList_[plot], 2)=="topMass") histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(0,500);
	    if(getStringEntry(plotList_[plot], 2)=="lepEta") histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(-2.0,2.0);
	    if(!(plotList_[plot].Contains("xSec"))&&(getStringEntry(plotList_[plot], 2)=="topY"   ||
						     getStringEntry(plotList_[plot], 2)=="topYHad"|| getStringEntry(plotList_[plot], 2)=="topYLep")){
	      histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(-3,3);
	    }
	    if(getStringEntry(plotList_[plot], 2)=="lepPt" ) histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(30,190);
	    if(getStringEntry(plotList_[plot], 2)=="ttbarY") histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(-1,1);
	    if(getStringEntry(plotList_[plot], 2)=="PartonJetDRall") histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(0,4);
	    // draw efficiency plots as line
	    if(getStringEntry(plotList_[plot], 1)=="efficiency") histo_[plotList_[plot]][sample]->Draw("p e1");
	    // others as histo (stack)
	    else histo_[plotList_[plot]][sample]->Draw("hist");
	    histo_[plotList_[plot]][42] = (TH1F*)(histo_[plotList_[plot]][sample]->Clone());
	  }
	  // draw other plots into same canvas
	  else{ 
	    // draw data as points
	    if(sample==kData){ 
	      histo_[plotList_[plot]][sample]->Draw("p e1 X0 same");
	    }
	    else{
	      // draw efficiencies as points
	      if(getStringEntry(plotList_[plot], 1)=="efficiency") histo_[plotList_[plot]][sample]->Draw("p e1 same");
	      // draw others as histo (stack)
	      else histo_[plotList_[plot]][sample]->Draw("hist same");
	    }
	  }
	  first=false;
	  // draw legend for recoYield plots
	  TString title=plotCanvas_[canvasNumber]->GetTitle();
	  if(title.Contains("analyzeTopRecoKinematicsKinFit")){
	    leg->SetX1NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.20);
	    leg->SetY1NDC(1.0 - gStyle->GetPadTopMargin()   - gStyle->GetTickLength() - 0.24);
	    leg->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength());
	    leg->SetY2NDC(1.0 - gStyle->GetPadTopMargin()   - gStyle->GetTickLength());
	    leg->Draw("SAME");
	  }
	  // redraw axis at the end
	  if((histo_.count(plotList_[plot])>0)&&(sample==kData)) histo_[plotList_[plot]][42]->Draw("axis same");	 
	  TString plotType=getStringEntry(plotList_[plot], 1);
	  if(plotType.Contains("xSec")||plotType.Contains("Reco")){
	    if (decayChannel=="muon") DrawDecayChLabel("#mu + Jets");
	    else if(decayChannel=="electron") DrawDecayChLabel("e + Jets");
	    DrawCMSLabels(true,luminosity);
	  }
	}
      }	     
      
      // a3) for 2D plots
      if((plot>=N1Dplots)&&(histo2_.count(plotList_[plot])>0)&&(histo2_[plotList_[plot]].count(sample)>0)){
	// new Canvas for every plot
	plotCanvas_[canvasNumber]->cd(0);
	plotCanvas_[canvasNumber]->SetRightMargin(myStyle.GetPadRightMargin());
	plotCanvas_[canvasNumber]->SetTitle(getStringEntry(plotList_[plot],2)+getStringEntry(plotList_[plot],1)+getTStringFromInt(sample)); 
	if(verbose>1){
	  std::cout << "plotting " << plotList_[plot];
	  std::cout << " from sample " << sampleLabel(sample,decayChannel);
	  std::cout << " to canvas " << canvasNumber  << " ( ";
	  std::cout << plotCanvas_[canvasNumber]->GetTitle() << " )"  << std::endl;
	}
	// draw histo
	plotCanvas_[canvasNumber]->SetRightMargin(0.15);
	histo2_[plotList_[plot]][sample]->GetXaxis()->SetNoExponent(true);
	histo2_[plotList_[plot]][sample]->GetYaxis()->SetNoExponent(true);
	histo2_[plotList_[plot]][sample]->Draw("colz");
	// print correlation factor
	double d = histo2_[plotList_[plot]][sample]->GetCorrelationFactor();
	char correlation[20];
	sprintf(correlation, "%f", d);
	TString corr = (TString)correlation;
	DrawLabel("correlation: "+corr, 0.35, 0.92, 0.75, 0.99, 0.7);
	++canvasNumber;
      }
    }    

    // for 1D hists and efficiency hists: next canvas
    if(((plot<N1Dplots)||(plot>=N1Dplots+N2Dplots))&&(histo_.count(plotList_[plot])>0)) ++canvasNumber;
  }
  // b) for legends
  plotCanvas_[canvasNumber]->cd(0);
  plotCanvas_[canvasNumber]->SetTitle("legend");
  leg0->Draw("");
  ++canvasNumber;

  // ---
  // saving
  // ---
  if(save){
    // pdf and eps only for standard analysis without variation
    if(systematicVariation==sysNo){
      if(verbose<=1) gErrorIgnoreLevel=kWarning;
      // a) as pdf
      saveCanvas(plotCanvas_, outputFolder, pdfName, true, false);
      // b) as eps and png
      for(unsigned int idx=0; idx<plotCanvas_.size(); idx++){
	TString saveToFolder=outputFolder;
	TString title=(plotCanvas_[idx])->GetTitle();
	if(title.Contains("efficiency"                     )) saveToFolder+="effAndAcc/";
	if(title.Contains("analyzeTopPartonLevelKinematics")) saveToFolder+="partonLevel/";
	if(title.Contains("analyzeHypoKinFit"              )) saveToFolder+="kinFitPerformance/";
	if(title.Contains("xSec"                           )) saveToFolder+="xSec/";

	if(title.Contains("analyzeTopRecoKinematicsKinFit" )) saveToFolder+="recoYield/";
	if(title.Contains("0")                              ) saveToFolder=outputFolder+"genRecoCorrPlots/";
	if(!title.Contains("canv")){
	  plotCanvas_[idx]->Print(saveToFolder+(TString)(plotCanvas_[idx]->GetTitle())+".eps"); 
	  plotCanvas_[idx]->Print(saveToFolder+(TString)(plotCanvas_[idx]->GetTitle())+".png");
	}
      }
    }
    // c) root file
    if(verbose>0) std::cout << "will save all plots in outputfile named " << outputFileName << std::endl;
    for(unsigned int idx=0; idx<plotCanvas_.size(); ++idx){
      // get correct folder
      TString title=(plotCanvas_[idx])->GetTitle();
      TString outputfolder="";
      TString possibleFolderNames[5]={"efficiency", "analyzeTopRecoKinematicsKinFit", "xSec", "analyzeTopPartonLevelKinematics", "analyzeHypoKinFit"};
      for(unsigned int name=0; name<sizeof(possibleFolderNames)/sizeof(TString); ++name){
	// add another subfoder indicating the systematic variation
	if(title.Contains(possibleFolderNames[name])){ 
	  outputfolder=possibleFolderNames[name]+"/"+sysLabel(systematicVariation);
	  plotCanvas_[idx]->SetTitle(title.ReplaceAll(possibleFolderNames[name],""));
	}
      }
      // save only canvas from selected subfolders or legends
      if(outputfolder!=""||title.Contains("legend")) saveToRootFile(outputFileName, plotCanvas_[idx], true, verbose-1,outputfolder);
    }
  }
  
  // delete pointer
  delete leg0;
}

