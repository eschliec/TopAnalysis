#ifndef basicFunctions_h
#define basicFunctions_h

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <map>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <typeinfo>

#include <TROOT.h>
#include <TSystem.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TKey.h>
#include <TGraphAsymmErrors.h>
#include <TString.h>

#include <TLine.h>
#include <TBox.h>
#include <TPaveLabel.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TF1.h>
#include <TBox.h>
#include <TGaxis.h>
#include <TError.h>
#include <TMath.h>

#include "HHStyle.h"

#ifdef DILEPTON_MACRO
namespace semileptonic {
#endif

  // ==========================================
  //  Cross-section variables and labels
  // ==========================================

  // basic variables
 
  TString xSecVariables[]     = {"topPt",     "topY",     "ttbarPt",     "ttbarY",     "ttbarMass",     "lepPt" ,    "lepEta",     "bqPt",     "bqEta"    };
  TString xSecVariablesNorm[] = {"topPtNorm", "topYNorm", "ttbarPtNorm", "ttbarYNorm", "ttbarMassNorm", "lepPtNorm" ,"lepEtaNorm", "bqPtNorm", "bqEtaNorm"}; 
  TString xSecVariablesIncl[] = {"inclusive"};

  TString xSecLabel[] = {"p_{T}^{t}/[GeV]" ,"y^{t}/ ",
			 "p_{T}^{t#bar{t}}/[GeV]", "y^{t#bar{t}}/ ", "m^{t#bar{t}}/[GeV]", 
			 "p_{T}^{l}/[GeV]", "#eta^{l}/ ", 
			 "p_{T}^{b}/[GeV]", "#eta^{b}/ "};

  // cross-check variables
 
  TString xSecVariablesCCVar[]     = {"topPtPlus",    "topPtMinus",    "topYPlus",    "topYMinus",    "lepEtaPlus",    "lepEtaMinus"    };
  TString xSecVariablesCCVarNorm[] = {"topPtPlusNorm","topPtMinusNorm","topYPlusNorm","topYMinusNorm","lepEtaPlusNorm","lepEtaMinusNorm"};

  TString xSecLabelCCVar[] = {"p_{T}^{t}/[GeV]", "p_{T}^{#bar{t}}/[GeV]", "y^{t}/ ", "y^{#bar{t}}/ ", "#eta^{#mu^{+}}/ ", "#eta^{#mu^{-}}/ "};

  // =================================================
  //  Parameters for samples and style options 
  // =================================================
  
  // Enumerator to assign an unique integer value to each sample

                   /*0:*/    /*1:*/    /*2:*/    /*3:*/    
  enum samples    {kSig    , kBkg    , kZjets  , kWjets  , 
		   /*4:*/    /*5:*/    /*6:*/    /*7:*/  
		   kQCD    , kSTop   , kDiBos  , kData,   
		   /*8*/     /*9*/     /*10*/    /*11*/    /*12*/    /*13*/
		   kQCDEM1 , kQCDEM2 , kQCDEM3 , kQCDBCE1, kQCDBCE2, kQCDBCE3,  
		   /*14*/    /*15*/    /*16*/
		   kWW     , kWZ     , kZZ     , 
		   /*17*/    /*18*/    /*19*/    /*20*/    /*21*/    /*22*/
		   kSTops  , kSATops , kSTopt  , kSATopt , kSToptW , kSAToptW,
		   /*23*/    /*24*/    /*25*/    /*26*/
		   kSigPow , kBkgPow,  kSigMca , kBkgMca,
                   /*27*/
		   ENDOFSAMPLEENUM};

  // Colors for event samples

  int color_[] ={ kRed    , kRed-7  , kAzure-2, kGreen-3, 
		  kYellow , kMagenta, 10      , kBlack  , 
		  kYellow , kYellow , kYellow , kYellow , kYellow , kYellow ,
		  10      , 10      , 10      , 
		  kMagenta, kMagenta, kMagenta, kMagenta, kMagenta, kMagenta };

  // Colors for data and theory curves

  const unsigned int constDataColor     = kBlack;
  const unsigned int constMadgraphColor = kRed;
  const unsigned int constMcatnloColor  = kBlue;
  const unsigned int constPowhegColor   = kGreen+1;
  const unsigned int constNnloColor     = kOrange-3;

  // Line style for theory curves

  const unsigned int constPowhegStyle  = 7;
  const unsigned int constNnloStyle    = 5;

  // Marker style

  int marker_[] = {20, 22, 29, 23, 
		   21, 27, 28, 20, 
		   21, 21, 21, 21, 21, 21,
		   28, 28, 28,
		   27, 27, 27, 27, 27, 27};

  // Enumerator to assign an unique integer value to each systematic variation

  enum systematicVariation { /* 0:*/ sysNo,
			     /* 1:*/ sysLumiUp,                  /* 2:*/ sysLumiDown,                
			     /* 3:*/ sysPUUp,                    /* 4:*/ sysPUDown,                  
			     /* 5:*/ sysJESUp,                   /* 6:*/ sysJESDown,                 
			     /* 7:*/ sysJERUp,                   /* 8:*/ sysJERDown,                 
			     /* 9:*/ sysLepEffSFNormUp,          /*10:*/ sysLepEffSFNormDown,        
			     /*11:*/ sysLepEffSFShapeUpEta,      /*12:*/ sysLepEffSFShapeDownEta,
			     /*13:*/ sysLepEffSFShapeUpPt,       /*14:*/ sysLepEffSFShapeDownPt,  
			     //      sysLepEffSFShapeUpPt40,             sysLepEffSFShapeDownPt40,
		             //      sysLepEffSFNormUpStat,              sysLepEffSFNormUpStat,
			     /*15:*/ sysTriggerEffSFJetNormUp,   /*16:*/ sysTriggerEffSFJetNormDown,
			     /*17:*/ sysTriggerEffSFJetShapeUp,  /*18:*/ sysTriggerEffSFJetShapeDown,    
			     /*19:*/ sysBtagSFUp,                /*20:*/ sysBtagSFDown, 
			     // old: sysBtagSFShapeUpPt100,              sysBtagSFShapeDownPt100,
			     //      sysBtagSFShapeUpEta1p2,             sysBtagSFShapeDownEta1p2,
			     /*21:*/ sysBtagSFShapeUpPt65,       /*22:*/ sysBtagSFShapeDownPt65,
			     /*23:*/ sysBtagSFShapeUpEta0p7,     /*24:*/ sysBtagSFShapeDownEta0p7,
			     /*25:*/ sysMisTagSFUp,              /*26:*/ sysMisTagSFDown,
			     /*27:*/ sysTopScaleUp,              /*28:*/ sysTopScaleDown,            
			     /*29:*/ sysVBosonScaleUp,           /*30:*/ sysVBosonScaleDown,          
			     /*31:*/ sysSingleTopScaleUp,        /*32:*/ sysSingleTopScaleDown,     
			     /*33:*/ sysTopMatchUp,              /*34:*/ sysTopMatchDown,            
			     /*35:*/ sysVBosonMatchUp,           /*36:*/ sysVBosonMatchDown,         
			     /*37:*/ sysTopMassUp,               /*38:*/ sysTopMassDown,            
			     /*39:*/ sysQCDUp,                   /*40:*/ sysQCDDown,                 
			     /*41:*/ sysSTopUp,                  /*42:*/ sysSTopDown,               
			     /*43:*/ sysDiBosUp,                 /*44:*/ sysDiBosDown,              
			     /*45:*/ sysPDFUp,                   /*46:*/ sysPDFDown,
			     /*47:*/ sysHadUp,                   /*48:*/ sysHadDown,
			     /*49:*/ sysGenMCatNLO,              /*50:*/ sysGenPowheg,
			     /*51:*/ sysShapeUp,                 /*52:*/ sysShapeDown,
			     /*53:*/ ENDOFSYSENUM,               /*54:*/ sysTest };

  // ============================
  //  Numerical Constants
  // ============================

  const double ttbarCrossSection=165.6;                      // NNNLO Kidonakis, recalculated for mtop=172.5 GeV (cf. TOP-11-008)
  const double ttbarCrossSectionError=sqrt(6.2*6.2+9.1*9.1); // Scale and PDF uncertainties on NNLO value
                                                             // --> the scale contributions has been symetrized

  const double SF_TopMassDownUncertainty=0.9/11.0; // scale factors for top mass uncertainty
  const double SF_TopMassUpUncertainty  =0.9/12.0; // --> world average is presently known at +/-0.9 GeV (arXiv:1107.5255v3 [hep-ex])
                                                   // --> systematic samples are varied by -11.0/+12.0 GeV 
                                                   // --> linearily rescale uncertainty on top mass in combineTopDiffXSecUncertainties.C

  const double constHadUncertainty   = 0.050; // relative uncertainty // outdated and only used as placeholder for bquark quantities
  const double globalLumiUncertainty = 0.022; // relative uncertainty
	
  const double constLumiElec = 4980.0;	
  const double constLumiMuon = 4955.0;
  
  const double BRPDG=0.145888;

  // ============== 
  //  Functions
  // ============== 

  TString sysLabel(unsigned int sys)
  {
    // this function returns a TString that corresponds
    // to the systematic variation "sys" of the enumerator "systematicVariation"
    // modified quantities: none
    // used functions: none
    // used enumerators: none (label correspond to systematicVariation)

    switch (sys)
    {
      case sysNo                       : return "sysNo";
      case sysLumiUp                   : return "sysLumiUp";
      case sysLumiDown                 : return "sysLumiDown";                
      case sysPUUp                     : return "sysPUUp";                  
      case sysPUDown                   : return "sysPUDown";                  
      case sysJESUp                    : return "sysJESUp";
      case sysJESDown                  : return "sysJESDown";
      case sysJERUp                    : return "sysJERUp";                 
      case sysJERDown                  : return "sysJERDown";                
      case sysLepEffSFNormUp           : return "sysLepEffSFNormUp";
      case sysLepEffSFNormDown         : return "sysLepEffSFNormDown";
      case sysLepEffSFShapeUpEta       : return "sysLepEffSFShapeUpEta";
      case sysLepEffSFShapeDownEta     : return "sysLepEffSFShapeDownEta";
      case sysLepEffSFShapeUpPt        : return "sysLepEffSFShapeUpPt";
      case sysLepEffSFShapeDownPt      : return "sysLepEffSFShapeDownPt"; 
   // case sysLepEffSFShapeUpPt40      : return "sysLepEffSFShapeUpPt40";
   // case sysLepEffSFShapeDownPt40    : return "sysLepEffSFShapeDownPt40";
   // case sysLepEffSFNormUpStat       : return "sysLepEffSFNormUpStat";
   // case sysLepEffSFNormDownStat     : return "sysLepEffSFNormDownStat";
      case sysTriggerEffSFJetNormUp    : return "sysTriggerEffSFJetNormUp";
      case sysTriggerEffSFJetNormDown  : return "sysTriggerEffSFJetNormDown";
      case sysTriggerEffSFJetShapeUp   : return "sysTriggerEffSFJetShapeUp";
      case sysTriggerEffSFJetShapeDown : return "sysTriggerEffSFJetShapeDown"; 
      case sysBtagSFUp                 : return "sysBtagSFUp";
      case sysBtagSFDown               : return "sysBtagSFDown";
   // case sysBtagSFShapeUpPt100       : return "sysBtagSFShapeUpPt100";
   // case sysBtagSFShapeDownPt100     : return "sysBtagSFShapeDownPt100";
   // case sysBtagSFShapeUpEta1p2      : return "sysBtagSFShapeUpEta1p2";
   // case sysBtagSFShapeDownEta1p2    : return "sysBtagSFShapeDownEta1p2";
      case sysBtagSFShapeUpPt65        : return "sysBtagSFShapeUpPt65";
      case sysBtagSFShapeDownPt65      : return "sysBtagSFShapeDownPt65";
      case sysBtagSFShapeUpEta0p7      : return "sysBtagSFShapeUpEta0p7";
      case sysBtagSFShapeDownEta0p7    : return "sysBtagSFShapeDownEta0p7"; 
   // case sysBTagSFFullShapeUpPt65    : return "sysBTagSFFullShapeUpPt65";
   // case sysBTagSFFullShapeDownPt65  : return "sysBTagSFFullShapeDownPt65";
   // case sysBTagSFFullShapeUpEta0p7  : return "sysBTagSFFullShapeUpEta0p7";
   // case sysBTagSFFullShapeDownEta0p7: return "sysBTagSFFullShapeDownEta0p7"; 
      case sysMisTagSFUp               : return "sysMisTagSFUp";
      case sysMisTagSFDown             : return "sysMisTagSFDown";
      case sysTopScaleUp               : return "sysTopScaleUp";
      case sysTopScaleDown             : return "sysTopScaleDown";
      case sysVBosonScaleUp            : return "sysVBosonScaleUp";
      case sysVBosonScaleDown          : return "sysVBosonScaleDown";
      case sysSingleTopScaleUp         : return "sysSingleTopScaleUp";
      case sysSingleTopScaleDown       : return "sysSingleTopScaleDown";
      case sysTopMatchUp               : return "sysTopMatchUp";
      case sysTopMatchDown             : return "sysTopMatchDown";            
      case sysVBosonMatchUp            : return "sysVBosonMatchUp";          
      case sysVBosonMatchDown          : return "sysVBosonMatchDown";
      case sysTopMassUp                : return "sysTopMassUp";  
      case sysTopMassDown              : return "sysTopMassDown";            
      case sysQCDUp                    : return "sysQCDUp";           
      case sysQCDDown                  : return "sysQCDDown";
      case sysSTopUp                   : return "sysSTopUp";               
      case sysSTopDown                 : return "sysSTopDown";
      case sysDiBosUp                  : return "sysDiBosUp";              
      case sysDiBosDown                : return "sysDiBosDown";
      case sysPDFUp                    : return "sysPDFUp";              
      case sysPDFDown                  : return "sysPDFDown";
      case sysHadUp                    : return "sysHadronizationUp";
      case sysHadDown                  : return "sysHadronizationDown";
      case sysGenMCatNLO               : return "sysGenMCatNLO";
      case sysGenPowheg                : return "sysGenPowheg";  
      case sysShapeUp                  : return "sysShapeUp";              
      case sysShapeDown                : return "sysShapeDown"; 
      default                          : std::cout << "ERROR: the chosen input for function sysLabel is not valid" << std::endl;
                                         std::cout << "chosen systematic variation:  " << sys            << std::endl;
                                         std::cout << "maximum systematic variation: " << ENDOFSYSENUM-1 << std::endl;
                                         return "";
    }
  }

  TString sampleLabel(unsigned int sample, const std::string decayChannel)
  {
    // this function returns the name of the entered MC process
    // corresponding to the enumerator entry "sample" as defined in samples
    // modified quantities: NONE
    // used functions:      NONE
    // used enumerators:    sample
    
    // list all MC process/sample names
    if(sample==kSig    ) return "t#bar{t} Signal";
    if(sample==kSigPow ) return "t#bar{t} Signal POWHEG";
    if(sample==kSigMca ) return "t#bar{t} Signal MC@NLO";
    if(sample==kBkg    ) return "t#bar{t} Other";
    if(sample==kBkgPow ) return "t#bar{t} Other POWHEG";
    if(sample==kBkgMca ) return "t#bar{t} Other MC@NLO";
    if(sample==kSTop   ) return "Single Top";
    if(sample==kSToptW ) return "Single Top tW";
    if(sample==kSTops  ) return "Single Top s";
    if(sample==kSTopt  ) return "Single Top t";
    if(sample==kSAToptW) return "Single Antitop tW";
    if(sample==kSATops ) return "Single Antitop s";
    if(sample==kSATopt ) return "Single Antitop t";
    if(sample==kWjets  ) return "W+Jets";
    if(sample==kZjets  ) return "Z+Jets";
    if(sample==kDiBos  ) return "Diboson";
    if(sample==kWW     ) return "WW";
    if(sample==kWZ     ) return "WZ";
    if(sample==kZZ     ) return "ZZ";
    if(sample==kQCD    ) return "QCD Multijet";
    if(sample==kQCDEM1 ) return "QCDEM1";
    if(sample==kQCDEM2 ) return "QCDEM2";
    if(sample==kQCDEM3 ) return "QCDEM3";
    if(sample==kQCDBCE1) return "QCDBCE1"; 
    if(sample==kQCDBCE2) return "QCDBCE2";
    if(sample==kQCDBCE3) return "QCDBCE3"; 
    if(sample==kData   ) return "Data";
    else return "Default "+decayChannel;
  }

  double effSFAB(int sys=sysNo, std::string decayChannel="unset")
  {
    // this function returns the muon eff SF
    // as derived from Z->mumu tag and probe method
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE
    // "sys": if sysMuEffSFup/sysMuEffSFdown is
    // chosen, the corresponding systematic shifted
    // SF is returned

    if(sys>=ENDOFSYSENUM&&sys!=sysTest){
      std::cout << "ERROR in effSFAB:" << std::endl;
      std::cout << "sys must be smaller than " << ENDOFSYSENUM << std::endl;
    }

    // SF is now applied as event weight in analyzer
    //return 1.0;

    // combined single muon SF Run A+B from tag and probe
    double result = -1.;
    if (decayChannel.compare("muon")==0) {
      result = 1.0; // SF is now applied as event weight in analyzer
    }
    else if (decayChannel.compare("electron")==0) {
      // result = 0.9658; 
      result = 1.0; // SF is now applied as event weight in analyzer
    }
    else if (decayChannel.compare("combined")==0) result = 1.0;
    // errors for the derived SF
    //double errorUp   = 0.03*result;
    //double errorDown = 0.03*result;
    // return requestet SF
    // SF uncertainty is now applied as event weight in analyzer
    //if(sys==sysMuEffSFUp  ) result+=errorUp;
    //if(sys==sysMuEffSFDown) result-=errorDown;
    return result;
  }
  
  // BR correction for ttbar->lnuqq'bb'
  // used for ttbar SG and BG (which is mainly from tau decays)
  double BRcorrectionSemileptonic = 0.985608;
  
  void histogramStyle(TH1& hist, int sampleType, bool filled=true, double markersize=1.2, unsigned int color=0)
  {
    // this function configures the style of a TH1 histogram "hist"
    // using "sampleTyp" to identify the corresponding sample from
    // the enumerator samples definition and assign the corresponding
    // marker and color
    // modified quantities: hist
    // used functions: marker_, color_
    // used enumerators: samples
    // "filled": 0 - line only
    //           1 - area under plot filled

    hist.SetStats(kFALSE);
    if(sampleType==kData || !filled){
      if(!filled)hist.SetLineWidth(2);
      hist.SetLineColor(color_[sampleType]);
      hist.SetMarkerColor(color_[sampleType]);
      if(sampleType==kQCD){
	hist.SetLineColor(kYellow+1);
	hist.SetMarkerColor(kYellow+1);
      }
      if(sampleType==kDiBos){
	hist.SetLineColor(kGray);
	hist.SetMarkerColor(kGray);
      }
      hist.SetMarkerStyle(marker_[sampleType]);
      hist.SetMarkerSize(markersize);
      hist.SetFillStyle(0);
    }
    else{
      hist.SetFillColor(color_[sampleType]);
    }
    if(color!=0){
      hist.SetLineColor(color);
      hist.SetFillColor(color);
      hist.SetMarkerColor(color);
    }
  }

  void histStyle2D(TH2& hist, const TString titleHisto, const TString titleX, const TString titleY)
  {
    // this function configures the style of a TH2 histogram "hist"
    // modified quantities: hist
    // used functions: NONE
    // used enumerators: NONE
    // "titleHisto": set tite of the histo (drawn in upper left)
    // "titleX": title of the x-axis
    // "titleY": title of the y-axis

    hist.SetTitle(titleHisto);
    if(titleX!="") hist.GetXaxis()->SetTitle(titleX);
    if(titleY!="") hist.GetYaxis()->SetTitle(titleY);
    hist.SetStats(kFALSE);
  }

  double readLineFromFile(int line, TString file="crossSectionCalculation.txt")
  {
    // this function to reads and returns a double value
    // from a specific line "line" of the .txt-like file "file"
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE

    // define variables
    std::ifstream finDouble (file);
    std::string readIn;
    // check if file exists
    if (!finDouble){
      std::cout << "can not open file " << file << std::endl;
      return -1;
    }
    // loop lines of the file
    for(int l=1; !finDouble.eof(); ++l){
      // save line content in readIn
      getline(finDouble, readIn);
      // convert your chosen line into double and return it
      if(l==line) return atof(readIn.c_str());
    }
    // if line is not found
    std::cout << "can not find line" << line << std::endl;
    return -1.;
  }

  //void canvasStyle(TCanvas& canv) re-definition of header to avoid warnings when compiling, must be changed once function body is not empty anymore
  void canvasStyle()
  {
    // function is not called anymore - just kept if indiviudal configurations might be required
    // could be considered to be deleted
  }

  void axesStyle(TH1& hist, const char* titleX, const char* titleY, float yMin=-123, float yMax=-123)
  {
    // this function configures the axis style of a TH1 histogram "hist"
    // modified quantities: hist
    // used functions: NONE
    // used enumerators: NONE
    // "titleX": title of the x-axis
    // "titleY": title of the y-axis
    // "yMin": minimum of the y axis
    // "yMax": maximum of the y axis

    hist.SetTitle("");
    hist.GetXaxis()->SetTitle(titleX);
    hist.GetYaxis()->SetTitle(titleY);
    if(yMin!=-123) hist.SetMinimum(yMin);
    if(yMax!=-123) hist.SetMaximum(yMax);
  }

  void legendStyle(TLegend& leg, TString header="", double textsize=0.03, unsigned int aligncode=12)
  {
    // this function configures the style of legends
    // modified quantities: leg
    // used functions: NONE
    // used enumerators: NONE
    // "header": header of the legend

    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.SetTextSize(textsize);
    leg.SetTextAlign(aligncode);
    if (header!="") leg.SetHeader(header); // setting the header adds additional row even if contents is empty (""), important when using TLegend::GetNRows()
  }

  template <class T>
    void writeToFile(T output, TString file, bool append)
    {
      // this function writes the content of "output"
      // into "file" which is of a type like .txt
      // "output" might be of type TString or double
      // modified quantities: file
      // used functions: NONE
      // used enumerators: NONE
      // "append": if false, "file" will be recreated and
      // the existing "file" will be deleted

      // a) write into file
      if(!append){
	std::ofstream fout(file);
	fout << output << std::endl;
	fout.close();
      }
      // b) add output to the end of the file
      if(append){
	std::ofstream fapp(file, std::ios::app);
	fapp << output << std::endl;;
	fapp.close();
      }
    }

  void drawLine(const double xmin, const double ymin, const double xmax, const double ymax, const unsigned int color=kBlack, const double lineWidth=2, const unsigned int lineStyle=1)
  {
    // this function draws a line withe the chosen coordinates,
    // color and width into the active canvas
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE
    TLine *line = new TLine();
    line->SetLineWidth(lineWidth);
    line->SetLineStyle(lineStyle);
    line->SetLineColor(color);
    line->DrawLine(xmin, ymin, xmax, ymax);
  }

  int roundToInt(double value, bool roundDown=false)
  {
    // function to round a double "value"
    // to an int and return this one
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE
    // "roundDown": choose if you always want to round down

    // create output
    int outputInt=0;
    // take care of negative numbers
    if(value<0){
      std::cout << "no negative numbers allowed in roundToInt" << std::endl;
      return 0;
    }
    // get number before the digit
    for(int x=0; value>x; ++x){
      outputInt=x;
    }
    // see if the number behind the digit is > 0.5
    if((roundDown==false)&&(value-outputInt) >=0.5) return (outputInt+1);
    // return result
    return outputInt;
  }

  Double_t round(Double_t value, Int_t k)
  {
    // function to round a double "value"
    // to a double with "k" digits

    // no need for rounding if the value is 0
    if(value==0.)return value;
    Double_t result = value;
    // search window: 10^20 to 10^-20
    for(Int_t i = 20; i > -20; i--){
      Double_t exp = pow(10.,i);
      // searching for the 'k'th non-zero digit
      if(round(value/(exp*pow(10.,k-1))) != 0){
	// round to the found digit
	result = round(value/(exp))*exp;
	break;
      }
    }
    return result;
  }
  
  TString getTStringFromInt(int i)
  {
    // function to convert an int "i" to
    // a TString and return this one
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE
    char result[20];
    sprintf(result, "%i", i);
    return (TString)result;
  }

  TString getTStringFromDouble(double d, int precision=2, bool output=false)
  {
    // function to convert an double "d" to
    // a TString and return this one
    // precision: number of digits to be considered
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE
    TString conv="%.";
    conv+=getTStringFromInt(precision);
    conv+="f";
    char result[30];
    if(output) std::cout << "sprintf(result, conv, d)=sprintf(" << result << ", " << conv << ", " << d << ")" << std::endl;
    sprintf(result, conv, d);
    return (TString)result;
  }

  TH1F* divideByBinwidth(TH1F* histo, bool calculateError=true)
    {
      // function divides the #entries in every bin of the input plot "histo"
      // by its binwidth and returns the result
      // the errors are recalculated if "calculateError"=1 is choosen
      // careful: not done for underflow/overflow
      // modified quantities: NONE
      // used functions: NONE
      // used enumerators: NONE

      // create output histo
      TH1F* outputHisto = (TH1F*)histo->Clone();
      // loop bins, excluding underflow and overflow
      for(int i=1; i<= histo->GetNbinsX(); i++){
	double binwidth=(double)(histo->GetBinWidth(i));
	// take care of bins width width 0
	// (these could not be filled and
	// exist because of technical reasons only)
	if(binwidth==0)binwidth=1;
	// recalculate errors
	double error = (double)(histo->GetBinError(i)) / binwidth;
	// recalculate bin content
	outputHisto->SetBinContent(i, ( (double)(histo->GetBinContent(i)) / binwidth) );
	if(calculateError) outputHisto->SetBinError(i, error);
      }
      // return result
      return outputHisto;
    }

  void DrawDecayChLabel(TString decaychannel="", double textSize=0.04)
  {
    // Draw label for Decay Channel in upper left corner of plot

    TPaveText *decch = new TPaveText();

    decch -> AddText(decaychannel);

    decch -> SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
    decch -> SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );
    decch -> SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
    decch -> SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength()        );

    decch -> SetFillStyle(0);
    decch -> SetBorderSize(0);
    decch -> SetTextSize(textSize);
    decch -> SetTextAlign(12);
    decch -> Draw("same");
  }

  void DrawCMSLabels(bool cmsprelim=true, double luminosity=0.0, double textSize=0.04)
  {
    // Draw official labels (CMS Preliminary, luminosity and CM energy) above plot

    TPaveText *label = new TPaveText();

    label -> SetX1NDC(gStyle->GetPadLeftMargin());
    label -> SetY1NDC(1.0-gStyle->GetPadTopMargin());
    label -> SetX2NDC(1.0-gStyle->GetPadRightMargin());
    label -> SetY2NDC(1.0);
    label -> SetTextFont(42);

    if (cmsprelim)
      {
	label -> AddText(Form("CMS Preliminary, %2.1f fb^{-1} at #sqrt{s} = 7 TeV",luminosity/1000));
      }
    else
      {
	label -> AddText(Form("CMS, %2.1f fb^{-1} at #sqrt{s} = 7 TeV",luminosity/1000));
      }

    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->SetTextSize(textSize);
    label->SetTextAlign(32);
    label->Draw("same");
  }

  void DrawLabel(TString text, const double x1, const double y1, const double x2, const double y2, int centering=12, double textSize=0.04)
  {
    // function to directly draw a label into the active canvas
    // the coordinates of the window are: "x1", "y1", "x2", "y2"
    // the label will contain the text "text" and the textsize is "textSize"
    // modified quantities: the active canvas
    // used functions: NONE
    // used enumerators: NONE

    TPaveText *label = new TPaveText(x1, y1, x2, y2, "br NDC");
    label->AddText(text);
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    if(textSize!=0) label->SetTextSize(textSize);
    label->SetTextAlign(centering);
    label->Draw("same");
  }

  void scaleByFactor(TH1F*& histo, const double scaleValue)
  {
    // function to scale the TH1F plot "histo" and its errors
    // by an entered "scaleValue" correctly for each bin
    // N'=N/L, sN'=sN/L
    // modified quantities: histo
    // used functions: NONE
    // used enumerators: NONE

    // scaling factor = 0?
    if(scaleValue==0) std::cout << "WARNING: scaling histogram by 0" << std::endl;
    // loop bins (including underflow and overflow)
    for(int bin=0; bin<=histo->GetNbinsX()+1; ++bin){
      double oldError = histo->GetBinError(bin);
      // scale bin value and error
      histo->SetBinContent(bin, histo->GetBinContent(bin)/scaleValue);
      histo->SetBinError  (bin, oldError                 /scaleValue);
    }
  }

  double lumiweight(unsigned int sample, double luminosity, unsigned int kSys, const std::string decayChannel)
  {
    // this function derives the lumiweight for every standard MC
    // sample "sample" based on the theoretical cross section, the
    // number of generated events and the chosen "luminosity"
    // Furthermore, the BR correction is considered for ttbar signal
    // NOTE: enter luminosity IN / pb!!!!
    // modified quantities: NONE
    // used functions: BRcorrectionSemileptonic
    // used enumerators: samples, systematicVariation
    TString calledfunction="lumiweight(unsigned int sample=";
    calledfunction+=sample;
    calledfunction+=", double luminosity=";
    calledfunction+=luminosity;
    calledfunction+=", unsigned int kSys=";
    calledfunction+=kSys;
    calledfunction+=", const std::string decayChannel=";
    calledfunction+=decayChannel;
    calledfunction+=")";

    // function internal detail level of text output
    int verbose=0;
    // a) check if input is valid
    // sample existing?
    if(sample>=ENDOFSAMPLEENUM){
      std::cout << "ERROR: invalid sample label for lumiweight calculation, no scaling will be done" << std::endl;
      std::cout << "       ( " << calledfunction << " )" << std::endl;
      return 1.;
    }
    // valid luminosity entered?
    if(luminosity<=0){
      std::cout << "ERROR: chosen luminosity for lumiweight calculation is <= 0, no scaling will be done" << std::endl;
      std::cout << "       ( " << calledfunction << " )" << std::endl;
      return 1.;
    }
    // valid decayChannel entered?
    if(!(decayChannel.compare("muon")==0 || decayChannel.compare("electron")==0)){
      // for combined channel: no weighting but also no text output
      if(!(decayChannel.compare("combined")==0)){ 
	std::cout << "ERROR: chosen decaychannel does not correspond to electron or muon, no scaling will be done" << std::endl;
	std::cout << "       ( " << calledfunction << " )" << std::endl;
      }
      return 1.;
    }
    // b) define in/output for weight calculation
    double crossSection=1;
    double Nevents=0;
    double weight =0;

    // c) fill #events in sample and cross sections (in / pb)
    //
    // The second number always indicates number without kinematics filter for non-Pythia samples
    //
    //  MadGraph: ttbar signal + background
    if((sample==kSig)||(sample==kBkg)){
      crossSection=ttbarCrossSection; 
      // Fall11
      Nevents = 3697693;//3701947;
      // Systematic samples
      if(kSys==sysTopScaleUp  ) Nevents=3691845;//3696269;  
      if(kSys==sysTopScaleDown) Nevents=4000585;//4004587;  
      if(kSys==sysTopMatchUp  ) Nevents=4025110;//4029823; 
      if(kSys==sysTopMatchDown) Nevents=1544084;//1545688;
      if(kSys==sysTopMassUp   ) Nevents=1669928;//1671859; 
      if(kSys==sysTopMassDown ) Nevents=1618341;//1620072; 
    }
    else if((sample==kSigPow)||(sample==kBkgPow)){
      crossSection=ttbarCrossSection; 
      // Fall11
      Nevents = 16420479;//16439970;
    }
    else if((sample==kSigMca)||(sample==kBkgMca)){
      crossSection=ttbarCrossSection; 
      // Fall11
      // MC@NLO: total events is sum of weights
      //    21745196 total     weights
      //    19277872 positive  weights (88.65%)
      //     2467324 negative  weights (11.35%)
      //    16810548 effective weights
      // each |weight| = 190.41256 -> number directly from LHE
      Nevents = 3200939480 ; // 16810548 * 190.41256
    }
    // MadGraph: W->lnu+jets
    else if(sample==kWjets){
      crossSection=31314.0;
      // Fall11
      Nevents = 81171633;//81345381;
      // Systematic samples:
      if(kSys==sysVBosonScaleUp  ) Nevents= 9761537;// 9784907; 
      if(kSys==sysVBosonScaleDown) Nevents=10075524;//10092532;  
      if(kSys==sysVBosonMatchUp  ) Nevents=10438789;//10461655;
      if(kSys==sysVBosonMatchDown) Nevents= 9936639;// 9956679;  
    }
    // MadGraph: DY->ll+jets
    else if(sample==kZjets){
      crossSection=3048.0;
      // Fall11
      Nevents = 36138495;//36209629;
      // Summer11 systematic samples:
      if(kSys==sysVBosonScaleUp  ) Nevents=1589496;//1593052;
      if(kSys==sysVBosonScaleDown) Nevents=1656321;//1658995;
      if(kSys==sysVBosonMatchUp  ) Nevents=1638100;//1641367;
      if(kSys==sysVBosonMatchDown) Nevents=1612036;//1615032;
    }
    // Pythia6: QCD Mu enriched
    else if(sample==kQCD&&decayChannel.compare("muon")==0){
      // Fall11
      crossSection=296600000.*0.0002855; // generator crossSection * prefilter efficiency
      Nevents     =25080241;
    }
    // QCD e+jets channel 
    // Fall11
    // a) subsamples
    else if(sample==kQCDEM1){
      crossSection= 0.0106*236100000;  // generator prefilter efficiency * LO PYTHIA crossSection
      Nevents     = 35721833;
    }
    else if(sample==kQCDEM2){
      crossSection= 0.061*59440000;    // generator prefilter efficiency * LO PYTHIA crossSection
      Nevents     = 70392060;
    }
    else if(sample==kQCDEM3){
      crossSection= 0.159*898200;      // generator prefilter efficiency * LO PYTHIA crossSection
      Nevents     = 8150672;
    }
    else if(sample==kQCDBCE1){
      crossSection= 0.00059*236100000; // generator prefilter efficiency * LO PYTHIA crossSection
      Nevents     = 2071133;
    }
    else if(sample==kQCDBCE2){
      crossSection= 0.00242*59440000;  // generator prefilter efficiency * LO PYTHIA crossSection
      Nevents     = 2030033; 
    }
    else if(sample==kQCDBCE3){
      crossSection= 0.0105*898200;     // generator prefilter efficiency * LO PYTHIA crossSection
      Nevents     = 1082691; 
    }
    // b)combined
    else if(sample==kQCD&&decayChannel.compare("electron")==0){
      // Fall11: already added in combineMCsamples.C
      // with cross section as weight, only
      // lumi normalization is done here
      Nevents     =1.;
      crossSection=1;
    }
    // Powheg: singleTop/Antitop Samples
    // a) subsamples
    else if(sample==kSATops){
      crossSection= 1.44;
      Nevents     = 137662;//137980;
      // scale variation
      if(kSys==sysTopScaleUp  ) Nevents =153584;//153981;
      if(kSys==sysTopScaleDown) Nevents =153650;//153971;
    }
    else if(sample==kSTops){
      crossSection= 3.19;
      Nevents     = 259595;//259971;
      // scale variation
      if(kSys==sysTopScaleUp  ) Nevents =285513;//285972;
      if(kSys==sysTopScaleDown) Nevents =285182;//285602;
    }
    else if(sample==kSATopt){
      crossSection= 22.65;
      Nevents     = 1939703;//1944826;
      // scale variation
      if(kSys==sysTopScaleUp  ) Nevents =563929;//565520;
      if(kSys==sysTopScaleDown) Nevents =564169;//565454;
    }
    else if(sample==kSTopt){
      crossSection= 41.92;
      Nevents     = 3891841;//3900171;
      // scale variation
      // FIXME: t channel single top scale samples after kinematics filter
      if(kSys==sysTopScaleUp  ) Nevents = 1029915;//1032197;
      if(kSys==sysTopScaleDown) Nevents = 1037460;//1039406;
    }
    else if(sample==kSAToptW){
      crossSection= 7.87;
      Nevents     = 808200;//809984;
      // scale variation
      if(kSys==sysTopScaleUp  ) Nevents =436750;//437798;
      if(kSys==sysTopScaleDown) Nevents =437007;//437863;
    }
    else if(sample==kSToptW){
      crossSection= 7.87;
      Nevents     = 812600;//814390;
      // scale variation
      if(kSys==sysTopScaleUp  ) Nevents =436656;//437736;
      if(kSys==sysTopScaleDown) Nevents =436986;//437819;
    }
    // b) combined single top sample
    else if(sample==kSTop){
      // Fall11: already added in combineMCsamples.C
      // with cross section as weight,
      // lumi normalization is done here
      Nevents=1;
      crossSection=1.;
    }
    // Pythia6: DiBoson Samples
    // Fall11
    // a) subsamples
    else if(sample==kWW){
      crossSection= 43.0;
      Nevents     = 4225916;
    }
    else if(sample==kWZ){
      crossSection= 18.2;
      Nevents     = 4265243;
    }
    else if(sample==kZZ){
      crossSection= 5.9;
      Nevents     = 4191045;
    }
    // b) combined diboson: WW,WZ,ZZ
    else if(sample==kDiBos){
      // Fall11: already added in combineMCsamples.C
      // with cross section as weight,
      // lumi normalization is done here
      crossSection=1.;
      Nevents=1;
    }
    // Data
    else if(sample==kData){
      return 1.;
      // this leads to a weight
      // of 1.0 as data needs no weight
    }
    // unknown input
    else{
      std::cout << "ERROR: unknown configuration sample=" << sample;
      std::cout << ", decayChannel=" << decayChannel << std::endl;
      exit(0);
    }
    // d) calculate weight
    weight = effSFAB(kSys,decayChannel) * luminosity / ( Nevents / crossSection );
    if(verbose>0){
      std::cout << "sample: " << sampleLabel(sample,decayChannel) << std::endl;
      std::cout << "systematic var.: " << sysLabel(kSys) << std::endl;
    }
    // e) systematic shifts
    double weight2=weight;
    // e0) ttbar data basd reweighting
    if(kSys==sysTest){
      if(sample==kSig&&decayChannel.compare("electron")==0) weight2*=549.162/557.161;
      if(sample==kSig&&decayChannel.compare("muon"    )==0) weight2*=546.5/567.8;
      if(sample==kBkg&&decayChannel.compare("electron")==0) weight2*=3148.5/3193.74;
      if(sample==kBkg&&decayChannel.compare("muon"    )==0) weight2*=3151.23/3196.51;
    }
    if(verbose>1) std::cout << "weight before scaling: " << weight2 << std::endl;
    // e1) for ttbar->lnu: BR correction
    if((sample==kSig)||(sample==kSigPow)||(sample==kSigMca)) weight *= BRcorrectionSemileptonic;
    // e2) systematic higher/lower BG
    double scale=0;
    // (i) more/less DiBoson
    if(sample==kWW||sample==kWZ||sample==kZZ||sample==kDiBos){
      scale=0.3;
      if(kSys==sysDiBosUp  ) weight*=(1.0+scale);
      if(kSys==sysDiBosDown) weight*=(1.0-scale);
    }
    // (ii) more/less QCD
    if(sample==kQCD){
      scale=0.5;
      if(kSys==sysQCDUp  ) weight*=(1.0+scale);
      if(kSys==sysQCDDown) weight*=(1.0-scale);
    }
    // (iii) more/less single top
    if(sample==kSATops||sample==kSATopt||sample==kSAToptW||sample==kSTops||sample==kSTopt||sample==kSToptW||sample==kSTop){
      scale=0.3;
      if(kSys==sysSTopUp  ) weight*=(1.0+scale);
      if(kSys==sysSTopDown) weight*=(1.0-scale);
    }
    if(scale!=0&&verbose>0) std::cout << "possible scale factor: " << scale << std::endl;
    // printout for systematic variation
    if(verbose==1){
      std::cout << "weight";
      if(verbose>1){
	if(weight!=weight2) std::cout << "(scaled)";
	if(weight==weight2) std::cout << "(not scaled)";
      }
      std::cout << ": " << weight << std::endl;
      if(verbose>1) std::cout << "ratio: " << weight/weight2 << std::endl;
      if(weight!=weight2&&((sample==kSig)||(sample==kSigPow)||(sample==kSigMca))) std::cout << "(BR correction applied)" << std::endl;
    }
    // return result
    if(sample!=kData&&weight==1){
      std::cout << "WARNING: function lumiweight";
      std::cout << " gives result 1 for sample:" << std::endl;
      std::cout << sampleLabel(sample,decayChannel) << ", " << decayChannel << " channel" << std::endl;
      std::cout << "used xSec:       " << crossSection << std::endl;
      std::cout << "used Nevents:    " << Nevents      << std::endl;
      std::cout << "used luminosity: " << luminosity   << std::endl;
    }
    return weight;
  }

  TString getStringEntry(const TString inputTString, unsigned int entry=42, const TString seperateBy="/")
  {
    // this function devides "inputTString" using "seperateBy"
    // as seperator and returns the "entry"-th element
    // 42: default, means last element after seperator
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE
    char *path, *element;
    // enable output for debugging
    bool verbose=false;
    // save inputstring in char* path
    std::string inputString=(std::string)inputTString;
    path = new char [inputString.size()+1];
    strcpy (path, inputString.c_str());
    //  path = (char*)inputTString.Data();
    // get pointer element to elements of path using seperateBy as seperator
    element=strtok (path,(char*)(seperateBy.Data()));
    // get vector with splittet elements
    std::vector<TString> result_;
    if(verbose){
      std::cout << "input: " << inputTString << std::endl;
      std::cout << "seperator: " << seperateBy << std::endl;
      std::cout << "splitted in the following elements:" << std::endl;
    }
    while(element!=NULL){
      if(verbose) std::cout << element << std::endl;
      std::string output="";
      output=output+element;
      result_.push_back((TString)output);
      element=strtok(NULL,(char*)(seperateBy.Data()));
    }
    delete[] path;
    delete[] element;
    // default configuration: get last entry
    if(entry==42) entry = (result_.size());
    // get requestet element
    if(verbose) std::cout << "requested " << "element #" << entry << " of " << result_.size() << std::endl;
    TString result="";
    // check if element exists
    if(entry<=result_.size()) result=result_[entry-1];
    else{
      std::cout << "warning - the requested element (" << entry <<") does not exist" << std::endl;
      std::cout << "the string "  << inputTString << " with seperator " << seperateBy << " gives only " <<  result_.size() << "elements!" << std::endl;
    }
    if(verbose) std::cout << "will return: " <<  result_[entry-1] << std::endl;
    // return output
    return result_[entry-1];
  }

  TString TopFilename(unsigned int sample, unsigned int sys, const std::string decayChannel)
  {
    // this function contains the basic convention for the MC
    // .root files and returns the correct names for choosen samplesample"
    // and systematic variation "sys"
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: samples, systematicVariation

    TString fileName = "decayChannel_undefined" ;
    if      (decayChannel.compare("electron")==0) fileName ="elecDiffXSec";
    else if (decayChannel.compare("muon")==0)     fileName ="muonDiffXSec";
    else if (decayChannel.compare("combined")==0) fileName ="combinedDiffXSec"; 
    else fileName ="muonDiffXSec"; // default
    // name of data file is given directly in the .C file
    if(sample==kData) return "";
    // standard MC filenames
    if(sample==kSig    )fileName += "Sig";
    if(sample==kBkg    )fileName += "Bkg";
    if(sample==kSigPow )fileName += "SigPowheg";
    if(sample==kBkgPow )fileName += "BkgPowheg";
    if(sample==kSigMca )fileName += "SigMcatnlo";
    if(sample==kBkgMca )fileName += "BkgMcatnlo";
    if(sample==kWjets  )fileName += "Wjets";
    if(sample==kZjets  )fileName += "Zjets";
    if(sample==kDiBos  )fileName += "VV";
    if(sample==kQCD    )fileName += "QCD";
    if(sample==kSTop   )fileName += "SingleTop";
    // subsamples are all located in subfolder called "MergedFiles"
    if(sample==kWW     )fileName = "MergedFiles/"+fileName+"WW";
    if(sample==kWZ     )fileName = "MergedFiles/"+fileName+"WZ";
    if(sample==kZZ     )fileName = "MergedFiles/"+fileName+"ZZ";
    if(sample==kQCDEM1 )fileName = "MergedFiles/"+fileName+"QCDEM1";
    if(sample==kQCDEM2 )fileName = "MergedFiles/"+fileName+"QCDEM2";
    if(sample==kQCDEM3 )fileName = "MergedFiles/"+fileName+"QCDEM3";
    if(sample==kQCDBCE1)fileName = "MergedFiles/"+fileName+"QCDBCE1"; 
    if(sample==kQCDBCE2)fileName = "MergedFiles/"+fileName+"QCDBCE2";
    if(sample==kQCDBCE3)fileName = "MergedFiles/"+fileName+"QCDBCE3";  
    if(sample==kSToptW )fileName = "MergedFiles/"+fileName+"SingleTopTW";
    if(sample==kSTops  )fileName = "MergedFiles/"+fileName+"SingleTopS";
    if(sample==kSTopt  )fileName = "MergedFiles/"+fileName+"SingleTopT";
    if(sample==kSAToptW)fileName = "MergedFiles/"+fileName+"SingleAntiTopTW";
    if(sample==kSATops )fileName = "MergedFiles/"+fileName+"SingleAntiTopS";
    if(sample==kSATopt )fileName = "MergedFiles/"+fileName+"SingleAntiTopT";
    // take care of systematic variations
    // they are located in dedicated subfolders
    // JES
    bool MCTagSummer11 = 0;
    if(sys==sysJESUp  ) fileName = "JESUp/"+fileName+"JESUp";
    if(sys==sysJESDown) fileName = "JESDown/"+fileName+"JESDown";
    // JER
    if(sys==sysJERUp  ) fileName = "JERUp/"+fileName+"JERUp";
    if(sys==sysJERDown) fileName = "JERDown/"+fileName+"JERDown";
    // Shape variation
    // only for new MC and ttbar signal
    if(sample==kSig){
      if(sys==sysShapeUp  ) fileName = "MCShapeUp/"+fileName+"MCShapeVarUp";
      if(sys==sysShapeDown) fileName = "MCShapeDown/"+fileName+"MCShapeVarDown";
    }
    // Shape variation
    // only for new MC and ttbar signal
    if(sample==kSig){
      if(sys==sysPDFUp  ) fileName = "PDFUp/"+fileName+"PdfVarUp";
      if(sys==sysPDFDown) fileName = "PDFDown/"+fileName+"PdfVarDown";
    }
    // data based pt ttbar reweighting
    if(sample==kSig||sample==kBkg){
      if(sys==sysTest) fileName = "ttbarReweight/"+fileName+"SysDistortdata";                                                                                                     }  
    // Q2-Scale
    // (a) top
    if((sys==sysTopScaleUp  )&&((sample==kSig)||(sample==kBkg)||(sample==kSTop)||(sample==kSToptW)||(sample==kSTops)||(sample==kSTopt)||(sample==kSAToptW)||(sample==kSATops)||(sample==kSATopt))) fileName = "ScaleUp/"+fileName+"ScaleUp";    
    if((sys==sysTopScaleDown)&&((sample==kSig)||(sample==kBkg)||(sample==kSTop)||(sample==kSToptW)||(sample==kSTops)||(sample==kSTopt)||(sample==kSAToptW)||(sample==kSATops)||(sample==kSATopt))) fileName = "ScaleDown/"+fileName+"ScaleDown";
    // (b) V+jets
    if((sys==sysVBosonScaleUp  )&&(sample==kWjets)) fileName = "ScaleUp/"+fileName+"ScaleUp";    
    if((sys==sysVBosonScaleDown)&&(sample==kWjets)) fileName = "ScaleDown/"+fileName+"ScaleDown";
    if((sys==sysVBosonScaleUp  )&&(sample==kZjets)) {fileName = "ScaleUp/"+fileName+"ScaleUp";     MCTagSummer11=1;}
    if((sys==sysVBosonScaleDown)&&(sample==kZjets)) {fileName = "ScaleDown/"+fileName+"ScaleDown"; MCTagSummer11=1;}
    // (c) SingleTop
    if((sys==sysSingleTopScaleUp)  &&(sample==kSTop)) fileName = "ScaleUp/"+fileName+"ScaleUp";   
    if((sys==sysSingleTopScaleDown)&&(sample==kSTop)) fileName = "ScaleDown/"+fileName+"ScaleDown";
    // Matching Scale
    // (a) top   
    if((sys==sysTopMatchUp  )&&((sample==kSig)||(sample==kBkg))) fileName = "MatchUp/"+fileName+"MatchUp";   
    if((sys==sysTopMatchDown)&&((sample==kSig)||(sample==kBkg))) fileName = "MatchDown/"+fileName+"MatchDown";
    // (b) V+jets
    if((sys==sysVBosonMatchUp  )&&(sample==kWjets)) fileName = "MatchUp/"+fileName+"MatchUp";
    if((sys==sysVBosonMatchDown)&&(sample==kWjets)) fileName = "MatchDown/"+fileName+"MatchDown";
    if((sys==sysVBosonMatchUp  )&&(sample==kZjets)) {fileName = "MatchUp/"+fileName+"MatchUp";     MCTagSummer11=1;}
    if((sys==sysVBosonMatchDown)&&(sample==kZjets)) {fileName = "MatchDown/"+fileName+"MatchDown"; MCTagSummer11=1;}
    // Top Mass
    if((sys==sysTopMassUp  )&&((sample==kSig)||(sample==kBkg))) fileName = "TopMassUp/"+fileName+"TopMassUp";    
    if((sys==sysTopMassDown)&&((sample==kSig)||(sample==kBkg))) fileName = "TopMassDown/"+fileName+"TopMassDown";
    // label for MC production cycle
    fileName += (MCTagSummer11) ? "Summer11" : "Fall11";
    fileName += "PF.root";
    // return output
    return fileName;
  }

  void saveCanvas(const std::vector<TCanvas*> MyCanvas, const TString outputFolder, const TString pdfName, const bool savePdf=true, const bool saveEps=true )
  {
    // introduce function that saves every single canvas in
    // MyCanvas as ./outputFolder/CanvasTitle.eps and in addition
    // all in one rootFile: ./outputFolder/pdfName.pdf
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE

    // a) save all plots in one pdf
    if(savePdf){
      MyCanvas[0]->Print(outputFolder+pdfName+".pdf(", "pdf");
      for(unsigned int idx=1; idx<MyCanvas.size()-1; idx++){
	MyCanvas[idx]->Print(outputFolder+pdfName+".pdf", "pdf");
      }
      MyCanvas[MyCanvas.size()-1]->Print(outputFolder+pdfName+".pdf)", "pdf");
    }
    // b) save every plot as pdf with title as name
    if(saveEps){
      for(unsigned int idx=0; idx<MyCanvas.size(); idx++){
	MyCanvas[idx]->Print(outputFolder+(TString)(MyCanvas[idx]->GetTitle())+".eps");
      }
    }
  }

  std::map<unsigned int, TFile*> getStdTopAnalysisFiles( const TString inputFolder, const unsigned int systematicVariation, const TString dataFile, const std::string decayChannel, TString ttbarMC="Madgraph")
    {
      // this function returns a map containing all existing .root in "inputFolder"
      // corresponding to the choosen systematic variation "systematicVariation"
      // and the data file "dataFile" for which the direct path/name.root is needed
      // The MC .root file names correspond to the standard names as defined in the
      // function TopFilename
      // modified quantities: NONE
      // used functions: TopFilename
      // used enumerators: samples
      // ttbarMC: use "Madgraph", "Mcatnlo" or "Powheg" for the respective samples


      // open our standard analysis files and save them in a map
      std::map<unsigned int, TFile*> files_;
      // loop samples
      for(int sample = kSig; sample<=kSAToptW; sample++){
	// there are no QCD subsamples for muon channel
	if(!(sample>=kQCDEM1&&sample<=kQCDBCE3&&decayChannel.compare("muon")==0)){
	  TString fileName;
	  if(sample!=kData){
	    if(((sample==kSig)||(sample==kBkg))&&ttbarMC=="Powheg") {
	      if(sample==kSig) fileName = inputFolder+"/"+TopFilename(kSigPow, systematicVariation, decayChannel);
	      if(sample==kBkg) fileName = inputFolder+"/"+TopFilename(kBkgPow, systematicVariation, decayChannel);
	    }
	    else if(((sample==kSig)||(sample==kBkg))&&ttbarMC=="Mcatnlo") {
	      if(sample==kSig) fileName = inputFolder+"/"+TopFilename(kSigMca, systematicVariation, decayChannel);
	      if(sample==kBkg) fileName = inputFolder+"/"+TopFilename(kBkgMca, systematicVariation, decayChannel);
	    }
	    else fileName = inputFolder+"/"+TopFilename(sample, systematicVariation, decayChannel);
	  }
	  if(sample==kData) fileName = dataFile;
	  // if file exists - save in map
	  if((fileName!="no")&&(fileName!="")){
	    TFile* file = TFile::Open(fileName);
	    if(file&&!(file->IsZombie())) files_[sample]= file;
	  }
	}
      }
      return files_;
    }

  void getAllPlots( std::map<unsigned int, TFile*> files_, const std::vector<TString> plotList_,  std::map< TString, std::map <unsigned int, TH1F*> >& histo_, std::map< TString, std::map <unsigned int, TH2F*> >& histo2_, const unsigned int N1Dplots, int& Nplots, const int verbose=0, const std::string decayChannel = "unset", std::vector<TString> *vecRedundantPartOfNameInData = 0, bool SSV=false)
  {
    // this function searches for every plot listed in "plotList_" in all files listed in "files_",
    // saves all 1D histos into "histo_" and all 2D histos into "histo2_"
    // empty plots will be neglected
    // modified quantities: "histo_", "histo2_", "Nplots"
    // used functions: sampleLabel
    // used enumerators: samples
    // "N1Dplots": the #1D plots is needed as input to distinguish between 1D and 2D plots
    // "Nplots": the total # of plots is calculated
    // "verbose": set detail level of output ( 0: no output, 1: std output 2: output for debugging )
    // "redundantPartOfNameInData": delete this part in the (folder)name when loading the plots from data
    //                              (needed to handle systematic variations where foldername in data and MC is different)
    // "SSV": for all btagging plots SSV control plots are used instead of the default (CSV) plots

    // loop plots
    for(unsigned int plot=0; plot<plotList_.size(); ++plot){
      // check if plot exists in any sample
      bool existsInAnySample=false;
      // loop samples
      for(unsigned int sample=kSig; sample<=kSAToptW; ++sample){
	// check existence of sample
	if(files_.count(sample)!=0){
	  // create plot container
	  TH1* targetPlot;
	  // delete additional part of MC foldername
	  // that does not exist in data 
	  TString plotname=plotList_[plot];	  
	  if(sample==kData){
	    if (vecRedundantPartOfNameInData != 0 && vecRedundantPartOfNameInData->size() != 0){
	      std::vector<TString>::iterator iter;
	      for (iter = (vecRedundantPartOfNameInData->begin()); iter != (vecRedundantPartOfNameInData->end()); iter++){
		plotname.ReplaceAll((*iter), ""); 
	      }
	    }
	  }
	  files_[sample]->GetObject(plotList_[plot], targetPlot);
	  // Check existence of plot
	  if(targetPlot){ 
	    existsInAnySample=true;
	    // go to end of loop
	    sample=kSAToptW;
	  }
	}
      }
      // end program and draw error if plot does not exist at all
      if(!existsInAnySample){
	std::cout << "ERROR: plot "+plotList_[plot]+" does not exist in any file!" << std::endl;
	exit(1);
      }
      else{
	// otherwise: get plots from sample
	// loop samples
	for(unsigned int sample=kSig; sample<=kSAToptW; ++sample){
	  // delete additional part of MC foldername
	  // that does not exist in data 
	  TString plotname=plotList_[plot];
	  // replace CSV with SSV plots
	  if(SSV&&!plotname.Contains("bottomJet")) plotname.ReplaceAll("Tagged", "SSV");
	  if(sample==kData){
	    if (vecRedundantPartOfNameInData != 0 && vecRedundantPartOfNameInData->size() != 0){
	      std::vector<TString>::iterator iter;
	      for (iter = (vecRedundantPartOfNameInData->begin()); iter != (vecRedundantPartOfNameInData->end()); iter++){
		if (plotname.Contains((*iter))) plotname.ReplaceAll((*iter), "");
	      }
	    }
	  }
	  // check if file exists
	  // give warning if file does not exist
	  if((files_.count(sample)==0)&&(plot==0)&&(verbose>0)) std::cout << "file for " << sampleLabel(sample,decayChannel) << " does not exist- continue and neglect this sample" << std::endl;
	  if(files_.count(sample)>0){
	    // create plot container
	    TH1* targetPlot;
	    if(verbose>1){
	      std::cout << std::endl << "sample: " << sample << ", " << files_[sample]->GetName() << std::endl;
	      std::cout << "plot: " << plot << ", " << plotname << std::endl;
	    }
	    files_[sample]->GetObject(plotname, targetPlot);
	    // Check if plot exits
	    // give warning if plot does not exist
	    if((!targetPlot)&&(verbose>0)) std::cout << "can not find plot "+plotname << " in file "+(TString)(files_[sample]->GetName()) << " - continue and neglect this plot" << std::endl;
	    if(targetPlot){
	      // check if plot is empty
	      bool emptyPlot=false;
	      if((plot<N1Dplots )&&((((TH1*)(files_[sample]->Get(plotname)))->GetEntries())==0.)) emptyPlot=true;
	      if((plot>=N1Dplots)&&((((TH2*)(files_[sample]->Get(plotname)))->GetEntries())==0.)) emptyPlot=true;
	      if(emptyPlot && verbose>1) std::cout << "plot "+plotList_[plot] << " in file "+(TString)(files_[sample]->GetName()) << " is empty- continue and neglect this plot" << std::endl;
	      // to avoid problems with samples where no event is passing the selection we will drop this requirement by now
	      emptyPlot=false;
	      if(!emptyPlot){
		// save plot in corresponding map
		if(plot<N1Dplots ) histo_ [plotList_[plot]][sample] = (TH1F*)(files_[sample]->Get(plotname));
		if(plot>=N1Dplots) histo2_[plotList_[plot]][sample] = (TH2F*)(files_[sample]->Get(plotname));
		// count every existing 2D plot (every sample is counted separetly as it will be drawn into an own canvas)
		if(plot>=N1Dplots) Nplots++;
	      }
	    }
	  }
	}
	// count every existing type of 1D plots (neglect #samples here as they well be drawn into the same canvas)
	if((plot<N1Dplots)&&(histo_.count(plotList_[plot])>0)) Nplots++;
      }
    }
  }

  void scaleByLuminosity(const std::vector<TString> plotList_,  std::map< TString, std::map <unsigned int, TH1F*> >& histo_, std::map< TString, std::map <unsigned int, TH2F*> >& histo2_, const unsigned int N1Dplots, const double luminosity, const int verbose=1, const int systematicVariation=sysNo, const std::string decayChannel = "unset", TString ttbarMC="Madgraph")
  {
    // this function scales every histo in histo_ and histo2_ to the corresponding luminosity
    // Additionally the mu eff SF is applied
    // modified quantities: "histo_", "histo2_"
    // used functions: sampleLabel, effSFAB
    // used enumerators: samples
    // "N1Dplots": the #1D plots is needed as input to destinguish between 1D and 2D plots
    // "verbose": set detail level of output ( 0: no output, 1: std output 2: output for debugging )
    // "luminosity": [/pb]
    // "systematicVariation": specify systematic variation corresponding to enum systematicVariation
    // ttbarMC: use "Madgraph", "Mcatnlo" or "Powheg" for the respective samples

    // loop samples
    for(unsigned int sample=kSig; sample<=kSAToptW; ++sample) {
      unsigned int sample2=sample;
      if(ttbarMC=="Powheg"){
	if(sample==kSig) sample2=kSigPow;
	if(sample==kBkg) sample2=kBkgPow;
      }
     if(ttbarMC=="Mcatnlo"){
	if(sample==kSig) sample2=kSigMca;
	if(sample==kBkg) sample2=kBkgMca;
      }
      // loop plots
      for(unsigned int plot=0; plot<plotList_.size(); ++plot){
	// a) 1D
	// check if 1D plot exists
	if((plot<N1Dplots)&&(histo_.count(plotList_[plot])>0)&&(histo_[plotList_[plot]].count(sample)>0)){
	  if(verbose>1) std::cout << std::endl << "plot: "+plotList_[plot] << " for sample " << sampleLabel(sample2,decayChannel) << ":" << std::endl;
	  if(verbose>1) std::cout << "#events before weighting: " << histo_[plotList_[plot]][sample]->Integral(0, histo_[plotList_[plot]][sample]->GetNbinsX()+1) << std::endl;
	  // scale MC samples to same luminosity
	  double weight = lumiweight(sample2, luminosity, systematicVariation, decayChannel);
	  histo_[plotList_[plot]][sample]->Scale(weight);
	  if(verbose>1) std::cout << "weight: " << weight << std::endl;
	  if(verbose>1) std::cout << "#events after weighting: " << histo_[plotList_[plot]][sample]->Integral(0, histo_[plotList_[plot]][sample]->GetNbinsX()+1) << std::endl;
	  // apply eff. SF for MC
	  //if(sample!=kData) histo_[plotList_[plot]][sample]->Scale(effSFAB(systematicVariation,decayChannel));
	}
	// b) 2D
	// check if 2D plot exists
	if((plot>=N1Dplots)&&(histo2_.count(plotList_[plot])>0)&&(histo2_[plotList_[plot]].count(sample)>0)){
	  if(verbose>1) std::cout << std::endl << "plot: "+plotList_[plot] << " for sample " << sampleLabel(sample2,decayChannel) << ":" << std::endl;
	  if(verbose>1) std::cout << "#events before weighting: " << histo2_[plotList_[plot]][sample]->Integral(0, histo2_[plotList_[plot]][sample]->GetNbinsX()+1, 0, histo2_[plotList_[plot]][sample]->GetNbinsY()+1) << std::endl;
	  // scale MC samples to same luminosity
	  double weight = lumiweight(sample2, luminosity, systematicVariation, decayChannel);
	  histo2_[plotList_[plot]][sample]->Scale(weight);
	  if(verbose>1) std::cout << "weight: " << weight << std::endl;
	  if(verbose>1) std::cout << "#events after weighting: " << histo2_[plotList_[plot]][sample]->Integral(0, histo2_[plotList_[plot]][sample]->GetNbinsX()+1, 0, histo2_[plotList_[plot]][sample]->GetNbinsY()+1) << std::endl;
	  // apply eff. SF for MC
	  //if(sample!=kData) histo2_[plotList_[plot]][sample]->Scale(effSFAB(systematicVariation,decayChannel));
	}
      }
    }
  }

  void AddSingleTopAndDiBoson(const std::vector<TString> plotList_,  std::map< TString, std::map <unsigned int, TH1F*> >& histo_, std::map< TString, std::map <unsigned int, TH2F*> >& histo2_, const unsigned int N1Dplots, const int verbose=1, bool reCreate=false, const std::string decayChannel = "unset")
  {
    // this function creates plots for all diboson and all single
    // top samples combined if the combined SingleTop and DiBoson
    // samples do not exist as .root file.
    // Every plot in "plotList_" existing in each of the three
    // sTop/ diBoson files will be combined and saved in the histo_
    // and histo2_ maps
    // NOTE: all plots from the samples are considered to be weighted
    // to the same luminosity before
    // modified quantities: "histo_", "histo2_"
    // used functions: sampleLabel
    // used enumerators: samples
    // "N1Dplots": the #1D plots is needed as input to distinguish between 1D and 2D plots
    // "verbose": set detail level of output ( 0: no output, 1: std output 2: output for debugging )
    // "reCreate": choose if you want to create the combined plot from the
    // single plots if it alredy exists. Careful: the old one will be deleted

    if(verbose>0) std::cout << std::endl;
    // loop plots
    for(unsigned int plot=0; plot<plotList_.size(); ++plot){
      // loop STop and DiBoson
      for(unsigned int sample=kSTop; sample<=kDiBos; ++sample){
	if(verbose>1) std::cout << "plot " << plotList_[plot] << ": " << std::endl;
	// mark first plot found
	bool first=true;
	// check that combined plot does not already exist
	// if existing and reCreate==1, it will be overwritten
	bool combinedplotExists=false;
	if((plot<N1Dplots )&&(histo_ [plotList_[plot]].count(sample)>0)) combinedplotExists=true;
	if((plot>=N1Dplots)&&(histo2_[plotList_[plot]].count(sample)>0)) combinedplotExists=true;
	if(reCreate==1) combinedplotExists=false;
	if((combinedplotExists==true)&&(verbose>0)){
	  std::cout << "combined plot " << plotList_[plot];
	  std::cout << " for sample " << sampleLabel(sample,decayChannel);
	  std::cout << " already exists, will keep this one" << std::endl;
	}
	if(!combinedplotExists){
	  // loop all three subchannels
	  unsigned int firstSubChannel=kSTops;
	  if(sample==kDiBos) firstSubChannel=kWW;
	  for(unsigned int subSample=firstSubChannel; subSample<=firstSubChannel+5; ++subSample){
	    // single top has 6, diboson only 2 subsamples
	    if((sample==kDiBos&&subSample<firstSubChannel+3)||sample==kSTop){
	      // check if plot exists for the subchannel
	      bool subPlotExists=false;
	      if((plot<N1Dplots )&&(histo_ [plotList_[plot]].count(subSample)>0)) subPlotExists=true;
	      if((plot>=N1Dplots)&&(histo2_[plotList_[plot]].count(subSample)>0)) subPlotExists=true;
	      if((subPlotExists==false)&&(verbose>0)){
		std::cout << "plot " << plotList_[plot];
		std::cout << " does not exist for subSample ";
		std::cout << sampleLabel(subSample,decayChannel) << std::endl;
	      }
	      if(subPlotExists){
		// add histo
		// a) 1D
		if(plot<N1Dplots){
		  if(first ) histo_[plotList_[plot]][sample]   =  (TH1F*)histo_[plotList_[plot]][subSample]->Clone();
		  if(!first) histo_[plotList_[plot]][sample]->Add((TH1F*)histo_[plotList_[plot]][subSample]->Clone());
		}
		// b) 2D
		if(plot>=N1Dplots){
		  if(first ) histo2_[plotList_[plot]][sample]   =  (TH2F*)histo2_[plotList_[plot]][subSample]->Clone();
		  if(!first) histo2_[plotList_[plot]][sample]->Add((TH2F*)histo2_[plotList_[plot]][subSample]->Clone());
		}
		// indicate that already one plot is found
		first=0;
		// print out information
		if(verbose>0){
		  std::cout << "will add " << plotList_[plot];
		  std::cout << " from " << sampleLabel(subSample,decayChannel) << std::endl;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  void createStackPlot(const std::vector<TString> plotList_, std::map< TString, std::map <unsigned int, TH1F*> >& histo_, const unsigned int plot, const unsigned int N1Dplots, const int verbose=1, const std::string decayChannel="muon")
  {
    // this function will transform the histogram "plotList_"["plot"]
    // within "histo_" into stacked plots.
    // modified quantities: "histo_"
    // used functions: sampleLabel
    // used enumerators: samples
    // "verbose": set detail level of output ( 0: no output, 1: std output 2: output for debugging )
    // "N1Dplots": the #1D plots is needed as input to destinguish between 1D and 2D plots

    // loop samples backwards
    for(int sample=kDiBos-1; sample>=kSig; --sample){
      // check if plot exists
      if((plot<N1Dplots)&&(histo_.count(plotList_[plot])>0)&&(histo_[plotList_[plot]].count(sample)>0)){
	// check if previous plot also exists
	if(histo_[plotList_[plot]].count(sample+1)>0){
	  histo_[plotList_[plot]][sample]->Add((TH1F*)histo_[plotList_[plot]][sample+1]->Clone());
	  if(verbose>1) std::cout << "add "+plotList_[plot] << " plot for " << sampleLabel(sample,decayChannel) << " and " << sampleLabel(sample+1,decayChannel) << std::endl;
	  if(verbose>1) std::cout << "this new stacked plot contains now " <<  histo_[plotList_[plot]][sample]->Integral(0, histo_[plotList_[plot]][sample]->GetNbinsX()+1) << " events" << std::endl;
	}
      }
    }
  }

  float modulo(const float a, const float b)
  {
    // this function calculates a modulo b
    // where a and b are float
    // modified quantities: none
    // used functions: none
    // used enumerators: none
    // "a": dividend
    // "b": divisor

    // round a and b to the 3rd decimal place
    // this should prevent errors due to rounding effects
    float astar=ceil(a*1000.-0.5)/1000.;
    float bstar=ceil(b*1000.-0.5)/1000.;
    //  std::cout << std::endl << astar << " modulo ";
    //  std::cout << bstar << std::endl;
    float rest=astar-bstar;
    if(rest<0){
      std::cout << "can not compute " << astar << " modulo " << bstar << std::endl;
      std::cout << "because " << astar << " < " << bstar << std::endl;
      exit(1);
    }
    for(float n=1.; rest>0.; n++){
      rest=astar-n*bstar;
      if(rest<0){
	// allow some precision tolerance because of rounding effects
	if(rest>-0.000001) rest=0;
	else rest = astar-(n-1.)*bstar;
	break;
      }
    }
    return rest;
  }

  void reBinTH1F(TH1F& histoUnbinned, const std::vector<double> &binlowerEdges_, const int verbose=0)
  {
    // this function rebins an histogram using a variable binning
    // modified quantities: "histoUnbinned"
    // ATTENTION: histoUnbinned needs to have an initial equidistant binning
    // used functions: modulo
    // used enumerators: none
    // "histoUnbinned": plot to be binned
    // "binlowerEdges": vector containing all lower bin edges starting at xaxis.min, ending with xaxis.max
    // "verbose": set detail level of output ( 0: no output, 1: std output 2: output for debugging )

    unsigned int NinitialBins=histoUnbinned.GetNbinsX();
    double xMin=histoUnbinned.GetBinLowEdge(1);
    double xMax=histoUnbinned.GetBinLowEdge(NinitialBins+1);
    if(verbose>1){
      std::cout << "plot:            " << histoUnbinned.GetName() << std::endl;
      std::cout << "min:             " << xMin << std::endl;
      std::cout << "max:             " << xMax << std::endl;
      std::cout << "N(initial bins): " << NinitialBins << std::endl;
    }
    // check if chosen binning is valid
    // 1) N(bins)
    if(NinitialBins<=binlowerEdges_.size()){
      std::cout << "histo " << histoUnbinned.GetName() << " cannot be rebinned" << std::endl;
      std::cout << "N(chosen bins) < N(initial bins)!" << std::endl;
      exit(1);
    }
    // 2) coarseness of chosen binning
    //  double initialBinWidth=(xMax-xMin)/(double)NinitialBins;
    double initialBinWidth = (double)(roundToInt(10000.*histoUnbinned.GetBinWidth(1)))/10000.;

    if(verbose>1) std::cout << "initial binwidth: " << initialBinWidth << std::endl;
    if(binlowerEdges_.size()>1){
      for(unsigned int finalBin=1; finalBin<binlowerEdges_.size()-1; ++finalBin){
	double finalBinWidth=binlowerEdges_[finalBin]-binlowerEdges_[finalBin-1];
	if(verbose>1){
	  std::cout << "bin #" << finalBin << ": ";
	  std::cout << std::setprecision(20) << std::fixed << finalBinWidth << " modulo ";
	  std::cout << std::setprecision(20) << std::fixed << initialBinWidth;
	  std::cout << " = " << modulo(finalBinWidth, initialBinWidth) << std::endl;
	}
	if(modulo(finalBinWidth, initialBinWidth)!=0&&verbose>-1){
	  std::cout << "WARNING: histo " << histoUnbinned.GetName() << " may not be rebinned correctly!" << std::endl;
	  std::cout << "the ininital binning is to coarse for the chosen binning!" << std::endl;
	  //std::cout << "attention: probably error in modulo function," << std::endl;
	  std::cout << "bin #" << finalBin << " of the choosen binning has ";
	  std::cout << "finalBinWidth modulo initialBinWidth of:" <<std::endl;
	  std::cout << std::setprecision(20) << std::fixed << finalBinWidth;
	  std::cout << " modulo " << std::setprecision(20) << std::fixed << initialBinWidth << " = ";
	  std::cout << std::setprecision(20) << std::fixed << modulo(finalBinWidth,initialBinWidth) << std::endl;
	  //exit(1);
	}
      }
    }
    // fill vector entries in array
    // because TH1F constructor needs an array
    const unsigned int arraySize=binlowerEdges_.size();
    double *binLowerEdgeArray = new double[arraySize];
 
    for(unsigned int arrayEntry=0; arrayEntry<arraySize; ++arrayEntry){
      binLowerEdgeArray[arrayEntry]=binlowerEdges_[arrayEntry];
      if(verbose>1) std::cout << "array entry #" << arrayEntry << ": " << binLowerEdgeArray[arrayEntry] << " ";
    }
    if(verbose>1) std::cout << std::endl;
   
    if(verbose>1) std::cout << "N(bins) before rebinning: " << histoUnbinned.GetNbinsX() << std::endl;
    histoUnbinned = *(TH1F*)histoUnbinned.Rebin(arraySize-1, histoUnbinned.GetName(), binLowerEdgeArray);
    if(verbose>1){
      std::cout << "N(bins) after rebinning: " << histoUnbinned.GetNbinsX() << " (should be " << arraySize-1 << ")" << std::endl;
      std::cout << "from which only " << binlowerEdges_.size() << " bins are used" << std::endl;
      std::cout << "the rest will be ignored and exists only for technical reasons" << std::endl;
    }
    delete binLowerEdgeArray;
  }

    // ===============================================================
    // ===============================================================

    template <class histoType, class varType> 
      histoType* reBinTH1FIrregularNewBinning(histoType *histoOldBinning, const std::vector<varType> &vecBinning, TString plotname, bool rescale=1, int verbose=-1)
    {
	//  This function rebins a histogram using a variable binning
	// 
	//  (1) It is not required to have an equidistant binning.
	//  (2) Any type of ROOT-histgramme can be used, so the template 
	//      arguments should be 
	//      (a) histoT = TH1D,   TH1F,  ....
	//      (b) varT   = double, float, ....
	//  
	//  modified quantities: none
	//  used functions:      none
	//  used enumerators:    none
	//  
	//  "histoOldBinning":   plot to be re-binned
	//  "vecBinning":        vector containing all bin edges 
	//                       from xaxis.min to xaxis.max
	//  "rescale":           rescale the rebinned histogramme
	//                       (applicable for cross-section e.g.) 
        //  "verbose":           level of output 
      
        unsigned int vecIndex=0;
	
	// fill vector into array to use appropriate constructor of TH1-classes
	const varType *binArray = vecBinning.data();
	
	// create histo with new binning
	TString name = (TString)histoOldBinning->GetName();

	// suppress error massage
	// FIXME: source of error should be revised at some point
        int initialIgnoreLevel=gErrorIgnoreLevel;
	if(verbose==-1) gErrorIgnoreLevel=kError;
	histoType histoNewBinning = histoType ("histoNewBinning"+plotname+name,"histoNewBinning"+plotname+name,vecBinning.size()-1,binArray);
	gErrorIgnoreLevel=initialIgnoreLevel;

	// fill contents of histoOldBinning into histoNewBinning and rescale if applicable
	for (vecIndex = 0; vecIndex < vecBinning.size()-1; vecIndex++){
	    
	    varType lowEdge      = vecBinning[vecIndex]; 
	    //if (plotname=="topPt"&&vecIndex==0&&lowEdge==0.0) lowEdge+=10;  // adhoc fix to compensate for minimum top-Pt cut in NNLO curve
	    varType highEdge     = vecBinning[vecIndex+1];
	    varType newBinWidth  = highEdge - lowEdge;
	    varType newBinCenter = 0.5*(highEdge+lowEdge);
	    varType binSum       = 0.0;	    	  
	    
	    for (int j=0; j<histoOldBinning->GetNbinsX(); j++){
		
		varType oldBin = histoOldBinning->GetBinLowEdge(j+1);
		
		if( (oldBin>lowEdge) && (oldBin<=highEdge) ){		   
		    if (rescale) binSum+=histoOldBinning->GetBinContent(j+1) * histoOldBinning->GetBinWidth(j+1);
		    else         binSum+=histoOldBinning->GetBinContent(j+1);
		}
	    }

	    if (rescale) histoNewBinning.Fill(newBinCenter,binSum/newBinWidth);
	    else histoNewBinning.Fill(newBinCenter,binSum);
	}
	histoType * result=(histoType *)histoNewBinning.Clone();
	return result;
    }

    // ===============================================================
    // ===============================================================


  bool plotExists(std::map< TString, std::map <unsigned int, TH1F*> > histo_, const TString plotName, const unsigned int sample)
  {
    // this function checks the existence of an specific
    // entry "histo_[plotName][sample]" in the map "histo_"
    // that contains all 1D plots
    // modified quantities: none
    // used functions: none
    // used enumerators: none

    bool result=false;
    if((histo_.count(plotName)>0)&&(histo_[plotName].count(sample)>0)){
      result=true;
    }
    return result;
  }

  bool plotExists(std::map< TString, std::map <unsigned int, TH2F*> > histo_, const TString plotName, const unsigned int sample)
  {
    // this function checks the existence of an specific
    // entry "histo_[plotName][sample]" in the map "histo_"
    // that contains all 1D plots
    // modified quantities: none
    // used functions: none
    // used enumerators: none

    bool result=false;
    if((histo_.count(plotName)>0)&&(histo_[plotName].count(sample)>0)){
      result=true;
    }
    return result;
  }

  void equalReBinTH1(const int reBinFactor, std::map< TString, std::map <unsigned int, TH1F*> > histo_, const TString plotName, const unsigned int sample)
  {
    // this uses an equal rebinning (factor "reBinFactor")
    // for an specific entry "histo_[plotName][sample]"
    // in the map "histo_" that contains all 1D plots
    // modified quantities: histo_[plotName][sample]
    // used functions: plotExists
    // used enumerators: none

    if(plotExists(histo_, plotName, sample)){
      (histo_[plotName][sample])->Rebin(reBinFactor);
    }

  }
  
  std::map<TString, std::vector<double> > makeVariableBinning(bool addCrossCheckVariables=false)
    {
      // this function creates a map with the hard coded
      // bin values for variable binning
      // NOTE: it is important to quote the overflow bin
      // of the initial binning as last bin here!!!
      // otherwise dividing per binwidth function
      // might later give nan and histo wouldn't be drawn
      // modified quantities: none
      // used functions: none
      // used enumerators: none

      std::map<TString, std::vector<double> > result;
      std::vector<double> bins_;

      // pt(top)
      double topPtBins[]={0.0, 60.0, 100.0, 150.0, 200.0 , 260.0, 320.0, 400.0};  
      // PAS binning: double topPtBins[]={0., 60., 120., 200., 280., 400., 800.};
      bins_.insert( bins_.begin(), topPtBins, topPtBins + sizeof(topPtBins)/sizeof(double) );
      result["topPt"]      = bins_;
      if (addCrossCheckVariables){
	result["topPtPlus"]  = bins_;
	result["topPtMinus"] = bins_;
      }
      bins_.clear();

      // y(top)
      double topYBins[]={-2.5, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.5}; 
      // PAS binning: double topYBins[]={-5., -2.5, -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 2.5, 5.};
      bins_.insert( bins_.begin(), topYBins, topYBins + sizeof(topYBins)/sizeof(double) );
      result["topY"]      = bins_;
      if (addCrossCheckVariables){
	result["topYPlus"]  = bins_;
	result["topYMinus"] = bins_;
      }
      bins_.clear();

      // pt(ttbar)
      double ttbarPtBins[]={0.0, 20.0, 45.0, 75.0, 120.0, 190.0, 300.0};
      // PAS binning: double ttbarPtBins[]={0., 20., 60., 110., 200., 300.}; // PAS Binning
      bins_.insert( bins_.begin(), ttbarPtBins, ttbarPtBins + sizeof(ttbarPtBins)/sizeof(double) );
      result["ttbarPt"]=bins_;
      bins_.clear();

      // y(ttbar)
      // old: double ttbarYBins[]={-5., -1.3, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.3, 5.};
      double ttbarYBins[]={-2.5, -1.3, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.3, 2.5};
      // PAS binning: double ttbarYBins[]={-5., -1.3, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.3, 5.};
      bins_.insert( bins_.begin(), ttbarYBins, ttbarYBins + sizeof(ttbarYBins)/sizeof(double) );
      result["ttbarY"]=bins_;
      bins_.clear();

      // m(ttbar)
      double ttbarMassBins[]={345.0, 400.0, 470.0, 550.0, 650.0, 800.0, 1100.0, 1600.0};
      // First option: double ttbarMassBins[]={0.0, 345.0, 400.0, 470.0, 550.0, 650.0, 800.0, 1200.0};  
      // Korea:        double ttbarMassBins[]={0.0, 345.0, 400.0, 450.0, 500.0, 550.0, 600.0, 700.0, 800.0, 1200.0}; 
      // PAS binning:  double ttbarMassBins[]={0., 345., 410., 480., 580., 750., 1200.};
      bins_.insert( bins_.begin(), ttbarMassBins, ttbarMassBins + sizeof(ttbarMassBins)/sizeof(double) );
      result["ttbarMass"]=bins_;
      bins_.clear();

      // pt(lepton)
      double lepPtBins[]={30., 35., 40., 45., 50., 60., 70., 80., 100., 120., 150., 200.};
      bins_.insert( bins_.begin(), lepPtBins, lepPtBins + sizeof(lepPtBins)/sizeof(double) );
      result["lepPt"]=bins_;
      bins_.clear();

      // eta(lepton)
      double lepEtaBins[]={-2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1};
      bins_.insert( bins_.begin(), lepEtaBins, lepEtaBins + sizeof(lepEtaBins)/sizeof(double) );
      result["lepEta"]      = bins_;
      if (addCrossCheckVariables){
	result["lepEtaPlus"]  = bins_;
	result["lepEtaMinus"] = bins_;
      }
      bins_.clear();

      // pt(bquark)
      double bqPtBins[]={30., 60., 95., 140., 200., 400.};
      bins_.insert( bins_.begin(), bqPtBins, bqPtBins + sizeof(bqPtBins)/sizeof(double) );
      result["bqPt"]=bins_;
      bins_.clear();

      // eta(bquark)
      double bqEtaBins[]={-2.4, -1.5, -1, -0.5, 0., 0.5, 1., 1.5, 2.4};
      bins_.insert( bins_.begin(), bqEtaBins, bqEtaBins + sizeof(bqEtaBins)/sizeof(double) );
      result["bqEta"]=bins_;
      bins_.clear();
      return result;
    }

  void setXAxisRange(TH1* his, TString variable)
  {
    // this function restricts the drawn range of the x axis in the differential variable
    // modified quantities: xAxis of his
    // used functions: none
    // used enumerators: none
      
    // restrict axis
    if(variable.Contains("topPt"))          his->GetXaxis()->SetRangeUser(0.  , 399. ); 
    else if(variable.Contains("topY"     )) his->GetXaxis()->SetRangeUser(-2.5, 2.49 );
    else if(variable.Contains("ttbarY"   )) his->GetXaxis()->SetRangeUser(-2.5, 2.49 );
    else if(variable.Contains("ttbarMass")) his->GetXaxis()->SetRangeUser(346., 1599.);
    else if(variable.Contains("ttbarPt"  )) his->GetXaxis()->SetRangeUser(0.  , 299. );
    else if(variable.Contains("lepPt"    )) his->GetXaxis()->SetRangeUser(30  , 199. );
    else if(variable.Contains("lepEta"   )) his->GetXaxis()->SetRangeUser(-2.1, 2.09 );
    else if(variable.Contains("bqPt"     )) his->GetXaxis()->SetRangeUser(30. , 399. );
    else if(variable.Contains("bqEta"    )) his->GetXaxis()->SetRangeUser(-2.4, 2.39 );
  }
    
  template <class T>
    unsigned int positionInVector(std::vector<T> vec_, T object)
    {
      // this function returns the position
      // of element object in vector vec_
      // modified quantities: none
      // used functions: none
      // used enumerators: none

      // loop vector elements
      for(unsigned int element=0; element<vec_.size(); ++element)
	{
	  // ask for requested element
	  if(vec_[element]==object){
	    return element;
	  }
	}
      // return -1 if requested element is not found
      std::cout << "requested element " << object << " is not found in  vector" << std::endl;
      return -1;
    }

  void DivideYieldByEfficiencyAndLumi(TH1F* yield, TH1F* efficiency, double luminosity, bool includeEffError)
  {
    // this function divides yield by efficiency
    // and calculates the correct error in each bin based
    // on gaussian error propagation
    // modified quantities: yield
    // used functions: none
    // used enumerators: none

    // check if #bins are the same
    if(yield->GetNbinsX()!=efficiency->GetNbinsX()){
      std::cout << "#bins in yield and efficiency plots are not the same!" << std::endl;
      exit(1);
    }
    // loop all bins
    for(int bin=1; bin<=yield->GetNbinsX(); ++bin){
      double content  = yield     ->GetBinContent(bin);
      double eff      = efficiency->GetBinContent(bin);
      double binwidth = yield     ->GetBinWidth  (bin);
      if(binwidth==0) binwidth=1;
      double events   = binwidth*content;
      // value (events/binwidth/eff/lumi)
      double xSec = content/(eff*luminosity);
      // error
      double effError=efficiency->GetBinError(bin);
      if(!includeEffError) effError=0;
      double xSecError = 1/(binwidth*luminosity)*sqrt( ((events)/(eff*eff)) + ((events*effError)/(eff*eff))* ((events*effError)/(eff*eff)) );
      // check if efficiency is 0
      if(eff==0){
	xSec=0;
	xSecError=0;
      }
      // set value and error
      yield->SetBinContent(bin, xSec     );
      yield->SetBinError  (bin, xSecError);
    }
  }

  template <class T>
    void saveToRootFile(const TString& outputFile, const T& object, const bool& overwrite=false, const int& verbose=1, const TString& folder="")
    {
      // this function saves objects of class T
      // (such as TH1) in an output rootfile with name outputFile
      // modified quantities: outputFile
      // used functions: none
      // used enumerators: none

      bool saveObject=true;
      // check if file exist
      TFile* file = TFile::Open(outputFile, "UPDATE");
      // if not exist: create
      if(!file){
	if(verbose>1) std::cout << "file " << outputFile << " does not exist, will be created" << std::endl;
	file = new TFile(outputFile, "RECREATE");
      }
      // check if file is broken
      if(file->IsZombie()){
	std:: cout << "file " << outputFile << " is broken" << std::endl;
	// if broken: don't save object
	saveObject=false;
      }
      else{
	// if not broken: open file
	file->cd();
	// create folder if it does not exist
	if(folder!=""){
	  if(verbose>1) std::cout << "check existence of directory " << folder << std::endl;
	  // see how many subfolders are within folder
	  unsigned int Nfolders=folder.CountChar(*"/")+1;
	  if(verbose>1) std::cout << "these are " << Nfolders << " subdirectories" << std::endl;
	  // loop subdirectories
	  for(unsigned int subfolderNumber=1; subfolderNumber<=Nfolders; ++subfolderNumber){
	    TString subfolder= getStringEntry(folder, subfolderNumber, "/");
	    if(gDirectory->GetDirectory(subfolder)==0){
	      if(verbose>1) std::cout << "subfolder " << subfolder  << " not existing - will create it" << std::endl;
	      gDirectory->mkdir(subfolder);
	    }
	    else if(verbose>1) std::cout << "subfolder " << subfolder  << " is already existing" << std::endl;
	    // go to directory
	    gDirectory->cd(subfolder);
	  }
	}
	// if you don't want to overwrite, check if object exists
	int count=-1;
	if(!overwrite){
	  TString saveObjectName=(TString)object->GetTitle();
	  if(verbose>1) std::cout << "searching for object " << saveObjectName << std::endl;
	  // loop all objects in file
	  TList * list = gDirectory->GetListOfKeys();
	  while( list->At( count+1 ) != list->Last()){
	    ++count;
	    TObject *folderObject = list->At(count);
	    TString folderObjectName = (TString)folderObject->GetName();
	    if(verbose>1) std::cout << "candidate #" << count+1 << ": " << folderObjectName << std::endl;
	    // check if object you want to save is already existing
	    // by comparing the names
	    if(folderObjectName==saveObjectName){
	      if(verbose>1){
		std::cout << "already exists in file " << outputFile << std::endl;
		std::cout << "will keep the old one!" << std::endl << std::endl;
	      }
	      saveObject=false;
	      break;
	    }
	  }
	}
      }
      // save object
      if(saveObject){
	if(verbose>0){
	  std::cout << "saving object " << (TString)object->GetTitle();
	  std::cout << " to file " << outputFile;
	  if(folder!="") std::cout << " ( folder " << folder << " )";
	  std::cout << std::endl << std::endl;
	}
	// overwrite existing
	object->Write(object->GetTitle(), TObject::kOverwrite);
      }
      // close file
      file->Close();
    }

  template <class T>
    void saveToRootFileAll(const TString& outputFile, const std::vector< std::pair< T*,TString > > &objects_, const bool& overwrite=false, const int& verbose=1)
    {
      // this function saves all objects of class T in saveCanvas_
      // (such as TH1) in an output rootfile with name outputFile
      // modified quantities: outputFile
      // used functions: none
      // used enumerators: none

      bool saveObject=true;
      // check if file exist
      TFile* file = TFile::Open(outputFile, "UPDATE");
      // if not exist: create
      if(!file){
	if(verbose>1) std::cout << "file " << outputFile << " does not exist, will be created" << std::endl;
	file = new TFile(outputFile, "RECREATE");
      }
      // check if file is broken
      if(file->IsZombie()){
	std:: cout << "file " << outputFile << " is broken" << std::endl;
	// if broken: don't save object
	saveObject=false;
      }
      else{
	// if not broken: open file
	file->cd();
	// loop list of plots to be saved
	for(unsigned int idx=0; idx<objects_.size(); ++idx){  
	  //get object and folder
	  TString folder=(TString)(objects_[idx].second);
	  T* object=(T*)(objects_[idx].first);
	  // create folder if it does not exist
	  if(folder!=""){
	    if(verbose>1) std::cout << "check existence of directory " << folder << std::endl;
	    // see how many subfolders are within folder
	    unsigned int Nfolders=folder.CountChar(*"/")+1;
	    if(verbose>1) std::cout << "these are " << Nfolders << " subdirectories" << std::endl;
	    // loop subdirectories
	    for(unsigned int subfolderNumber=1; subfolderNumber<=Nfolders; ++subfolderNumber){
	      TString subfolder= getStringEntry(folder, subfolderNumber, "/");
	      if(gDirectory->GetDirectory(subfolder)==0){
		if(verbose>1) std::cout << "subfolder " << subfolder  << " not existing - will create it" << std::endl;
		gDirectory->mkdir(subfolder);
	      }
	      else if(verbose>1) std::cout << "subfolder " << subfolder  << " is already existing" << std::endl;
	      // go to directory
	      gDirectory->cd(subfolder);
	    }
	  }
	  // if you don't want to overwrite, check if object exists
	  int count=-1;
	  if(!overwrite){
	    TString saveObjectName=(TString)object->GetTitle();
	    if(verbose>1) std::cout << "searching for object " << saveObjectName << std::endl;
	    // loop all objects in file
	    TList * list = gDirectory->GetListOfKeys();
	    while( list->At( count+1 ) != list->Last()){
	      ++count;
	      TObject *folderObject = list->At(count);
	      TString folderObjectName = (TString)folderObject->GetName();
	      if(verbose>1) std::cout << "candidate #" << count+1 << ": " << folderObjectName << std::endl;
	      // check if object you want to save is already existing
	      // by comparing the names
	      if(folderObjectName==saveObjectName){
		if(verbose>1){
		  std::cout << "already exists in file " << outputFile << std::endl;
		  std::cout << "will keep the old one!" << std::endl << std::endl;
		}
		saveObject=false;
		break;
	      }
	    }
	  }
	  // save object
	  if(saveObject){
	    if(verbose>0){
	      std::cout << "saving object " << (TString)object->GetTitle();
	      std::cout << " to file " << outputFile;
	      if(folder!="") std::cout << " ( folder " << folder << " )";
	      std::cout << std::endl << std::endl;
	    }
	    // overwrite existing
	    object->Write(object->GetTitle(), TObject::kOverwrite);
	  }
	  // return to main folder
	  for(int nfolder=1; nfolder<=folder.CountChar(*"/")+1; ++nfolder){
	    gDirectory->cd("..");
	  }
	}
      }
      // close file
      file->Close();
    }

  int drawRatio(const TH1* histNumerator, TH1* histDenominator, const Double_t& ratioMin, const Double_t& ratioMax, TStyle myStyle, int verbose=0, const std::vector<double> err_=std::vector<double>(0), TString ratioLabelNominator="N_{data}", TString ratioLabelDenominator="N_{MC}", TString ratioDrawOption="p e X0", int ratioDrawColor=kBlack)
  {
    // this function draws a pad with the ratio of 'histNumerator' and 'histDenominator'
    // the range of the ratio is 'ratioMin' to 'ratioMax'
    // to the systematic variation "sys" of the enumerator "systematicVariation"
    // per default only the gaussian error of the 'histNumerator' is considered:
    // (error(bin i) = sqrt(histNumerator->GetBinContent(i))/histDenominator->GetBinContent(i))
    // if 'err_' is present and its size equals the number of bins in the histos,
    // its valus are considered as error for the ratio
    // NOTE: x Axis is transferred from histDenominator to the bottom of the canvas
    // modified quantities: none
    // used functions: none
    // used enumerators: none

    // check that histos have the same binning
    if(histNumerator->GetNbinsX()!=histDenominator->GetNbinsX()){
      std::cout << " ERROR when calling drawRatio - histos have different number of bins" << std::endl;
      return -1;
    }
    
    if(verbose>1){
      std::cout << "building ratio plot of " << histNumerator->GetName();
      std::cout << " and " << histDenominator->GetName() << std::endl;
    }
    // create ratio
    TH1F* ratio = (TH1F*)histNumerator->Clone();
    ratio->Divide(histDenominator);
    // calculate error for ratio
    // a) from err_
    if(err_.size()==(unsigned int)histNumerator->GetNbinsX()){
      if(verbose>0) std::cout << "ratio error from vector" << std::endl;
      for(int bin=1; bin<=histNumerator->GetNbinsX(); bin++){
	ratio->SetBinError(bin, err_[bin-1]);
      }
    }
    else{
      // b) default: only gaussian error of histNumerator
      if(verbose>0) std::cout << "ratio error from statistical error of " << histNumerator->GetName() << " only" << std::endl;
      for(int bin=1; bin<=histNumerator->GetNbinsX(); bin++){
	ratio->SetBinError(bin, sqrt(histNumerator->GetBinContent(bin))/histDenominator->GetBinContent(bin));
      }
    }
    // get some values from old pad
    //Int_t    logx = gPad->GetLogx();
    //Double_t left = gPad->GetLeftMargin();
    //Double_t right = gPad->GetRightMargin();
      
    Int_t    logx  = myStyle.GetOptLogx();
    Double_t left  = myStyle.GetPadLeftMargin();
    Double_t right = myStyle.GetPadRightMargin();
    
    // y:x size ratio for canvas
    double canvAsym = 4./3.;
    // ratio size of pad with plot and pad with ratio
    double ratioSize = 0.36;
    // change old pad
    gPad->SetBottomMargin(ratioSize);
    gPad->SetRightMargin(right);
    gPad->SetLeftMargin(left);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(0);
    gPad->SetFillColor(10);
    // create new pad for ratio plot
    TPad *rPad;
    rPad = new TPad("rPad","",0,0,1,ratioSize+0.001);
#ifdef DILEPTON_MACRO
    rPad->SetFillColor(10);
#else
    rPad->SetFillStyle(0);
    rPad->SetFillColor(0);
#endif
    rPad->SetBorderSize(0);
    rPad->SetBorderMode(0);
    rPad->Draw();
    // configure ratio plot
    rPad->cd();
    rPad->SetLogy(0);
    rPad->SetLogx(logx);
    rPad->SetTicky(1);
    double scaleFactor = 1./(canvAsym*ratioSize);
    ratio->SetStats(kFALSE);
    ratio->SetTitle("");
    ratio->SetMaximum(ratioMax);
    ratio->SetMinimum(ratioMin);
    ratio->SetLineWidth(2);
    // configure axis of ratio plot
    ratio->GetXaxis()->SetTitleSize(histDenominator->GetXaxis()->GetTitleSize()*scaleFactor*1.3);
    ratio->GetXaxis()->SetTitleOffset(histDenominator->GetXaxis()->GetTitleOffset()*0.9);
    ratio->GetXaxis()->SetLabelSize(histDenominator->GetXaxis()->GetLabelSize()*scaleFactor*1.4);
    ratio->GetXaxis()->SetTitle(histDenominator->GetXaxis()->GetTitle());
    ratio->GetXaxis()->SetNdivisions(histDenominator->GetNdivisions());
    ratio->GetYaxis()->CenterTitle();
    ratio->GetYaxis()->SetTitle("#frac{"+ratioLabelNominator+"}{"+ratioLabelDenominator+"}");
    ratio->GetYaxis()->SetTitleSize(histDenominator->GetYaxis()->GetTitleSize()*scaleFactor);
    ratio->GetYaxis()->SetTitleOffset(histDenominator->GetYaxis()->GetTitleOffset()/scaleFactor);
    ratio->GetYaxis()->SetLabelSize(histDenominator->GetYaxis()->GetLabelSize()*scaleFactor);
    ratio->GetYaxis()->SetLabelOffset(histDenominator->GetYaxis()->GetLabelOffset()*3.3);
    ratio->GetYaxis()->SetTickLength(0.03);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetXaxis()->SetRange(histDenominator->GetXaxis()->GetFirst(), histDenominator->GetXaxis()->GetLast());
    ratio->GetXaxis()->SetNoExponent(true);
    // delete axis of initial plot
    histDenominator->GetXaxis()->SetLabelSize(0);
    histDenominator->GetXaxis()->SetTitleSize(0);
    // draw ratio plot
    ratio->SetLineColor(ratioDrawColor);
    ratio->DrawClone(ratioDrawOption);
    ratio->SetMarkerSize(1.2);
    ratio->SetMarkerColor(ratioDrawColor);
    ratio->DrawClone(ratioDrawOption+" same");
    rPad->SetTopMargin(0.0);
    rPad->SetBottomMargin(0.15*scaleFactor);
    rPad->SetRightMargin(right);
    gPad->SetLeftMargin(left);
    gPad->RedrawAxis();
    // draw grid
    rPad->SetGrid(1,1);

    // draw a horizontal lines on a given histogram
    // a) at 1
    Double_t xmin = ratio->GetXaxis()->GetXmin();
    Double_t xmax = ratio->GetXaxis()->GetXmax();
    TString height = ""; height += 1;
    TF1 *f = new TF1("f", height, xmin, xmax);
    f->SetLineStyle(1);
    f->SetLineWidth(1);
    f->SetLineColor(kBlack);
    f->Draw("L same");
    // b) at upper end of ratio pad
    TString height2 = ""; height2 += ratioMax;
    TF1 *f2 = new TF1("f2", height2, xmin, xmax);
    f2->SetLineStyle(1);
    f2->SetLineWidth(1);
    f2->SetLineColor(kBlack);
    f2->Draw("L same");

    return 0;    
  }

  unsigned int theoryColor(TString theo="madgraph"){
    // this function returns the default color for a given theory specified by "theo"
    // modified quantities: none
    // used functions: none
    // used enumerators: none
    if(theo.Contains("mcatnlo")||theo.Contains("MC@NLO")||theo.Contains("mc@nlo")||theo.Contains("McAtNlo")||theo.Contains("Mc@Nlo")) return constMcatnloColor; 
    if(theo.Contains("powheg" )||theo.Contains("Powheg")||theo.Contains("POWHEG")||theo.Contains("PowHeg")) return constPowhegColor; 
    if(theo.Contains("nnlo"   )||theo.Contains("kidonakis")) return constNnloColor;
    if(theo.Contains("data")) return kBlack;
    return constMadgraphColor;
  }

  TString xSecLabelName(TString variable="", bool noUnit=false){
    // this function returns the x axis label for a given quantity specified by "variable"
    // if noUnit==true -> skip unit
    // modified quantities: none
    // used functions: none
    // used enumerators: none
    
    TString strUnitGeV = " #left[GeV#right]";
    if(noUnit) strUnitGeV="";

    if     (variable == "topPt"      ) return "p_{T}^{t}"+strUnitGeV;
    else if(variable == "topPtPlus"  ) return "p_{T}^{t}"+strUnitGeV;
    else if(variable == "topPtMinus" ) return "p_{T}^{#bar{t}}"+strUnitGeV;
    else if(variable == "topY"       ) return "y^{t}";
    else if(variable == "topYPlus"   ) return "y^{t}";
    else if(variable == "topYMinus"  ) return "y^{#bar{t}}";
    else if(variable == "ttbarPt"    ) return "p_{T}^{t#bar{t}}"+strUnitGeV;
    else if(variable == "ttbarY"     ) return "y^{t#bar{t}}";
    else if(variable == "ttbarMass"  ) return "m^{t#bar{t}}"+strUnitGeV;
    else if(variable == "lepPt"      ) return "p_{T}^{l}"+strUnitGeV;
    else if(variable == "lepEta"     ) return "#eta^{l}"; 
    else if(variable == "lepEtaPlus" ) return "#eta^{l^{+}}";
    else if(variable == "lepEtaMinus") return "#eta^{l^{-}}";
    else if(variable == "bqPt"       ) return "p_{T}^{b}"+strUnitGeV;
    else if(variable == "bqEta"      ) return "#eta^{b}";
    else return "Default Label for variable "+variable;
  }

  TCanvas* drawFinalResultRatio(TH1F* histNumeratorData, const Double_t& ratioMin, const Double_t& ratioMax, TStyle myStyle, int verbose=0, std::vector<TH1F*> histDenominatorTheory_=std::vector<TH1F*>(0), TCanvas* canv=0, double rangeMin=-1., double rangeMax=-1.)
  {
    // this function draws a pad with the ratio "histNumeratorData" over "histDenominatorTheoryX" 
    // for up to five specified theory curves, using "histNumeratorDataDown" and "histNumeratorDataUp"
    // the range of the ratio is 'ratioMin' to 'ratioMax'
    // to the systematic variation "sys" of the enumerator "systematicVariation"
    // the uncertainty of "histNumeratorData" only is considered for the uncertainty of the ratio:
    // ( error(bin i) = sqrt(1/histNumerator->GetBinContent(i)) / histDenominator->GetBinContent(i) )
    // NOTE: x Axis is transferred from histDenominator to the bottom of the canvas
    // modified quantities: none
    // used functions: setXAxisRange, theoryColor
    // used enumerators: none


    // some basic printout
    if(histDenominatorTheory_.size()>0&&histNumeratorData){
      if(verbose>0) std::cout << "calling drawFinalResultRatio for " << histNumeratorData->GetName() << " and " << histDenominatorTheory_.size() << " theory curves" << std::endl;
      // check that histos have the same binning
      for(unsigned int nTheory=0; nTheory<histDenominatorTheory_.size(); ++nTheory){
	if(!histDenominatorTheory_[nTheory]||histNumeratorData->GetNbinsX()!=histDenominatorTheory_[nTheory]->GetNbinsX()){
	  std::cout << "error when calling drawRatio - data and theory" << nTheory+1 << " histo have different number of bins!" << std::endl;
	  std::cout << histNumeratorData->GetNbinsX() << " (data) vs. " << histDenominatorTheory_[nTheory]->GetNbinsX() << " (theory)" << std::endl;
	  std::cout << "set this ratio to 1!" << std::endl;
	  // if not: theory is set to data -> ratio=1
	  histDenominatorTheory_[nTheory]=(TH1F*)histNumeratorData->Clone("WARNINGerrorInRatio");
	}
      }
      if(verbose>1) std::cout << "range: " << histNumeratorData->GetNbinsX() << " bins from " <<  histNumeratorData->GetBinLowEdge(histNumeratorData->GetXaxis()->GetFirst()) << " to " << histNumeratorData->GetBinLowEdge(histNumeratorData->GetXaxis()->GetLast()+1) << std::endl;

      // enter canvas
      if(!canv){
	std::cout << "canvas handed over in drawFinalResultRatio as canv is empty" << std::endl;
	return canv;
      }	       
      canv->cd();
      canv->Draw();
      
      // create ratios
      std::vector<TH1*> ratio_;
      for(unsigned int nTheory=0; nTheory<histDenominatorTheory_.size(); ++nTheory){
	if(verbose>1){
	  std::cout << "ratio for theory #" << nTheory;
	  std::cout << "  (" << histDenominatorTheory_[nTheory]->GetName();
	  std::cout << " / " << histNumeratorData->GetName() << ")" << std::endl;
	}
	// clone theory curve
	TH1F* tempRatio = (TH1F*)histDenominatorTheory_[nTheory]->Clone();
	tempRatio->Reset("icms");
	// loop bins
	for(int bin=1; bin<=histDenominatorTheory_[nTheory]->GetNbinsX(); ++bin){
	  //std::cout << "bin: " << bin << std::endl;
	  if(histNumeratorData->GetBinWidth(bin)!=0.&& histNumeratorData->GetBinContent(bin)!=0){
	    // set content
	    //std::cout << "histDenominatorTheory_[nTheory]->GetBinContent(bin): "<< histDenominatorTheory_[nTheory]->GetBinContent(bin) << std::endl;
	    //std::cout << "histNumeratorData->GetBinContent(bin): " << histNumeratorData->GetBinContent(bin) << std::endl;
	    tempRatio->SetBinContent(bin, histDenominatorTheory_[nTheory]->GetBinContent(bin) / histNumeratorData->GetBinContent(bin));
	    // set error: simple gaussian error propagation neglecting the theory error
	    tempRatio->SetBinError  (bin, ( histDenominatorTheory_[nTheory]->GetBinContent(bin) / ( histNumeratorData->GetBinContent(bin)* histNumeratorData->GetBinContent(bin)) ) * histNumeratorData->GetBinError(bin));	    
	  }
	  // set technical bins with width==0 to zero
	  else{
	    tempRatio->SetBinError  (bin, 0.5);
	    tempRatio->SetBinContent(bin, 1  );
	  }
	  if(verbose>1) {
	    std::cout << histDenominatorTheory_[nTheory]->GetName() << " bin " << bin << ":" << tempRatio->GetBinContent(bin) << "+-" << tempRatio->GetBinError(bin) << std::endl;
	    std::cout << "rel err ratio: " << tempRatio->GetBinError(bin)/tempRatio->GetBinContent(bin) << std::endl;
	    std::cout << "rel err data:  " << histNumeratorData->GetBinError(bin)/histNumeratorData->GetBinContent(bin) << std::endl;
	  }
	}
	ratio_.push_back((TH1F*)tempRatio->Clone());
      }
      // get some values from style
      Int_t    logx  = myStyle.GetOptLogx();
      Double_t left  = myStyle.GetPadLeftMargin();
      Double_t right = myStyle.GetPadRightMargin();
      setXAxisRange(histNumeratorData, (TString)histNumeratorData->GetName());
      if(rangeMin!=-1.&&rangeMax!=-1.){
	histNumeratorData->GetXaxis()->SetRange(histNumeratorData->FindBin(rangeMin),histNumeratorData->FindBin(rangeMax));
      }
      //Double_t xmin = histNumeratorData->GetXaxis()->GetXmin();
      //Double_t xmax = histNumeratorData->GetXaxis()->GetXmax();
      // y:x size ratio for canvas
      double canvAsym = 4./3.;
      // ratio size of pad with plot and pad with ratio
      double ratioSize = 0.36;
      double scaleFactor = 1./(canvAsym*ratioSize);
      // change old pad
      gPad->SetBottomMargin(ratioSize);
      gPad->SetRightMargin(right);
      gPad->SetLeftMargin(left);
      gPad->SetBorderMode(0);
      gPad->SetBorderSize(0);
      gPad->SetFillColor(10);
      // create pad to hide old axis
      TPad *whitebox = new TPad("whitebox","",0.97*left,ratioSize-0.2,1,ratioSize+0.001);
      whitebox->SetFillColor(10);
      whitebox->SetFillStyle(1001);
      whitebox->SetBorderSize(0);
      whitebox->SetBorderMode(0);
      //whitebox->Draw("");
      // create pad to hide the outer left part of the
      // rapidity x-axis without masking the 0 of the y-axis
      TPad *whitebox2 = new TPad("whitebox2","",0.176332,0.144974,0.982758,0.346332);
      whitebox2->SetFillColor(10);
      whitebox2->SetFillStyle(1001);
      whitebox2->SetBorderSize(0);
      whitebox2->SetBorderMode(0);
      //if(((TString)histNumeratorData->GetName()).Contains("Y"))	whitebox2->Draw("");
      // create new pad for ratio plot
      TPad *rPad = new TPad("rPad","",0,0,1,ratioSize+0.001);
#ifdef DILEPTON_MACRO
      rPad->SetFillColor(10);
#else
      rPad->SetFillStyle(0);
      rPad->SetFillColor(0);
#endif
      rPad->SetBorderSize(0);
      rPad->SetBorderMode(0);
      rPad->SetLogy(0);
      rPad->SetLogx(logx);
      rPad->SetTicky(1);
      // draw grid
      TPad *grid =(TPad*)rPad->Clone("grid");
      grid->SetGrid(1,1);
      grid->Draw("");
      rPad->SetGrid(1,1);
      rPad->Draw("");
      rPad->cd();
      for(unsigned int nTheory=0; nTheory<ratio_.size(); ++nTheory){
	ratio_[nTheory]->SetStats(kFALSE);
	ratio_[nTheory]->SetTitle("");
	ratio_[nTheory]->SetMaximum(ratioMax);
	ratio_[nTheory]->SetMinimum(ratioMin);
	unsigned int color = histDenominatorTheory_[nTheory]->GetMarkerColor();
	if(color==kBlack) color=theoryColor((TString)histDenominatorTheory_[nTheory]->GetName());
	ratio_[nTheory]->SetLineColor(color);
	double linewidth=histDenominatorTheory_[nTheory]->GetLineWidth();
	linewidth=2;
	ratio_[nTheory]->SetLineWidth(linewidth);
	ratio_[nTheory]->SetMarkerColor(color);
	ratio_[nTheory]->SetFillColor(color-4);
	ratio_[nTheory]->SetFillStyle(0);
	ratio_[nTheory]->SetMarkerSize(0.2);
	ratio_[nTheory]->SetLineStyle(histDenominatorTheory_[nTheory]->GetLineStyle());
	// configure axis of ratio_[nTheory] plot
	ratio_[nTheory]->GetXaxis()->SetTitleSize(histDenominatorTheory_[nTheory]->GetXaxis()->GetTitleSize()*scaleFactor*1.3);
	ratio_[nTheory]->GetXaxis()->SetTitleOffset(histDenominatorTheory_[nTheory]->GetXaxis()->GetTitleOffset()*0.9);
	ratio_[nTheory]->GetXaxis()->SetLabelSize(histDenominatorTheory_[nTheory]->GetXaxis()->GetLabelSize()*scaleFactor*1.4);
	ratio_[nTheory]->GetXaxis()->SetTitle(histDenominatorTheory_[nTheory]->GetXaxis()->GetTitle());
	ratio_[nTheory]->GetXaxis()->SetNdivisions(histDenominatorTheory_[nTheory]->GetNdivisions());
	ratio_[nTheory]->GetYaxis()->CenterTitle();
	ratio_[nTheory]->GetYaxis()->SetTitle("#frac{theory}{data}");
	ratio_[nTheory]->GetYaxis()->SetTitleSize(histDenominatorTheory_[nTheory]->GetYaxis()->GetTitleSize()*scaleFactor);
	ratio_[nTheory]->GetYaxis()->SetTitleOffset(histDenominatorTheory_[nTheory]->GetYaxis()->GetTitleOffset()/scaleFactor);
	ratio_[nTheory]->GetYaxis()->SetLabelSize(histDenominatorTheory_[nTheory]->GetYaxis()->GetLabelSize()*scaleFactor);
	ratio_[nTheory]->GetYaxis()->SetLabelOffset(histDenominatorTheory_[nTheory]->GetYaxis()->GetLabelOffset()*3.3);
	ratio_[nTheory]->GetYaxis()->SetTickLength(0.03);
	ratio_[nTheory]->GetYaxis()->SetNdivisions(505);
	ratio_[nTheory]->GetXaxis()->SetRange(histNumeratorData->GetXaxis()->GetFirst(), histNumeratorData->GetXaxis()->GetLast());
	if(nTheory==0){
	  setXAxisRange(ratio_[nTheory], (TString)ratio_[nTheory]->GetName());
	  if(rangeMin!=-1.&&rangeMax!=-1.){
	    if(verbose>0) std::cout << "manual range chosen: " << ratio_[nTheory]->FindBin(rangeMin) << " .. " << ratio_[nTheory]->FindBin(rangeMax) << std::endl;
	    ratio_[nTheory]->GetXaxis()->SetRange(ratio_[nTheory]->FindBin(rangeMin),ratio_[nTheory]->FindBin(rangeMax));
	  }
	}
	TString titleX=histNumeratorData->GetXaxis()->GetTitle();
	if(titleX!="") ratio_[nTheory]->GetXaxis()->SetTitle(titleX);
	titleX=(TString)(histNumeratorData->GetXaxis()->GetTitle());
	if(titleX!="") ratio_[nTheory]->GetXaxis()->SetTitle(titleX);
	ratio_[nTheory]->GetXaxis()->SetNoExponent(true);
	// delete axis of initial plot
	histDenominatorTheory_[nTheory]->GetXaxis()->SetLabelSize(0);
	histDenominatorTheory_[nTheory]->GetXaxis()->SetTitleSize(0);
	// draw ratio_[nTheory] plot
	if(nTheory==0){
	  // line at 1
	  TH1F* one=(TH1F*)ratio_[nTheory]->Clone();
	  one->Divide(one);
	  one->SetLineWidth(1);
	  one->SetLineColor(kBlack);
	  one->SetFillStyle(0);
	  one->DrawClone("hist");
	  // line at upper border
	  TH1F* up=(TH1F*)one->Clone("up");
	  up->SetLineWidth(3);
	  up->Scale(ratioMax);
	  up->DrawClone("hist same");
	  // errorband for data error around 1
	  TH1F* errorband=(TH1F*)one->Clone("errorband");
	  for(int bin=1; bin<=errorband->GetNbinsX(); ++bin){
	    errorband->SetBinContent(bin, 1);
	    errorband->SetBinError  (bin, histNumeratorData->GetBinError(bin)/histNumeratorData->GetBinContent(bin));
	    errorband->SetMarkerColor(kYellow);
	    errorband->SetFillColor(kYellow);
	    errorband->SetFillStyle(1001);
	  }
	  errorband->DrawClone("e2 p same");
	TLegend *leg  = new TLegend(); 
	//leg->SetX1NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.20);
	//leg->SetY1NDC(1.0 - gStyle->GetPadTopMargin()   - gStyle->GetTickLength() - 0.03*leg->GetNRows());
	//leg->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength()+0.20);
	//leg->SetY2NDC(1.0 - gStyle->GetPadTopMargin()   - gStyle->GetTickLength()-0.20);
	leg->SetX1NDC(0.2);
	leg->SetY1NDC(1.0);
	leg->SetX2NDC(0.4);
	leg->SetY2NDC(0.8);
	leg ->SetFillStyle(0);
	leg ->SetBorderSize(0);
	leg ->SetTextSize(0.1);
	leg ->SetTextAlign(12);
	leg ->AddEntry(errorband, "data stat+sys error", "F");
	leg ->DrawClone("same");
	}
	//if(nTheory==0) ratio_[nTheory]->DrawClone("hist");
	ratio_[nTheory]->DrawClone("hist same");
	//ratio_[nTheory]->DrawClone("p e X0 same");
	//ratio_[nTheory]->Print("./"+(TString)(ratio_[nTheory]->GetName())+"ratio.png");
      }
      rPad->SetTopMargin(0.0);
      rPad->SetBottomMargin(0.15*scaleFactor);
      rPad->SetRightMargin(right);
      gPad->SetLeftMargin(left);
      gPad->RedrawAxis();
      //rPad->Print("./"+(TString)(histNumeratorData->GetName())+".png");
      gPad->cd();
      rPad->Draw();
      return canv;
    }
    else{
      if(verbose>0){ 
	std::cout << "calling drawFinalResultRatio for ";
	if(histNumeratorData) std::cout << histNumeratorData->GetName();
	else  std::cout << "missing data curve";
	std::cout << " and " << histDenominatorTheory_.size() << " theory curves" << std::endl;
      }
    }
    return (TCanvas*)0;
  }

  TH1F* killEmptyBins(TH1F* histo, int verbose=0){
    // this function removes all bins with width==0 from "histo"
    // verbose: level of text output
    // modified quantities: none
    // used functions: none
    // used enumerators: none
    
    // collect nonemty bins
    std::vector<double> binEdges_;
    std::vector<double> binValues_;
    std::vector<double> blockList_;
    // loop original histogram
    if(verbose>0) std::cout << "scanning " << histo->GetName() << " for technical bins with width 0" << std::endl;
    for(int bin =1; bin<=histo->GetNbinsX(); ++bin){
       if(verbose>1) std::cout << "bin #" << bin << " ( lowedge " << histo->GetBinLowEdge(bin) << ")" << std::endl;
      if(histo->GetBinWidth(bin)!=0.){
	if(verbose>1) std::cout << "nonzero width!"<< std::endl;
	binValues_.push_back(histo->GetBinContent(bin));
	if(binEdges_.size()==0||binEdges_[binEdges_.size()-1]!=histo->GetBinLowEdge(bin)){
	  if(verbose>1) std::cout << "add "<< histo->GetBinLowEdge(bin) << std::endl;
	  binEdges_.push_back(histo->GetBinLowEdge(bin));
	}
	if(binEdges_[binEdges_.size()-1]!=histo->GetBinLowEdge(bin+1)){
	  binEdges_.push_back(histo->GetBinLowEdge(bin+1));
	  if(verbose>1) std::cout << "add "<< histo->GetBinLowEdge(bin+1) << std::endl;
	}
      }
      else{
	blockList_.push_back(bin);
	if(verbose>1) std::cout << "blocked!" << std::endl;
      }
    }
    //if(binEdges_[binEdges_.size()-1]=!histo->GetBinLowEdge(histo->GetNbinsX()+1)) binEdges_.push_back(histo->GetBinLowEdge(histo->GetNbinsX()+1));
    if(verbose>0) std::cout << "scanning 2" << std::endl;
    // create histo without empty bins
    int Nbins=binValues_.size();
    double bins[binEdges_.size()];
    if(verbose>0) std::cout << Nbins << " bins: ";
    for(unsigned int entry=0; entry<binEdges_.size(); ++entry){
      bins[entry]=binEdges_[entry];
      if(verbose>0) std::cout << bins[entry] << ", ";
    }
    if(verbose>0) std::cout << "end" << std::endl;
    if(verbose>0) std::cout << "scanning 3" << std::endl;
    TString lab="mod"+(TString)histo->GetName();
    TH1F* result = new TH1F(lab, lab, Nbins, bins);
    // fill new histo
    int resultbin=1;
    for(int bin =1; bin<=histo->GetNbinsX(); ++bin){
      // check if bin has to be refilled
      bool fill=true;
      for(unsigned int blocked=0; blocked<blockList_.size(); ++blocked){
	if(blockList_[blocked]==bin) fill=false;
      }
      if(fill){
	result->SetBinContent(resultbin,histo->GetBinContent(bin));
	result->SetBinError  (resultbin,histo->GetBinError  (bin));
	resultbin++;
      }
    }
    // first bin
    result->SetBinContent(0,histo->GetBinContent(0));
    result->SetBinError  (0,histo->GetBinError  (0));
    // last bin
    result->SetBinContent(result->GetNbinsX()+1,histo->GetBinContent(histo->GetNbinsX()+1));
    result->SetBinError  (result->GetNbinsX()+1,histo->GetBinError  (histo->GetNbinsX()+1));
    // print out final plot
    if(verbose>1){
      for(int bin =0; bin<=result->GetNbinsX()+1; ++bin){
	std::cout << "bin #" << bin << ": " <<  result->GetBinContent(bin) << " +- " << result->GetBinError(bin);
	std::cout << "(" << result->GetBinLowEdge(bin) << ")" << std::endl;
      }
    }
    // return cleaned plot
    return result;
  }

  double getInclusiveXSec(TH1* hist, int verbose=0)
  {
    // this function integrates a given
    // histogramm 'hist' by summing up the product
    // of binwidth and bincontent for all bins.
    // NOTE: binwidth==0 is treated as binwidth==1
    // modified quantities: none
    // used functions: none
    // used enumerators: none
    double value=0;
    if(verbose>0) std::cout << "histo: " << hist->GetTitle() << std::endl;
    for(int bin=0; bin<=hist->GetNbinsX()+1; ++bin){

      double binWidth=hist->GetBinWidth(bin);
      if(binWidth==0.||bin==0||bin==hist->GetNbinsX()+1) binWidth=1.;
      double thisBin=binWidth*hist->GetBinContent(bin);
      value+=thisBin;
      if(verbose>1){ 
	std::cout << "xSec bin" << bin << ": " << thisBin << std::endl;
	std::cout << "( with lower edge" << hist->GetBinLowEdge(bin) << ")" << std::endl;
      }
    }
    if(verbose>0) std::cout << "integrated xSec: " << value << std::endl;
    return value;
  }


  void whipEmptyBinsAway(TGraphAsymmErrors* hist, int verbose=0)
  {
    // this function shifts all unset points
    // (at x,y=0,0 with error 0,0)
    // out of sight to 0, -1000
    // modified quantities: none
    // used functions: none
    // used enumerators: none
    if(verbose>1) std::cout << "checking " << hist->GetName() << " for empty bins" << std::endl;
    // loop points
    for(int bin=0; bin<=hist->GetN()-1; ++bin){
      // search for unset points
      if( hist->GetX()[bin]==0 &&
	  hist->GetY()[bin]==0 &&
	  hist->GetErrorXhigh(bin)==0 &&
	  hist->GetErrorXlow (bin)==0 &&
	  hist->GetErrorYhigh(bin)==0 &&
	  hist->GetErrorYlow (bin)==0 ){
	hist->SetPoint(bin, 0, -1000);
	if(verbose>1) std::cout << "bin #" << bin << " set to 0,-1000" << std::endl;
      }
    }
  }

  TH1F* getTheoryPrediction(TString plotName="", TString fileName="")
  {

    // check if input is valid
    if(plotName==""||fileName==""){
      std::cout << "ERROR: no plotName or filename chosen in function getTheoryPrediction" << std::endl;
      exit(0);
    }
    // load file
    TFile* file = new (TFile)(fileName);
    // check if file is valid
    if(file->IsZombie()){
      std::cout << "ERROR: input file " << fileName << " in function getTheoryPrediction is zombie" << std::endl;
      exit(0);
    }
    // get plot
    TH1* targetPlot;
    file->GetObject(plotName, targetPlot);
    // check if plot exists
    if(!targetPlot){
      std::cout << "ERROR in function getTheoryPrediction: plot " << plotName;
      std::cout << " does not exist in file " << fileName << std::endl;
      exit(0);
    }
    gROOT->cd();
    TH1F* result = (TH1F*)(targetPlot->Clone());
    file->Close();
    delete file;
    return result;
  }

  double signum(double value){
    // this function is a
    // simple signum function
    // modified quantities: none
    // used functions: none
    // used enumerators: none

    if(value>0) return  1.0;
    if(value<0) return -1.0;
    return 0;
  }

  void closeStdTopAnalysisFiles(std::map<unsigned int, TFile*> files_)
    {
      // this function closes and deletes the files saved in "files_"
      // modified quantities: files_
      // used functions: NONE
      // used enumerators: NONE

      // close and delete input files
      for(std::map<unsigned int, TFile*>::const_iterator file=files_.begin(); file!=files_.end(); ++file){
	if(file->second){
	  file->second->Close();	  delete file->second;
	}
      }
    }

  // ===========================================================================================
  // ===========================================================================================

  void makeUncertaintyBands(std::map< TString, std::map <unsigned int, TH1F*> >& histo_,
			    std::map< TString, TH1F* >& histoErrorBand_,
			    std::vector<TString>& plotList_,
			    unsigned int& Nplots)
  {
    // this function generates set of errorbands ("histoErrorBand_") for all samples 
    // and histogrammes in "histo_", using the "Nplots" plots listed in "plotList_" 
    // taking into account luminosity and cross section error
    // 
    // modified quantities: histoErrorBand_
    // used functions:      none
    // used enumerators:    samples

    for(unsigned int plot=0; plot<Nplots; plot++){
      
      TString plotName = plotList_[plot];
      
      // Initialize and reset histograms
      
      TH1F* histoSumRef       = (TH1F*)histo_[plotName][kSig]->Clone();
      TH1F* histoSumTTbarOnly = (TH1F*)histo_[plotName][kSig]->Clone();
      
      histoSumRef       -> Reset("ICESM");
      histoSumTTbarOnly -> Reset("ICESM");
      
      // Integral over all samples before accessing the differences
      
      for(unsigned int sample=kSig; sample<kData; ++sample){

	if((plot<Nplots)&&(plotExists(histo_, plotName, sample))){

	  histoSumRef -> Add(histo_[plotName][sample]);

	  if (sample == kSig || sample == kBkg) 
	    histoSumTTbarOnly -> Add(histo_[plotName][sample]);
	}
      }
      
      // Compare summed histograms, symmetrize deviations and store relative deviation to error histogram after adding constant contributions
      
      double relXSecError = ttbarCrossSectionError/ttbarCrossSection;
      
      for (int bin = 0; bin < histoSumRef->GetNbinsX(); bin++){

	double xSecError       = (histoSumTTbarOnly->GetBinContent(bin+1))*relXSecError;
	double luminosityError = (histoSumRef->GetBinContent(bin+1))*globalLumiUncertainty;
	
	double totalError      = sqrt(xSecError*xSecError+luminosityError*luminosityError);
	
	histoSumRef->SetBinError(bin+1,totalError);
      }
      
      histoErrorBand_[plotName] = (TH1F*)histoSumRef->Clone();
    }
  }

  // ===========================================================================================
  // ===========================================================================================
  
  void addCanvas(std::vector<TCanvas*> &plotCanvas_){
    // this function adds a new canvas to the vector of canvas "plotCanvas_"
    // modified quantities: plotCanvas_
    // used functions: NONE
    // used enumerators: NONE
    char canvname[10];
    int number =plotCanvas_.size()+1;
    sprintf(canvname,"canv%i",number);    
    plotCanvas_.push_back( new TCanvas( canvname, canvname, 600, 600) );
    //canvasStyle(*plotCanvas_[number-1]);
  }

  TH1F* useFittedFunctions(TH1F* input, TString model="", TString plotname="", int verbose=0){
    // this function modifys "input" in a way that wiggly parts 
    // with large uncertainties are replaced by fitted functions
    // modified quantities: input
    // used functions: NONE
    // used enumerators: NONE
    // input: fine binned theory curve
    // model: indicates theory (madgraph, powheg or mcatnalo)
    // verbose: level of output
    int initialIgnoreLevel=gErrorIgnoreLevel;
    // check chosen PS
    TString PS="";
    if(!plotname.Contains("PhaseSpace")) PS+="full";
    else if(plotname.Contains("Hadron")) PS+="had";
    else if(plotname.Contains("Parton")) PS+="part";
    else std::cout << "error in useFittedFunctions, can not get PS from plotname " << plotname << std::endl;
    // if large sample:
    double largeSampleSF=1;
    bool largeSample=false;
    if(model=="madgraph"&&((PS=="full"&&input->GetEntries()>10000000)||(PS!="full"&&input->GetEntries()>3000000))){
      largeSampleSF=17.;
      largeSample=true;
    }
    if(verbose==0) gErrorIgnoreLevel=kFatal;
    if(verbose>1) std::cout << plotname << "(" << model << ")" << std::endl;
    // clone histo for output
    TH1F* result=(TH1F*)input->Clone();
    // define fit function A
    double fitLowEdge=0;
    double fitHighEdge=0;
    TString def="";
    double a=0.;
    double b=0.;
    double c=0.;
    double d=0.;
    double e=0.;
    TString addOpt="";
    // define fit function B
    double fitLowEdgeB=0;
    double fitHighEdgeB=0;
    TString defB="";
    double aB=0.;
    double bB=0.;
    double cB=0.;
    double dB=0.;    double eB=0.;
    TString addOptB="";
    // define fit function C
    double fitLowEdgeC=0;
    double fitHighEdgeC=0;
    TString defC="";
    double aC=0.;
    double bC=0.;
    double cC=0.;
    double dC=0.;
    double eC=0.;
    TString addOptC="";
    // configure parameter for all specific plots
    if(model=="mcatnlo"){
      if(plotname.Contains("TTbarM")){
	//tail:
	fitLowEdge=430.;
	fitHighEdge=1600.;
	def="[0]*exp([1]*x)+[2]";
	a=0.11;
	b=-0.0074;
	c=0.00001;
	addOpt="LL";
	// start:
	fitLowEdgeB=345.;
	fitHighEdgeB=400.;
	defB="TMath::GammaDist(x,[0],[1],[2])";
	aB=1.53;
	bB=345.;
	cB=96.2;
	addOptB="LL";
	// very early start:
	fitLowEdgeC=317.;
	fitHighEdgeC=336.;
	defC="[0]*exp([1]*x)+[2]";
	aC=0.00000373757;
	bC=0.0125775;
	cC=-0.000194989;
	addOptC="LL";
      }
      else if(plotname.Contains("TTbarPt")){
	//tail:
	fitLowEdge=50.;
	fitHighEdge=300.;
	def="[0]*exp([1]*x)+[2]*x+[3]*x*x+[4]";
	a= 0.0214186;
	b=-0.0274519;
	c=-0.0000135818;
	d= 0.0000000213453;
	e= 0.00227380;
	addOpt="LL";
      }
     else if(plotname.Contains("TTbarY")){
	fitLowEdge=0.5;
	fitHighEdge=-0.75;
	def="TMath::Exp(x*x*[0]+x*x*x*x*[1]+x*x*x*x*x*x*[2])*[3]";
	a=-0.52;
	b=-0.047;
	c=-0.094;
	d=812149;
	addOpt="LL";
     }
     else if(plotname.Contains("TopY")){
       //tail:
       fitLowEdge=-1.;
       fitHighEdge=1.;
       def="[0]*TMath::Gaus(x,0,[1])";
       a=0.00402075;
       b=1.14871;
       addOpt="LL";
     }
     //else if(plotname.Contains("LepPt")){
     //  fitLowEdge=45.;
     //  fitHighEdge=210.;
     //  def="TMath::Exp(x*[0]+x*x*[1])*[2]";
     //  a=-0.0414389;
     //  b=0.0000334359;
     //  c=0.133893;
     //}
    }
    else if(model=="madgraph"){
      if(plotname.Contains("ttbarMass")){
	//tail:
	fitLowEdge=470.0;
	fitHighEdge=720.0;
	def="[0]*exp([1]*x)+[2]";
	a= (PS=="full" ? 164748.0 : 28165.)*largeSampleSF;
	b= (PS=="full" ? -0.00810181 : -0.00756);
	c= (PS=="full" ? 3.59495 : 2.77)*largeSampleSF;
	addOpt="LL";
        fitLowEdgeB=680.0;
        fitHighEdgeB=920.0;
        defB="[0]*exp([1]*x)+[2]";
        aB= (PS=="full" ? 164748.0 : 28165.)*largeSampleSF;
        bB= (PS=="full" ? -0.00810181 : -0.00756);
        cB= (PS=="full" ? 3.59495 : 2.77)*largeSampleSF;
        addOptB="LL";  
	//end of tail:
	fitLowEdgeC=880.0;
	fitHighEdgeC=1600.0;
	defC="[0]*exp([1]*x)+[2]";
	aC= 54699.0*largeSampleSF;  
	bC=-0.0068;
	cC= 0.979;
	addOptC="LL";
	// start:
	//if(!largeSample){
	//  fitLowEdgeC=345.;
	//  fitHighEdgeC=400.;
	//  defC="[3]*TMath::GammaDist(x,[0],[1],[2])";
	//  aC=(PS=="full" ? 1.50 : 1.56 );
	//  bC=(PS=="full" ? 345. : 344.8);
	//  cC=(PS=="full" ? 79.3 : 80.8 );
	//  dC=(PS=="full" ? 1010990 : 1000000.);
	//  addOptC="LL";
	//}
      }
/*       else if(plotname.Contains("ttbarPt")){ */
/* 	//tail: */
/* 	fitLowEdge=50.; */
/* 	fitHighEdge=300.; */
/* 	def="[0]*exp([1]*x)+[2]*x+[3]*x*x+[4]"; */
/* 	a= 5180.*largeSampleSF; */
/* 	b=-0.03; */
/* 	c=-4.48; */
/* 	d= 0.007; */
/* 	e= 745.5; */
/* 	addOpt="LL"; */
/*       } */
/*       else if(plotname.Contains("topY")){ */
/* 	//tail: */
/* 	if(!largeSample){ */
/* 	  fitLowEdge=-1.; */
/* 	  fitHighEdge=1.; */
/* 	  def="[0]*TMath::Gaus(x,0,[1])"; */
/* 	  a=19845.2; */
/* 	  b=1.09186; */
/* 	  addOpt="LL"; */
/* 	} */
/*       } */
      else if(plotname.Contains("TTbarY")){
	fitLowEdge=-0.5;
	fitHighEdge=0.75;
	def="TMath::Exp(x*x*[0]+x*x*x*x*[1]+x*x*x*x*x*x*[2])*[3]";
	a=-0.52;
	b=-0.047;
	c=-0.094;
	d=812149;
	addOpt="LL";
      }
      else if(plotname.Contains("lepPt")){ 
/* 	// head */
/* 	fitLowEdge = 30.; */
/* 	fitHighEdge= 58.; */
/* 	def="TMath::Exp(x*[0]+x*x*[1]+x*x*x*[2])*[3]"; */
/* 	a=(largeSample ?  5.53215e-02 :  4.26320e-02); */
/* 	b=(largeSample ? -1.36559e-03 : -1.11187e-03); */
/* 	c=(largeSample ?  7.65236e-06 :  6.03412e-06); */
/* 	d=(largeSample ?  6.78192e+04 :  4.84434e+03); */
 	// tail 
 	fitLowEdgeB =50.; 
 	fitHighEdgeB=200.; 
 	defB="TMath::Exp(x*[0]+x*x*[1]+x*x*x*[2])*[3]"; 
 	aB=(largeSample ? -2.03019e-02 : -1.39410e-02); 
 	bB=(largeSample ? -9.19906e-05 : -1.49428e-04); 
 	cB=(largeSample ?  2.63030e-07 :  4.11490e-07); 
 	dB=(largeSample ?  3.18800e+05 :  1.51488e+04); 
      } 
      else if(plotname.Contains("bqPt")){
 	// head
 	fitLowEdge =30.0;
 	fitHighEdge=93.0;
 	def="TMath::Exp(x*[0]+x*x*[1]+x*x*x*[2])*[3]";
 	a= 6.96386e-02;
 	b=-1.01480e-03;
 	c= 3.86642e-06;
 	d= 1.63953e+03;
 	// tail
 	fitLowEdgeB = 92.0;
 	fitHighEdgeB=410.0;
 	defB="TMath::Exp(x*[0]+x*x*[1]+x*x*x*[2])*[3]";
 	aB=-1.87272e-02;
 	bB=-3.43912e-05;
 	cB= 7.16568e-08;
 	dB= 2.67318e+04;
      }
/*       else if(plotname.Contains("bqPt")&&plotname.Contains("Hadron")){ */
/* 	// tail */
/* 	fitLowEdge =100.0; */
/* 	fitHighEdge=400.0; */
/* 	def="TMath::Exp(x*[0]+x*x*[1]+x*x*x*[2])*[3]"; */
/* 	a=-2.08896e-02; */
/* 	b=-2.63966e-05; */
/* 	c= 6.52234e-08; */
/* 	d= 5.32653e+05; */
/*       } */
    }
    else if(model=="powheg"){
      if(plotname.Contains("ttbarMass")){
      fitLowEdge=480.;
      fitHighEdge=1600.;
      def="[0]*exp([1]*x)+[2]";
      a=40556.;
      b=-0.007;
      c=3.1;
      addOpt="LL";
      // start:
      fitLowEdgeB=345.;
      fitHighEdgeB=500.;
      defB="[3]*TMath::GammaDist(x,[0],[1],[2])";
      aB=1.58;
      bB=344.;
      cB=86.6;
      dB=234473.;
      addOptB="LL";
      // very early start:
      fitLowEdgeC=295.;
      fitHighEdgeC=337.;
      defC="[0]*exp([1]*x)+[2]";
      aC=0.000000000530501;
      bC=0.075;
      cC=3.81;
      addOptC="LL";
      }
      else if(plotname.Contains("ttbarPt")){
	//tail:
	fitLowEdge=50.;
	fitHighEdge=300.;
	def="[0]*exp([1]*x)+[2]*x+[3]*x*x+[4]";
	a= 5180.;
	b=-0.03;
	c=-4.48;
	d= 0.007;
	e= 745.5;
	addOpt="LL";
      }
     else if(plotname.Contains("ttbarY")){
	fitLowEdge=-1.5;
	fitHighEdge=1.5;
	def="TMath::Exp(x*x*[0]+x*x*x*x*[1]+x*x*x*x*x*x*[2])*[3]";
	a=-0.598570;
	b=-0.120246;
	c=-0.00770785;
	d=12987.6;
	addOpt="LL";
      }
     //else if(plotname.Contains("lepPt")){
     //  fitLowEdge=45.;
     //  fitHighEdge=210.;
     //  def="TMath::Exp(x*[0]+x*x*[1])*[2]";
     //  a=-0.0414389;
     //  b=0.0000334359;
     //  c=0.133893;
     //}
    }
    // use fit function A
    if(def!=""&&fitLowEdge!=fitHighEdge){
      if(verbose>0) std::cout << "will perform fit" << std::endl;
      // empty bins within fit range for output histo
      bool lastBin=false;
      for(int bin=1; bin<=result->GetNbinsX()+1; ++bin){
	bool thisBin=false;
	double binX=result->GetBinCenter(bin);
	if(binX>=fitLowEdge&&binX<=fitHighEdge){ 
	  result->SetBinContent(bin, 0.);
	  thisBin=true;
	}
	// adapt fit range to binning
	// low
	if(lastBin==false&&thisBin==true) fitLowEdge=result->GetBinLowEdge(bin);
	// high
	if(lastBin==true&&thisBin==false) fitHighEdge=result->GetBinLowEdge(bin);	
	lastBin=thisBin;
      }
      if(verbose>0){
	std::cout << "fitting " << def << " in: " << fitLowEdge << ".." << fitHighEdge << std::endl;
	if(a!=0) std::cout << "start value [0]:" << a << std::endl;
	if(b!=0) std::cout << "start value [1]:" << b << std::endl;
	if(c!=0) std::cout << "start value [2]:" << c << std::endl;
	if(d!=0) std::cout << "start value [3]:" << d << std::endl;
	if(e!=0) std::cout << "start value [4]:" << e << std::endl;
      }
      // perform fit
      TF1* function=new TF1("function",def,fitLowEdge, fitHighEdge);
      if(a!=0.){
	function->SetParameter(0,a);
	//if(def.Contains("GammaDist")) function->SetParLimits(0,0.6*a,1.4*a);
      }
      if(b!=0.){
	function->SetParameter(1,b);
	if(def.Contains("GammaDist")) function->SetParLimits(1,0.85*b,1.15*b);
      }
      if(c!=0.){
	function->SetParameter(2,c);
	//if(def.Contains("GammaDist")) function->SetParLimits(2,0.6*c,1.4*c);
      }
      if(d!=0.){
	function->SetParameter(3,d);
	//if(def.Contains("GammaDist")) function->SetParLimits(3,0.6*d,1.4*d);
      }
      if(e!=0.){
	function->SetParameter(4,e);
	//if(def.Contains("GammaDist")) function->SetParLimits(4,0.6*e,1.4*e);
      }
      TString option="Q";
      if(verbose>0)option="";
      if(verbose>1)option="V";
      option+="N";
      option+=addOpt;
      option+="R";
      //option+="M";
      input->Fit(function,option,"same",fitLowEdge, fitHighEdge);
      // add part from fit
      result->Add(function);
    }
    // use fit function B
    if(defB!=""&&fitLowEdgeB!=fitHighEdgeB){
      if(verbose>0) std::cout << "will perform fit" << std::endl;
      // empty bins within fit range for output histo
      bool lastBin=false;
      for(int bin=0; bin<=result->GetNbinsX()+1; ++bin){
	bool thisBin=false;
	double binX=result->GetBinCenter(bin);
	if(binX>=fitLowEdgeB&&binX<=fitHighEdgeB){ 
	  result->SetBinContent(bin, 0.);
	  thisBin=true;
	}
	// adapt fit range to binning
	// low
	if(lastBin==false&&thisBin==true) fitLowEdgeB=result->GetBinLowEdge(bin);
	// high
	if(lastBin==true&&thisBin==false) fitHighEdgeB=result->GetBinLowEdge(bin);	
	lastBin=thisBin;
      }
      if(verbose>0){
	std::cout << "fitting " << defB << " in: " << fitLowEdgeB << ".." << fitHighEdgeB << std::endl;
	if(aB!=0) std::cout << "start value [0]:" << aB << std::endl;
	if(bB!=0) std::cout << "start value [1]:" << bB << std::endl;
	if(cB!=0) std::cout << "start value [2]:" << cB << std::endl;
	if(dB!=0) std::cout << "start value [3]:" << dB << std::endl;
	if(eB!=0) std::cout << "start value [4]:" << eB << std::endl;
      }
      // perform fit
      TF1* functionB=new TF1("functionB",defB,fitLowEdgeB, fitHighEdgeB);
      if(aB!=0.){
	functionB->SetParameter(0,aB);
	if(defB.Contains("GammaDist")){
	  if(aB<0.)functionB->SetParLimits(0,1.4*aB,0.6*aB);
	  if(aB>0.)functionB->SetParLimits(0,0.6*aB,1.4*aB);	  
	}
      }
      if(bB!=0.){
	functionB->SetParameter(1,bB);
	if(defB.Contains("GammaDist")){
	  if(bB<0.) functionB->SetParLimits(1,1.15*bB,0.85*bB);
	  if(bB>0.) functionB->SetParLimits(1,0.85*bB,1.15*bB);
	}
      }
      if(cB!=0.){
	functionB->SetParameter(2,cB);
	if(defB.Contains("GammaDist")){
	  if(cB<0.) functionB->SetParLimits(2,1.4*cB,0.6*cB);
	  if(cB>0.) functionB->SetParLimits(2,0.6*cB,1.4*cB);
	}
      }
      if(dB!=0.){
	functionB->SetParameter(3,dB);
	if(defB.Contains("GammaDist")){
	  if(dB<0.) functionB->SetParLimits(3,1.4*dB,0.6*dB);
	  if(dB>0.) functionB->SetParLimits(3,0.6*dB,1.4*dB);
	}
      }
      if(eB!=0.){
	functionB->SetParameter(4,eB);
	if(defB.Contains("GammaDist")){
	  if(eB<0.) functionB->SetParLimits(4,1.4*eB,0.6*eB);
	  if(eB>0.) functionB->SetParLimits(4,0.6*eB,1.4*eB);
	}
      }
      TString option="Q";
      if(verbose>0)option="";
      if(verbose>1)option="V";
      option+="N";
      option+=addOptB;
      option+="R";
      //option+="M";
      input->Fit(functionB,option,"same",fitLowEdgeB, fitHighEdgeB);
      // add part from fit
      result->Add(functionB);
    }
   // use fit function C
    if(defC!=""&&fitLowEdgeC!=fitHighEdgeC){
      if(verbose>0) std::cout << "will perform fit" << std::endl;
      // empty bins within fit range for output histo
      bool lastBin=false;
      for(int bin=1; bin<=result->GetNbinsX()+1; ++bin){
	bool thisBin=false;
	double binX=result->GetBinCenter(bin);
	if(binX>=fitLowEdgeC&&binX<=fitHighEdgeC){ 
	  result->SetBinContent(bin, 0.);
	  thisBin=true;
	}
	// adapt fit range to binning
	// low
	if(lastBin==false&&thisBin==true) fitLowEdgeC=result->GetBinLowEdge(bin);
	// high
	if(lastBin==true&&thisBin==false) fitHighEdgeC=result->GetBinLowEdge(bin);	
	lastBin=thisBin;
      }
      if(verbose>0){
	std::cout << "fitting " << defC << " in: " << fitLowEdgeC << ".." << fitHighEdgeC << std::endl;
	if(aC!=0) std::cout << "start value [0]:" << aC << std::endl;
	if(bC!=0) std::cout << "start value [1]:" << bC << std::endl;
	if(cC!=0) std::cout << "start value [2]:" << cC << std::endl;
	if(dC!=0) std::cout << "start value [3]:" << dC << std::endl;
	if(eC!=0) std::cout << "start value [4]:" << eC << std::endl;
      }
      // perform fit
      TF1* functionC=new TF1("functionC",defC,fitLowEdgeC, fitHighEdgeC);
      if(aC!=0.){
	functionC->SetParameter(0,aC);
	if(defC.Contains("GammaDist")) functionC->SetParLimits(0,0.6*aC,1.4*aC);
      }
      if(bC!=0.){
	functionC->SetParameter(1,bC);
	if(defC.Contains("GammaDist")) functionC->SetParLimits(1,0.85*bC,1.15*bC);
      }
      if(cC!=0.){
	functionC->SetParameter(2,cC);
	if(defC.Contains("GammaDist")) functionC->SetParLimits(2,0.6*cC,1.4*cC);
      }
      if(dC!=0.){
	functionC->SetParameter(3,dC);
	if(defC.Contains("GammaDist")) functionC->SetParLimits(3,0.6*dC,1.4*dC);
      }
      if(eC!=0.){
	functionC->SetParameter(4,eC);
	if(defC.Contains("GammaDist")) functionC->SetParLimits(4,0.6*eC,1.4*eC);
      }
      TString option="Q";
      if(verbose>0)option="";
      if(verbose>1)option="V";
      option+="N";
      option+=addOptC;
      option+="R";
      //option+="M";
      input->Fit(functionC,option,"same",fitLowEdgeC, fitHighEdgeC);
      // add part from fit
      result->Add(functionC);
    }
    gErrorIgnoreLevel=initialIgnoreLevel;
    // return result
    return result;
  }

  void DrawTheoryCurve(TString filename="", TString plotname="", bool normalize=true, int smoothFactor=0, int rebinFactor=0, int color=kBlack, int linestyle=1, double rangeLow=-1., double rangeHigh=-1., bool errorbands=false, int errorRebinFactor=0, int errorSmoothFactor=0, int verbose=0, bool drawOnlyErrors=false, bool drawRawPlot=false, TString model="", bool smoothcurves=true, TString PS="Parton")
  {
    // this function draws "plot" from "file" into the active canvas
    // modified quantities: NONE
    // used functions: getTheoryPrediction, histogramStyle, useFittedFunctions, makeVariableBinning
    // used enumerators: NONE
    // filename: name of input file
    // plotname: name of histogram
    // normalize: indicates whether absolute or normalized plot is shown
    // smoothFactor: argument of ->Smooth()
    // rebinFactor:  argument of ->Rebin()
    // color: color of histogram
    // rangeLow:  lower plot range of histo
    // rangeHigh: upper plot range of histo
    // errorbands: additional errorbands are calculated and drawn
    // INFORMATION: if errorbands==true is chosen, the varied up/down plots are expected
    //              to have the same name with an additional "_Up"/"_Down" in the end
    // errorRebinFactor: argument of ->Rebin()
    // errorSmoothFactor: argument of ->Smooth()
    // verbose: level of output
    // drawOnlyErrors: draw only error bands but no central curve
    // drawRawPlot: draw in addition the unbinned and unsmoothed plot 
    //              to see whether rebinning and smoothing has changed the shape 
    // model: indicates theory (madgraph, powheg or mcatnalo)
    // smoothcurves: indicates wheter smooth or binned curve is drawn
    // PS: "Parton", "Hadron"
    
/*     if(filename=="/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun/combinedDiffXSecSigMcatnloFall11PF.root"){ */
/*       std::cout << std::endl << "NEW IMPLEMENTATION!!!! " << std::endl <<  plotname << ": "; */
/*       std::cout << " normalize? " << normalize << " , errorbands? " << errorbands << " , drawOnlyErrors? " << drawOnlyErrors << " , model? " << model << " , smoothcurves? " << smoothcurves << std::endl; */
/*     } */

    // get unified variable name from specific file dependend plotnames
    TString plotname2=plotname;
    if(plotname2.BeginsWith("h")) plotname2.Replace(0,1,"");
    plotname2.ReplaceAll("Vis"    ,"");
    plotname2.ReplaceAll("MC@NLO" ,"");
    plotname2.ReplaceAll("MC@NLO2","");
    plotname2.ReplaceAll("POWHEG","" );
    plotname2.ReplaceAll("Top","top" );
    plotname2.ReplaceAll("TTbarM","ttbarMass" );
    plotname2.ReplaceAll("TTbar","ttbar");
    plotname2.ReplaceAll("Lep","lep");
    plotname2.ReplaceAll("Bottom","bq");
    plotname2.ReplaceAll("Gen","");
    if(plotname2.Contains("/")){
      plotname2.ReplaceAll(getStringEntry(plotname2,1)+"/","");
    }
    //if(plotname2.Contains("bq")&&model=="madgraph"&&smoothcurves) verbose=1;
    // create variable bin edges for non smooth curves
    std::map<TString, std::vector<double> > binning_ = makeVariableBinning();
    // output
    if(verbose>1){
      std::cout << std::endl;
      std::cout << "variable:           " << plotname     << std::endl;
      std::cout << "variable (modifed): " << plotname2    << std::endl;
      if(rangeLow!=-1.||rangeHigh!=-1.) std::cout << "range:              " << rangeLow << " .... " << rangeHigh << std::endl;
      std::cout << "file:               " << filename     << std::endl;
      std::cout << "model:              " << model        << std::endl;
      std::cout << "normalize:          " << normalize    << std::endl;
      std::cout << "smoothFactor:       " << smoothFactor << std::endl;
      std::cout << "rebinFactor:        " << rebinFactor  << std::endl;
      std::cout << "color:              " << color        << std::endl;
      std::cout << "draw errorbands:    " << errorbands   << std::endl;
      std::cout << "draw smooth curves: " << smoothcurves << std::endl;
      if(errorbands){
	std::cout << "errorRebinFactor:   " << errorRebinFactor << std::endl;
	std::cout << "errorSmoothFactor:  " << errorSmoothFactor << std::endl;
      }
    }
    // --- 
    // absolute normalization
    // --- 
    // scale to inclusive xSec
    double weight=ttbarCrossSection;
    // consider efficiency in phase space
    //if(model=="powheg" ) weight*=255743./1227541.;
    //if(model=="mcatnlo") weight*=119134./546459; // MaGraph value for the moment
    // BR correction = <BR W->e/mu> * <BR W->qq> (from pdg)
    weight*=0.5*(0.1075+0.1057)*0.676;
    // N-> <e+mu>
    weight*=0.5;
    // output
    if(!normalize&&(model=="powheg"||model=="mcatnlo")&&verbose>0){
      std::cout << "incl xSec(PS): " << 2*weight << " (" << model << ")"<< std::endl;
    }
    // --- 
    // draw raw plot to estimate effect of rebinning and smoothing
    // --- 
    if(drawRawPlot){
      TH1F* raw=getTheoryPrediction(plotname, filename);
      raw->SetMarkerColor(kBlack);
      raw->SetLineColor(kBlack);
      //raw->SetMarkerSize(2);
      //raw->SetMarkerStyle(29);
      if(normalize) raw->Scale(1./(raw->Integral(0,raw->GetNbinsX()+1)));
      // absolute normalization for POWHEG and MCatNlo curve
      if(!normalize&&(model=="powheg"||model=="mcatnlo")) raw->Scale(weight/getInclusiveXSec(raw));
      raw->Scale(1./(raw->GetBinWidth(1)));
      raw->Draw("c same"); 
    }
    // --- 
    // central prediction
    // --- 
    // get plot
    TH1F* load=getTheoryPrediction(plotname, filename);
    TH1F* result=0;
    // replace low statistic parts with fitted curve
    if(smoothcurves) result=useFittedFunctions(load, model, plotname, verbose);
    // variable rebinning for binned curve
    TString plotname3 = plotname2;
    plotname3.ReplaceAll("Gen","");
    if(!smoothcurves){
      // check if maxima of x-axes of loaded plot and new binning is the same
      double newMax=binning_[plotname3][binning_[plotname3].size()-1];
      double oldMax=load->GetBinLowEdge(load->GetNbinsX()+1);
      if(newMax<=oldMax){
	result=(TH1F*)(load->Clone());
	reBinTH1F(*result, binning_[plotname3], verbose);
      }
      // if not: use different rebinning function and ignore bin content in the missing range
      else{
	if(verbose>0){
	  std::cout << "use alternative rebinning function" << std::endl;
	  std::cout << "WARNING: in DrawTheoryCurve: maximum of loaded "<< model;
	  std::cout << " theory curve is smaller then the maximum of the requested binning! Will ignore the missing bin content!" << std::endl;
	}
	result= reBinTH1FIrregularNewBinning(load, binning_[plotname3] , plotname3, false);
      }
    }
    // rename
    TString name=plotname2;
    if(model=="powheg") name+="POWHEG";
    else if(model=="mcatnlo"){
      if(filename.Contains("catnlo" )) name+="MC@NLO";
      else if(filename.Contains("NtupleCteq6m" )) name+="MC@NLO2";
    }
    else if(model=="madgraph") name+="MADGRAPH";
    if(verbose>0) std::cout << "name of theory curve: " << name << std::endl;
    result->SetName(name);
    result->SetTitle(name);
     // configure style
    histogramStyle(*result, kSig, false, 1.2, color);
    result->SetLineStyle(linestyle);
    // smoothing 1
    if(smoothcurves&&smoothFactor){ 
      if(!plotname.Contains("bqEta")) result->Smooth(smoothFactor);
      else{
	TH1F* duplicate=(TH1F*)result->Clone();
	result->Smooth(smoothFactor); 
	// no smoothing for |eta|>1.6
	for(int bin=0; bin<=result->GetNbinsX(); ++bin){
	  if(std::abs(result->GetBinCenter(bin))>1.6) result->SetBinContent(bin, duplicate->GetBinContent(bin));	  
	}
      }
    }  
    // equal rebinning before
    if(smoothcurves&&rebinFactor) result->Rebin(rebinFactor);
    // normalize to area
    if(normalize) result->Scale(1./(result->Integral(0,result->GetNbinsX()+1)));
    // absolute normalization for POWHEG and MCatNlo curve
    if(!normalize&&(model=="powheg"||model=="mcatnlo")) result->Scale(weight/getInclusiveXSec(result));
    // divide by binwidth
    if(smoothcurves) result->Scale(1.0/result->GetBinWidth(1));
    else result=divideByBinwidth(result, false);
    // smoothing 2
    if(smoothcurves&&smoothFactor){
      if(!plotname.Contains("bqEta")) result->Smooth(smoothFactor);
      else{
        TH1F* duplicate=(TH1F*)result->Clone(); 
        result->Smooth(smoothFactor);
        // no smoothing for |eta|>1.6
        for(int bin=0; bin<=result->GetNbinsX(); ++bin){
          if(std::abs(result->GetBinCenter(bin))>1.6) result->SetBinContent(bin, duplicate->GetBinContent(bin));
        }
      }                    
    }   
    // set range
    if(rangeLow!=-1.&&rangeHigh!=-1.) result->GetXaxis()->SetRangeUser(rangeLow, rangeHigh);
    // --- 
    // error bands
    // --- 
    if(errorbands){
      // get values   
      TH1F* loadcentral   =getTheoryPrediction(plotname        , filename);
      TH1F* loadErrorUp   =getTheoryPrediction(plotname+"_Up"  , filename);
      TH1F* loadErrorDown =getTheoryPrediction(plotname+"_Down", filename);
      // for bq and lep error bands: use mcatnlo central curve for correct PS from default running
      // and transfer relative error from parton level studies
      TString plotname4="analyzeTop"+PS+"LevelKinematics";     
      if(PS=="Hadron"){
	if     (plotname.Contains("Lep"   )) plotname4+="Lepton";
	else if(plotname.Contains("Bottom")) plotname4+="Bjets";
      }
      if(plotname.Contains("Vis")) plotname4+="PhaseSpace"; 
      plotname4+="/"+plotname2;

      if(PS=="Hadron"&&(plotname.Contains("Lep")||plotname.Contains("Bottom"))) plotname4+="Gen";
      if(verbose>1) std::cout << plotname << "->" << plotname4 << std::endl;
      TH1F* loadcentral2=0;
      if(verbose>0){
	std::cout << "transfer relative uncertainties to central value from: " << plotname4 << " in " 
		  << "/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun/combinedDiffXSecSigMcatnloFall11PF.root" << std::endl;
      }       
      loadcentral2=getTheoryPrediction(plotname4,"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun/combinedDiffXSecSigMcatnloFall11PF.root");

      TH1F* central=0;    
      TH1F* central2=0;
      TH1F* ErrorUp=0;
      TH1F* ErrorDown=0;
      if(smoothcurves){
	// replace low statistic parts with fitted curve
	central  =useFittedFunctions(loadcentral , model,plotname,        verbose);
	central2 =useFittedFunctions(loadcentral2, model,plotname,        verbose);      
	ErrorUp  =useFittedFunctions(loadErrorUp,  model,plotname+"_Up",  verbose);
	ErrorDown=useFittedFunctions(loadErrorDown,model,plotname+"_Down",verbose);
	// smoothing 1
	if(errorSmoothFactor){
	  central  ->Smooth(errorSmoothFactor);
	  central2 ->Smooth(errorSmoothFactor);
	  ErrorUp  ->Smooth(errorSmoothFactor);
	  ErrorDown->Smooth(errorSmoothFactor);
	}
	// rebinning
	if(errorRebinFactor){ 
	  central  ->Rebin(errorRebinFactor);
	  central2 ->Rebin(errorRebinFactor);
	  ErrorUp  ->Rebin(errorRebinFactor);
	  ErrorDown->Rebin(errorRebinFactor);
	}
      }
      // variable rebinning for binned curve
      if(!smoothcurves){
	// check if maxima of x-axes of loaded plots and new binning is the same
	double newMax         =binning_[plotname3][binning_[plotname3].size()-1];
	double oldMaxcentral  =loadcentral  ->GetBinLowEdge(load->GetNbinsX()+1);
	double oldMaxcentral2 =loadcentral2 ->GetBinLowEdge(loadcentral2->GetNbinsX()+1);
	double oldMaxErrorUp  =loadErrorUp  ->GetBinLowEdge(load->GetNbinsX()+1);
	double oldMaxErrorDown=loadErrorDown->GetBinLowEdge(load->GetNbinsX()+1);
	if(newMax<=oldMaxcentral){
	  central=(TH1F*)(loadcentral->Clone());
	  reBinTH1F(*central, binning_[plotname3], verbose);
	}
	// central curve- if not: use different rebinning function and ignore bin content in the missing range
	else{
	  if(verbose>0){
	    std::cout << "WARNING: in DrawTheoryCurve: maximum of loaded "<< model;
	    std::cout << " central theory curve is smaller then the maximum of the requested binning! Will ignore the missing bin content!" << std::endl;
	  }
	  central =reBinTH1FIrregularNewBinning(loadcentral , binning_[plotname3] , plotname3, false);
	}
	if(newMax<=oldMaxcentral2){
	  central2=(TH1F*)(loadcentral2->Clone());
	  reBinTH1F(*central2, binning_[plotname3], verbose);
	}
	// central curve- if not: use different rebinning function and ignore bin content in the missing range
	else{
	  if(verbose>0){
	    std::cout << "WARNING: in DrawTheoryCurve: maximum of loaded "<< model;
	    std::cout << " central2 theory curve is smaller then the maximum of the requested binning! Will ignore the missing bin content!" << std::endl;
	  }
	  central2=reBinTH1FIrregularNewBinning(loadcentral2, binning_[plotname3] , plotname3, false);
	}	
	if(newMax<=oldMaxErrorUp){
	  ErrorUp=(TH1F*)(loadErrorUp->Clone());
	  reBinTH1F(*ErrorUp, binning_[plotname3], verbose);
	}
	// error up- if not: use different rebinning function and ignore bin content in the missing range
	else{
	  if(verbose>0){
	    std::cout << "WARNING: in DrawTheoryCurve: maximum of loaded "<< model;
	    std::cout << " ErrorUp theory curve is smaller then the maximum of the requested binning! Will ignore the missing bin content!" << std::endl;
	  }
	  ErrorUp=reBinTH1FIrregularNewBinning(loadErrorUp, binning_[plotname3] , plotname3, false);
	}
	if(newMax<=oldMaxErrorDown){
	  ErrorDown=(TH1F*)(loadErrorDown->Clone());
	  reBinTH1F(*ErrorDown, binning_[plotname3], verbose);
	}
	// error down- if not: use different rebinning function and ignore bin content in the missing range
	else{
	  if(verbose>0){
	    std::cout << "WARNING: in DrawTheoryCurve: maximum of loaded "<< model;
	    std::cout << " ErrorDown theory curve is smaller then the maximum of the requested binning! Will ignore the missing bin content!" << std::endl;
	  }
	  ErrorDown=reBinTH1FIrregularNewBinning(loadErrorDown, binning_[plotname3] , plotname3, false);
	}
      }
      // normalize to area
      if(normalize){
	central  ->Scale(1./(central  ->Integral(0,central  ->GetNbinsX()+1)));
	central2 ->Scale(1./(central2 ->Integral(0,central2 ->GetNbinsX()+1)));
	ErrorUp  ->Scale(1./(central  ->Integral(0,central  ->GetNbinsX()+1)));
	ErrorDown->Scale(1./(central  ->Integral(0,central  ->GetNbinsX()+1)));
      }
      // divide by binwidth
      if(smoothcurves){
	central  ->Scale(1.0/central  ->GetBinWidth(1));
	central2 ->Scale(1.0/central  ->GetBinWidth(1));
	ErrorUp  ->Scale(1.0/ErrorUp  ->GetBinWidth(1));
	ErrorDown->Scale(1.0/ErrorDown->GetBinWidth(1));
      }
      else{
	central  =divideByBinwidth(central  , false);
	central2 =divideByBinwidth(central2 , false);
	ErrorUp  =divideByBinwidth(ErrorUp  , false);
	ErrorDown=divideByBinwidth(ErrorDown, false);
      }
      // smoothing 2
      if(smoothcurves&&errorSmoothFactor){
	central  ->Smooth(errorSmoothFactor);
	central2 ->Smooth(errorSmoothFactor);
	ErrorUp  ->Smooth(errorSmoothFactor);
	ErrorDown->Smooth(errorSmoothFactor);
      }
      // create errorbands
      const int constNBINS = central2->GetNbinsX();
      TGraphAsymmErrors *errorBands = new TGraphAsymmErrors(constNBINS);      
      TString errorBandName=name;
      errorBandName.ReplaceAll("MC@NLO2","MC@NLO");
      errorBandName+="errorBand";
      if(verbose>0) std::cout << "name of error bands: " << errorBandName << std::endl;
      errorBands->SetName(errorBandName);
      errorBands->SetTitle(errorBandName);
      // loop bins
      for(Int_t iBin=1; iBin<=constNBINS; iBin++){
	Double_t centralVal = central  ->GetBinContent(iBin);
	Double_t maxVal     = ErrorUp  ->GetBinContent(iBin);
	Double_t minVal     = ErrorDown->GetBinContent(iBin);
	// transfer relative uncertainty from parton level studies to central value for specific phase space
	Double_t centralValue = central2 ->GetBinContent(iBin);
	Double_t maxValue     = centralValue * (1+(maxVal-centralVal)/centralVal);
	Double_t minValue     = centralValue * (1+(minVal-centralVal)/centralVal);
	if(verbose>1&&(plotname.Contains("Lep")||plotname.Contains("Bottom"))){
	  std::cout << std::endl << plotname << ", bin " << iBin << "/" << central->GetNbinsX() << std::endl;
	  std::cout << "central: " << centralValue << std::endl;
	  std::cout << "up     : " << maxValue     << std::endl;
	  std::cout << "down   : " << minValue     << std::endl;
	  std::cout << "c: " << centralVal << std::endl;
	  std::cout << "u: " << maxVal << " ( = " << 1+(maxVal-centralVal)/centralVal << " )" << std::endl;
	  std::cout << "d: " << minVal << " ( = " << 1+(minVal-centralVal)/centralVal << " )" << std::endl;
	}
	// take care of points outside x-plotrange
	if((rangeLow!=-1.&&central2->GetBinCenter(iBin)<rangeLow)||(rangeHigh!=-1.&&central2->GetBinCenter(iBin)>rangeHigh)){
	  errorBands->SetPoint      (iBin, central2->GetBinCenter(iBin), 0.);
	  //std::cout << "errorBands->SetPoint      (iBin, central2->GetBinCenter(iBin), 0.)" << std::endl;
	  errorBands->SetPointEXlow (iBin, 0.);
	  errorBands->SetPointEXhigh(iBin, 0.);
	}
	else {
	  //std::cout << "errorBands->SetPoint      (iBin, central2->GetBinCenter(iBin), centralValue  )" << std::endl;
	  errorBands->SetPoint(iBin, central2->GetBinCenter(iBin), centralValue  );
	  if(smoothcurves){
	    errorBands->SetPointEXlow (iBin, central2->GetXaxis()->GetBinLowEdge(iBin) );
	    errorBands->SetPointEXhigh(iBin, central2->GetXaxis()->GetBinUpEdge (iBin) );
	  }
	  else{
	    errorBands->SetPointEXlow (iBin, 0.5*central2->GetBinWidth(iBin));
	    errorBands->SetPointEXhigh(iBin, 0.5*central2->GetBinWidth(iBin));
	  }
	}
	// symmetrize errors 
	//double difference=std::abs((maxValue-minValue)*0.5);
	//if((maxValue>centralValue&&minValue>centralValue)||(maxValue<centralValue&&minValue<centralValue)){
	//  difference=(std::abs((maxValue-centralValue))+std::abs((minValue-centralValue)))*0.5;
	//}
	if((maxValue>centralValue&&minValue>centralValue)||(maxValue<centralValue&&minValue<centralValue)){
	  double larger = (minValue > maxValue) ? minValue : maxValue;
	  maxValue=larger;
	  minValue=centralValue - (larger - centralValue); 
	}
	errorBands->SetPointEYhigh(iBin, maxValue-centralValue);
	errorBands->SetPointEYlow (iBin, centralValue-minValue);
	// per bin output
	if(verbose>1){
	  std::cout << "bin " << iBin << std::endl;
	  std::cout << "central TH1: " << std::endl;
	  std::cout << "(<x>,y)= (" << central2->GetBinCenter(iBin) << "," << central2->GetBinContent(iBin) << ")" << std::endl;
	  std::cout << "error input:" << std::endl;
	  std::cout << "central, min, max: " << central2->GetBinContent(iBin);
	  std::cout << ", " << ErrorDown->GetBinContent(iBin) << ", ";
	  std::cout <<  ErrorUp->GetBinContent(iBin) << std::endl;
	  std::cout << "error band" << std::endl;
	  std::cout << "(x,y): (" << errorBands->GetX()[iBin] << "," << errorBands->GetY()[iBin] << ")" << std::endl;
	  std::cout << "error xlow:  " << errorBands->GetErrorXlow(iBin)  << std::endl;
	  std::cout << "error xhigh: " << errorBands->GetErrorXhigh(iBin) << std::endl;
	  std::cout << "error ylow:  " << errorBands->GetErrorYlow(iBin)  << std::endl;
	  std::cout << "error yhigh: " << errorBands->GetErrorYhigh(iBin) << std::endl;
	}
      }
      // plotstyle for error bands
      errorBands->SetFillColor(kGray);
      errorBands->SetFillStyle(1001); // NB: explicitly needed, otherwise filling invisible due to default "0"
      errorBands->SetLineColor  (result->GetLineColor());
      errorBands->SetLineWidth  (result->GetLineWidth());
      errorBands->SetMarkerColor(result->GetMarkerColor());
      // draw error bands
      if(smoothcurves) errorBands->Draw("e3 same");
      else errorBands->Draw("e2 same");
    }    
    // --- 
    //draw central value
    // --- 
    if(!drawOnlyErrors){ 
      if(smoothcurves) result->Draw("hist c same"); 
      else result->Draw("hist same");
      if(verbose>0) std::cout << "theory curve drawn" << std::endl;
    }
  }

  void errorWeightedMeanCombination(const TH1& plotMu, const TH1& plotEl, TH1& plotCombination, unsigned int verbose=0){
    // this function does a simple error weighted mean combination 
    // for all bins of "plotMu" and "plotEl" seperately, saving the
    // result in "plotCombination"
    // verbose: level of output
    // modified quantities: plotCombination
    
    // loop bins
    for(int bin=1; bin<=plotCombination.GetNbinsX(); ++bin){
      if(verbose>1){
	std::cout << ", bin" << bin << ":" << std::endl;
	std::cout << "(" << plotCombination.GetBinLowEdge(bin) << " .... " << plotCombination.GetBinLowEdge(bin+1) << ")" << std::endl;
      }
      // consider only non-empty bins
      if(plotMu.GetBinContent(bin)!=0&&plotEl.GetBinContent(bin)!=0){
	double xSecMu     =plotMu.GetBinContent(bin);
	double xSecEl     =plotEl.GetBinContent(bin);
	double xSecErrorMu=plotMu.GetBinError  (bin);
	double xSecErrorEl=plotEl.GetBinError  (bin);
	double combinedxSec=(xSecMu/(xSecErrorMu*xSecErrorMu)+xSecEl/(xSecErrorEl*xSecErrorEl));
	combinedxSec/=(1/(xSecErrorMu*xSecErrorMu)+1/(xSecErrorEl*xSecErrorEl));
	double combinedxSecError=sqrt(1/(1/(xSecErrorMu*xSecErrorMu)+1/(xSecErrorEl*xSecErrorEl)));
	plotCombination.SetBinContent(bin, combinedxSec     );
	plotCombination.SetBinError  (bin, combinedxSecError);
	if(verbose>1){
	  std::cout << "muon:     " << xSecMu       << " +/- " << xSecErrorMu       << std::endl;
	  std::cout << "elec:     " << xSecEl       << " +/- " << xSecErrorEl       << std::endl;
	  std::cout << "combined: " << combinedxSec << " +/- " << combinedxSecError << std::endl;
	}
      }
    }
  }


    
    // ===============================================================
    // ===============================================================
    
    double regParameter(TString variable, TString decayChannel, int verbose=0, bool fullPS=false, bool tau=false, bool hadronPS=false, TString closureTestSpecifier="")
    {
	// this function returns k/value for SVD 
	// unfolding for the corresponding variable
	// modified quantities: NONE
	// used functions: NONE
	// used enumerators: NONE
	// variable: name of variable
	// verbose: level of output
	
	// NB: k-value should be independent from decay channel
	// NB: at the moment k=N(bins) is used!

	double k=-1.0;

	// Unfolding with optimal tau
	
	if(tau){
	    
	    if(decayChannel.Contains("muon")){

	        if(closureTestSpecifier==""){
		// STANDARD data
		    // PAS Binning
		    //if(variable == "topPt")          k = (fullPS) ? 2.07 : 4.17523;
		    //else if(variable == "topY")      k = (fullPS) ? 2.86 : 2.96168;
		    //else if(variable == "ttbarPt")   k = (fullPS) ? 2.04 : 2.19367; 
		    //else if(variable == "ttbarY")    k = (fullPS) ? 1.50 : 2.11527;
		    //else if(variable == "ttbarMass") k = (fullPS) ? 1.32 : 2.49955;
		    // New Binning Revision
		    if     (variable.Contains("topPt")    ) k = (fullPS) ? 3.67 : 3.67; // hadron level not important
		    else if(variable.Contains("topY" )    ) k = (fullPS) ? 3.70 : 3.70; // hadron level not important
		    else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.11 : 2.11; // hadron level not important
		    else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 3.09 : 3.09; // hadron level not important
		    else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.26 : 1.26; // hadron level not important
		    else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.88 : ((hadronPS) ? 2.30  : 2.88);
		    else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.53 : ((hadronPS) ? 2e-6  : 1.53); 
		    else if(variable.Contains("bqPt")     ) k = (fullPS) ? 3.64 : ((hadronPS) ? 4.30  : 3.64); 
		    else if(variable.Contains("bqEta")    ) k = (fullPS) ? 4.24 : ((hadronPS) ? 3.45  : 4.24);
		}
		else if(closureTestSpecifier=="Up"){
		  // so far only redetermined for parton lv PS, hadron are copied from standard
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 2.47 : 2.47; // hadron level not important
		  else if(variable.Contains("topY" )    ) k = (fullPS) ? 3.54 : 3.54; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.33 : 2.33; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 3.04 : 3.04; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 0.77 : 0.77; // hadron level not important
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.61 : (hadronPS) ? 2.31  : 2.61;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.47 : (hadronPS) ? 1e-5  : 1.47; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 2.57 : (hadronPS) ? 3.06  : 2.57; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 4.12 : (hadronPS) ? 3.94  : 4.12;
		}
		else if(closureTestSpecifier=="Down"){
		    // so far only redetermined for parton lv PS, hadron are copied from standard
		    if     (variable.Contains("topPt")    ) k = (fullPS) ? 4.54 : 4.54; // hadron level not important
		    else if(variable.Contains("topY" )    ) k = (fullPS) ? 4.07 : 4.07; // hadron level not important
		    else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.56 : 2.56; // hadron level not important
		    else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 3.26 : 3.26; // hadron level not important
		    else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.67 : 1.67; // hadron level not important
		    else if(variable.Contains("lepPt")    ) k = (fullPS) ? 3.26 : (hadronPS) ? 2.31  : 3.26;
		    else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.72 : (hadronPS) ? 1e-5  :1.72; 
		    else if(variable.Contains("bqPt")     ) k = (fullPS) ? 4.51 : (hadronPS) ? 3.06  : 4.51; 
		    else if(variable.Contains("bqEta")    ) k = (fullPS) ? 4.69 : (hadronPS) ? 3.94  : 4.69;
		}
		else if(closureTestSpecifier=="Up0p015"){
		  // so far only redetermined for parton lv PS, hadron are copied from standard
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 2.59 : 2.59; // hadron level not important
		  else if(variable.Contains("topY" )    ) k = (fullPS) ? 3.59 : 3.59; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.37 : 2.37; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 3.06 : 3.06; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 0.82 : 0.82; // hadron level not important
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.66 : (hadronPS) ? 2.31  : 2.66;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.49 : (hadronPS) ? 1e-5  : 1.49; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 2.69 : (hadronPS) ? 3.06  : 2.69; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 4.17 : (hadronPS) ? 3.94  : 4.17;
		}
		else if(closureTestSpecifier=="Down0p015"){
		    // so far only redetermined for parton lv PS, hadron are copied from standard
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 5.13 : 5.13; // hadron level not important
		  else if(variable.Contains("topY" )    ) k = (fullPS) ? 4.18 : 4.18; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.84 : 2.84; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 3.31 : 3.31; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.91 : 1.91; // hadron level not important
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 3.39 : (hadronPS) ? 2.31  : 3.39;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.76 : (hadronPS) ? 1e-5  :1.76; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 5.06 : (hadronPS) ? 3.06  : 5.06; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 4.79 : (hadronPS) ? 3.94  : 4.79;
		}
		else if(closureTestSpecifier=="500"){
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 3.55 : 3.61; // hadron level not important
		  else if(variable.Contains("topY" )    ) k = (fullPS) ? 3.66 : 3.82; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.46 : 2.53; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 2.99 : 3.13; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.38 : 1.38; // hadron level not important
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.85 : (hadronPS) ? 2.31  : 2.95;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.54 : (hadronPS) ? 1e-5  : 1.61; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 3.47 : (hadronPS) ? 3.06  : 3.59; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 4.23 : (hadronPS) ? 3.94  : 4.38;
		}
		else if(closureTestSpecifier=="750"){
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 3.38 : 3.61; // hadron level not important
		  else if(variable.Contains("topY" )    ) k = (fullPS) ? 3.78 : 3.82; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.51 : 2.53; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 3.10 : 3.13; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.34 : 1.38; // hadron level not important
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.91 : (hadronPS) ? 2.31  : 2.95;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.59 : (hadronPS) ? 1e-5  : 1.61; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 3.41 : (hadronPS) ? 3.06  : 3.59; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 4.34 : (hadronPS) ? 3.94  : 4.38;
		}
		else if(closureTestSpecifier=="750x0p03" || closureTestSpecifier=="750x0p1"){
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 3.59 : 3.61; // hadron level not important
		  else if(variable.Contains("topY" )    ) k = (fullPS) ? 3.81 : 3.82; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.53 : 2.53; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 3.12 : 3.13; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.37 : 1.38; // hadron level not important
		  //old:
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.91 : (hadronPS) ? 2.31  : 2.95;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.59 : (hadronPS) ? 1e-5  : 1.61; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 3.41 : (hadronPS) ? 3.06  : 3.59; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 4.34 : (hadronPS) ? 3.94  : 4.38;
		}
		else if(closureTestSpecifier=="750x0p25"){
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 3.55 : 3.61; // hadron level not important
		  else if(variable.Contains("topY" )    ) k = (fullPS) ? 3.81 : 3.82; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.53 : 2.53; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 3.12 : 3.13; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.37 : 1.38; // hadron level not important
		  //old:
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.91 : (hadronPS) ? 2.31  : 2.95;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.59 : (hadronPS) ? 1e-5  : 1.61; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 3.41 : (hadronPS) ? 3.06  : 3.59; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 4.34 : (hadronPS) ? 3.94  : 4.38;
		}
		else if(closureTestSpecifier=="750x0p5"){
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 3.50 : 3.61; // hadron level not important
		  else if(variable.Contains("topY" )    ) k = (fullPS) ? 3.80 : 3.82; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.52 : 2.53; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 3.11 : 3.13; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.36 : 1.38; // hadron level not important
		  //old:
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.91 : (hadronPS) ? 2.31  : 2.95;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.59 : (hadronPS) ? 1e-5  : 1.61; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 3.41 : (hadronPS) ? 3.06  : 3.59; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 4.34 : (hadronPS) ? 3.94  : 4.38;
		}
		else if(closureTestSpecifier=="750x2"){
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 3.19 : 3.61; // hadron level not important
		  else if(variable.Contains("topY" )    ) k = (fullPS) ? 3.73 : 3.82; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.48 : 2.53; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 3.05 : 3.13; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.31 : 1.38; // hadron level not important
		  // old:
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.91 : (hadronPS) ? 2.31  : 2.95;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.59 : (hadronPS) ? 1e-5  : 1.61; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 3.41 : (hadronPS) ? 3.06  : 3.59; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 4.34 : (hadronPS) ? 3.94  : 4.38;
		}
		else if(closureTestSpecifier=="750x4"){
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 2.66 : 3.61; // hadron level not important
		  else if(variable.Contains("topY" )    ) k = (fullPS) ? 3.65 : 3.82; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.43 : 2.53; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 2.99 : 3.13; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.25 : 1.38; // hadron level not important
		  // old:
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.91 : (hadronPS) ? 2.31  : 2.95;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.59 : (hadronPS) ? 1e-5  : 1.61; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 3.41 : (hadronPS) ? 3.06  : 3.59; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 4.34 : (hadronPS) ? 3.94  : 4.38;
		}
		else if(closureTestSpecifier.Contains("NoDistort")){
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 3.61 : 3.61; // hadron level not important
		  else if(variable.Contains("topY" )    ) k = (fullPS) ? 3.82 : 3.82; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.53 : 2.53; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 3.13 : 3.13; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.38 : 1.38; // hadron level not important
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.95 : (hadronPS) ? 2.31  : 2.95;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.61 : (hadronPS) ? 1e-5  : 1.61; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 3.59 : (hadronPS) ? 3.06  : 3.59; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 4.38 : (hadronPS) ? 3.94  : 4.38;
		}

            }
	    else if (decayChannel.Contains("electron")){

	        if(closureTestSpecifier==""){
		// STANDARD data
		    // PAS Binning
		    //if(variable == "topPt")          k = (fullPS) ? 3.30 : 4.17523;
		    //else if(variable == "topY")      k = (fullPS) ? 4.23 : 2.96168;
		    //else if(variable == "ttbarPt")   k = (fullPS) ? 2.98 : 2.19367; 
		    //else if(variable == "ttbarY")    k = (fullPS) ? 2.71 : 2.11527;
		    //else if(variable == "ttbarMass") k = (fullPS) ? 3.40 : 2.49955;
		    // New Binning Revision
		    if     (variable.Contains("topPt")    ) k = (fullPS) ? 3.57 : 3.57; // hadron level not important
		    else if(variable.Contains("topY")     ) k = (fullPS) ? 3.24 : 3.24; // hadron level not important
		    else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.18 : 2.18; // hadron level not important
		    else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 2.72 : 2.72; // hadron level not important
		    else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.31 : 1.31; // hadron level not important
		    else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.66 : ((hadronPS) ? 2.33   : 2.66);
		    else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.07 : ((hadronPS) ? 2e-6   : 1.07); 
		    else if(variable.Contains("bqPt")     ) k = (fullPS) ? 3.52 : ((hadronPS) ? 4.07   : 3.52); 
		    else if(variable.Contains("bqEta")    ) k = (fullPS) ? 3.78 : ((hadronPS) ? 3.04   : 3.78);
		}
		else if(closureTestSpecifier=="Up"){
		  // so far only redetermined for parton lv PS, hadron are copied from standard
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 2.36 : 2.36; // hadron level not important
		  else if(variable.Contains("topY")     ) k = (fullPS) ? 3.12 : 3.12; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.06 : 2.06; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 2.72 : 2.72; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 0.78 : 0.78; // hadron level not important
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.52 : (hadronPS) ? 2.32   : 2.52;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 0.91 : (hadronPS) ? 1e-5   : 0.91; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 2.38 : (hadronPS) ? 2.86   : 2.38; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 3.66 : (hadronPS) ? 3.44   : 3.66;
		}
		else if(closureTestSpecifier=="Down"){
		  // so far only redetermined for parton lv PS, hadron are copied from standard
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 4.26 : 4.26; // hadron level not important
		  else if(variable.Contains("topY")     ) k = (fullPS) ? 3.62 : 3.62; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.43 : 2.43; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 2.93 : 2.93; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.78 : 1.78; // hadron level not important
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 3.01 : (hadronPS) ? 2.32   : 3.01;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.05 : (hadronPS) ? 1e-5   : 1.05; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 4.14 : (hadronPS) ? 2.86   : 4.14; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 4.21 : (hadronPS) ? 3.44   : 4.21;
		}
		else if(closureTestSpecifier=="Up0p015"){
		  // so far only redetermined for parton lv PS, hadron are copied from standard
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 2.47 : 2.47; // hadron level not important
		  else if(variable.Contains("topY")     ) k = (fullPS) ? 3.18 : 3.18; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.09 : 2.09; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 2.74 : 2.74; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 0.84 : 0.84; // hadron level not important
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.56 : (hadronPS) ? 2.32   : 2.56;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.05 : (hadronPS) ? 1e-5   : 1.05; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 2.48 : (hadronPS) ? 2.86   : 2.48; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 3.75 : (hadronPS) ? 3.44   : 3.75;
		}
		else if(closureTestSpecifier=="Down0p015"){
		  // so far only redetermined for parton lv PS, hadron are copied from standard
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 4.75 : 4.75; // hadron level not important
		  else if(variable.Contains("topY")     ) k = (fullPS) ? 3.75 : 3.75; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.49 : 2.49; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 2.97 : 2.97; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.99 : 1.99; // hadron level not important
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 3.10 : (hadronPS) ? 2.32   : 3.10;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.27 : (hadronPS) ? 1e-5   : 1.27; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 4.66 : (hadronPS) ? 2.86   : 4.66; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 4.35 : (hadronPS) ? 3.44   : 4.35;
		}
		if(closureTestSpecifier=="500"){
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 3.34 : 3.40; // hadron level not important
		  else if(variable.Contains("topY")     ) k = (fullPS) ? 3.27 : 3.36; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.17 : 2.23; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 2.68 : 2.78; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.42 : 1.43; // hadron level not important
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.67 : (hadronPS) ? 2.32   : 2.76;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.10 : (hadronPS) ? 1e-5   : 0.99; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 3.22 : (hadronPS) ? 2.86   : 3.32; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 3.83 : (hadronPS) ? 3.44   : 3.91;
		}
		if(closureTestSpecifier=="750"){
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 3.22 : 3.40; // hadron level not important
		  else if(variable.Contains("topY")     ) k = (fullPS) ? 3.36 : 3.36; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.20 : 2.23; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 2.77 : 2.78; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.39 : 1.43; // hadron level not important
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.73 : (hadronPS) ? 2.32   : 2.76;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.14 : (hadronPS) ? 1e-5   : 0.99; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 3.17 : (hadronPS) ? 2.86   : 3.32; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 3.91 : (hadronPS) ? 3.44   : 3.91;
		}
		if(closureTestSpecifier=="750x0p03" || closureTestSpecifier=="750x0p1"){
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 3.39 : 3.40; // hadron level not important
		  else if(variable.Contains("topY")     ) k = (fullPS) ? 3.39 : 3.36; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.23 : 2.23; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 2.80 : 2.78; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.42 : 1.43; // hadron level not important
		  //old:
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.73 : (hadronPS) ? 2.32   : 2.76;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.14 : (hadronPS) ? 1e-5   : 0.99; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 3.17 : (hadronPS) ? 2.86   : 3.32; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 3.91 : (hadronPS) ? 3.44   : 3.91;
		}
		if(closureTestSpecifier=="750x0p25"){
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 3.36 : 3.40; // hadron level not important
		  else if(variable.Contains("topY")     ) k = (fullPS) ? 3.39 : 3.36; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.23 : 2.23; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 2.79 : 2.78; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.42 : 1.43; // hadron level not important
		  //old:
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.73 : (hadronPS) ? 2.32   : 2.76;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.14 : (hadronPS) ? 1e-5   : 0.99; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 3.17 : (hadronPS) ? 2.86   : 3.32; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 3.91 : (hadronPS) ? 3.44   : 3.91;
		}
		if(closureTestSpecifier=="750x0p5"){
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 3.31 : 3.40; // hadron level not important
		  else if(variable.Contains("topY")     ) k = (fullPS) ? 3.38 : 3.36; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.20 : 2.23; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 2.78 : 2.78; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.40 : 1.43; // hadron level not important
		  //old:
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.73 : (hadronPS) ? 2.32   : 2.76;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.14 : (hadronPS) ? 1e-5   : 0.99; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 3.17 : (hadronPS) ? 2.86   : 3.32; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 3.91 : (hadronPS) ? 3.44   : 3.91;
		}
		if(closureTestSpecifier=="750x2"){
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 3.06 : 3.40; // hadron level not important
		  else if(variable.Contains("topY")     ) k = (fullPS) ? 3.32 : 3.36; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.18 : 2.23; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 2.73 : 2.78; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.35 : 1.43; // hadron level not important
		  // old:
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.73 : (hadronPS) ? 2.32   : 2.76;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.14 : (hadronPS) ? 1e-5   : 0.99; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 3.17 : (hadronPS) ? 2.86   : 3.32; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 3.91 : (hadronPS) ? 3.44   : 3.91;
		}
		if(closureTestSpecifier=="750x4"){
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 2.79 : 3.40; // hadron level not important
		  else if(variable.Contains("topY")     ) k = (fullPS) ? 3.24 : 3.36; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.13 : 2.23; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 2.67 : 2.78; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.29 : 1.43; // hadron level not important
		  // old:
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.73 : (hadronPS) ? 2.32   : 2.76;
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.14 : (hadronPS) ? 1e-5   : 0.99; 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 3.17 : (hadronPS) ? 2.86   : 3.32; 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 3.91 : (hadronPS) ? 3.44   : 3.91;
		}
		else if(closureTestSpecifier.Contains("NoDistort")){
		  if     (variable.Contains("topPt")    ) k = (fullPS) ? 3.40 : 3.40; // hadron level not important
		  else if(variable.Contains("topY")     ) k = (fullPS) ? 3.40 : 3.40; // hadron level not important
		  else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 2.23 : 2.23; // hadron level not important
		  else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 2.80 : 2.80; // hadron level not important
		  else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.42 : 1.42; // hadron level not important
		  else if(variable.Contains("lepPt")    ) k = (fullPS) ? 2.77 : ((hadronPS) ? 2.32   : 2.77);
		  else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.15 : ((hadronPS) ? 1e-5   : 1.15); 
		  else if(variable.Contains("bqPt")     ) k = (fullPS) ? 3.32 : ((hadronPS) ? 2.86   : 3.32); 
		  else if(variable.Contains("bqEta")    ) k = (fullPS) ? 3.97 : ((hadronPS) ? 3.44   : 3.97);
		}

            }
	    else if(decayChannel.Contains("combined")){  
	      // New Binning Revision
	      if     (variable.Contains("topPt")    ) k = (fullPS) ? 5.13 : 5.13; // hadron level not important
	      else if(variable.Contains("topY")     ) k = (fullPS) ? 4.92 : 4.92; // hadron level not important
	      else if(variable.Contains("ttbarPt")  ) k = (fullPS) ? 3.23 : 3.23; // hadron level not important 
	      else if(variable.Contains("ttbarY")   ) k = (fullPS) ? 4.13 : 4.13; // hadron level not important
	      else if(variable.Contains("ttbarMass")) k = (fullPS) ? 1.81 : 1.81; // hadron level not important
	      else if(variable.Contains("lepPt")    ) k = (fullPS) ? 3.93 : ((hadronPS) ? 3.26  : 3.93);
	      else if(variable.Contains("lepEta")   ) k = (fullPS) ? 1.87 : ((hadronPS) ? 0.0006: 1.87); 
	      else if(variable.Contains("bqPt")     ) k = (fullPS) ? 5.05 : ((hadronPS) ? 5.9   : 5.06); 
	      else if(variable.Contains("bqEta")    ) k = (fullPS) ? 5.68 : ((hadronPS) ? 4.64  : 5.68);
	    }
	}
	else{
	    
	    // Default unfolding with number of bins
	    // PAS Binning
	    //if(variable == "topPt")          k =  5;
	    //else if(variable == "topY")      k =  8;
	    //else if(variable == "ttbarPt")   k =  5;
	    //else if(variable == "ttbarY")    k =  8;
	    //else if(variable == "ttbarMass") k =  5;
	    // New Binning
	    if     (variable.Contains("topPt")    ) k =  7;
	    else if(variable.Contains("topY")     ) k = 10;
	    else if(variable.Contains("ttbarPt")  ) k =  6;
	    else if(variable.Contains("ttbarY")   ) k = 10;
	    else if(variable.Contains("ttbarMass")) k =  7;
	    else if(variable.Contains("lepPt")    ) k = 11;
	    else if(variable.Contains("lepEta")   ) k = 14;
	    else if(variable.Contains("bqPt")     ) k =  5;
	    else if(variable.Contains("bqEta")    ) k =  8;
	}      
	
	// output
	if(verbose>1){
	    if(tau) std::cout << "tau";
	    else std::cout << "k";
	    std::cout << "(" << variable << ") = " << k << std::endl;
	}
	// check result
	if(k<0){
	    std::cout << "ERROR in regParameter:" << std::endl; 
	    std::cout << "invalid k=" << k << " for variable " << variable << std::endl;
	    std::cout << "unknown variable " << variable << " or decayChannel " << decayChannel << std::endl;
	}
	return k;
    }
    
    double getValue(TString fileName, TString canvName, TString plotName, int bin){
      // get file
      TFile* file = TFile::Open(fileName);
      if(!file||(file->IsZombie())){ 
	std::cout << "WARNING in getUncFromRootFile: file " << fileName << " not found" << std::endl;
	return 0.;
      }
      // get canvas
      TCanvas* canv=(TCanvas*)file->Get(canvName);
      if(!canv){ 
	std::cout << "WARNING in getUncFromRootFile: canvas " << canvName << " in file " << fileName << " not found" << std::endl;
	return 0.;
      }
      // get plot
      TH1F* plot  = (TH1F*)canv->GetPrimitive(plotName);
      if(!plot){ 
	std::cout << "WARNING in getUncFromRootFile: plot " << plotName << " within canvas " << canvName << " in file " << fileName << " not found" << std::endl;
	return 0.;
      }
      // get uncertainty value
      if(plot->GetNbinsX()+1<bin||bin<0){
	std::cout << "WARNING in getUncFromRootFile: bin " << bin << " in plot " << plotName << " within canvas " << canvName << " in file " << fileName << " not found (" << plot->GetNbinsX() << " bins)" << std::endl;
	return 0.;
      }
      if(plot->GetBinContent(bin)<0){
	std::cout << "WARNING in getUncFromRootFile: bin " << bin << " in plot " << plotName << " within canvas " << canvName << " in file " << fileName << " has a negative value" << std::endl;
      }
      std::cout << plot->GetBinContent(bin) << std::endl;
      return plot->GetBinContent(bin);
    }


    // ===============================================================
    // ===============================================================

    template <class T> bool SortVectorPrescription(const T& a, const T& b){
	// This function is used to defined
	// a sorting presciption to be used 
	// with std::sort
	return a < b;
    }	

    // ===============================================================
    // ===============================================================
  
    
    

#ifdef DILEPTON_MACRO
}
#endif

#endif
