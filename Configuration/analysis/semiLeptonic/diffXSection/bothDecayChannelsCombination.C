#include "basicFunctions.h"

void bothDecayChannelsCombination(double luminosity=1143, bool save=true, unsigned int verbose=0, TString inputFolderName="TOP2011/110819_AnalysisRun"){
	
  // ---
  //    Setup
  // ---
  // set root style
	
  TStyle myStyle("HHStyle","HHStyle");
  setHHStyle(myStyle);
  myStyle.cd();
  TGaxis::SetMaxDigits(2);
  gROOT->SetStyle("HHStyle");

  // say if closure test with reweighted m(ttbar) on parton level is done
  // will plot additionally the modified diff. norm. xSec on parton level
  bool reweightClosure=false;
  // Zprime pseudo data test
  // choose "", "500" or "750"
  TString zprime="";
  // adjust data label
  TString dataLabel= (reweightClosure||zprime!="") ? "Pseudo-Data": "Data"; 
  // choose if xSec are extrapolated to whole phase space
  bool extrapolate=false;
  TString PS="";
  if(!extrapolate)PS="PhaseSpace";
  // GOSSIE quick fix: cut of m(ttbar) below 354 GeV
  bool cutTtbarMass=true;
  // decay channels
  enum channel {kMuon, kElectron};
  // define folders where XSec plots are stored
  TString xSecFolder = "xSec";
  // get right file by luminosity
  TString dataSample="";
  if(luminosity<50) dataSample="2010";
  if(luminosity>50) dataSample="2011";
  // define muon and electron input rootfile
  std::map<unsigned int, TFile*> files_;
  files_[kMuon    ] = new (TFile)("diffXSecTopSemiMu"  +dataSample+".root");
  files_[kElectron] = new (TFile)("diffXSecTopSemiElec"+dataSample+".root");
  // define variables for combination
  // NOTE: these must be included in xSecVariables_ 
  // in analyzeHypothesisKinFit.C and combineTopDiffXSecUncertainties.C
  std::vector<TString> xSecVariables_;
  TString xSecVariables[] ={"topPt", "topY", "ttbarPt", "ttbarMass", "ttbarY", "lepPt" ,"lepEta", "topPtNorm", "topYNorm", "ttbarPtNorm", "ttbarMassNorm", "ttbarYNorm", "lepPtNorm" ,"lepEtaNorm"};
  xSecVariables_.insert( xSecVariables_.begin(), xSecVariables, xSecVariables + sizeof(xSecVariables)/sizeof(TString) );
  // create plot container for combined e+#mu plots
  std::map< TString, std::map <unsigned int, TH1F*> > histo_;
  std::map< TString, std::map <unsigned int, TCanvas*> > canvas_;
  // ---
  //    Combination for all kinematic variables x systematic variations available
  // ---
  // loop systematic variations
  for(unsigned int sys=sysNo; sys<=sysDiBosDown; ++sys){
    TString subfolder=sysLabel(sys);
    // loop variables
    for(unsigned int i=0; i<xSecVariables_.size(); ++i){
      // get canvas
      TCanvas* canvasMu   = (TCanvas*)(files_[kMuon    ]->Get(xSecFolder+"/"+subfolder+"/"+xSecVariables_[i]));
      TCanvas* canvasTheo = (TCanvas*)(files_[kMuon    ]->Get(xSecFolder+"/"+sysLabel(sysNo)+"/"+xSecVariables_[i]));
      TCanvas* canvasEl   = (TCanvas*)(files_[kElectron]->Get(xSecFolder+"/"+subfolder+"/"+xSecVariables_[i]));
      if(canvasMu&&canvasEl){
	// get data plots for all systematics
	TString plotName = xSecVariables_[i];
	plotName.ReplaceAll("Norm","");
	TH1F* plotMu   = (TH1F*)canvasMu  ->GetPrimitive(plotName+"kData");
	TH1F* plotTheo = (TH1F*)canvasTheo->GetPrimitive(plotName        );
	TH1F* plotEl   = (TH1F*)canvasEl  ->GetPrimitive(plotName+"kData");
	if(plotMu&&plotEl){ 
	  if(verbose>1){
	    std::cout << "plot "+plotName+"kData in "+xSecFolder+"/"+subfolder+"/"+xSecVariables_[i];
	    std::cout << " for both channels found!" << std::endl;
	  }
	  // combine the results
	  TH1F* plotCombination=(TH1F*)(plotMu->Clone());
	  plotCombination->Scale(0.0);
	  // loopins
	  for(int bin=1; bin<=plotCombination->GetNbinsX(); ++bin){
	    // consider only non-empty bins
	    if(plotMu->GetBinContent(bin)!=0&&plotEl->GetBinContent(bin)!=0){
	      double xSecMu     =plotMu->GetBinContent(bin);
	      double xSecEl     =plotEl->GetBinContent(bin);
	      double xSecErrorMu=plotMu->GetBinError  (bin);
	      double xSecErrorEl=plotEl->GetBinError  (bin);
	      double combinedxSec=(xSecMu/(xSecErrorMu*xSecErrorMu)+xSecEl/(xSecErrorEl*xSecErrorEl));
	      combinedxSec/=(1/(xSecErrorMu*xSecErrorMu)+1/(xSecErrorEl*xSecErrorEl));
	      double combinedxSecError=sqrt(1/(1/(xSecErrorMu*xSecErrorMu)+1/(xSecErrorEl*xSecErrorEl)));
	      plotCombination->SetBinContent(bin, combinedxSec     );
	      plotCombination->SetBinError  (bin, combinedxSecError);
	      if(verbose>1){
		std::cout << xSecVariables_[i] << ", bin" << bin << ":" << std::endl;
		std::cout << "mu: " << xSecMu << " +/- " << xSecErrorMu << std::endl;
		std::cout << "el: " << xSecEl << " +/- " << xSecErrorEl << std::endl;
		std::cout << "combined: " << combinedxSec << " +/- " << combinedxSecError << std::endl;
	      }
	    }
	  }
	  // ---
	  //    additional histos for reweighting closure test
	  // ---
	  if(reweightClosure&&sys==sysNo){
	    // get reweighted samples
	    TString muReweighted="/afs/naf.desy.de/group/cms/scratch/tophh/TOP2011/110819_AnalysisRun/muonDiffXSecSigMadD6TSummer11ReweightedttbarMassUpPF.root";
	    TString elReweighted="/afs/naf.desy.de/group/cms/scratch/tophh/TOP2011/110819_AnalysisRun/elecDiffXSecSigMadD6TSummer11ReweightedttbarMassUpPF.root";
	    TFile* mufile = new (TFile)(muReweighted);
	    TFile* elfile = new (TFile)(elReweighted);
	    // get plot
	    TString partonPlot="analyzeTopPartonLevelKinematicsPhaseSpace/"+plotName;
	    histo_["reweighted"+plotName     ][kSig] = (TH1F*)(mufile->Get(partonPlot)->Clone("mu"+plotName));
	    histo_["reweighted"+plotName+"El"][kSig] = (TH1F*)(elfile->Get(partonPlot)->Clone("el"+plotName));
	    histo_["reweighted"+plotName     ][kSig]->Add(histo_["reweighted"+plotName+"El"][kSig]);
	    // apply standard rebinning
	    std::map<TString, std::vector<double> > binning_ = makeVariableBinning();
	    reBinTH1F(*histo_["reweighted"+plotName][kSig], binning_[plotName], verbose-1);
	    // scale to unit area
	    histo_["reweighted"+plotName][kSig]->Scale(1/histo_["reweighted"+plotName][kSig]->Integral(0,histo_["reweighted"+plotName][kSig]->GetNbinsX()+1));
	    // divide by binwidth
	    histo_["reweighted"+plotName][kSig]=divideByBinwidth(histo_["reweighted"+plotName][kSig], verbose-1);
	    // delete files
	    //delete mufile;
	    //mufile=NULL;
	    //delete elfile;
	    //elfile=NULL;
	    //set style
	    histogramStyle(*histo_["reweighted"+plotName][kSig], kSig, false, 1.2, kRed+1);
	    histo_["reweighted"+plotName][kSig]->SetLineColor(kMagenta);
	  }
	  // ---
	  //    additional histos for z prime closure test
	  // ---
	  if(zprime!=""&&sys==sysNo){
	    TString muSig="/afs/naf.desy.de/group/cms/scratch/tophh/TOP2011/110819_AnalysisRun/muonDiffXSecSigMadD6TSummer11PF.root";
	    TString elSig="/afs/naf.desy.de/group/cms/scratch/tophh/TOP2011/110819_AnalysisRun/elecDiffXSecSigMadD6TSummer11PF.root";
	    TString muZprime="/afs/naf.desy.de/group/cms/scratch/tophh/TOP2011/110819_AnalysisRun/muonDiffXSecZPrime_M"+zprime+"_W"+zprime+"0_MadSummer11PF.root";
	    TString elZprime="/afs/naf.desy.de/group/cms/scratch/tophh/TOP2011/110819_AnalysisRun/elecDiffXSecZPrime_M"+zprime+"_W"+zprime+"0_MadSummer11PF.root";
	    TFile* muSigfile = new (TFile)(muSig);
	    TFile* elSigfile = new (TFile)(elSig);
	    TFile* muZprimefile = new (TFile)(muZprime);
	    TFile* elZprimefile = new (TFile)(elZprime);
	     // get plot
	    TString partonPlot="analyzeTopPartonLevelKinematicsPhaseSpace/"+plotName;
	    histo_["modified"+plotName           ][kSig] = (TH1F*)(muSigfile   ->Get(partonPlot)->Clone("modified"+plotName));
	    histo_["modified"+plotName+"ElSig"   ][kSig] = (TH1F*)(elSigfile   ->Get(partonPlot)->Clone("elSig"+plotName   ));
	    histo_["modified"+plotName+"muZprime"][kSig] = (TH1F*)(muZprimefile->Get(partonPlot)->Clone("zprime"+plotName  ));
	    histo_["modified"+plotName+"ElZprime"][kSig] = (TH1F*)(elZprimefile->Get(partonPlot)->Clone("zprime"+plotName  ));
	    // relative apply lumiweight- needed because these are different samples
	    histo_["modified"+plotName           ][kSig]->Scale(lumiweight(kSig, luminosity, 0, "muon"    ));
	    histo_["modified"+plotName+"ElSig"   ][kSig]->Scale(lumiweight(kSig, luminosity, 0, "electron"));
	    double zPrimeLumiWeight=1;
	    if     (zprime=="500") zPrimeLumiWeight=(10*16.2208794979645*luminosity)/232074;
	    else if(zprime=="750") zPrimeLumiWeight=(10*3.16951400706147*luminosity)/206525;
	    histo_["modified"+plotName+"muZprime"][kSig]->Scale(zPrimeLumiWeight);
	    histo_["modified"+plotName+"ElZprime"][kSig]->Scale(zPrimeLumiWeight);
	    // add plots
	    histo_["modified"+plotName][kSig]->Add(histo_["modified"+plotName+"ElSig"   ][kSig]);
	    histo_["modified"+plotName][kSig]->Add(histo_["modified"+plotName+"muZprime"][kSig]);
	    histo_["modified"+plotName][kSig]->Add(histo_["modified"+plotName+"ElZprime"][kSig]);
	    // apply standard rebinning
	    std::map<TString, std::vector<double> > binning_ = makeVariableBinning();
	    reBinTH1F(*histo_["modified"+plotName][kSig], binning_[plotName], verbose-1);
	    // scale to unit area
	    histo_["modified"+plotName][kSig]->Scale(1/histo_["modified"+plotName][kSig]->Integral(0,histo_["modified"+plotName][kSig]->GetNbinsX()+1));
	    // divide by binwidth
	    histo_["modified"+plotName][kSig]=divideByBinwidth(histo_["modified"+plotName][kSig], verbose-1);
	    // set style
	    histogramStyle(*histo_["modified"+plotName][kSig], kSig, false, 1.2, kRed+1);
	    histo_["modified"+plotName][kSig]->SetLineColor(kMagenta);
	  }
	  // adapt plot style
	  double max = plotTheo->GetMaximum();
	  max*=1.5;
	  plotTheo->GetXaxis()->SetNoExponent(true);
	  if(max>1&&max<100) plotTheo->GetYaxis()->SetNoExponent(true);
	  else plotTheo->GetYaxis()->SetNoExponent(false);
	  TString yTitle=plotTheo->GetYaxis()->GetTitle();
	  TString xTitle=plotTheo->GetXaxis()->GetTitle();
	  if(verbose>1)std::cout << "original title: " << yTitle  << std::endl;
	  yTitle.ReplaceAll("#mu","l");
	  xTitle.ReplaceAll("#mu","l");
	  if(verbose>1)std::cout << "new title: " << yTitle  << std::endl;
	  plotTheo->GetYaxis()->SetTitle(yTitle);
	  plotTheo->GetXaxis()->SetTitle(xTitle);
	  // create Canvas
	  TCanvas* combicanvas = new TCanvas("combicanvas", "combicanvas", 600, 600);
	  // get style from old canvas
	  combicanvas->SetLogy(canvasMu->GetLogy());
	  // plot into canvas
	  combicanvas->cd(0);
	  combicanvas->SetBottomMargin(myStyle.GetPadBottomMargin());
	  combicanvas->SetLeftMargin(myStyle.GetPadLeftMargin());
	  gStyle->SetEndErrorSize(8);
	  plotTheo->Draw("hist");
	  // get unbinned Madgraph theory curve
	  // add muon and electron channel to
	  // minimize statistical fluctuations
	  TH1F* unbinnedTheory = getTheoryPrediction("analyzeTopPartonLevelKinematics"+PS+"/"+plotName,"/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName+"/"+TopFilename(kSig, 0, "muon"));
	  unbinnedTheory->Add(getTheoryPrediction("analyzeTopPartonLevelKinematics"+PS+"/"+plotName,"/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName+"/"+TopFilename(kSig, 0, "electron")));
	  // get MC@NLO curve
	  TString plotName2="";
	  if     (plotName=="topPt"    ) plotName2="hVisTopPt";
	  else if(plotName=="topY"     ) plotName2="hVisTopY";
	  else if(plotName=="ttbarPt"  ) plotName2="hVisTTbarPt";
	  else if(plotName=="ttbarY"   ) plotName2="hVisTTbarY";
	  else if(plotName=="ttbarMass") plotName2="hVisTTbarM";
	  else if(plotName=="lepPt"    ) plotName2="hVisLepPt";
	  else if(plotName=="lepEta"   ) plotName2="hVisLepEta";
	  else{ 
	    std::cout << "no valid name for input variable ";
	    std::cout << plotName << " for MC@Nlo chosen" << std::endl;
	    exit(0);
	  }
	  TH1F* unbinnedTheoryMCAtNLO = getTheoryPrediction(plotName2,"/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName+"/ttbarNtupleCteq6m.root");
	  TH1F* unbinnedTheory2 = (TH1F*)unbinnedTheory->Clone();
	  TH1F* unbinnedTheoryMCAtNLO2 = (TH1F*)unbinnedTheoryMCAtNLO->Clone();
          // get histogram central value & up/down variations for MC@NLO error bands
          TString plotName2_Up   = plotName2; plotName2_Up  .Append("_Up"  ) ;
          TString plotName2_Down = plotName2; plotName2_Down.Append("_Down") ;
	  TH1F* errorTheoryMCatNLOcentral = getTheoryPrediction(plotName2     ,"/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName+"/ttbarNtupleCteq6m.root");
	  TH1F* errorTheoryMCatNLOup      = getTheoryPrediction(plotName2_Up  ,"/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName+"/ttbarNtupleCteq6m.root");
	  TH1F* errorTheoryMCatNLOdown    = getTheoryPrediction(plotName2_Down,"/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName+"/ttbarNtupleCteq6m.root");
	  // normalize to unsit area for diff. norm. plots
	  if(xSecVariables_[i].Contains("Norm")){
	    if(!xSecVariables_[i].Contains("lep")&&!xSecVariables_[i].Contains("Y")&&!xSecVariables_[i].Contains("ttbarPt")&&!xSecVariables_[i].Contains("ttbarMass")) unbinnedTheory->Rebin(10);
	    if(xSecVariables_[i].Contains("ttbarMass")) unbinnedTheory->Rebin(6);
	    unbinnedTheory ->Scale(1.0/(unbinnedTheory ->Integral(0,unbinnedTheory->GetNbinsX()+1)));
	    unbinnedTheory ->Scale(1.0/(unbinnedTheory ->GetBinWidth(1)));
	    unbinnedTheory2->Scale(1.0/(unbinnedTheory2->Integral(0,unbinnedTheory2->GetNbinsX()+1)));
	    unbinnedTheory2->Scale(1.0/(unbinnedTheory2->GetBinWidth(1)));
	    if(xSecVariables_[i].Contains("Y")) unbinnedTheoryMCAtNLO->Rebin(20);
	    else if(xSecVariables_[i].Contains("ttbarPt"))unbinnedTheoryMCAtNLO->Rebin(2);
	    else if(xSecVariables_[i].Contains("Mass"))unbinnedTheoryMCAtNLO->Rebin(4);
	    else if(xSecVariables_[i].Contains("Mass"))unbinnedTheoryMCAtNLO->Rebin(2);
	    else if(!xSecVariables_[i].Contains("lep"))unbinnedTheoryMCAtNLO->Rebin(2);
	    else if(xSecVariables_[i].Contains("lepPt")) unbinnedTheoryMCAtNLO->Rebin(1);
	    else unbinnedTheoryMCAtNLO->Rebin(10);
	    unbinnedTheoryMCAtNLO->Scale(1.0/(unbinnedTheoryMCAtNLO->Integral(0,unbinnedTheoryMCAtNLO->GetNbinsX()+1)));
	    unbinnedTheoryMCAtNLO->Scale(1.0/(unbinnedTheoryMCAtNLO->GetBinWidth(1)));
	    unbinnedTheoryMCAtNLO2->Scale(1.0/(unbinnedTheoryMCAtNLO2->Integral(0,unbinnedTheoryMCAtNLO2->GetNbinsX()+1)));
	    unbinnedTheoryMCAtNLO2->Scale(1.0/(unbinnedTheoryMCAtNLO2->GetBinWidth(1)));
	  }
	  // take into account that e and mu channel were added 
	  else{
	    unbinnedTheory->Scale(0.5);
	    unbinnedTheory->Scale(lumiweight(kSig, luminosity, 0, "muon"));
	    unbinnedTheory->Scale(1/luminosity);
	    unbinnedTheory->Scale(1.0/(unbinnedTheory->GetBinWidth(1)));
	    unbinnedTheory2->Scale(0.5);
	    unbinnedTheory2->Scale(lumiweight(kSig, luminosity, 0, "muon"));
	    unbinnedTheory2->Scale(1.0/(unbinnedTheory2->GetBinWidth(1)));
	    unbinnedTheoryMCAtNLO->Scale(1.0/(unbinnedTheoryMCAtNLO->GetBinWidth(1)));
	    // cross section / Ngen
	    unbinnedTheoryMCAtNLO->Scale(ttbarCrossSection/10000);
	    unbinnedTheoryMCAtNLO2->Scale(1.0/(unbinnedTheoryMCAtNLO2->GetBinWidth(1)));
	    // cross section / Ngen
	    unbinnedTheoryMCAtNLO2->Scale(ttbarCrossSection/10000);

	  }

          // uncertainty bands MC@NLO
          Int_t errorRebin  = 0 ;
          Int_t errorSmooth = 0 ;
          if(xSecVariables_[i].Contains("lepPt"  )) { errorRebin =  0 ; errorSmooth =   10 ; }
          if(xSecVariables_[i].Contains("lepEta" )) { errorRebin =  5 ; errorSmooth =   10 ; }
          if(xSecVariables_[i].Contains("topPt"  )) { errorRebin = 10 ; errorSmooth =   10 ; }
          if(xSecVariables_[i].Contains("topY"   )) { errorRebin =  5 ; errorSmooth =   10 ; }
          if(xSecVariables_[i].Contains("ttbarPt")) { errorRebin =  2 ; errorSmooth =    5 ; }
          if(xSecVariables_[i].Contains("ttbarY" )) { errorRebin =  5 ; errorSmooth =   10 ; }
          if(xSecVariables_[i].Contains("ttbarM" )) { errorRebin =  5 ; errorSmooth =   10 ; }
 
          if (errorRebin) {
            errorTheoryMCatNLOcentral->Rebin(errorRebin) ;
            errorTheoryMCatNLOup     ->Rebin(errorRebin) ;
            errorTheoryMCatNLOdown   ->Rebin(errorRebin) ;
          }
        
          errorTheoryMCatNLOcentral->Scale(1./errorTheoryMCatNLOcentral->GetBinWidth(1)) ;
          errorTheoryMCatNLOup     ->Scale(1./errorTheoryMCatNLOup     ->GetBinWidth(1)) ;
          errorTheoryMCatNLOdown   ->Scale(1./errorTheoryMCatNLOdown   ->GetBinWidth(1)) ;
    
          // these are used to check if smoothing does not distort distributions too much  
	  //TH1F* errorTheoryMCatNLOcentralClean = (TH1F*) errorTheoryMCatNLOcentral->Clone() ;
	  //TH1F* errorTheoryMCatNLOupClean      = (TH1F*) errorTheoryMCatNLOup     ->Clone() ;
	  //TH1F* errorTheoryMCatNLOdownClean    = (TH1F*) errorTheoryMCatNLOdown   ->Clone() ;

          if (errorSmooth) {
            errorTheoryMCatNLOcentral->Smooth(errorSmooth);
            errorTheoryMCatNLOup     ->Smooth(errorSmooth);
            errorTheoryMCatNLOdown   ->Smooth(errorSmooth);
          }

          // with a TGraphAsymmErrors we can do the nice error band plotting without running in trouble when overlaying multiple histograms
          TGraphAsymmErrors * errorBandsMCatNLO = new TGraphAsymmErrors(errorTheoryMCatNLOcentral->GetNbinsX()-1);
          for (Int_t iBin = 1; iBin < errorTheoryMCatNLOcentral->GetNbinsX() ; iBin++ ) {

            Double_t centralValue = errorTheoryMCatNLOcentral->GetBinContent(iBin) ;
            Double_t     maxValue = errorTheoryMCatNLOup     ->GetBinContent(iBin) ;
            Double_t     minValue = errorTheoryMCatNLOdown   ->GetBinContent(iBin) ;

            //errorBandsMCatNLO->SetPoint        ( iBin, errorTheoryMCatNLOcentral->GetBinCenter(iBin), centralValue  );
            //errorBandsMCatNLO->SetPointEXlow   ( iBin, errorTheoryMCatNLOcentral->GetXaxis()->GetBinLowEdge(iBin)   );
            //errorBandsMCatNLO->SetPointEXhigh  ( iBin, errorTheoryMCatNLOcentral->GetXaxis()->GetBinUpEdge (iBin)   );
	    // GOSSIE quick fix
	    if (cutTtbarMass&&xSecVariables_[i].Contains("ttbarM" ) &&
		errorTheoryMCatNLOcentral->GetBinCenter(iBin) < 345.) {
	      errorBandsMCatNLO->SetPoint        ( iBin, errorTheoryMCatNLOcentral->GetBinCenter(iBin), 0.            );
	      errorBandsMCatNLO->SetPointEXlow   ( iBin, 0.);
	      errorBandsMCatNLO->SetPointEXhigh  ( iBin, 0.);
	    } 
	    else {
	      errorBandsMCatNLO->SetPoint        ( iBin, errorTheoryMCatNLOcentral->GetBinCenter(iBin), centralValue  );
	      errorBandsMCatNLO->SetPointEXlow   ( iBin, errorTheoryMCatNLOcentral->GetXaxis()->GetBinLowEdge(iBin)   );
	      errorBandsMCatNLO->SetPointEXhigh  ( iBin, errorTheoryMCatNLOcentral->GetXaxis()->GetBinUpEdge (iBin)   );
            }

            if (maxValue > minValue) {
              errorBandsMCatNLO->SetPointEYhigh( iBin, maxValue     - centralValue                                  );
              errorBandsMCatNLO->SetPointEYlow ( iBin, centralValue - minValue                                      );
            } else {
              errorBandsMCatNLO->SetPointEYhigh( iBin, minValue     - centralValue                                  );
              errorBandsMCatNLO->SetPointEYlow ( iBin, centralValue - maxValue                                      );
            }
          }

	  histogramStyle(*unbinnedTheory        , kSig, false, 1.2, kRed+1); 
	  histogramStyle(*unbinnedTheory2       , kSig, false, 1.2, kRed+1);
	  histogramStyle(*unbinnedTheoryMCAtNLO , kSig, false, 1.2, kAzure);
	  histogramStyle(*unbinnedTheoryMCAtNLO2, kSig, false, 1.2, kAzure);
	  if(xSecVariables_[i].Contains("Y"      )) unbinnedTheory->Smooth(2);
	  else if(xSecVariables_[i].Contains("ttbarPt")) unbinnedTheory->Smooth(10);
	  else unbinnedTheory->Smooth(10);
	  if(xSecVariables_[i].Contains("Mass")) unbinnedTheoryMCAtNLO->Smooth(50);
	  else if(!xSecVariables_[i].Contains("Y")) unbinnedTheoryMCAtNLO->Smooth(10);
	  else unbinnedTheoryMCAtNLO->Smooth(3);
  if(xSecVariables_[i].Contains("Norm")){
            // adapt style for error bands to central curve mc@nlo
            errorBandsMCatNLO->SetFillColor  (kGray)    ;
            errorBandsMCatNLO->SetFillStyle  (1001)     ; // NB: explicitly needed, otherwise filling invisible due to default "0"
            errorBandsMCatNLO->SetLineColor  (unbinnedTheoryMCAtNLO->GetLineColor()) ;
            errorBandsMCatNLO->SetLineWidth  (unbinnedTheoryMCAtNLO->GetLineWidth()) ;
            errorBandsMCatNLO->SetMarkerColor(unbinnedTheoryMCAtNLO->GetLineColor()) ;
            errorBandsMCatNLO->Draw("e3 same")          ; 
  	    // draw unsmoothed theory curves
	    // to see if rebinning/smoothing changes the shape
	    //unbinnedTheoryMCAtNLO2->Draw("same");
	    //unbinnedTheory2->Draw("same");
	    //errorTheoryMCatNLOcentralClean->Draw("hist c same");
	    //errorTheoryMCatNLOupClean     ->Draw("hist c same");
	    //errorTheoryMCatNLOdownClean   ->Draw("hist c same");
	    // draw smoothed theory curves
	    plotTheo->Draw("hist same")                 ;
	    unbinnedTheoryMCAtNLO->Draw("hist c same")  ;
	    // draw reweighted histo for closure test
	    if(reweightClosure&&sys==sysNo){ 
	      histo_["reweighted"+plotName][kSig]->Draw("hist same");
	    }
	    // draw distorted parton truth histo including zprime
	    if(zprime!=""&&sys==sysNo){
	      histo_["modified"+plotName][kSig]->Draw("hist same");
	    }
	    // GOSSIE quick fix
	    if(cutTtbarMass&&xSecVariables_[i].Contains("ttbarMass")) {
	      unbinnedTheory->GetXaxis()->SetRangeUser(325.,1200.);
	    }
	    unbinnedTheory       ->Draw("c same")       ;
	    TLegend *leg = new TLegend();
	    leg->SetX1NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength()-0.25);
	    leg->SetY1NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength()-0.15);
	    leg->SetX2NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength());
	    leg->SetY2NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength());
	    leg->SetTextFont(42);
	    leg->SetTextSize(0.04);
	    leg->SetFillStyle(0);
	    leg->SetBorderSize(0);
	    leg->SetTextAlign(12);
	    //leg->SetHeader("");
	    leg->AddEntry(plotCombination  , dataLabel , "P" );
	    leg->AddEntry(unbinnedTheory   , "MadGraph", "L" );
	    // add additional labels for closure test(s)
	    if(reweightClosure&&sys==sysNo) leg->AddEntry(histo_["reweighted"+plotName][kSig], "MadGraph reweighted", "L");
	    if(zprime!=""&&sys==sysNo)  leg->AddEntry(histo_["modified"+plotName][kSig], "t#bar{t} & "+zprime+" GeV Z'", "L");
	    leg->AddEntry(errorBandsMCatNLO, "MC@NLO  ", "FL");
	    leg->Draw("same");
	  }
	  plotCombination->Draw("e1 same");
          gPad->RedrawAxis(); 
	  DrawCMSLabels(true,luminosity);
	  DrawDecayChLabel("e/#mu + Jets Combined");
	  histo_[xSecVariables_[i]][sys]=(TH1F*)(plotCombination->Clone());
	  // save combined e+mu plot for systematic error calculation afterwards
	  canvas_[xSecVariables_[i]][sys]=(TCanvas*)(combicanvas->Clone());
	  canvas_[xSecVariables_[i]][sys]->SetTitle(xSecVariables_[i]);
	  canvas_[xSecVariables_[i]][sys]->SetName(xSecVariables_[i]);
	  if(save&&sys==sysNo){
	    int initialIgnoreLevel=gErrorIgnoreLevel;
	    if(verbose==0) gErrorIgnoreLevel=kWarning;
	    combicanvas->Print("./diffXSecFromSignal/plots/combined/"+dataSample+"/xSec/xSec"+xSecVariables_[i]+".eps"); 
	    combicanvas->Print("./diffXSecFromSignal/plots/combined/"+dataSample+"/xSec/xSec"+xSecVariables_[i]+".png");
	    gErrorIgnoreLevel=initialIgnoreLevel;
	  }
	  // close Canvas
	  combicanvas->Close();
	  //delete combicanvas;
	  delete combicanvas;
	}
	if(!(plotMu||plotEl)&&verbose>1){ 
	  std::cout << "ERROR: plot " << xSecVariables_[i]+"kData" << " not found in ";
	  std::cout << xSecFolder+"/"+subfolder+"/"+xSecVariables_[i];
	  std::cout << " for decay channel ";
	  if(!plotMu) std::cout << "muon";
	  if(!plotEl) std::cout << " electron";
	  std::cout << std::endl;
	}
      }
      if(!(canvasMu||canvasEl)&&verbose>1) std::cout << "ERROR: canvas " << xSecFolder+"/"+subfolder+"/"+xSecVariables_[i] << " not found in either one or both decay channels!" << std::endl;
    }
  }
  // ---
  //    save combined #mu+e plots to an .root file
  // ---
  // loop variables
  for(unsigned int i=0; i<xSecVariables_.size(); ++i){
    // loop systematic variations
    for(unsigned int sys=sysNo; sys<=sysDiBosDown; ++sys){
      if(histo_[xSecVariables_[i]][sys]){
	saveToRootFile("diffXSecTopSemiLep.root", canvas_[xSecVariables_[i]][sys], true, verbose,"xSec/"+sysLabel(sys));
      }
    }
  }
}

