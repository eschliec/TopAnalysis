#include "basicFunctions.h"
#include "BCC.h"

void combineTopDiffXSecUncertainties(double luminosity=1143, bool save=true, unsigned int verbose=1, TString inputFolderName="TOP2011/110819_AnalysisRun", TString decayChannel="muon", bool adpatOldUncertainties=true){
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
  
  TStyle myStyle("HHStyle","HHStyle");
  setHHStyle(myStyle);
  TGaxis::SetMaxDigits(2);
  myStyle.cd();
  gROOT->SetStyle("HHStyle");

  // ---
  //    parameter Configuration
  // ---
  // set detail level of ouStput 
  // 0: no output, 1: std output 2: output for debugging
  // save: enable saving
  // dataSample: see if its "2010" or "2011" data
  TString dataSample="2011";
  if(luminosity<50.) dataSample="2010";
  // adpatOldUncertainties: choose if you want to apapt 
  // uncertainties from 
  // the 2010 mu+jets analysis if they are not available 
  // for 2011 analysis, NOTE: these uncertainties 
  // must exist in ./diffXSecTopSemiMu2010.root
  // create container to save relative shifted plots
  std::map< TString, std::map <unsigned int, TCanvas*> > adoptedUncertaintyPlot_;
  // outputFile: target rootfile
  // NOTE: this must be identical with TString outputFileName 
  // in analyzeHypothesisKinFit.C
  TString outputFile="diffXSecTopSemi";
  if(decayChannel=="muon"    ) outputFile+="Mu"+dataSample;
  if(decayChannel=="electron") outputFile+="Elec"+dataSample;
  if(decayChannel=="combined"){
    outputFile+="Lep";
    adpatOldUncertainties=false;
  }
  outputFile+=".root";
  // define folder where XSec plots are stored
  TString xSecFolder = "xSec";
  // save all plots into the following folder
  TString outputFolder = "./diffXSecFromSignal/plots/"+decayChannel+"/";
  if(dataSample!="") outputFolder+=dataSample;
  
  // ---
  //    basic printout and variable definitions
  // ---
  // define variables for combination
  // NOTE: these must be identical to those defined in 
  // xSecVariables_ in analyzeHypothesisKinFit.C
  std::vector<TString> xSecVariables_;
  TString xSecVariables[] ={"topPt", "topY", "ttbarPt", "ttbarMass", "ttbarY", "lepPt", "lepEta", "topPtNorm", "topYNorm", "ttbarPtNorm", "ttbarMassNorm", "ttbarYNorm", "lepPtNorm", "lepEtaNorm", "inclusive"};
  xSecVariables_.insert( xSecVariables_.begin(), xSecVariables, xSecVariables + sizeof(xSecVariables)/sizeof(TString) );
  // chose min/max value[%] for relative uncertainty plots
  double errMax=40.;
  double errMin=-1*errMax;
  // create plot container
  std::map< TString, std::map <unsigned int, TH1F*> > histo_;
  std::map< TString, std::map <unsigned int, TH1F*> > relativeUncertainties_;
  // create container to indicate that plots have been found 
  // and therefor systematics will be calculated
  std::map<TString, std::map<unsigned int, bool> > calculateError_;
  // create container for combined Errors
  std::map<TString, TGraphAsymmErrors*> totalErrors_;
  // setup for bin center corrected values
  gSystem->Load("BCC_C.so");
  bool mergeLepAndHadTop=true;
  std::vector<TString> xSecVariableBranchNames_;
  for(unsigned int i=0; i<xSecVariables_.size(); ++i){
    if(!xSecVariables_[i].Contains("Norm")&&xSecVariables_[i]!="inclusive"){
      if(xSecVariables_[i].Contains("top")){
	xSecVariableBranchNames_.push_back(xSecVariables_[i]+"Had");
	xSecVariableBranchNames_.push_back(xSecVariables_[i]+"Lep");
      }
      else xSecVariableBranchNames_.push_back(xSecVariables_[i]);
    }
  }
  // parameter printout
  if(verbose>0) std::cout << std::endl << "executing combineTopDiffXSecUncertainties with " << dataSample << " data" << std::endl << std::endl;
  if(verbose>0) std::cout << "target file: " << outputFile << std::endl;
  if(verbose>0) std::cout << "target variables for cross sections:" << std::endl;
  // loop variables
  for(unsigned int i=0; i<xSecVariables_.size(); ++i){
    if(verbose>0){
      std::cout << xSecVariables_[i]; 
      if(i<xSecVariables_.size()-1) std::cout << ", ";
      else std::cout << std::endl;
    }
  }
  if(verbose>0) std::cout << "target folder containing cross section plots: " << xSecFolder << std::endl;
  if(verbose>0&&save) std::cout << "final plots will be saved in " << outputFile << " and as .eps in " << outputFolder << std::endl;
  // loading bin center corrections
  std::cout << "loading bin center corrections" << std::endl;
  BCC b("/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName+"/"+TopFilename(kSig, 0, (std::string)decayChannel),"analyzeTopPartonLevelKinematicsPhaseSpace",xSecVariableBranchNames_,mergeLepAndHadTop);
  b.runBCCCalculation();
  std::map<TString, std::vector<double> > correctedCenters_ = b.getMapWithCorrectedCentersInX();
  std::map<TString, std::vector<double> > corrCenterErrors_ = b.getMapWithCenterErrorsInX();
 // output results
  if(verbose>1){
    for (std::map<TString, std::vector<double> >::iterator iter1 = correctedCenters_.begin(); iter1 !=  correctedCenters_.end(); iter1++ )
      {
      std::cout << iter1->first << ": ";
      for (std::vector<double>::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++)
	{
	  std::cout << (*iter2) << " ";
	}
      std::cout << std::endl;
    }
    
    for (std::map<TString, std::vector<double> >::iterator iter1 = corrCenterErrors_.begin(); iter1 != corrCenterErrors_.end(); iter1++ )
      {
	std::cout << iter1->first << ": ";
	
	for (std::vector<double>::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++)
	  {
	    std::cout << (*iter2) << " ";
	  }
	std::cout << std::endl;
      }
  }

  // ---
  //    open rootfile
  // ---
  TFile* file = TFile::Open(outputFile, "UPDATE");
  if(verbose>1) std::cout << "opened file "+outputFile << std::endl;
  // open mu+jets 2010 analysis file
  TFile* oldAnalysisfile=0;
  if(adpatOldUncertainties){ 
    oldAnalysisfile = TFile::Open("/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName+"/diffXSecTopSemiMu2010.root", "READ");
    if(verbose>1) std::cout << "opened old file diffXSecTopSemiMu2010.root" << std::endl;
  }
  // check if file exist
  if(file&&!file->IsZombie()){
    if(verbose>0) std::cout << "target file found!" << std::endl;
    // ---
    //    part A: get all data xSec with all systematic variations
    // ---
    if(verbose>0) std::cout << std::endl << "part A: loading plots" << std::endl;
    // loop systematic variations
    for(unsigned int sys=sysNo; sys<=sysDiBosDown; ++sys){
      TString subfolder=sysLabel(sys);
      // loop variables
      for(unsigned int i=0; i<xSecVariables_.size(); ++i){
	// get canvas       
	TCanvas* canvas = (TCanvas*)file->Get(xSecFolder+"/"+subfolder+"/"+xSecVariables_[i]);
	if(canvas){
	  // get data plot for all systematics
	  TString plotName = xSecVariables_[i];
	  plotName.ReplaceAll("Norm","");
	  TH1F* plot= (TH1F*)canvas->GetPrimitive(plotName+"kData");
	  if(plot){ 
	    if(verbose>1) std::cout << "plot "+plotName+"kData in "+xSecFolder+"/"+subfolder+"/"+xSecVariables_[i] << " found!" << std::endl;
	    // go to root directory, keep plot when closing rootfile
	    gROOT->cd();
	    histo_[xSecVariables_[i]][sys]=(TH1F*)(plot->Clone());
	    calculateError_[xSecVariables_[i]][sys]=true;
	    // exclude "old" PU error
	    if(sys==sysPileUp) calculateError_[xSecVariables_[i]][sys]=false;
	    // exclude flat SF error
	    if(sys==sysflatTrigSF) calculateError_[xSecVariables_[i]][sys]=false;
	  }
	  if(!plot&&verbose>1){ 
	    std::cout << "ERROR: plot " << xSecVariables_[i]+"kData" << " not found in ";
	    std::cout << xSecFolder+"/"+subfolder+"/"+xSecVariables_[i] << std::endl;
	  }
	}	
	if(adpatOldUncertainties&&!canvas&&sys!=sysNo&&sys!=sysPileUp&&!(sys>=sysBtagSFUp&&sys<=sysMisTagSFdown)&&!xSecVariables_[i].Contains("lep")){
	  if(verbose>1) std::cout << "use 2010 uncertainties for variation " << sysLabel(sys)+", "+xSecVariables_[i] << std::endl;
	  // adapt 2010 UNCERTAINTIES
	  // if systematic variation does not exist for 2011 data (indicated by luminosity), 
	  // take the binwise relative uncertainty from 2010 mu+jets analysis (no PF2PAT)
	  // and save relative shifted plot for 2011 result
	  // exclude PU uncertainty
	  // exclude new uncertainties, that does not exist in 2010 analysis
	  if(luminosity>40.0){
	    if(verbose>1) std::cout << "luminosity>40/pb" << std::endl;
	    // copy std analysis plot to scale bins later with relative uncertainties 
	    TString plotName = xSecVariables_[i];
	    plotName.ReplaceAll("Norm","");
	    TH1F* currentShiftedCrossSectionPlot=(TH1F*)(histo_[xSecVariables_[i]][sysNo]->Clone(plotName+"kData"));
	    if(verbose>1) std::cout << "std analysis file for this variation copied" << std::endl;
	    // loop all bins of current 2011 analysis histogram without variations
	    for(int bin=1; bin<=histo_[xSecVariables_[i]][sysNo]->GetNbinsX(); ++bin){
	      if (xSecVariables_[i].Contains("ttbarMass") && bin== 1) continue;   // Ad-hoc bug fix for non-empty first bin of ttbarmass (Thomas)
	      if (xSecVariables_[i].Contains("topY")      && bin== 1) continue;   // Ad-hoc bug fix for non-empty first bin of ttbarmass (Martijn)
	      if (xSecVariables_[i].Contains("topY")      && bin==10) continue;   // Ad-hoc bug fix for non-empty first bin of ttbarmass (Martijn)
	      //ensure to have non-empty bins
	      if(histo_[xSecVariables_[i]][sysNo]->GetBinContent(bin)>0){
		if(verbose>1) std::cout << "bin " << bin << std::endl;
		// get plot with relative uncertainties for given quantity and bin
		TCanvas* oldCanvas = (TCanvas*)oldAnalysisfile->Get("relativeUncertainties/"+xSecVariables_[i]+"/relSysPlotBin"+getTStringFromInt(bin)+xSecVariables_[i]);
		if(verbose>1) std::cout << "got canvas" << std::endl; 
		TH1F* plot= (TH1F*)oldCanvas->GetPrimitive("relSysPlotBin"+getTStringFromInt(bin)+xSecVariables_[i]);
		if(plot){ 
		  if(verbose>1) std::cout << "got plot" << std::endl;
		  // take into account that some systematics in 2010
		  // were not present and therefore the last (diboson)
		  // bins are shifted
		  int sysold=sys;
		  if(sys==sysDiBosDown||sys==sysDiBosUp){
		    int NbinsOld=plot->GetNbinsX();
		    int NbinsNew=sysDiBosDown+5;
		    sysold-=(NbinsNew-NbinsOld);
		  }
		  double relSysErrorBin = plot->GetBinContent(sysold);
		  if(verbose>1) std::cout << "bin rel uncertainty: " << relSysErrorBin << std::endl;
		  // apply binwise SF for the relative uncertainty 
		  // from 2010 mu+jets analysis to the copy of the 
		  // std 2011 analysis plot
		  double binSFForSystUncertainty=(1.0+relSysErrorBin/100.);
		  if(verbose>1) std::cout << "binSF: " << binSFForSystUncertainty << std::endl;
		  currentShiftedCrossSectionPlot->SetBinContent(bin, binSFForSystUncertainty*currentShiftedCrossSectionPlot->GetBinContent(bin));
		  currentShiftedCrossSectionPlot->SetBinError  (bin, binSFForSystUncertainty*currentShiftedCrossSectionPlot->GetBinError(bin)  );
		  if(verbose>1) std::cout << "scaled xSec: " << currentShiftedCrossSectionPlot->GetBinContent(bin) << std::endl;
		  histo_[xSecVariables_[i]][sys]=(TH1F*)(currentShiftedCrossSectionPlot->Clone());
		  // save plot in canvas 
		  TCanvas* currentUncertainty = new TCanvas( xSecVariables_[i], xSecVariables_[i], 600, 600);
		  //canvasStyle(*currentUncertainty);
		  currentUncertainty->cd();
		  currentShiftedCrossSectionPlot->Draw("");
		  // draw label to indicate that sys error is adapted from 2010 mu+jets
		  DrawLabel("!adopted from 2010 #mu+jets!", 0.2, 0.7, 0.7, 0.8, 0.07);
		  // 
		  adoptedUncertaintyPlot_[xSecVariables_[i]][sys]=(TCanvas*)(currentUncertainty->Clone());
		  delete currentUncertainty;
		  // enable error calculation
		  if(!calculateError_[xSecVariables_[i]].count(sys)>0) calculateError_[xSecVariables_[i]][sys]=true;
		}
	      }
	    }
	    if(verbose>1) std::cout << "scaled xSec first bin: " << histo_[xSecVariables_[i]][sys]->GetBinContent(1) << std::endl;
	  }
	}
	if(!canvas&&verbose>1) std::cout << "ERROR: canvas " << xSecFolder+"/"+subfolder+"/"+xSecVariables_[i] << " not found!" << std::endl;
      }
    }
    
    // ---
    //    part B:calculate combined systematic+statistic error
    // ---
    if(verbose>0) std::cout << std::endl << "part B: calculate systematic errors for all variables and bins" << std::endl;
    // loop variables
    for(unsigned int i=0; i<xSecVariables_.size(); ++i){
      if(verbose>0) std::cout << std::endl << "variable: " << xSecVariables_[i] << std::endl;
      // check if any plot of the chosen variable has been found
      if(calculateError_.count(xSecVariables_[i])>0&&calculateError_[xSecVariables_[i]][sysNo]==true){
	// define object to save asymmetric errors
	int Nbins=histo_[xSecVariables_[i]][sysNo]->GetNbinsX();
	TGraphAsymmErrors* combinedErrors= new TGraphAsymmErrors(Nbins);
	// get plot without variation
	TH1F* noSysPlot=(TH1F*)histo_[xSecVariables_[i]][sysNo]->Clone();
	// loop bins
	for(int bin=1; bin<=noSysPlot->GetNbinsX(); ++bin){
	  double totalSystematicErrorUp  =0;
	  double totalSystematicErrorDown=0;
	  double stdBinXSecValue=histo_[xSecVariables_[i]][sysNo]->GetBinContent(bin);
	  // jump to bin without vanishing bin content
	  while(stdBinXSecValue==0&&bin<Nbins){
	    if(verbose>1) std::cout << "bin #" << bin << "/" << Nbins << " has bin content ==0, will be skipped!" << std::endl;
	    ++bin;
	    stdBinXSecValue=histo_[xSecVariables_[i]][sysNo]->GetBinContent(bin);
	  }
	  // if last bin is empty -> stop looping, next variable
	  if(bin==histo_[xSecVariables_[i]][sysNo]->GetNbinsX()&&stdBinXSecValue==0){
	    if(verbose>1){
	      std::cout << "bin #" << bin << "/" << Nbins << " has bin content ==0, will be skipped!" << std::endl;
	      std::cout << "last bin reached" << std::endl;
	    }
	  }
	  // else: continue with calculation
	  else{
	    if(verbose>0){ 
	      std::cout << std::endl << xSecVariables_[i] << " bin #" << bin << "/" << Nbins;
	      std::cout << " ( std value[pb/binwidth]: " << stdBinXSecValue << " )" << std::endl;
	    }
	    // create plot that indicates all relative uncertainties
	    // one plot for every variable in every bin
	    TH1F* relSysPlot = new TH1F("relSysPlotBin"+getTStringFromInt(bin)+xSecVariables_[i], "relSysPlotBin"+getTStringFromInt(bin)+xSecVariables_[i], sysDiBosDown+5, sysNo, sysDiBosDown+5);
	    relSysPlot->GetXaxis()->SetLabelSize(0.025);
	    // loop systematic variations
	    for(unsigned int sys=sysNo+1; sys<=sysDiBosDown; ++sys){
	      // set labels for relative uncertainties for bin & variable in map relativeUncertainties_
	      relSysPlot->GetXaxis()->SetBinLabel(sys, sysLabel(sys).ReplaceAll("sys",""));
	      // create plot that indicates the relative systematic uncertainty
	      double sysDiff=0;
	      if(verbose>0) std::cout << sysLabel(sys);
	      // check if chosen systematic variation of chosen variable has been found and std value is nonzero
	      if(calculateError_[xSecVariables_[i]].count(sys)<=0||stdBinXSecValue==0){
		if(verbose>0){
		  if(stdBinXSecValue!=0) std::cout << ": not found" << std::endl;
		  else std::cout << ": std value is 0" << std::endl;
		}
	      }
	      else{
		// check if variable has to be considerered for total systematic error
		if(calculateError_[xSecVariables_[i]][sys]==true){
		  if(verbose>0) std::cout << "(considered): ";
		  double sysBinXSecValue=histo_[xSecVariables_[i]][sys]->GetBinContent(bin);
		  sysDiff=sysBinXSecValue-stdBinXSecValue;
		  // set btagging error to 0 as long as it is not there
		  //if(sys==sysBtagDown||sys==sysBtagUp) sysDiff=0;
		  // asymmetric errors
		  //if(sysDiff>0) totalSystematicErrorUp  +=sysDiff*sysDiff;
		  //if(sysDiff<0) totalSystematicErrorDown+=sysDiff*sysDiff;
		  totalSystematicErrorUp  +=0.5*sysDiff*0.5*sysDiff;
		  totalSystematicErrorDown+=0.5*sysDiff*0.5*sysDiff;
		}
		else{ 
		  if(verbose>0) std::cout << "(not considered): ";
		}
		// print single systematic uncertainty absolut and relative for bin & variable
		if(verbose>0) std::cout << sysDiff << " ( = " << 100*sysDiff/stdBinXSecValue << "% )" << std::endl;
		// save relative systematic uncertainties for bin & variable
		relSysPlot->SetBinContent(sys, 100*sysDiff/stdBinXSecValue);
	      }
	      // for last systematic 
	      if(sys==sysDiBosDown){
		// go to root directory keep plot when closing rootfile
		gROOT->cd();
		// finally save relative uncertainties in map relativeUncertainties_
		relativeUncertainties_[xSecVariables_[i]][bin]=(TH1F*)relSysPlot->Clone();
		// print total systematic uncertainty absolut and relative for bin & variable
		if(verbose>0){
		  std::cout << "total systematic uncertainty: " << std::endl;
		  std::cout << " +" << sqrt(totalSystematicErrorUp)   << " ( =";
		  std::cout << 100*sqrt(totalSystematicErrorUp  )/stdBinXSecValue << "% )" << std::endl;
		  std::cout << " -" << sqrt(totalSystematicErrorDown) << " ( =";
		  std::cout << 100*sqrt(totalSystematicErrorDown)/stdBinXSecValue << "% )" << std::endl;
		}
		// get statistical error
		double statErrorBinVar = noSysPlot->GetBinError(bin);
		double combinedErrorUpBinVar   = sqrt(totalSystematicErrorUp  +statErrorBinVar*statErrorBinVar);
		double combinedErrorDownBinVar = sqrt(totalSystematicErrorDown+statErrorBinVar*statErrorBinVar);
		// print statistical uncertainty absolut and relative for bin & variable
		if(verbose>0){
		  std::cout << "total statistic uncertainty: " << std::endl;
		  std::cout << " +/- " << statErrorBinVar << " ( =";
		  std::cout << 100*statErrorBinVar/stdBinXSecValue << "% )" << std::endl;
		}
		if(verbose>0){
		  std::cout << "total uncertainty: " << std::endl;
		  std::cout << " +" <<  combinedErrorUpBinVar   << " ( =";
		  std::cout << 100*combinedErrorUpBinVar/stdBinXSecValue   << "% )" << std::endl;
		  std::cout << " -" <<  combinedErrorDownBinVar << " ( =";
		  std::cout << 100*combinedErrorDownBinVar/stdBinXSecValue << "% )" << std::endl;
		}
		// save relative uncertainties for bin & variable in map relativeUncertainties_
		// a) statistic uncertainty
		relativeUncertainties_[xSecVariables_[i]][bin]->SetBinContent(sysDiBosDown+1, 100*statErrorBinVar/stdBinXSecValue);
		relativeUncertainties_[xSecVariables_[i]][bin]->GetXaxis()->SetBinLabel(sysDiBosDown+1, "statistical");
		// b) total systematic uncertainty up
		relativeUncertainties_[xSecVariables_[i]][bin]->SetBinContent(sysDiBosDown+2, 100*sqrt(totalSystematicErrorUp  )/stdBinXSecValue);
		relativeUncertainties_[xSecVariables_[i]][bin]->GetXaxis()->SetBinLabel(sysDiBosDown+2, "total syst. Up");
		// c) total systematic uncertainty down
		relativeUncertainties_[xSecVariables_[i]][bin]->SetBinContent(sysDiBosDown+3, -100*sqrt(totalSystematicErrorDown)/stdBinXSecValue);
		relativeUncertainties_[xSecVariables_[i]][bin]->GetXaxis()->SetBinLabel(sysDiBosDown+3, "total syst. Down");
		// d) total stat+systematic uncertainty up
		relativeUncertainties_[xSecVariables_[i]][bin]->SetBinContent(sysDiBosDown+4, 100*combinedErrorUpBinVar/stdBinXSecValue);
		relativeUncertainties_[xSecVariables_[i]][bin]->GetXaxis()->SetBinLabel(sysDiBosDown+4, "total Up");
		// e) total stat+systematic uncertainty down
		relativeUncertainties_[xSecVariables_[i]][bin]->SetBinContent(sysDiBosDown+5, -100*combinedErrorDownBinVar/stdBinXSecValue);
		relativeUncertainties_[xSecVariables_[i]][bin]->GetXaxis()->SetBinLabel(sysDiBosDown+5, "total Down");
		// f) set 2010 PU uncertainty for 2011 analysis to 0
		if(calculateError_[xSecVariables_[i]][sysPileUp]==false) relativeUncertainties_[xSecVariables_[i]][bin]->SetBinContent(sysPileUp, 0);
		// set combined errors for final xSec plot
		double pointXValue = histo_[xSecVariables_[i]][sysNo]->GetBinCenter(bin);
		double pointXError = 0;
		if(decayChannel=="combined"&&xSecVariables_[i]!="inclusive"){
		  TString plotName=xSecVariables_[i];
		  plotName.ReplaceAll("Norm","");
		  pointXValue = correctedCenters_[plotName].at(bin-1);
		  pointXError = corrCenterErrors_[plotName].at(bin-1);
		  pointXError=0;
		}
		combinedErrors->SetPoint(0, 0, -1000);
		combinedErrors->SetPoint(bin, pointXValue, histo_[xSecVariables_[i]][sysNo]->GetBinContent(bin));
		combinedErrors->SetPointError(bin, pointXError, pointXError, combinedErrorDownBinVar, combinedErrorUpBinVar);
		// define style for relative error plots
		histogramStyle(*relativeUncertainties_[xSecVariables_[i]][bin], kSig, true, 2.0, kBlack); 
		relativeUncertainties_[xSecVariables_[i]][bin]->GetXaxis()->LabelsOption("v");
		relativeUncertainties_[xSecVariables_[i]][bin]->GetXaxis()->SetLabelSize(0.05);
		relativeUncertainties_[xSecVariables_[i]][bin]->SetMaximum(errMax);
		relativeUncertainties_[xSecVariables_[i]][bin]->SetMinimum(errMin);
		double histMax = relativeUncertainties_[xSecVariables_[i]][bin]->GetBinContent(relativeUncertainties_[xSecVariables_[i]][bin]->GetMaximumBin());
		double histMin = relativeUncertainties_[xSecVariables_[i]][bin]->GetBinContent(relativeUncertainties_[xSecVariables_[i]][bin]->GetMinimumBin());
		if(histMax>errMax||histMin<errMin){
		  double newRange=histMax;
		  if(fabs(histMax)<fabs(histMin)) newRange=-1*histMin;
		  relativeUncertainties_[xSecVariables_[i]][bin]->SetMaximum(1.2*newRange);
		  relativeUncertainties_[xSecVariables_[i]][bin]->SetMinimum(-1.2*newRange);
		}
		relativeUncertainties_[xSecVariables_[i]][bin]->SetStats(kFALSE);
		relativeUncertainties_[xSecVariables_[i]][bin]->GetYaxis()->SetTitle("rel. uncertainty (%)");
		relativeUncertainties_[xSecVariables_[i]][bin]->GetYaxis()->SetTitleOffset(1.2);
		// save asymmetric errors in map totalErrors_
		combinedErrors->SetLineWidth(1.0); // histo_[xSecVariables_[i]][sysNo]->GetLineWidth());
		combinedErrors->SetLineColor(histo_[xSecVariables_[i]][sysNo]->GetLineColor());
		totalErrors_[xSecVariables_[i]]=(TGraphAsymmErrors*)(combinedErrors->Clone());
		whipEmptyBinsAway(totalErrors_[xSecVariables_[i]], verbose);
		// remove ugly first point in m(ttbar)
		if(xSecVariables_[i].Contains("ttbarMass")){
		  totalErrors_[xSecVariables_[i]]->SetPoint(1, 0, -1000);
		}
	      }
	    }
	    delete relSysPlot;
	  }
	}
      }
    }
    // close file
    // needed to be able to use the saveToRootFile function 
    // which also opens the file
    file->Close();
    
    // ---
    //    part C:draw plot with combined systematic+statistic error
    //           in new cross section plot
    // ---
    if(save){
      if(verbose>0) std::cout << std::endl << "part C: save xSec plots with stat.+sys. errors and plots with relative uncertainties" << std::endl;
      // loop variables
      for(unsigned int i=0; i<xSecVariables_.size(); ++i){
	// check if total combined errors are available
	if(totalErrors_.count(xSecVariables_[i])>0){
	  if(verbose>0) std::cout << std::endl << "variable: " << xSecVariables_[i] << std::endl;
	  // save rel uncertainties plot for all variables and bins
	  // loop bins
	  for(int bin=0; bin<=histo_[xSecVariables_[i]][sysNo]->GetNbinsX(); ++bin){
	    // make sure that plot exists
	    if(relativeUncertainties_.count(xSecVariables_[i])>0&&relativeUncertainties_[xSecVariables_[i]].count(bin)>0){
	      // create canvas for plot
	      TCanvas* relUnCertaintyCanvas = new TCanvas(relativeUncertainties_[xSecVariables_[i]][bin]->GetName() ,relativeUncertainties_[xSecVariables_[i]][bin]->GetTitle(), 600, 600);
	      //canvasStyle(*relUnCertaintyCanvas);
	      relUnCertaintyCanvas->SetBottomMargin(0.30);
	      relUnCertaintyCanvas->SetLeftMargin(0.10);
	      relUnCertaintyCanvas->SetRightMargin(0.10);
	      relUnCertaintyCanvas->SetTopMargin(0.07);
	      relUnCertaintyCanvas->cd(0);
	      relUnCertaintyCanvas->SetGrid(1,1);
	      // draw plot into canvas
	      relativeUncertainties_[xSecVariables_[i]][bin]->Draw("hist");
	      // draw axis also on the right side of canvas
	      int xPosition=sysDiBosDown+5;
	      double histMax = relativeUncertainties_[xSecVariables_[i]][bin]->GetBinContent(relativeUncertainties_[xSecVariables_[i]][bin]->GetMaximumBin());
	      double histMin = relativeUncertainties_[xSecVariables_[i]][bin]->GetBinContent(relativeUncertainties_[xSecVariables_[i]][bin]->GetMinimumBin());
	      double max = errMax;
	      double min = errMin;
	      if(histMax>errMax||histMin<errMin){
		double newRange=1.2*histMax;
		if(fabs(histMax)<fabs(histMin)) newRange=-1.2*histMin;
		max=newRange;
		min=-1*newRange;
	      }
	      TGaxis *axis = new TGaxis(xPosition,-1*max,xPosition,max,min,max,relativeUncertainties_[xSecVariables_[i]][bin]->GetYaxis()->GetNdivisions(),"+L");
	      axis->SetLabelSize(myStyle.GetLabelSize());
	      axis->SetLabelFont(myStyle.GetLabelFont());
	      axis->SetLabelOffset(myStyle.GetLabelOffset());
	      axis->Draw("same");
	      // redraw to have statistical error as +/-
	      TH1F* relUnCertaintyCopy = (TH1F*)relativeUncertainties_[xSecVariables_[i]][bin]->Clone();
	      relUnCertaintyCopy->SetBinContent(sysDiBosDown+1, (-1.)*(relUnCertaintyCopy->GetBinContent(sysDiBosDown+1)));
	      relUnCertaintyCopy->DrawClone("hist same");
	      // draw every systematic variation with different color
	      int colourCounter=1;
	      unsigned int colour=kBlack;
	      for(unsigned int sys=1; sys<=(unsigned int)relativeUncertainties_[xSecVariables_[i]][bin]->GetNbinsX(); ++sys, ++colourCounter){
		if(sys==sysPileUp  ) ++colourCounter;
		if(colourCounter==5) colourCounter=1;
		if(colourCounter==1||colourCounter==2) colour=kBlack;
		if(colourCounter==3||colourCounter==4) colour=kRed;
		relUnCertaintyCopy->Scale(0.);
		relUnCertaintyCopy->SetBinContent(sys, relativeUncertainties_[xSecVariables_[i]][bin]->GetBinContent(sys));
		if(sys<=sysDiBosDown){
		  relUnCertaintyCopy->SetFillColor(colour);
		  relUnCertaintyCopy->DrawCopy("hist same");
		}
		if(sys>=sysDiBosDown+2&&sys<=sysDiBosDown+3){
		  relUnCertaintyCopy->SetFillColor(kOrange);
		  relUnCertaintyCopy->DrawCopy("hist same");
		}
		if(sys>=sysDiBosDown+4){
		  relUnCertaintyCopy->SetFillColor(kBlue);
		  relUnCertaintyCopy->DrawCopy("hist same");
		}
	      }
	      // save canvas to file
	      if(save) saveToRootFile( outputFile, relUnCertaintyCanvas, true, verbose, "relativeUncertainties/"+xSecVariables_[i]);
	      // save canvas as eps
	      int initialIgnoreLevel=gErrorIgnoreLevel;
	      if(verbose==0) gErrorIgnoreLevel=kWarning;
	      if(save){
		TString saveName=outputFolder+"/uncertainties/relativeUncertainties"+xSecVariables_[i]+"Bin"+getTStringFromInt(bin);
		if(decayChannel=="combined") saveName+="Combined";
		relUnCertaintyCanvas->Print(saveName+".eps");
		relUnCertaintyCanvas->Print(saveName+".png");
	      }
	      gErrorIgnoreLevel=initialIgnoreLevel;
	      // delete canvas
	      delete relUnCertaintyCanvas;
	      delete relUnCertaintyCopy;
	      delete axis;
	    }
	  }
	  // open file to get xSec canvas from it
	  file = TFile::Open(outputFile, "UPDATE");
	  // get Canvas
	  TCanvas* canvas = (TCanvas*)(file->Get(xSecFolder+"/sysNo/"+xSecVariables_[i])->Clone());
	  // close file
	  // needed to be able to use the saveToRootFile function 
	  // which also opens the file
	  file->Close();
	  if(canvas){
	    // for inclusive xSecs
	    if(xSecVariables_[i].Contains("inclusive")){
	      TCanvas* canvas2 = new TCanvas(xSecVariables_[i], xSecVariables_[i], 600, 600);
	      //canvasStyle(*canvas2);
	      canvas2->SetLeftMargin(0.05);
	      canvas2->cd();
	      canvas2->SetTitle(xSecVariables_[i]);
	      //canvas2->SetGrid(1,0);
	      // draw xSec 
	      int NxSecBins=3000;
	      double xSecMax=300.0;
	      double xSecMin=0.;
	      TH1F* finalInclusiveXSec = new TH1F( "finalInclusive", "finalInclusive", NxSecBins, xSecMin, xSecMax);
	      // x axis as xSec[pb] indicator
	      int xSecBin=roundToInt(histo_["inclusive"][sysNo]->GetBinContent(1)*((double)NxSecBins)/xSecMax)+1;
	      finalInclusiveXSec->SetBinContent(xSecBin, 1);
	      finalInclusiveXSec->SetBinError(xSecBin, 0);
	      histogramStyle(*finalInclusiveXSec, kData , false);
	      axesStyle(*finalInclusiveXSec, "#sigma(t#bar{t}#rightarrowX) [pb]", " ", 0, 4);
	      finalInclusiveXSec->GetYaxis()->SetNdivisions(0);
	      finalInclusiveXSec->GetXaxis()->SetNdivisions(510);
	      finalInclusiveXSec->GetXaxis()->SetNoExponent(true);
	      finalInclusiveXSec->GetYaxis()->SetTitleSize(0.0);							
	      finalInclusiveXSec->Draw("axis");
	      // draw Theory expectation     
	      double theoryXSecNLL=165;
	      double theoryXSecNLO=157.5;
	      double theoryErrorNLOUp  =23.2;
	      double theoryErrorNLODown=24.4;
	      double theoryErrorNLLUp  =sqrt(5*5+9*9);
	      double theoryErrorNLLDown=sqrt(9*9+9*9);
	      TBox* TheoryError = new TBox(theoryXSecNLO-theoryErrorNLODown, 0.1, theoryXSecNLO+theoryErrorNLOUp, 3.9);
	      TBox* TheoryError2= new TBox(theoryXSecNLL-theoryErrorNLLDown, 0.2 ,theoryXSecNLL+theoryErrorNLLUp, 3.8);
	      TheoryError->SetFillColor(kGray);
	      TheoryError2->SetFillColor(kGray+1);
	      TheoryError->Draw ("same");
	      TheoryError2->Draw("same");
	      // our Analysis result
	      double xSecError=histo_["inclusive"][sysNo]->GetBinError(1);
	      double xSecValue=histo_["inclusive"][sysNo]->GetBinContent(1);
	      double totalErrorUp  = totalErrors_[xSecVariables_[i]]->GetErrorYhigh(1);
	      double totalErrorDown= totalErrors_[xSecVariables_[i]]->GetErrorYlow(1);
	      // 	      if(verbose>0){
	      // 		std::cout << "xSec: "        << xSecValue << std::endl;
	      // 		std::cout << "stat. error: +/-" << xSecError << std::endl;
	      // 		std::cout << "tot. error up  : " << totalErrorUp   << std::endl;	      
	      // 		std::cout << "tot. error down: " << totalErrorDown << std::endl;
	      // 	      }
	      totalErrors_[xSecVariables_[i]+"ANinclusive"] = new TGraphAsymmErrors(1);
	      totalErrors_[xSecVariables_[i]+"ANinclusive"]->SetLineWidth(3);
	      totalErrors_[xSecVariables_[i]+"ANinclusive"]->SetMarkerSize(2.2);
	      totalErrors_[xSecVariables_[i]+"ANinclusive"]->SetMarkerStyle(20);
	      totalErrors_[xSecVariables_[i]+"ANinclusive"]->SetMarkerColor(kRed);
	      totalErrors_[xSecVariables_[i]+"ANinclusive"]->SetLineColor  (kRed);
	      totalErrors_[xSecVariables_[i]+"ANinclusive"]->SetPoint(0, xSecValue, 1);
	      totalErrors_[xSecVariables_[i]+"ANinclusive"]->SetPointError(0, totalErrorDown, totalErrorUp, 0, 0);
	      totalErrors_[xSecVariables_[i]+"ANinclusive"]->Draw("p same");
	      drawLine(xSecValue-xSecError, 0.99, xSecValue+xSecError, 0.99, kRed, 7, 1);
	      // 2010 CMS result
	      double cmsxSecValue=158.;
	      double cmsErrorUp  =sqrt(10.*10.+6.*6.+15.*15.);
	      double cmsErrorDown=sqrt(10.*10.+6.*6.+15.*15.);
	      double cmsStatError=10;
	      totalErrors_[xSecVariables_[i]+"CMSinclusive"] = new TGraphAsymmErrors(1);
	      totalErrors_[xSecVariables_[i]+"CMSinclusive"]->SetLineWidth(3);
	      totalErrors_[xSecVariables_[i]+"CMSinclusive"]->SetMarkerSize(2.2);
	      totalErrors_[xSecVariables_[i]+"CMSinclusive"]->SetMarkerStyle(22);
	      totalErrors_[xSecVariables_[i]+"CMSinclusive"]->SetMarkerColor(kRed);
	      totalErrors_[xSecVariables_[i]+"CMSinclusive"]->SetLineColor  (kRed);
	      totalErrors_[xSecVariables_[i]+"CMSinclusive"]->SetPoint(0, cmsxSecValue, 3);
	      totalErrors_[xSecVariables_[i]+"CMSinclusive"]->SetPointError(0, cmsErrorDown, cmsErrorUp, 0, 0);
	      totalErrors_[xSecVariables_[i]+"CMSinclusive"]->Draw("p same");
	      drawLine(cmsxSecValue-cmsStatError, 2.99, cmsxSecValue+cmsStatError, 2.99, kRed, 7, 1);
	      finalInclusiveXSec->Draw("axis same");
	      // label 
	      DrawLabel("CMS 2010 combined"     , 0.07, 0.70, 0.4, 0.75); 
	      DrawLabel("2010 data, 36 pb^{-1}" , 0.07, 0.65, 0.4, 0.70);
	      DrawLabel("(TOP-11-001)"          , 0.07, 0.60, 0.4, 0.65);
	      TString channelLabel="unknown";
	      TString dataLabel=Form(dataSample+" data, %2.1f fb^{-1}",luminosity/1000);
	      if(decayChannel.Contains("mu")) channelLabel="#mu + Jets";
	      if(decayChannel.Contains("el")) channelLabel="e + Jets";
	      DrawLabel(channelLabel, 0.07, 0.30, 0.4, 0.35);
	      DrawLabel(dataLabel,    0.07, 0.25, 0.4, 0.30);
	      canvas=(TCanvas*)canvas2->Clone();
	    }
	    // for differential xSecs
	    else{
	      canvas->cd();
	      canvas->SetTitle(xSecVariables_[i]);
	      // unset exp xaxis
	      totalErrors_[xSecVariables_[i]]->GetXaxis()->SetNoExponent(true);
	      double max=totalErrors_[xSecVariables_[i]]->GetMaximum();
	      max*=1.3;
	      if(max>1&&max<100) totalErrors_[xSecVariables_[i]]->GetYaxis()->SetNoExponent(false);
	      else totalErrors_[xSecVariables_[i]]->GetYaxis()->SetNoExponent(true);
	      // for combined cross sections:
	      // change TH1F with statistical errors into 
	      // TGraphAsymmErrors with bin center corrections
	      if(decayChannel=="combined"&&xSecVariables_[i]!="inclusive"){
		TString plotName = xSecVariables_[i];
		std::cout << "shifted stat value for: variable " << plotName << std::endl;
		plotName.ReplaceAll("Norm","");
		// get data plot
		TH1F* dataStat = (TH1F*)canvas->GetPrimitive(plotName+"kData");
		// convert to TGraphAsymmErrors
		TGraphAsymmErrors* statErrors= new TGraphAsymmErrors(dataStat->GetNbinsX());
		statErrors->SetLineWidth(1.0) ; // totalErrors_[xSecVariables_[i]]->GetLineWidth());
		statErrors->SetLineColor(totalErrors_[xSecVariables_[i]]->GetLineColor());
		for(int bin=1; bin<=dataStat->GetNbinsX(); ++bin){
		  if(dataStat->GetBinWidth(bin)!=0){
		    double pointXValue = correctedCenters_[plotName].at(bin-1);
		    double pointXError = corrCenterErrors_[plotName].at(bin-1);
		    double pointYError = dataStat->GetBinError(bin);
		    statErrors->SetPoint(bin, pointXValue, dataStat->GetBinContent(bin));
		    statErrors->SetPoint(0, 0, -1000);
		    pointXError=0;
		    statErrors->SetPointError(bin, pointXError, pointXError, pointYError, pointYError);
		  }
		}
		whipEmptyBinsAway(statErrors, verbose);
		// print new object
		statErrors->Draw("p same");
		// delete (old) TH1F plot from canvas
		canvas->GetListOfPrimitives()->Remove(canvas->GetPrimitive(plotName+"kData"));
	      }
	      // Draw errors into Canvas
	      totalErrors_[xSecVariables_[i]]->Draw("p z same");
	      canvas->SetName (xSecVariables_[i]);
	    }
	    // save Canvas
	    int initialIgnoreLevel=gErrorIgnoreLevel;
	    if(verbose==0) gErrorIgnoreLevel=kWarning;
	    // a) within rootFile
	    if(save) saveToRootFile(outputFile, canvas, true, verbose, "finalXSec");
	    // b) as eps and png
	    if(save){
	      TString saveName=outputFolder+"/xSec/finalXSec"+xSecVariables_[i];
	      if(decayChannel=="combined") saveName+="Combined";
	      canvas->Print(saveName+".eps");
	      canvas->Print(saveName+".png");
	    }
	    gErrorIgnoreLevel=initialIgnoreLevel;
	    //delete statErrors;
	  }
	}
      }
    }
  }
  else std::cout << std::endl << "ERROR:target file does not exist or is broken!" << std::endl;
  // close files
  file->Close();
  if(adpatOldUncertainties) oldAnalysisfile->Close();
  // save relative uncertainty shifted canvas to file
  if(adpatOldUncertainties&&save){ 
    // loop variables
    for(unsigned int i=0; i<xSecVariables_.size(); ++i){
      // loop systematic variations
      for(unsigned int sys=sysNo+1; sys<=sysDiBosDown; ++sys){
	if(adoptedUncertaintyPlot_[xSecVariables_[i]][sys]){
	  saveToRootFile(outputFile, adoptedUncertaintyPlot_[xSecVariables_[i]][sys], false, 0,"xSec/"+sysLabel(sys));
	}
      }
    }
  }
}
