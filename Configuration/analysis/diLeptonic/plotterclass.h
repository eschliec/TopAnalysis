#ifndef plotterclass_h
#define plotterclass_h

#include "THStack.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH1F.h"
#include <vector>
#include <map>
#include <iostream>
#include <cstdio>
#include <fstream>
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TExec.h"
#include "TStyle.h"
#include "TMath.h"
#include "TROOT.h"
#include <sstream>
#include "basicFunctions.h"
#include "HHStyle.h"
#include "TGraphAsymmErrors.h"

// DAVID
#include "DilepSVDFunctions.h"
#include "DilepSVDFunctions.C"

#define TOPXSEC (225.197)
#define LUMI (5100)

class Plotter {

 public:
  Plotter();
  Plotter(TString, TString, TString, double, double);
  ~Plotter();
  void   setOptions(TString, TString, TString, TString, int, bool, bool, bool, double, double, double, double, int, std::vector<double>, std::vector<double>);
  void   setDataSet(std::vector<TString>, std::vector<double>, std::vector<TString>, std::vector<int>, TString);
  void   setDataSet(TString);
  void   setSystDataSet(TString, TString, TString);
  void   fillSystHisto();
  void   fillHisto();
  void setStyle(TH1D&, unsigned int);
  void setStyle(TH1&, unsigned int);
  void   write();
  double CalcXSec(std::vector<TString> datasetVec, double InclusiveXsectionVec[4],double InclusiveXsectionStatErrorVec[4], TString Systematic, TString Shift);
  void MakeTable();
  void PlotXSec();
  void CalcDiffXSec(TH1* varhists[], TH1* RecoPlot, TH1* GenPlot, TH2* genReco2d, double DiffXSecVec[4][10], double DiffXSecStatErrorVec[4][10], double GenDiffXSecVec[4][10]);
  void PlotDiffXSec();
  void DYScaleFactor();
  void CalcInclSystematics(TString Systematic, TString systSampleUp, TString systSampleDown, int syst_number);
  void CalcDiffSystematics(TString Systematic, TString systSampleUp, TString systSampleDown, int syst_number);
  void InclFlatSystematics(int syst_number);
  void DiffFlatSystematics(int syst_number,  int nbins);
  TLegend* getNewLegend();
  TLegend* getNewLegendpre();
  TH1* GetNloCurve(const char *particle, const char *quantity, const char *generator);
  TH1* GetNloCurve(TString NewName, TString Generator);
  TH1F* ConvertGraphToHisto(TGraphErrors *pGraph);
  TH1F* reBinTH1FIrregularNewBinning(TH1F *histoOldBinning, TString plotname, bool rescale);
  // DAVID
  void UnfoldingOptions(bool doSVD);
  void SetOutpath(TString path); 
  void ApplyMCATNLOWeight(TH1* hist, TString Systematic, TString Shift, TString Sample);
  TLegend* ControlLegend(int HistsSize, TH1* drawhists[], std::vector<TString> legends, TLegend *leg);
  TLegend* ControlLegend(int HistsSize, TH1D* drawhists[], std::vector<TString> Legends, TLegend *leg);
  void DrawLabel(TString text, const double x1, const double y1, const double x2, const double y2, int centering, double textSize);
 private:
  TString name;
  TString specialComment;
  int bins, datafiles, rebin;
  double rangemin, rangemax, ymin, ymax;
  std::vector<TFile> files, filesUp, filesDown;
  std::vector<TString> dataset, datasetUp, datasetDown;
  std::vector<double> scales;
  std::vector<TString> legends, legendsUp, legendsDown, legendsSyst;
  std::vector<int> colors, colorsUp, colorsDown;
  std::vector<double> XAxisbins, XAxisbinCenters;
  double DYScale[4];
  TString DYEntry;
  TString YAxis;
  TString XAxis;
  TString channel;
  int channelType; //0=ee 1==mumu 2==emu 3==combined  
  TH1D onehist;
  std::vector<TH1D> hists;
  std::vector<TH1D> systhistsUp;
  std::vector<TH1D> systhistsDown;
  bool initialized, logX, logY, doDYScale;
  int lumi, signalHist;
  TString channelLabel[4];
  double InclusiveXsection[4];
  double InclusiveXsectionStatError[4];
  double InclusiveXsectionTotalError[4];
  double InclusiveXsectionSysError[4];
  double InclusiveXsectionSysErrorBySyst[4][20];
  double InclusiveXsectionBySystNom[4][20];
  double InclusiveXsectionBySystUp[4][20];
  double InclusiveXsectionBySystDown[4][20];
  double InclusiveXsectionStatErrorBySystNom[4][20];
  double InclusiveXsectionStatErrorBySystUp[4][20];
  double InclusiveXsectionStatErrorBySystDown[4][20];
  double VisXsection[4][10];//needed for bin-by-bin corrections to combination [channel][bin]
  double VisXsectionError[4][10];//needed for bin-by-bin corrections to combination [channel][bin]
  double GenDiffXSec[4][10];//the differential X-section per channel by bin [channel][bin]
  double GenDiffXSecError[4][10];//the differential X-section Error per channel by bin [channel][bin]
  double TotalVisXSection[4];
  double TotalVisXSectionNom[4];
  double TotalVisXSectionUp[4];
  double TotalVisXSectionDown[4];
  
  double DiffXSec[4][10];//the differential X-section per channel by bin [channel][bin]
  double DiffXSecStatError[4][10];//the differential X-section Error per channel by bin [channel][bin]
  double DiffXSecSysError[4][10];//the differential X-section Error per channel by bin [channel][bin]
  double DiffXSecSysErrorBySyst[4][10][20];//the differential X-section Error per channel by bin [channel][bin][systematic]
  double DiffXSecTotalError[4][10];//the differential X-section Error per channel by bin [channel][bin]
  double DiffXSecBySystNom[4][10][20];//the differential X-section Error per channel by bin [channel][bin][systematic]
  double DiffXSecBySystUp[4][10][20];//the differential X-section Error per channel by bin [channel][bin][systematic]
  double DiffXSecBySystDown[4][10][20];//the differential X-section Error per channel by bin [channel][bin][systematic]
  double DiffXSecStatErrorBySystNom[4][10][20];//the differential X-section Error per channel by bin [channel][bin][systematic]
  double DiffXSecStatErrorBySystUp[4][10][20];//the differential X-section Error per channel by bin [channel][bin][systematic]
  double DiffXSecStatErrorBySystDown[4][10][20];//the differential X-section Error per channel by bin [channel][bin][systematic]
  double SignalEventswithWeight;

  // DAVID
  bool doUnfolding; 
  bool doSystematics;
  bool drawNLOCurves;
  TString outpath;
  TString outpathPlots;
  TString subfolderChannel;
  TString subfolderSpecial;
};


// DAVID
void Plotter::UnfoldingOptions(bool doSVD)
{
  doUnfolding = doSVD;
  doSystematics = false;
  drawNLOCurves = false;
}


// DAVID
void Plotter::SetOutpath(TString path)
{
  outpath = path;
}



void Plotter::DYScaleFactor(){
  TH1::AddDirectory(kFALSE);
  ifstream FileList("FileLists/HistoFileList_Nominal_combined.txt");
  TString filename;
  
  double NoutEEDYMC=0, NinEEDYMC=0, NoutMuMuDYMC=0, NinMuMuDYMC=0;//Number of events in/out of z-veto region for the DY MC
  double NinEE=0, NinMuMu=0, NinEMu=0;//Number of events in z-veto region for data
  double NinEEloose=0, NinMuMuloose=0;//Number of data events in Z-Veto region with MET cut
  double NinEEMC=0, NinMuMuMC=0;//All other MC events

  while(!FileList.eof()){
    FileList>>filename;
    if(filename!=""){
      TFile *ftemp = TFile::Open(filename);
      if(filename.Contains("ee")){
	if(filename.Contains("run")){
	  TH1D *htemp = (TH1D*)ftemp->Get("Zh1"); NinEE+=htemp->Integral();
	  TH1D *htemp1 = (TH1D*)ftemp->Get("Looseh1"); NinEEloose+=htemp1->Integral();	  
	}
	else if(filename.Contains("dy")){
	  if(filename.Contains("50inf")){
	    TH1D *htemp = (TH1D*)ftemp->Get("Zh1"); NinEEDYMC+=htemp->Integral();
	    TH1D *htemp1 = (TH1D*)ftemp->Get("TTh1"); NoutEEDYMC+=htemp1->Integral();
	  }
	  else{TH1D *htemp = (TH1D*)ftemp->Get("TTh1"); NoutEEDYMC+=htemp->Integral();}
	}	
	else{
	  TH1D *htemp = (TH1D*)ftemp->Get("Zh1"); NinEEMC+=htemp->Integral();
	}
      }
      
      if(filename.Contains("emu")){
	if(filename.Contains("run")){TH1D *htemp = (TH1D*)ftemp->Get("Zh1"); NinEMu+=htemp->Integral();}
      }
	
      if(filename.Contains("mumu")){
	if(filename.Contains("run")){
	  TH1D *htemp = (TH1D*)ftemp->Get("Zh1"); NinMuMu+=htemp->Integral();
	  TH1D *htemp1 = (TH1D*)ftemp->Get("Looseh1"); NinMuMuloose+=htemp1->Integral();
	}
	else if(filename.Contains("dy")){
	  if(filename.Contains("50inf")){
	    TH1D *htemp = (TH1D*)ftemp->Get("Zh1"); NinMuMuDYMC+=htemp->Integral();
	    TH1D *htemp1 = (TH1D*)ftemp->Get("TTh1"); NoutMuMuDYMC+=htemp1->Integral();}
	  else{TH1D *htemp = (TH1D*)ftemp->Get("TTh1"); NoutMuMuDYMC+=htemp->Integral();}
	}	
	else{
	  TH1D *htemp = (TH1D*)ftemp->Get("Zh1"); NinMuMuMC+=htemp->Integral();
	}
      }      
      delete ftemp;
    }
  }
  double NoutMCEE = (NoutEEDYMC/NinEEDYMC)*(NinEE - 0.5*NinEMu*sqrt(NinEEloose/NinMuMuloose));
  double NoutMCMuMu = (NoutMuMuDYMC/NinMuMuDYMC)*(NinMuMu - 0.5*NinEMu*sqrt(NinMuMuloose/NinEEloose));
  
  double DYSFEE = NoutMCEE/NoutEEDYMC;
  double DYSFMuMu = NoutMCMuMu/NoutMuMuDYMC;

  cout << endl;
  cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
  cout << "Calculation of DY Scale Factors for '" << name << "'  " << endl;
   
  cout<<"DYSFEE:                 "<<DYSFEE<<endl;
  cout<<"DYSFMuMu:               "<<DYSFMuMu<<endl;

  cout<<"NinEEloose:             "<<NinEEloose<<endl;
  cout<<"NinMMloose:             "<<NinMuMuloose<<endl;

  cout<<"kee:                    "<<sqrt(NinEEloose/NinMuMuloose)<<" +- "<<0.5*TMath::Sqrt(1./NinMuMuloose + 1./NinEEloose)<<endl;
  cout<<"kmumu:                  "<<sqrt(NinMuMuloose/NinEEloose)<<" +- "<<0.5*TMath::Sqrt(1./NinMuMuloose + 1./NinEEloose)<<endl;

  cout<<"Rout/Rin Mumu:          "<<(NoutMuMuDYMC/NinMuMuDYMC)<<endl;
  cout<<"Rout/Rin ee:            "<<(NoutEEDYMC/NinEEDYMC)<<endl;

  cout<<"Est. From Data(ee):     "<<NoutMCEE<<endl;
  cout<<"Est. From Data(mumu):   "<<NoutMCMuMu<<endl;

  cout<<"Est. From MC(ee):       "<<NoutEEDYMC<<endl;
  cout<<"Est. From MC(mumu):     "<<NoutMuMuDYMC<<endl;
 
  cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
  cout << endl;

  if(doDYScale==true){
    DYScale[0]=DYSFEE;
    DYScale[1]=DYSFMuMu;
    DYScale[2]=1.;
    DYScale[3]=(DYSFEE+DYSFMuMu)/2;//not correct, but close, fix later
  }else{
    DYScale[0]=1.;//need to make a switch for control plots that don't want DYScale
    DYScale[1]=1.;
    DYScale[2]=1.;
    DYScale[3]=1.;
  }
}

void Plotter::InclFlatSystematics(int syst_number){
  
  //trigger uncertainties
  /*
  if (channelType==0){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .025;}//ee 
  if (channelType==1){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .015;}//mumu  
  if (channelType==2){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .02;}//emu  
  if (channelType==3){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = 
      1/(sqrt((1/(InclusiveXsectionSysErrorBySyst[0][syst_number]*InclusiveXsectionSysErrorBySyst[0][syst_number]) + 
  	       1/(InclusiveXsectionSysErrorBySyst[1][syst_number]*InclusiveXsectionSysErrorBySyst[1][syst_number]) +
  	       1/(InclusiveXsectionSysErrorBySyst[2][syst_number]*InclusiveXsectionSysErrorBySyst[2][syst_number]))));}//combined  
  */

  if (channelType==0){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = 0.015;}//ee 
  if (channelType==1){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = 0.015;}//mumu  
  if (channelType==2){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = 0.015;}//emu  
  if (channelType==3){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = 0.015;}//combined  

  syst_number++;

  //Lepton selection

  InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .04;//all 
  syst_number++;

  //B-tagging (for now)
  /*  if (channelType==0){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .034;}//ee 
  if (channelType==1){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .034;}//mumu  
  if (channelType==2){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .034;}//emu  
  if (channelType==3){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .034;}//combined  
  syst_number++;
  */
  //KinFit 
  if (channelType==0){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .02;}//ee 
  if (channelType==1){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .02;}//mumu  
  if (channelType==2){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .017;}//emu  
  if (channelType==3){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .017;}//combined 
  syst_number++;

//  //Hadronization (for now)
//  if (channelType==0){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .049;}//ee 
//  if (channelType==1){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .043;}//mumu  
//  if (channelType==2){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .054;}//emu  
//  if (channelType==3){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = 
//      1/(sqrt((1/(InclusiveXsectionSysErrorBySyst[0][syst_number]*InclusiveXsectionSysErrorBySyst[0][syst_number]) + 
//  	       1/(InclusiveXsectionSysErrorBySyst[1][syst_number]*InclusiveXsectionSysErrorBySyst[1][syst_number]) +
//  	       1/(InclusiveXsectionSysErrorBySyst[2][syst_number]*InclusiveXsectionSysErrorBySyst[2][syst_number]))));}//combined  
//  syst_number++;

  //PDF (for now). Taken rfom l+jet people

  if (channelType==0){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .015;}//ee
  if (channelType==1){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .015;}//mumu
  if (channelType==2){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .015;}//emu
  if (channelType==3){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .015;}//combined
  syst_number++;


  //Branching Ratio
  if (channelType==0){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .015;}//ee 
  if (channelType==1){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .015;}//mumu  
  if (channelType==2){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .015;}//emu  
  if (channelType==3){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .015;}//combined  
  syst_number++;

  //Luminosity
  if (channelType==0){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .022;}//ee 
  if (channelType==1){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .022;}//mumu  
  if (channelType==2){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .022;}//emu  
  if (channelType==3){InclusiveXsectionSysErrorBySyst[channelType][syst_number] = .022;}//combined  
  syst_number++;

}
void Plotter::DiffFlatSystematics(int syst_number, int nbins){
  
  for(int bin = 0; bin<nbins; bin++){
    int syst = syst_number;
    //trigger uncertainties
    legendsSyst.push_back("trigger");
    /*
    if (channelType==0){DiffXSecSysErrorBySyst[channelType][bin][syst] = .025;}//ee 
    if (channelType==1){DiffXSecSysErrorBySyst[channelType][bin][syst] = .015;}//mumu  
    if (channelType==2){DiffXSecSysErrorBySyst[channelType][bin][syst] = .02;}//emu  
    if (channelType==3){DiffXSecSysErrorBySyst[channelType][bin][syst] = 
	1/(sqrt((1/(DiffXSecSysErrorBySyst[0][bin][syst]*DiffXSecSysErrorBySyst[0][bin][syst]) + 
		 1/(DiffXSecSysErrorBySyst[1][bin][syst]*DiffXSecSysErrorBySyst[1][bin][syst]) +
		 1/(DiffXSecSysErrorBySyst[2][bin][syst]*DiffXSecSysErrorBySyst[2][bin][syst]))));}//combined  
    */
    if (channelType==0){DiffXSecSysErrorBySyst[channelType][bin][syst] = .015;}//ee 
    if (channelType==1){DiffXSecSysErrorBySyst[channelType][bin][syst] = .015;}//mumu  
    if (channelType==2){DiffXSecSysErrorBySyst[channelType][bin][syst] = .015;}//emu  
    if (channelType==3){DiffXSecSysErrorBySyst[channelType][bin][syst] = .015;}//combined
    syst++;
    
    //Lepton selection
    legendsSyst.push_back("lepton");
    
    DiffXSecSysErrorBySyst[channelType][bin][syst] = .02;//all (different from inclusive?) 
    syst++;

    //B-tagging (for now)
    /*    legendsSyst.push_back("b-tagging");
    if (channelType==0){DiffXSecSysErrorBySyst[channelType][bin][syst] = .017;}//ee 
    if (channelType==1){DiffXSecSysErrorBySyst[channelType][bin][syst] = .017;}//mumu  
    if (channelType==2){DiffXSecSysErrorBySyst[channelType][bin][syst] = .017;}//emu  
    if (channelType==3){DiffXSecSysErrorBySyst[channelType][bin][syst] = .017;}
	//	1/(sqrt((1/(DiffXSecSysErrorBySyst[0][bin][syst]*DiffXSecSysErrorBySyst[0][bin][syst]) + 
	//	 1/(DiffXSecSysErrorBySyst[1][bin][syst]*DiffXSecSysErrorBySyst[1][bin][syst]) +
	//	 1/(DiffXSecSysErrorBySyst[2][bin][syst]*DiffXSecSysErrorBySyst[2][bin][syst]))));}//combined  
    syst++;
    */
    //KinFit 
    legendsSyst.push_back("kin fit");
    if (channelType==0){DiffXSecSysErrorBySyst[channelType][bin][syst] = .005;}//ee 
    if (channelType==1){DiffXSecSysErrorBySyst[channelType][bin][syst] = .005;}//mumu  
    if (channelType==2){DiffXSecSysErrorBySyst[channelType][bin][syst] = .005;}//emu  
    if (channelType==3){DiffXSecSysErrorBySyst[channelType][bin][syst] = .005;}
    //	1/(sqrt((1/(DiffXSecSysErrorBySyst[0][bin][syst]*DiffXSecSysErrorBySyst[0][bin][syst]) + 
    //		 1/(DiffXSecSysErrorBySyst[1][bin][syst]*DiffXSecSysErrorBySyst[1][bin][syst]) +
    //		 1/(DiffXSecSysErrorBySyst[2][bin][syst]*DiffXSecSysErrorBySyst[2][bin][syst]))));}//combined  
    syst++;
    
//    //Hadronization (for now)
//    legendsSyst.push_back("Hadr.");
//    if (channelType==0){DiffXSecSysErrorBySyst[channelType][bin][syst] = .006;}//ee 
//    if (channelType==1){DiffXSecSysErrorBySyst[channelType][bin][syst] = .01;}//mumu  
//    if (channelType==2){DiffXSecSysErrorBySyst[channelType][bin][syst] = .003;}//emu  
//    if (channelType==3){DiffXSecSysErrorBySyst[channelType][bin][syst] = 
//	1/(sqrt(     1/(DiffXSecSysErrorBySyst[0][bin][syst]*DiffXSecSysErrorBySyst[0][bin][syst]) + 
//		         1/(DiffXSecSysErrorBySyst[1][bin][syst]*DiffXSecSysErrorBySyst[1][bin][syst]) +
//		         1/(DiffXSecSysErrorBySyst[2][bin][syst]*DiffXSecSysErrorBySyst[2][bin][syst])    ));}//combined  
//    syst++;

    //PDF (for now). Taken from l+jet people
    legendsSyst.push_back("PDF");
    if (channelType==0){DiffXSecSysErrorBySyst[channelType][bin][syst] = .015;}//ee 
    if (channelType==1){DiffXSecSysErrorBySyst[channelType][bin][syst] = .015;}//mumu  
    if (channelType==2){DiffXSecSysErrorBySyst[channelType][bin][syst] = .015;}//emu  
    if (channelType==3){DiffXSecSysErrorBySyst[channelType][bin][syst] = .015;}//combined  
    syst++;
    
  }
   
}

void Plotter::CalcInclSystematics(TString Systematic, TString systSampleUp, TString systSampleDown, int syst_number){
  setSystDataSet(Systematic, systSampleUp, systSampleDown);
  fillSystHisto();

  double InclusiveXsectionNom[4], InclusiveXsectionUp[4], InclusiveXsectionDown[4];
  double InclusiveXsectionStatErrorNom[4], InclusiveXsectionStatErrorUp[4], InclusiveXsectionStatErrorDown[4];

  double Sys_Error_Up, Sys_Error_Down, Sys_Error;
  double scale = 1.;

  if(channelType!=3){//set everything to zero
      InclusiveXsectionDown[channelType]=0;
      InclusiveXsectionStatErrorDown[channelType]=0;
      InclusiveXsectionUp[channelType]=0;
      InclusiveXsectionStatErrorUp[channelType]=0;
      InclusiveXsectionNom[channelType]=0;
      InclusiveXsectionStatErrorNom[channelType]=0;
  }else{
    for(int i = 0; i<3; i++){//fill the siffxsec vector with the values for the previous systematics
	  InclusiveXsectionDown[i]=InclusiveXsectionBySystDown[i][syst_number];
	  InclusiveXsectionStatErrorDown[i]=InclusiveXsectionStatErrorBySystDown[i][syst_number];
	  InclusiveXsectionUp[i]=InclusiveXsectionBySystUp[i][syst_number];
	  InclusiveXsectionStatErrorUp[i]=InclusiveXsectionStatErrorBySystUp[i][syst_number];
	  InclusiveXsectionNom[i]=InclusiveXsectionBySystNom[i][syst_number];
	  InclusiveXsectionStatErrorNom[i]=InclusiveXsectionStatErrorBySystNom[i][syst_number];
    }    
  }   
  CalcXSec(dataset, InclusiveXsectionNom,InclusiveXsectionStatErrorNom, Systematic, "None"); 
  CalcXSec(datasetUp, InclusiveXsectionUp,InclusiveXsectionStatErrorUp, Systematic, "Up");  
  CalcXSec(datasetDown, InclusiveXsectionDown,InclusiveXsectionStatErrorDown, Systematic, "Down");  
  
  InclusiveXsectionBySystDown[channelType][syst_number]=InclusiveXsectionDown[channelType];
  InclusiveXsectionStatErrorBySystDown[channelType][syst_number]=InclusiveXsectionStatErrorDown[channelType];
  InclusiveXsectionBySystUp[channelType][syst_number]=InclusiveXsectionUp[channelType];
  InclusiveXsectionStatErrorBySystUp[channelType][syst_number]=InclusiveXsectionStatErrorUp[channelType];
  InclusiveXsectionBySystNom[channelType][syst_number]=InclusiveXsectionNom[channelType];
  InclusiveXsectionStatErrorBySystNom[channelType][syst_number]=InclusiveXsectionStatErrorNom[channelType];

  Sys_Error_Up   = abs(InclusiveXsectionNom[channelType] - InclusiveXsectionUp[channelType])/InclusiveXsectionNom[channelType];
  Sys_Error_Down = abs(InclusiveXsectionNom[channelType] - InclusiveXsectionDown[channelType])/InclusiveXsectionNom[channelType];
  Sys_Error  = (Sys_Error_Up+Sys_Error_Down)/(2.*scale);
  if(Systematic=="MASS"){
    Sys_Error=Sys_Error/12.;
  }

  //  Sum_Errors += Sys_Error;
  InclusiveXsectionSysErrorBySyst[channelType][syst_number] = Sys_Error; 
  cout<<"Inclusive Systematic for '"<<Systematic<<"':";
  cout<<" " << Sys_Error*100 << " % " << endl;
    
}

void Plotter::CalcDiffSystematics(TString Systematic, TString systSampleUp, TString systSampleDown, int syst_number){

  //    TH2 *genReco2d=0; 
  // TString newname = name;
  // if(name.Contains("Hyp")){//Histogram naming convention has to be smarter
  //   newname.ReplaceAll("Hyp",3,"",0);
  // }
  cout << endl;
  cout << endl;
  cout << "    Preparing to Calculate " << Systematic << "-Uncertainty ... " << endl;

  if(Systematic=="BTAG_PT_" || Systematic=="BTAG_ETA_"){  legendsSyst.push_back("BTAG");}
  else {legendsSyst.push_back(Systematic);}
  double Xbins[XAxisbins.size()];
  //  double binWidth[XAxisbinCenters.size()];
  for(unsigned int i = 0; i<XAxisbins.size();i++){Xbins[i]=XAxisbins[i];} 
  setSystDataSet(Systematic, systSampleUp, systSampleDown);   
  fillSystHisto();  
  TH1D* stacksum = (TH1D*)hists[datafiles].Rebin(bins,"stack",Xbins);   
  TH1D* stacksumUp = (TH1D*)systhistsUp[datafiles].Rebin(bins,"stackup",Xbins);  
  TH1D* stacksumDown = (TH1D*)systhistsDown[datafiles].Rebin(bins,"stackdown",Xbins); 

  TH2 *genReco2d=NULL;
  TH2 *genReco2dUp=NULL;
  TH2 *genReco2dDown=NULL;

  TH1 *RecoPlotFineBins = NULL;
  TH1 *RecoPlotFineBinsUp = NULL;
  TH1 *RecoPlotFineBinsDown = NULL;

  TH1 *GenPlotTheory = NULL;
  TH1 *GenPlotTheoryUp = NULL;
  TH1 *GenPlotTheoryDown = NULL;

  TH1 *GenPlot = NULL;
  TH1 *GenPlotUp = NULL;
  TH1 *GenPlotDown = NULL;

  TH1 *RecoPlot = NULL;
  TH1 *RecoPlotUp = NULL;
  TH1 *RecoPlotDown = NULL;

  bool init = false;
  bool initUp = false;
  bool initDown = false;
 
  TH1 *varhists[hists.size()];
  TH1 *varhistsUp[hists.size()];
  TH1 *varhistsDown[hists.size()];
  
  double Sys_Error_Up, Sys_Error_Down, Sys_Error;//, Sum_Errors;

  TString newname = name;
  if(name.Contains("Hyp")){//Histogram naming convention has to be smarter
    newname.ReplaceAll("Hyp",3,"",0);
  } 
   
  if ( doUnfolding == false ) {
    // Prepare Histograms//DY!!
    //Nominal
    for (unsigned int i =0; i<hists.size(); i++){
      varhists[i]=hists[i].Rebin(bins,"varhists",Xbins);            
      setStyle(*varhists[i], i);
      if((legends[i] == DYEntry) && channelType!=2){
	varhists[i]->Scale(DYScale[channelType]);
      }

      if(legends[i] == "t#bar{t} Signal"){
	TFile *ftemp = TFile::Open(dataset[i]);
	if(init==false){
	  RecoPlotFineBins =  (TH1D*)ftemp->Get("Reco"+newname)->Clone();
	  genReco2d = (TH2*)ftemp->Get("GenReco"+newname)->Clone();
	  GenPlotTheory=(TH1D*)ftemp->Get("VisGen"+newname)->Clone();
	  if(newname.Contains("Lepton")||newname.Contains("Top")||newname.Contains("BJet")){
	    RecoPlotFineBins->Add((TH1D*)ftemp->Get("RecoAnti"+newname)->Clone());
	    genReco2d->Add((TH2*)ftemp->Get("GenRecoAnti"+newname)->Clone());
	    GenPlotTheory->Add((TH1D*)ftemp->Get("VisGenAnti"+newname)->Clone());
	  }
	  init =true;
	} else {//account for more than one signal histogram
	  RecoPlotFineBins->Add((TH1D*)ftemp->Get("Reco"+newname)->Clone());
	  genReco2d->Add((TH2*)ftemp->Get("GenReco"+newname)->Clone());
	  GenPlotTheory->Add((TH1D*)ftemp->Get("VisGen"+newname)->Clone());
	  if(newname.Contains("Lepton")||newname.Contains("Top")||newname.Contains("BJet")){
	    GenPlotTheory->Add((TH1D*)ftemp->Get("VisGenAnti"+newname)->Clone());
	    genReco2d->Add((TH2*)ftemp->Get("GenRecoAnti"+newname)->Clone());
	    RecoPlotFineBins->Add((TH1D*)ftemp->Get("RecoAnti"+newname)->Clone());
	  }	  
	}
	delete ftemp;
      }           
    }
    GenPlot = GenPlotTheory->Rebin(bins,"genplot",Xbins);	
    RecoPlot = (TH1D*)RecoPlotFineBins->Rebin(bins,"recohists",Xbins);
  
    double DiffXSecNom[4][10];
    double DiffXSecStatErrorNom[4][10];
    double GenDiffXSecNom[4][10];
    if(channelType!=3){//set everything to zero
      for(int j = 0; j<10; j++){
	DiffXSecNom[channelType][j]=0;
	DiffXSecStatErrorNom[channelType][j]=0;
      }  
    }else{
      for(int i = 0; i<3; i++){//fill the siffxsec vector with the values for the previous systematics
	for(int j = 0; j<10; j++){
	  DiffXSecNom[i][j]=DiffXSecBySystNom[i][j][syst_number];
	  DiffXSecStatErrorNom[i][j]=DiffXSecStatErrorBySystNom[i][j][syst_number];
	}  
      }    
    }
    CalcDiffXSec(varhists, RecoPlot, GenPlot, genReco2d, DiffXSecNom,DiffXSecStatErrorNom,GenDiffXSecNom);//calculate the differential cross-section
    for(int j = 0; j<10; j++){
      DiffXSecBySystNom[channelType][j][syst_number]=DiffXSecNom[channelType][j];
      DiffXSecStatErrorBySystNom[channelType][j][syst_number]=DiffXSecStatErrorNom[channelType][j];
    }  
  
    //Shift Up
    for (unsigned int i =0; i<systhistsUp.size(); i++){
      varhistsUp[i]=systhistsUp[i].Rebin(bins,"varhists",Xbins);            
      setStyle(*varhistsUp[i], i);
      if((legends[i] == DYEntry) && channelType!=2){
	varhistsUp[i]->Scale(DYScale[channelType]);
	//if(Systematic == "DY_"){
	//varhistsUp[i]->Scale(1.3);
	//}
      }else if(legends[i] == "t#bar{t} Signal"){
	TFile *ftempUp = TFile::Open(datasetUp[i]);
	if(initUp==false){
	  RecoPlotFineBinsUp =  (TH1D*)ftempUp->Get("Reco"+newname)->Clone();
	  genReco2dUp = (TH2*)ftempUp->Get("GenReco"+newname)->Clone();
	  GenPlotTheoryUp=(TH1D*)ftempUp->Get("VisGen"+newname)->Clone();
	  if(newname.Contains("Lepton")||newname.Contains("Top")||newname.Contains("BJet")){
	    RecoPlotFineBinsUp->Add((TH1D*)ftempUp->Get("RecoAnti"+newname)->Clone());
	    genReco2dUp->Add((TH2*)ftempUp->Get("GenRecoAnti"+newname)->Clone());
	    GenPlotTheoryUp->Add((TH1D*)ftempUp->Get("VisGenAnti"+newname)->Clone());
	  }
	  initUp =true;
	} else {//account for more than one signal histogram
	  RecoPlotFineBinsUp->Add((TH1D*)ftempUp->Get("Reco"+newname)->Clone());
	  genReco2dUp->Add((TH2*)ftempUp->Get("GenReco"+newname)->Clone());
	  GenPlotTheoryUp->Add((TH1D*)ftempUp->Get("VisGen"+newname)->Clone());
	  if(newname.Contains("Lepton")||newname.Contains("Top")||newname.Contains("BJet")){
	    GenPlotTheoryUp->Add((TH1D*)ftempUp->Get("VisGenAnti"+newname)->Clone());
	    genReco2dUp->Add((TH2*)ftempUp->Get("GenRecoAnti"+newname)->Clone());
	    RecoPlotFineBinsUp->Add((TH1D*)ftempUp->Get("RecoAnti"+newname)->Clone());
	  }	  
	}
	delete ftempUp;
      }else{
	if(Systematic == "BG_" || Systematic == "DY_"){
	  varhistsUp[i]->Scale(1.3);
	}
      }
    }
    GenPlotUp = GenPlotTheoryUp->Rebin(bins,"genplot",Xbins);	
    RecoPlotUp = (TH1D*)RecoPlotFineBinsUp->Rebin(bins,"recohists",Xbins);

    double DiffXSecUp[4][10];
    double DiffXSecStatErrorUp[4][10];
    double GenDiffXSecUp[4][10];
    if(channelType!=3){//set everything to zero
      for(int j = 0; j<10; j++){
	DiffXSecUp[channelType][j]=0;
	DiffXSecStatErrorUp[channelType][j]=0;
      }  
    }else{
      for(int i = 0; i<3; i++){//fill the siffxsec vector with the values for the previous systematics
	for(int j = 0; j<10; j++){
	  DiffXSecUp[i][j]=DiffXSecBySystUp[i][j][syst_number];
	  DiffXSecStatErrorUp[i][j]=DiffXSecStatErrorBySystUp[i][j][syst_number];
	}  
      }    
    }
    CalcDiffXSec(varhistsUp, RecoPlotUp, GenPlotUp, genReco2d, DiffXSecUp,DiffXSecStatErrorUp,GenDiffXSecUp);//calculate the differential cross-section
    for(int j = 0; j<10; j++){
      DiffXSecBySystUp[channelType][j][syst_number]=DiffXSecUp[channelType][j];
      DiffXSecStatErrorBySystUp[channelType][j][syst_number]=DiffXSecStatErrorUp[channelType][j];
    }  

    //Shift Down
    for (unsigned int i =0; i<systhistsDown.size(); i++){
      varhistsDown[i]=systhistsDown[i].Rebin(bins,"varhists",Xbins);            
      setStyle(*varhistsDown[i], i);
      if((legends[i] == DYEntry) && channelType!=2){
	varhistsDown[i]->Scale(DYScale[channelType]);
	//if(Systematic == "DY_"){
	//varhistsDown[i]->Scale(0.7);
	//}
      }else if(legends[i] == "t#bar{t} Signal"){
	TFile *ftempDown = TFile::Open(datasetDown[i]);
	if(initDown==false){
	  RecoPlotFineBinsDown =  (TH1D*)ftempDown->Get("Reco"+newname)->Clone();
	  genReco2dDown = (TH2*)ftempDown->Get("GenReco"+newname)->Clone();
	  GenPlotTheoryDown=(TH1D*)ftempDown->Get("VisGen"+newname)->Clone();
	  if(newname.Contains("Lepton")||newname.Contains("Top")||newname.Contains("BJet")){
	    RecoPlotFineBinsDown->Add((TH1D*)ftempDown->Get("RecoAnti"+newname)->Clone());
	    genReco2dDown->Add((TH2*)ftempDown->Get("GenRecoAnti"+newname)->Clone());
	    GenPlotTheoryDown->Add((TH1D*)ftempDown->Get("VisGenAnti"+newname)->Clone());
	  }
	  initDown =true;
	} else {//account for more than one signal histogram
	  RecoPlotFineBinsDown->Add((TH1D*)ftempDown->Get("Reco"+newname)->Clone());
	  genReco2dDown->Add((TH2*)ftempDown->Get("GenReco"+newname)->Clone());
	  GenPlotTheoryDown->Add((TH1D*)ftempDown->Get("VisGen"+newname)->Clone());
	  if(newname.Contains("Lepton")||newname.Contains("Top")||newname.Contains("BJet")){
	    GenPlotTheoryDown->Add((TH1D*)ftempDown->Get("VisGenAnti"+newname)->Clone());
	    genReco2dDown->Add((TH2*)ftempDown->Get("GenRecoAnti"+newname)->Clone());
	    RecoPlotFineBinsDown->Add((TH1D*)ftempDown->Get("RecoAnti"+newname)->Clone());
	  }	  
	}
	delete ftempDown;
      }else{
	if(Systematic == "BG_" || Systematic == "DY_"){
	  varhistsDown[i]->Scale(0.7);
	}
      }
    }
    GenPlotDown = GenPlotTheoryDown->Rebin(bins,"genplot",Xbins);	
    RecoPlotDown = (TH1D*)RecoPlotFineBinsDown->Rebin(bins,"recohists",Xbins);

    double DiffXSecDown[4][10];
    double DiffXSecStatErrorDown[4][10];
    double GenDiffXSecDown[4][10];
    if(channelType!=3){//set everything to zero
      for(int j = 0; j<10; j++){
	DiffXSecDown[channelType][j]=0;
	DiffXSecStatErrorDown[channelType][j]=0;
      }  
    }else{
      for(int i = 0; i<3; i++){//fill the siffxsec vector with the values for the previous systematics
	for(int j = 0; j<10; j++){
	  DiffXSecDown[i][j]=DiffXSecBySystDown[i][j][syst_number];
	  DiffXSecStatErrorDown[i][j]=DiffXSecStatErrorBySystDown[i][j][syst_number];
	}  
      }    
    }
    CalcDiffXSec(varhistsDown, RecoPlotDown, GenPlotDown, genReco2d, DiffXSecDown,DiffXSecStatErrorDown,GenDiffXSecDown);//calculate the differential cross-section
    for(int j = 0; j<10; j++){
      DiffXSecBySystDown[channelType][j][syst_number]=DiffXSecDown[channelType][j];
      DiffXSecStatErrorBySystDown[channelType][j][syst_number]=DiffXSecStatErrorDown[channelType][j];
    } 

    //DYScale Factor...
    /*  for(unsigned int i=datafiles+1; i<hists.size() ; i++){ // prepare histos and leg
	TFile *ftemp = TFile::Open(dataset[i]);
	if(genReco2d==NULL){
	genReco2d = (TH2*)ftemp->Get("GenReco"+newname)->Clone();
	GenPlot = (TH1*)ftemp->Get("VisGen"+newname)->Clone();
	} else{
	genReco2d->Add((TH2*)ftemp->Get("GenReco"+newname)->Clone());
	GenPlot->Add((TH1*)ftemp->Get("VisGen"+newname)->Clone());
	}
	TH1 *htemp = (TH1D*)hists[i].Rebin(bins,"htemp",Xbins);
	if((legends[i] == DYEntry) && channelType!=2 ){
	htemp->Scale(DYScale[channelType]);
	}
	stacksum->Add(htemp);
	delete ftemp;
	} 
	for(unsigned int i=datafiles+1; i<systhistsUp.size() ; i++){ // prepare histos and leg
	TFile *ftempUp = TFile::Open(datasetUp[i]);
	if(genReco2dUp==NULL){
	genReco2dUp = (TH2*)ftempUp->Get("GenReco"+newname)->Clone();
	} else{
	genReco2dUp->Add((TH2*)ftempUp->Get("GenReco"+newname)->Clone());
	}
	TH1 *htemp = (TH1D*)systhistsUp[i].Rebin(bins,"htempup",Xbins);
	if((legends[i] == DYEntry) && channelType!=2 ){
	htemp->Scale(DYScale[channelType]);
	}
	stacksumUp->Add(htemp);
	if((legends[i] == DYEntry) && Systematic == "DY_"){
	htemp->Scale(1.5);
	}
	if((legends[i] != DYEntry)&& (legends[i] != "t#bar{t} Signal") && Systematic == "BG_"){
	htemp->Scale(1.3);
	}

	delete ftempUp;
	} 
	for(unsigned int i=datafiles+1; i<systhistsDown.size() ; i++){ // prepare histos and leg
	TFile *ftempDown = TFile::Open(datasetDown[i]);
	if(genReco2dDown==NULL){
	genReco2dDown = (TH2*)ftempDown->Get("GenReco"+newname)->Clone();
	} else{
	genReco2dDown->Add((TH2*)ftempDown->Get("GenReco"+newname)->Clone());
	}
	TH1 *htemp = (TH1D*)systhistsDown[i].Rebin(bins,"htempdown",Xbins);
	if((legends[i] == DYEntry) && channelType!=2 ){
	htemp->Scale(DYScale[channelType]);
	}
	if((legends[i] == DYEntry) && Systematic == "DY_"){
	htemp->Scale(0.5);
	}
	if((legends[i] != DYEntry)&& (legends[i] != "t#bar{t} Signal") && Systematic == "BG_"){
	htemp->Scale(0.7);
	}
	stacksumDown->Add(htemp);
	delete ftempDown;
	}*/
   
    double scale = 1.;
    TH1 *h_DiffXSec = (TH1D*)varhists[0]->Clone();
    TH1 *h_DiffXSecUp = (TH1D*)varhistsUp[0]->Clone();
    TH1 *h_DiffXSecDown = (TH1D*)varhistsDown[0]->Clone();
    h_DiffXSec->Reset();
    h_DiffXSecUp->Reset();
    h_DiffXSecDown->Reset();
    for (Int_t i=0; i<bins; ++i) {
      h_DiffXSec->SetBinContent(i+1,DiffXSecNom[channelType][i]);
      h_DiffXSecUp->SetBinContent(i+1,DiffXSecUp[channelType][i]);
      h_DiffXSecDown->SetBinContent(i+1,DiffXSecDown[channelType][i]);
    }
    double datascale=0;
    double datascaleUp=0;
    double datascaleDown=0;
    datascale = h_DiffXSec->Integral("width");
    datascaleUp = h_DiffXSecUp->Integral("width");
    datascaleDown = h_DiffXSecDown->Integral("width");
    //if(name.Contains("LeptonEta")){
    TotalVisXSectionNom[channelType] = datascale;
    TotalVisXSectionUp[channelType] = datascaleUp;
    TotalVisXSectionDown[channelType] = datascaleDown;
    //}
    //cout<<"VISIBLE CROSS-SECTION: "<<TotalVisXSection[channelType]<<endl;
    for (Int_t i=0; i<bins; ++i) {
      DiffXSecNom[channelType][i] = DiffXSecNom[channelType][i]/TotalVisXSectionNom[channelType];
      DiffXSecUp[channelType][i] = DiffXSecUp[channelType][i]/TotalVisXSectionUp[channelType];
      DiffXSecDown[channelType][i] = DiffXSecDown[channelType][i]/TotalVisXSectionDown[channelType];
    }

    /*  for(Int_t bin = 0; bin <= stacksum->GetNbinsX(); ++bin) {
	if( stacksum->GetBinContent(bin+1) == 0 )  continue;
	Sys_Error_Up   = abs(stacksum->GetBinContent(bin+1) - stacksumUp->GetBinContent(bin+1))/stacksum->GetBinContent(bin+1);
	Sys_Error_Down = abs(stacksum->GetBinContent(bin+1) - stacksumDown->GetBinContent(bin+1))/stacksum->GetBinContent(bin+1);
	Sys_Error  = (Sys_Error_Up+Sys_Error_Down)/(2.*scale);
	//    Sum_Errors += Sys_Error;
	DiffXSecSysErrorBySyst[channelType][bin][syst_number] = Sys_Error;//the differential X-section Error per channel by bin [channel][bin][systematic]
	} */
   
    for (Int_t i=0; i<bins; ++i) {
      if( DiffXSecNom[channelType][i]== 0 )  continue;
      Sys_Error_Up   = abs(DiffXSecNom[channelType][i] - DiffXSecUp[channelType][i])/DiffXSecNom[channelType][i];
      Sys_Error_Down   = abs(DiffXSecNom[channelType][i] - DiffXSecDown[channelType][i])/DiffXSecNom[channelType][i];
      Sys_Error  = (Sys_Error_Up+Sys_Error_Down)/(2.*scale);
      //    Sum_Errors += Sys_Error;
      DiffXSecSysErrorBySyst[channelType][i][syst_number] = Sys_Error;//the differential X-section Error per channel by bin [channel][bin][systematic]
      cout<<"Differential Systematic: "<<Systematic<<" NAME: "<<name<<" channel: "<<channelType<<endl;
      cout<<"DiffXSecNom[channelType][i]: "<<DiffXSecNom[channelType][i]<<endl;
      cout<<"DiffXSecUp[channelType][i]: "<<DiffXSecUp[channelType][i]<<endl;
      cout<<"DiffXSecDown[channelType][i]: "<<DiffXSecDown[channelType][i]<<endl;
      cout<<"SYS ERROR: "<<Sys_Error<<endl;
      cout<<endl;
    } 
  }  
  // DAVID. Guckst du hier! 
  if ( doUnfolding == true ) {
  
    // SVD Helper Class
    DilepSVDFunctions mySVDFunctions; 
    mySVDFunctions.SetOutputPath(outpath);  
		
    // Variables for the needed histograms
    TH1D* theDataHist = NULL;
    TH1D* theBgrHist = NULL; 
    TH1D* theBgrHistUp = NULL;
    TH1D* theBgrHistDown = NULL;
    TH1D* theTtBgrHist = NULL; 
    TH1D* theTtBgrHistUp = NULL;
    TH1D* theTtBgrHistDown = NULL;
    TH1D* theRecHist = NULL;
    TH1D* theRecHistUp = NULL;
    TH1D* theRecHistDown = NULL;
    TH1D* theGenHist = NULL;
    TH1D* theGenHistUp = NULL;
    TH1D* theGenHistDown = NULL; 
    TH2D* theRespHist = NULL;
    TH2D* theRespHistUp = NULL;
    TH2D* theRespHistDown = NULL;   
		 
  
    // DAVID:
    // Data, Signal and Background
    // can be obtained from vectors that Tyler fills.
    // These are the vectors
    // "hists", "systhistsUp" amd "systhistsDown"
    // Notice, that the DY Background will be scaled with
    // the DYScale.
    for ( size_t i = 0; i < hists.size() ; i++ ) {
      TString histsIntegral = TString::Format("%20.2f", ((hists[i])).Integral());
      cout << "        " << histsIntegral << " Events in Sample " <<    legends[i] << " (Nom)" << endl;
      if ( legends[i] == "Data" ) { 
		if ( theDataHist == NULL ) {
		  theDataHist = (TH1D*) (hists[i]).Clone("theDataHist");
		} else {
		  theDataHist->Add(&(hists[i]));
		}
      } else if ( legends[i] == "t#bar{t} Signal") {  
		if ( theRecHist == NULL ) {
		  theRecHist = (TH1D*) (hists[i]).Clone("theRecHist");
		} else {
		  theRecHist->Add(&(hists[i]));
		}   
      } else if ( legends[i] == "t#bar{t} Other") {  
		if ( theTtBgrHist == NULL ) {
		  theTtBgrHist = (TH1D*) (hists[i]).Clone("theTtBgrHist");
		} else {
		  theTtBgrHist->Add(&(hists[i]));
		}
		if ( theBgrHist == NULL ) {
		  theBgrHist = (TH1D*) (hists[i]).Clone("theBgrHist");
		} else {
		  theBgrHist->Add(&(hists[i]));
		}  
      } else {  
		// Calculate DYScale fist
		double theDYScale = 1.;
		if ( (legends[i] == DYEntry) && (channelType != 2 ) ) {
		  theDYScale = DYScale[channelType];
		}
					
		// Get Background HIst, thereby applying DY Scale
		if ( theBgrHist == NULL ) {
		  theBgrHist = (TH1D*) (hists[i]).Clone("theBgrHist");
		  theBgrHist->Scale(theDYScale);
		} else {
		  theBgrHist->Add(&(hists[i]), theDYScale);
		}  
      }
    }
    
    for ( size_t i = 0; i < systhistsUp.size() ; i++ ) {  
      TString histsIntegral = TString::Format("%20.2f", ((systhistsUp[i])).Integral());
      cout << "        " << histsIntegral << " Events in Sample " <<    legendsUp[i] << " (Up)" << endl;  
      if ( legendsUp[i] == "Data") {
      } 
      else if ( legendsUp[i] == "t#bar{t} Signal") {   
		if ( theRecHistUp == NULL ) {  
		  theRecHistUp = (TH1D*) (systhistsUp[i]).Clone("theRecHistUp");  
		} else {
		  theRecHistUp->Add(&(systhistsUp[i]));
		}  
      } 
      else if ( legendsUp[i] == "t#bar{t} Other") {  
		if ( theTtBgrHistUp == NULL ) { 
		  theTtBgrHistUp = (TH1D*) (systhistsUp[i]).Clone("theTtBgrHistUp");  
		} else {
		  theTtBgrHistUp->Add(&(systhistsUp[i]));
		}  
		if ( theBgrHistUp == NULL ) {
		  theBgrHistUp = (TH1D*) (systhistsUp[i]).Clone("theBgrHistUp");
		} else {
		  theBgrHistUp->Add(&(systhistsUp[i]));
		} 
      } else {  

    	// Calculate DYScale fist
	    double theDYScale = 1.;
    	if ( (legendsUp[i] == DYEntry) && (channelType != 2 ) ) {
	      theDYScale = DYScale[channelType];
    	}

        // Background Scale Factor
        double theBgrScaleFactor = 1.;
        if ( legendsUp[i] == DYEntry ) {
            theBgrScaleFactor = theDYScale;
        }
        if(Systematic == "BG_") {
          if ( legendsUp[i] != "t#bar{t} Other" && legendsUp[i] != DYEntry) {
             theBgrScaleFactor=theBgrScaleFactor*1.3;
          }
        }
        if(Systematic == "DY_") {
          if ( legendsUp[i] == DYEntry) {
             theBgrScaleFactor=theBgrScaleFactor*1.3;
          }
        }
        
        // Fetching and Scaling the Background Histos 
        if ( theBgrHistUp == NULL ) {
          theBgrHistUp = (TH1D*) (systhistsUp[i]).Clone("theBgrHistUp");
          theBgrHistUp->Scale(theBgrScaleFactor);
        } else {
          theBgrHistUp->Add(&(systhistsUp[i]), theBgrScaleFactor);
        }

      }
    }
    
    for ( size_t i = 0; i < systhistsDown.size() ; i++ ) { 
      TString histsIntegral = TString::Format("%20.2f", ((systhistsDown[i])).Integral());
      cout << "        " << histsIntegral << " Events in Sample " <<    legendsDown[i] << " (Down)" << endl;   
      if ( legendsDown[i] == "Data") {
      } else if ( legendsDown[i] == "t#bar{t} Signal") {   
		if ( theRecHistDown == NULL ) { 
		  theRecHistDown = (TH1D*) (systhistsDown[i]).Clone("theRecHistDown"); 
		} else {
		  theRecHistDown->Add(&(systhistsDown[i]));
		}  
      } else if ( legendsDown[i] == "t#bar{t} Other") {   
		if ( theTtBgrHistDown == NULL ) { 
		  theTtBgrHistDown = (TH1D*) (systhistsDown[i]).Clone("theTtBgrHistDown"); 
		} else {
		  theTtBgrHistDown->Add(&(systhistsDown[i]));
		}  
		if ( theBgrHistDown == NULL ) {
		  theBgrHistDown = (TH1D*) (systhistsDown[i]).Clone("theBgrHistDown");
		} else {
		  theBgrHistDown->Add(&(systhistsDown[i]));
		}
      } else {  
        
        // Calculate DYScale fist
        double theDYScale = 1.;
        if ( (legendsDown[i] == DYEntry) && (channelType != 2 ) ) {
          theDYScale = DYScale[channelType];
        }

        // Background Scale Factor
        double theBgrScaleFactor = 1.;
        if ( legendsDown[i] == DYEntry ) {
            theBgrScaleFactor = theDYScale;
        }
        if(Systematic == "BG_") {
          if ( legendsDown[i] != "t#bar{t} Other" && legendsDown[i] != DYEntry) {
             theBgrScaleFactor=theBgrScaleFactor*0.7;
          }
        }
        if(Systematic == "DY_") {
          if ( legendsDown[i] == DYEntry) {
             theBgrScaleFactor=theBgrScaleFactor*0.7;
          }
        }

        // Fetching and Scaling the Background Histos
        if ( theBgrHistDown == NULL ) {
          theBgrHistDown = (TH1D*) (systhistsDown[i]).Clone("theBgrHistDown");
          theBgrHistDown->Scale(theBgrScaleFactor);
        } else {
          theBgrHistDown->Add(&(systhistsDown[i]), theBgrScaleFactor);
        }
 
      } 
    }
        
         
		
    // Do we have to take care of antiparticles?
    bool doAntiParticles = false;
    if ( newname.Contains("Lepton")) doAntiParticles = true;
    if ( newname.Contains("BJet")) doAntiParticles = true;
    if ( newname.Contains("Top")) doAntiParticles = true;
		 
		 
	 
    // DAVID
    // Response Matrix and visible Gen Dist
    // have to be obtained from the files directly
    for ( size_t i = 0; i < hists.size() ; i++ ) {
      if(legends[i] == "t#bar{t} Signal"){
		 
		// Histo Names
		TString respHistName("GenReco"+newname);
		TString genHistName("Gen"+newname);
		 
		    	
		// Nominal Sample
		TFile* ftemp = TFile::Open(dataset[i]);
		TH2D* resptemp = (TH2D*) ftemp->Get("GenReco"+newname);
		TH1D* gentemp = (TH1D*) ftemp->Get("VisGen"+newname); 
		if ( theRespHist == NULL ) {
		  theRespHist = (TH2D*) (resptemp)->Clone("theRespHist");
		} else {
		  theRespHist->Add(resptemp);
		}
		if ( theGenHist == NULL ) {
		  theGenHist = (TH1D*) (gentemp)->Clone("theGenHist");
		} else {
		  theGenHist->Add(gentemp);
		}
					 
		// Antiparticles 
		if ( doAntiParticles == true ) {
		  resptemp =  (TH2D*) ftemp->Get("GenRecoAnti"+newname);
		  gentemp = (TH1D*) ftemp->Get("VisGenAnti"+newname); 
		  theRespHist->Add(resptemp);
		  theGenHist->Add(gentemp);
		} 
					
		delete ftemp;
      }
    }
    
    for ( size_t i = 0; i < systhistsUp.size() ; i++ ) { 
      if(legendsUp[i] == "t#bar{t} Signal"){
			    	
		// Shift Up
		TFile* ftempUp = TFile::Open(datasetUp[i]);
		TH2D* resptempUp = (TH2D*) ftempUp->Get("GenReco"+newname);
		TH1D* gentempUp = (TH1D*) ftempUp->Get("VisGen"+newname); 
		if ( theRespHistUp == NULL ) {
		  theRespHistUp = (TH2D*) (resptempUp)->Clone("theRespHistUp");
		} else {
		  theRespHistUp->Add(resptempUp);
		}
		if ( theGenHistUp == NULL ) {
		  theGenHistUp = (TH1D*) (gentempUp)->Clone("theGenHistUp");
		} else {
		  theGenHistUp->Add(gentempUp);
		} 
					 
		// Antiparticles 
		if ( doAntiParticles == true ) {
		  resptempUp =  (TH2D*) ftempUp->Get("GenRecoAnti"+newname);
		  gentempUp = (TH1D*) ftempUp->Get("VisGenAnti"+newname); 
		  theRespHistUp->Add(resptempUp);
		  theGenHistUp->Add(gentempUp);
		}
					
		delete ftempUp;
		

      }
    }
    
  
    
    for ( size_t i = 0; i < systhistsDown.size() ; i++ ) { 
      if(legendsDown[i] == "t#bar{t} Signal"){
					
		// Shift Down
		TFile* ftempDown = TFile::Open(datasetDown[i]);
		TH2D* resptempDown = (TH2D*) ftempDown->Get("GenReco"+newname);
		TH1D* gentempDown = (TH1D*) ftempDown->Get("VisGen"+newname); 
		if ( theRespHistDown == NULL ) {
		  theRespHistDown = (TH2D*) (resptempDown)->Clone("theRespHistDown");
		} else {
		  theRespHistDown->Add(resptempDown);
		}
		if ( theGenHistDown == NULL ) {
		  theGenHistDown = (TH1D*) (gentempDown)->Clone("theGenHistDown");
		} else {
		  theGenHistDown->Add(gentempDown);
		} 
					 
		// Antiparticles 
		if ( doAntiParticles == true ) {
		  resptempDown =  (TH2D*) ftempDown->Get("GenRecoAnti"+newname);
		  gentempDown = (TH1D*) ftempDown->Get("VisGenAnti"+newname); 
		  theRespHistDown->Add(resptempDown);
		  theGenHistDown->Add(gentempDown);
		}
					
		delete ftempDown;
		 
		 
      }
    } 
     
      
    // Apply Scale Factor for MC@NLO 
	if ( legendsSyst.size() > 0 ) {  
        ApplyMCATNLOWeight(theRespHistUp, legendsSyst.back(), "Up", "ttbarsignal");  
        ApplyMCATNLOWeight(theRespHistDown, legendsSyst.back(), "Down", "ttbarsignal");  
        ApplyMCATNLOWeight(theGenHistUp, legendsSyst.back(), "Up", "ttbarsignal");  
        ApplyMCATNLOWeight(theGenHistDown, legendsSyst.back(), "Down", "ttbarsignal");  
    }
		 
    // Get the binning 
    double* theBins = Xbins;
    int numberBins = bins;
		
		
    // Names and Labels
    TString channelLabelStr(channelLabel[channelType]);
    TString theChannelName = channelLabelStr; 		
    if ( channelLabelStr.Contains("#mu#mu")  ) theChannelName = "mumu";
    if ( channelLabelStr.Contains("e#mu")    ) theChannelName = "emu";
    if ( channelLabelStr.Contains("ee")      ) theChannelName = "ee";
    if ( channelLabelStr.Contains("Dilepton Combined")    ) theChannelName = "combined";
    TString theParticleName = "";
    if ( name.Contains("Lepton") ) theParticleName = "Leptons";
    if ( name.Contains("LLBar")   ) theParticleName = "LepPair";
    if ( name.Contains("Top")     ) theParticleName = "TopQuarks";
    if ( name.Contains("TTBar")   ) theParticleName = "TtBar";
    if ( name.Contains("BJet")    ) theParticleName = "BJets";
    TString theQuantityName = "";
    if ( name.Contains("pT")      ) theQuantityName = "Pt";
    if ( name.Contains("Eta")     ) theQuantityName = "Eta";
    if ( name.Contains("Rapidity")) theQuantityName = "Rapidity";
    if ( name.Contains("Mass")    ) theQuantityName = "Mass";
    TString theSpecialPostfix = "";
    if ( specialComment.CompareTo("Standard") != 0 ) {
    	//theSpecialPostfix = specialComment;
    } 
    TString theSystematicName = Systematic; 
 
    cout << endl;
    cout << endl;
   
    // Get the integrals for the normalization
    double totalDataEventsNom  = TopSVDFunctions::SVD_Integral1D((TH1D*)theDataHist, 0, false);
    double totalBgrEventsNom   = TopSVDFunctions::SVD_Integral1D((TH1D*)theBgrHist, 0, false);
    double totalBgrEventsUp    = TopSVDFunctions::SVD_Integral1D((TH1D*)theBgrHistUp, 0, false);
    double totalBgrEventsDown  = TopSVDFunctions::SVD_Integral1D((TH1D*)theBgrHistDown, 0, false);
    double totalTtBgrEventsNom   = TopSVDFunctions::SVD_Integral1D((TH1D*)theTtBgrHist, 0, false);
    double totalTtBgrEventsUp    = TopSVDFunctions::SVD_Integral1D((TH1D*)theTtBgrHistUp, 0, false);
    double totalTtBgrEventsDown  = TopSVDFunctions::SVD_Integral1D((TH1D*)theTtBgrHistDown, 0, false);
    double totalRecEventsNom   = TopSVDFunctions::SVD_Integral1D((TH1D*)theRecHist, 0, false);
    double totalRecEventsUp    = TopSVDFunctions::SVD_Integral1D((TH1D*)theRecHistUp, 0, false);
    double totalRecEventsDown  = TopSVDFunctions::SVD_Integral1D((TH1D*)theRecHistDown, 0, false);
    double totalGenEventsNom   = TopSVDFunctions::SVD_Integral1D((TH1D*)theGenHist, 0, false);
    double totalGenEventsUp    = TopSVDFunctions::SVD_Integral1D((TH1D*)theGenHistUp, 0, false);
    double totalGenEventsDown  = TopSVDFunctions::SVD_Integral1D((TH1D*)theGenHistDown, 0, false);

 
    // UNFOLDING OF SYSTEMATICS
    // Retrieve histograms with the unfolded quantities.
    // Note: The unfolded histograms have additional side bins!
    // Keep this in mind when accessing bin content via indices 
    TH1D* symmSysErrors = NULL; 
    mySVDFunctions.SVD_DoUnfoldSys(
				   theDataHist,
				   theBgrHist, theBgrHistUp, theBgrHistDown, 
				   theTtBgrHist, theTtBgrHistUp, theTtBgrHistDown, 
				   theGenHist, theGenHistUp, theGenHistDown, 
				   theRecHist, theRecHistUp, theRecHistDown, 
				   theRespHist, theRespHistUp, theRespHistDown, 
                   totalDataEventsNom, 
                   totalBgrEventsNom,  totalBgrEventsUp,  totalBgrEventsDown, 
                   totalTtBgrEventsNom,  totalTtBgrEventsUp,  totalTtBgrEventsDown, 
                   totalRecEventsNom,  totalRecEventsUp,  totalRecEventsDown, 
                   totalGenEventsNom,  totalGenEventsUp,  totalGenEventsDown,  
				   theBins, numberBins,
				   symmSysErrors,  
				   theChannelName, theParticleName, theQuantityName, theSpecialPostfix, theSystematicName
				   ); 
    
    
//     // Save the shifts in Tyler's triple-matrix ...
//     for(Int_t bin = 0; bin <= stacksum->GetNbinsX(); ++bin) {
//         Sys_Error = symmSysErrors->GetBinContent(bin+2); // Keep in mind the extra layer of OF bins
//         if(Systematic == "MASS"){
//         Sys_Error = Sys_Error/12.;
//         }
//         // Save it
//         DiffXSecSysErrorBySyst[channelType][bin][syst_number] = Sys_Error;//the differential X-section Error per channel by bin [channel][bin][systematic]
//     }
    
    
    
    //Symetrize Eta and Rapidity distributions
    if (theQuantityName == "Eta" || theQuantityName == "Rapidity" ){
        for(int j=0; j<(int) symmSysErrors->GetNbinsX(); ++j){
            cout<<"In bin "<<j<<" binCenter "<<symmSysErrors->GetBinCenter(j+1)<<" Content "<<symmSysErrors->GetBinContent(j+1)<<endl;
        }

        int Nbins = stacksum->GetNbinsX();
        cout<<"Nbins in "<<symmSysErrors->GetName()<<" = "<<symmSysErrors->GetNbinsX()<<endl;
        //There are 2 extra bins coming from the unfolding ==>  skip the underflow+1 bin from left and and overflow+1 bin from right
        for(int i=0; i<Nbins; ++i){
            cout<<"(2nd loop) In bin "<<i<<" binCenter "<<symmSysErrors->GetBinCenter(i+2)<<" Content "<<symmSysErrors->GetBinContent(i+2)<<endl;
            cout<<"                     binCenter "<<symmSysErrors->GetBinCenter(Nbins-i+1)<<" Content "<<symmSysErrors->GetBinContent(Nbins-i+1)<<endl;
            Sys_Error = 0.5*(symmSysErrors->GetBinContent(i+2)+symmSysErrors->GetBinContent(Nbins+1-i));
            cout<<"Symetrized error "<<Sys_Error<<endl;
            if(Systematic == "MASS"){
                Sys_Error = Sys_Error/12.;
            }
            // Save it
            DiffXSecSysErrorBySyst[channelType][i][syst_number] = Sys_Error;//the differential X-section Error per channel by bin [channel][bin][systematic]
        }
    }
    else{
        // Save the shifts in Tyler's triple-matrix ...
        for(Int_t bin = 0; bin < stacksum->GetNbinsX(); ++bin) {
            Sys_Error = symmSysErrors->GetBinContent(bin+2); // Keep in mind the extra layer of OF bins
            if(Systematic == "MASS"){
                Sys_Error = Sys_Error/12.;
            }
            // Save it
            DiffXSecSysErrorBySyst[channelType][bin][syst_number] = Sys_Error;//the differential X-section Error per channel by bin [channel][bin][systematic]
        }
    }
  }
     
  // Old pre SVD code
  if ( doUnfolding == false) {
    double scale = 1.0;
    for(Int_t bin = 0; bin <= stacksum->GetNbinsX(); ++bin) {
      if( stacksum->GetBinContent(bin+1) == 0 )  continue;
      Sys_Error_Up   = abs(stacksum->GetBinContent(bin+1) - stacksumUp->GetBinContent(bin+1))/stacksum->GetBinContent(bin+1);
      Sys_Error_Down = abs(stacksum->GetBinContent(bin+1) - stacksumDown->GetBinContent(bin+1))/stacksum->GetBinContent(bin+1);
      Sys_Error  = (Sys_Error_Up+Sys_Error_Down)/(2.*scale);
      //    Sum_Errors += Sys_Error;
      DiffXSecSysErrorBySyst[channelType][bin][syst_number] = Sys_Error;//the differential X-section Error per channel by bin [channel][bin][systematic]
    } 
  }
}

Plotter::Plotter()
{
  name="defaultName";
  specialComment="Standard";
  rangemin=0;
  rangemax=3;
  YAxis="N_{events}";
  initialized=false;
  datafiles = 0;
  
  // DAVID  
  outpath = "";
  outpathPlots = "Plots";
  subfolderChannel = "";
  subfolderSpecial = "";
    doSystematics = false;
}

Plotter::Plotter(TString name_, TString XAxis_,TString YAxis_, double rangemin_, double rangemax_)
{
  name=name_;
  rangemin=rangemin_;
  rangemax=rangemax_;
  XAxis=XAxis_;
  YAxis=YAxis_;
  initialized=false;
}

Plotter::~Plotter()
{
}

void Plotter::setOptions(TString name_, TString specialComment_, TString YAxis_, TString XAxis_, int rebin_, bool doDYScale_, bool logX_, bool logY_, double ymin_, double ymax_, double rangemin_, double rangemax_, int bins_, std::vector<double>XAxisbins_, std::vector<double>XAxisbinCenters_)
{
  XAxisbins.clear();
  XAxisbinCenters.clear();
  XAxisbinCenters = XAxisbinCenters_;
  XAxisbins = XAxisbins_;
  rebin=rebin_;
  name=name_;
  specialComment=specialComment_;
  bins=bins_;
  doDYScale = doDYScale_;
  logX = logX_;
  logY = logY_;
  ymin = ymin_;
  ymax = ymax_;
  rangemin=rangemin_;
  rangemax=rangemax_;
  YAxis=YAxis_;
  XAxis=XAxis_;
  if(XAxis.Contains("band#bar{b}")){//Histogram naming convention has to be smarter
    XAxis.ReplaceAll("band#bar{b}",11,"b and #bar{b}",13);
  }
  if(XAxis.Contains("tand#bar{t}")){//Histogram naming convention has to be smarter
    XAxis.ReplaceAll("tand#bar{t}",11,"t and #bar{t}",13);
  }
  if(XAxis.Contains("l^{+}andl^{-}")){//Histogram naming convention has to be smarter
    XAxis.ReplaceAll("l^{+}andl^{-}",13,"l^{+} and l^{-}",15);
  }
  if(YAxis.Contains("Toppairs")){
    YAxis.ReplaceAll("Toppairs",8,"Top-quark pairs",15);
  }
  if(YAxis.Contains("Topquarks")){
    YAxis.ReplaceAll("Topquarks",9, "Top quarks",10);
  }
  if(YAxis.Contains("Numberof")){
      YAxis.ReplaceAll("Numberof", 8, "Number of ",10);
  }

  DYScale[0]=1.;
  DYScale[1]=1.;
  DYScale[2]=1.;
  DYScale[3]=1.;
  for(int i=0; i<4; i++){
    for(int j=0; j<10; j++){
      DiffXSec[i][j]=0;
      DiffXSecSysError[i][j]=0;
      DiffXSecStatError[i][j]=0;
      DiffXSecTotalError[i][j]=0;
      for(int k=0; k<20; k++){
	DiffXSecSysErrorBySyst[i][j][k]=0;
	DiffXSecStatErrorBySystNom[i][j][k]=0;
	DiffXSecStatErrorBySystUp[i][j][k]=0;
	DiffXSecStatErrorBySystDown[i][j][k]=0;
	DiffXSecBySystNom[i][j][k]=0;
	DiffXSecBySystUp[i][j][k]=0;
	DiffXSecBySystDown[i][j][k]=0;
      }
    }
  }
}


void Plotter::setDataSet(std::vector<TString> dataset_, std::vector<double> scales_, std::vector<TString> legends_, std::vector<int> colors_, TString DYEntry_)
{
  dataset.clear();
  scales.clear();
  legends.clear();
  legendsSyst.clear();
  colors.clear();
  dataset=dataset_;
  scales=scales_;
  legends=legends_;
  colors=colors_;
  DYEntry=DYEntry_;
  
}

void Plotter::setDataSet(TString mode)
{
  initialized=false;
  legendsSyst.clear();

  channel=mode;
  if(channel =="ee"){channelType=0;channelLabel[0]="ee";}
  if(channel =="mumu"){channelType=1;channelLabel[1]="#mu#mu";}
  if(channel =="emu"){channelType=2;channelLabel[2]="e#mu";}
  if(channel =="combined"){channelType=3;channelLabel[3]="Dilepton Combined";}


  // Set dataset specific subfolders
  outpathPlots = "./Plots";
  subfolderChannel = channel; 
  subfolderChannel.Prepend("/");
  subfolderSpecial = "";
  if ( specialComment.CompareTo("Standard") != 0 ) {
  	subfolderSpecial = specialComment;
  	subfolderSpecial.Prepend("/");
  }
  

  lumi=LUMI;
  
  DYEntry = "Z / #gamma* #rightarrow ee/#mu#mu";

  if(channel=="ee" || channel=="emu" ||channel=="mumu"){  
    ifstream FileList("FileLists/HistoFileList_Nominal_"+channel+".txt");
    TString filename;
    datafiles=0;

    dataset.clear();
    legends.clear();
    colors.clear();
  
    while(!FileList.eof()){
      FileList>>filename;

      if(filename!=""){
	dataset.push_back(filename);
	if(filename.Contains("run")){legends.push_back("Data"); colors.push_back(kBlack);datafiles++;}
	else if(filename.Contains("ttbarsignal")){legends.push_back("t#bar{t} Signal"); colors.push_back(kRed+1);}
	else if(filename.Contains("ttbarbg")){legends.push_back("t#bar{t} Other"); colors.push_back(kRed-7);}
	else if(filename.Contains("single")){legends.push_back("Single Top"); colors.push_back(kMagenta);}
	else if(filename.Contains("ww") ||filename.Contains("wz")||filename.Contains("zz")){legends.push_back("Diboson"); colors.push_back(10);}
	else if(filename.Contains("dytautau")){legends.push_back("Z / #gamma* #rightarrow #tau#tau"); colors.push_back(kAzure+8);}
	else if(filename.Contains("dymumu")||filename.Contains("dyee")){legends.push_back("Z / #gamma* #rightarrow ee/#mu#mu"); colors.push_back(kAzure-2);}
//	else if(filename.Contains("dyee")){legends.push_back("Z / #gamma* #rightarrow ee"); colors.push_back(kAzure-2);}
	else if(filename.Contains("wtolnu")){legends.push_back("W+Jets"); colors.push_back(kGreen-3);}
	else if(filename.Contains("qcd")){legends.push_back("QCD Multijet"); colors.push_back(kYellow);}
      }
    }
  }

  if(channel=="combined"){
    datafiles=0;
    dataset.clear();
    legends.clear();
    colors.clear();
    ifstream FileList("FileLists/HistoFileList_Nominal_combined.txt");
    TString filename;

    while(!FileList.eof()){
      FileList>>filename;

      if(filename!=""){
	dataset.push_back(filename);
	
	if(filename.Contains("run")){legends.push_back("Data"); colors.push_back(kBlack);datafiles++;}
	else if(filename.Contains("ttbarsignal")){legends.push_back("t#bar{t} Signal"); colors.push_back(kRed+1);}
	else if(filename.Contains("ttbarbg")){legends.push_back("t#bar{t} Other"); colors.push_back(kRed-7);}
	else if(filename.Contains("single")){legends.push_back("Single Top"); colors.push_back(kMagenta);}
	else if(filename.Contains("ww") ||filename.Contains("wz")||filename.Contains("zz")){legends.push_back("Diboson"); colors.push_back(10);}
	else if(filename.Contains("dytautau")){legends.push_back("Z / #gamma* #rightarrow #tau#tau"); colors.push_back(kAzure+8);}
	else if(filename.Contains("dymumu")||filename.Contains("dyee")){legends.push_back("Z / #gamma* #rightarrow ee/#mu#mu"); colors.push_back(kAzure-2);}
	else if(filename.Contains("wtolnu")){legends.push_back("W+Jets"); colors.push_back(kGreen-3);}
	else if(filename.Contains("qcd")){legends.push_back("QCD Multijet"); colors.push_back(kYellow);}
      }
    }
  }
}

void Plotter::setSystDataSet(TString systematic, TString systSampleUp, TString systSampleDown)
{
  initialized=false;

  DYEntry = "Z / #gamma* #rightarrow ee/#mu#mu";
  lumi=4966;
  if(channel=="ee" || channel=="emu" || channel=="mumu"){  
    TString HistoFileUp,HistoFileDown; 
    if(systematic=="DY_" || systematic=="BG_"){
      HistoFileUp = "FileLists/HistoFileList_Nominal_"+channel+".txt";
      HistoFileDown = "FileLists/HistoFileList_Nominal_"+channel+".txt";
    }else{
      HistoFileUp = "FileLists/HistoFileList_"+systSampleUp+"_"+channel+".txt";
      HistoFileDown = "FileLists/HistoFileList_"+systSampleDown+"_"+channel+".txt";
    }
    ifstream FileListUp(HistoFileUp);
    ifstream FileListDown(HistoFileDown);
    TString filenameUp, filenameDown;

    datasetUp.clear();
    datasetDown.clear();
    legendsUp.clear();
    legendsDown.clear();
    colorsUp.clear();
    colorsDown.clear(); 
    while(!FileListUp.eof()){
      FileListUp>>filenameUp; 

      if(filenameUp!=""){
	datasetUp.push_back(filenameUp);
	if(filenameUp.Contains("run")){legendsUp.push_back("Data"); colorsUp.push_back(kBlack);}
	else if(filenameUp.Contains("ttbarsignal")){legendsUp.push_back("t#bar{t} Signal"); colorsUp.push_back(kRed+1);}
	else if(filenameUp.Contains("ttbarbg")){legendsUp.push_back("t#bar{t} Other"); colorsUp.push_back(kRed-7);}
	else if(filenameUp.Contains("single")){legendsUp.push_back("Single Top"); colorsUp.push_back(kMagenta);}
	else if(filenameUp.Contains("ww") ||filenameUp.Contains("wz")||filenameUp.Contains("zz")){legendsUp.push_back("Diboson"); colorsUp.push_back(10);}
	else if(filenameUp.Contains("dytautau")){legendsUp.push_back("Z / #gamma* #rightarrow #tau#tau"); colorsUp.push_back(kAzure+8);}
	else if(filenameUp.Contains("dymumu")||filenameUp.Contains("dyee")){legendsUp.push_back("Z / #gamma* #rightarrow ee/#mu#mu"); colorsUp.push_back(kAzure-2);}
	else if(filenameUp.Contains("wtolnu")){legendsUp.push_back("W+Jets"); colorsUp.push_back(kGreen-3);}
	else if(filenameUp.Contains("qcd")){legendsUp.push_back("QCD Multijet"); colorsUp.push_back(kYellow);}
      }
    }
    while(!FileListDown.eof()){
      FileListDown>>filenameDown;
      if(filenameDown!=""){
	datasetDown.push_back(filenameDown);
	if(filenameDown.Contains("run")){legendsDown.push_back("Data"); colorsDown.push_back(kBlack);}
	else if(filenameDown.Contains("ttbarsignal")){legendsDown.push_back("t#bar{t} Signal"); colorsDown.push_back(kRed+1);}
	else if(filenameDown.Contains("ttbarbg")){legendsDown.push_back("t#bar{t} Other"); colorsDown.push_back(kRed-7);}
	else if(filenameDown.Contains("single")){legendsDown.push_back("Single Top"); colorsDown.push_back(kMagenta);}
	else if(filenameDown.Contains("ww") ||filenameDown.Contains("wz")||filenameDown.Contains("zz")){legendsDown.push_back("Diboson"); colorsDown.push_back(10);}
	else if(filenameDown.Contains("dytautau")){legendsDown.push_back("Z / #gamma* #rightarrow #tau#tau"); colorsDown.push_back(kAzure+8);}
	else if(filenameDown.Contains("dymumu")||filenameDown.Contains("dyee")){legendsDown.push_back("Z / #gamma* #rightarrow ee/#mu#mu"); colorsDown.push_back(kAzure-2);}
	else if(filenameDown.Contains("wtolnu")){legendsDown.push_back("W+Jets"); colorsDown.push_back(kGreen-3);}
	else if(filenameDown.Contains("qcd")){legendsDown.push_back("QCD Multijet"); colorsDown.push_back(kYellow);}
      }
    }
  }

  if(channel=="combined"){  
    TString HistoFileUp,HistoFileDown; 
    if(systematic=="DY_" || systematic=="BG_"){
      HistoFileUp = "FileLists/HistoFileList_Nominal_"+channel+".txt";
      HistoFileDown = "FileLists/HistoFileList_Nominal_"+channel+".txt";
    }else{
      HistoFileUp = "FileLists/HistoFileList_"+systSampleUp+"_"+channel+".txt";
      HistoFileDown = "FileLists/HistoFileList_"+systSampleDown+"_"+channel+".txt";
    }
    ifstream FileListUp(HistoFileUp);
    ifstream FileListDown(HistoFileDown);
    TString filenameUp, filenameDown;

    datasetUp.clear();
    datasetDown.clear();
    legendsUp.clear();
    legendsDown.clear();
    colorsUp.clear();
    colorsDown.clear();
  
    while(!FileListUp.eof()){
      FileListUp>>filenameUp;

      if(filenameUp!=""){
	datasetUp.push_back(filenameUp);
	if(filenameUp.Contains("run")){legendsUp.push_back("Data"); colorsUp.push_back(kBlack);}
	else if(filenameUp.Contains("ttbarsignal")){legendsUp.push_back("t#bar{t} Signal"); colorsUp.push_back(kRed+1);}
	else if(filenameUp.Contains("ttbarbg")){legendsUp.push_back("t#bar{t} Other"); colorsUp.push_back(kRed-7);}
	else if(filenameUp.Contains("single")){legendsUp.push_back("Single Top"); colorsUp.push_back(kMagenta);}
	else if(filenameUp.Contains("ww") ||filenameUp.Contains("wz")||filenameUp.Contains("zz")){legendsUp.push_back("Diboson"); colorsUp.push_back(10);}
	else if(filenameUp.Contains("dytautau")){legendsUp.push_back("Z / #gamma* #rightarrow #tau#tau"); colorsUp.push_back(kAzure+8);}
	else if(filenameUp.Contains("dymumu")||filenameUp.Contains("dyee")){legendsUp.push_back("Z / #gamma* #rightarrow ee/#mu#mu"); colorsUp.push_back(kAzure-2);}
	else if(filenameUp.Contains("wtolnu")){legendsUp.push_back("W+Jets"); colorsUp.push_back(kGreen-3);}
	else if(filenameUp.Contains("qcd")){legendsUp.push_back("QCD Multijet"); colorsUp.push_back(kYellow);}
      }
    }
    while(!FileListDown.eof()){
      FileListDown>>filenameDown;
      if(filenameDown!=""){
	datasetDown.push_back(filenameDown);
	if(filenameDown.Contains("run")){legendsDown.push_back("Data"); colorsDown.push_back(kBlack);}
	else if(filenameDown.Contains("ttbarsignal")){legendsDown.push_back("t#bar{t} Signal"); colorsDown.push_back(kRed+1);}
	else if(filenameDown.Contains("ttbarbg")){legendsDown.push_back("t#bar{t} Other"); colorsDown.push_back(kRed-7);}
	else if(filenameDown.Contains("single")){legendsDown.push_back("Single Top"); colorsDown.push_back(kMagenta);}
	else if(filenameDown.Contains("ww") ||filenameDown.Contains("wz")||filenameDown.Contains("zz")){legendsDown.push_back("Diboson"); colorsDown.push_back(10);}
	else if(filenameDown.Contains("dytautau")){legendsDown.push_back("Z / #gamma* #rightarrow #tau#tau"); colorsDown.push_back(kAzure+8);}
	else if(filenameDown.Contains("dymumu")||filenameDown.Contains("dyee")){legendsDown.push_back("Z / #gamma* #rightarrow ee/#mu#mu"); colorsDown.push_back(kAzure-2);}
	else if(filenameDown.Contains("wtolnu")){legendsDown.push_back("W+Jets"); colorsDown.push_back(kGreen-3);}
	else if(filenameDown.Contains("qcd")){legendsDown.push_back("QCD Multijet"); colorsDown.push_back(kYellow);}
      }
    }
  }
}


void Plotter::fillHisto()
{
  TH1::AddDirectory(kFALSE);
  if(!initialized){ 
    hists.clear();
    for(unsigned int i=0; i<dataset.size(); i++){
      TFile *ftemp = TFile::Open(dataset[i]);
      TH1D *hist = (TH1D*)ftemp->Get(name)->Clone();
      if(name.Contains("Lepton")){
      	TString stemp = name;
	    stemp.ReplaceAll("Lepton",6,"AntiLepton",10);
	    TH1D *hist2 = (TH1D*)ftemp->Get(stemp)->Clone();     
      	hist->Add(hist2);
      }
      if(name.Contains("BJet")){
      	TString stemp = name;
	    stemp.ReplaceAll("BJet",4,"AntiBJet",8);
	    TH1D *hist2 = (TH1D*)ftemp->Get(stemp)->Clone();     
      	hist->Add(hist2);
      }
      if(name.Contains("Top")){
      	TString stemp = name;
	    stemp.ReplaceAll("Top",3,"AntiTop",7);
	    TH1D *hist2 = (TH1D*)ftemp->Get(stemp)->Clone();     
      	hist->Add(hist2);
      }
      
      
      
      
      setHHStyle(*gStyle);
      //gStyle->SetErrorX(0);

      //cout << "[Plotter::fillHisto] datasetname = " << dataset[i] << endl;

      hists.push_back(*hist);
      //ftemp->Close("R");
      //delete hist;
      delete ftemp;
    }
    initialized=true;
  }
}


void Plotter::ApplyMCATNLOWeight(TH1* hist, TString Systematic, TString Shift, TString Sample)
{
    // This Scale Factor is needed, because the lumiWeight
    // which is stored on the Ntupels does not account
    // for the MC weights that MC@NLO uses.
    
    
    // First, we test if we need to apply the weight in the first place
    bool doApplyMCATNLOWeight = false; 
    if ( (Systematic == "HAD") && (Shift == "Up") && Sample.Contains("ttbar") ) doApplyMCATNLOWeight = true;
    
    // Exit, if nothing needs to be done
    if ( doApplyMCATNLOWeight == false ) return;
    
    
    // Here, we calculate the factor.
    
    // Take the number of non-weighted generated events which comes
    // from the file unmerged/MCATNLO/ee_ttbarsignalplustau.txt
    // which is
    double MCATNLO_Events = 21745199.;
    //
    // The absolut value of the weight for all these is
    double MCATNLO_Weight = 190.41256;
    //
    // There is a fraction with positive weights
    double MCATNLO_posWeights = 0.8865;
    //
    // And a fraction with negative weights
    double MCATNLO_negWeights = 0.1135;
    //
    // Such that the weighted number of events is 
    double MCATNLO_Weights = MCATNLO_Weight * MCATNLO_Events * (MCATNLO_posWeights - MCATNLO_negWeights);
    // 
    // Therefore, the scale factor to be applied is
    double MCATNLO_ScaleFactor = MCATNLO_Weights / MCATNLO_Events;
    double MCATNLO_ScaleFactorInv = MCATNLO_Events / MCATNLO_Weights;
    
    
    // Output
    cout << endl; 
    cout << endl; 
    cout << "ATTENTION!" << endl;
    cout << "Applying a scale factor to the MC@NLO Sample to account for MC Weights" << endl;  
    cout << "    Histo Name:           " << hist->GetName() << endl;
    cout << "    Histo Title:          " << hist->GetTitle() << endl;
    cout << "    Systematic:           " << Systematic << endl;
    cout << "    Shift:                " << Shift << endl;
    cout << "    Sample:               " << Sample << endl;
    cout << "    Channel:              " << channel << endl;
    cout << "    Quantity:             " << name << endl; 
    cout << "    The factor is:        " << MCATNLO_ScaleFactor << endl;
    cout << "    The inverse of it is: " << MCATNLO_ScaleFactorInv << endl;
    cout << endl;
  
    // Apply the weight
    hist->Scale(MCATNLO_ScaleFactorInv);
  
}

void Plotter::fillSystHisto()
{
	 
  TH1::AddDirectory(kFALSE);
  if(!initialized){ 
    hists.clear();
    systhistsUp.clear();
    systhistsDown.clear();
    
    // Fill nominal sample
    for(unsigned int i=0; i<dataset.size(); i++){
      TFile *ftemp = TFile::Open(dataset[i]);
      TH1D *hist = (TH1D*)ftemp->Get(name)->Clone();
      
 
      
      if(name.Contains("Lepton")){
      	TString stemp = name;
	    stemp.ReplaceAll("Lepton",6,"AntiLepton",10);
	    TH1D *hist2 = (TH1D*)ftemp->Get(stemp)->Clone();     
      	hist->Add(hist2);
      }
      if(name.Contains("BJet")){
      	TString stemp = name;
	    stemp.ReplaceAll("BJet",4,"AntiBJet",8);
	    TH1D *hist2 = (TH1D*)ftemp->Get(stemp)->Clone();     
      	hist->Add(hist2);
      }
      if(name.Contains("Top")){
      	TString stemp = name;
	    stemp.ReplaceAll("Top",3,"AntiTop",7);
	    TH1D *hist2 = (TH1D*)ftemp->Get(stemp)->Clone();     
      	hist->Add(hist2);
      }
                 
     
      setHHStyle(*gStyle);
      //gStyle->SetErrorX(0);
      hists.push_back(*hist);
      //ftemp->Close("R");
      //delete hist;
      delete ftemp;
    }
    
    
    // Fill up-shifted sample 
    for(unsigned int i=0; i<datasetUp.size(); i++){ 
      TFile *ftemp = TFile::Open(datasetUp[i]);
      TH1D *hist = (TH1D*)ftemp->Get(name)->Clone();     
      if(name.Contains("BJet")){
      	TString stemp = name;
	    stemp.ReplaceAll("BJet",4,"AntiBJet",8);
	    TH1D *hist2 = (TH1D*)ftemp->Get(stemp)->Clone();     
      	hist->Add(hist2);
      }
      if(name.Contains("Lepton")){
      	TString stemp = name;
	    stemp.ReplaceAll("Lepton",6,"AntiLepton",10);
	    TH1D *hist2 = (TH1D*)ftemp->Get(stemp)->Clone();     
      	hist->Add(hist2);
      }
      if(name.Contains("Top")){
      	TString stemp = name;
	    stemp.ReplaceAll("Top",3,"AntiTop",7);
	    TH1D *hist2 = (TH1D*)ftemp->Get(stemp)->Clone();     
      	hist->Add(hist2);
      }
      
      
      
      // Apply Scale Factor for MC@NLO
      if ( legendsSyst.size() > 0 ) { 
      	ApplyMCATNLOWeight(hist, legendsSyst.back(), "Up", datasetUp[i]);
      } 
      
      setHHStyle(*gStyle);
      //gStyle->SetErrorX(0);
      systhistsUp.push_back(*hist);
      delete ftemp;
    }
    
    // Fill down-shifted sample
    for(unsigned int i=0; i<datasetDown.size(); i++){
      TFile *ftemp = TFile::Open(datasetDown[i]);
      TH1D *hist = (TH1D*)ftemp->Get(name)->Clone();     
      if(name.Contains("BJet")){
      	TString stemp = name;
	    stemp.ReplaceAll("BJet",4,"AntiBJet",8);
	    TH1D *hist2 = (TH1D*)ftemp->Get(stemp)->Clone();     
      	hist->Add(hist2);
      }
      if(name.Contains("Lepton")){
      	TString stemp = name;
	    stemp.ReplaceAll("Lepton",6,"AntiLepton",10);
	    TH1D *hist2 = (TH1D*)ftemp->Get(stemp)->Clone();     
       	hist->Add(hist2);
      }
      if(name.Contains("Top")){
      	TString stemp = name;
	    stemp.ReplaceAll("Top",3,"AntiTop",7);
	    TH1D *hist2 = (TH1D*)ftemp->Get(stemp)->Clone();     
      	hist->Add(hist2);
      }
      
      // Apply Scale Factor for MC@NLO
      if ( legendsSyst.size() > 0 ) {
          ApplyMCATNLOWeight(hist, legendsSyst.back(), "Down", datasetDown[i]);
      }

      
      setHHStyle(*gStyle);
      //gStyle->SetErrorX(0);
      systhistsDown.push_back(*hist);
      delete ftemp;
    }
    initialized=true;
  }
}



void Plotter::write() // do scaling, stacking, legending, and write in file 
{
  if(initialized){

  TCanvas * c = new TCanvas("","");

  THStack * stack=  new THStack("def", "def");
  TLegend * leg =  new TLegend(0.70,0.55,0.98,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetX1NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.25);
  leg->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 - leg->GetNRows()*0.04);
  leg->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength());
  leg->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength());
  TH1D *drawhists[hists.size()];



  std::stringstream ss;
  ss << DYScale[channelType];
  TString scale;
  scale=(TString)ss.str();
  int legchange=0;
  leg->Clear();
  c->Clear();
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  c->SetName("");
  c->SetTitle("");
  for(unsigned int i=0; i<hists.size() ; i++){ // prepare histos and leg
    drawhists[i]=(TH1D*) hists[i].Clone();
    if(rebin>1) drawhists[i]->Rebin(rebin);
    setStyle(*drawhists[i], i);
    if(legends[i] != "Data"){
      if(legends[i] == legends[0]){
	    drawhists[0]->Add(drawhists[i]);
      }

      if(legends[i] == "t#bar{t} Signal"){signalHist = i;}
      if((legends[i] == DYEntry) && channelType!=2 ){
	    drawhists[i]->Scale(DYScale[channelType]);
      }
      if(i > 1){
	    if(legends[i] != legends[i-1]){
	      legchange = i; 
	      if((legends[i] == DYEntry)&& DYScale[channelType] != 1) leg->AddEntry(drawhists[i], legends[i],"f");
	      else leg->AddEntry(drawhists[i], legends[i] ,"f");
	    }else{
	      drawhists[legchange]->Add(drawhists[i]);	      
	    }
      }

      if(i!=(hists.size()-1)){
	if(legends[i]!=legends[i+1]){
	  //cout<<legends[i]<<endl;
	  drawhists[i]->SetLineColor(1);
	}
      }else{
	drawhists[i]->SetLineColor(1);
      }


      if(legends[i] != legends[i-1]){
	drawhists[i]->SetLineColor(1);
	stack->Add(drawhists[i]); 
      }
    }
    else{
      if(i==0) leg->AddEntry(drawhists[i], legends[i] ,"pe");
      if(i>0){
	    if(legends[i] != legends[i-1]){
	      leg->AddEntry(drawhists[i], legends[i] ,"pe");
	    }
	    if(legends[i] == legends[0]){
	      drawhists[0]->Add(drawhists[i]);
	    }
      }
    }
  } 
  
  leg = ControlLegend(hists.size(), drawhists, legends, leg);
  
  TFile *f0 = new TFile("SigBackground.root","UPDATE");

  TList* l = stack->GetHists(); 
  TH1D* stacksum = (TH1D*) l->At(0)->Clone();
 
  TString aaa = "a";
  for (int i = 1; i < l->GetEntries(); ++i) {
    aaa=aaa+"a"; 
    stacksum->Add((TH1D*)l->At(i));
    if(legends[datafiles+i] == "t#bar{t} Other") {
      stacksum->Write(name+"_"+channel+aaa+"_Background");
    }
    if(legends[datafiles+i] == "t#bar{t} Signal") {
	stacksum->Write(name+"_"+channel+aaa+"_Signal");
    }

  } 
  f0->Close();
  //stat uncertainty::make a function 
  TH1D* syshist =0;
  syshist = (TH1D*)stacksum->Clone();
  double lumierr = 0.045; 
  for(Int_t i=0; i<=syshist->GetNbinsX(); ++i){
    
    Double_t binc = 0;
    binc += stacksum->GetBinContent(i);
    syshist->SetBinContent(i, binc);
    // calculate uncertainty: lumi uncertainty
    Double_t binerr2 = binc*binc*lumierr*lumierr;
    Double_t topunc = 0; // uncertainty on top xsec
    
    //Kidonakis
    double topxsec = TOPXSEC;
    double topxsecErr2 = 2.2*2.2 + 4.4*4.4 + 5.5*5.5; //topxsecErr2 = lumiErr*lumiErr + topxsecScaleErr*topxsecScaleErr + topxsecPDFErr*topxsecPDFErr

    double topRelUnc =  TMath::Sqrt(topxsecErr2)/topxsec;
    topunc += drawhists[signalHist]->GetBinContent(i)*topRelUnc;
    binerr2 += (topunc*topunc);
    syshist->SetLineColor(1);
    syshist->SetBinError(i, TMath::Sqrt(binerr2));
  }    
  
  if(logY)c->SetLogy();
  syshist->SetFillStyle(3004);
  syshist->SetFillColor(kBlack);
  leg->AddEntry( syshist, "Uncertainty", "f" );
  
  drawhists[0]->SetMinimum(ymin);
  if(rangemin!=0 || rangemax!=0)drawhists[0]->SetAxisRange(rangemin, rangemax, "X");
  if(ymax==0){
    if(logY){  
      drawhists[0]->SetMaximum(18*drawhists[0]->GetBinContent(drawhists[0]->GetMaximumBin()));
    }
    else{drawhists[0]->SetMaximum(1.5*drawhists[0]->GetBinContent(drawhists[0]->GetMaximumBin()));}
  }else{
    drawhists[0]->SetMaximum(ymax);
  }
  drawhists[0]->GetXaxis()->SetNoExponent(kTRUE);
  TGaxis::SetMaxDigits(2);
  //Removal of extra ticks in JetMult plots
  if(name.Contains("jet") && name.Contains("Multi")){ drawhists[0]->GetXaxis()->SetNdivisions(drawhists[0]->GetNbinsX(),0,0, 1);}
  //Add the binwidth to the yaxis in yield plots
  TString ytitle = TString(drawhists[0]->GetYaxis()->GetTitle()).Copy();
  double binwidth = drawhists[0]->GetXaxis()->GetBinWidth(1);
  std::ostringstream width;
  width<<binwidth;
  if(name.Contains("Rapidity") || name.Contains("Eta")){ytitle.Append(" / ").Append(width.str());}
  else if(name.Contains("pT") || name.Contains("Mass") || name.Contains("mass") || name.Contains("MET") || name.Contains("HT")){ytitle.Append(" / ").Append(width.str()).Append(" GeV");};
  drawhists[0]->GetYaxis()->SetTitle(ytitle);
  
  drawhists[0]->Draw("e1"); //############## 
  //drawhists[0]->Draw("e"); //############## 
  
  stack->Draw("same HIST");
  gPad->RedrawAxis();
  TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0.5)");//this is frustrating and stupid but apparently necessary...
  setex1->Draw();
  syshist->SetMarkerStyle(0);//<===================
  syshist->Draw("same,E2");
  TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.)");
  setex2->Draw();
  drawhists[0]->Draw("same,e1"); //#############
  //drawhists[0]->Draw("same,e"); //#############
  
  DrawCMSLabels(false, lumi);
  
  DrawDecayChLabel(channelLabel[channelType]);    
  leg->Draw("SAME");  
  drawRatio(drawhists[0], stacksum, 0.5, 1.9, *gStyle);

    
  // Create Directory for Output Plots 
  gSystem->MakeDirectory(outpathPlots);
  int cntMkdirTries = 0;
  gSystem->MakeDirectory(outpathPlots+"/"+subfolderChannel+"/"+subfolderSpecial);  
  c->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/"+name+".eps");  
  c->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/"+name+".C");  
  c->Clear();  
  leg->Clear();  
  delete c;  
  delete leg;
  
  //delete stack;
  }
  else std::cout << "Histogram " << name << " not filled during the process." << std::endl;
}

void Plotter::setStyle(TH1 &hist, unsigned int i)
{
  hist.SetFillColor(colors[i]);
  hist.SetLineColor(colors[i]);
  hist.SetLineWidth(1);
  

  if(legends[i] == "Data"){
    hist.SetFillColor(0);
    hist.SetMarkerStyle(20); 
    hist.SetMarkerSize(1.);
    hist.SetLineWidth(1);
    hist.GetXaxis()->SetLabelFont(42);
    hist.GetYaxis()->SetLabelFont(42);
    hist.GetXaxis()->SetTitleFont(42);
    hist.GetYaxis()->SetTitleFont(42);
    hist.GetYaxis()->SetTitleOffset(1.7);
    hist.GetXaxis()->SetTitleOffset(1.25);
    if(name.Contains("pT") ||name.Contains("Mass") ){
      hist.GetXaxis()->SetTitle(XAxis+" #left[GeV#right]");
      if(name.Contains("Rapidity")) hist.GetXaxis()->SetTitle(XAxis);
    }else  hist.GetXaxis()->SetTitle(XAxis);
    
    if(name.Contains("pT") ||name.Contains("Mass")){
      hist.GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d"+XAxis+"}"+" #left[GeV^{-1}#right]"); 
      if(name.Contains("Rapidity")) hist.GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d"+XAxis+"}");     
    }else{hist.GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d"+XAxis+"}");
    }
  }
}

void Plotter::setStyle(TH1D &hist, unsigned int i)
{
  hist.SetFillColor(colors[i]);
  hist.SetLineColor(colors[i]);
  

  if(legends[i] == "Data"){
    hist.SetMarkerStyle(20); 
    hist.SetMarkerSize(1.);
    hist.SetLineWidth(1);
    hist.GetXaxis()->SetLabelFont(42);
    hist.GetYaxis()->SetLabelFont(42);
    hist.GetXaxis()->SetTitleSize(0.04);
    hist.GetYaxis()->SetTitleSize(0.04);
    hist.GetXaxis()->SetTitleFont(42);
    hist.GetYaxis()->SetTitleFont(42);
    if(name.Contains("pT") || name.Contains("Mass")){
      hist.GetXaxis()->SetTitle(XAxis+" #left[GeV#right]");
      if(name.Contains("Rapidity")) hist.GetXaxis()->SetTitle(XAxis);
    }else    hist.GetXaxis()->SetTitle(XAxis);
    //hist.GetXaxis()->SetTitle(XAxis);
    hist.GetYaxis()->SetTitle(YAxis);
    hist.GetYaxis()->SetTitleOffset(1.7);
    hist.GetXaxis()->SetTitleOffset(1.25);
  }
}

void Plotter::PlotXSec(){

  TH1::AddDirectory(kFALSE);
  for(int i =0; i<16; i++){
//  for(int i =0; i<15; i++){
    InclusiveXsectionSysErrorBySyst[channelType][i] = 0.;
  }
  
  
  if(doSystematics){
    cout << endl;
    cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
    cout << "Starting Calculation of Inclusive Systematics for '" << name << "' in Channel '" << channel << "':" << endl;
    cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
    CalcInclSystematics("JES","JESUP", "JESDOWN", 0);
    CalcInclSystematics("RES", "RESUP", "RESDOWN", 1);
    CalcInclSystematics("PU_", "PU_UP", "PU_DOWN", 2);
    CalcInclSystematics("SCALE", "SCALEUP", "SCALEDOWN", 3);
    CalcInclSystematics("MATCH", "MATCHUP", "MATCHDOWN", 4);
    CalcInclSystematics("MASS", "MASSUP", "MASSDOWN", 5);
    CalcInclSystematics("DY_", "DY_UP", "DY_DOWN", 6);
    CalcInclSystematics("BG_","BG_UP", "BG_DOWN", 7);
    CalcInclSystematics("HAD", "MCATNLO", "POWHEG", 8);
    CalcInclSystematics("BTAG_", "BTAG_UP", "BTAG_DOWN", 9);
    InclFlatSystematics(10);
      
    cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
    cout << "Finished Calculation of Inclusive Systematics for '" << name << "' in Channel '" << channel << "':" << endl;
    cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
    cout << endl;
  }
  
  CalcXSec(dataset, InclusiveXsection, InclusiveXsectionStatError, "","");

  double syst_square=0;


  for(int i =0; i<15; i++){
    syst_square += InclusiveXsectionSysErrorBySyst[channelType][i]*InclusiveXsectionSysErrorBySyst[channelType][i];
  }
  InclusiveXsectionSysError[channelType] = sqrt(syst_square);
  cout<<"&^&^&^&^&^&^^&^&^ InclusiveXsectionSysError[channelType]: "<<InclusiveXsectionSysError[channelType]<<endl;

  if(channelType==3){

     // measured results with statistical error
   Double_t mx[]   = {      0.50,       1.50,       2.50,       3.50};
   Double_t mexl[] = {      0.00,       0.00,       0.00,       0.00};
   Double_t mexh[] = {      0.00,       0.00,       0.00,       0.00};

   TGraphAsymmErrors *mplot = new TGraphAsymmErrors(4, mx, InclusiveXsection, mexl, mexh,InclusiveXsectionStatError, InclusiveXsectionStatError);
   mplot->SetMarkerStyle(20);
   mplot->SetMarkerColor(kBlack);
   mplot->SetMarkerSize(1.5);
   mplot->SetLineColor(kBlack);
   
   for(int i=0; i<4; i++){
     InclusiveXsectionTotalError[i] = sqrt(InclusiveXsectionStatError[i]*InclusiveXsectionStatError[i] +InclusiveXsection[i]*InclusiveXsectionSysError[i]*InclusiveXsection[i]*InclusiveXsectionSysError[i]);
   }

   TGraphAsymmErrors *mplotwithsys = new TGraphAsymmErrors(4, mx, InclusiveXsection, mexl, mexh,InclusiveXsectionTotalError, InclusiveXsectionTotalError);
   mplotwithsys->SetMarkerStyle(20);
   mplotwithsys->SetMarkerColor(kBlack);
   mplotwithsys->SetMarkerSize(1.5);
   mplotwithsys->SetLineColor(kBlack);

   // mstw
   Double_t mstwmean = 157.5;
   Double_t mstwx[]   = {    -0.5,     0.5,	1.5,	 2.5,	  3.5,     4.5};
   Double_t mstwy[]   = {mstwmean,mstwmean,mstwmean,mstwmean,mstwmean,mstwmean};
   Double_t mstwexl[] = {      .4,	.4,	 .5,	  .5,	   .5,      .5};
   Double_t mstwexh[] = {      .5,	.5,	 .5,	  .5,	   .4,      .4};
   Double_t mstweyl[] = {    24.4,    24.4,    24.4,	24.4,	 24.4,    24.4};
   Double_t mstweyh[] = {    23.2,    23.2,    23.2,	23.2,	 23.2,    23.2};

   TGraphAsymmErrors *mstwplot = new TGraphAsymmErrors(6, mstwx, mstwy, mstwexl, mstwexh, mstweyl, mstweyh);
   mstwplot->SetLineColor(kGreen+1);
   mstwplot->SetLineWidth(4);
   mstwplot->SetFillColor(kGreen+1);
   mstwplot->SetFillStyle(3004);

   // herapdf
   Double_t heramean = 144.156;
   Double_t herapdfx[]   = {	-0.5,	  0.5,     1.5,     2.5,     3.5,     4.5};
   Double_t herapdfy[]   = {heramean,heramean,heramean,heramean,heramean,heramean};
   Double_t herapdfexl[] = {	  .4,	   .4,      .5,      .5,      .5,      .5};
   Double_t herapdfexh[] = {	  .5,	   .5,      .5,      .5,      .4,      .4};
   Double_t herapdfeyl[] = {  13.849,  13.849,  13.849,  13.849,  13.849,  13.849};
   Double_t herapdfeyh[] = {   5.475,	5.475,   5.475,   5.475,   5.475,   5.475};

   TGraphAsymmErrors *herapdfplot = new TGraphAsymmErrors(6, herapdfx, herapdfy, herapdfexl, herapdfexh, herapdfeyl, herapdfeyh);
   herapdfplot->SetLineColor(kBlue+1);
   herapdfplot->SetLineWidth(4);
   herapdfplot->SetFillColor(kBlue+1);
   herapdfplot->SetFillStyle(3005);

   TH1F* framehist = new TH1F("framehist","",4,0.,4.);
   framehist->SetMinimum(0);
   framehist->SetMaximum(310);
   framehist->GetXaxis()->SetTickLength(0);
   framehist->GetXaxis()->SetBinLabel(1,"");
   framehist->GetXaxis()->SetBinLabel(2,"");
   framehist->GetXaxis()->SetBinLabel(3,"");
   framehist->GetYaxis()->SetTitle("#sigma [pb]");
   framehist->GetYaxis()->CenterTitle(kTRUE);

   TPaveText* box1 = new TPaveText(0.25,0.33,0.33,0.43,"NDC");
   box1->SetFillColor(10);
   box1->SetTextSize(0.04);
   box1->AddText("ee");

   TPaveText* box2 = new TPaveText(0.44,0.33,0.52,0.43,"NDC");
   box2->SetFillColor(10);
   box2->SetTextSize(0.04);
   box2->AddText("#mu#mu");

   TPaveText* box3 = new TPaveText(0.62,0.33,0.72,0.43,"NDC");
   box3->SetFillColor(10);
   box3->SetTextSize(0.04);
   box3->AddText("e#mu");

   TPaveText* box4 = new TPaveText(0.82,0.33,0.90,0.43,"NDC");
   box4->SetFillColor(10);
   box4->SetTextSize(0.04);
   box4->AddText("combined");

   TLegend* leg = getNewLegend(); // new TLegend( 0.56, 0.18, 0.89, 0.33 );
   leg->SetBorderSize( 0 );
   leg->SetFillColor( 0 );
   leg->SetTextFont(62);
   leg->SetTextSize(0.03);
   leg->AddEntry( mplot,       "Measurements",            "p"  );
   leg->AddEntry( mstwplot,    "MCFM #otimes MSTW08",     "lf" );
   leg->AddEntry( herapdfplot, "MCFM #otimes HERAPDF1.0", "lf" );

   TCanvas* c = new TCanvas("plot", "plot", 1200, 800);
   framehist->Draw();
   herapdfplot->Draw("C,2,SAME");
   mstwplot->Draw("C,2,SAME");
   gStyle->SetEndErrorSize(8);
   mplot->Draw("p,SAME");
   mplotwithsys->Draw("p,SAME,Z");
   leg ->Draw("SAME");
   box1->Draw("SAME");
   box2->Draw("SAME");
   box3->Draw("SAME");
   box4->Draw("SAME");
   gSystem->MakeDirectory(outpathPlots);
   gSystem->MakeDirectory(outpathPlots+subfolderChannel);
   gSystem->MakeDirectory(outpathPlots+subfolderChannel+subfolderSpecial);
   c->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/"+"InclusiveXSec.eps");
   c->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/"+"InclusiveXSec.C");
   c->Clear();
   delete c;

   cout<<"!!!!!!!!!!!!!!!!!!!!ee Cross Section: "<<InclusiveXsection[0]<<" +/- "<<InclusiveXsectionStatError[0]<<"(stat) +/- "<<InclusiveXsection[0]*InclusiveXsectionSysError[0]<<"(sys)"<<endl;
   cout<<"!!!!!!!!!!!!!!!!!!!!mumu Cross Section: "<<InclusiveXsection[1]<<" +/- "<<InclusiveXsectionStatError[1]<<"(stat) +/- "<<InclusiveXsection[0]*InclusiveXsectionSysError[1]<<"(sys)"<<endl;
   cout<<"!!!!!!!!!!!!!!!!!!!!emu Cross Section: "<<InclusiveXsection[2]<<" +/- "<<InclusiveXsectionStatError[2]<<"(stat) +/- "<<InclusiveXsection[0]*InclusiveXsectionSysError[2]<<"(sys)"<<endl;
   cout<<"!!!!!!!!!!!!!!!!!!!!Combined Cross Section: "<<InclusiveXsection[3]<<" +/- "<<InclusiveXsectionStatError[3]<<"(stat) +/- "<<InclusiveXsection[0]*InclusiveXsectionSysError[3]<<"(sys)"<<endl;


  }
}

void Plotter::MakeTable(){

  TH1D *numhists5[hists.size()];
  TH1D *numhists6[hists.size()];
  TH1D *numhists7[hists.size()];
  TH1D *numhists8[hists.size()];
  TH1D *numhists9[hists.size()];

  for(unsigned int i=0; i<dataset.size(); i++){
    TFile *ftemp = TFile::Open(dataset[i]);
    TH1D *temp_hist5 = (TH1D*)ftemp->Get("step5")->Clone();     
    numhists5[i]=temp_hist5;
    TH1D *temp_hist6 = (TH1D*)ftemp->Get("step6")->Clone();     
    numhists6[i]=temp_hist6;
    TH1D *temp_hist7 = (TH1D*)ftemp->Get("step7")->Clone();     
    numhists7[i]=temp_hist7;
    TH1D *temp_hist8 = (TH1D*)ftemp->Get("step8")->Clone();     
    numhists8[i]=temp_hist8;
    TH1D *temp_hist9 = (TH1D*)ftemp->Get("step9")->Clone();     
    numhists9[i]=temp_hist9;
    delete ftemp;
  }

  

  for(unsigned int i=0; i<hists.size() ; i++){ // prepare histos and leg
    if((legends[i] == DYEntry) && channelType!=2){
      //numhists5[i]->Scale(DYScale[channelType]);//DYscale not applied in step5 and 6?
      //numhists6[i]->Scale(DYScale[channelType]);
      numhists7[i]->Scale(DYScale[channelType]);
      numhists8[i]->Scale(DYScale[channelType]);
      numhists9[i]->Scale(DYScale[channelType]);
    }
  }  

  ////////////////////////////Make output for tables
  double tmp_num5 = 0;
  double tmp_num6 = 0;
  double tmp_num7 = 0;
  double tmp_num8 = 0;
  double tmp_num9 = 0;
  
  ofstream EventFile5;
  ofstream EventFile6;
  ofstream EventFile7;
  ofstream EventFile8;
  ofstream EventFile9; 
  string EventFilestring = outpathPlots.Data();
  EventFilestring.append(subfolderChannel.Data());
  EventFilestring.append(subfolderSpecial.Data());
  gSystem->MakeDirectory(outpathPlots);
  gSystem->MakeDirectory(outpathPlots+"/"+subfolderChannel+"/"+subfolderSpecial);  
  string EventFilestring5;
  string EventFilestring6;
  string EventFilestring7;
  string EventFilestring8;
  string EventFilestring9; 
  EventFilestring5 =EventFilestring;EventFilestring5.append("/Events5.txt");
  EventFilestring6 =EventFilestring;EventFilestring6.append("/Events6.txt");
  EventFilestring7 =EventFilestring;EventFilestring7.append("/Events7.txt");
  EventFilestring8 =EventFilestring;EventFilestring8.append("/Events8.txt");
  EventFilestring9 =EventFilestring;EventFilestring9.append("/Events9.txt");
  EventFile5.open(EventFilestring5.c_str());
  EventFile6.open(EventFilestring6.c_str());
  EventFile7.open(EventFilestring7.c_str());
  EventFile8.open(EventFilestring8.c_str());
  EventFile9.open(EventFilestring9.c_str());
  double bg_num5 = 0;
  double bg_num6 = 0;
  double bg_num7 = 0;
  double bg_num8 = 0;
  double bg_num9 = 0;
  for(unsigned int i=0; i<hists.size() ; i++){ 
    tmp_num5+=numhists5[i]->Integral();
    tmp_num6+=numhists6[i]->Integral();
    tmp_num7+=numhists7[i]->Integral();
    tmp_num8+=numhists8[i]->Integral();
    tmp_num9+=numhists9[i]->Integral();

    if(i==(hists.size()-1)){
      EventFile5<<legends[i]<<": "<<tmp_num5<<endl;
      EventFile6<<legends[i]<<": "<<tmp_num6<<endl;
      EventFile7<<legends[i]<<": "<<tmp_num7<<endl;
      EventFile8<<legends[i]<<": "<<tmp_num8<<endl;
      EventFile9<<legends[i]<<": "<<tmp_num9<<endl;
      bg_num5+=tmp_num5;
      bg_num6+=tmp_num6;
      bg_num7+=tmp_num7;
      bg_num8+=tmp_num8;
      bg_num9+=tmp_num9;
      tmp_num5=0;
      tmp_num6=0;
      tmp_num7=0;
      tmp_num8=0;
      tmp_num9=0;
    }else if(legends[i]!=legends[i+1]){
      EventFile5<<legends[i]<<": "<<tmp_num5<<endl;
      EventFile6<<legends[i]<<": "<<tmp_num6<<endl;
      EventFile7<<legends[i]<<": "<<tmp_num7<<endl;
      EventFile8<<legends[i]<<": "<<tmp_num8<<endl;
      EventFile9<<legends[i]<<": "<<tmp_num9<<endl;
      if(legends[i]!="Data"){
    	bg_num5+=tmp_num5;
    	bg_num6+=tmp_num6;
    	bg_num7+=tmp_num7;
    	bg_num8+=tmp_num8;
	    bg_num9+=tmp_num9;
      }
      tmp_num5=0;
      tmp_num6=0;
      tmp_num7=0;
      tmp_num8=0;
      tmp_num9=0;
    }

  }
  EventFile5<<"Total background: "<<bg_num5<<endl;
  EventFile5.close();  
  EventFile6<<"Total background: "<<bg_num6<<endl;
  EventFile6.close();  
  EventFile7<<"Total background: "<<bg_num7<<endl;
  EventFile7.close();
  EventFile8<<"Total background: "<<bg_num8<<endl;
  EventFile8.close();  
  EventFile9<<"Total background: "<<bg_num9<<endl;
  EventFile9.close();  
}

double Plotter::CalcXSec(std::vector<TString> datasetVec, double InclusiveXsectionVec[4],double InclusiveXsectionStatErrorVec[4], TString Systematic, TString Shift){
  //  double BranchingFraction[4]={0.0167, 0.0162, 0.0328, 0.06569};//[ee, mumu, emu, combined] including tau
  //  double BranchingFraction[4]={0.011556, 0.011724, 0.022725, 0.04524};//[ee, mumu, emu, combined] not including tau

  
  double BranchingFraction[4]={0.01166, 0.01166, 0.02332, 0.04666};//[ee, mumu, emu, combined] not including tau
  //lumi = 4966;
  lumi = LUMI;

  TH1D *numhists[hists.size()];
  double numbers[4]={0., 0., 0., 0.};//[0]=data, [1]=Signal, [2]Signal(only lumi & PU weights), [3]background (non-ttbar)
  double TTbarBGnum =0;

  for(unsigned int i=0; i<datasetVec.size(); i++){
    TFile *ftemp = TFile::Open(datasetVec[i]);
    TH1D *hist = (TH1D*)ftemp->Get("step9")->Clone();   
    
    // Apply Scale Factor for MC@NLO
    ApplyMCATNLOWeight(hist, Systematic, Shift,  datasetVec[i]); 
    
      
    numhists[i]=hist;
    delete ftemp;
  }
 
  for(unsigned int i=0; i<hists.size() ; i++){ // prepare histos and leg 
    if(legends[i] == "Data"){
      numbers[0]+=numhists[i]->Integral();
    }
    else if(legends[i] == "t#bar{t} Signal"){ 
      TFile *ftemp2 = TFile::Open(datasetVec[i]);  
      TH1D *NoPUPlot = (TH1D*)ftemp2->Get("step9")->Clone(); 
    
      // Apply Scale Factor for MC@NLO
      ApplyMCATNLOWeight(NoPUPlot, Systematic, Shift,  datasetVec[i]);  
      
      numbers[1]+=NoPUPlot->Integral(); 
      delete ftemp2;
       
      TFile *ftemp = TFile::Open(datasetVec[i]);   
      //TH1D *GenPlot = (TH1D*)ftemp->Get("VisGenAll")->Clone();  
      TH1D *GenPlot = (TH1D*)ftemp->Get("GenAll")->Clone();  
      
      // Apply Scale Factor for MC@NLO
      ApplyMCATNLOWeight(GenPlot, Systematic, Shift,  datasetVec[i]); 
      
      numbers[2]+=GenPlot->Integral();  
      delete ftemp; 
      
    }  else if(legends[i] == "t#bar{t} Other"){
      //TFile *ftemp2 = TFile::Open(datasetVec[i]);
      // TH1D *NoPUPlot = (TH1D*)ftemp2->Get("step9")->Clone();
      // TTbarBGnum+=NoPUPlot->Integral();
      //delete ftemp2;
      //}  
    
      
      numbers[3]+=numhists[i]->Integral();
      
    } else {      
      if((legends[i] == DYEntry) && channelType!=2){
	    numhists[i]->Scale(DYScale[channelType]);
      }
      if((legends[i] == DYEntry) && Systematic == "DY_" && Shift == "Up"){
	    numhists[i]->Scale(1.3);
      }
      if((legends[i] == DYEntry) && Systematic == "DY_" && Shift == "Down"){
	    numhists[i]->Scale(0.7);
      }
      if(Systematic == "BG_" && Shift=="Up" && legends[i]!= "t#bar{t} Other" && legends[i] != DYEntry){
	    numhists[i]->Scale(1.3);
      }
      if(Systematic == "BG_" && Shift=="Down" && legends[i]!= "t#bar{t} Other" && legends[i] != DYEntry){
	     numhists[i]->Scale(0.7);
      }  
      
      numbers[3]+=numhists[i]->Integral();
    }   
  }  
 
  ////////////////////////////Make output for tables
  
  double tmp_num = 0;

  double signalFraction = 0; 
  
  signalFraction = numbers[1]/(numbers[1]+TTbarBGnum); // is 1 right now, since TTbarBGnum is 0

  ofstream EventFile;  
  string EventFilestring = outpathPlots.Data();
  EventFilestring.append(subfolderChannel.Data());
  EventFilestring.append(subfolderSpecial.Data());  
  EventFilestring.append("/Events.txt");
  EventFile.open(EventFilestring.c_str());
  double bg_num = 0;
  for(unsigned int i=0; i<hists.size() ; i++){ 
    tmp_num+=numhists[i]->Integral();

    if(i==(hists.size()-1)){
      EventFile<<legends[i]<<": "<<tmp_num<<endl;
      bg_num+=tmp_num;
      tmp_num=0;
    }else if(legends[i]!=legends[i+1]){
      EventFile<<legends[i]<<": "<<tmp_num<<endl;
      if(legends[i]!="Data")bg_num+=tmp_num;
      tmp_num=0;
    }

  }
  EventFile<<"Total background: "<<bg_num<<endl;
  EventFile.close();
 
  //  double xsec = (signalFraction*(numbers[0]-numbers[3]))/((numbers[1]/numbers[2])*BranchingFraction[channelType]*lumi);
  double xsec = ((numbers[0]-numbers[3]))/((numbers[1]/numbers[2])*BranchingFraction[channelType]*lumi);
  double xsecstaterror = TMath::Sqrt(numbers[0])/((numbers[1]/numbers[2])*BranchingFraction[channelType]*lumi);

  if(Systematic == "HAD") { 
    cout<<"numbers[0] (All data)     : "<<numbers[0]<<endl;
    cout<<"numbers[1] (Rec Level MC) : "<<numbers[1]<<endl;
    cout<<"numbers[2] (Gen Level MC) : "<<numbers[2]<<endl;
    cout<<"numbers[3] (Background)   : "<<numbers[3]<<endl;
    cout<<"Global Efficiency: "<<(numbers[1]/numbers[2])<<endl;      
  }
  if(channelType!=3){
    InclusiveXsectionVec[channelType] = xsec;
    InclusiveXsectionStatErrorVec[channelType] = xsecstaterror;
  }else{
    InclusiveXsectionVec[channelType] =( InclusiveXsectionVec[0]/(InclusiveXsectionStatErrorVec[0]*InclusiveXsectionStatErrorVec[0])
				                        +InclusiveXsectionVec[1]/(InclusiveXsectionStatErrorVec[1]*InclusiveXsectionStatErrorVec[1])			
		   		                        +InclusiveXsectionVec[2]/(InclusiveXsectionStatErrorVec[2]*InclusiveXsectionStatErrorVec[2]) )/
				                         ( 1/(InclusiveXsectionStatErrorVec[0]*InclusiveXsectionStatErrorVec[0])
				                        +  1/(InclusiveXsectionStatErrorVec[1]*InclusiveXsectionStatErrorVec[1])			
				                        +  1/(InclusiveXsectionStatErrorVec[2]*InclusiveXsectionStatErrorVec[2])   );			

    InclusiveXsectionStatErrorVec[channelType] =1/(TMath::Sqrt(
                                      (1/(InclusiveXsectionStatErrorVec[0]*InclusiveXsectionStatErrorVec[0]))
							         +(1/(InclusiveXsectionStatErrorVec[1]*InclusiveXsectionStatErrorVec[1]))			
    				                 +(1/(InclusiveXsectionStatErrorVec[2]*InclusiveXsectionStatErrorVec[2]))      ));	 
  } 
  return xsec;
}
void Plotter::CalcDiffXSec(TH1 *varhists[], TH1* RecoPlot, TH1* GenPlot, TH2* genReco2d, double DiffXSecVec[4][10], double DiffXSecStatErrorVec[4][10], double GenDiffXSecVec[4][10]){

  double SignalEvents = 63244696.0;
  double Xbins[XAxisbins.size()];
  double binWidth[XAxisbinCenters.size()];
  double topxsec = TOPXSEC;
  for(unsigned int i = 0; i<XAxisbins.size();i++){Xbins[i]=XAxisbins[i];}
  double DataSum[XAxisbinCenters.size()];
  double GenSignalSum[XAxisbinCenters.size()];
  double BGSum[XAxisbinCenters.size()];
  bool init = false;
  double efficiencies[XAxisbinCenters.size()];
  
  for (unsigned int hist =0; hist<hists.size(); hist++){
    if(legends[hist] == "Data"){
      for (Int_t bin=0; bin<bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
	DataSum[bin]+=varhists[hist]->GetBinContent(bin+1);
      }
    }
    else if((legends[hist] == "t#bar{t} Signal")&&init==false){
      init=true;
      for (Int_t bin=0; bin<bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
	efficiencies[bin] = (RecoPlot->GetBinContent(bin+1)) / (GenPlot->GetBinContent(bin+1));
	GenSignalSum[bin] = GenPlot->GetBinContent(bin+1);
	cout<<"GenSignalSum[bin]: "<<GenSignalSum[bin]<<endl;
	cout<<"efficiencies[bin]: "<<efficiencies[bin]<<endl;
      }      
    }
    else{
      for (Int_t bin=0; bin<bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
	BGSum[bin]+=varhists[hist]->GetBinContent(bin+1);
      }
    }      
  }
    
  // DAVID 
  if ( doUnfolding == true ) {
    
    // SVD Helper Class
    DilepSVDFunctions mySVDFunctions; 
    mySVDFunctions.SetOutputPath(outpath); 

		
    // Getting the histogram 
    TH1* theDataHist = NULL;
    TH1* theBgrHist = NULL;
    TH1* theTtBgrHist = NULL;
    TH1* theRecHist = NULL;
    TH1* theGenHist = GenPlot; 
    TH1* theRespHist = genReco2d;
    for ( size_t i = 0; i < hists.size() ; i++ ) {
      if ( legends[i] == "Data" ) {
		if ( theDataHist == NULL ) {
		  theDataHist = (TH1*) (varhists[i])->Clone("theDataHist");
		} else {
		  theDataHist->Add(varhists[i]);
		}
      } else if ( legends[i] == "t#bar{t} Signal") {
		if ( theRecHist == NULL ) {
		  theRecHist = (TH1*) (varhists[i])->Clone("theRecHist");
		} else {
		  theRecHist->Add(varhists[i]);
		}
      } else if ( legends[i] == "t#bar{t} Other") {
      	// Attention: The TtBar-Bgr must be added to the regular background
		if ( theTtBgrHist == NULL ) {
		  theTtBgrHist = (TH1*) (varhists[i])->Clone("theRecHist");
		} else {
		  theTtBgrHist->Add(varhists[i]);
		}		
		if ( theBgrHist == NULL ) {
		  theBgrHist = (TH1*) (varhists[i])->Clone("theBgrHist");
		} else {
		  theBgrHist->Add(varhists[i]);
		}
      } else {
		if ( theBgrHist == NULL ) {
		  theBgrHist = (TH1*) (varhists[i])->Clone("theBgrHist");
		} else {
		  theBgrHist->Add(varhists[i]);
		}
      }
    }
	
    // Binning
    double* theBins = Xbins;
    int numberBins = bins;
		
		
    // Names and Labels
    TString channelLabelStr(channelLabel[channelType]);
    TString theChannelName = channelLabelStr; 		
    if ( channelLabelStr.Contains("#mu#mu")  ) theChannelName = "mumu";
    if ( channelLabelStr.Contains("e#mu")    ) theChannelName = "emu";
    if ( channelLabelStr.Contains("ee")      ) theChannelName = "ee";
    if ( channelLabelStr.Contains("Dilepton Combined")    ) theChannelName = "combined";
    TString theParticleName = "";
    if ( name.Contains("Lepton") ) theParticleName = "Leptons";
    if ( name.Contains("LLBar")   ) theParticleName = "LepPair";
    if ( name.Contains("Top")     ) theParticleName = "TopQuarks";
    if ( name.Contains("TTBar")   ) theParticleName = "TtBar";
    if ( name.Contains("BJet")    ) theParticleName = "BJets";
    TString theQuantityName = "";
    if ( name.Contains("pT")      ) theQuantityName = "Pt";
    if ( name.Contains("Eta")     ) theQuantityName = "Eta";
    if ( name.Contains("Rapidity")) theQuantityName = "Rapidity";
    if ( name.Contains("Mass")    ) theQuantityName = "Mass";
    TString theSpecialPostfix = "";
    if ( specialComment.CompareTo("Standard") != 0 ) {
    	//theSpecialPostfix = specialComment;
    } 
    
    
    double totalDataEventsNom[1]  = {TopSVDFunctions::SVD_Integral1D((TH1D*)theDataHist, 0, false)}; 
    double totalBgrEventsNom[1]   = {TopSVDFunctions::SVD_Integral1D((TH1D*)theBgrHist, 0, false)};
    double totalTtBgrEventsNom[1]   = {TopSVDFunctions::SVD_Integral1D((TH1D*)theTtBgrHist, 0, false)};
    double totalRecEventsNom[1]   = {TopSVDFunctions::SVD_Integral1D((TH1D*)theRecHist, 0, false)};
    double totalGenEventsNom[1]  = {TopSVDFunctions::SVD_Integral1D((TH1D*)theGenHist, 0, true)}; 
    
    
    // UNFOLDING 
    // Retrieve a histogram with the unfolded quantities.
    // Note: The unfolded histogram has additional side bins!
    // Keep this in mind when accessing bin content via indices
    TH1D* unfoldedDistribution = NULL;
    TH1D* unfoldedDistributionNormalized = NULL;
    int numSystematics = 0;
    mySVDFunctions.SVD_DoUnfold(
				(TH1D*) theDataHist, 
				(TH1D*) theBgrHist, 
				(TH1D*) theTtBgrHist, 
				(TH1D*) theGenHist, 
				(TH1D*) theRecHist, 
				(TH2D*) theRespHist, 
                totalDataEventsNom, 
                totalBgrEventsNom, 
                totalTtBgrEventsNom,    
                totalRecEventsNom,  
                totalGenEventsNom,   
				theBins, numberBins,  
				unfoldedDistribution, 
				unfoldedDistributionNormalized,
				numSystematics,
				theChannelName, theParticleName, theQuantityName, theSpecialPostfix, "");
 	
 		
		
    // Make a vector from the result
    double UnfoldingResult[XAxisbinCenters.size()];
    double UnfoldingError[XAxisbinCenters.size()];
    for ( size_t i = 0; i < XAxisbinCenters.size() ; i++ ) {
      UnfoldingResult[i] = unfoldedDistributionNormalized->GetBinContent(i+2);//account for extra row in SVD unfolding
      UnfoldingError[i] = unfoldedDistributionNormalized->GetBinError(i+2);
      //UnfoldingResult[i] = unfoldedDistribution->GetBinContent(i+2);//account for extra row in SVD unfolding
      //UnfoldingError[i] = unfoldedDistribution->GetBinError(i+2);
    }
		
		 	
    SignalEventswithWeight=0;
     // CROSS SECTION CALCULATION
    for (Int_t i=0; i<bins; ++i) {
      SignalEventswithWeight+=GenSignalSum[i];
    }

    for (Int_t i=0; i<bins; ++i) {
      if(channelType!=3){
	binWidth[i] = Xbins[i+1]-Xbins[i];       
	DiffXSecVec[channelType][i] = UnfoldingResult[i]/(binWidth[i]);
	DiffXSecStatErrorVec[channelType][i] = UnfoldingError[i]/(binWidth[i]); // statistical error 
	GenDiffXSecVec[channelType][i] = (GenSignalSum[i]*topxsec)/(SignalEventswithWeight*binWidth[i]);//DIRTY (signal*topxsec)/(total events*binwidth)
	cout<<"GenSignalSum[i]: "<<GenSignalSum[i]<<endl;
	cout<<"topxsec: "<<topxsec<<endl;
	cout<<"SignalEventswithWeight: "<<SignalEventswithWeight<<endl;
	cout<<"binWidth: "<<binWidth[i]<<endl;
	cout<<"GenDiffXSecVec[channelType][i]: "<<GenDiffXSecVec[channelType][i]<<endl;
	GenDiffXSecError[channelType][i] = TMath::Sqrt(DataSum[i])/(efficiencies[i]*lumi*binWidth[i]); // statistical error
				
	if(name.Contains("Lepton")||name.Contains("Top")||name.Contains("BJet")){
	  //	  DiffXSecVec[channelType][i]=DiffXSecVec[channelType][i]/2.;
	  DiffXSecVec[channelType][i]=DiffXSecVec[channelType][i];
	  GenDiffXSecVec[channelType][i]=GenDiffXSecVec[channelType][i]/2.;
	  DiffXSecStatErrorVec[channelType][i]=DiffXSecStatErrorVec[channelType][i];
	}
      }else{//For the combination
	binWidth[i] = Xbins[i+1]-Xbins[i];      
	DiffXSecVec[channelType][i] =(DiffXSecVec[0][i]/(DiffXSecStatErrorVec[0][i]*DiffXSecStatErrorVec[0][i])
				   +DiffXSecVec[1][i]/(DiffXSecStatErrorVec[1][i]*DiffXSecStatErrorVec[1][i])			
				   +DiffXSecVec[2][i]/(DiffXSecStatErrorVec[2][i]*DiffXSecStatErrorVec[2][i]))/
	  (1/(DiffXSecStatErrorVec[0][i]*DiffXSecStatErrorVec[0][i])
	   +(1/(DiffXSecStatErrorVec[1][i]*DiffXSecStatErrorVec[1][i]))			
	   +(1/(DiffXSecStatErrorVec[2][i]*DiffXSecStatErrorVec[2][i])));			
			
	DiffXSecStatErrorVec[channelType][i]=1/(TMath::Sqrt((1/(DiffXSecStatErrorVec[0][i]*DiffXSecStatErrorVec[0][i]))
							 +(1/(DiffXSecStatErrorVec[1][i]*DiffXSecStatErrorVec[1][i]))			
							 +(1/(DiffXSecStatErrorVec[2][i]*DiffXSecStatErrorVec[2][i]))));			 
	GenDiffXSecVec[channelType][i] = (GenSignalSum[i]*topxsec)/(SignalEventswithWeight*binWidth[i]);//DIRTY (signal*topxsec)/(total events*binwidth)
	GenDiffXSecError[channelType][i] = TMath::Sqrt(DataSum[i])/(efficiencies[i]*lumi*binWidth[i]); // statistical error
	if(name.Contains("Lepton")||name.Contains("Top")||name.Contains("BJet")){
	  DiffXSecVec[channelType][i]=DiffXSecVec[channelType][i];
	  GenDiffXSecVec[channelType][i]=GenDiffXSecVec[channelType][i]/2.;
	  DiffXSecStatErrorVec[channelType][i]=DiffXSecStatErrorVec[channelType][i];
	}
      }
      //		    h_DiffXSec->SetBinContent(i+1,DiffXSec[channelType][i]);
      //   h_DiffXSec->SetBinError(i+1,DiffXSecStatError[channelType][i]);
      //   h_GenDiffXSec->SetBinContent(i+1,GenDiffXSec[channelType][i]);	
    }	      	
  } else { // BBB Unfolding
    
    for (Int_t i=0; i<bins; ++i) {
      if(channelType!=3){
	binWidth[i] = Xbins[i+1]-Xbins[i];      
	//	cout<<"Datasum[i]: "<<DataSum[i]<<" BGSum[i]: "<<BGSum[i]<<" efficiencies[i]: "<<efficiencies[i]<<" binWidth[i]: "<<binWidth[i]<<" lumi: "<<lumi<<endl; 
	DiffXSecVec[channelType][i] = (DataSum[i]-BGSum[i])/(efficiencies[i]*binWidth[i]*lumi);
	DiffXSecStatErrorVec[channelType][i] = TMath::Sqrt(DataSum[i])/(efficiencies[i]*lumi*binWidth[i]); // statistical error
	//	GenDiffXSec[channelType][i] = (GenSignalSum[i]*topxsec)/(SignalEvents*BranchingFraction[channelType]*binWidth[i]);//DIRTY (signal*topxsec)/(total events*bf*binwidth)
	GenDiffXSecVec[channelType][i] = (GenSignalSum[i]*topxsec)/(SignalEventswithWeight*binWidth[i]);//DIRTY (signal*topxsec)/(total events*binwidth)
	GenDiffXSecError[channelType][i] = TMath::Sqrt(DataSum[i])/(efficiencies[i]*lumi*binWidth[i]); // statistical error
	
	if(name.Contains("Lepton")||name.Contains("Top")||name.Contains("BJet")){
	  DiffXSecVec[channelType][i]=DiffXSecVec[channelType][i]/2.;
	  GenDiffXSecVec[channelType][i]=GenDiffXSecVec[channelType][i]/2.;
	  DiffXSecStatErrorVec[channelType][i]=DiffXSecStatErrorVec[channelType][i]/2.;
	}
	if (efficiencies[i] == 0) { //cannot divide by zero
	  cout << "WARNING in PlotDifferentialCrossSection: Efficieny is zero in bin " << i << " while creating " << name << endl;
	}
	//cout<<endl;
      }else{//For the combination
	binWidth[i] = Xbins[i+1]-Xbins[i];      
	DiffXSecVec[channelType][i] =(DiffXSecVec[0][i]/(DiffXSecStatErrorVec[0][i]*DiffXSecStatErrorVec[0][i])
				      +DiffXSecVec[1][i]/(DiffXSecStatErrorVec[1][i]*DiffXSecStatErrorVec[1][i])			
				      +DiffXSecVec[2][i]/(DiffXSecStatErrorVec[2][i]*DiffXSecStatErrorVec[2][i]))/
	  (1/(DiffXSecStatErrorVec[0][i]*DiffXSecStatErrorVec[0][i])
	   +(1/(DiffXSecStatErrorVec[1][i]*DiffXSecStatErrorVec[1][i]))			
	   +(1/(DiffXSecStatErrorVec[2][i]*DiffXSecStatErrorVec[2][i])));			
	
	DiffXSecStatErrorVec[channelType][i]=1/(TMath::Sqrt((1/(DiffXSecStatErrorVec[0][i]*DiffXSecStatErrorVec[0][i]))
							    +(1/(DiffXSecStatErrorVec[1][i]*DiffXSecStatErrorVec[1][i]))			
							    +(1/(DiffXSecStatErrorVec[2][i]*DiffXSecStatErrorVec[2][i]))));			
	
	//cout<<"&&&&&&&&&&&&&&&!!!!!!!!ee DiffCross Sec: "<<DiffXSec[0][i]<<" +/- "<<DiffXSecStatError[0][i]<<endl;
	//cout<<"&&&&&&&&&&&&&&&!!!!!!!!mumu DiffCross Sec: "<<DiffXSec[1][i]<<" +/- "<<DiffXSecStatError[1][i]<<endl;
	//cout<<"&&&&&&&&&&&&&&&!!!!!!!!emu DiffCross Sec: "<<DiffXSec[2][i]<<" +/- "<<DiffXSecStatError[2][i]<<endl;
	//cout<<"&&&&&&&&&&&&&&&!!!!!!!!Combined DiffCross Sec: "<<DiffXSec[3][i]<<" +/- "<<DiffXSecStatError[3][i]<<endl;
	GenDiffXSecVec[channelType][i] = (GenSignalSum[i]*topxsec)/(SignalEventswithWeight*binWidth[i]);//DIRTY (signal*topxsec)/(total events*binwidth)
	GenDiffXSecError[channelType][i] = TMath::Sqrt(DataSum[i])/(efficiencies[i]*lumi*binWidth[i]); // statistical error
	if(name.Contains("Lepton")||name.Contains("Top")||name.Contains("BJet")){
	  DiffXSecVec[channelType][i]=DiffXSecVec[channelType][i]/2.;
	  GenDiffXSecVec[channelType][i]=GenDiffXSecVec[channelType][i]/2.;
	  DiffXSecStatErrorVec[channelType][i]=DiffXSecStatErrorVec[channelType][i]/2.;
	}
      }
    }
  }
}
void Plotter::PlotDiffXSec(){
    TH1::AddDirectory(kFALSE); 
    TGaxis::SetMaxDigits(2);
    doSystematics = false;
    if(doSystematics){//############### Syst ################
      cout << endl;
      cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
      cout << "Starting Calculation of Differential Systematics for '" << name << "' in Channel '" << channel << "':" << endl;  
      cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
      CalcDiffSystematics("JES", "JESUP", "JESDOWN", 0);
      CalcDiffSystematics("RES", "RESUP", "RESDOWN", 1);
      CalcDiffSystematics("PU_", "PU_UP", "PU_DOWN", 2);
      CalcDiffSystematics("SCALE", "SCALEUP", "SCALEDOWN", 3);
      CalcDiffSystematics("MATCH", "MATCHUP", "MATCHDOWN", 4);
      CalcDiffSystematics("MASS", "MASSUP", "MASSDOWN", 5);
      CalcDiffSystematics("DY_", "DY_UP", "DY_DOWN", 6);
      CalcDiffSystematics("BG_", "BG_UP", "BG_DOWN", 7);  
      CalcDiffSystematics("HAD", "MCATNLO", "POWHEG", 8); 
      CalcDiffSystematics("BTAG_PT_", "BTAG_PT_UP", "BTAG_PT_DOWN", 9); 
      CalcDiffSystematics("BTAG_ETA_", "BTAG_ETA_UP", "BTAG_ETA_DOWN", 10); 
      DiffFlatSystematics(11,bins); 
      cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
      cout << "Finished Calculation of Differential Systematics for '" << name << "' in Channel '" << channel << "':" << endl;  
      cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
      cout << endl; 
    }
    double topxsec = TOPXSEC;
    //double BranchingFraction[4]={0.0167, 0.0162, 0.0328, 0.06569};//[ee, mumu, emu]
    double SignalEvents = 63244696.0;
    double Xbins[XAxisbins.size()];
    for(unsigned int i = 0; i<XAxisbins.size();i++){Xbins[i]=XAxisbins[i];}
    double binCenters[XAxisbinCenters.size()];
    for(unsigned int i = 0; i<XAxisbinCenters.size();i++){binCenters[i]=XAxisbinCenters[i];}

    TH1 *RecoPlot = new TH1D;
    TH1 *RecoPlotFineBins = new TH1D;
    TH1 *GenPlot =new TH1D;
    TH1 *GenPlotTheory =new TH1D;
    //double DataSum[nbins]={0};
    //double GenSignalSum[nbins]={0};
    //double BGSum[nbins]={0};
    double DataSum[XAxisbinCenters.size()];
    double GenSignalSum[XAxisbinCenters.size()];
    double BGSum[XAxisbinCenters.size()];
    bool init = false;
    TH1 *varhists[hists.size()];
    TH2 *genReco2d=0; 
    TString newname = name;
    if(name.Contains("Hyp")){//Histogram naming convention has to be smarter
      newname.ReplaceAll("Hyp",3,"",0);
    }
    for (unsigned int i =0; i<hists.size(); i++){
      varhists[i]=hists[i].Rebin(bins,"varhists",Xbins);  
      setStyle(*varhists[i], i);
      if(legends[i] == "t#bar{t} Signal"){
	TFile *ftemp = TFile::Open(dataset[i]);
	if(init==false){
           std::cout << ftemp->Get("Reco"+newname) << " from " << dataset[i] << " Reco"+newname+"\n";
	  RecoPlotFineBins =  (TH1D*)ftemp->Get("Reco"+newname)->Clone();
	  genReco2d = (TH2*)ftemp->Get("GenReco"+newname)->Clone();
	  GenPlotTheory=(TH1D*)ftemp->Get("VisGen"+newname)->Clone();
	  
	  if(newname.Contains("Lepton")||newname.Contains("Top")||newname.Contains("BJet")){
	    RecoPlotFineBins->Add((TH1D*)ftemp->Get("RecoAnti"+newname)->Clone());
	    genReco2d->Add((TH2*)ftemp->Get("GenRecoAnti"+newname)->Clone());
	    GenPlotTheory->Add((TH1D*)ftemp->Get("VisGenAnti"+newname)->Clone());
	  }
	  init =true;
	} else {//account for more than one signal histogram
	  RecoPlotFineBins->Add((TH1D*)ftemp->Get("Reco"+newname)->Clone());
	  genReco2d->Add((TH2*)ftemp->Get("GenReco"+newname)->Clone());
	  GenPlotTheory->Add((TH1D*)ftemp->Get("VisGen"+newname)->Clone());
	  if(newname.Contains("Lepton")||newname.Contains("Top")||newname.Contains("BJet")){
	    GenPlotTheory->Add((TH1D*)ftemp->Get("VisGenAnti"+newname)->Clone());
	    genReco2d->Add((TH2*)ftemp->Get("GenRecoAnti"+newname)->Clone());
	    RecoPlotFineBins->Add((TH1D*)ftemp->Get("RecoAnti"+newname)->Clone());
	  }	  
	}
	delete ftemp;
      }           
    }
    GenPlot = GenPlotTheory->Rebin(bins,"genplot",Xbins);	
    RecoPlot = (TH1D*)RecoPlotFineBins->Rebin(bins,"recohists",Xbins);
    THStack * stack=  new THStack("def", "def");
    TLegend *leg = getNewLegendpre();
    int legchange = 0;
    TH1 *varhistsPlotting[hists.size()];
    
    for(unsigned int i=0; i<hists.size() ; i++){ // prepare histos and leg
      setStyle(*varhists[i], i);
      varhistsPlotting[i]=(TH1*)varhists[i]->Clone();
      if(legends[i] != "Data"){
	if((legends[i] == DYEntry) && channelType!=2){
	  varhists[i]->Scale(DYScale[channelType]);
	  varhistsPlotting[i]->Scale(DYScale[channelType]);
	}

	if(i!=(hists.size()-1)){
	  if(legends[i]!=legends[i+1]){
	    //cout<<legends[i]<<endl;
	    varhistsPlotting[i]->SetLineColor(1);
	  }
	}else{
	  varhistsPlotting[i]->SetLineColor(1);
	}
	
	if(legends[i] != legends[i-1]){
	  varhistsPlotting[i]->SetLineColor(1);
	  stack->Add(varhistsPlotting[i]);  
	}
	if(i > 1){
	  if(legends[i] != legends[i-1]){
	    legchange = i;
	    if( (legends[i] == DYEntry) && DYScale[channelType]!= 1){
	      leg->AddEntry(varhistsPlotting[i], legends[i], "f");
	    }else leg->AddEntry(varhistsPlotting[i], legends[i] ,"f");
	  }else{
	    varhistsPlotting[legchange]->Add(varhistsPlotting[i]);	      
	  }
	}	
      }
      else{
	if(i==0) leg->AddEntry(varhistsPlotting[i], legends[i] ,"pe");
	if(i>0){
	  if(legends[i] != legends[i-1]){//check
	    //	    leg->AddEntry(varhists[i], legends[i] ,"pe");
	  }
	}
      }
    }
    ///////////////////////////////////
    //purity and stability plots as taken from CombinedCrossSection... ...
    
    TH1* genHist = (TH1*)GenPlot->Clone();
    TH1* genRecHist = new TH1D("","",bins,Xbins);
    int intbinsX[XAxisbins.size()];
    int intbinsY[XAxisbins.size()];
    // fill the elements of the main diagonal of the 2d hist into binned 1D histogram
    for (unsigned int i=0; i<XAxisbins.size(); ++i) {
        intbinsX[i] = genReco2d->GetXaxis()->FindBin(Xbins[i]+0.001);
        intbinsY[i] = genReco2d->GetYaxis()->FindBin(Xbins[i]+0.001);
//         std::cout << "PSE_bins " << intbinsX[i] << " " << intbinsY[i] << "\n";
        
//         if (intbinsX[i] != intbinsY[i]) std::cout << "ARGH!!!!!\n";
        if (i>0) {
            genRecHist->SetBinContent(i,
                        ((TH2D*)genReco2d)->Integral( intbinsX[i-1],intbinsX[i]-1,intbinsY[i-1],intbinsY[i]-1));
        }
    }

    TH1* genPseHist = ((TH2D*)genReco2d)->ProjectionY();
    TH1* recPseHist = ((TH2D*)genReco2d)->ProjectionX();
    
    TH1* genBinHist    = genPseHist->Rebin(bins,"genBinHist", Xbins);
    TH1* recBinHist    = recPseHist->Rebin(bins,"recBinHist", Xbins);

    genRecHist->SetBinContent(0,      0);
    genRecHist->SetBinContent(bins+1,0);
    genBinHist->SetBinContent(0,      0);
    genBinHist->SetBinContent(bins+1,0);
    recBinHist->SetBinContent(0,      0);
    recBinHist->SetBinContent(bins+1,0);
    genHist   ->SetBinContent(0,      0);
    genHist   ->SetBinContent(bins+1,0);

    // this is realy ugly but necessary:
    // As it seems, somewhere a double is tranformed into a float so that
    // efficiencies can be larger than 1.
    for(Int_t i=1; i<=genRecHist->GetNbinsX(); ++i){
      if(genRecHist->GetBinContent(i) > recBinHist->GetBinContent(i)){
        genRecHist->SetBinContent(i,recBinHist->GetBinContent(i));
        cout << "WARNING in PlotDifferentialCrossSections: number of events generated and reconstructed in bin" << i
	    << " = " << genRecHist->GetBinContent(i) << " is larger than number of reconstructed events in that bin"
	    << " = " << recBinHist->GetBinContent(i) << endl;
      }
      if(genRecHist->GetBinContent(i) > genBinHist->GetBinContent(i)){
        genRecHist->SetBinContent(i,genBinHist->GetBinContent(i));
	    cout << "WARNING in PlotDifferentialCrossSections: number of events generated and reconstructed in bin " << i
	    << " is larger than number of genrated events in that bin" << endl;
      }
    }

    // efficiency, purity, stability
    TGraphAsymmErrors* grE; // for efficiency
    TGraphAsymmErrors* grP; // for purity
    TGraphAsymmErrors* grS; // for stability

    // efficiency
    grE = new TGraphAsymmErrors(recBinHist, genHist);
    grE->SetMinimum(0);
    grE->SetMaximum(1);
    grE->SetLineColor(8);
    grE->SetLineWidth(2);
    grE->SetMarkerSize(2);
    grE->SetMarkerStyle(21);
    grE->SetMarkerColor(8);

    // purity
    grP = new TGraphAsymmErrors(genRecHist, recBinHist);
    grP->SetLineColor(4);
    grP->SetLineWidth(2);
    grP->SetMarkerSize(2);
    grP->SetMarkerStyle(23);
    grP->SetMarkerColor(4);

    // stability
    grS = new TGraphAsymmErrors(genRecHist, genBinHist);
    grS->SetLineColor(2);
    grS->SetLineWidth(2);
    grS->SetMarkerSize(2);
    grS->SetMarkerStyle(22);
    grS->SetMarkerColor(2);


    grE->GetXaxis()->SetTitle(XAxis);
    TCanvas * cESP = new TCanvas("ESP","ESP");

    // this is a dummy to get the x axis range corrct
    recBinHist->Reset();
    recBinHist->Draw();
    recBinHist->SetMaximum(1.);
    //    FormatHisto(recBinHist);
    recBinHist->GetXaxis()->SetTitle(TString("Reconstructed ").Copy().Append(XAxis));
    recBinHist->GetXaxis()->SetNoExponent(kTRUE);
    grE->GetXaxis()->SetNoExponent(kTRUE);
    grE->Draw("P,SAME");
    grP->Draw("P,SAME");
    grS->Draw("P,SAME");
    TLegend* leg3 = getNewLegend(); // new TLegend(0.60,0.73,0.95,0.83);
    leg3->SetFillStyle(0);
    leg3->SetBorderSize(0);
    leg3->AddEntry(grE, "Efficiency", "p" );
    leg3->AddEntry(grP, "Purity",    "p" );
    leg3->AddEntry(grS, "Stability", "p" );
    leg3->Draw("SAME");


    cESP->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/ESP_"+name+".eps");
    cESP->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/ESP_"+name+".C");
    cESP->Clear();
    delete cESP;
    double efficiencies[XAxisbinCenters.size()];

    init = false;
    for (unsigned int hist =0; hist<hists.size(); hist++){
      if(legends[hist] == "Data"){
	    for (Int_t bin=0; bin<bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
	      DataSum[bin]+=varhists[hist]->GetBinContent(bin+1);
	    }
      }
      else if((legends[hist] == "t#bar{t} Signal")&&init==false){
	signalHist=hist;
	init=true;
	    for (Int_t bin=0; bin<bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
	      efficiencies[bin] = (RecoPlot->GetBinContent(bin+1)) / (GenPlot->GetBinContent(bin+1));
	      GenSignalSum[bin] += GenPlot->GetBinContent(bin+1);
	      //	  cout<<"efficiencies[bin]: "<<efficiencies[bin]<<endl;
	    }      
      }
      else{
		for (Int_t bin=0; bin<bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
		  BGSum[bin]+=varhists[hist]->GetBinContent(bin+1);
		}
      }      
    }


    double totalDataSum = 0;
    for (Int_t bin=0; bin<bins; ++bin) {
      totalDataSum+=DataSum[bin];
    }
    TH1 *h_DiffXSec = (TH1D*)varhists[0]->Clone();
    TH1 *h_GenDiffXSec = (TH1D*)varhists[0]->Clone();
    h_DiffXSec->Reset();
    h_GenDiffXSec->Reset();

    CalcDiffXSec(varhists, RecoPlot, GenPlot, genReco2d, DiffXSec,DiffXSecStatError,GenDiffXSec);

    for (Int_t i=0; i<bins; ++i) {
      cout<<GenDiffXSec[channelType][i]<<endl;
      h_DiffXSec->SetBinContent(i+1,DiffXSec[channelType][i]);
      h_DiffXSec->SetBinError(i+1,DiffXSecStatError[channelType][i]);
      h_GenDiffXSec->SetBinContent(i+1,GenDiffXSec[channelType][i]);	
    }

    //data normalization
    double datascale;
    datascale = h_DiffXSec->Integral("width");
    //    if(name.Contains("LeptonEta")){
    //TotalVisXSection[channelType] = datascale;

      TotalVisXSection[channelType] = 1.;
      //}

    //cout<<"VISIBLE CROSS-SECTION: "<<TotalVisXSection[channelType]<<endl;
    h_DiffXSec->Scale(1/TotalVisXSection[channelType]);
    //theory curves

    double DiffXSecPlot[XAxisbinCenters.size()];
    double DiffXSecStatErrorPlot[XAxisbinCenters.size()];
    double DiffXSecTotalErrorPlot[XAxisbinCenters.size()];

    double ModelSysPlot[XAxisbinCenters.size()];
    double ExpSysPlot[XAxisbinCenters.size()];;

    for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
      double syst_square = 0;
      ExpSysPlot[bin]=0.;
      ModelSysPlot[bin]=0.;
      if(doSystematics){
	for(int syst =0; syst<15; syst++){ 
	  syst_square += DiffXSecSysErrorBySyst[channelType][bin][syst]*DiffXSecSysErrorBySyst[channelType][bin][syst];
	  if(legendsSyst[syst]=="RES" ||legendsSyst[syst]=="JES" ||legendsSyst[syst]=="PU_" ||legendsSyst[syst]=="DY_" ||legendsSyst[syst]=="BG_" ||legendsSyst[syst]=="trigger" ||legendsSyst[syst]=="lepton" ||legendsSyst[syst]=="BTAG" ||legendsSyst[syst]=="kin fit"){
	    ExpSysPlot[bin]+=DiffXSecSysErrorBySyst[channelType][bin][syst]*DiffXSecSysErrorBySyst[channelType][bin][syst];
	  }
	  else{
	    ModelSysPlot[bin]+=DiffXSecSysErrorBySyst[channelType][bin][syst]*DiffXSecSysErrorBySyst[channelType][bin][syst];
	  }
	}
      }
      ExpSysPlot[bin]=sqrt(ExpSysPlot[bin]);
      ModelSysPlot[bin]=sqrt(ModelSysPlot[bin]);
      DiffXSecSysError[channelType][bin] = sqrt(syst_square)*DiffXSec[channelType][bin];
      DiffXSecStatError[channelType][bin] = DiffXSecStatError[channelType][bin]/TotalVisXSection[channelType];
      DiffXSecTotalError[channelType][bin] = sqrt(DiffXSecSysError[channelType][bin]*DiffXSecSysError[channelType][bin] + DiffXSecStatError[channelType][bin]*DiffXSecStatError[channelType][bin]);
      DiffXSecPlot[bin]=DiffXSec[channelType][bin]/TotalVisXSection[channelType];
      DiffXSecStatErrorPlot[bin]=DiffXSecStatError[channelType][bin];//TotalVisXSection[channelType];
      DiffXSecTotalErrorPlot[bin]=DiffXSecTotalError[channelType][bin];//TotalVisXSection[channelType];
    } 

    //create a file for Results!!
    ofstream ResultsFile, ResultsLateX;  
    string ResultsFilestring = outpathPlots.Data();
    ResultsFilestring.append(subfolderChannel.Data());
    ResultsFilestring.append(subfolderSpecial.Data());   
    ResultsFilestring.append("/"); 
    ResultsFilestring.append(newname); 
    ResultsFilestring.append("Results.txt");
    ResultsFile.open(ResultsFilestring.c_str());
    
    
    string ResultsFilestringLatex = outpathPlots.Data();
    ResultsFilestringLatex.append(subfolderChannel.Data());
    ResultsFilestringLatex.append(subfolderSpecial.Data()); 
    ResultsFilestringLatex.append("/"); 
    ResultsFilestringLatex.append(newname); 
    ResultsFilestringLatex.append("ResultsLaTeX.txt");
    ResultsLateX.open(ResultsFilestringLatex.c_str());
    ResultsLateX<<"Bin Center & Bin & 1/#sigma d#sigma/dX & stat(\%) & syst(\%) & total(\%)"<<endl;
    for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
      ResultsFile<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[bin]<<" bin: "<<Xbins[bin]<<" to "<<Xbins[bin+1]<<" DiffXsec: "<<DiffXSecPlot[bin]<<" StatError(percent): "<<DiffXSecStatErrorPlot[bin]/DiffXSecPlot[bin]<<" SysError: "<<DiffXSecSysError[channelType][bin]/DiffXSecPlot[bin]<<" TotalError: "<<DiffXSecTotalErrorPlot[bin]/DiffXSecPlot[bin]<<endl;
      ResultsLateX<<"$"<<h_DiffXSec->GetBinCenter(bin+1)<<"$ & $" <<h_DiffXSec->GetBinLowEdge(bin+1)<<"$ to $"<<h_DiffXSec->GetBinLowEdge(bin+2)<<"$ & ";
      ResultsLateX<<DiffXSecPlot[bin]<<" & "<<setprecision(3)<<DiffXSecStatErrorPlot[bin]*100./DiffXSecPlot[bin]<<" & "<<setprecision(3)<<100.*DiffXSecSysError[channelType][bin]/DiffXSecPlot[bin]<<" & "<<setprecision(3)<<100.*DiffXSecTotalErrorPlot[bin]/DiffXSecPlot[bin]<< "\\\\" << endl;
    }
    ResultsFile.close();
    ResultsLateX.close();
    
    if(doSystematics){
      //The Markus plots
      TCanvas * c10 = new TCanvas("Markus","Markus");
      THStack* SystHists = new THStack("MSTACK","MSTACK");
      TLegend * leg10 =  new TLegend(0.20,0.65,0.45,0.90);
      std::map<int, int> FillOrder;
      FillOrder[14] = 0;   //JES
      FillOrder[13] = 1;   //RES
      FillOrder[12] = 2;   //PU
      FillOrder[11] = 6;   //DY
      FillOrder[10] = 7;   //BG
      FillOrder[9] = 11;   //Trigg
      FillOrder[8] = 12;   //Lep
      FillOrder[7] = 10;  //Btag
      FillOrder[6] = 9;  //Btag
      FillOrder[5] = 13;  //KinFit
      FillOrder[4] = 3;   //SCALE
      FillOrder[3] = 5;  //MASS
      FillOrder[2] = 4;  //MATCH
      FillOrder[1] = 8; //HAD
      FillOrder[0] = 14; //PDF
      
      ofstream ResultsSysFilestring; 
      string ResultsSystLaTeX = outpathPlots.Data();
      ResultsSystLaTeX.append(subfolderChannel.Data());
      ResultsSystLaTeX.append(subfolderSpecial.Data());    
      ResultsSystLaTeX.append("/"); 
      ResultsSystLaTeX.append(newname); 
      ResultsSystLaTeX.append("SystematicsLaTeX.txt");
      FILE *systfile;
      systfile = fopen(ResultsSystLaTeX.c_str(), "w");
      
      for(int systs =0; systs<15; systs++){
	int syst = FillOrder[systs];
	if (syst==10) {continue;}//Skip the BTAG_ETA systematic because it's added in quadrature to BTAG_PT
	TH1D* systtemp = (TH1D*)varhists[0]->Clone();
	systtemp->Reset();
	double TotalSyst=0.0, TotalSqSyst=0.0;
	double AvgSyst= 0.0, SqAvgSys=0.0;
	
	for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
	  if(syst==9){//sum up in quadrature the pT and Eta variations of btagging
	    DiffXSecSysErrorBySyst[channelType][bin][syst]= TMath::Sqrt((DiffXSecSysErrorBySyst[channelType][bin][syst]*DiffXSecSysErrorBySyst[channelType][bin][syst])+(DiffXSecSysErrorBySyst[channelType][bin][syst+1]*DiffXSecSysErrorBySyst[channelType][bin][syst+1]));
	  }
	  //systtemp->SetBinContent(bin+1,(DiffXSecSysErrorBySyst[channelType][bin][syst]/DiffXSec[channelType][bin])*(DiffXSecSysErrorBySyst[channelType][bin][syst]/DiffXSec[channelType][bin]));
	  systtemp->SetBinContent(bin+1,(DiffXSecSysErrorBySyst[channelType][bin][syst]*DiffXSecSysErrorBySyst[channelType][bin][syst]));
	  if(bin==0){  
	    if(syst==0) fprintf(systfile, "JES    ");
	    if(syst==1) fprintf(systfile, "RES    ");
	    if(syst==2) fprintf(systfile, "PU     ");
	    if(syst==3) fprintf(systfile, "SCALE  ");
	    if(syst==4) fprintf(systfile, "MATCH  ");
	    if(syst==5) fprintf(systfile, "MASS   ");
	    if(syst==6) fprintf(systfile, "DY     ");
	    if(syst==7) fprintf(systfile, "BG     ");
	    if(syst==8) fprintf(systfile, "HAD    ");
	    if(syst==9) fprintf(systfile, "Btag   ");
	    //ssyt==10 BTAG_ETA, it summed in quadrature with BTAG_PT syst==9
	    if(syst==11) fprintf(systfile, "Trigg  ");
	    if(syst==12) fprintf(systfile, "Lep.   ");
	    if(syst==13) fprintf(systfile, "KinFit ");
	    if(syst==14) fprintf(systfile, "PDF    ");
	  }
	  fprintf(systfile, "%2.5f ", TMath::Sqrt(systtemp->GetBinContent(bin+1))*100);
	  if(bin>0 && bin<bins-1){//Exclude the 2 side bins
	    TotalSyst=TotalSyst+TMath::Sqrt(systtemp->GetBinContent(bin+1));
	    TotalSqSyst=TotalSqSyst+systtemp->GetBinContent(bin+1);
	  }
	}
	AvgSyst=TotalSyst/(bins-2);
	SqAvgSys=TMath::Sqrt(TotalSqSyst/(bins-2));
	fprintf(systfile, "Lin.Avg.(%)= %.5f  Quad.Avg.(%)=%.5f\n", 100*AvgSyst, 100*SqAvgSys);
	systtemp->SetFillColor(15-systs);
	SystHists->Add((TH1D*)systtemp->Clone());
	leg10->AddEntry(systtemp->Clone(), legendsSyst[syst], "f");
	delete systtemp;
      }
      SystHists->Draw();
      fclose(systfile);
      
      if(name.Contains("pT") ||name.Contains("Mass") ){
	SystHists->GetHistogram()->GetXaxis()->SetTitle(XAxis.Copy().Append(" #left[GeV#right]"));
	if(name.Contains("Rapidity")) SystHists->GetHistogram()->GetXaxis()->SetTitle(XAxis);
      }
      else  SystHists->GetHistogram()->GetXaxis()->SetTitle(XAxis);
      SystHists->GetHistogram()->GetYaxis()->SetTitle("#sum #left( #frac{#Delta #sigma}{#sigma} #right)^{2}");
      SystHists->GetXaxis()->SetNoExponent(kTRUE);
      
      
      leg10->SetFillColor(0);
      leg10->Draw("SAME");
      c10->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/MSP_"+name+".eps");
      c10->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/MSP_"+name+".C");
      c10->Clear();
      delete leg10;
      delete c10;
      
      //The Experimental/Model/Statistical plot
      TCanvas * c11 = new TCanvas("EMS","EMS");
      TH1D* ExpHist = (TH1D*)varhists[0]->Clone();
      TH1D* ModelHist = (TH1D*)varhists[0]->Clone();
      TH1D* StatHist = (TH1D*)varhists[0]->Clone();
      TH1D* TotalHist = (TH1D*)varhists[0]->Clone();
      TLegend * leg11 =  new TLegend(0.65,0.60,0.90,0.85);
      ExpHist->Reset();ModelHist->Reset();StatHist->Reset();TotalHist->Reset();
      for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
	ExpHist->SetBinContent(bin+1,100*ExpSysPlot[bin]);
	ModelHist->SetBinContent(bin+1,100*ModelSysPlot[bin]);
	StatHist->SetBinContent(bin+1,100*DiffXSecStatError[channelType][bin]/DiffXSec[channelType][bin]);
	TotalHist->SetBinContent(bin+1,100*DiffXSecTotalError[channelType][bin]/DiffXSec[channelType][bin]);
      }
      TotalHist->SetMinimum(0.);
      TotalHist->GetYaxis()->SetTitle("#frac{#Delta#sigma}{#sigma} [%]");
      TotalHist->SetLineColor(1);
      ExpHist->SetLineColor(kRed);
      StatHist->SetLineColor(kGreen);
      ModelHist->SetLineColor(kBlue);
      leg11->SetFillColor(0);
      leg11->AddEntry(ExpHist->Clone(), "Experimental Uncertainty", "l");
      leg11->AddEntry(StatHist->Clone(), "Statistical Uncertainty", "l");
      leg11->AddEntry(ModelHist->Clone(), "Model Uncertainty", "l");
      leg11->AddEntry(TotalHist->Clone(), "Total Uncertainty", "l");
      TotalHist->Draw();ModelHist->Draw("SAME");ExpHist->Draw("SAME");StatHist->Draw("SAME");
      leg11->Draw("SAME");
      TotalHist->GetXaxis()->SetNoExponent(kTRUE);
      c11->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/SEM_"+name+".eps");
      c11->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/SEM_"+name+".C");
      c11->Clear();
      delete ExpHist;delete StatHist;delete ModelHist;delete TotalHist;
      delete leg11;
      delete c11;
    }
    Double_t mexl[XAxisbinCenters.size()];
    Double_t mexh[XAxisbinCenters.size()];
    for (unsigned int j=0; j<XAxisbinCenters.size();j++){mexl[j]=0;mexh[j]=0;}
    TGraphAsymmErrors *tga_DiffXSecPlot = new TGraphAsymmErrors(bins, binCenters, DiffXSecPlot, mexl, mexh, DiffXSecStatErrorPlot, DiffXSecStatErrorPlot);
    tga_DiffXSecPlot->SetMarkerStyle(1);
    tga_DiffXSecPlot->SetMarkerColor(kBlack);
    tga_DiffXSecPlot->SetMarkerSize(1);
    tga_DiffXSecPlot->SetLineColor(kBlack);
   
    TGraphAsymmErrors *tga_DiffXSecPlotwithSys = new TGraphAsymmErrors(bins, binCenters, DiffXSecPlot, mexl, mexh, DiffXSecTotalErrorPlot, DiffXSecTotalErrorPlot);
    tga_DiffXSecPlotwithSys->SetMarkerStyle(20);
    tga_DiffXSecPlotwithSys->SetMarkerColor(kBlack);
    tga_DiffXSecPlotwithSys->SetMarkerSize(1);
    tga_DiffXSecPlotwithSys->SetLineColor(kBlack);

    //    for( Int_t j = 0; j< GenPlotTheory->GetNbinsX(); j++ ){
    // GenPlotTheory->SetBinContent(j,(GenPlotTheory->GetBinContent(j)*topxsec)/(SignalEvents*GenPlotTheory->GetBinWidth(j)));
    //}
    GenPlotTheory->Scale(topxsec/(SignalEventswithWeight*GenPlotTheory->GetBinWidth(1)));
    GenPlot->Scale(topxsec/(SignalEventswithWeight*GenPlot->GetBinWidth(1)));
    if(name.Contains("Lepton")||name.Contains("Top")||name.Contains("BJet")){
      GenPlotTheory->Scale(1./2.);
      //h_GenDiffXSec->Scale(1./2.);
    }
    double genscale = 1./GenPlotTheory->Integral("width");
    //double genscale = 1.;
    
    //    genscale = 1./ h_GenDiffXSec->Integral("width");
    GenPlotTheory->Scale(genscale);
    h_GenDiffXSec->Scale(genscale);

    bool binned_theory=true; //############
  
    
  
    TH1* mcnlohist=0, *mcnlohistup=0, *mcnlohistdown=0, *powheghist=0;
    TH1* mcnlohistnorm=0;
    TGraph *mcatnloBand=0;
    
    TH1* mcnlohistnormBinned, *mcnlohistupBinned, *mcnlohistdownBinned, *mcnlohistBinned;
    TH1* powheghistBinned;
    
    TH1F *Kidoth1(nullptr);
    TH1F *Kidoth1_Binned(nullptr);
    TFile *KidoFile(nullptr);
    
    if (drawNLOCurves) {
        mcnlohist = GetNloCurve(newname,"MCATNLO");
        double mcnloscale = 1./mcnlohist->Integral("width");
        if (binned_theory==false) mcnlohist->Rebin(2);mcnlohist->Scale(0.5); //#####
        mcnlohist->Scale(mcnloscale);


        if(name.Contains("LeptonpT")){mcnlohistnorm = GetNloCurve("Leptons","Pt","MCatNLO");}//temprorary until I change the naming convention in the root file
        else if(name.Contains("LeptonEta")){mcnlohistnorm = GetNloCurve("Leptons","Eta","MCatNLO");}
        else if(name.Contains("LLBarpT")){mcnlohistnorm = GetNloCurve("LepPair","Pt","MCatNLO");}
        else if(name.Contains("LLBarMass")){mcnlohistnorm = GetNloCurve("LepPair","Mass","MCatNLO");}
        else if(name.Contains("ToppT")){mcnlohistnorm = GetNloCurve("TopQuarks","Pt","MCatNLO");}
        else if(name.Contains("TopRapidity")){mcnlohistnorm = GetNloCurve("TopQuarks","Rapidity","MCatNLO");}
        else if(name.Contains("TTBarpT")){mcnlohistnorm = GetNloCurve("TtBar","Pt","MCatNLO");}
        else if(name.Contains("TTBarRapidity")){mcnlohistnorm = GetNloCurve("TtBar","Rapidity","MCatNLO");}
        else if(name.Contains("TTBarMass")){mcnlohistnorm = GetNloCurve("TtBar","Mass","MCatNLO");}
        else if(name.Contains("BJetpT")){mcnlohistnorm = GetNloCurve("Jets","Pt","MCatNLO");}
        else if(name.Contains("BJetEta")){mcnlohistnorm = GetNloCurve("Jets","Eta","MCatNLO");}
        else {mcnlohistnorm = new TH1();}
        //    if (binned_theory==false) mcnlohistnorm->Rebin(5);mcnlohistnorm->Scale(0.2);
        mcnlohistnormBinned    = mcnlohistnorm->Rebin(bins,"genBinHist", Xbins);

        if(name.Contains("LeptonpT")){mcnlohistup = GetNloCurve("Leptons","Pt","MCNLOup");}//temprorary until I change the naming convention in the root file
        else if(name.Contains("LeptonEta")){mcnlohistup = GetNloCurve("Leptons","Eta","MCNLOup");}
        else if(name.Contains("LLBarpT")){mcnlohistup = GetNloCurve("LepPair","Pt","MCNLOup");}
        else if(name.Contains("LLBarMass")){mcnlohistup = GetNloCurve("LepPair","Mass","MCNLOup");}
        else if(name.Contains("ToppT")){mcnlohistup = GetNloCurve("TopQuarks","Pt","MCNLOup");}
        else if(name.Contains("TopRapidity")){mcnlohistup = GetNloCurve("TopQuarks","Rapidity","MCNLOup");}
        else if(name.Contains("TTBarpT")){mcnlohistup = GetNloCurve("TtBar","Pt","MCNLOup");}
        else if(name.Contains("TTBarRapidity")){mcnlohistup = GetNloCurve("TtBar","Rapidity","MCNLOup");}
        else if(name.Contains("TTBarMass")){mcnlohistup = GetNloCurve("TtBar","Mass","MCNLOup");}
        else if(name.Contains("BJetpT")){mcnlohistup = GetNloCurve("Jets","Pt","MCNLOup");}
        else if(name.Contains("BJetEta")){mcnlohistup = GetNloCurve("Jets","Eta","MCNLOup");}
        else {mcnlohistup = new TH1();}
        //    if (binned_theory==false) mcnlohistup->Rebin(5);mcnlohistup->Scale(0.2);
        mcnlohistupBinned    = mcnlohistup->Rebin(bins,"genBinHist", Xbins);

        //    mcnlohistupBinned->Scale(1./mcnlohistupBinned->Integral("width"));

        if(name.Contains("LeptonpT")){mcnlohistdown = GetNloCurve("Leptons","Pt","MCNLOdown");}//temprorary until I change the naming convention in the root file
        else if(name.Contains("LeptonEta")){mcnlohistdown = GetNloCurve("Leptons","Eta","MCNLOdown");}
        else if(name.Contains("LLBarpT")){mcnlohistdown = GetNloCurve("LepPair","Pt","MCNLOdown");}
        else if(name.Contains("LLBarMass")){mcnlohistdown = GetNloCurve("LepPair","Mass","MCNLOdown");}
        else if(name.Contains("ToppT")){mcnlohistdown = GetNloCurve("TopQuarks","Pt","MCNLOdown");}
        else if(name.Contains("TopRapidity")){mcnlohistdown = GetNloCurve("TopQuarks","Rapidity","MCNLOdown");}
        else if(name.Contains("TTBarpT")){mcnlohistdown = GetNloCurve("TtBar","Pt","MCNLOdown");}
        else if(name.Contains("TTBarRapidity")){mcnlohistdown = GetNloCurve("TtBar","Rapidity","MCNLOdown");}
        else if(name.Contains("TTBarMass")){mcnlohistdown = GetNloCurve("TtBar","Mass","MCNLOdown");}
        else if(name.Contains("BJetpT")){mcnlohistdown = GetNloCurve("Jets","Pt","MCNLOdown");}
        else if(name.Contains("BJetEta")){mcnlohistdown = GetNloCurve("Jets","Eta","MCNLOdown");}
        else {mcnlohistdown = new TH1();}
        //    if (binned_theory==false) mcnlohistdown->Rebin(5);mcnlohistdown->Scale(0.2);
        mcnlohistdownBinned    = mcnlohistdown->Rebin(bins,"genBinHist", Xbins);
        //mcnlohistdown->Scale(mcnloscale);
        //    mcnlohistdownBinned->Scale(1./mcnlohistdownBinned->Integral("width"));

        powheghist = GetNloCurve(newname, "POWHEG");
        double powhegscale = 1./powheghist->Integral("width");
        if (binned_theory==false) powheghist->Rebin(2);powheghist->Scale(0.5);
        powheghist->Scale(powhegscale);
            
        powheghistBinned = powheghist->Rebin(bins,"powhegplot",Xbins);	
        for (Int_t bin=0; bin<bins; bin++){
        powheghistBinned->SetBinContent(bin+1,powheghistBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/powheghist->GetBinWidth(1)));
        }
        powheghistBinned->Scale(1./powheghistBinned->Integral("width"));

        mcnlohistBinned = mcnlohist->Rebin(bins,"mcnloplot",Xbins);	
        for (Int_t bin=0; bin<bins; bin++){
        mcnlohistBinned->SetBinContent(bin+1,mcnlohistBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/mcnlohist->GetBinWidth(1)));
        mcnlohistupBinned->SetBinContent(bin+1,mcnlohistupBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/mcnlohistup->GetBinWidth(1)));
        mcnlohistdownBinned->SetBinContent(bin+1,mcnlohistdownBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/mcnlohistdown->GetBinWidth(1)));
        mcnlohistnormBinned->SetBinContent(bin+1,mcnlohistnormBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/mcnlohistnorm->GetBinWidth(1)));
        }
        mcnlohistBinned->Scale(1./mcnlohistBinned->Integral("width"));
        mcnlohistupBinned->Scale(1./mcnlohistnormBinned->Integral("width"));
        mcnlohistdownBinned->Scale(1./mcnlohistnormBinned->Integral("width"));
        mcnlohistnormBinned->Scale(1./mcnlohistnormBinned->Integral("width"));

        for (Int_t bin=0; bin<bins; bin++){
        mcnlohistupBinned->SetBinContent(bin+1,(mcnlohistupBinned->GetBinContent(bin+1)/mcnlohistnormBinned->GetBinContent(bin+1))*mcnlohistBinned->GetBinContent(bin+1));
        mcnlohistdownBinned->SetBinContent(bin+1,(mcnlohistdownBinned->GetBinContent(bin+1)/mcnlohistnormBinned->GetBinContent(bin+1))*mcnlohistBinned->GetBinContent(bin+1));
        }

        //Uncertainty band for MC@NLO
        const Int_t nMCNLOBins = mcnlohistup->GetNbinsX();
        Double_t x[bins];
        Double_t xband[2*bins];
        Double_t errup[bins];
        Double_t errdn[bins];
        Double_t errorband[2*bins];
        
        for( Int_t j = 0; j< bins; j++ ){
        x[j]=mcnlohistBinned->GetBinCenter(j+1);
        errup[j]=(mcnlohistupBinned->GetBinContent(j+1)/mcnlohistnormBinned->GetBinContent(j+1))*mcnlohistBinned->GetBinContent(j+1);
        errdn[j]=(mcnlohistdownBinned->GetBinContent(j+1)/mcnlohistnormBinned->GetBinContent(j+1))*mcnlohistBinned->GetBinContent(j+1);
        
        xband[j] = x[j];
        errorband[j] = errdn[j]; //lower band
        xband[2*bins-j-1] = x[j];
        errorband[2*bins-j-1] = errup[j]; //upper band
        
        }
        
        mcatnloBand = new TGraph(2*bins, xband, errorband);
        mcatnloBand->SetFillColor(kGray);
        mcatnloBand->SetFillStyle(1001);
        mcatnloBand->SetLineColor(kBlue);
        mcatnloBand->SetLineWidth(2);
        mcatnloBand->SetLineStyle(5);
        
        if(name.Contains("ToppT") || name.Contains("TopRapidity")){
        KidoFile=TFile::Open("dilepton_kidonakisNNLO.root");
        if(name.Contains("ToppT")){
            Kidoth1_Binned = (TH1F*)KidoFile->Get("topPt");
        }
        else if(name.Contains("TopRapidity")){
            Kidoth1_Binned = (TH1F*)KidoFile->Get("topY");
        }
        }
        
    //    TH1 *MCFMHist;
    //    TFile* MCFMfile = new TFile("diffCrossSections_normalized_tt_bbl_todk_MSTW200_172_172_ful_central.root","READ");
    //
    //    if(name.Contains("LeptonpT")){MCFMfile->GetObject<TH1>("pt_l", MCFMHist);}
    //    else if(name.Contains("LeptonEta")){MCFMfile->GetObject<TH1>("eta_l", MCFMHist);}
    //    else if(name.Contains("LLBarpT")){MCFMfile->GetObject<TH1>("pt_ll", MCFMHist);}
    //    else if(name.Contains("LLBarMass")){MCFMfile->GetObject<TH1>("m_ll", MCFMHist);}
    //    else if(name.Contains("ToppT")){MCFMfile->GetObject<TH1>("pt_t", MCFMHist);}
    //    else if(name.Contains("TopRapidity")){MCFMfile->GetObject<TH1>("y_t", MCFMHist);}
    //    else if(name.Contains("TTBarpT")){MCFMfile->GetObject<TH1>("pt_tt", MCFMHist);}
    //    else if(name.Contains("TTBarRapidity")){MCFMfile->GetObject<TH1>("y_tt", MCFMHist);}
    //    else if(name.Contains("TTBarMass")){MCFMfile->GetObject<TH1>("m_tt", MCFMHist);}
    //    else{cout<<"probably going to crash soon"<<endl;}
    }
    
    TCanvas * c = new TCanvas("DiffXS","DiffXS");
    
    if(logY){
      c->SetLogy();
    }
    h_DiffXSec->SetMarkerStyle(20);
    //MCFMHist->SetMarkerStyle(2);
    h_GenDiffXSec->SetLineWidth(2);
    h_GenDiffXSec->SetMinimum(ymin);
    if(ymax!=0){
      
      if(logY){  
	h_GenDiffXSec->SetMaximum(18*h_GenDiffXSec->GetBinContent(h_GenDiffXSec->GetMaximumBin()));
      }
      else{ h_GenDiffXSec->SetMaximum(1.5*h_GenDiffXSec->GetBinContent(h_GenDiffXSec->GetMaximumBin()));}
    }
    h_GenDiffXSec->GetXaxis()->SetNoExponent(kTRUE);
    if (name.Contains("Rapidity") || name.Contains("Eta")){h_GenDiffXSec->GetYaxis()->SetNoExponent(kTRUE);}
    h_GenDiffXSec->Draw();
    if (ymax!=0) h_GenDiffXSec->SetMaximum(ymax);
    //    h_DiffXSec->Draw("SAME, EP0");
    gStyle->SetEndErrorSize(8);
    if (drawNLOCurves) {
    //    mcatnloBand->Draw("same, F");
        mcnlohistupBinned->SetFillColor(kGray);
        mcnlohistupBinned->SetLineColor(kGray);
        mcnlohistupBinned->Draw("same");
        mcnlohistdownBinned->SetLineColor(10);
        mcnlohistdownBinned->SetFillColor(10);
        mcnlohistdownBinned->Draw("same");
    }
    GenPlotTheory->SetLineColor(kRed+1);
    GenPlotTheory->SetLineWidth(2);
    GenPlotTheory->SetLineStyle(1);
//    GenPlotTheory->Rebin(2);GenPlotTheory->Scale(1./2.);
    h_GenDiffXSec->SetLineColor(kRed+1);
    h_GenDiffXSec->SetLineStyle(1);

    //bool binned_theory=true; 
    if (drawNLOCurves) {
        mcnlohist->SetLineColor(kBlue); //#####################
        mcnlohist->SetLineStyle(5);
        mcnlohistBinned->SetLineColor(kBlue); //#####################
        mcnlohistBinned->SetLineWidth(2);
        mcnlohistBinned->SetLineStyle(5);
        powheghist->SetLineColor(kGreen+1); //#####################
        powheghist->SetLineStyle(7);
        powheghistBinned->SetLineColor(kGreen+1); //#####################
        powheghistBinned->SetLineWidth(2);
        powheghistBinned->SetLineStyle(7);
        
        if(binned_theory==false){
        mcnlohist->Draw("SAME,C");
        powheghist->Draw("SAME,C");
        }else{
        mcnlohistBinned->Draw("SAME");
        powheghistBinned->Draw("SAME");
        }

        if(name.Contains("ToppT") || name.Contains("TopRapidity")){
        Kidoth1_Binned->SetLineWidth(2);
        Kidoth1_Binned->SetLineColor(kOrange-3); //########################
        Kidoth1_Binned->SetLineStyle(1);
        Kidoth1_Binned->Draw("SAME");
        }
        //MCFMHist->Draw("SAME");
        //h_DiffXSec->Draw("SAME, EP0");
    }

    if(!name.Contains("HypLLBarpT") && !name.Contains("HypTTBarpT") && !name.Contains("HypLeptonpT") && !name.Contains("HypBJetpT")){
        TH1D *SmoothMadgraph =(TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
        SmoothMadgraph->Smooth(10);
        SmoothMadgraph->Draw("SAME, L");
    }
    else if( !name.Contains("HypTTBarpT") && !name.Contains("HypLeptonpT")){
        TH1D *SmoothMadgraph =(TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
        SmoothMadgraph->Smooth(4);
        SmoothMadgraph->Draw("SAME, L");
    }
    else if( !name.Contains("HypTTBarpT")){
        TH1D *SmoothMadgraph =(TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
        SmoothMadgraph->Smooth(2);
        SmoothMadgraph->Draw("SAME, L");
    }
    else {GenPlotTheory->Draw("SAME,C");} //### 150512 ###    

    h_GenDiffXSec->Draw("SAME"); //### 150512 ###

    DrawCMSLabels(false, lumi);
    
    DrawDecayChLabel(channelLabel[channelType]);    
    
    TLegend leg2 = *getNewLegend();
    leg2.AddEntry(h_DiffXSec, "Data",    "p");
    leg2.AddEntry(GenPlotTheory,            "MadGraph","l");
    if (drawNLOCurves) {
        if (mcnlohistup->GetEntries() && mcnlohistdown->GetEntries()) leg2.AddEntry(mcatnloBand,      "MC@NLO",  "fl");
        //else if (mcnlohist->GetEntries()) leg2.AddEntry(mcnlohist,      "MC@NLO",  "l");
        //if (mcnlohistupBinned->GetEntries()) leg2.AddEntry(mcnlohistupBinned,      "MC@NLO",  "fl");
        //if (mcnlohist->GetEntries()) leg2.AddEntry(mcnlohistBinned,      "MC@NLO",  "l");
        if (powheghist->GetEntries())  leg2.AddEntry(powheghistBinned,       "POWHEG",  "l");        
        if (name.Contains("ToppT") || name.Contains("TopRapidity")) leg2.AddEntry(Kidoth1_Binned,       "Approx. NNLO",  "l");
    }
    
    //if (MCFMHist->GetEntries())  leg2.AddEntry(MCFMHist,       "MCFM",  "p");        
    leg2.SetFillStyle(0);
    leg2.SetBorderSize(0);
    leg2.SetX1NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength()-0.30);
    leg2.SetY1NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength()-0.05*(double)leg2.GetNRows());
    leg2.SetX2NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength());
    leg2.SetY2NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength());
    leg2.SetTextSize(0.04);
    leg2.Draw("same");
    if (drawNLOCurves) {
        if (name.Contains("ToppT"))        DrawLabel("(arXiv:1009.4935)", leg2.GetX1NDC()+0.06, leg2.GetY1NDC()-0.025, leg2.GetX2NDC(), leg2.GetY1NDC(), 12, 0.025);
        if (name.Contains("TopRapidity"))  DrawLabel("(arXiv:1105.5167)", leg2.GetX1NDC()+0.06, leg2.GetY1NDC()-0.025, leg2.GetX2NDC(), leg2.GetY1NDC(), 12, 0.025);
    }
    
    h_GenDiffXSec->Draw("SAME");
    gStyle->SetEndErrorSize(10);
    tga_DiffXSecPlot->Draw("p, SAME");
    tga_DiffXSecPlotwithSys->Draw("p, SAME, Z");
    gPad->RedrawAxis();
    
    c->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/DiffXS_"+name+".eps"); 
    c->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/DiffXS_"+name+".C"); 
    c->Clear();
    delete c;
    gStyle->SetEndErrorSize(0);
    
    TCanvas * c1 = new TCanvas("DiffXS","DiffXS");
    TList* l = stack->GetHists();
    TH1D* stacksum = (TH1D*) l->At(0)->Clone();
    for (int i = 1; i < l->GetEntries(); ++i) {
      stacksum->Add((TH1D*)l->At(i));
    } 
    for(unsigned int i=1; i<hists.size() ; i++){ // sum all data plots to first histogram
      if(legends[i] == legends[0]){
	varhists[0]->Add(varhists[i]);
      }
    }
    TH1D* syshist =0;
    syshist = (TH1D*)stacksum->Clone();
    double lumierr = 0.045; 
    //stat uncertainty::make a function 
    for(Int_t i=0; i<=syshist->GetNbinsX(); ++i){
      
      Double_t binc = 0;
      binc += stacksum->GetBinContent(i);
      syshist->SetBinContent(i, binc);
      // calculate uncertainty: lumi uncertainty
      Double_t binerr2 = binc*binc*lumierr*lumierr;
      Double_t topunc = 0; // uncertainty on top xsec
      
      double topxsec = TOPXSEC;
      double topxsecErr2 = 2.2*2.2 + 11.6*11.6;
      
      double topRelUnc =  TMath::Sqrt(topxsecErr2)/topxsec;
      //Functionality for multiple signal histograms
      topunc += varhists[signalHist]->GetBinContent(i)*topRelUnc;
      binerr2 += (topunc*topunc);
      syshist->SetLineColor(1);
      syshist->SetBinError(i, TMath::Sqrt(binerr2));
    }    

    leg = ControlLegend(hists.size(), varhistsPlotting, legends, leg);
    syshist->SetFillStyle(3004);
    syshist->SetFillColor(kBlack);
    leg->AddEntry( syshist, "Uncertainty", "f" );


    varhists[0]->SetMaximum(1.5*varhists[0]->GetBinContent(varhists[0]->GetMaximumBin()));

    varhists[0]->SetMinimum(0);
    varhists[0]->GetYaxis()->SetTitle("events");
    varhists[0]->GetXaxis()->SetNoExponent(kTRUE);
    //varhists[0]->GetYaxis()->SetNoExponent(kTRUE);
    //varhists[0]->Draw("e");  //#########
    varhists[0]->Draw("e"); 
    
    //Add the binwidth to the yaxis in yield plots
    TString ytitle = TString(varhists[0]->GetYaxis()->GetTitle()).Copy();
    double binwidth = varhists[0]->GetXaxis()->GetBinWidth(1);
    std::ostringstream width;
    width<<binwidth;
    if(name.Contains("Rapidity") || name.Contains("Eta")){ytitle.Append(" / ").Append(width.str());}
    else if(name.Contains("pT") || name.Contains("Mass") || name.Contains("mass") || name.Contains("MET") || name.Contains("HT")){ytitle.Append(" / ").Append(width.str()).Append(" GeV");};
    varhists[0]->GetYaxis()->SetTitle(ytitle);

    stack->Draw("same HIST");

    /*    TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0.5)");//this is frustrating and stupid but apparently necessary...
    setex1->Draw();
    syshist->Draw("same,E2");
    TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.)");
    setex2->Draw();*/
    varhists[0]->Draw("same, e1"); //############
    //varhists[0]->Draw("same, e"); 
    DrawCMSLabels(false, lumi);
    DrawDecayChLabel(channelLabel[channelType]);    
    leg->Draw("SAME");
    gPad->RedrawAxis();
    TFile *f1 = new TFile("KinFitPlots.root","UPDATE");
    stacksum->Write(name+"_"+channel+"_MC");
    varhists[0]->Write(name+"_"+channel+"_Data");
    f1->Close();
    c1->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/preunfolded_"+name+".eps");
    c1->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/preunfolded_"+name+".C");
    c1->Clear();
    delete c1; 	



}

// get generator cross section curve for NLO prediction
TH1* Plotter::GetNloCurve(const char *particle, const char *quantity, const char *generator){

  TH1::AddDirectory(kFALSE);
  TString histname;
  if(strcmp(particle, "TopQuarks")==0||strcmp(particle, "TtBar")==0){
    histname="total_";
  }else{
    histname="visible_";
  }
  histname.Append(particle);
  histname.Append(quantity);
  histname.Append("_");
  histname.Append(generator);
  
  TH1* hist;
  
  TFile* file = new TFile;
  
  if(strcmp(generator, "Powheg")==0){file = TFile::Open("selectionRoot/Nominal/emu/ttbarsignalplustau_powheg.root","READ");}
  else if(strcmp(generator, "MCatNLO")==0){file = TFile::Open("MCatNLO_status3_v20120729.root","READ");}
  else if(strcmp(generator, "MCNLOup")==0){file = TFile::Open("MCatNLO_Uncert_Up_status3_v20120729.root","READ");}
  else if(strcmp(generator, "MCNLOdown")==0){file = TFile::Open("MCatNLO_Uncert_Down_status3_v20120729.root","READ");}
  
  if (file && !file->IsZombie()) {
    file->GetObject<TH1>(histname, hist);

    if(!hist){
      std::cerr << "WARNING in GetNloCurve: input histogram '" << histname << "' could not been opened! Returning dummy!" << endl;
      hist = new TH1();
      return hist;
    }
    
    TH1D* rethist = (TH1D*)hist->Clone();
    TH1D* weight = (TH1D*)file->Get(TString("total_LeptonsPt_").Append(generator));
    
    Double_t wgt = 1.;
    if(!weight){
      std::cerr << "WARNING in GetNloCurve: histogram to extract original number of events could not be opened! No weighting applied!" << endl;
    } else{
      Double_t nevents = weight->GetEntries();
      //
      Double_t crosssection = TOPXSEC; //######
      Double_t binw = hist->GetBinWidth(1);
      wgt = crosssection/nevents/binw;
    }
    //rethist->Scale(wgt);
    return rethist;
  }
  
  std::cerr << "WARNING in GetNloCurve: input file could not been opened! Returning dummy!" << endl;
  hist = new TH1D();
  delete file;
  return hist;
}

TH1* Plotter::GetNloCurve(TString NewName, TString Generator){

  TH1::AddDirectory(kFALSE);
  TString histname("VisGen"+NewName);
  
  TH1* hist;
  TH1* hist1;
  TH1* hist2;
  
  TFile* file = new TFile;
  TFile* file1 = new TFile;
  TFile* file2 = new TFile;

    
  if(Generator=="MCATNLO"){
    if(channelType == 0)file = TFile::Open("selectionRoot/"+Generator+"/ee/ttbarsignalplustau_mcatnlo.root","READ");
    else if(channelType == 1)file = TFile::Open("selectionRoot/"+Generator+"/mumu/ttbarsignalplustau_mcatnlo.root","READ");
    else if(channelType == 2)file = TFile::Open("selectionRoot/"+Generator+"/emu/ttbarsignalplustau_mcatnlo.root","READ");
    else {
      file = TFile::Open("selectionRoot/"+Generator+"/emu/ttbarsignalplustau_mcatnlo.root","READ");
      file1 = TFile::Open("selectionRoot/"+Generator+"/ee/ttbarsignalplustau_mcatnlo.root","READ");
      file2 = TFile::Open("selectionRoot/"+Generator+"/mumu/ttbarsignalplustau_mcatnlo.root","READ");
    }
  }else{
    if(channelType == 0)file = TFile::Open("selectionRoot/"+Generator+"/ee/ttbarsignalplustau_powheg.root","READ");
    else if(channelType == 1)file = TFile::Open("selectionRoot/"+Generator+"/mumu/ttbarsignalplustau_powheg.root","READ");
    else if(channelType == 2)file = TFile::Open("selectionRoot/"+Generator+"/emu/ttbarsignalplustau_powheg.root","READ");
    else {
      file = TFile::Open("selectionRoot/"+Generator+"/emu/ttbarsignalplustau_powheg.root","READ");
      file1 = TFile::Open("selectionRoot/"+Generator+"/ee/ttbarsignalplustau_powheg.root","READ");
      file2 = TFile::Open("selectionRoot/"+Generator+"/mumu/ttbarsignalplustau_powheg.root","READ");
    }
  }

  if (file && !file->IsZombie()) {
    if (channelType<3)hist=(TH1*)file->Get("VisGen"+NewName)->Clone();
    else {
      hist=(TH1*)file->Get("VisGen"+NewName)->Clone();
      hist->Add((TH1*)file1->Get("VisGen"+NewName)->Clone());
      hist->Add((TH1*)file2->Get("VisGen"+NewName)->Clone());
    }
    if(NewName.Contains("Lepton")||NewName.Contains("Top")||NewName.Contains("BJet")){
      if(channelType<3)hist->Add((TH1*)file->Get("VisGenAnti"+NewName)->Clone());
      else{
	hist->Add((TH1*)file->Get("VisGenAnti"+NewName)->Clone());
	hist->Add((TH1*)file1->Get("VisGenAnti"+NewName)->Clone());
	hist->Add((TH1*)file2->Get("VisGenAnti"+NewName)->Clone());
      }
    }
    if(!hist){
      std::cerr << "WARNING in GetNloCurve: input histogram '" << histname << "' could not been opened! Returning dummy!" << endl;
      hist = new TH1();
      return hist;
    }
    
    TH1D* rethist = (TH1D*)hist->Clone();
    
    Double_t wgt = 1.;
    Double_t nevents = 16420479;//weight->GetEntries();
    Double_t crosssection = TOPXSEC;
    Double_t binw = hist->GetBinWidth(1);
    wgt = crosssection/nevents/binw;
    rethist->Scale(wgt);
    return rethist;
  }
  
  std::cerr << "WARNING in GetNloCurve: input file could not been opened! Returning dummy!" << endl;
  hist = new TH1D();
  delete file;  delete file1;  delete file2;
  return hist;
}

// get new legend
TLegend* Plotter::getNewLegend() {
  TLegend *leg = new TLegend();
  leg->SetX1NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength()-0.25);
  leg->SetY1NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength()-0.20);
  leg->SetX2NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength());
  leg->SetY2NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength());
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextAlign(12);
  return leg;
}

// get new legend
TLegend* Plotter::getNewLegendpre() {
  TLegend *leg = new TLegend();
  leg->SetX1NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength()-0.25);
  leg->SetY1NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength()-0.30);
  leg->SetX2NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength());
  leg->SetY2NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength());
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextAlign(12);
  return leg;
}

TH1F* Plotter::ConvertGraphToHisto(TGraphErrors *pGraph){
  // takes data from a graph, determines binning and fills data into histogram
  Int_t NPoints = pGraph->GetN();
  Double_t BinLimits[NPoints+1];
  // sort graph
  pGraph->Sort();
  // determine lower limit of histogram: half the distance to next point
  Double_t x0,x1,y;
  pGraph->GetPoint(0,x0,y);
  pGraph->GetPoint(1,x1,y);
  Double_t Distance = TMath::Abs(x0-x1);
  BinLimits[0] = x0 - Distance/2.;
  // now set upper limits for all the other points
  for (Int_t k = 0 ; k<NPoints-1;k++){
    pGraph->GetPoint(k,x0,y);
    pGraph->GetPoint(k+1,x1,y);
    Distance = TMath::Abs(x0-x1);
    BinLimits[k+1] = x0 + Distance/2.;}
  // for the last point set upper limit similar to first point:
  pGraph->GetPoint(NPoints-2,x0,y);
  pGraph->GetPoint(NPoints-1,x1,y);
  Distance = TMath::Abs(x0-x1);
  BinLimits[NPoints] = x1 + Distance/2.;
  // now we know the binning and can create the histogram:
  TString Name = "ConvertedHisto"; 
  // make name unique 
  Name+= rand();
  TH1F *ThisHist = new TH1F(Name,"Converted Histogram",NPoints,BinLimits);
  // now fill the histogram
  for (Int_t i = 0; i<pGraph->GetN();i++){
    Double_t x2,y2;
    pGraph->GetPoint(i,x2,y2);
    ThisHist->SetBinContent(i+1,y2);
  }
  return ThisHist;
}

//TH1F* Plotter::reBinTH1FIrregularNewBinning(TH1F *histoOldBinning, const std::vector<double> &vecBinning, TString plotname, bool rescale=1){
TH1F* Plotter::reBinTH1FIrregularNewBinning(TH1F *histoOldBinning, TString plotname, bool rescale){
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
  cout << endl;
  cout << endl;
  cout << "asdfasdfasdfasdfasdf hallo david " << plotname << " " << rescale << endl;
  cout << "histoOldBinning = ";
  for ( int i = 0 ; i < histoOldBinning->GetXaxis()->GetNbins() + 1; i++ ) cout << " " << histoOldBinning->GetXaxis()->GetBinLowEdge(i+1);
  cout << endl;
  cout << endl;
  cout << endl;
  
   
  unsigned int vecIndex=0;

  // fill vector into array to use appropriate constructor of TH1-classes
  const double *binArray = XAxisbins.data();
	
  // create histo with new binning
  TH1F *histoNewBinning = new TH1F("histoNewBinning"+plotname,"histoNewBinning"+plotname,XAxisbins.size()-1,binArray);
	
  // fill contents of histoOldBinning into histoNewBinning and rescale if applicable
  for (vecIndex = 0; vecIndex < XAxisbins.size()-1; vecIndex++){
	    
    double lowEdge      = XAxisbins[vecIndex]; 
    if (plotname=="topPt"&&vecIndex==0&&lowEdge==0.0) lowEdge+=10;  // adhoc fix to compensate for minimum top-Pt cut in NNLO curve
    double highEdge     = XAxisbins[vecIndex+1];
    double newBinWidth  = highEdge - lowEdge;
    double newBinCenter = 0.5*(highEdge+lowEdge);
    double binSum       = 0.0;	    	  
	    
    for (int j=1; j<histoOldBinning->GetNbinsX(); j++){
		
      double oldBin = histoOldBinning->GetBinCenter(j);
		
      if( (oldBin>=lowEdge) && (oldBin<highEdge) ){		   
	if (rescale) binSum+=histoOldBinning->GetBinContent(j) * histoOldBinning->GetBinWidth(j);
	else         binSum+=histoOldBinning->GetBinContent(j);
      }
    }

    if (rescale) histoNewBinning->Fill(newBinCenter,binSum/newBinWidth);
    else histoNewBinning->Fill(newBinCenter,binSum);
  }

  return (TH1F*)histoNewBinning->Clone();
}

TLegend* Plotter::ControlLegend(int HistsSize, TH1* drawhists[], std::vector<TString> Legends, TLegend *leg){
    //hardcoded ControlPlot legend
    std::vector<TString> OrderedLegends;    
    OrderedLegends.push_back("Data");
    OrderedLegends.push_back("t#bar{t} Signal");
    OrderedLegends.push_back("t#bar{t} Other");
    OrderedLegends.push_back("Single Top");
    OrderedLegends.push_back("W+Jets");
    OrderedLegends.push_back("Z / #gamma* #rightarrow ee/#mu#mu");
    OrderedLegends.push_back("Z / #gamma* #rightarrow #tau#tau");
    OrderedLegends.push_back("Diboson");
    OrderedLegends.push_back("QCD Multijet");
      
    leg->Clear();
    leg->SetX1NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength()-0.25);
    leg->SetY1NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength()-0.30);
    leg->SetX2NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength());
    leg->SetY2NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength());
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextAlign(12);
    for(int i=0; i<(int)OrderedLegends.size(); ++i){
        for(int j=0; j<HistsSize; ++j){
            if (OrderedLegends[i] == Legends[j]){
                if( OrderedLegends[i] == "Data"){
                    leg->AddEntry(drawhists[j], OrderedLegends[i], "pe");
                    break;
                }
                else{
                    leg->AddEntry(drawhists[j], OrderedLegends[i], "f");
                    break;
                }
            }
        }
    }
    return leg;
}

TLegend* Plotter::ControlLegend(int HistsSize, TH1D* drawhists[], std::vector<TString> Legends, TLegend *leg){
    
    //hardcoded ControlPlot legend
    std::vector<TString> OrderedLegends;    
    OrderedLegends.push_back("Data");
    OrderedLegends.push_back("t#bar{t} Signal");
    OrderedLegends.push_back("t#bar{t} Other");
    OrderedLegends.push_back("Single Top");
    OrderedLegends.push_back("W+Jets");
    OrderedLegends.push_back("Z / #gamma* #rightarrow ee/#mu#mu");
    OrderedLegends.push_back("Z / #gamma* #rightarrow #tau#tau");
    OrderedLegends.push_back("Diboson");
    OrderedLegends.push_back("QCD Multijet");
    
    leg->Clear();
    leg->SetX1NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength()-0.25);
    leg->SetY1NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength()-0.30);
    leg->SetX2NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength());
    leg->SetY2NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength());
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextAlign(12);
    for(int i=0; i<(int)OrderedLegends.size(); ++i){
        for(int j=0; j<HistsSize; ++j){
            if (OrderedLegends[i] == Legends[j]){
                if( OrderedLegends[i] == "Data"){
                    leg->AddEntry(drawhists[j], OrderedLegends[i], "pe");
                    break;
                }
                else{
                    leg->AddEntry(drawhists[j], OrderedLegends[i], "f");
                    break;
                }
            }
        }
    }
    return leg;
}

void Plotter::DrawLabel(TString text, const double x1, const double y1, const double x2, const double y2, int centering, double textSize){
    //Function to add Kidonakis references to DiffXSection plots' legends
    TPaveText *label = new TPaveText(x1, y1, x2, y2, "br NDC");
    label->AddText(text);
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    
    if (textSize!=0) label->SetTextSize(textSize);
    label->SetTextAlign(centering);
    label->Draw("same");
}



#endif