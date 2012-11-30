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
#include <set>

// DAVID
#include "DilepSVDFunctions.h"
#include "DilepSVDFunctions.C"

using namespace std;


class Plotter {

 public:
  Plotter();
  Plotter(TString, TString, TString, double, double);
  ~Plotter();
  void   setOptions(TString name_, TString specialComment_, TString YAxis_, TString XAxis_, int rebin_, bool doDYScale_, bool logX_, bool logY_, double ymin_, double ymax_, double rangemin_, double rangemax_, int bins_, std::vector<double> XAxisbins_, std::vector<double> XAxisbinCenters_);
  void   setDataSet(std::vector<TString>, std::vector<double>, std::vector<TString>, std::vector<int>, TString);
  void   setDataSet(TString, TString);
  bool   fillHisto();
  void   setStyle(TH1*, unsigned int, bool = false);
  void   unfolding();
  void   preunfolding(TString Chan="", TString System="");
  void   write(TString, TString);
  void   setLumi(double);
  
  double CalcXSec(std::vector<TString> datasetVec, double InclusiveXsectionVec[4],double InclusiveXsectionStatErrorVec[4], TString Systematic, TString Shift);
  void MakeTable();


  void PlotXSec();
  //  void CalcDiffXSec(TH1* varhists[], TH1* RecoPlot, TH1* GenPlot, TH2* genReco2d, double DiffXSecVec[4][10], double DiffXSecStatErrorVec[4][10]); 
  void CalcDiffXSec(TString, TString);
  void CalcDiffSystematics(TString, TString, TString, TString, double);
  void PlotDiffXSec(TString);

  void DYScaleFactor();


  TLegend* getNewLegend();
  TLegend* getNewLegendpre();

  TH1* GetNloCurve(const char *particle, const char *quantity, const char *generator);
  TH1* GetNloCurve(TString NewName, TString Generator);
  TH1F* ConvertGraphToHisto(TGraphErrors *pGraph);
  TH1F* reBinTH1FIrregularNewBinning(TH1F *histoOldBinning, TString plotname, bool rescale);

  //IVAN's Scaling Code
  void ApplyFlatWeights(TH1* varhists,   const double weight);
  void ApplyFlatWeights(TH1* varhists[], const double weight);
  double SampleXSection(TString filename);
  double CalcLumiWeight(TString WhichSample);
  std::vector<TString> InputFileList(TString mode, TString Systematic);

  // DAVID
  void UnfoldingOptions(bool doSVD);
  void SetOutpath(TString path); 
  void ApplyMCATNLOWeight(TH1* hist, TString Systematic, TString Shift, TString Sample);
  TLegend* ControlLegend(int HistsSize, TH1* drawhists[], std::vector<TString> legends, TLegend *leg);
  TLegend* ControlLegend(int HistsSize, TH1D* drawhists[], std::vector<TString> Legends, TLegend *leg);
  void DrawLabel(TString text, const double x1, const double y1, const double x2, const double y2, int centering, double textSize);
  
  void ListOfSystematics(set<TString>);
  set<TString> ListOfSyst;

 private:

  TString name;

  TString specialComment;
  int bins, datafiles, rebin;
  double rangemin, rangemax, ymin, ymax;

  std::vector<TString> dataset, datasetUp, datasetDown;
  std::vector<double> scales;
  std::vector<TString> legends, legendsSyst;
  std::vector<int> colors;
  std::vector<double> XAxisbins, XAxisbinCenters;

  vector<double> DYScale;
  TString DYEntry;
  TString YAxis;
  TString XAxis;
  TString channel;
  int channelType; //0=ee 1==mumu 2==emu 3==combined  

  std::vector<TH1D> hists;
  std::vector<TH1D> systhistsUp;
  std::vector<TH1D> systhistsDown;

  bool initialized, logX, logY, doDYScale;
  int lumi, signalHist;

  vector<TString> channelLabel;

  double SignalEventswithWeight;
  //static constexpr double topxsec = 220.0; //again changes with normalization; //Use this option when compiling via g++
  static const double topxsec = 213.192; //again changes with normalization;       //Use this option when running it via 'root -l -b -q 'Histo.C+g...'
  
  // DAVID
  bool doUnfolding; 
  bool doSystematics;
  bool drawNLOCurves, drawMadSpinCorr, drawMCATNLO, drawKidonakis, drawPOWHEG;
  TString outpath;
  TString outpathPlots;
  TString outpathResults;
  TString subfolderChannel;
  TString subfolderSpecial;

};


void Plotter::setLumi(double newLumi)
{
    this->lumi = newLumi;
}


// DAVID
void Plotter::UnfoldingOptions(bool doSVD)
{
  doUnfolding = doSVD;
  doSystematics = true;
  drawNLOCurves = true; // boolean to draw/not-draw extra theory curves in the Diff.XSection plots

  drawMadSpinCorr  = true;
  drawMCATNLO      = false;
  drawKidonakis    = false;
  drawPOWHEG       = true;

}


// DAVID
void Plotter::SetOutpath(TString path)
{
  outpath = path;
}


void Plotter::unfolding()
{

    TString sys_array[] = {"DY_","BG_","PU_", "JER", "JES", "TRIG_","MASS", "MATCH", "SCALE","BTAG_ETA_","BTAG_PT_"};
    double sys_array_flat_value[]={0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0};
    TString channel_array[] = {"ee","mumu","emu","combined"};

    vector<TString> vec_systematic (sys_array, sys_array + sizeof(sys_array)/sizeof(sys_array[0]));
    vector<double>  vec_flat_value (sys_array_flat_value, sys_array_flat_value +sizeof(sys_array_flat_value)/sizeof(sys_array_flat_value[0]));
    vector<TString> vec_channel (channel_array, channel_array + sizeof(channel_array)/sizeof(channel_array[0]));

    if(vec_systematic.size() != vec_flat_value.size()){
        cout<<"WARNING (in unfolding)!!!\n You are using 2 different size vectors: 'vec_systematic' and 'vec_flat_value'.\n Fix this and rerun again"<<endl;
        return;
    }

    for(int chan = 0; chan < (int)vec_channel.size(); chan++){ //loop over channels
    
        cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
        cout << "Starting Calculation of Diff. X Section for '" << name << "' in Channel '" << vec_channel.at(chan) << "':" << endl;
        cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
        
        CalcDiffXSec(vec_channel.at(chan),"Nominal");
        
        cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
        cout << "Finished Calculation of Diff. X Section for '" << name << "' in Channel '" << vec_channel.at(chan) << "':" << endl;
        cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;

        if(doSystematics){//############### Syst ################
            for(int sys = 0; sys < (int)vec_systematic.size(); sys++){ //loop over systematics
                cout << endl;
                cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
                cout << "Starting Calculation of Differential Systematics for '" << name << "' in Channel '" << vec_channel.at(chan) << "':" << endl;
                cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
                
                CalcDiffSystematics(vec_channel.at(chan), vec_systematic.at(sys), vec_systematic.at(sys)+"UP",vec_systematic.at(sys)+"DOWN",vec_flat_value.at(sys));
                
                cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
                cout << "Finished Calculation of Differential Systematics for '" << name << "' in Channel '" << vec_channel.at(chan) << endl;
                cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
                cout << endl;
            }//systematic loop
        }
    }//channel loop

}

void Plotter::preunfolding(TString Chan, TString System)
{
    
    //The preunfolding funtion can run now in a single channel and/or single systematic sample
    //for that in just run the Histo.C with the correct parameters. Check the 'Histo' function in Histo.C file
    
    TString sys_array[] = {"Nominal","DY_UP","DY_DOWN","BG_UP","BG_DOWN","PU_UP","PU_DOWN", "JERUP", "JERDOWN", "JESUP", "JESDOWN", "TRIG_UP", "TRIG_DOWN", "MASSUP", "MASSDOWN",
                           "MATCHUP", "MATCHDOWN", "SCALEUP", "SCALEDOWN", "BTAG_ETA_UP", "BTAG_ETA_DOWN", "BTAG_PT_UP", "BTAG_PT_DOWN"};
    TString channel_array[] = {"ee","mumu","emu","combined"};
    vector<TString> vec_systematic (sys_array, sys_array + sizeof(sys_array)/sizeof(sys_array[0]));
    vector<TString> vec_channel (channel_array, channel_array + sizeof(channel_array)/sizeof(channel_array[0]));

    //Check if we run on a single channel or in all of them
    if(Chan!=""){
        vec_channel.clear();
        vec_channel.push_back(Chan);
    }
    if(System!=""){
        vec_systematic.clear();
        vec_systematic.push_back(System);
    }

    for(int i=0; i<(int) vec_channel.size(); i++){//loop over channels
        for(int j=0; j<(int)vec_systematic.size(); j++){//loop over systematics
            write(vec_channel.at(i),vec_systematic.at(j));
        }//systematic loop

    }//channel loop
}

void Plotter::DYScaleFactor(){

    if(!doDYScale){
        DYScale.at(0)=1.;//need to make a switch for control plots that don't want DYScale
        DYScale.at(1)=1.;
        DYScale.at(2)=1.;
        DYScale.at(3)=1.;
        return;
    }

    if(DYScale.size() <4){DYScale.insert(DYScale.begin(), 4, 1.);} //If DYScale has not the appropriate size initialize to (1., 1., 1., 1.)

    cout<<"\n\nBegin DYSCALE FACTOR calculation"<<endl;

    TH1::AddDirectory(kFALSE);

    vector<TString> Vec_Files = InputFileList("combined", "Nominal");//Read the hardcoded list of files
    if(Vec_Files.size()<1) {cout<<"WARNING(in DYScaleFactor)!!! No datasets available to calculate DY SF. EXITING!!"<<endl; return;}
    
    double NoutEEDYMC=0, NinEEDYMC=0, NoutMuMuDYMC=0, NinMuMuDYMC=0;//Number of events in/out of z-veto region for the DY MC
    double NinEE=0, NinMuMu=0, NinEMu=0;//Number of events in z-veto region for data
    double NinEEloose=0, NinMuMuloose=0;//Number of data events in Z-Veto region with MET cut
    double NinEEMC=0, NinMuMuMC=0;//All other MC events

    for(int i=0; i<(int) Vec_Files.size(); i++){
        double LumiWeight = CalcLumiWeight(Vec_Files.at(i));
        double allWeights=LumiWeight;//calculate here all the flat-weights we apply: Lumi*others*...
        TFile *ftemp = TFile::Open(Vec_Files.at(i));
        if(Vec_Files.at(i).Contains("ee") || Vec_Files.at(i).Contains("mumu")){
            if(Vec_Files.at(i).Contains("run")){
                TH1D *htemp = (TH1D*)ftemp->Get("Zh1");
                TH1D *htemp1 = (TH1D*)ftemp->Get("Looseh1");
                ApplyFlatWeights(htemp, allWeights);
                ApplyFlatWeights(htemp1, allWeights);
                if(Vec_Files.at(i).Contains("ee")){
                    NinEE+=htemp->Integral();
                    NinEEloose+=htemp1->Integral();
                }
                if(Vec_Files.at(i).Contains("mumu")){
                    NinMuMu+=htemp->Integral();
                    NinMuMuloose+=htemp1->Integral();
                }
            }
            else if(Vec_Files.at(i).Contains("dy")){
                if(Vec_Files.at(i).Contains("50inf")){
                    TH1D *htemp = (TH1D*)ftemp->Get("Zh1");
                    TH1D *htemp1 = (TH1D*)ftemp->Get("TTh1");
                    ApplyFlatWeights(htemp, LumiWeight);
                    ApplyFlatWeights(htemp1, LumiWeight);
                    if(Vec_Files.at(i).Contains("ee")){
                        NinEEDYMC+=htemp->Integral();
                        NoutEEDYMC+=htemp1->Integral();
                    }
                    if(Vec_Files.at(i).Contains("mumu")){
                        NinMuMuDYMC+=htemp->Integral();
                        NoutMuMuDYMC+=htemp1->Integral();
                    }
                }
                else{
                    TH1D *htemp = (TH1D*)ftemp->Get("TTh1");
                    ApplyFlatWeights(htemp, LumiWeight);
                    if(Vec_Files.at(i).Contains("ee")){   NoutEEDYMC+=htemp->Integral();}
                    if(Vec_Files.at(i).Contains("mumu")){ NoutMuMuDYMC+=htemp->Integral();}
                }
            }
            else{
                TH1D *htemp = (TH1D*)ftemp->Get("Zh1");
                ApplyFlatWeights(htemp, LumiWeight);
                if(Vec_Files.at(i).Contains("ee")){   NinEEMC+=htemp->Integral();   }
                if(Vec_Files.at(i).Contains("mumu")){ NinMuMuMC+=htemp->Integral(); }
            }
        }
        
        if(Vec_Files.at(i).Contains("emu") && Vec_Files.at(i).Contains("run")){
            TH1D *htemp = (TH1D*)ftemp->Get("Zh1");
            ApplyFlatWeights(htemp, LumiWeight);
            NinEMu+=htemp->Integral();
        }

        delete ftemp;
    }
    
    double kee = sqrt(NinEEloose/NinMuMuloose);
    double kmumu = sqrt(NinMuMuloose/NinEEloose);
    
    double RoutinEE = NoutEEDYMC/NinEEDYMC;
    double RoutinMuMu = NoutMuMuDYMC/NinMuMuDYMC;
    
    double NoutMCEE = RoutinEE*(NinEE - 0.5*NinEMu*kee);
    double NoutMCMuMu = RoutinMuMu*(NinMuMu - 0.5*NinEMu*kmumu);

    double DYSFEE = NoutMCEE/NoutEEDYMC;
    double DYSFMuMu = NoutMCMuMu/NoutMuMuDYMC;

    cout << endl;
    cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
    cout << "Calculation of DY Scale Factors for '" << name << "'  " << endl;

    cout<<"DYSFEE:                 "<<DYSFEE<<endl;
    cout<<"DYSFMuMu:               "<<DYSFMuMu<<endl;

    cout<<"NinEEloose:             "<<NinEEloose<<endl;
    cout<<"NinMMloose:             "<<NinMuMuloose<<endl;

    cout<<"kee:                    "<<kee<<" +- "<<0.5*TMath::Sqrt(1./NinMuMuloose + 1./NinEEloose)<<endl;
    cout<<"kmumu:                  "<<kmumu<<" +- "<<0.5*TMath::Sqrt(1./NinMuMuloose + 1./NinEEloose)<<endl;

    cout<<"Rout/Rin ee:            "<<RoutinEE<<endl;
    cout<<"Rout/Rin Mumu:          "<<RoutinMuMu<<endl;

    cout<<"Est. From Data(ee):     "<<NoutMCEE<<endl;
    cout<<"Est. From Data(mumu):   "<<NoutMCMuMu<<endl;

    cout<<"Est. From MC(ee):       "<<NoutEEDYMC<<endl;
    cout<<"Est. From MC(mumu):     "<<NoutMuMuDYMC<<endl;

    cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
    cout << endl;

    
    DYScale.at(0)=DYSFEE;
    DYScale.at(1)=DYSFMuMu;
    DYScale.at(2)=1.;
    DYScale.at(3)=(DYSFEE+DYSFMuMu)/2;//not correct, but close, fix later

    cout<<"End DYSCALE FACTOR calculation\n"<<endl;

}



void Plotter::CalcDiffSystematics(TString Channel, TString Systematic, TString SystematicUp, TString SystematicDown, double flat_Syst){

    cout << endl;
    cout << endl;
    cout << "    Preparing to Calculate " << Systematic << "-Uncertainty ... " << endl;

    ofstream ResultsFile;
    gSystem->MakeDirectory("UnfoldingResults");
    gSystem->MakeDirectory("UnfoldingResults/"+Systematic);
    gSystem->MakeDirectory("UnfoldingResults/"+Systematic+"/"+Channel);
    
    string ResultsFilestring = outpathResults.Data();
    ResultsFilestring.append(subfolderSpecial.Data());
    ResultsFilestring.append("/");
    ResultsFilestring.append(Systematic);
    ResultsFilestring.append("/");
    ResultsFilestring.append(Channel);
    ResultsFilestring.append("/");
    ResultsFilestring.append(name);
    ResultsFilestring.append("Results.txt");
    cout<<ResultsFilestring<<endl;
    ResultsFile.open(ResultsFilestring.c_str());
        
    double Xbins[XAxisbins.size()];
    for(unsigned int i = 0; i<XAxisbins.size();i++){Xbins[i]=XAxisbins[i];} 

    if(flat_Syst > 0.0){// flat systematics will be: flat_Syst != 0.0
        for(int bin = 0; bin<(int)XAxisbinCenters.size(); bin++){
            ResultsFile<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[bin]<<" bin: "<<Xbins[bin]<<" to "<<Xbins[bin+1]<<" SystematicError: "<<flat_Syst<<endl;
        }
    }

    double Sys_Error;//, Sum_Errors;

    TString newname = name;
    if(name.Contains("Hyp")){//Histogram naming convention has to be smarter
        newname.ReplaceAll("Hyp",3,"",0);
    }
    
    if(Channel=="ee"){channelType=0;}
    if(Channel=="mumu"){channelType=1;}
    if(Channel=="emu"){channelType=2;}
    if(Channel=="combined"){channelType=3;}
    
    // DAVID. Guckst du hier! 
    if ( doUnfolding == true && flat_Syst == 0.0 && channelType != 3) {//Only do it in non combined channel
    
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


        TFile *ftemp = TFile::Open("preunfolded/Nominal/"+Channel+"/"+name+"_UnfoldingHistos.root");
        theDataHist =  (TH1D*)ftemp->Get("aDataHist")->Clone();
        theBgrHist =  (TH1D*)ftemp->Get("aBgrHist")->Clone();
        theTtBgrHist =  (TH1D*)ftemp->Get("aTtBgrHist")->Clone();
        theRecHist =  (TH1D*)ftemp->Get("aRecHist")->Clone();
        theGenHist =  (TH1D*)ftemp->Get("aGenHist")->Clone();
        theRespHist =  (TH2D*)ftemp->Get("aRespHist")->Clone();
        delete ftemp;

        TFile *ftempUp = TFile::Open("preunfolded/"+SystematicUp+"/"+Channel+"/"+name+"_UnfoldingHistos.root");
        theBgrHistUp =  (TH1D*)ftempUp->Get("aBgrHist")->Clone();
        theTtBgrHistUp =  (TH1D*)ftempUp->Get("aTtBgrHist")->Clone();
        theRecHistUp =  (TH1D*)ftempUp->Get("aRecHist")->Clone();
        theGenHistUp =  (TH1D*)ftempUp->Get("aGenHist")->Clone();
        theRespHistUp =  (TH2D*)ftempUp->Get("aRespHist")->Clone();
        delete ftempUp;

        TFile *ftempDown = TFile::Open("preunfolded/"+SystematicDown+"/"+Channel+"/"+name+"_UnfoldingHistos.root");
        theBgrHistDown =  (TH1D*)ftempDown->Get("aBgrHist")->Clone();
        theTtBgrHistDown =  (TH1D*)ftempDown->Get("aTtBgrHist")->Clone();
        theRecHistDown =  (TH1D*)ftempDown->Get("aRecHist")->Clone();
        theGenHistDown =  (TH1D*)ftempDown->Get("aGenHist")->Clone();
        theRespHistDown =  (TH2D*)ftempDown->Get("aRespHist")->Clone();
        delete ftempDown;

        // Apply Scale Factor for MC@NLO 
        /*    if ( legendsSyst.size() > 0 ) {  
        ApplyMCATNLOWeight(theRespHistUp, legendsSyst.back(), "Up", "ttbarsignal");  
        ApplyMCATNLOWeight(theRespHistDown, legendsSyst.back(), "Down", "ttbarsignal");  
        ApplyMCATNLOWeight(theGenHistUp, legendsSyst.back(), "Up", "ttbarsignal");  
        ApplyMCATNLOWeight(theGenHistDown, legendsSyst.back(), "Down", "ttbarsignal");  
        }*/
                    
        // Get the binning 
        double* theBins = Xbins;
        int numberBins = bins;
                    
        // Names and Labels
        TString channelLabelStr(channelLabel[channelType]);
        TString theChannelName = Channel;
        //if ( channelLabelStr.Contains("#mu#mu")  ) theChannelName = "mumu";
        //if ( channelLabelStr.Contains("e#mu")    ) theChannelName = "emu";
        //if ( channelLabelStr.Contains("ee")      ) theChannelName = "ee";
        //if ( channelLabelStr.Contains("Dilepton Combined")    ) theChannelName = "combined";
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
        
        
        
        
        //Symetrize Eta and Rapidity distributions
        if (theQuantityName == "Eta" || theQuantityName == "Rapidity" ){
            for(int j=0; j<(int) symmSysErrors->GetNbinsX(); ++j){
                cout<<"In bin "<<j<<" binCenter "<<symmSysErrors->GetBinCenter(j+1)<<" Content "<<symmSysErrors->GetBinContent(j+1)<<endl;
            }

            int Nbins = theDataHist->GetNbinsX();
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

                //cout<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[i]<<" bin: "<<Xbins[i]<<" to "<<Xbins[i+1]<<" SystematicError: "<<Sys_Error<<endl;
                ResultsFile<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[i]<<" bin: "<<Xbins[i]<<" to "<<Xbins[i+1]<<" SystematicRelError: "<<Sys_Error<<endl;
            }
        }
        else{
            // Save the shifts in Tyler's triple-matrix ...
            for(Int_t bin = 0; bin < theDataHist->GetNbinsX(); ++bin) {
                Sys_Error = symmSysErrors->GetBinContent(bin+2); // Keep in mind the extra layer of OF bins
                if(Systematic == "MASS"){
                    Sys_Error = Sys_Error/12.;
                }
                // Save it
                //cout<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[bin]<<" bin: "<<Xbins[bin]<<" to "<<Xbins[bin+1]<<" SystematicError: "<<Sys_Error<<endl;
                ResultsFile<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[bin]<<" bin: "<<Xbins[bin]<<" to "<<Xbins[bin+1]<<" SystematicRelError: "<<Sys_Error<<endl;
            }
        }
    }
    
    
    
    if (channelType==3){//for 'combined' channel: do an statistical combination of the the 3 independent channels using the central value's statistical error as weight
        
        TString eefilename="UnfoldingResults/Nominal/ee/"+name+"Results.txt";
        TString emufilename="UnfoldingResults/Nominal/emu/"+name+"Results.txt";
        TString mumufilename="UnfoldingResults/Nominal/mumu/"+name+"Results.txt";
        TString combfilename="UnfoldingResults/Nominal/combined/"+name+"Results.txt";
        TString eeErrfilename="UnfoldingResults/"+Systematic+"/ee/"+name+"Results.txt";
        TString emuErrfilename="UnfoldingResults/"+Systematic+"/emu/"+name+"Results.txt";
        TString mumuErrfilename="UnfoldingResults/"+Systematic+"/mumu/"+name+"Results.txt";
        
        //check the existence of the file
        if(gSystem->AccessPathName(eefilename) || gSystem->AccessPathName(emufilename) || gSystem->AccessPathName(mumufilename) || gSystem->AccessPathName(eeErrfilename) || gSystem->AccessPathName(emuErrfilename) || gSystem->AccessPathName(mumuErrfilename)){
            cout<<"WARNING (in CalcDiffSystematics)!!"<<endl;
            cout<<"One of the input files you use for the combined XSection measurement doesn't exist!! Exiting!!"<<endl;
            return;
        }
        
        ifstream ResultsEE(eefilename);
        ifstream ResultsEMu(emufilename);
        ifstream ResultsMuMu(mumufilename);
        ifstream ResultsComb(combfilename);
        ifstream SystErrEE(eeErrfilename);
        ifstream SystErrEMu(emuErrfilename);
        ifstream SystErrMuMu(mumuErrfilename);
        
        double perChannelDiffXSecPlot[4][bins];      //perChannelDiffXSecPlot[channel][bin]
        double perChannelDiffXSecSystError[3][bins];      //perChannelDiffXSecPlot[channel][bin]
        double perChannelDiffXSecStatError[4][bins]; //perChannelDiffXSecStatError[channel][bin]
        double perChannelGenDiffXSec[4][bins];       //perChannelGenDiffXSec[channel][bin]
        TString Dummy="";
        
        for (Int_t bin=0; bin<bins; bin++){//Retrieve arrays for plotting
            ResultsEE>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecPlot[0][bin]>>Dummy>>perChannelDiffXSecStatError[0][bin]>>Dummy>>perChannelGenDiffXSec[0][bin];
            ResultsEMu>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecPlot[1][bin]>>Dummy>>perChannelDiffXSecStatError[1][bin]>>Dummy>>perChannelGenDiffXSec[1][bin];
            ResultsMuMu>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecPlot[2][bin]>>Dummy>>perChannelDiffXSecStatError[2][bin]>>Dummy>>perChannelGenDiffXSec[2][bin];
            ResultsComb>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecPlot[3][bin]>>Dummy>>perChannelDiffXSecStatError[3][bin]>>Dummy>>perChannelGenDiffXSec[3][bin];
            SystErrEE>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecSystError[0][bin];
            SystErrEMu>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecSystError[1][bin];
            SystErrMuMu>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecSystError[2][bin];
        }
        ResultsEE.close(); ResultsEMu.close(); ResultsMuMu.close();

        double binWidth[bins];
        double DiffXSecVecVaried[bins];
        //do the actual combined Diff.XSection calculation
        for(int i=0; i<bins; i++){
            binWidth[i] = Xbins[i+1]-Xbins[i];
            for(int j=0; j<3; j++){//check if any stat error is 0, in this case set their contribution to 0!!
                if(perChannelDiffXSecStatError[j][i] == 0){
                    perChannelDiffXSecStatError[j][i] = 1e100;
                    perChannelDiffXSecPlot[j][i] = 0;
                }
            }
            DiffXSecVecVaried[i] =((perChannelDiffXSecSystError[0][i] * perChannelDiffXSecPlot[0][i]/(perChannelDiffXSecStatError[0][i]*perChannelDiffXSecStatError[0][i]))
                                  +(perChannelDiffXSecSystError[1][i] * perChannelDiffXSecPlot[1][i]/(perChannelDiffXSecStatError[1][i]*perChannelDiffXSecStatError[1][i]))
                                  +(perChannelDiffXSecSystError[2][i] * perChannelDiffXSecPlot[2][i]/(perChannelDiffXSecStatError[2][i]*perChannelDiffXSecStatError[2][i])))/
                                  ((1/(perChannelDiffXSecStatError[0][i]*perChannelDiffXSecStatError[0][i]))+
                                   (1/(perChannelDiffXSecStatError[1][i]*perChannelDiffXSecStatError[1][i]))+
                                   (1/(perChannelDiffXSecStatError[2][i]*perChannelDiffXSecStatError[2][i])));
            
            double Sys_Error = 1.-(TMath::Abs(DiffXSecVecVaried[i]-perChannelDiffXSecPlot[3][i])/perChannelDiffXSecPlot[3][i]);
            ResultsFile<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[i]<<" bin: "<<Xbins[i]<<" to "<<Xbins[i+1]<<" SystematicRelError: "<<Sys_Error<<endl;
        }
    }
    
    ResultsFile.close();
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

    channelLabel.insert(channelLabel.begin(),4, "");

    // DAVID
    outpath = "";
    outpathPlots = "Plots";
    outpathPlots = "UnfoldingResults";
    subfolderChannel = "";
    subfolderSpecial = "";
    
    
    //Ivan: Initialize list of systematics

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

void Plotter::setOptions(TString name_, TString specialComment_, TString YAxis_, TString XAxis_, int rebin_, bool doDYScale_, bool logX_, bool logY_, double ymin_, double ymax_, double rangemin_, double rangemax_, int bins_, std::vector<double> XAxisbins_, std::vector<double> XAxisbinCenters_)
{
    name=name_; //Variable name you want to plot
    specialComment=specialComment_; 
    YAxis=YAxis_; //Y-axis title
    XAxis=XAxis_; //X-axis title
    rebin=rebin_; //Nr. of bins to be merged together
    doDYScale = doDYScale_; //Apply DY scale factor?
    logX = logX_; //Draw X-axis in Log scale
    logY = logY_; //Draw Y-axis in Log scale
    ymin = ymin_; //Min. value in Y-axis
    ymax = ymax_; //Max. value in Y-axis
    rangemin=rangemin_; //Min. value in X-axis
    rangemax=rangemax_; //Max. value in X-axis
    bins=bins_; //Number of bins to plot
    XAxisbins.clear(); 
    XAxisbins = XAxisbins_; // Bins edges=bins+1
    XAxisbinCenters.clear();
    XAxisbinCenters = XAxisbinCenters_; //Central point for BinCenterCorrection=bins

    //Modify the X/Y-axis labels
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

    DYScale.insert(DYScale.begin(), 4, 1.);//Initialize the DY scale-factor to (1., 1., 1., 1.)
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

void Plotter::setDataSet(TString mode, TString Systematic)
{
    initialized=false;
    legendsSyst.clear();

    if(channelLabel.size()<4){channelLabel.insert(channelLabel.begin(), 4, "");}

    if(mode =="ee"){channelType=0;channelLabel.at(0)="ee";}
    if(mode =="mumu"){channelType=1;channelLabel.at(1)="#mu#mu";}
    if(mode =="emu"){channelType=2;channelLabel.at(2)="e#mu";}
    if(mode =="combined"){channelType=3;channelLabel.at(3)="Dilepton Combined";}

    // Set dataset specific subfolders
    outpathPlots = "./Plots";
    outpathResults = "./UnfoldingResults";
    subfolderChannel = mode;
    subfolderChannel.Prepend("/");
    subfolderSpecial = "";
    if ( specialComment.CompareTo("Standard") != 0 ) { subfolderSpecial = specialComment.Prepend("/");}

    DYEntry = "Z / #gamma* #rightarrow ee/#mu#mu";

    if(Systematic.Contains("DY_") || Systematic.Contains("BG_")){Systematic = "Nominal";}//We just need to vary the nominal DY and BG systematics

    ifstream FileList("FileLists/HistoFileList_"+Systematic+"_"+mode+".txt");
    TString filename;
    datafiles=0;

    dataset.clear();
    legends.clear();
    colors.clear();

    while(!FileList.eof()){
        FileList>>filename;
        if(filename==""){continue;}//Skip empty lines
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


void Plotter::ListOfSystematics(set<TString> listofsystematics){
    
    //List of systematics valid for running. Please check (or add any new systematic) the function 'ListOfSystematics' in Histo.C
    
    ListOfSyst = listofsystematics;
    
}


std::vector<TString> Plotter::InputFileList(TString mode, TString Systematic)
{
    //Hard code the input file list. This functions handles correctly the data files. So no symlinks are anymore needed
    //The file structure must be: selectionRoot/SYSTEMATIC/channel/XXX.root
    //This function will take the data (aka: run2012A, run2012B, run2012C, ...) from the Nominal systematic for ALL the systematic
    //This function will take every sample, except the ttbar, from Nominal for any signal systematic: POWHEG, MCATNLO, SPINCORR, MATCH, SCALE, MASS
    
    std::vector<TString> FileVector;
    FileVector.clear();
    
    if(ListOfSyst.find(Systematic) == ListOfSyst.end()){
        cout<<"WARNING(in InputFileList)!!! Using a non valid systematic: "<<Systematic<<endl;
        cout<<"Please use one among the supported ones:"<<endl;
        for(set<TString>::iterator iter = ListOfSyst.begin(); iter!=ListOfSyst.end(); ++iter){ cout<<*iter<<endl;}
        return FileVector;
    }
    
    if( mode.CompareTo("combined") && mode.CompareTo("ee") && mode.CompareTo("emu") && mode.CompareTo("mumu")){
        cout<<"The decay channel you provided is not supported."<<endl;
        cout<<"Please use: ee, emu, mumu, combined"<<endl;
        return FileVector;
    }
    
    if(!mode.CompareTo("combined")){
        vector<TString> eemode   = Plotter::InputFileList(TString("ee"), Systematic);
        vector<TString> emumode  = Plotter::InputFileList(TString("emu"), Systematic);
        vector<TString> mumumode = Plotter::InputFileList(TString("mumu"), Systematic);
        FileVector.insert(FileVector.end(), eemode.begin(), eemode.end());
        FileVector.insert(FileVector.end(), emumode.begin(), emumode.end());
        FileVector.insert(FileVector.end(), mumumode.begin(), mumumode.end());
        return FileVector;
    }

    //data is only stored in the Nominal directory
    FileVector.push_back(TString("selectionRoot/Nominal/").Append(mode.Copy()).Append("/").Append(mode.Copy()).Append("_run2012A.root"));
    FileVector.push_back(TString("selectionRoot/Nominal/").Append(mode.Copy()).Append("/").Append(mode.Copy()).Append("_run2012Arecover.root"));
    FileVector.push_back(TString("selectionRoot/Nominal/").Append(mode.Copy()).Append("/").Append(mode.Copy()).Append("_run2012B.root"));
    FileVector.push_back(TString("selectionRoot/Nominal/").Append(mode.Copy()).Append("/").Append(mode.Copy()).Append("_run2012C_24Aug.root"));
    FileVector.push_back(TString("selectionRoot/Nominal/").Append(mode.Copy()).Append("/").Append(mode.Copy()).Append("_run2012C_PromptReco.root"));
    
    //MC depends on the specific Systematic: Signal systematics only use different signal samples
    TString tempName;
    if(!Systematic.CompareTo("JERUP") || !Systematic.CompareTo("JERDOWN") || !Systematic.CompareTo("JESUP") || !Systematic.CompareTo("JESDOWN") ||
       !Systematic.CompareTo("PU_UP") || !Systematic.CompareTo("PU_DOWN") || !Systematic.CompareTo("BTAG_UP") || !Systematic.CompareTo("BTAG_DOWN") ||
       !Systematic.CompareTo("BTAG_PT_UP") || !Systematic.CompareTo("BTAG_PT_DOWN") || !Systematic.CompareTo("BTAG_ETA_UP") || !Systematic.CompareTo("BTAG_ETA_DOWN")
    ){
        tempName = TString("selectionRoot/") + Systematic + "/" + mode + "/" + mode;
    }
    else{
        tempName = TString("selectionRoot/Nominal/") + mode + "/" + mode;
    }
    
    FileVector.push_back(tempName + "_dyee1050.root");
    FileVector.push_back(tempName + "_dyee50inf.root");
    FileVector.push_back(tempName + "_dymumu1050.root");
    FileVector.push_back(tempName + "_dymumu50inf.root");
    FileVector.push_back(tempName + "_singleantitop_tw.root");
    FileVector.push_back(tempName + "_singletop_tw.root");
    FileVector.push_back(tempName + "_wtolnu.root");
    FileVector.push_back(tempName + "_wwtoall.root");
    FileVector.push_back(tempName + "_wztoall.root");
    FileVector.push_back(tempName + "_zztoall.root");
    FileVector.push_back(tempName + "_dytautau1050.root");
    FileVector.push_back(tempName + "_dytautau50inf.root");
    FileVector.push_back(tempName + "_qcdem2030.root");
    FileVector.push_back(tempName + "_qcdem3080.root");
    FileVector.push_back(tempName + "_qcdem80170.root");
    FileVector.push_back(tempName + "_qcdmu15.root");
    
    TString ttbarsignalplustau = TString("selectionRoot/") + Systematic + "/" + mode + "/" + mode + "_ttbarsignalplustau.root";
    TString ttbarbgviatau      = TString("selectionRoot/") + Systematic + "/" + mode + "/" + mode + "_ttbarbgviatau.root";
    TString ttbarbg            = TString("selectionRoot/") + Systematic + "/" + mode + "/" + mode + "_ttbarbg.root";

    //Add extra filename for signal systematic filenames
    if(!Systematic.CompareTo("POWHEG") || !Systematic.CompareTo("MCATNLO") || 
       !Systematic.CompareTo("MASSUP") || !Systematic.CompareTo("MASSDOWN") ||
       !Systematic.CompareTo("MATCHINGUP") || !Systematic.CompareTo("MATCHINGDOWN") ||
       !Systematic.CompareTo("MATCHUP") || !Systematic.CompareTo("MATCHDOWN") ||
       !Systematic.CompareTo("SCALEUP") || !Systematic.CompareTo("SCALEDOWN")){
        
        if(!Systematic.CompareTo("MATCHINGUP") || !Systematic.CompareTo("MATCHINGDOWN")){ Systematic.Remove(5, 3);}
        TString tmpSyst = Systematic;
        tmpSyst.Prepend("_").ToLower();
        
        ttbarsignalplustau.ReplaceAll(".root", tmpSyst + ".root");
        ttbarbg.ReplaceAll(".root", tmpSyst + ".root");
        ttbarbgviatau.ReplaceAll(".root", tmpSyst + ".root");
    }
    
    if(!Systematic.CompareTo("SPINCORR")){
        ttbarsignalplustau.ReplaceAll(".root", "_FullLeptMadgraphWithSpinCorrelation.root");
        ttbarbgviatau.ReplaceAll(".root", "_FullLeptMadgraphWithSpinCorrelation.root");
        ttbarbg.ReplaceAll(".root", "_HadronicMadgraphWithSpinCorrelation.root");
        FileVector.push_back(ttbarbg);
        ttbarbg.ReplaceAll("Hadronic", "SemiLept");
    }
    
    FileVector.push_back(ttbarsignalplustau);
    FileVector.push_back(ttbarbg);
    FileVector.push_back(ttbarbgviatau);
    
    return FileVector;
    
}


bool Plotter::fillHisto()
{   
    if (initialized) { return true; }
    TH1::AddDirectory(kFALSE);
    hists.clear();
    for(unsigned int i=0; i<dataset.size(); i++){
        TFile *ftemp = TFile::Open(dataset[i]);
        if (!ftemp) {
            cerr << "*** Could not open file " << dataset[i] << "!" << endl;
            return false;
        }
        TH1D *hist = dynamic_cast<TH1D*>(ftemp->Get(name));
        if (!hist) {
            cerr << "*** Could not find histogram " << name << " in file " << dataset[i] << "!" << endl;
            return false;
        }
        hist = static_cast<TH1D*>(hist->Clone());
        if (!name.Contains("bcp_") && !name.Contains("Lead") 
            && (name.Contains("Lepton") || name.Contains("BJet") || name.Contains("Top")))
        {
            TString stemp = name;
            if(name.Contains("Lepton"))    {stemp.ReplaceAll("Lepton",6,"AntiLepton",10);}
            else if(name.Contains("BJet")) {stemp.ReplaceAll("BJet",4,"AntiBJet",8);}
            else if(name.Contains("Top"))  {stemp.ReplaceAll("Top",3,"AntiTop",7);}
            TH1D *other = dynamic_cast<TH1D*>(ftemp->Get(stemp));
            if (other) hist->Add(other); 
            else cerr << "Cannot find corresponding anti-quantity histogram: "<<stemp<<endl;
        }

        //Rescaling to the data luminosity
        double LumiWeight = CalcLumiWeight(dataset[i]);
        ApplyFlatWeights(hist, LumiWeight);

//         //Apply any other flat weight (still to do: define the weights somewhere else)
//         ApplyFlatWeights(hist, lepSF);
//         ApplyFlatWeights(hist, KinRecSF);
    
        setHHStyle(*gStyle);

        hists.push_back(*hist);
        delete ftemp;
    }
    initialized=true;
    return true;
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
    //    cout << "    Channel:              " << mode << endl;
    cout << "    Quantity:             " << name << endl; 
    cout << "    The factor is:        " << MCATNLO_ScaleFactor << endl;
    cout << "    The inverse of it is: " << MCATNLO_ScaleFactorInv << endl;
    cout << endl;
  
    // Apply the weight
    hist->Scale(MCATNLO_ScaleFactorInv);
  
}


void Plotter::write(TString Channel, TString Systematic) // do scaling, stacking, legending, and write in file 
{
    setDataSet(Channel,Systematic);
    if (!fillHisto()) return;

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
    TH1* aDataHist = NULL;
    TH1* aBgrHist = NULL;
    TH1* aTtBgrHist = NULL;
    TH1* aRecHist = NULL;
    TH1* aGenHist = NULL; 
    TH1* aRespHist = NULL;


    double Xbins[XAxisbins.size()];
    TString newname = name;

    if(name.Contains("Hyp")){//Histogram naming convention has to be smarter
        newname.ReplaceAll("Hyp",3,"",0);
    }

    bool init=false;

    for(unsigned int i = 0; i<XAxisbins.size();i++){Xbins[i]=XAxisbins[i];}

    for(unsigned int i=0; i<hists.size() ; i++){ // prepare histos and leg
        drawhists[i]=(TH1D*) hists[i].Clone();//rebin and scale the histograms
        if(rebin>1) drawhists[i]->Rebin(rebin);
        setStyle(drawhists[i], i, true);
    }


    gSystem->MakeDirectory("preunfolded");
    gSystem->MakeDirectory("preunfolded/"+Systematic);
    gSystem->MakeDirectory("preunfolded/"+Systematic+"/"+Channel);


    for(unsigned int i=0; i<hists.size() ; i++){ // prepare histos and leg
    if(legends[i] != "Data"){
        //      drawhists[i]->Scale(12.1/5.1);

        if(XAxisbins.size()>1){//only distributions we want to unfold will have a binning vector
            if(legends[i] == "t#bar{t} Signal"){
                TFile *ftemp = TFile::Open(dataset[i]);
                if(init==false){
                    aRespHist = (TH2*)ftemp->Get("GenReco"+newname)->Clone();
                    aGenHist=(TH1D*)ftemp->Get("VisGen"+newname)->Clone();
                    //Rebin(bins,"aTtBgrHist",Xbins);
                    if(!newname.Contains("Lead") && (newname.Contains("Lepton")||newname.Contains("Top")||newname.Contains("BJet"))){
                        aRespHist->Add((TH2*)ftemp->Get("GenRecoAnti"+newname)->Clone());
                        aGenHist->Add((TH1D*)ftemp->Get("VisGenAnti"+newname)->Clone());
                    }
                    init =true;
                } else {//account for more than one signal histogram
                    aRespHist->Add((TH2*)ftemp->Get("GenReco"+newname)->Clone());
                    aGenHist->Add((TH1D*)ftemp->Get("VisGen"+newname)->Clone());
                    if(!newname.Contains("Lead") && (newname.Contains("Lepton")||newname.Contains("Top")||newname.Contains("BJet"))){
                        aGenHist->Add((TH1D*)ftemp->Get("VisGenAnti"+newname)->Clone());
                        aRespHist->Add((TH2*)ftemp->Get("GenRecoAnti"+newname)->Clone());
                    }
                }
                delete ftemp;
            }
            //cout<<"Legend: "<<legends[i]<<endl;
            if(legends[i] == "t#bar{t} Signal"){
                signalHist = i; //What is happening for the combined channel?? only 1 integer value??
                aRecHist = drawhists[i]->Rebin(bins,"aRecHist",Xbins); //What happens in the combined channel? Are the 3 channels summed up?? NO!!
                //cout<<"Added "<<legends[i]<<" to aRecHist"<<endl;
            }else if(legends[i] == "t#bar{t} Other"){//IMPORTANT: TTbar Other are added to the ttbarbackground histogram AND the Background Hist gram
                if(aTtBgrHist == NULL){//check to see if this has been created
                    aTtBgrHist = drawhists[i]->Rebin(bins,"aTtBgrHist",Xbins);
                    //cout<<"Created "<<legends[i]<<" to aTtBgrHist"<<endl;
                }else{
                    aTtBgrHist->Add(drawhists[i]->Rebin(bins,"aTtBgrHist",Xbins));
                    //cout<<"Added "<<legends[i]<<" to aTtBgrHist"<<endl;
                }
                if(aBgrHist == NULL){
                    aBgrHist = drawhists[i]->Rebin(bins,"aBgrHist",Xbins);
                    //cout<<"Created "<<legends[i]<<" to aBgrHist"<<endl;
                }else{
                    aBgrHist->Add(drawhists[i]->Rebin(bins,"aBgrHist",Xbins));
                    //cout<<"Added "<<legends[i]<<" to aBgrHist"<<endl;
                }
            }else if((legends[i] == DYEntry)){
                if (channelType!=2) drawhists[i]->Scale(DYScale[channelType]);

                //Here we take into account the systematic shifts needed for DY systematic because it only modifies the nominal dataset
                if(Systematic == "DY_UP"){
                    drawhists[i]->Scale(1.3);
                }
                if(Systematic == "DY_DOWN"){
                    drawhists[i]->Scale(0.7);
                }
                if(aBgrHist == NULL){
                    aBgrHist = drawhists[i]->Rebin(bins,"aBgrHist",Xbins);
                    //cout<<"Created "<<legends[i]<<" to aBgrHist"<<endl;
                }else{
                    aBgrHist->Add(drawhists[i]->Rebin(bins,"aBgrHist",Xbins));
                    //cout<<"Added "<<legends[i]<<" to aBgrHist"<<endl;
                }
            }else{

                //Here we take into account the systematic shifts needed for BG systematic because it only modifies the nominal dataset
                if(Systematic == "BG_UP"){
                    drawhists[i]->Scale(1.5);
                }
                if(Systematic == "BG_DOWN"){
                    drawhists[i]->Scale(0.5);
                }

                if(aBgrHist == NULL){
                    aBgrHist = drawhists[i]->Rebin(bins,"aBgrHist",Xbins);
                    //cout<<"Created "<<legends[i]<<" to aBgrHist"<<endl;
                }else{
                    aBgrHist->Add(drawhists[i]->Rebin(bins,"aBgrHist",Xbins));
                    //cout<<"Added "<<legends[i]<<" to aBgrHist"<<endl;
                    //cout<<"integral: "<<aBgrHist->Integral()<<endl;
                }
            }
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


    if(XAxisbins.size()>1){//only distributions we want to unfold will have a binning vector
        aDataHist = drawhists[0]->Rebin(bins,"aDataHist",Xbins);
        //  cout<<"Added data to aDataHist"<<endl;

        TFile *f15 = new TFile("preunfolded/"+Systematic+"/"+Channel+"/"+name+"_UnfoldingHistos.root","RECREATE");
        aDataHist->Write("aDataHist");
        aTtBgrHist->Write("aTtBgrHist");
        aBgrHist->Write("aBgrHist");
        aGenHist->Write("aGenHist");
        aRespHist->Write("aRespHist");
        aRecHist->Write("aRecHist");

        f15->Close();
    }

    leg = ControlLegend(hists.size(), drawhists, legends, leg);

    if(name.Contains("HypjetMultiXSec")){

        double InclusiveXsectionWrite[4], InclusiveXsectionStatErrorWrite[4];
        CalcXSec(dataset, InclusiveXsectionWrite, InclusiveXsectionStatErrorWrite, "","");

        cout<<"&&&&&&&&&&&%%%%%%%%%%^^^^^^^^^^^^"<<endl;
        cout<<"&&&&&&&&&&&%%%%%%%%%%^^^^^^^^^^^^"<<endl;
        cout<<"&&&&&&&&&&&%%%%%%%%%%^^^^^^^^^^^^"<<endl;
        cout<<"Look ma! I calculated an inclusive cross-section! It is "<<InclusiveXsectionWrite[channelType]<<endl;
        cout<<"With a stat uncertainty of: "<<InclusiveXsectionStatErrorWrite[channelType]<<endl;

        ofstream InclusiveXSecResult("preunfolded/"+Systematic+"/"+Channel+"/InclusiveXSection.txt");  
        InclusiveXSecResult<<"InclusiveXSection: "<<InclusiveXsectionWrite[channelType]<<" StatError: "<<InclusiveXsectionStatErrorWrite[channelType]<<endl;
        InclusiveXSecResult.close();
    }


    TList* l = stack->GetHists(); 
    TH1D* stacksum = (TH1D*) l->At(0)->Clone();

    for (int i = 1; i < l->GetEntries(); ++i) { stacksum->Add((TH1D*)l->At(i));}
    //  f0->Close();
    //stat uncertainty::make a function 
    TH1D* syshist =0;
    syshist = (TH1D*)stacksum->Clone();
//     double lumierr = 0.045; 
    for(Int_t i=0; i<=syshist->GetNbinsX(); ++i){

        Double_t binc = 0;
        binc += stacksum->GetBinContent(i);
        syshist->SetBinContent(i, binc);
/*        // calculate uncertainty: lumi uncertainty
         Double_t binerr2 = binc*binc*lumierr*lumierr;
         Double_t topunc = 0; // uncertainty on top xsec

        //Kidonakis
         double topxsecErr2 = 2.2*2.2 + 4.4*4.4 + 5.5*5.5; //topxsecErr2 = lumiErr*lumiErr + topxsecScaleErr*topxsecScaleErr + topxsecPDFErr*topxsecPDFErr

         double topRelUnc =  TMath::Sqrt(topxsecErr2)/topxsec;
        //    topunc += drawhists[signalHist]->GetBinContent(i)*topRelUnc;
        // binerr2 += (topunc*topunc);
        // syshist->SetLineColor(1);
        // syshist->SetBinError(i, TMath::Sqrt(binerr2));
*/
    }

    if(logY)c->SetLogy();

    syshist->SetFillStyle(3004);
    syshist->SetFillColor(kBlack);

    //leg->AddEntry( syshist, "Uncertainty", "f" );

    drawhists[0]->SetMinimum(ymin);

    if(rangemin!=0 || rangemax!=0) {drawhists[0]->SetAxisRange(rangemin, rangemax, "X");}

    if(ymax==0){
        if(logY){drawhists[0]->SetMaximum(18  * drawhists[0]->GetBinContent(drawhists[0]->GetMaximumBin()));}
        else    {drawhists[0]->SetMaximum(1.5 * drawhists[0]->GetBinContent(drawhists[0]->GetMaximumBin()));}
    }
    else{drawhists[0]->SetMaximum(ymax);}

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
    drawhists[0]->Draw("e1");

    stack->Draw("same HIST");
    gPad->RedrawAxis();
    TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0.5)");//this is frustrating and stupid but apparently necessary...
    setex1->Draw();
    syshist->SetMarkerStyle(0);
    //syshist->Draw("same,E2");
    TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.)");
    setex2->Draw();
    drawhists[0]->Draw("same,e1");

    DrawCMSLabels(2, lumi, 8);
    DrawDecayChLabel(channelLabel[channelType]);
    leg->Draw("SAME");
    //drawRatio(drawhists[0], stacksum, 0.75, 1.5, *gStyle);

    // Create Directory for Output Plots 
    gSystem->MakeDirectory(outpathPlots);
    gSystem->MakeDirectory(outpathPlots+"/"+subfolderChannel);
    gSystem->MakeDirectory(outpathPlots+"/"+subfolderChannel+"/"+Systematic);
    c->Print(outpathPlots+subfolderChannel+"/"+Systematic+"/"+name+".eps");
    //c->Print(outpathPlots+subfolderChannel+"/"+Systematic+"/"+name+".C");
    c->Clear();
    leg->Clear();
    delete c;
    delete leg;
}

void Plotter::setStyle(TH1 *hist, unsigned int i, bool isControlPlot)
{
    hist->SetFillColor(colors[i]);
    hist->SetLineColor(colors[i]);
    hist->SetLineWidth(1);

    if(legends[i] == "Data"){
        hist->SetFillColor(0);
        hist->SetMarkerStyle(20);
        hist->SetMarkerSize(1.);
        hist->SetLineWidth(1);
        hist->GetXaxis()->SetLabelFont(42);
        hist->GetYaxis()->SetLabelFont(42);
        hist->GetXaxis()->SetTitleFont(42);
        hist->GetYaxis()->SetTitleFont(42);
        hist->GetYaxis()->SetTitleOffset(1.7);
        hist->GetXaxis()->SetTitleOffset(1.25);
        if ((name.Contains("pT") || name.Contains("Mass")) && !name.Contains("Rapidity")) {
            hist->GetXaxis()->SetTitle(XAxis+" #left[GeV#right]");
            hist->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d"+XAxis+"}"+" #left[GeV^{-1}#right]"); 
        } else {
            hist->GetXaxis()->SetTitle(XAxis);
            hist->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d"+XAxis+"}");     
        }
        if (isControlPlot) hist->GetYaxis()->SetTitle(YAxis);
    }
}


void Plotter::PlotXSec(){

  TH1::AddDirectory(kFALSE);
  
  //CalcXSec(dataset, InclusiveXsection, InclusiveXsectionStatError, "","");

  double syst_square=0;

  //Here we'll just pull the values from the txt file like it was done for the diff case.

  for(int i =0; i<15; i++){
    //    syst_square += InclusiveXsectionSysErrorBySyst[channelType][i]*InclusiveXsectionSysErrorBySyst[channelType][i];
  }
  //  InclusiveXsectionSysError[channelType] = sqrt(syst_square);
  //cout<<"&^&^&^&^&^&^^&^&^ InclusiveXsectionSysError[channelType]: "<<InclusiveXsectionSysError[channelType]<<endl;

  double InclusiveXsectionStatErrorPlot[4], InclusiveXsectionSysErrorPlot[4], InclusiveXsectionPlot[4], InclusiveXsectionTotalErrorPlot[4];;

  if(channelType==3){

     // measured results with statistical error
   Double_t mx[]   = {      0.50,       1.50,       2.50,       3.50};
   Double_t mexl[] = {      0.00,       0.00,       0.00,       0.00};
   Double_t mexh[] = {      0.00,       0.00,       0.00,       0.00};

   TGraphAsymmErrors *mplot = new TGraphAsymmErrors(4, mx, InclusiveXsectionPlot, mexl, mexh,InclusiveXsectionStatErrorPlot, InclusiveXsectionStatErrorPlot);
   mplot->SetMarkerStyle(20);
   mplot->SetMarkerColor(kBlack);
   mplot->SetMarkerSize(1.5);
   mplot->SetLineColor(kBlack);
   
   for(int i=0; i<4; i++){
     InclusiveXsectionTotalErrorPlot[i] = sqrt(InclusiveXsectionStatErrorPlot[i]*InclusiveXsectionStatErrorPlot[i] +InclusiveXsectionPlot[i]*InclusiveXsectionSysErrorPlot[i]*InclusiveXsectionPlot[i]*InclusiveXsectionSysErrorPlot[i]);
   }

   TGraphAsymmErrors *mplotwithsys = new TGraphAsymmErrors(4, mx, InclusiveXsectionPlot, mexl, mexh,InclusiveXsectionTotalErrorPlot, InclusiveXsectionTotalErrorPlot);
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

   cout<<"!!!!!!!!!!!!!!!!!!!!ee Cross Section: "<<InclusiveXsectionPlot[0]<<" +/- "<<InclusiveXsectionStatErrorPlot[0]<<"(stat) +/- "<<InclusiveXsectionPlot[0]*InclusiveXsectionSysErrorPlot[0]<<"(sys)"<<endl;
   cout<<"!!!!!!!!!!!!!!!!!!!!mumu Cross Section: "<<InclusiveXsectionPlot[1]<<" +/- "<<InclusiveXsectionStatErrorPlot[1]<<"(stat) +/- "<<InclusiveXsectionPlot[0]*InclusiveXsectionSysErrorPlot[1]<<"(sys)"<<endl;
   cout<<"!!!!!!!!!!!!!!!!!!!!emu Cross Section: "<<InclusiveXsectionPlot[2]<<" +/- "<<InclusiveXsectionStatErrorPlot[2]<<"(stat) +/- "<<InclusiveXsectionPlot[0]*InclusiveXsectionSysErrorPlot[2]<<"(sys)"<<endl;
   cout<<"!!!!!!!!!!!!!!!!!!!!Combined Cross Section: "<<InclusiveXsectionPlot[3]<<" +/- "<<InclusiveXsectionStatErrorPlot[3]<<"(stat) +/- "<<InclusiveXsectionPlot[0]*InclusiveXsectionSysErrorPlot[3]<<"(sys)"<<endl;


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
    gSystem->MakeDirectory(outpathPlots+"/"+subfolderChannel);
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

  double BranchingFraction[4]={0.01166, 0.01166, 0.02332, 0.04666};//[ee, mumu, emu, combined] not including tau

  TH1D *numhists[hists.size()];
  double numbers[4]={0., 0., 0., 0.};//[0]=data, [1]=Signal, [2]Signal(only lumi & PU weights), [3]background (non-ttbar)
  double TTbarBGnum =0;

  for(unsigned int i=0; i<datasetVec.size(); i++){
    TFile *ftemp = TFile::Open(datasetVec[i]);
    TH1D *hist = (TH1D*)ftemp->Get("step9")->Clone();   

    double LumiWeight = CalcLumiWeight(datasetVec[i]);
    cout << "LumiWeight for " << dataset[i] << " = " << LumiWeight << endl;
    ApplyFlatWeights(hist, LumiWeight);
    
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
      TH1D *GenPlot = (TH1D*)ftemp->Get("GenAll")->Clone();  
      
      // Apply Scale Factor for MC@NLO
      ApplyMCATNLOWeight(GenPlot, Systematic, Shift,  datasetVec[i]); 
      
      numbers[2]+=GenPlot->Integral();  
      delete ftemp; 
      
    }  else if(legends[i] == "t#bar{t} Other"){

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
void Plotter::CalcDiffXSec(TString Channel, TString Systematic){

    
    cout<<"\n\n Preparing to calculate the central value of DiffXSection in channel '"<<Channel<<"' for systematic '"<<Systematic<<"'"<<endl;
    
    double Xbins[XAxisbins.size()];
    double binWidth[XAxisbinCenters.size()];
    for(unsigned int i = 0; i<XAxisbins.size();i++){Xbins[i]=XAxisbins[i];}
    double GenSignalSum[XAxisbinCenters.size()];

    // Getting the histogram
    TH1* theDataHist = NULL;
    TH1* theBgrHist = NULL;
    TH1* theTtBgrHist = NULL;
    TH1* theRecHist = NULL;
    //    TH1* theGenHist = GenPlot;
    //TH1* theRespHist = genReco2d;
    TH1* theGenHist = NULL; 
    TH1* theRespHist = NULL;


    TFile *ftemp = TFile::Open("preunfolded/"+Systematic+"/"+Channel+"/"+name+"_UnfoldingHistos.root");
    theDataHist =  (TH1D*)ftemp->Get("aDataHist")->Clone();
    theBgrHist =  (TH1D*)ftemp->Get("aBgrHist")->Clone();
    theTtBgrHist =  (TH1D*)ftemp->Get("aTtBgrHist")->Clone();
    theRecHist =  (TH1D*)ftemp->Get("aRecHist")->Clone();
    theGenHist =  (TH1D*)ftemp->Get("aGenHist")->Clone();
    theRespHist =  (TH2D*)ftemp->Get("aRespHist")->Clone();
    delete ftemp;
    //  aDataHist = drawhists[0]

    for (Int_t bin=0; bin<bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
        GenSignalSum[bin] = theGenHist->Rebin(bins,"aDataHist",Xbins)->GetBinContent(bin+1);
    }


    double GenDiffXSecVec[4][bins];
    double DiffXSecVec[4][bins];
    double DiffXSecStatErrorVec[4][bins];

    if(Channel=="ee"){channelType=0;}
    if(Channel=="mumu"){channelType=1;}
    if(Channel=="emu"){channelType=2;}
    if(Channel=="combined"){channelType=3;}
    
    // DAVID 
    if ( doUnfolding == true && channelType !=3) {//do the unfolding only in the individual channels: ee, emu, mumu

        // SVD Helper Class
        DilepSVDFunctions mySVDFunctions;
        mySVDFunctions.SetOutputPath(outpath);

        // Binning
        double* theBins = Xbins;
        int numberBins = bins;

        // Names and Labels
        TString channelLabelStr(channelLabel[channelType]);
        TString theChannelName = Channel;
        TString theParticleName = "";
        if ( name.Contains("Lepton")  ) theParticleName = "Leptons";
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
            UnfoldingError[i] = unfoldedDistributionNormalized->GetBinError(i+2); //absolute bin error
            //UnfoldingResult[i] = unfoldedDistribution->GetBinContent(i+2);//account for extra row in SVD unfolding
            //UnfoldingError[i] = unfoldedDistribution->GetBinError(i+2);
        }

        SignalEventswithWeight=0;
        // CROSS SECTION CALCULATION
        for (Int_t i=0; i<bins; ++i) {
            SignalEventswithWeight+=GenSignalSum[i];
        }

        for (Int_t i=0; i<bins; ++i) {
            //      if(Channel!="combined"){
            binWidth[i] = Xbins[i+1]-Xbins[i];
            DiffXSecVec[channelType][i] = UnfoldingResult[i]/(binWidth[i]);
            DiffXSecStatErrorVec[channelType][i] = UnfoldingError[i]/(binWidth[i]); // absolute statistical error
            GenDiffXSecVec[channelType][i] = (GenSignalSum[i]*topxsec)/(SignalEventswithWeight*binWidth[i]);//DIRTY (signal*topxsec)/(total events*binwidth)

            if(name.Contains("Lepton")||name.Contains("Top")||name.Contains("BJet")){
                GenDiffXSecVec[channelType][i] = GenDiffXSecVec[channelType][i]/2.;
            }
        }
    }
    
    if (doUnfolding && channelType==3){//for 'combined' channel: do an statistical combination of the the 3 independent channels
        
        TString eefilename="UnfoldingResults/Nominal/ee/"+name+"Results.txt";
        TString emufilename="UnfoldingResults/Nominal/emu/"+name+"Results.txt";
        TString mumufilename="UnfoldingResults/Nominal/mumu/"+name+"Results.txt";
        
        //check the existence of the file
        if(gSystem->AccessPathName(eefilename) || gSystem->AccessPathName(emufilename) || gSystem->AccessPathName(mumufilename) ){
            cout<<"WARNING (in CalcDiffXSec)!!"<<endl;
            cout<<"One of the input files you use for the combined XSection measurement doesn't exist!! Exiting!!"<<endl;
            return;
        }
        
        ifstream ResultsEE(eefilename);
        ifstream ResultsEMu(emufilename);
        ifstream ResultsMuMu(mumufilename);
        double perChannelDiffXSecPlot[3][bins];      //perChannelDiffXSecPlot[channel][bin]
        double perChannelDiffXSecStatError[3][bins]; //perChannelDiffXSecStatError[channel][bin]
        double perChannelGenDiffXSec[3][bins];       //perChannelGenDiffXSec[channel][bin]
        TString Dummy="";
        for (Int_t bin=0; bin<bins; bin++){//Retrieve arrays for plotting
            ResultsEE>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecPlot[0][bin]>>Dummy>>perChannelDiffXSecStatError[0][bin]>>Dummy>>perChannelGenDiffXSec[0][bin];
            ResultsEMu>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecPlot[1][bin]>>Dummy>>perChannelDiffXSecStatError[1][bin]>>Dummy>>perChannelGenDiffXSec[1][bin];
            ResultsMuMu>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecPlot[2][bin]>>Dummy>>perChannelDiffXSecStatError[2][bin]>>Dummy>>perChannelGenDiffXSec[2][bin];
            //cout<<"DiffXSecPlot[bin]: "<<DiffXSecPlot[bin]<<endl;
            //cout<<"GenDiffXSecPlot[bin]: "<<GenDiffXSecPlot[bin]<<endl;
        }
        ResultsEE.close(); ResultsEMu.close(); ResultsMuMu.close();
        
        //for gen level distribution
        for (Int_t i=0; i<bins; ++i) {
            SignalEventswithWeight+=GenSignalSum[i];
        }
        
        //do the actual combined Diff.XSection calculation
        for(int i=0; i<bins; i++){
            binWidth[i] = Xbins[i+1]-Xbins[i];
            for(int j=0; j<3; j++){//check if any stat error is 0, in this case set their contribution to 0!!
                if(perChannelDiffXSecStatError[j][i] == 0){
                    perChannelDiffXSecStatError[j][i] = 1e100;
                    perChannelDiffXSecPlot[j][i] = 0;
                }
            }
            DiffXSecVec[channelType][i] =(perChannelDiffXSecPlot[0][i]/(perChannelDiffXSecStatError[0][i]*perChannelDiffXSecStatError[0][i])
                                         +perChannelDiffXSecPlot[1][i]/(perChannelDiffXSecStatError[1][i]*perChannelDiffXSecStatError[1][i])
                                         +perChannelDiffXSecPlot[2][i]/(perChannelDiffXSecStatError[2][i]*perChannelDiffXSecStatError[2][i]))/
                                         ((1/(perChannelDiffXSecStatError[0][i]*perChannelDiffXSecStatError[0][i]))+
                                          (1/(perChannelDiffXSecStatError[1][i]*perChannelDiffXSecStatError[1][i]))+
                                          (1/(perChannelDiffXSecStatError[2][i]*perChannelDiffXSecStatError[2][i])));
                        
            DiffXSecStatErrorVec[channelType][i]=1/(TMath::Sqrt((1/(perChannelDiffXSecStatError[0][i]*perChannelDiffXSecStatError[0][i]))+
                                                                (1/(perChannelDiffXSecStatError[1][i]*perChannelDiffXSecStatError[1][i]))+
                                                                (1/(perChannelDiffXSecStatError[2][i]*perChannelDiffXSecStatError[2][i]))));
            GenDiffXSecVec[channelType][i] = (GenSignalSum[i]*topxsec)/(SignalEventswithWeight*binWidth[i]);//DIRTY (signal*topxsec)/(total events*binwidth)
        }
    }
    
    ofstream ResultsFile, ResultsLateX;

    gSystem->MakeDirectory("UnfoldingResults");
    gSystem->MakeDirectory("UnfoldingResults/"+Systematic);
    gSystem->MakeDirectory("UnfoldingResults/"+Systematic+"/"+Channel);

    string ResultsFilestring = outpathResults.Data();
    ResultsFilestring.append(subfolderSpecial.Data());
    ResultsFilestring.append("/");
    ResultsFilestring.append(Systematic);
    ResultsFilestring.append("/");
    ResultsFilestring.append(Channel);
    ResultsFilestring.append("/");
    ResultsFilestring.append(name);
    ResultsFilestring.append("Results.txt");
    ResultsFile.open(ResultsFilestring.c_str());


    string ResultsFilestringLatex = outpathPlots.Data();
    ResultsFilestringLatex.append(subfolderChannel.Data());
    ResultsFilestringLatex.append(subfolderSpecial.Data());
    ResultsFilestringLatex.append("/");
    ResultsFilestringLatex.append(name);
    ResultsFilestringLatex.append("ResultsLaTeX.txt");
    ResultsLateX.open(ResultsFilestringLatex.c_str());
    ResultsLateX<<"Bin Center & Bin & 1/#sigma d#sigma/dX & stat(\%) & syst(\%) & total(\%)"<<endl;
    for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
        //cout<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[bin]<<" bin: "<<Xbins[bin]<<" to "<<Xbins[bin+1]<<" DiffXsec: "<<DiffXSecVec[channelType][bin]<<" StatError: "<<DiffXSecStatErrorVec[channelType][bin]<<" GenDiffXsec: "<<GenDiffXSecVec[channelType][bin]<<endl;
        ResultsFile<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[bin]<<" bin: "<<Xbins[bin]<<" to "<<Xbins[bin+1]<<" DiffXsec: "<<DiffXSecVec[channelType][bin]<<" StatError: "<<DiffXSecStatErrorVec[channelType][bin]<<" GenDiffXsec: "<<GenDiffXSecVec[channelType][bin]<<endl;
    }
    ResultsFile.close();
    ResultsLateX.close();
}

void Plotter::PlotDiffXSec(TString Channel){

    if(Channel=="ee"){channelType=0;}
    if(Channel=="mumu"){channelType=1;}
    if(Channel=="emu"){channelType=2;}
    if(Channel=="combined"){channelType=3;}
    
    subfolderChannel = Channel;
    subfolderChannel.Prepend("/");

    TH1::AddDirectory(kFALSE);
    TGaxis::SetMaxDigits(2);

    double Xbins[XAxisbins.size()];
    for(unsigned int i = 0; i<XAxisbins.size();i++){Xbins[i]=XAxisbins[i];}
    double binCenters[XAxisbinCenters.size()];
    for(unsigned int i = 0; i<XAxisbinCenters.size();i++){binCenters[i]=XAxisbinCenters[i];}

    TH1 *RecoPlot = new TH1D;
    TH1 *GenPlot =new TH1D;
    TH1 *GenPlotTheory =new TH1D;

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

    TFile *ftemp = TFile::Open("preunfolded/Nominal/"+Channel+"/"+name+"_UnfoldingHistos.root");

    //    theDataHist =  (TH1D*)ftemp->Get("aDataHist")->Clone();
    //theBgrHist =  (TH1D*)ftemp->Get("aBgrHist")->Clone();
    //theTtBgrHist =  (TH1D*)ftemp->Get("aTtBgrHist")->Clone();
    RecoPlot =  (TH1D*)ftemp->Get("aRecHist")->Clone();
    GenPlotTheory =  (TH1D*)ftemp->Get("aGenHist")->Clone();
    genReco2d =  (TH2D*)ftemp->Get("aRespHist")->Clone();
    delete ftemp;


    for (unsigned int i =0; i<hists.size(); i++){
      varhists[i]=hists[i].Rebin(bins,"varhists",Xbins);
      setStyle(varhists[i], i);
    }

    GenPlot = GenPlotTheory->Rebin(bins,"genplot",Xbins);

    THStack * stack=  new THStack("def", "def");
    TLegend *leg = getNewLegendpre();
    int legchange = 0;
    TH1 *varhistsPlotting[hists.size()];


    for(unsigned int i=0; i<hists.size() ; i++){ // prepare histos and leg
        setStyle(varhists[i], i);
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
                    }
                    else leg->AddEntry(varhistsPlotting[i], legends[i] ,"f");
                }
                else{varhistsPlotting[legchange]->Add(varhistsPlotting[i]);}
            }
        }
        else{ if(i==0) leg->AddEntry(varhistsPlotting[i], legends[i] ,"pe");}
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

        if (i>0) {
            genRecHist->SetBinContent(i,((TH2D*)genReco2d)->Integral( intbinsX[i-1],intbinsX[i]-1,intbinsY[i-1],intbinsY[i]-1));
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
                 << " is larger than number of generated events in that bin" << endl;
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
    //cESP->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/ESP_"+name+".C");
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
                //cout<<"efficiencies[bin]: "<<efficiencies[bin]<<endl;
            }
        }
        else{
            for (Int_t bin=0; bin<bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
                BGSum[bin]+=varhists[hist]->GetBinContent(bin+1);
            }
        }
    }

    double totalDataSum = 0;
    double GenDiffXSecPlot[XAxisbinCenters.size()];
    for (Int_t bin=0; bin<bins; ++bin) {
        totalDataSum+=DataSum[bin];
    }

    TH1 *h_DiffXSec    = (TH1D*)varhists[0]->Clone();   h_DiffXSec->Reset();
    TH1 *h_GenDiffXSec = (TH1D*)varhists[0]->Clone();   h_GenDiffXSec->Reset();

    //The systematic array is filled in the order in which the Stack is filled
    TString sys_array[] = {"MATCH", "MASS", "SCALE", "BTAG_ETA_", "BTAG_PT_", "TRIG_", "BG_", "DY_", "PU_", "JER", "JES"};//For the time being uintil all systematics are finalished
    vector<TString> vec_systematic (sys_array, sys_array + sizeof(sys_array)/sizeof(sys_array[0]));

    double DiffXSecPlot[XAxisbinCenters.size()];
    double DiffXSecStatErrorPlot[XAxisbinCenters.size()];
    double DiffXSecSysErrorbySysPlot[XAxisbinCenters.size()][(int)vec_systematic.size()];
    double DiffXSecSysErrorPlot[XAxisbinCenters.size()];
    double DiffXSecTotalErrorPlot[XAxisbinCenters.size()];

    double ModelSysPlot[XAxisbinCenters.size()];
    double ExpSysPlot[XAxisbinCenters.size()];

    //Read absolute systematic uncertainty of each bin from file
    for(int Syst=0; Syst<(int)vec_systematic.size(); Syst++){
        ifstream SysResultsList("UnfoldingResults/"+vec_systematic.at(Syst)+"/"+Channel+"/"+name+"Results.txt");
        for (Int_t bin=0; bin<bins; ++bin) {
            TString DUMMY;
            SysResultsList>>DUMMY>>XAxisbinCenters[bin]>>DUMMY>>Xbins[bin]>>DUMMY>>Xbins[bin+1]>>DUMMY>>DiffXSecSysErrorbySysPlot[bin][Syst];
        }
        SysResultsList.close();
    }

    TString Dummy;
    //Read central and absolute statistical uncertainty values from Nominal
    ifstream ResultsList("UnfoldingResults/Nominal/"+Channel+"/"+name+"Results.txt");
    for (Int_t bin=0; bin<bins; bin++){//Retrieve arrays for plotting
        ResultsList>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>DiffXSecPlot[bin]>>Dummy>>DiffXSecStatErrorPlot[bin]>>Dummy>>GenDiffXSecPlot[bin];
        //cout<<"DiffXSecPlot[bin]: "<<DiffXSecPlot[bin]<<endl;
        //cout<<"GenDiffXSecPlot[bin]: "<<GenDiffXSecPlot[bin]<<endl;
    }

    for (Int_t i=0; i<bins; ++i) {
        h_DiffXSec->SetBinContent(i+1,DiffXSecPlot[i]);
        h_DiffXSec->SetBinError(i+1,DiffXSecStatErrorPlot[i]);
        h_GenDiffXSec->SetBinContent(i+1,GenDiffXSecPlot[i]);
    }

    double TotalVisXSection = 1.; //this can currently be set to 1. because David's code takes care of the normalization, but just incase we need it

    h_DiffXSec->Scale(1/TotalVisXSection);

    for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
        double syst_square = 0;
        ExpSysPlot[bin]=0.;
        ModelSysPlot[bin]=0.;
        if(doSystematics){
            for(int syst =0; syst<(int)vec_systematic.size(); syst++){
                syst_square += DiffXSecSysErrorbySysPlot[bin][syst]*DiffXSecSysErrorbySysPlot[bin][syst];
                if(vec_systematic.at(syst).Contains("JER") || vec_systematic.at(syst).Contains("JES") || vec_systematic.at(syst).Contains("PU_") ||
                   vec_systematic.at(syst).Contains("DY_") || vec_systematic.at(syst).Contains("BG_") || vec_systematic.at(syst).Contains("TRIG") ||
                   vec_systematic.at(syst).Contains("lepton") || vec_systematic.at(syst).Contains("BTAG_") || vec_systematic.at(syst).Contains("kin fit")){
                    ExpSysPlot[bin]+=DiffXSecSysErrorbySysPlot[bin][syst]*DiffXSecSysErrorbySysPlot[bin][syst];
                }
                else{
                    ModelSysPlot[bin]+=DiffXSecSysErrorbySysPlot[bin][syst]*DiffXSecSysErrorbySysPlot[bin][syst];
                }
            }
            DiffXSecSysErrorPlot[bin]+=syst_square;
        }
        ExpSysPlot[bin]=sqrt(ExpSysPlot[bin]);
        ModelSysPlot[bin]=sqrt(ModelSysPlot[bin]);
        DiffXSecStatErrorPlot[bin] = DiffXSecStatErrorPlot[bin]/TotalVisXSection;
        DiffXSecPlot[bin]=DiffXSecPlot[bin]/TotalVisXSection;
        DiffXSecSysErrorPlot[bin]=sqrt(DiffXSecSysErrorPlot[bin])*DiffXSecPlot[bin]; //absolute systematic error in bin 'bin'
        DiffXSecTotalErrorPlot[bin]=sqrt(DiffXSecSysErrorPlot[bin]*DiffXSecSysErrorPlot[bin] + DiffXSecStatErrorPlot[bin]*DiffXSecStatErrorPlot[bin]);//absolute total error
    }

    //create a file for Results!!
    //Right now this is a check to make sure what is calculated in the unfolding step is still what we see!
    ofstream ResultsFile, ResultsLateX;
    string ResultsFilestring = outpathPlots.Data();
    ResultsFilestring.append(subfolderChannel.Data());
    ResultsFilestring.append(subfolderSpecial.Data());
    ResultsFilestring.append("/");
    ResultsFilestring.append(newname);
    ResultsFilestring.append("ResultsAfter.txt");
    ResultsFile.open(ResultsFilestring.c_str());
    
    string ResultsFilestringLatex = outpathPlots.Data();
    ResultsFilestringLatex.append(subfolderChannel.Data());
    ResultsFilestringLatex.append(subfolderSpecial.Data());
    ResultsFilestringLatex.append("/");
    ResultsFilestringLatex.append(newname);
    ResultsFilestringLatex.append("ResultsLaTeXAfter.txt");
    ResultsLateX.open(ResultsFilestringLatex.c_str());
    ResultsLateX<<"Bin Center & Bin & 1/#sigma d#sigma/dX & stat(\%) & syst(\%) & total(\%)"<<endl;

    for (Int_t bin=0; bin<bins; bin++){
      ResultsFile<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[bin]<<" bin: "<<Xbins[bin]<<" to "<<Xbins[bin+1]<<" DiffXsec: "<<DiffXSecPlot[bin]<<" StatError(percent): "<<DiffXSecStatErrorPlot[bin]/DiffXSecPlot[bin]<<" SysError: "<<DiffXSecSysErrorPlot[bin]<<" TotalError: "<<DiffXSecTotalErrorPlot[bin]/DiffXSecPlot[bin]<<" GenDiffXSec: "<<GenDiffXSecPlot[bin]<<endl;
      ResultsLateX<<"$"<<h_DiffXSec->GetBinCenter(bin+1)<<"$ & $" <<h_DiffXSec->GetBinLowEdge(bin+1)<<"$ to $"<<h_DiffXSec->GetBinLowEdge(bin+2)<<"$ & ";
      ResultsLateX<<DiffXSecPlot[bin]<<" & "<<setprecision(3)<<DiffXSecStatErrorPlot[bin]*100./DiffXSecPlot[bin]<<" & "<<setprecision(3)<<100.*DiffXSecSysErrorPlot[bin]<<" & "<<setprecision(3)<<100.*DiffXSecTotalErrorPlot[bin]/DiffXSecPlot[bin]<< "\\\\" << endl;
    }
    ResultsFile.close();
    ResultsLateX.close();
    
    if(doSystematics){

        //The Markus plots
        TCanvas * c10 = new TCanvas("Markus","Markus");
        THStack* SystHists = new THStack("MSTACK","MSTACK");
        TLegend * leg10 =  new TLegend(0.20,0.65,0.45,0.90);

        ofstream ResultsSysFilestring; 
        string ResultsSystLaTeX = outpathPlots.Data();
        ResultsSystLaTeX.append(subfolderChannel.Data());
        ResultsSystLaTeX.append(subfolderSpecial.Data());
        ResultsSystLaTeX.append("/");
        ResultsSystLaTeX.append(newname);
        ResultsSystLaTeX.append("SystematicsLaTeX.txt");
        FILE *systfile;
        systfile = fopen(ResultsSystLaTeX.c_str(), "w");
        
        for(int systs =0; systs<(int)vec_systematic.size(); systs++){
            if (vec_systematic.at(systs) == "BTAG_ETA_") {continue;}//Skip the BTAG_ETA systematic because it's added in quadrature to BTAG_PT
            TH1D* systtemp = (TH1D*)varhists[0]->Clone();
            systtemp->Reset();
            double TotalSyst=0.0, TotalSqSyst=0.0;
            double AvgSyst= 0.0, SqAvgSys=0.0;
            
            for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
                if(vec_systematic.at(systs) == "BTAG_PT_"){
                    DiffXSecSysErrorbySysPlot[bin][systs]= TMath::Sqrt(
                        (DiffXSecSysErrorbySysPlot[bin][systs]*DiffXSecSysErrorbySysPlot[bin][systs])
                        +(DiffXSecSysErrorbySysPlot[bin][systs+1]*DiffXSecSysErrorbySysPlot[bin][systs+1])
                    );
                }

                systtemp->SetBinContent(bin+1,(DiffXSecSysErrorbySysPlot[bin][systs]*DiffXSecSysErrorbySysPlot[bin][systs]));
                if(bin==0){
                    fprintf(systfile, vec_systematic.at(systs)+" ");
                }
                fprintf(systfile, "%2.5f ", TMath::Sqrt(systtemp->GetBinContent(bin+1))*100);
                if(bin>0 && bin<bins-1){//Exclude the 2 side bins
                    TotalSyst=TotalSyst+TMath::Sqrt(systtemp->GetBinContent(bin+1));
                    TotalSqSyst=TotalSqSyst+systtemp->GetBinContent(bin+1);
                }
            }
            AvgSyst=TotalSyst/(bins-2);
            SqAvgSys=TMath::Sqrt(TotalSqSyst/(bins-2));
            fprintf(systfile, "Lin.Avg.(%%)= %.5f  Quad.Avg.(%%)=%.5f\n", 100*AvgSyst, 100*SqAvgSys);
            systtemp->SetFillColor((int)vec_systematic.size()-systs);
            SystHists->Add((TH1D*)systtemp->Clone());
            leg10->AddEntry(systtemp->Clone(), vec_systematic.at(systs), "f");
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
        //c10->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/MSP_"+name+".C");
        c10->Clear();
        delete leg10;
        delete c10;
        
        //The Experimental/Model/Statistical plot
        TCanvas * c11 = new TCanvas("EMS","EMS");
        TH1D* ExpHist = (TH1D*)varhists[0]->Clone();    ExpHist->Reset();
        TH1D* ModelHist = (TH1D*)varhists[0]->Clone();  ModelHist->Reset();
        TH1D* StatHist = (TH1D*)varhists[0]->Clone();   StatHist->Reset();
        TH1D* TotalHist = (TH1D*)varhists[0]->Clone();  TotalHist->Reset();
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
            ExpHist->SetBinContent(bin+1,100*ExpSysPlot[bin]);
            ModelHist->SetBinContent(bin+1,100*ModelSysPlot[bin]);
            StatHist->SetBinContent(bin+1,100*DiffXSecStatErrorPlot[bin]/DiffXSecPlot[bin]);
            TotalHist->SetBinContent(bin+1,100*DiffXSecTotalErrorPlot[bin]/DiffXSecPlot[bin]);
        }
        TotalHist->SetMinimum(0.);
        TotalHist->GetYaxis()->SetTitle("#frac{#Delta#sigma}{#sigma} [%]");
        TotalHist->SetLineColor(1);
        ExpHist->SetLineColor(kRed);
        StatHist->SetLineColor(kGreen);
        ModelHist->SetLineColor(kBlue);
        TLegend * leg11 =  new TLegend(0.65,0.60,0.90,0.85);
        leg11->SetFillColor(0);
        leg11->AddEntry(ExpHist->Clone(), "Experimental Uncertainty", "l");
        leg11->AddEntry(StatHist->Clone(), "Statistical Uncertainty", "l");
        leg11->AddEntry(ModelHist->Clone(), "Model Uncertainty", "l");
        leg11->AddEntry(TotalHist->Clone(), "Total Uncertainty", "l");
        TotalHist->Draw();ModelHist->Draw("SAME");ExpHist->Draw("SAME");StatHist->Draw("SAME");
        leg11->Draw("SAME");
        TotalHist->GetXaxis()->SetNoExponent(kTRUE);
        c11->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/SEM_"+name+".eps");
        //c11->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/SEM_"+name+".C");
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

    GenPlotTheory->Scale(topxsec/(SignalEventswithWeight*GenPlotTheory->GetBinWidth(1)));
    GenPlot->Scale(topxsec/(SignalEventswithWeight*GenPlot->GetBinWidth(1)));
    if(name.Contains("Lepton")||name.Contains("Top")||name.Contains("BJet")){
      GenPlotTheory->Scale(1./2.);
    }
    double genscale = 1./GenPlotTheory->Integral("width");
    GenPlotTheory->Scale(genscale);
    
    genscale = 1./ h_GenDiffXSec->Integral("width");
    h_GenDiffXSec->Scale(genscale);
    bool binned_theory=true; //############
  
    
  
    TH1* mcnlohist=0, *mcnlohistup=0, *mcnlohistdown=0, *powheghist=0, *spincorrhist=0;
    TH1* mcnlohistnorm=0;
    TGraph *mcatnloBand=0;
    
    TH1* mcnlohistnormBinned = 0, *mcnlohistupBinned = 0;
    TH1 *mcnlohistdownBinned = 0, *mcnlohistBinned = 0;
    TH1* powheghistBinned = 0;
    TH1* spincorrhistBinned = 0;
    
    TH1F *Kidoth1_Binned = 0;
    TFile *KidoFile = 0;
    
    if (drawNLOCurves && drawMCATNLO) {
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
        else {mcnlohistnorm = new TH1D();}
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
        else {mcnlohistup = new TH1D();}
        //    if (binned_theory==false) mcnlohistup->Rebin(5);mcnlohistup->Scale(0.2);
        mcnlohistupBinned    = mcnlohistup->Rebin(bins,"genBinHist", Xbins);


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
        else {mcnlohistdown = new TH1D();}
        //    if (binned_theory==false) mcnlohistdown->Rebin(5);mcnlohistdown->Scale(0.2);
        mcnlohistdownBinned    = mcnlohistdown->Rebin(bins,"genBinHist", Xbins);

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
        
    }
    if(drawNLOCurves && drawPOWHEG){
        powheghist = GetNloCurve(newname, "POWHEG");
        double powhegscale = 1./powheghist->Integral("width");
        if (binned_theory==false) powheghist->Rebin(2);powheghist->Scale(0.5);
        powheghist->Scale(powhegscale);
            
        powheghistBinned = powheghist->Rebin(bins,"powhegplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){
            powheghistBinned->SetBinContent(bin+1,powheghistBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/powheghist->GetBinWidth(1)));
        }
        powheghistBinned->Scale(1./powheghistBinned->Integral("width"));
    }
    
    if(drawNLOCurves && drawMadSpinCorr){
        spincorrhist = GetNloCurve(newname, "SPINCORR");
        double spincorrhistscale = 1./spincorrhist->Integral("width");
        if (binned_theory==false) spincorrhist->Rebin(2);spincorrhist->Scale(0.5);
        spincorrhist->Scale(spincorrhistscale);
            
        spincorrhistBinned = spincorrhist->Rebin(bins,"spincorrplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){
            spincorrhistBinned->SetBinContent(bin+1,spincorrhistBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/spincorrhist->GetBinWidth(1)));
        }
        spincorrhistBinned->Scale(1./spincorrhistBinned->Integral("width"));
    }
    
    
    
    if(drawNLOCurves && drawKidonakis && (name.Contains("ToppT") || name.Contains("TopRapidity"))){
        KidoFile=TFile::Open("dilepton_kidonakisNNLO.root");
        if(name.Contains("ToppT")){
            Kidoth1_Binned = (TH1F*)KidoFile->Get("topPt");
        }
        else if(name.Contains("TopRapidity")){
            Kidoth1_Binned = (TH1F*)KidoFile->Get("topY");
        }
    }

    if(drawNLOCurves){
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
    if (drawNLOCurves && drawMCATNLO) {
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

    h_GenDiffXSec->SetLineColor(kRed+1);
    h_GenDiffXSec->SetLineStyle(1);

    if (drawNLOCurves && drawMCATNLO) {
        mcnlohist->SetLineColor(kBlue); //#####################
        mcnlohist->SetLineStyle(5);
        mcnlohistBinned->SetLineColor(kBlue); //#####################
        mcnlohistBinned->SetLineWidth(2);
        mcnlohistBinned->SetLineStyle(5);
        if(binned_theory==false){mcnlohist->Draw("SAME,C");}
        else{mcnlohistBinned->Draw("SAME");}
    }
    if(drawNLOCurves && drawPOWHEG){
        powheghist->SetLineColor(kGreen+1); //#####################
        powheghist->SetLineStyle(7);
        powheghistBinned->SetLineColor(kGreen+1); //#####################
        powheghistBinned->SetLineWidth(2);
        powheghistBinned->SetLineStyle(7);
        
        if(binned_theory==false){powheghist->Draw("SAME,C");}
        else{powheghistBinned->Draw("SAME");}
    }
    if(drawNLOCurves && drawMadSpinCorr){
        spincorrhist->SetLineColor(kOrange+3); //#####################
        spincorrhist->SetLineStyle(2);
        spincorrhistBinned->SetLineColor(kOrange+3); //#####################
        spincorrhistBinned->SetLineWidth(2);
        spincorrhistBinned->SetLineStyle(2);
        
        if(binned_theory==false){spincorrhist->Draw("SAME,C");}
        else{spincorrhistBinned->Draw("SAME");}
    }
    
    if(drawNLOCurves && drawKidonakis){
        if(!name.Contains("ToppT") && !name.Contains("TopRapidity")){}
        Kidoth1_Binned->SetLineWidth(2);
        Kidoth1_Binned->SetLineColor(kOrange-3); //########################
        Kidoth1_Binned->SetLineStyle(1);
        Kidoth1_Binned->Draw("SAME");
    }
    /*
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
*/
    h_GenDiffXSec->Draw("SAME"); //### 150512 ###

    DrawCMSLabels(2, lumi, 8);
    
    DrawDecayChLabel(channelLabel[channelType]);
    
    TLegend leg2 = *getNewLegend();
    leg2.AddEntry(h_DiffXSec, "Data",    "p");
    leg2.AddEntry(GenPlotTheory,            "MadGraph","l");
    if (drawNLOCurves) {
        if (drawMCATNLO && mcnlohistup->GetEntries() && mcnlohistdown->GetEntries())   leg2.AddEntry(mcatnloBand,      "MC@NLO",  "fl");
        if (drawPOWHEG && powheghist->GetEntries())                                    leg2.AddEntry(powheghistBinned, "POWHEG",  "l");
        if (drawMadSpinCorr && spincorrhist->GetEntries())                             leg2.AddEntry(spincorrhistBinned, "MadGraph SC",  "l");
        if (drawKidonakis && (name.Contains("ToppT") || name.Contains("TopRapidity"))) leg2.AddEntry(Kidoth1_Binned,   "Approx. NNLO",  "l");
    }
    
    leg2.SetFillStyle(0);
    leg2.SetBorderSize(0);
    leg2.SetX1NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength()-0.30);
    leg2.SetY1NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength()-0.05*(double)leg2.GetNRows());
    leg2.SetX2NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength());
    leg2.SetY2NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength());
    leg2.SetTextSize(0.04);
    leg2.Draw("same");
     if (drawNLOCurves && drawKidonakis &&  (name.Contains("ToppT") || name.Contains("TopRapidity"))){
        DrawLabel("(arXiv:1210.7813v1)", leg2.GetX1NDC()+0.06, leg2.GetY1NDC()-0.025, leg2.GetX2NDC(), leg2.GetY1NDC(), 12, 0.025);
    }
    
    h_GenDiffXSec->Draw("SAME");
    gStyle->SetEndErrorSize(10);
    tga_DiffXSecPlot->Draw("p, SAME");
    tga_DiffXSecPlotwithSys->Draw("p, SAME, Z");
    gPad->RedrawAxis();
    
    c->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/DiffXS_"+name+".eps"); 
    //c->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/DiffXS_"+name+".C"); 
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
    //leg->AddEntry( syshist, "Uncertainty", "f" );


    varhists[0]->SetMaximum(1.5*varhists[0]->GetBinContent(varhists[0]->GetMaximumBin()));

    varhists[0]->SetMinimum(0);
    varhists[0]->GetYaxis()->SetTitle("events");
    varhists[0]->GetXaxis()->SetNoExponent(kTRUE);
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

    //Only necessary if we want error bands

    /*    TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0.5)");//this is frustrating and stupid but apparently necessary...
    setex1->Draw();
    syshist->Draw("same,E2");
    TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.)");
    setex2->Draw();*/
    varhists[0]->Draw("same, e1"); //############
    //varhists[0]->Draw("same, e"); 
    DrawCMSLabels(2, lumi, 8);
    DrawDecayChLabel(channelLabel[channelType]);    
    leg->Draw("SAME");
    gPad->RedrawAxis();
    //    TFile *f1 = new TFile("KinFitPlots.root","UPDATE");
    //stacksum->Write(name+"_"+channel+"_MC");
    //varhists[0]->Write(name+"_"+channel+"_Data");
    //f1->Close();
    c1->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/preunfolded_"+name+".eps");
    //c1->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/preunfolded_"+name+".C");
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
  
  TFile* file = 0;
  
  if(strcmp(generator, "Powheg")==0){file = TFile::Open("selectionRoot/Nominal/emu/ttbarsignalplustau_powheg.root","READ");}
  else if(strcmp(generator, "MCatNLO")==0){file = TFile::Open("MCatNLO_status3_v20120729.root","READ");}
  else if(strcmp(generator, "MCNLOup")==0){file = TFile::Open("MCatNLO_Uncert_Up_status3_v20120729.root","READ");}
  else if(strcmp(generator, "MCNLOdown")==0){file = TFile::Open("MCatNLO_Uncert_Down_status3_v20120729.root","READ");}
  
  if (file && !file->IsZombie()) {
    file->GetObject<TH1>(histname, hist);

    if(!hist){
      std::cerr << "WARNING in GetNloCurve: input histogram '" << histname << "' could not been opened! Returning dummy!" << endl;
      hist = new TH1D();
      return hist;
    }
    
    TH1D* rethist = (TH1D*)hist->Clone();
    TH1D* weight = (TH1D*)file->Get(TString("total_LeptonsPt_").Append(generator));
    
    Double_t wgt = 1.;
    if(!weight){
      std::cerr << "WARNING in GetNloCurve: histogram to extract original number of events could not be opened! No weighting applied!" << endl;
    }

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
    
    TFile* file = 0;
    TFile* file1 = 0;
    TFile* file2 = 0;

        
    if(Generator=="MCATNLO"){
        if(channelType == 0)file = TFile::Open("selectionRoot/"+Generator+"/ee/ttbarsignalplustau_mcatnlo.root","READ");
        else if(channelType == 1)file = TFile::Open("selectionRoot/"+Generator+"/mumu/ttbarsignalplustau_mcatnlo.root","READ");
        else if(channelType == 2)file = TFile::Open("selectionRoot/"+Generator+"/emu/ttbarsignalplustau_mcatnlo.root","READ");
        else {
            file = TFile::Open("selectionRoot/"+Generator+"/emu/ttbarsignalplustau_mcatnlo.root","READ");
            file1 = TFile::Open("selectionRoot/"+Generator+"/ee/ttbarsignalplustau_mcatnlo.root","READ");
            file2 = TFile::Open("selectionRoot/"+Generator+"/mumu/ttbarsignalplustau_mcatnlo.root","READ");
        }
    }
    else if(Generator="POWHEG"){
        if(channelType == 0)file = TFile::Open("selectionRoot/"+Generator+"/ee/ee_ttbarsignalplustau_powheg.root","READ");
        else if(channelType == 1)file = TFile::Open("selectionRoot/"+Generator+"/mumu/mumu_ttbarsignalplustau_powheg.root","READ");
        else if(channelType == 2)file = TFile::Open("selectionRoot/"+Generator+"/emu/emu_ttbarsignalplustau_powheg.root","READ");
        else {
            file = TFile::Open("selectionRoot/"+Generator+"/emu/emu_ttbarsignalplustau_powheg.root","READ");
            file1 = TFile::Open("selectionRoot/"+Generator+"/ee/ee_ttbarsignalplustau_powheg.root","READ");
            file2 = TFile::Open("selectionRoot/"+Generator+"/mumu/mumu_ttbarsignalplustau_powheg.root","READ");
        }
    }
    else{
        if(channelType == 0)file = TFile::Open("selectionRoot/"+Generator+"/ee/ee_ttbarsignalplustau_ttbarsignalplustau_FullLeptMadgraphWithSpinCorrelation.root","READ");
        else if(channelType == 1)file = TFile::Open("selectionRoot/"+Generator+"/mumu/mumu_ttbarsignalplustau_ttbarsignalplustau_FullLeptMadgraphWithSpinCorrelation.root","READ");
        else if(channelType == 2)file = TFile::Open("selectionRoot/"+Generator+"/emu/emu_ttbarsignalplustau_ttbarsignalplustau_FullLeptMadgraphWithSpinCorrelation.root","READ");
        else {
            file = TFile::Open("selectionRoot/"+Generator+"/emu/emu_ttbarsignalplustau_ttbarsignalplustau_FullLeptMadgraphWithSpinCorrelation.root","READ");
            file1 = TFile::Open("selectionRoot/"+Generator+"/ee/ee_ttbarsignalplustau_ttbarsignalplustau_FullLeptMadgraphWithSpinCorrelation.root","READ");
            file2 = TFile::Open("selectionRoot/"+Generator+"/mumu/mumu_ttbarsignalplustau_ttbarsignalplustau_FullLeptMadgraphWithSpinCorrelation.root","READ");
        }

    }

    if (file && !file->IsZombie()) {
        if (channelType<3)hist=(TH1*)file->Get("VisGen"+NewName)->Clone();
        else {
        hist=(TH1*)file->Get("VisGen"+NewName)->Clone();
        hist->Add((TH1*)file1->Get("VisGen"+NewName)->Clone());
        hist->Add((TH1*)file2->Get("VisGen"+NewName)->Clone());
        }
        if(!NewName.Contains("Lead") && (NewName.Contains("Lepton")||NewName.Contains("Top")||NewName.Contains("BJet"))){
            if(channelType<3)hist->Add((TH1*)file->Get("VisGenAnti"+NewName)->Clone());
            else{
                hist->Add((TH1*)file->Get("VisGenAnti"+NewName)->Clone());
                hist->Add((TH1*)file1->Get("VisGenAnti"+NewName)->Clone());
                hist->Add((TH1*)file2->Get("VisGenAnti"+NewName)->Clone());
            }
        }
        if(!hist){
            std::cerr << "WARNING in GetNloCurve: input histogram '" << histname << "' could not been opened! Returning dummy!" << endl;
            hist = new TH1D();
            return hist;
        }

        TH1D* rethist = (TH1D*)hist->Clone();

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



void Plotter::ApplyFlatWeights(TH1* varhists, const double weight){

    if(weight == 0) {cout<<"Warning: the weight your applying is 0. This will remove your distribution."<<endl;}
    //if(weight >=1e3){cout<<"Warning: the weight your applying is >= 1e3. This will enlarge too much your distribution."<<endl;}
    varhists->Scale(weight);
}


void Plotter::ApplyFlatWeights(TH1* varhists[], const double weight){

    if(weight == 0) {cout<<"Warning: the weight your applying is 0. This will remove your distribution."<<endl;}
    if(weight >=1e3){cout<<"Warning: the weight your applying is >= 1e3. This will enlarge too much your distribution."<<endl;}
    for(size_t i=0; i<hists.size(); i++){
        varhists[i]->Scale(weight);
    }
}

double Plotter::CalcLumiWeight(TString WhichSample){
    if (WhichSample.Contains("run")) return 1;
    double lumiWeight=0;
    if(WhichSample!=""){
        double XSection = SampleXSection(WhichSample);
        if(XSection <= 0.){
            cout<<"Sample XSection is <0. Can't calculate luminosity weight!! returning"<<endl;
            return 0;
        };
        
        //From 'filename' get the number of weighted (MC weights) event processed.
        TFile *f = TFile::Open(WhichSample);
        if(!f || f->IsZombie()){
            cout<<"The file you requested is not valid"<<endl;
            return 0;
        };
        TH1 *h_NrOfEvts = dynamic_cast<TH1*>(f->Get("weightedEvents"));
        if(!h_NrOfEvts){
            cout<<"The histogram you requested is not valid"<<endl;
            return 0;
        };
        double NrOfEvts = h_NrOfEvts->GetBinContent(1);
       
        lumiWeight = lumi*XSection/NrOfEvts;

        f->Close();
        h_NrOfEvts->Delete();
    }
    
    if (lumiWeight == 0) {
        cout << WhichSample << " has lumi weight 0\n";
    }

    return lumiWeight;
}


double Plotter::SampleXSection(TString filename){
    
    //MC cross sections taken from:
    //  https://twiki.cern.ch/twiki/bin/view/CMS/StandardModelCrossSectionsat8TeV
    //  AN-12/194    AN-12/228
    
    double XSec=-1.;
    if(filename.Contains("run"))              {XSec = 1;}
    else if(filename.Contains("ttbar"))       {XSec = 213.192;}//225.197;}
    else if(filename.Contains("single"))      {XSec = 22.2;}
    else if(filename.Contains("ww"))          {XSec = 54.838;}
    else if(filename.Contains("wz"))          {XSec = 33.21;}
    else if(filename.Contains("zz"))          {XSec = 17.654;}
    else if(filename.Contains("1050"))        {XSec = 860.5;}
    else if(filename.Contains("50inf"))       {XSec = 3532.8;}
    else if(filename.Contains("wtolnu"))      {XSec = 36257.2;}
    else if(filename.Contains("qcdmu15"))     {XSec = 3.640E8*3.7E-4;}
    else if(filename.Contains("qcdmu2030"))   {XSec = 2.870E8*6.500E-3;}
    else if(filename.Contains("qcdmu3050"))   {XSec = 6.609E7*12.20E-3;}
    else if(filename.Contains("qcdmu5080"))   {XSec = 8.802E6*21.80E-3;}
    else if(filename.Contains("qcdmu80120"))  {XSec = 1.024E6*39.50E-3;}
    else if(filename.Contains("qcdmu120170")) {XSec = 1.578E5*47.30E-3;}
    else if(filename.Contains("qcdem2030"))   {XSec = 2.886E8*10.10E-3;}
    else if(filename.Contains("qcdem3080"))   {XSec = 7.433E7*62.10E-3;}
    else if(filename.Contains("qcdem80170"))  {XSec = 1.191E6*153.9E-3;}
    else if(filename.Contains("qcdbcem2030")) {XSec = 2.886E8*5.800E-4;}
    else if(filename.Contains("qcdbcem3080")) {XSec = 7.424E7*2.250E-3;}
    else if(filename.Contains("qcdbcem80170")){XSec = 1.191E6*10.90E-3;}
    
    return XSec;
}





#endif
