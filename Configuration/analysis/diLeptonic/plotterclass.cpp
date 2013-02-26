#include "plotterclass.h"

#include <fstream>
#include <iostream>
#include <cstdio>
#include <sstream>
#include <cmath>

#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TExec.h>
#include <TStyle.h>
#include <TMath.h>
#include <TROOT.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>

#include "utils.h"
  
  // DAVID
#include "DilepSVDFunctions.h"

using namespace std;

const double Plotter::topxsec = 244.849; //again changes with normalization, must be set outside of the class
const bool Plotter::doClosureTest = 0; //Signal is MC - so add BG to it and dont do DY scaling

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

  drawSmoothMadgraph = true;
  drawMadSpinCorr  = false;
  drawMCATNLO      = true;
  drawKidonakis    = true;
  drawPOWHEG       = true;

}


// DAVID
void Plotter::SetOutpath(TString path)
{
  outpath = path;
}


void Plotter::unfolding()
{
    TString sys_array[] = {"HAD_", "LEPT_", "KIN_", "DY_","BG_","PU_", "TRIG_","MASS_", "MATCH_", "SCALE_", "BTAG_ETA_","BTAG_PT_", "BTAG_LJET_ETA_", "BTAG_LJET_PT_", "JER_", "JES_"};

    double sys_array_flat_value[]={0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0.0025, 0.0025, 0.0025, 0.0025, 0, 0};
    //Right now BTAG = 0.5% = Sqrt(BTAG_ETA ** 2 + BTAG_PT ** 2 + BTAG_LJET_ETA ** 2 + BTAG_LJET_PT ** 2)
    //assumed  BTAG_ETA = BTAG_PT = BTAG_LJET_ETA = BTAG_LJET_PT = 0.25% so a total error of 0.5% comes up in BTag
    
    
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

        //############### Syst ################
        if (doSystematics) { 
            for(int sys = 0; sys < (int)vec_systematic.size(); sys++){ //loop over systematics
                cout << endl;
                cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
                cout << "Starting Calculation of Differential Systematics for '" << name << "' in Channel '" << vec_channel.at(chan) << "':" << endl;
                cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;

                if ( !vec_systematic.at(sys).BeginsWith("HAD_") ) {
                     CalcDiffSystematics(vec_channel.at(chan), 
                                        vec_systematic.at(sys), 
                                        vec_systematic.at(sys)+"UP",
                                        vec_systematic.at(sys)+"DOWN",
                                        vec_flat_value.at(sys));
                } else {
                    if(vec_channel.at(chan) != "combined"){
		        write(vec_channel.at(chan),"POWHEG");//make sure the preunfolded files exist
		        write(vec_channel.at(chan),"MCATNLO");
                        CalcDiffXSec(vec_channel.at(chan), "POWHEG");
                        CalcDiffXSec(vec_channel.at(chan), "MCATNLO");
                        GetDiffToNominal(vec_channel.at(chan), TString("POWHEG"), name);
                        GetDiffToNominal(vec_channel.at(chan), TString("MCATNLO"), name);
                        CalcUpDownDifference(vec_channel.at(chan), "POWHEG", "MCATNLO", name);
                    }
                    else if (vec_channel.at(chan) == "combined"){
                        CalcDiffSystematics(vec_channel.at(chan), vec_systematic.at(sys), vec_systematic.at(sys), vec_systematic.at(sys),vec_flat_value.at(sys));
                    }
                }

                cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
                cout << "Finished Calculation of Differential Systematics for '" << name << "' in Channel '" << vec_channel.at(chan) << endl;
                cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
                cout << endl;
            }//systematic loop
        }
    }//channel loop

}

void Plotter::preunfolding(TString Channel, TString Systematic)
{
    write(Channel, Systematic);
}


void Plotter::DYScaleFactor(TString SpecialComment){

    DYScale = {1,1,1,1}; 

    if(!doDYScale || doClosureTest) return; //need to make a switch for control plots that don't want DYScale

    TString nameAppendix = "";
    if ( !SpecialComment.BeginsWith("_post") &&  SpecialComment != "Standard" ){
        cout<<"\n\n*******************************************************************"<<endl;
        cout<<"ERROR: When calculating the DY Scale factor you must specify in which step you want to calculate the DY SF:"<<endl;
        cout<<" '_postZcut', '_post2jets', '_postMET', '_post1btag', '_postKinReco' or 'Standard' = _postKinReco"<<endl;
        cout<<"*******************************************************************\n\n"<<endl;
        exit(444);
    }
    if (SpecialComment.BeginsWith("_post")){
        nameAppendix = SpecialComment;
    } else if ( SpecialComment == "Standard") {
        nameAppendix = "_postKinReco";
    }

    cout<<"\n\nBegin DYSCALE FACTOR calculation at selection step "<<nameAppendix<<endl;
    
    vector<TString> Vec_Files = InputFileList("combined", "Nominal");//Read the hardcoded list of files
    if(Vec_Files.size()<1) {cout<<"WARNING(in DYScaleFactor)!!! No datasets available to calculate DY SF. EXITING!!"<<endl; return;}
    
    double NoutEEDYMC=0, NinEEDYMC=0, NoutMuMuDYMC=0, NinMuMuDYMC=0;//Number of events in/out of z-veto region for the DY MC
    double NinEE=0, NinMuMu=0, NinEMu=0;//Number of events in z-veto region for data
    double NinEEloose=0, NinMuMuloose=0;//Number of data events in Z-Veto region with MET cut
    double NinEEMC=0, NinMuMuMC=0;//All other MC events

    for(size_t i=0; i < Vec_Files.size(); i++){
        double LumiWeight = CalcLumiWeight(Vec_Files.at(i));
        double allWeights=LumiWeight;//calculate here all the flat-weights we apply: Lumi*others*...
        if(Vec_Files.at(i).Contains("ee") || Vec_Files.at(i).Contains("mumu")){
            if(Vec_Files.at(i).Contains("run")){
                TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString("Zh1").Append(nameAppendix));
                TH1D *htemp1 = fileReader->GetClone<TH1D>(Vec_Files.at(i), "Looseh1");
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
                    TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString("Zh1").Append(nameAppendix));
                    TH1D *htemp1 = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString("TTh1").Append(nameAppendix));
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
                    delete htemp; delete htemp1;
                }
                else{
                    TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString("TTh1").Append(nameAppendix));
                    ApplyFlatWeights(htemp, LumiWeight);
                    if(Vec_Files.at(i).Contains("ee")){   NoutEEDYMC+=htemp->Integral();}
                    if(Vec_Files.at(i).Contains("mumu")){ NoutMuMuDYMC+=htemp->Integral();}
                    delete htemp;
                }
            }
            else{
                TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString("Zh1").Append(nameAppendix));
                ApplyFlatWeights(htemp, LumiWeight);
                if(Vec_Files.at(i).Contains("ee")){   NinEEMC+=htemp->Integral();   }
                if(Vec_Files.at(i).Contains("mumu")){ NinMuMuMC+=htemp->Integral(); }
                delete htemp;
            }
        }
        
        if(Vec_Files.at(i).Contains("emu") && Vec_Files.at(i).Contains("run")){
            TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString("Zh1").Append(nameAppendix));
            ApplyFlatWeights(htemp, LumiWeight);
            NinEMu+=htemp->Integral();
            delete htemp;
        }

    }
    
    
    double kee = sqrt(NinEEloose/NinMuMuloose);
//     printf("kee = sqrt(%.2f/%.2f) = %.2f\n", NinEEloose, NinMuMuloose, kee);
    
    double kmumu = sqrt(NinMuMuloose/NinEEloose);
//     printf("kmumu = sqrt(%.2f/%.2f) = %.2f\n", NinMuMuloose, NinEEloose, kmumu);
    
    double RoutinEE = NoutEEDYMC/NinEEDYMC;
//     printf("RoutinEE = %.2f/%.2f = %.2f\n", NoutEEDYMC, NinEEDYMC, RoutinEE);
    
    double RoutinMuMu = NoutMuMuDYMC/NinMuMuDYMC;
//     printf("RoutinMuMu = %.2f/%.2f = %.2f\n", NoutMuMuDYMC, NinMuMuDYMC, RoutinMuMu);
    
    double NoutMCEE = RoutinEE*(NinEE - 0.5*NinEMu*kee);
    double NoutMCMuMu = RoutinMuMu*(NinMuMu - 0.5*NinEMu*kmumu);

    double DYSFEE = NoutMCEE/NoutEEDYMC;
    double DYSFMuMu = NoutMCMuMu/NoutMuMuDYMC;

    cout << endl;
    cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
    cout << "Calculation of DY Scale Factors for '" << name << "'  at selection step "<<nameAppendix << endl;

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

    if (DYSFEE < 0.2 || DYSFMuMu < 0.2) {
        cout << "The DY SF is too low (below 0.2). Something is probably wrong.\n";
        std::exit(1);
    }
    
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

    string ResultsFilestring = outpathResults.Data();
    ResultsFilestring.append(subfolderSpecial.Data());
    ResultsFilestring.append("/");
    ResultsFilestring.append(Systematic);
    ResultsFilestring.append("/");
    ResultsFilestring.append(Channel);
    gSystem->mkdir((char*)ResultsFilestring.c_str(), true);
    ResultsFilestring.append("/");
    ResultsFilestring.append(name);
    ResultsFilestring.append("Results.txt");
    cout<<ResultsFilestring<<endl;
    ResultsFile.open(ResultsFilestring.c_str());
        
    double Xbins[XAxisbins.size()];
    for(size_t i = 0; i<XAxisbins.size();i++){Xbins[i]=XAxisbins[i];} 

    if(flat_Syst > 0.0){// flat systematics will be: flat_Syst != 0.0
        for(size_t bin = 0; bin < XAxisbinCenters.size(); bin++){
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

        TString filename = "preunfolded/Nominal/"+Channel+"/"+name+"_UnfoldingHistos.root";
        theDataHist = fileReader->GetClone<TH1D>(filename, "aDataHist");
        theBgrHist =  fileReader->GetClone<TH1D>(filename, "aBgrHist");
        theTtBgrHist =  fileReader->GetClone<TH1D>(filename, "aTtBgrHist");
        theRecHist =  fileReader->GetClone<TH1D>(filename, "aRecHist");
        theGenHist =  fileReader->GetClone<TH1D>(filename, "aGenHist");
        theRespHist = fileReader->GetClone<TH2D>(filename, "aRespHist");

        TString filenameUp = "preunfolded/"+SystematicUp+"/"+Channel+"/"+name+"_UnfoldingHistos.root";
        theBgrHistUp =  fileReader->GetClone<TH1D>(filenameUp, "aBgrHist");
        theTtBgrHistUp =  fileReader->GetClone<TH1D>(filenameUp, "aTtBgrHist");
        theRecHistUp =  fileReader->GetClone<TH1D>(filenameUp, "aRecHist");
        theGenHistUp =  fileReader->GetClone<TH1D>(filenameUp, "aGenHist");
        theRespHistUp =  fileReader->GetClone<TH2D>(filenameUp, "aRespHist");

        TString filenameDown = "preunfolded/"+SystematicDown+"/"+Channel+"/"+name+"_UnfoldingHistos.root";
        theBgrHistDown =  fileReader->GetClone<TH1D>(filenameDown, "aBgrHist");
        theTtBgrHistDown =  fileReader->GetClone<TH1D>(filenameDown, "aTtBgrHist");
        theRecHistDown =  fileReader->GetClone<TH1D>(filenameDown, "aRecHist");
        theGenHistDown =  fileReader->GetClone<TH1D>(filenameDown, "aGenHist");
        theRespHistDown =  fileReader->GetClone<TH2D>(filenameDown, "aRespHist");

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
        theSpecialPostfix = name;
        if ( specialComment.CompareTo("Standard") != 0 ) {
            //theSpecialPostfix = specialComment;
        } 
        TString theSystematicName = Systematic; 

        cout << endl;
        cout << endl;
    
        // Get the integrals for the normalization
        double totalDataEventsNom  = TopSVDFunctions::SVD_Integral1D(theDataHist, 0, false);
        double totalBgrEventsNom   = TopSVDFunctions::SVD_Integral1D(theBgrHist, 0, false);
        double totalBgrEventsUp    = TopSVDFunctions::SVD_Integral1D(theBgrHistUp, 0, false);
        double totalBgrEventsDown  = TopSVDFunctions::SVD_Integral1D(theBgrHistDown, 0, false);
        double totalTtBgrEventsNom = TopSVDFunctions::SVD_Integral1D(theTtBgrHist, 0, false);
        double totalTtBgrEventsUp  = TopSVDFunctions::SVD_Integral1D(theTtBgrHistUp, 0, false);
        double totalTtBgrEventsDown= TopSVDFunctions::SVD_Integral1D(theTtBgrHistDown, 0, false);
        double totalRecEventsNom   = TopSVDFunctions::SVD_Integral1D(theRecHist, 0, false);
        double totalRecEventsUp    = TopSVDFunctions::SVD_Integral1D(theRecHistUp, 0, false);
        double totalRecEventsDown  = TopSVDFunctions::SVD_Integral1D(theRecHistDown, 0, false);
        double totalGenEventsNom   = TopSVDFunctions::SVD_Integral1D(theGenHist, 0, false);
        double totalGenEventsUp    = TopSVDFunctions::SVD_Integral1D(theGenHistUp, 0, false);
        double totalGenEventsDown  = TopSVDFunctions::SVD_Integral1D(theGenHistDown, 0, false);

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
                if(Systematic == "MASS_"){
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
                if(Systematic == "MASS_"){
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
        if(gSystem->AccessPathName(eefilename) || gSystem->AccessPathName(emufilename) || gSystem->AccessPathName(mumufilename) ||
           gSystem->AccessPathName(eeErrfilename) || gSystem->AccessPathName(emuErrfilename) || gSystem->AccessPathName(mumuErrfilename)){
            cout<<"WARNING (in CalcDiffSystematics)!!"<<endl;
            cout<<"One of the input files you use for the combined XSection measurement doesn't exist!! Exiting!!"<<endl;
            exit(22);
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

//         double binWidth[bins];
        double DiffXSecVecVaried[bins];
        //do the actual combined Diff.XSection calculation
        for(int i=0; i<bins; i++){
//             binWidth[i] = Xbins[i+1]-Xbins[i];
            for(int j=0; j<3; j++){//check if any stat error is 0, in this case set their contribution to 0!!
                if(perChannelDiffXSecStatError[j][i] == 0){
                    perChannelDiffXSecStatError[j][i] = 1e100;
                    perChannelDiffXSecPlot[j][i] = 0;
                }
            }
            DiffXSecVecVaried[i] =((( 1 + perChannelDiffXSecSystError[0][i]) * perChannelDiffXSecPlot[0][i]/(perChannelDiffXSecStatError[0][i]*perChannelDiffXSecStatError[0][i]))
                                  +(( 1 + perChannelDiffXSecSystError[1][i]) * perChannelDiffXSecPlot[1][i]/(perChannelDiffXSecStatError[1][i]*perChannelDiffXSecStatError[1][i]))
                                  +(( 1 + perChannelDiffXSecSystError[2][i]) * perChannelDiffXSecPlot[2][i]/(perChannelDiffXSecStatError[2][i]*perChannelDiffXSecStatError[2][i])))/
                                  ((1/(perChannelDiffXSecStatError[0][i]*perChannelDiffXSecStatError[0][i]))+
                                   (1/(perChannelDiffXSecStatError[1][i]*perChannelDiffXSecStatError[1][i]))+
                                   (1/(perChannelDiffXSecStatError[2][i]*perChannelDiffXSecStatError[2][i])));
            
            double Sys_Error = fabs(DiffXSecVecVaried[i]-perChannelDiffXSecPlot[3][i])/perChannelDiffXSecPlot[3][i];
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
    fileReader = RootFileReader::getInstance();
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
    if ( specialComment.CompareTo("Standard") != 0 ) {
        //subfolderSpecial = specialComment.Prepend("/");
    }

    DYEntry = "Z / #gamma* #rightarrow ee/#mu#mu";

    if(Systematic.Contains("DY_") || Systematic.Contains("BG_")){Systematic = "Nominal";}//We just need to vary the nominal DY and BG systematics

    TString histoListName = "FileLists/HistoFileList_"+Systematic+"_"+mode+".txt";
    cout << "reading " << histoListName << endl;
    ifstream FileList(histoListName);
    if (FileList.fail()) {
        cerr << "Error reading " << histoListName << endl;
        exit(1);
    }
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
        //for (auto s: FileVector) std::cout << "FileVector = " << s << "\n";
        //shouldn't this be sorted? i.e. first runs of all channels, then ...
        //need to ask Ivan !!!FIXME
        return FileVector;
    }

    //data is only stored in the Nominal directory
    TString nominalPath = TString("selectionRoot/Nominal/") + mode + "/" + mode;
    if (!doClosureTest) {
        FileVector.push_back(nominalPath + "_run2012A.root");
        FileVector.push_back(nominalPath + "_run2012Arecover.root");
        FileVector.push_back(nominalPath + "_run2012B.root");
        FileVector.push_back(nominalPath + "_run2012C_24Aug.root");
        FileVector.push_back(nominalPath + "_run2012C_PromptReco.root");
    } else {
        FileVector.push_back(nominalPath + "_ttbarsignalplustau_fakerun_nominal.root");
    }
    
    //MC depends on the specific Systematic: Signal systematics only use different signal samples
    TString tempName;
    if(!Systematic.CompareTo("JER_UP") || !Systematic.CompareTo("JER_DOWN") || !Systematic.CompareTo("JES_UP") || !Systematic.CompareTo("JES_DOWN") ||
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
       !Systematic.CompareTo("MASS_UP") || !Systematic.CompareTo("MASS_DOWN") ||
       !Systematic.CompareTo("MATCHINGUP") || !Systematic.CompareTo("MATCHINGDOWN") ||
       !Systematic.CompareTo("MATCH_UP") || !Systematic.CompareTo("MATCH_DOWN") ||
       !Systematic.CompareTo("SCALE_UP") || !Systematic.CompareTo("SCALE_DOWN")){
        
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
        TH1D *hist = fileReader->GetClone<TH1D>(dataset.at(i), name, true);
//         cout << i << ": " << dataset.at(i) << endl;
        if (!hist) return false;
        if (!name.Contains("_step") && !name.Contains("bcp_") && !name.Contains("Lead") && !name.EndsWith("bkr") && !name.EndsWith("akr")
            && (name.Contains("Lepton") || name.Contains("BJet") || name.Contains("Top")))
        {
            TString stemp = name;
            if(name.Contains("Lepton"))    {stemp.ReplaceAll("Lepton",6,"AntiLepton",10);}
            else if(name.Contains("BJet")) {stemp.ReplaceAll("BJet",4,"AntiBJet",8);}
            else if(name.Contains("Top"))  {stemp.ReplaceAll("Top",3,"AntiTop",7);}
            const TH1D* other = fileReader->Get<TH1D>(dataset.at(i), stemp);
            if (other) hist->Add(other); 
            else cerr << "Cannot find corresponding anti-quantity histogram: "<<stemp<<endl;
        }

        //Rescaling to the data luminosity
        double LumiWeight = CalcLumiWeight(dataset.at(i));
        //cout << "File " << dataset.at(i) << " has weight " << LumiWeight << "\n";

        ApplyFlatWeights(hist, LumiWeight);

        setHHStyle(*gStyle);

        hists.push_back(*hist);
    }
    if (doClosureTest) {
        for (size_t i = 1; i<dataset.size(); ++i) {
            if (!dataset.at(i).Contains("ttbarsignalplustau"))
                hists[0].Add(&hists[i]);
        }
    }
    initialized=true;
    return true;
}

void Plotter::addAndDelete_or_Assign(TH1*& addToThis, TH1* addThis)
{
    if (!addToThis) addToThis = addThis;
    else {
        addToThis->Add(addThis);
        delete addThis;
    }
}


void Plotter::write(TString Channel, TString Systematic) // do scaling, stacking, legending, and write in file 
{
    setDataSet(Channel,Systematic);
    if (!fillHisto()) return;

    if (hists.size() == 0) { 
        cerr << "***ERROR! No histograms available! " << Channel << "/" << Systematic << endl; 
        exit(11); 
    }
        
    TCanvas * c = new TCanvas("","");

    THStack * stack=  new THStack("def", "def");
    TLegend * leg =  new TLegend(0.70,0.55,0.98,0.85);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetX1NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.25);
    leg->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 - leg->GetNRows()*0.04);
    leg->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength());
    leg->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength());

    TH1 *drawhists[hists.size()];
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


    gSystem->mkdir("preunfolded/"+Systematic+"/"+Channel, true);


    for(unsigned int i=0; i<hists.size() ; ++i){ // prepare histos and leg
//         cout << "Legend ["<<i<<"] = " << legends.at(i) << endl;
    if(legends.at(i) != "Data"){
        //      drawhists[i]->Scale(12.1/5.1);

        if(XAxisbins.size()>1){//only distributions we want to unfold will have a binning vector
            if(legends.at(i) == "t#bar{t} Signal"){
                TString ftemp = dataset.at(i);
		double LumiWeight = CalcLumiWeight(dataset.at(i));
                if (!init) {
                    aRespHist = fileReader->GetClone<TH2>(ftemp, "GenReco"+newname);
                    aGenHist = fileReader->GetClone<TH1D>(ftemp, "VisGen"+newname);
                    //Rebin(bins,"aTtBgrHist",Xbins);
                    if(!newname.Contains("Lead") && (newname.Contains("Lepton")||newname.Contains("Top")||newname.Contains("BJet"))){
                        aRespHist->Add(fileReader->Get<TH2>(ftemp, "GenRecoAnti"+newname));
                        aGenHist->Add(fileReader->Get<TH1D>(ftemp, "VisGenAnti"+newname));
                    }
                    init =true;
                } else {//account for more than one signal histogram
                    aRespHist->Add(fileReader->Get<TH2>(ftemp, "GenReco"+newname));
                    aGenHist->Add(fileReader->Get<TH1D>(ftemp, "VisGen"+newname));
                    if(!newname.Contains("Lead") && (newname.Contains("Lepton")||newname.Contains("Top")||newname.Contains("BJet"))){
                        aGenHist->Add(fileReader->Get<TH1D>(ftemp, "VisGenAnti"+newname));
                        aRespHist->Add(fileReader->Get<TH2>(ftemp, "GenRecoAnti"+newname));
                    }
                }
		aRespHist->Scale(LumiWeight);
		aGenHist->Scale(LumiWeight);
	    }
            //cout<<"Legend: "<<legends.at(i)<<endl;
            if(legends.at(i) == "t#bar{t} Signal"){
                signalHist = i; //What is happening for the combined channel?? only 1 integer value??
                aRecHist = drawhists[i]->Rebin(bins,"aRecHist",Xbins); //What happens in the combined channel? Are the 3 channels summed up?? NO!!
                //cout<<"Added "<<legends.at(i)<<" to aRecHist"<<endl;
            }else if(legends.at(i) == "t#bar{t} Other"){//IMPORTANT: TTbar Other are added to the ttbarbackground histogram AND the Background Hist gram
                addAndDelete_or_Assign(aTtBgrHist, drawhists[i]->Rebin(bins,"aTtBgrHist",Xbins));
                addAndDelete_or_Assign(aBgrHist, drawhists[i]->Rebin(bins,"aTtBgrHist",Xbins));
            }else if((legends.at(i) == DYEntry)){
                if (channelType!=2) drawhists[i]->Scale(DYScale[channelType]);

                //Here we take into account the systematic shifts needed for DY systematic because it only modifies the nominal dataset
                if(Systematic == "DY_UP"){
                    drawhists[i]->Scale(1.3);
                }
                if(Systematic == "DY_DOWN"){
                    drawhists[i]->Scale(0.7);
                }
                addAndDelete_or_Assign(aBgrHist, drawhists[i]->Rebin(bins,"aBgrHist",Xbins));
            }else{
                //Here we take into account the systematic shifts needed for BG systematic because it only modifies the nominal dataset
                if(Systematic == "BG_UP"){
                    drawhists[i]->Scale(1.5);
                }
                if(Systematic == "BG_DOWN"){
                    drawhists[i]->Scale(0.5);
                }
                addAndDelete_or_Assign(aBgrHist, drawhists[i]->Rebin(bins,"aBgrHist",Xbins));
            }
        }

        if(i > 1){
            if(legends.at(i) != legends.at(i-1)){
                legchange = i; 
                if((legends.at(i) == DYEntry)&& DYScale[channelType] != 1) leg->AddEntry(drawhists[i], legends.at(i),"f");
                else leg->AddEntry(drawhists[i], legends.at(i) ,"f");
            }else{
                drawhists[legchange]->Add(drawhists[i]);
            }
        }

        if(i!=(hists.size()-1)){
            if(legends.at(i)!=legends.at(i+1)){
                drawhists[i]->SetLineColor(1);
            }
        }else{
            drawhists[i]->SetLineColor(1);
        }
        if(legends.at(i) != legends.at(i-1)){
            drawhists[i]->SetLineColor(1);
            stack->Add(drawhists[i]); 
        }
    }
    else{
        if(i==0) leg->AddEntry(drawhists[i], legends.at(i) ,"pe");
        if(i>0){
            if(legends.at(i) != legends.at(i-1)){
                leg->AddEntry(drawhists[i], legends.at(i) ,"pe");
            }
            if(legends.at(i) == legends.at(0)){
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
        delete f15;
    }

    leg = ControlLegend(hists.size(), drawhists, legends, leg);

    if(name.Contains("HypjetMultiXSec")){

        double InclusiveXsectionWrite[4], InclusiveXsectionStatErrorWrite[4];
        CalcXSec(dataset, InclusiveXsectionWrite, InclusiveXsectionStatErrorWrite, Systematic,"");

        MakeTable();
        if(channelType==3) PlotXSec();
    }


    TList* l = stack->GetHists(); 
    TH1D* stacksum = (TH1D*) l->At(0)->Clone();

    for (int i = 1; i < l->GetEntries(); ++i) { stacksum->Add((TH1D*)l->At(i));}

    TH1D* syshist =0;
    syshist = (TH1D*)stacksum->Clone();
    for(Int_t i=0; i<=syshist->GetNbinsX(); ++i){
        Double_t binc = 0;
        binc += stacksum->GetBinContent(i);
        syshist->SetBinContent(i, binc);
    }

    if(logY)c->SetLogy();
    syshist->SetFillStyle(3004);
    syshist->SetFillColor(kBlack);

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

    DrawCMSLabels(1, 8);
    DrawDecayChLabel(channelLabel[channelType]);
    leg->Draw("SAME");
    //drawRatio(drawhists[0], stacksum, 0.5, 1.7);

    // Create Directory for Output Plots 
    gSystem->mkdir(outpathPlots+"/"+subfolderChannel+"/"+Systematic, true);
    c->Print(outpathPlots+subfolderChannel+"/"+Systematic+"/"+name+".eps");
    //c->Print(outpathPlots+subfolderChannel+"/"+Systematic+"/"+name+".root");
    //c->Print(outpathPlots+subfolderChannel+"/"+Systematic+"/"+name+".C");
    
    TH1 *sumMC = 0; 
    TH1 *sumttbar = 0;
    for (size_t i = 0; i < hists.size(); ++i) {
        if (legends.at(i) != "Data") {
            if (sumMC) sumMC->Add(drawhists[i]);
            else sumMC = static_cast<TH1*>(drawhists[i]->Clone());
        }
        if (legends.at(i) == "t#bar{t} Signal") {
            if (sumttbar) sumttbar->Add(drawhists[i]);
            else sumttbar = static_cast<TH1*>(drawhists[i]->Clone());
        }
    }
    sumMC->SetName(name);
    
    //save Canvas AND sources in a root file
    TFile out_root(outpathPlots+subfolderChannel+"/"+Systematic+"/"+name+"_source.root", "RECREATE");
    drawhists[0]->Write(name+"_data");
    sumttbar->Write(name+"_signalmc");
    sumMC->Write(name+"_allmc");
    c->Write(name + "_canvas");
    out_root.Close();
    
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

    if(legends.at(i) == "Data"){
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

    TString channel_array[] = {"ee","mumu","emu","combined"};
    TString sys_array[] = {"PDF_", "HAD_", "MATCH_", "MASS_", "SCALE_", "BTAG_", "BTAG_LJET_", "KIN_", "LEPT_", "TRIG_", "BG_", "DY_", "PU_", "JER_", "JES_"};//For the time being uintil all systematics are finalished
    vector<TString> vec_systematic (sys_array, sys_array + sizeof(sys_array)/sizeof(sys_array[0]));
    vector<TString> vec_channel (channel_array, channel_array + sizeof(channel_array)/sizeof(channel_array[0]));

    double BR_Error = 0.015;
    double Lumi_Error = 0.044;

    double InclusiveXsectionPlot[4], InclusiveXsectionStatErrorPlot[4], InclusiveXsectionSysErrorPlot[4], InclusiveXsectionTotalErrorPlot[4];
    for (int j=0; j<(int)vec_channel.size(); j++){

        ifstream SysResultsList("Plots/Nominal/"+vec_channel.at(j)+"/InclXSec.txt");
        TString DUMMY;
        SysResultsList>>DUMMY>>DUMMY>>DUMMY>>DUMMY>>DUMMY>>InclusiveXsectionPlot[j]>>DUMMY>>InclusiveXsectionStatErrorPlot[j];
        SysResultsList.close();

        ofstream OutputFile(string("Plots/"+vec_channel.at(j)+"/InclXSecResultLateX.txt"));
        OutputFile<<"Inclusive XSection Numerical Results for channel "<<vec_channel.at(j)<<endl;

        double syst_square_for_channel=0.0;
        for (int i=0; i<(int) vec_systematic.size(); i++){

            ifstream SysUP, SysDOWN;

            if(vec_systematic.at(i) != "HAD_"){
                SysUP.open("Plots/"+vec_systematic.at(i)+"UP/"+vec_channel.at(j)+"/InclXSec.txt");
                SysDOWN.open("Plots/"+vec_systematic.at(i)+"DOWN/"+vec_channel.at(j)+"/InclXSec.txt");
            } else{
                SysUP.open("Plots/MCATNLO/"+vec_channel.at(j)+"/InclXSec.txt");
                SysDOWN.open("Plots/POWHEG/"+vec_channel.at(j)+"/InclXSec.txt");
            }
            double VarUp = 0, VarDown = 0, StatErrUp = 0, StatErrDown = 0;

            SysUP>>DUMMY>>DUMMY>>DUMMY>>DUMMY>>DUMMY>>VarUp>>DUMMY>>StatErrUp;
            SysDOWN>>DUMMY>>DUMMY>>DUMMY>>DUMMY>>DUMMY>>VarDown>>DUMMY>>StatErrDown;
            SysUP.close();
            SysDOWN.close();

            //systematic error in %
            double sys_err=(TMath::Abs(InclusiveXsectionPlot[j]-VarUp)+TMath::Abs(InclusiveXsectionPlot[j]-VarDown))*0.5/InclusiveXsectionPlot[j];
            if(vec_systematic.at(i).Contains("MASS_")) {sys_err = sys_err/12;}
            syst_square_for_channel+=sys_err*sys_err;
            OutputFile<<vec_systematic.at(i)<<" (%): "<<setprecision(3)<<sys_err*100<<endl;
        }
        OutputFile<<"BranchingRatio (%): "<<setprecision(3)<<BR_Error*100<<endl;
        OutputFile<<"Luminosity (%): "<<setprecision(3)<<Lumi_Error*100<<endl;

        InclusiveXsectionSysErrorPlot[j]=TMath::Sqrt(syst_square_for_channel + BR_Error * BR_Error + Lumi_Error * Lumi_Error);

        InclusiveXsectionTotalErrorPlot[j] = sqrt(InclusiveXsectionStatErrorPlot[j]*InclusiveXsectionStatErrorPlot[j] +
                                                  InclusiveXsectionPlot[j]*InclusiveXsectionSysErrorPlot[j]*InclusiveXsectionPlot[j]*InclusiveXsectionSysErrorPlot[j]
                                                 );

        OutputFile<<"\n\n*******************************************************************************\n\n";
        OutputFile<<" InclXsec[pb]     Stat.[pb]    Syst.[pb]   Total[pb]"<<endl;
        OutputFile<<setprecision(6)<<InclusiveXsectionPlot[j]<<" +- "<<setprecision(3)<<InclusiveXsectionStatErrorPlot[j]<<" +- "<<setprecision(4)<<InclusiveXsectionSysErrorPlot[j]*InclusiveXsectionPlot[j]<<" +- "<<setprecision(4)<<InclusiveXsectionTotalErrorPlot[j];
        OutputFile.close();
    }

    // measured results with statistical error
    Double_t mx[]   = {      0.50,       1.50,       2.50,       3.50};
    Double_t mexl[] = {      0.00,       0.00,       0.00,       0.00};
    Double_t mexh[] = {      0.00,       0.00,       0.00,       0.00};

    TGraphAsymmErrors *mplot = new TGraphAsymmErrors(4, mx, InclusiveXsectionPlot, mexl, mexh,InclusiveXsectionStatErrorPlot, InclusiveXsectionStatErrorPlot);
    mplot->SetMarkerStyle(20);
    mplot->GetYaxis()->SetNoExponent(kTRUE);
    mplot->SetMarkerColor(kBlack);
    mplot->SetMarkerSize(1.5);
    mplot->SetLineColor(kBlack);

    TGraphAsymmErrors *mplotwithsys = new TGraphAsymmErrors(4, mx, InclusiveXsectionPlot, mexl, mexh,InclusiveXsectionTotalErrorPlot, InclusiveXsectionTotalErrorPlot);
    mplotwithsys->SetMarkerStyle(20);
    mplotwithsys->SetMarkerColor(kBlack);
    mplotwithsys->SetMarkerSize(1.5);
    mplotwithsys->SetLineColor(kBlack);

    // kidonakis
    Double_t kidonmean = 234;
    Double_t kidonx[]   = {    -0.5,     0.5,   1.5,     2.5,     3.5,     4.5};
    Double_t kidony[]   = {kidonmean,kidonmean,kidonmean,kidonmean,kidonmean,kidonmean};
    Double_t kidonexl[] = {      .4,    .4,      .5,      .5,      .5,      .5};
    Double_t kidonexh[] = {      .5,    .5,      .5,      .5,      .4,      .4};
    Double_t kidoneyl[] = {    15.6,    15.6,    15.6,  15.6,    15.6,    15.6};
    Double_t kidoneyh[] = {    13.9,    13.9,    13.9,  13.9,    13.9,    13.9};

    TGraphAsymmErrors *kidonplot = new TGraphAsymmErrors(6, kidonx, kidony, kidonexl, kidonexh, kidoneyl, kidoneyh);
    kidonplot->SetLineColor(kGreen+1);
    kidonplot->SetLineWidth(4);
    kidonplot->SetFillColor(kGreen+1);
    kidonplot->SetFillStyle(3004);

    // mcfm
    Double_t mcfmmean = 225.197;
    Double_t mcfmx[]   = {      -0.5,     0.5,     1.5,     2.5,     3.5,     4.5};
    Double_t mcfmy[]   = {mcfmmean,mcfmmean,mcfmmean,mcfmmean,mcfmmean,mcfmmean};
    Double_t mcfmexl[] = {        .4,      .4,      .5,      .5,      .5,      .5};
    Double_t mcfmexh[] = {        .5,      .5,      .5,      .5,      .4,      .4};
    Double_t mcfmeyl[] = {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};
    Double_t mcfmeyh[] = {   0.0, 0.0,  0.0,  0.0,  0.0,  0.0};

    TGraphAsymmErrors *mcfmplot = new TGraphAsymmErrors(6, mcfmx, mcfmy, mcfmexl, mcfmexh, mcfmeyl, mcfmeyh);
    mcfmplot->SetLineColor(kBlue+1);
    mcfmplot->SetLineWidth(4);
    mcfmplot->SetFillColor(kBlue+1);
    mcfmplot->SetFillStyle(3005);

    TH1F* framehist = new TH1F("framehist","",4,0.,4.);
    framehist->SetMinimum(100);
    framehist->SetMaximum(380);
    framehist->GetXaxis()->SetTickLength(0);
    framehist->GetXaxis()->SetBinLabel(1,"");
    framehist->GetXaxis()->SetBinLabel(2,"");
    framehist->GetXaxis()->SetBinLabel(3,"");
    framehist->GetYaxis()->SetTitle("#sigma [pb]");
    framehist->GetYaxis()->CenterTitle(kTRUE);
    framehist->GetYaxis()->SetNoExponent(kTRUE);

    TPaveText* box1 = new TPaveText(0.25,0.20,0.33,0.30,"NDC");
    box1->SetFillColor(10);
    box1->SetTextSize(0.04);
    box1->AddText("ee");

    TPaveText* box2 = new TPaveText(0.44,0.20,0.52,0.30,"NDC");
    box2->SetFillColor(10);
    box2->SetTextSize(0.04);
    box2->AddText("#mu#mu");

    TPaveText* box3 = new TPaveText(0.62,0.20,0.72,0.30,"NDC");
    box3->SetFillColor(10);
    box3->SetTextSize(0.04);
    box3->AddText("e#mu");

    TPaveText* box4 = new TPaveText(0.82,0.20,0.90,0.30,"NDC");
    box4->SetFillColor(10);
    box4->SetTextSize(0.04);
    box4->AddText("combined");

    TLegend* leg =  new TLegend( 0.56, 0.65, 0.89, 0.85);
    leg->SetBorderSize( 0 );
    leg->SetFillColor( 0 );
    leg->SetTextFont(62);
    leg->SetTextSize(0.03);
    leg->AddEntry( mplot,       "Measurements",            "p"  );
    leg->AddEntry( mcfmplot, "MCFM #otimes CTQE66M", "lf" );
    leg->AddEntry( kidonplot,    "Kidonakis #otimes MSTW2008 NNLO",     "lf" );

    TCanvas* c = new TCanvas("plot", "plot", 1200, 800);
    framehist->Draw();
    mcfmplot->Draw("C,2,SAME");
    kidonplot->Draw("C,2,SAME");
    gStyle->SetEndErrorSize(8);
    mplot->Draw("p,SAME");
    mplotwithsys->Draw("p,SAME,Z");
    leg ->Draw("SAME");

    box1->Draw("SAME");
    box2->Draw("SAME");
    box3->Draw("SAME");
    box4->Draw("SAME");
   gSystem->mkdir(outpathPlots+subfolderChannel+subfolderSpecial, true);
    c->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/"+"InclusiveXSec.eps");
    c->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/"+"InclusiveXSec.C");
    c->Clear();
    delete c;

}


void Plotter::MakeTable(){

    TH1D *numhists5[hists.size()];
    TH1D *numhists6[hists.size()];
    TH1D *numhists7[hists.size()];
    TH1D *numhists8[hists.size()];
    TH1D *numhists9[hists.size()];

    for(unsigned int i=0; i<dataset.size(); i++){

        TH1D *temp_hist5 = fileReader->GetClone<TH1D>(dataset[i], "step5");
        TH1D *temp_hist6 = fileReader->GetClone<TH1D>(dataset[i], "step6");
        TH1D *temp_hist7 = fileReader->GetClone<TH1D>(dataset[i], "step7");
        TH1D *temp_hist8 = fileReader->GetClone<TH1D>(dataset[i], "step8");
        TH1D *temp_hist9 = fileReader->GetClone<TH1D>(dataset[i], "step9");
        
        double LumiWeight = CalcLumiWeight(dataset.at(i));
        ApplyFlatWeights(temp_hist5, LumiWeight);
        ApplyFlatWeights(temp_hist6, LumiWeight);
        ApplyFlatWeights(temp_hist7, LumiWeight);
        ApplyFlatWeights(temp_hist8, LumiWeight);
        ApplyFlatWeights(temp_hist9, LumiWeight);

        numhists5[i]=temp_hist5;
        numhists6[i]=temp_hist6;
        numhists7[i]=temp_hist7;
        numhists8[i]=temp_hist8;
        numhists9[i]=temp_hist9;

    }

    for(unsigned int i=0; i<hists.size() ; i++){ // prepare histos and leg
        if((legends.at(i) == DYEntry) && channelType!=2){
            //numhists5[i]->Scale(DYScale[channelType]);//DYscale not applied in step5 and 6?
            //numhists6[i]->Scale(DYScale[channelType]);
            numhists7[i]->Scale(DYScale.at(channelType));
            numhists8[i]->Scale(DYScale.at(channelType));
            numhists9[i]->Scale(DYScale.at(channelType));
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
    gSystem->mkdir(outpathPlots+"/"+subfolderChannel+"/"+subfolderSpecial, true);  
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
            EventFile5<<legends.at(i)<<": "<<tmp_num5<<endl;
            EventFile6<<legends.at(i)<<": "<<tmp_num6<<endl;
            EventFile7<<legends.at(i)<<": "<<tmp_num7<<endl;
            EventFile8<<legends.at(i)<<": "<<tmp_num8<<endl;
            EventFile9<<legends.at(i)<<": "<<tmp_num9<<endl;
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
        }
        else if(legends.at(i)!=legends.at(i+1)){
            EventFile5<<legends.at(i)<<": "<<tmp_num5<<endl;
            EventFile6<<legends.at(i)<<": "<<tmp_num6<<endl;
            EventFile7<<legends.at(i)<<": "<<tmp_num7<<endl;
            EventFile8<<legends.at(i)<<": "<<tmp_num8<<endl;
            EventFile9<<legends.at(i)<<": "<<tmp_num9<<endl;
            if(legends.at(i)!="Data"){
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
    cout<<"\nEvent yields saved in "<<EventFilestring5.c_str()<<"\n"<<endl;
}

double Plotter::CalcXSec(std::vector<TString> datasetVec, double InclusiveXsectionVec[4],double InclusiveXsectionStatErrorVec[4], TString Systematic, TString Shift){

    double BranchingFraction[4]={0.01166, 0.01166, 0.02332, 0.04666};//[ee, mumu, emu, combined] not including tau

    TH1D *numhists[hists.size()];
    double numbers[5]={0., 0., 0., 0., 0.};//[0]=data, [1]=Signal, [2]Signal(only lumi & PU weights), [3]ttbar background, [4]background(non-ttbar)
    double error_numbers[5]={0., 0., 0., 0., 0.};//Square of error: [0]=data, [1]=Signal, [2]Signal(only lumi & PU weights), [3]ttbar background, [4]background(non-ttbar)
//     double TTbarBGnum =0;

    if (Systematic.Contains("UP"))  { Shift="Up";}
    if (Systematic.Contains("DOWN")){ Shift="Down";}

    for(unsigned int i=0; i<datasetVec.size(); i++){
        TH1D *hist = fileReader->GetClone<TH1D>(datasetVec[i], "step9");
        //hist->Sumw2();

        double LumiWeight = CalcLumiWeight(datasetVec.at(i));
        ApplyFlatWeights(hist, LumiWeight);

        numhists[i]=hist;
    }

    for(unsigned int i=0; i<hists.size() ; i++){ // prepare histos and leg 
        if(legends.at(i) == "Data"){
            numbers[0]+=numhists[i]->Integral();
            error_numbers[0]+=numhists[i]->GetBinError(2) * numhists[i]->GetBinError(2); //This bin selection is hardcoded please change it if changes when filling in Analysis.C
        }
        else if(legends.at(i) == "t#bar{t} Signal"){
            TH1D *NoPUPlot = fileReader->GetClone<TH1D>(datasetVec.at(i), "step9");

            double LumiWeight = CalcLumiWeight(datasetVec.at(i));
            ApplyFlatWeights(NoPUPlot, LumiWeight);

            numbers[1]+=NoPUPlot->Integral();
            error_numbers[1]+=NoPUPlot->GetBinError(2) * NoPUPlot->GetBinError(2); //This bin selection is hardcoded please change it if changes when filling in Analysis.C

            TH1D *GenPlot = fileReader->GetClone<TH1D>(datasetVec.at(i), "GenAll");
            ApplyFlatWeights(GenPlot, LumiWeight);

            numbers[2]+=GenPlot->Integral();
            error_numbers[2]+=GenPlot->GetBinError(18) * GenPlot->GetBinError(18); //This bin selection is hardcoded please change it if changes when filling in Analysis.C

        }
        else if(legends.at(i) == "t#bar{t} Other"){
            numbers[3]+=numhists[i]->Integral();
            error_numbers[3]+=numhists[i]->GetBinError(2) * numhists[i]->GetBinError(2); //This bin selection is hardcoded please change it if changes when filling in Analysis.C
        }
        else {
            if((legends.at(i) == DYEntry) && channelType!=2){
                numhists[i]->Scale(DYScale[channelType]);
            }
            if((legends.at(i) == DYEntry) && Systematic.Contains("DY_") && Shift == "Up"){
                numhists[i]->Scale(1.3);
            }
            if((legends.at(i) == DYEntry) && Systematic.Contains("DY_") && Shift == "Down"){
                numhists[i]->Scale(0.7);
            }
            if(Systematic.Contains("BG_") && Shift=="Up" && legends.at(i)!= "t#bar{t} Other" && legends.at(i) != DYEntry){
                numhists[i]->Scale(1.3);
            }
            if(Systematic.Contains("BG_") && Shift=="Down" && legends.at(i)!= "t#bar{t} Other" && legends.at(i) != DYEntry){
                numhists[i]->Scale(0.7);
            }
            numbers[4]+=numhists[i]->Integral();
            error_numbers[4]+=numhists[i]->GetBinError(2) * numhists[i]->GetBinError(2); //This bin selection is hardcoded please change it if changes when filling in Analysis.C
        }
    }

    ////////////////////////////Make output for tables

    double tmp_num = 0;

//   double signalFraction = 0; 

//   signalFraction = numbers[1]/(numbers[1]+TTbarBGnum); // is 1 right now, since TTbarBGnum is 0

    ofstream EventFile, XSecFile;
    string EventFilestring = outpathPlots.Data();
    EventFilestring.append("/");
    EventFilestring.append(Systematic.Data());
    gSystem->MakeDirectory(EventFilestring.c_str());
    EventFilestring.append(subfolderChannel.Data());
    gSystem->MakeDirectory(EventFilestring.c_str());
    EventFilestring.append(subfolderSpecial.Data());
    gSystem->MakeDirectory(EventFilestring.c_str());
    string XSecFileString = EventFilestring;
    EventFilestring.append("/Events.txt");
    XSecFileString.append("/InclXSec.txt");// File to store the Inclusive XSection and stat error to be reused later to make the plots.
    EventFile.open(EventFilestring.c_str());
    XSecFile.open(XSecFileString.c_str());

    double bg_num = 0;
    for(unsigned int i=0; i<hists.size() ; i++){ 
        tmp_num+=numhists[i]->Integral();

        if(i==(hists.size()-1)){
            EventFile<<legends.at(i)<<": "<<tmp_num<<endl;
            bg_num+=tmp_num;
            tmp_num=0;
        }
        else if(legends.at(i)!=legends.at(i+1)){
            EventFile<<legends.at(i)<<": "<<tmp_num<<endl;
            if(legends.at(i)!="Data")bg_num+=tmp_num;
            tmp_num=0;
        }

    }
    EventFile<<"Total MCs: "<<bg_num<<endl;
    EventFile<<"\nDataEvents= "<<numbers[0]<<endl;
    EventFile<<"SignalReco= "<<numbers[1]<<endl;
    EventFile<<"SignalGen = "<<numbers[2]<<endl;
    EventFile<<"ttbar bags= "<<numbers[3]<<endl;
    EventFile<<"All Backgd= "<<numbers[4]<<endl;
    EventFile<<"Efficiency= "<<(numbers[1]/numbers[2])<<endl;
    EventFile<<"BrancRatio= "<<BranchingFraction[channelType]<<endl;

    double xsec = ( (numbers[0]-numbers[4]) * (numbers[1]/(numbers[1]+numbers[3])) ) / ( (numbers[1]/numbers[2])*BranchingFraction[channelType]*lumi);
    //    double xsecstaterror = xsec * ( TMath::Sqrt(error_numbers[0]) + TMath::Sqrt(error_numbers[4]) ) / (numbers[0] - numbers[4]);  //relative statistical error
    double xsecstaterror = TMath::Sqrt(error_numbers[0]) * (numbers[1]/(numbers[1]+numbers[3])) / ( (numbers[1]/numbers[2])*BranchingFraction[channelType]*lumi);

    if(channelType!=3){
        InclusiveXsectionVec[channelType] = xsec;
        InclusiveXsectionStatErrorVec[channelType] = xsecstaterror;
    }
    else{
        TString eefilename="Plots/"+Systematic+"/ee/InclXSec.txt";
        TString mumufilename="Plots/"+Systematic+"/mumu/InclXSec.txt";
        TString emufilename="Plots/"+Systematic+"/emu/InclXSec.txt";

        //check the existence of the file
        if( gSystem->AccessPathName(eefilename) || gSystem->AccessPathName(emufilename) || gSystem->AccessPathName(mumufilename)){
            cout<<"WARNING (in CalcXSec)!!"<<endl;
            cout<<"One of the input files you use for the combined XSection measurement doesn't exist!! Exiting!!"<<endl;
            return 0;
        }

        ifstream ResultsEE(eefilename);
        ifstream ResultsEMu(emufilename);
        ifstream ResultsMuMu(mumufilename);

        TString Dummy="";

        ResultsEE>>Dummy>>Dummy>>Dummy>>Dummy>>Dummy>>InclusiveXsectionVec[0]>>Dummy>>InclusiveXsectionStatErrorVec[0];
        ResultsMuMu>>Dummy>>Dummy>>Dummy>>Dummy>>Dummy>>InclusiveXsectionVec[1]>>Dummy>>InclusiveXsectionStatErrorVec[1];
        ResultsEMu>>Dummy>>Dummy>>Dummy>>Dummy>>Dummy>>InclusiveXsectionVec[2]>>Dummy>>InclusiveXsectionStatErrorVec[2];

        ResultsEE.close(); ResultsEMu.close(); ResultsMuMu.close();

        InclusiveXsectionVec[channelType] =( InclusiveXsectionVec[0]/(InclusiveXsectionStatErrorVec[0]*InclusiveXsectionStatErrorVec[0])
                                            +InclusiveXsectionVec[1]/(InclusiveXsectionStatErrorVec[1]*InclusiveXsectionStatErrorVec[1])
                                            +InclusiveXsectionVec[2]/(InclusiveXsectionStatErrorVec[2]*InclusiveXsectionStatErrorVec[2]) 
                                            )/
                                            ( 1/(InclusiveXsectionStatErrorVec[0]*InclusiveXsectionStatErrorVec[0])
                                            + 1/(InclusiveXsectionStatErrorVec[1]*InclusiveXsectionStatErrorVec[1])
                                            + 1/(InclusiveXsectionStatErrorVec[2]*InclusiveXsectionStatErrorVec[2])
                                            );

        InclusiveXsectionStatErrorVec[channelType] =1/(TMath::Sqrt(
                                            (1/(InclusiveXsectionStatErrorVec[0]*InclusiveXsectionStatErrorVec[0]))
                                            +(1/(InclusiveXsectionStatErrorVec[1]*InclusiveXsectionStatErrorVec[1]))
                                            +(1/(InclusiveXsectionStatErrorVec[2]*InclusiveXsectionStatErrorVec[2]))
                                                                    ));
    }
    EventFile<<"XSection  = "<<InclusiveXsectionVec[channelType]<<endl;
    EventFile<<"XSecStaErr= "<<InclusiveXsectionStatErrorVec[channelType]<<endl;
    EventFile.close();
    XSecFile<<"Systematic: "<<Systematic<<" Channel: "<<subfolderChannel<<" InclXSection: "<<InclusiveXsectionVec[channelType]<<" AbsStatError: "<<InclusiveXsectionStatErrorVec[channelType]<<endl;
    XSecFile.close();
    cout<<"Inclusive XSection information saved in: "<<XSecFileString.c_str()<<endl;
    return xsec;
}

void Plotter::CalcDiffXSec(TString Channel, TString Systematic){
    cout<<"\n\n Preparing to calculate the central value of DiffXSection in channel '"<<Channel<<"' for systematic '"<<Systematic<<"'"<<endl;
    
    double Xbins[XAxisbins.size()];
    double binWidth[XAxisbinCenters.size()];
    for(unsigned int i = 0; i<XAxisbins.size();i++){Xbins[i]=XAxisbins[i];}
    double GenSignalSum[XAxisbinCenters.size()];

    TString ftemp = "preunfolded/"+Systematic+"/"+Channel+"/"+name+"_UnfoldingHistos.root";
    TH1* theDataHist =  fileReader->GetClone<TH1D>(ftemp, "aDataHist");
    TH1* theBgrHist =  fileReader->GetClone<TH1D>(ftemp, "aBgrHist");
    TH1* theTtBgrHist =  fileReader->GetClone<TH1D>(ftemp, "aTtBgrHist");
    TH1* theRecHist =  fileReader->GetClone<TH1D>(ftemp, "aRecHist");
    TH1* theGenHist =  fileReader->GetClone<TH1D>(ftemp, "aGenHist");
    TH1* theRespHist =  fileReader->GetClone<TH2D>(ftemp, "aRespHist");

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
        //if (name.Contains("Lead")) theSpecialPostfix = name;
        theSpecialPostfix = name;
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
            ResultsMuMu>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecPlot[1][bin]>>Dummy>>perChannelDiffXSecStatError[1][bin]>>Dummy>>perChannelGenDiffXSec[1][bin];
            ResultsEMu>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecPlot[2][bin]>>Dummy>>perChannelDiffXSecStatError[2][bin]>>Dummy>>perChannelGenDiffXSec[2][bin];
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

    gSystem->mkdir("UnfoldingResults/"+Systematic+"/"+Channel, true);

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
    ResultsLateX<<"Bin Center & Bin & 1/#sigma d#sigma/dX & stat(%) & syst(%) & total(%)"<<endl;
    for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
        //cout<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[bin]<<" bin: "<<Xbins[bin]<<" to "<<Xbins[bin+1]<<" DiffXsec: "<<DiffXSecVec[channelType][bin]<<" StatError: "<<DiffXSecStatErrorVec[channelType][bin]<<" GenDiffXsec: "<<GenDiffXSecVec[channelType][bin]<<endl;
        ResultsFile<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[bin]<<" bin: "<<Xbins[bin]<<" to "<<Xbins[bin+1]<<" DiffXsec: "<<DiffXSecVec[channelType][bin]<<" StatError: "<<DiffXSecStatErrorVec[channelType][bin]<<" GenDiffXsec: "<<GenDiffXSecVec[channelType][bin]<<endl;
    }
    ResultsFile.close();
    ResultsLateX.close();
}

void Plotter::PlotDiffXSec(TString Channel){

    setDataSet(Channel,"Nominal");
    if (!fillHisto()) return;

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
    if ( drawSmoothMadgraph ) {
        for ( unsigned int i = 0; i<XAxisbinCenters.size();i++ ) {
            binCenters[i] = XAxisbinCenters[i];
        }
    } else {
        for ( unsigned int i = 0; i<XAxisbinCenters.size();i++ ) {
            binCenters[i] = XAxisbins[i] + (XAxisbins[i+1]-XAxisbins[i])/2;
        }
    }

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

    TString ftemp = "preunfolded/Nominal/"+Channel+"/"+name+"_UnfoldingHistos.root";

    //    theDataHist =  fileReader->GetClone<TH1D>(ftemp, "aDataHist");
    //theBgrHist =  fileReader->GetClone<TH1D>(ftemp, "aBgrHist");
    //theTtBgrHist =  fileReader->GetClone<TH1D>(ftemp, "aTtBgrHist");
//     RecoPlot =  fileReader->GetClone<TH1D>(ftemp, "aRecHist");
    TH1 *GenPlotTheory =  fileReader->GetClone<TH1D>(ftemp, "aGenHist");
    genReco2d =  fileReader->GetClone<TH2D>(ftemp, "aRespHist");


    for (unsigned int i =0; i<hists.size(); i++){
      varhists[i]=hists[i].Rebin(bins,"varhists",Xbins);
      setStyle(varhists[i], i);
    }

    TH1* GenPlot = GenPlotTheory->Rebin(bins,"genplot",Xbins);

    THStack * stack=  new THStack("def", "def");
    TLegend *leg = getNewLegendpre();
    int legchange = 0;
    TH1 *varhistsPlotting[hists.size()];


    for(unsigned int i=0; i<hists.size() ; i++){ // prepare histos and leg
        setStyle(varhists[i], i);
        varhistsPlotting[i]=(TH1*)varhists[i]->Clone();
        if(legends.at(i) != "Data"){
            if((legends.at(i) == DYEntry) && channelType!=2){
                varhists[i]->Scale(DYScale[channelType]);
                varhistsPlotting[i]->Scale(DYScale[channelType]);
            }

            if(i!=(hists.size()-1)){
                if(legends.at(i)!=legends.at(i+1)){
                //cout<<legends.at(i)<<endl;
                varhistsPlotting[i]->SetLineColor(1);
                }
            }else{
                varhistsPlotting[i]->SetLineColor(1);
            }
            
            if(legends.at(i) != legends.at(i-1)){
                varhistsPlotting[i]->SetLineColor(1);
                stack->Add(varhistsPlotting[i]);
            }
            if(i > 1){
                if(legends.at(i) != legends.at(i-1)){
                    legchange = i;
                    if( (legends.at(i) == DYEntry) && DYScale[channelType]!= 1){
                        leg->AddEntry(varhistsPlotting[i], legends.at(i), "f");
                    }
                    else leg->AddEntry(varhistsPlotting[i], legends.at(i) ,"f");
                }
                else{varhistsPlotting[legchange]->Add(varhistsPlotting[i]);}
            }
        }
        else{ if(i==0) leg->AddEntry(varhistsPlotting[i], legends.at(i) ,"pe");}
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
//     double efficiencies[XAxisbinCenters.size()];
    init = false;

    for (unsigned int hist =0; hist<hists.size(); hist++){
        if(legends.at(hist) == "Data"){
            for (Int_t bin=0; bin<bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
                DataSum[bin]+=varhists[hist]->GetBinContent(bin+1);
	    }
        }
        else if((legends.at(hist) == "t#bar{t} Signal")&&init==false){
	  signalHist=hist;
          init=true;
            for (Int_t bin=0; bin<bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
//                 efficiencies[bin] = (RecoPlot->GetBinContent(bin+1)) / (GenPlot->GetBinContent(bin+1));
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

    TH1 *h_DiffXSec    = (TH1D*)varhists[0]->Clone();   
    h_DiffXSec->Reset();
    TH1 *h_GenDiffXSec = (TH1D*)varhists[0]->Clone();   h_GenDiffXSec->Reset();

    //The systematic array is filled in the order in which the Stack is filled
    TString sys_array[] = {"PDF_", "HAD_", "MATCH_", "MASS_", "SCALE_", "BTAG_PT_", "BTAG_ETA_", "BTAG_LJET_PT_", "BTAG_LJET_ETA_", "KIN_", "LEPT_", "TRIG_", "BG_", "DY_", "PU_", "JER_", "JES_"};//For the time being uintil all systematics are finalished
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
                if(vec_systematic.at(syst).Contains("JER_") || vec_systematic.at(syst).Contains("JES_") || vec_systematic.at(syst).Contains("PU_") ||
                   vec_systematic.at(syst).Contains("DY_") || vec_systematic.at(syst).Contains("BG_") || vec_systematic.at(syst).Contains("TRIG_") ||
                   vec_systematic.at(syst).Contains("LEPT_") || vec_systematic.at(syst).Contains("BTAG_") || vec_systematic.at(syst).Contains("KIN_")){
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
            if (vec_systematic.at(systs) == "BTAG_ETA_" || vec_systematic.at(systs) == "BTAG_LJET_ETA_") {continue;}//Skip the BTAG_ETA systematic because it's added in quadrature to BTAG_PT
            TH1D* systtemp = (TH1D*)varhists[0]->Clone();
            systtemp->Reset();
            double TotalSyst=0.0, TotalSqSyst=0.0;
            double AvgSyst= 0.0, SqAvgSys=0.0;
            
            for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
                if(vec_systematic.at(systs) == "BTAG_PT_" || vec_systematic.at(systs) == "BTAG_LJET_PT_"){
                    DiffXSecSysErrorbySysPlot[bin][systs]= TMath::Sqrt(
                        (DiffXSecSysErrorbySysPlot[bin][systs]*DiffXSecSysErrorbySysPlot[bin][systs])
                        +(DiffXSecSysErrorbySysPlot[bin][systs+1]*DiffXSecSysErrorbySysPlot[bin][systs+1])
                    );
                }

                systtemp->SetBinContent(bin+1,(DiffXSecSysErrorbySysPlot[bin][systs]*DiffXSecSysErrorbySysPlot[bin][systs]));
                if(bin==0){
                    fprintf(systfile, "%s", (vec_systematic.at(systs)+" ").Data());
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

    //GenPlotTheory->Scale(topxsec/(SignalEventswithWeight*GenPlotTheory->GetBinWidth(1)));
    GenPlot->Scale(topxsec/(SignalEventswithWeight*GenPlot->GetBinWidth(1)));
    if( (name.Contains("Lepton")||name.Contains("Top")||name.Contains("BJet")) && !name.Contains("Lead")){
      GenPlotTheory->Scale(1./2.);
    }
//     GenPlotTheory->Rebin(2);
    double genscale = 1./GenPlotTheory->Integral("width");
    GenPlotTheory->Scale(genscale);
    
    genscale = 1./ h_GenDiffXSec->Integral("width");
    h_GenDiffXSec->Scale(genscale);
    
    cout << "+++ Generated cross sections";
    for (int i = 1; i <= h_GenDiffXSec->GetNbinsX(); ++i) cout << "\n" << i << "\t" << h_GenDiffXSec->GetBinContent(i);
    cout << "\n\n";
    
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
    
    bool canDrawMCATNLO = true;
    
    if (drawNLOCurves && drawMCATNLO) {
        mcnlohist = GetNloCurve(newname,"MCATNLO");
        double mcnloscale = 1./mcnlohist->Integral("width");
        if (binned_theory==false) mcnlohist->Rebin(2);mcnlohist->Scale(0.5); //#####
        mcnlohist->Scale(mcnloscale);
        
        mcnlohistBinned = mcnlohist->Rebin(bins,"mcnloplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){
            mcnlohistBinned->SetBinContent(bin+1,mcnlohistBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/mcnlohist->GetBinWidth(1)));
        }
        mcnlohistBinned->Scale(1./mcnlohistBinned->Integral("width"));
        
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
        //else if(name.Contains("LeptonBJetMass")){mcnlohistnorm = GetNloCurve("Jets","Eta","MCatNLO");}
        
        //    if (binned_theory==false) mcnlohistnorm->Rebin(5);mcnlohistnorm->Scale(0.2);
        if (mcnlohistnorm) {
            mcnlohistnormBinned = mcnlohistnorm->Rebin(bins,"genBinHistNorm", Xbins);

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
            //    if (binned_theory==false) mcnlohistdown->Rebin(5);mcnlohistdown->Scale(0.2);
            mcnlohistdownBinned    = mcnlohistdown->Rebin(bins,"genBinHist", Xbins);

            for (Int_t bin=0; bin<bins; bin++){
                mcnlohistupBinned->SetBinContent(bin+1,mcnlohistupBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/mcnlohistup->GetBinWidth(1)));
                mcnlohistdownBinned->SetBinContent(bin+1,mcnlohistdownBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/mcnlohistdown->GetBinWidth(1)));
                mcnlohistnormBinned->SetBinContent(bin+1,mcnlohistnormBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/mcnlohistnorm->GetBinWidth(1)));
            }
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
        } else {
            std::cout << "\n*************************\nMC@NLO Curve not available!\n**********************\n";
            canDrawMCATNLO = false;
        }
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
    
    if(drawNLOCurves && drawKidonakis && (name.Contains("ToppT") || name.Contains("TopRapidity")) && !name.Contains("Lead")){
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
    if (drawNLOCurves && drawMCATNLO && canDrawMCATNLO) {
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
    if(drawNLOCurves && drawKidonakis && (name.Contains("TopRapidity") || name.Contains("ToppT")) && !name.Contains("Lead")){
        Kidoth1_Binned->SetLineWidth(2);
        Kidoth1_Binned->SetLineColor(kOrange+4);
        Kidoth1_Binned->SetLineStyle(2);
        Kidoth1_Binned->Draw("SAME][");
    }
    
    if (drawSmoothMadgraph){
        if( name.Contains("HypTTBarMass")){
            GenPlotTheory->Rebin(15);
            GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
            TH1D *SmoothMadgraph =(TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(10);
            SmoothMadgraph->Draw("SAME, L");
        }
        else if(name.Contains("HypBJetEtaLead") || name.Contains("HypBJetEtaNLead") || name.Contains("HypBJetpTLead") || name.Contains("HypLeptonBjetMass")){
            GenPlotTheory->Rebin(2);
            GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
            TH1D *SmoothMadgraph =(TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(10);
            SmoothMadgraph->Draw("SAME, L");
        }
        else if(name.Contains("HypLeptonEtaNLead")){
            GenPlotTheory->Rebin(2);
            GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
            TH1D *SmoothMadgraph =(TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(5);
            SmoothMadgraph->Draw("SAME, L");
        }
        else if(name.Contains("HypLeptonEta") || name.Contains("HypTopRapidity") || 
                name.Contains("HypBJetEta") || name.Contains("HypBJetpTNLead") || name.Contains("HypTTBarRapidity")){
            TH1D *SmoothMadgraph =(TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(10);
            SmoothMadgraph->Draw("SAME, L");
        }
        else if (name.Contains("HypToppT") || name.Contains("HypLeptonBjetMass")){
            GenPlotTheory->Rebin(6);
            GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
            TH1D *SmoothMadgraph =(TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(10);
            SmoothMadgraph->Draw("SAME, L");
        }
        else if(name.Contains("HypLeptonpTLead")){
            GenPlotTheory->Rebin(2);
            GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
            TH1D *SmoothMadgraph =(TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(4);
            SmoothMadgraph->Draw("SAME, L");
        }
        else if(name.Contains("HypLeptonpTNLead")){
            GenPlotTheory->Rebin(2);
            GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
            GenPlotTheory->Draw("same,c");
        }
        else if(name.Contains("HypTTBarpT")){
            GenPlotTheory->Rebin(4);
            GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
            GenPlotTheory->Draw("same,c");
        }
        else if(name.Contains("HypLLBarMass")){
            GenPlotTheory->Rebin(4);
            GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
            TH1D *SmoothMadgraph =(TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(2);
            SmoothMadgraph->Draw("same,l");
        }
        else {GenPlotTheory->Draw("SAME,C");} //### 150512 ### 
    }

    h_GenDiffXSec->Draw("SAME"); //### 150512 ###
    DrawCMSLabels(1, 8);
    DrawDecayChLabel(channelLabel[channelType]);
    TLegend leg2 = *getNewLegend();
    leg2.AddEntry(h_DiffXSec, "Data", "p");
    leg2.AddEntry(GenPlotTheory, "MadGraph","l");
    if (drawNLOCurves) {
        if (drawMCATNLO && canDrawMCATNLO && mcnlohistup->GetEntries() && mcnlohistdown->GetEntries())   leg2.AddEntry(mcatnloBand,      "MC@NLO",  "fl");
        else if (drawMCATNLO && mcnlohist->GetEntries()) leg2.AddEntry(mcnlohist,      "MC@NLO",  "l");
        if (drawPOWHEG && powheghist->GetEntries())                                    leg2.AddEntry(powheghistBinned, "POWHEG",  "l");
        if (drawMadSpinCorr && spincorrhist->GetEntries())                             leg2.AddEntry(spincorrhistBinned, "MadGraph SC",  "l");
        if (drawKidonakis && !name.Contains("Lead") && (name.Contains("ToppT") || name.Contains("TopRapidity"))) leg2.AddEntry(Kidoth1_Binned,   "Approx. NNLO",  "l");
    }
    leg2.SetFillStyle(0);
    leg2.SetBorderSize(0);
    leg2.SetX1NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength()-0.30);
    leg2.SetY1NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength()-0.05*(double)leg2.GetNRows());
    leg2.SetX2NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength());
    leg2.SetY2NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength());
    leg2.SetTextSize(0.04);
    leg2.Draw("same");
    if (drawNLOCurves && drawKidonakis &&  (name.Contains("ToppT") || name.Contains("TopRapidity")) && !name.Contains("Lead")){
        DrawLabel("(arXiv:1210.7813)", leg2.GetX1NDC()+0.06, leg2.GetY1NDC()-0.025, leg2.GetX2NDC(), leg2.GetY1NDC(), 12, 0.025);
    }
    h_GenDiffXSec->Draw("SAME");
    gStyle->SetEndErrorSize(10);
    tga_DiffXSecPlot->Draw("p, SAME");
    tga_DiffXSecPlotwithSys->Draw("p, SAME, Z");
    gPad->RedrawAxis();
    c->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/DiffXS_"+name+".eps"); 
    TFile out_source(outpathPlots+subfolderChannel+subfolderSpecial+"/DiffXS_"+name+"_source.root", "RECREATE");
    c->Write("canvas");
    tga_DiffXSecPlot->Write("data_staterror_only");
    tga_DiffXSecPlotwithSys->Write("data");
    h_GenDiffXSec->Write("mc");
    out_source.Close();
    //c->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/DiffXS_"+name+".C"); 
    c->Clear();
    delete c;
    gStyle->SetEndErrorSize(0);

    TH1D *dataHisto = TGraphToTH1(tga_DiffXSecPlotwithSys, h_GenDiffXSec);
    cout<<"\nFor variable "<<name<<" in channel "<<subfolderChannel<< " the Chi2/ndof= "<<dataHisto->Chi2Test(h_GenDiffXSec, "WW,CHI2/NDF")<<"\n\n";
    
    PrintResultTotxtFile(Channel, binCenters, tga_DiffXSecPlot, tga_DiffXSecPlotwithSys);

    TCanvas * c1 = new TCanvas("DiffXS","DiffXS");
    TList* l = stack->GetHists();
    TH1D* stacksum = (TH1D*) l->At(0)->Clone();
    for (int i = 1; i < l->GetEntries(); ++i) {
      stacksum->Add((TH1D*)l->At(i));
    } 
    for(unsigned int i=1; i<hists.size() ; i++){ // sum all data plots to first histogram
      if(legends.at(i) == legends.at(0)){
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
    DrawCMSLabels(1, 8);
    DrawDecayChLabel(channelLabel[channelType]);    
    leg->Draw("SAME");
    gPad->RedrawAxis();
    //    TFile *f1 = new TFile("KinFitPlots.root","UPDATE");
    //stacksum->Write(name+"_"+channel+"_MC");
    //varhists[0]->Write(name+"_"+channel+"_Data");
    //f1->Close();
    c1->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/preunfolded_"+name+".eps");
    //c1->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/preunfolded_"+name+".C");
    //Draw
    TFile out_root(outpathPlots+subfolderChannel+subfolderSpecial+"/preunfolded_"+name+"_source.root", "RECREATE");
    varhists[0]->Write(name+"_data");
//     ->Write(name+"_signalmc");
    stacksum->Write(name+"_allmc");
    c1->Write(name + "_canvas");
    out_root.Close();
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

    TString filename;
    if (Generator == "MCATNLO") {
        filename = "_ttbarsignalplustau_mcatnlo.root";
    } else if (Generator == "POWHEG") {
        filename = "_ttbarsignalplustau_powheg.root";
    } else if (Generator == "SPINCORR") {
        filename = "_ttbarsignalplustau_FullLeptMadgraphWithSpinCorrelation.root";
    } else {
        std::cerr << "Unknown Generator!\n";
        std::exit(2);
    }
    
    const static std::vector<TString> channelName {"ee", "mumu", "emu"};
    std::vector<TString> files;
    assert(channelType >= 0);
    assert(channelType <= 3);
    for (int i = 0; i <= 2; ++i) {
        if (channelType == i || channelType == 3) {
            files.push_back("selectionRoot/"+Generator+"/"+channelName.at(i)+"/"+channelName.at(i) + filename);
            cout << "Getting NLO curve from: " << files.at(files.size()-1) << endl;
        }
    }
    
    TString histname("VisGen"+NewName);
    TH1* hist = fileReader->GetClone<TH1>(files.at(0), histname);
    for (size_t i = 1; i < files.size(); ++i) {
        hist->Add(fileReader->Get<TH1>(files.at(i), histname));
    }

    if (!NewName.Contains("Lead") 
        && (NewName.Contains("Lepton") || NewName.Contains("Top") || NewName.Contains("BJet")))
    {
        //loop over anti-particle histograms and add them
        TString antiName("VisGenAnti"+NewName);
        for (const auto& file : files) {
            hist->Add(fileReader->Get<TH1>(file, antiName));
        }
    }
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


TH1D* Plotter::TGraphToTH1(TGraphAsymmErrors *gr, TH1 *mc){
    
    if ( !gr || !mc){
        cout<<"Not TGraph"<<std::endl;;
        exit(34);
    }
    
    vector<double> xBinLow, xBinHigh;
    
    for (int i=1; i<=(int) mc->GetNbinsX(); ++i){
        xBinLow.push_back(mc->GetXaxis()->GetBinLowEdge(i));
        xBinHigh.push_back(mc->GetXaxis()->GetBinUpEdge(i));
    }

    vector<double> xpoints, ypoints;
    vector<double> yLError, yHError;
    vector<double> xLError, xHError;
    xpoints.clear(); xLError.clear(); xHError.clear();
    ypoints.clear(); yLError.clear(); yHError.clear();
    
    for ( int i=0; i<(int) gr->GetN(); ++i ) {
        xpoints.push_back(gr->GetX()[i]);
        xHError.push_back(gr->GetErrorXhigh(i));
        xLError.push_back(gr->GetErrorXlow(i));
        
        ypoints.push_back(gr->GetY()[i]);
        yHError.push_back(gr->GetErrorYhigh(i));
        yLError.push_back(gr->GetErrorYlow(i));
    }
    
    TH1D * dataHisto = (TH1D*)mc ->Clone();
    
    for ( int i=0; i<(int) dataHisto->GetNbinsX(); ++i ) {
        
        if ( TMath::Abs(yHError.at(i)-yLError.at(i)) > 1e-6 ){
            std::cout<<"Up and down error differents. Cannot proceed with the TGraphAsymmErrors -> TH1 conversion"<<std::endl;
            exit(22);
        }
        dataHisto->SetBinContent(i, ypoints.at(i));
        dataHisto->SetBinError(i, yHError.at(i));
    }

    return dataHisto;
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

double Plotter::CalcLumiWeight(const TString& WhichSample){
    if (WhichSample.Contains("run")) return 1;
    double lumiWeight=0;
    if(WhichSample!=""){
        double XSection = SampleXSection(WhichSample);
        if(XSection <= 0.){
            cout<<"Sample XSection is <0. Can't calculate luminosity weight!! returning"<<endl;
            return 0;
        }
        
        //From 'filename' get the number of weighted (MC weights) event processed.
        const TH1 *h_NrOfEvts = fileReader->Get<TH1>(WhichSample, "weightedEvents");
        double NrOfEvts = h_NrOfEvts->GetBinContent(1);
        lumiWeight = lumi*XSection/NrOfEvts;
    }
    
    if (lumiWeight == 0) {
        cout << WhichSample << " has lumi weight 0\n";
    }

    return lumiWeight;
}


double Plotter::SampleXSection(const TString& filename){
    
    //MC cross sections taken from:
    //  https://twiki.cern.ch/twiki/bin/view/CMS/StandardModelCrossSectionsat8TeV
    //  AN-12/194    AN-12/228
    
    if(filename.Contains("run"))              {return 1;}
    else if(filename.Contains("FullLept"))    {return topxsec * 0.1049;}
    else if(filename.Contains("SemiLept"))    {return topxsec * 0.4380;}
    else if(filename.Contains("Hadronic"))    {return topxsec * 0.4570;}
    else if(filename.Contains("ttbar"))       {return topxsec;}
    else if(filename.Contains("single"))      {return 11.1;}
    else if(filename.Contains("ww"))          {return 54.838;}
    else if(filename.Contains("wz"))          {return 33.21;}
    else if(filename.Contains("zz"))          {return 17.654;}
    else if(filename.Contains("1050"))        {return 860.5;}
    else if(filename.Contains("50inf"))       {return 3532.8;}
    else if(filename.Contains("wtolnu"))      {return 36257.2;}
    else if(filename.Contains("qcdmu15"))     {return 3.640E8*3.7E-4;}
    else if(filename.Contains("qcdmu2030"))   {return 2.870E8*6.500E-3;}
    else if(filename.Contains("qcdmu3050"))   {return 6.609E7*12.20E-3;}
    else if(filename.Contains("qcdmu5080"))   {return 8.802E6*21.80E-3;}
    else if(filename.Contains("qcdmu80120"))  {return 1.024E6*39.50E-3;}
    else if(filename.Contains("qcdmu120170")) {return 1.578E5*47.30E-3;}
    else if(filename.Contains("qcdem2030"))   {return 2.886E8*10.10E-3;}
    else if(filename.Contains("qcdem3080"))   {return 7.433E7*62.10E-3;}
    else if(filename.Contains("qcdem80170"))  {return 1.191E6*153.9E-3;}
    else if(filename.Contains("qcdbcem2030")) {return 2.886E8*5.800E-4;}
    else if(filename.Contains("qcdbcem3080")) {return 7.424E7*2.250E-3;}
    else if(filename.Contains("qcdbcem80170")){return 1.191E6*10.90E-3;}
    
    return -1;
}

// Draw label for Decay Channel in upper left corner of plot
void Plotter::DrawDecayChLabel(TString decaychannel, double textSize) {

    TPaveText *decch = new TPaveText();

    decch->AddText(decaychannel);

    decch->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
    decch->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );
    decch->SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
    decch->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength()        );

    decch->SetFillStyle(0);
    decch->SetBorderSize(0);
    if (textSize!=0) decch->SetTextSize(textSize);
    decch->SetTextAlign(12);
    decch->Draw("same");
}

// Draw official labels (CMS Preliminary, luminosity and CM energy) above plot
void Plotter::DrawCMSLabels(int cmsprelim, double energy, double textSize) {

    const char *text;
    if(cmsprelim ==2 ) {//Private work for PhDs students
        text = "Private Work, %2.1f fb^{-1} at #sqrt{s} = %2.f TeV";
    } else if (cmsprelim==1) {//CMS preliminary label
        text = "CMS Preliminary, %2.1f fb^{-1} at #sqrt{s} = %2.f TeV";
    } else {//CMS label
        text = "CMS, %2.1f fb^{-1} at #sqrt{s} = %2.f TeV";
    }
    
    TPaveText *label = new TPaveText();
    label->SetX1NDC(gStyle->GetPadLeftMargin());
    label->SetY1NDC(1.0-gStyle->GetPadTopMargin());
    label->SetX2NDC(1.0-gStyle->GetPadRightMargin());
    label->SetY2NDC(1.0);
    label->SetTextFont(42);
    label->AddText(Form(text, lumi/1000, energy));
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    if (textSize!=0) label->SetTextSize(textSize);
    label->SetTextAlign(32);
    label->Draw("same");
}


void Plotter::PrintResultTotxtFile (TString channel, double binCenters[], TGraphAsymmErrors *tga_DiffXSecPlot, TGraphAsymmErrors *tga_DiffXSecPlotwithSys){
    
    ofstream SavingFile;
    string filename = string("Plots/"+channel+"/"+name+"LaTeX.txt");
    SavingFile.open(filename.c_str(), ios_base::out);
    
    SavingFile<<"Variable: "<<name<<"   Channel: "<<channelLabel.at(channelType)<<endl;
    SavingFile<<"BinCenter & LowXbinEdge  &  HighXbinEdge  &   DiffXSec  &  StatError(\\%)  & SystError(\\%)  & TotalError(\\%) \\\\"<<endl;
    SavingFile<<"\\hline"<<endl;
    for(int i=0; i<(int)tga_DiffXSecPlot->GetN(); i++){
        double DiffXSec=tga_DiffXSecPlot->GetY()[i];
        double RelStatErr=0, RelSysErr=0, RelTotErr=0;
        if(DiffXSec!=0.0){
            RelStatErr = 100*(tga_DiffXSecPlot->GetErrorY(i))/DiffXSec;
            RelTotErr  = 100*(tga_DiffXSecPlotwithSys->GetErrorY(i))/DiffXSec;
            if(RelTotErr>=RelStatErr) RelSysErr = TMath::Sqrt(RelTotErr*RelTotErr - RelStatErr*RelStatErr);
        }
        SavingFile<<"$"<<setprecision(3)<<binCenters[i]<<"$ & $"<<XAxisbins.at(i)<<"$ to $"<<setprecision(3)<<XAxisbins.at(i+1)<<"$   &  ";
        SavingFile<<setprecision(5)<<DiffXSec<<"  &   "<<setprecision(3)<<RelStatErr<<" &    "<<setprecision(3)<<RelSysErr<<" &    ";
        SavingFile<<setprecision(3)<<RelTotErr<<" \\\\"<<endl;
    }
    SavingFile.close();
}



void Plotter::GetDiffToNominal(TString Channel, TString Systematic, TString Variable){
    
    //This function will read the unfolding result for Nominal and Variation1 (eg: MASS_UP) for a given variable
    // and it will return the relative error: abs(MASS_UP-Nominal)/Nominal
    

    ifstream NominalFile, SystematicFile;
    string nominalfilename = string ("UnfoldingResults/Nominal/"+Channel+"/"+Variable+"Results.txt");
    string systematicfilename = string ("UnfoldingResults/"+Systematic+"/"+Channel+"/"+Variable+"Results.txt");

    NominalFile.open(nominalfilename);
    SystematicFile.open(systematicfilename);
    
    
    vector <double> BinCenters, BinLowEdge, BinHigEdge;
    vector <double> NominalValue, SystematicValue;
    vector <double> RelativeError;
    TString Dummy = "";
    
    
    if (!NominalFile.is_open() || !SystematicFile.is_open()){
        cout<<"The input file cannot be opened. Exiting!!"<<endl;
        exit(443);}
    for (Int_t bin=0; bin<bins; bin++){
        double XAxisbinCenters_Nom = 0, XLowEdge_Nom = 0, XHigEdge_Nom = 0;
        double XAxisbinCenters_Sys = 0, XLowEdge_Sys = 0, XHigEdge_Sys = 0;
        double CentralValue_Nom = 0, StatError_Nom = 0, GenXSec_Nom = 0;
        double CentralValue_Sys = 0, StatError_Sys = 0, GenXSec_Sys = 0;
        NominalFile>>Dummy>>XAxisbinCenters_Nom>>Dummy>>XLowEdge_Nom>>Dummy>>XHigEdge_Nom>>Dummy>>CentralValue_Nom>>Dummy>>StatError_Nom>>Dummy>>GenXSec_Nom;
        SystematicFile>>Dummy>>XAxisbinCenters_Sys>>Dummy>>XLowEdge_Sys>>Dummy>>XHigEdge_Sys>>Dummy>>CentralValue_Sys>>Dummy>>StatError_Sys>>Dummy>>GenXSec_Sys;
        if( abs(XAxisbinCenters_Nom-XAxisbinCenters_Sys)>1e3 || 
            abs(XLowEdge_Nom-XLowEdge_Sys)>1e3 || abs(XHigEdge_Nom-XHigEdge_Sys)>1e3){
                cout<<"You are trying to read Unfolding results containing different bin centers and/or bin ranges."<<endl;
                cout<<"This is not a fatal error but you should check it"<<endl;
                cout<<"Variable '"<<Variable<<"' in channel '"<<Channel<<"'"<<endl;
            }
        BinCenters.push_back(XAxisbinCenters_Nom);
        BinLowEdge.push_back(XLowEdge_Nom);
        BinHigEdge.push_back(XHigEdge_Nom);
        NominalValue.push_back(CentralValue_Nom);
        SystematicValue.push_back(CentralValue_Sys);

        if(CentralValue_Nom != 0){
            RelativeError.push_back(abs(CentralValue_Nom-CentralValue_Sys)/CentralValue_Nom);
        }
    }
    NominalFile.close();
    SystematicFile.close();

    ofstream SystematicRelError;
    SystematicRelError.open(systematicfilename, ios::trunc);
    if(!SystematicRelError.is_open()){
        cout<<"The output file cannot be opened. Exiting!!"<<endl;
        exit(444);
    }
    for (int bin = 0; bin<(int)SystematicValue.size(); bin++){
        SystematicRelError<<"XAxisbinCenters[bin]: "<<BinCenters.at(bin)<<" bin: "<<BinLowEdge.at(bin)<<" to "<<BinHigEdge.at(bin)<<" SystematicRelError: "<<RelativeError.at(bin)<<endl;
    }
    SystematicRelError.close();

}


void Plotter::CalcUpDownDifference( TString Channel, TString Syst_Up, TString Syst_Down, TString Variable){
    
    //Function to get the error of a certain systematic: Err_Up - Err_Down
    
//This function will read the unfolding result for Nominal and Variation1 (eg: MASS_UP) for a given variable
    // and it will return the relative error: abs(MASS_UP-Nominal)/Nominal
    
    ofstream DownRelError;
    ifstream UpFile, DownFile;
    string upfilename = string ("UnfoldingResults/"+Syst_Up+"/"+Channel+"/"+Variable+"Results.txt");
    string downfilename = string ("UnfoldingResults/"+Syst_Down+"/"+Channel+"/"+Variable+"Results.txt");
    
    cout<<"\n\n\n"<<upfilename<<endl;
    cout<<downfilename<<endl;
    
    UpFile.open(upfilename);
    DownFile.open(downfilename);
    
    vector <double> BinCenters, BinLowEdge, BinHigEdge;
    vector <double> UpValue, DownValue;
    vector <double> RelativeError;
    TString Dummy = "";
    
    if (!UpFile.is_open() || !DownFile.is_open()){
        cout<<"The input file cannot be opened. Exiting!!"<<endl;
        exit(433);}
    for (Int_t bin=0; bin<bins; bin++){
        double XAxisbinCenters_Up = 0, XLowEdge_Up = 0, XHigEdge_Up = 0;
        double XAxisbinCenters_Down = 0, XLowEdge_Down = 0, XHigEdge_Down = 0;
        double CentralValue_Up = 0;
        double CentralValue_Down = 0;
        UpFile>>Dummy>>XAxisbinCenters_Up>>Dummy>>XLowEdge_Up>>Dummy>>XHigEdge_Up>>Dummy>>CentralValue_Up;
        DownFile>>Dummy>>XAxisbinCenters_Down>>Dummy>>XLowEdge_Down>>Dummy>>XHigEdge_Down>>Dummy>>CentralValue_Down;

        BinCenters.push_back(XAxisbinCenters_Up);
        BinLowEdge.push_back(XLowEdge_Up);
        BinHigEdge.push_back(XHigEdge_Up);
        UpValue.push_back(CentralValue_Up);
        DownValue.push_back(CentralValue_Down);

        if(CentralValue_Up != 0){
            RelativeError.push_back(abs(CentralValue_Up-CentralValue_Down));
        }
    }
    UpFile.close();
    DownFile.close();

    ofstream SystematicRelError;
    string systematicfilename = string ("UnfoldingResults/HAD_/"+Channel+"/"+Variable+"Results.txt");
    gSystem->mkdir("UnfoldingResults/HAD_/"+Channel, true);
    SystematicRelError.open(systematicfilename, ios::trunc);
    if(!SystematicRelError.is_open()){
        cout<<"The output file cannot be opened. Exiting!!"<<endl;
        exit(434);
    }
    for (int bin = 0; bin<(int)RelativeError.size(); bin++){
        SystematicRelError<<"XAxisbinCenters[bin]: "<<BinCenters.at(bin)<<" bin: "<<BinLowEdge.at(bin)<<" to "<<BinHigEdge.at(bin)<<" SystematicRelError: "<<RelativeError.at(bin)<<endl;
    }
    SystematicRelError.close();

}
