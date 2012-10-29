#include <iostream>
#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TProof.h>
#include <TSelector.h>
#include <TObjString.h>
#include "Analysis.h"
#include "PUReweighter.h"

void load_Analysis(TString validFilenamePattern){
   
    int filecounter = 0;

    ifstream infile ("selectionList.txt");
    TString filename;
    
    Analysis *selector = new Analysis();
    PUReweighter *pu = new PUReweighter();
    pu->setMCDistrSum12("S10");
    std::string pu_path(getenv("CMSSW_BASE"));
//     pu_path.append("/src/TopAnalysis/TopUtils/data/Data_PUDist_sysNo_69400_2012ABReReco.root");
    pu_path.append("/src/TopAnalysis/TopUtils/data/PU_Data_2012_5fbinv.root");
    pu->setDataTruePUInput(pu_path.c_str());
    selector->SetPUReweighter(pu);
    

    while(!infile.eof()){
        infile >> filename;
        if (filename == "" || filename[0] == '#') continue; //empty line? --> skip
        if (!filename.Contains(validFilenamePattern)) continue;

        //channel selection for later BTag eff
        std::cout << std::endl;
        std::cout << "PROCESSING File " << ++filecounter << " ("<< filename << ") from selectionList.txt" << std::endl;
        std::cout << std::endl;

        TFile file(filename);
        if (file.IsZombie()) { std::cerr << "Cannot open " << filename << std::endl; return; }

        TObjString *channel = dynamic_cast<TObjString*>(file.Get("writeNTuple/channelName"));
        TObjString *systematics = dynamic_cast<TObjString*>(file.Get("writeNTuple/systematicsName"));
        TObjString *samplename = dynamic_cast<TObjString*>(file.Get("writeNTuple/sampleName"));
        TObjString *o_isSignal = dynamic_cast<TObjString*>(file.Get("writeNTuple/isSignal"));
        TObjString *o_isMC = dynamic_cast<TObjString*>(file.Get("writeNTuple/isMC"));
        TH1* weightedEvents = dynamic_cast<TH1*>(file.Get("EventsBeforeSelection/weightedEvents"));
        if (!channel || !systematics || !o_isSignal || !o_isMC || !samplename) { 
            std::cout << "Error: info about sample missing!" << std::endl; 
            return;  
        }
        bool isSignal = o_isSignal->GetString() == "1";
        bool isMC = o_isMC->GetString() == "1";
        TString btagFile = "BTagEff/" + channel->GetString() + "/ttbarsignalplustau.root";
        
        selector->SetBTagFile(btagFile);
        selector->SetChannel(channel->GetString());
        selector->SetSignal(isSignal);
        selector->SetMC(isMC);
        selector->SetSystematic(systematics->GetString());
        selector->SetWeightedEvents(weightedEvents);
        selector->SetSamplename(samplename->GetString());
        selector->SetOutputfilename(filename);
        selector->SetRunViaTau(0);

        TTree *tree = dynamic_cast<TTree*>(file.Get("writeNTuple/NTuple"));
        if (! tree) { std::cerr << "Error: Tree not found!"; return; }
        
        TChain chain("writeNTuple/NTuple");
        chain.Add(filename);
//         chain.SetProof(); //will work from 5.34 onwards
        
        chain.Process(selector);
        
        if (isSignal) {
            selector->SetRunViaTau(1);
            filename.ReplaceAll("signalplustau", "bgviatau");
            selector->SetOutputfilename(filename);
            chain.Process(selector);
        }
        
        file.Close();
    }
}

int main(int argc, const char* argv[]) {
    TString validFilenamePattern = argc > 1 ? argv[1] : "";

//     TProof* p = TProof::Open(""); // not before ROOT 5.34
    load_Analysis(validFilenamePattern);
//     delete p;
}

