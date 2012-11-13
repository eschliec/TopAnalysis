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

void load_Analysis(TString validFilenamePattern, TString systematic){
   
    int filecounter = 0;

    ifstream infile ("selectionList.txt");
    TString filename;
    
    Analysis *selector = new Analysis();
    PUReweighter *pu = new PUReweighter();
    pu->setMCDistrSum12("S10");
    std::string pu_path(getenv("CMSSW_BASE"));
    if (systematic == "") {
        pu_path.append("/src/TopAnalysis/TopUtils/data/Data_PUDist_12fb.root");
    } else if (systematic == "PU_UP") {
        pu_path.append("/src/TopAnalysis/TopUtils/data/Data_PUDist_12fb_sysUp.root");
    } else if (systematic == "PU_DOWN") {
        pu_path.append("/src/TopAnalysis/TopUtils/data/Data_PUDist_12fb_sysDown.root");
    } else {
        cerr << "Systematics unknown!\n"; exit(3);
    }
    
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
        TObjString *systematics_from_file = dynamic_cast<TObjString*>(file.Get("writeNTuple/systematicsName"));
        TObjString *samplename = dynamic_cast<TObjString*>(file.Get("writeNTuple/sampleName"));
        TObjString *o_isSignal = dynamic_cast<TObjString*>(file.Get("writeNTuple/isSignal"));
        TObjString *o_isMC = dynamic_cast<TObjString*>(file.Get("writeNTuple/isMC"));
        TH1* weightedEvents = dynamic_cast<TH1*>(file.Get("EventsBeforeSelection/weightedEvents"));
        if (!channel || !systematics_from_file || !o_isSignal || !o_isMC || !samplename) { 
            std::cout << "Error: info about sample missing!" << std::endl; 
            return;  
        }
        bool isSignal = o_isSignal->GetString() == "1";
        bool isMC = o_isMC->GetString() == "1";
        TString btagFile = "BTagEff/Nominal/" + channel->GetString() + "/" 
            + channel->GetString() + "_ttbarsignalplustau.root";
        
        selector->SetBTagFile(btagFile);
        selector->SetChannel(channel->GetString());
        selector->SetSignal(isSignal);
        selector->SetMC(isMC);
        if (systematic == "") {
            selector->SetSystematic(systematics_from_file->GetString());
        } else {
            selector->SetSystematic(systematic);
        }
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

void syntaxError() {
    std::cerr << 
        "\nload_Analysis: Invalid syntax!\n---------------------------------\n" <<
        "Valid Parameters:\n" <<
        "-f pattern  -->  only process filenames containing the pattern\n" <<
        "-s [ PU_UP | PU_DOWN ]  -->  run with a systematic that runs on the normal ntuples\n" <<
        "\n";
    exit(1);
}

int main(int argc, char* const argv[]) {
    char opt;
    TString validFilenamePattern{""};
    TString syst{""};
    while ((opt = getopt(argc, argv, "f:s:")) != -1) {
        if (opt == 'f') {
            validFilenamePattern = optarg;
        } else if (opt == 's') {
            syst = optarg;
        } else {
            syntaxError();
        }
    }
    if (optind < argc) syntaxError();
//     TProof* p = TProof::Open(""); // not before ROOT 5.34
    load_Analysis(validFilenamePattern, syst);
//     delete p;
}
