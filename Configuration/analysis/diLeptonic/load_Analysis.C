#include <iostream>
#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TProof.h>
#include <TSelector.h>
#include <TObjString.h>
#include "Analysis.h"

void load_Analysis(){
   
    int filecounter = 0;

    ifstream infile ("selectionList.txt");
    TString filename;
    
    Analysis *selector = new Analysis(); //(Analysis*)TSelector::GetSelector("Analysis.C+");

    while(!infile.eof()){
        infile >> filename;
        if (filename == "" || filename[0] == '#') continue; //empty line? --> skip

        //channel selection for later BTag eff
        std::cout << std::endl;
        std::cout << "PROCESSING File " << ++filecounter << " ("<< filename << ") from selectionList.txt" << std::endl;
        std::cout << std::endl;

        TFile *f1 = TFile::Open(filename);
        if (!f1 || f1->IsZombie()) { std::cerr << "Cannot open " << filename << std::endl; return; }

        TObjString *channel = dynamic_cast<TObjString*>(f1->Get("writeNTuple/channelName"));
        TObjString *systematics = dynamic_cast<TObjString*>(f1->Get("writeNTuple/systematicsName"));
        TObjString *samplename = dynamic_cast<TObjString*>(f1->Get("writeNTuple/sampleName"));
        TObjString *o_isSignal = dynamic_cast<TObjString*>(f1->Get("writeNTuple/isSignal"));
        TObjString *o_isMC = dynamic_cast<TObjString*>(f1->Get("writeNTuple/isMC"));
        TH1* weightedEvents = dynamic_cast<TH1*>(f1->Get("EventsBeforeSelection/weightedEvents"));
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

        TTree *tree = dynamic_cast<TTree*>(f1->Get("writeNTuple/NTuple"));
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
        
        f1->Close();
        delete f1;
    }
}

int main() {
//     TProof* p = TProof::Open(""); // not before ROOT 5.34
    load_Analysis();
//     delete p;
}

