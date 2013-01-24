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
#include "CommandLineParameters.hh"

void load_Analysis(TString validFilenamePattern, TString givenChannel, TString systematic){
   
    ifstream infile ("selectionList.txt");
    if (!infile.good()) { 
        cerr << "Error! Please check the selectionList.txt file!\n" << endl; 
        exit(773); 
    }
    
    Analysis *selector = new Analysis();
    PUReweighter *pu = new PUReweighter();
    pu->setMCDistrSum12("S10");
    std::string pu_path(getenv("CMSSW_BASE"));
    if (systematic == "PU_UP") {
        pu_path.append("/src/TopAnalysis/TopUtils/data/Data_PUDist_12fb_sysUp.root");
        cout << "using pilup-up distribution\n";
    } else if (systematic == "PU_DOWN") {
        pu_path.append("/src/TopAnalysis/TopUtils/data/Data_PUDist_12fb_sysDown.root");
        cout << "using pilup-down distribution\n";
    } else {
        pu_path.append("/src/TopAnalysis/TopUtils/data/Data_PUDist_12fb.root");
        if (systematic != "") {
            cout << "Using Nominal PU distribution for " << systematic << " systematic!\n";
        }
    }
    
    pu->setDataTruePUInput(pu_path.c_str());
    selector->SetPUReweighter(pu);

    int filecounter = 0;
    while(!infile.eof()){
        TString filename;
        infile >> filename;
        if (filename == "" || filename[0] == '#') continue; //empty line? --> skip
        if (!filename.Contains(validFilenamePattern)) continue;

        //channel selection for later BTag eff
        std::cout << std::endl;
        std::cout << "PROCESSING File " << ++filecounter << " ("<< filename << ") from selectionList.txt" << std::endl;
        std::cout << std::endl;

        TFile file(filename);
        if (file.IsZombie()) { std::cerr << "Cannot open " << filename << std::endl; return; }

        TObjString *channel_file = dynamic_cast<TObjString*>(file.Get("writeNTuple/channelName"));
        TObjString *systematics_from_file = dynamic_cast<TObjString*>(file.Get("writeNTuple/systematicsName"));
        TObjString *samplename = dynamic_cast<TObjString*>(file.Get("writeNTuple/sampleName"));
        TObjString *o_isSignal = dynamic_cast<TObjString*>(file.Get("writeNTuple/isSignal"));
        TObjString *o_isMC = dynamic_cast<TObjString*>(file.Get("writeNTuple/isMC"));
        TH1* weightedEvents = dynamic_cast<TH1*>(file.Get("EventsBeforeSelection/weightedEvents"));
        if (!channel_file || !systematics_from_file || !o_isSignal || !o_isMC || !samplename) { 
            std::cout << "Error: info about sample missing!" << std::endl; 
            return;  
        }
        bool isSignal = o_isSignal->GetString() == "1";
        bool isMC = o_isMC->GetString() == "1";
        
        if (!isMC && systematic != "") {
            cout << "Sample is DATA, so not running again for systematic variation\n";
            continue;
        }
        
        if (systematic == "PDF" && (!isSignal || !(systematics_from_file->GetString() == "Nominal"))) {
            cout << "Skipping file: is not signal or not nominal -- and running PDFs\n";
            continue;
        }
        
        //determine the channels to run over
        std::vector<TString> channels;
        //is the "mode" (=channel) given in the file?
        if (channel_file->GetString() != "") {
            channels.push_back(channel_file->GetString());
        } else {
            if (givenChannel != "") {
                channels.push_back(givenChannel);
            } else {
                channels.push_back("emu");
                channels.push_back("ee");
                channels.push_back("mumu");
            }
        }
        
        for (const auto& channel : channels) {
            TString btagFile = "BTagEff/Nominal/" + channel + "/" 
                + channel + "_ttbarsignalplustau.root";
            
            selector->SetBTagFile(btagFile);
            selector->SetChannel(channel);
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
            if (! tree) { std::cerr << "Error: Tree not found!\n"; exit(854); }
            
            TChain chain("writeNTuple/NTuple");
            chain.Add(filename);
            // chain.SetProof(); //will work from 5.34 onwards
            
            if (systematic == "PDF") {
                TH1* pdfWeights = dynamic_cast<TH1*>(file.Get("EventsBeforeSelection/pdfEventWeights"));
                if (!pdfWeights) { std::cerr << "Error: pdfEventWeights histo missing!\n"; exit(831); }
                for (int pdf_no = 1; pdfWeights->GetBinContent(pdf_no) > 0; ++pdf_no) {
                    TString pdfName("PDF_");
                    pdfName += (pdf_no+1)/2;
                    pdfName += (pdf_no % 2 ? "_UP" : "_DOWN");
                    selector->SetSystematic(pdfName);
                    weightedEvents->SetBinContent(1, pdfWeights->GetBinContent(pdf_no));
                    selector->SetWeightedEvents(weightedEvents);
                    selector->SetPDF(pdf_no);
                    chain.Process(selector);
                }
            } else {
                chain.Process(selector);
                if (isSignal) {
                    selector->SetRunViaTau(1);
                    filename.ReplaceAll("signalplustau", "bgviatau");
                    selector->SetOutputfilename(filename);
                    chain.Process(selector);
                }
            }
        }
        file.Close();
    }
}

void syntaxError() {
    std::cerr << 
        "\nload_Analysis: Invalid syntax!\n---------------------------------\n" <<
        "Valid Parameters:\n" <<
        "-f pattern\n" <<
        "   only process filenames containing the pattern\n" <<
        "-s [ PU_UP | PU_DOWN | TRIG_UP | TRIG_DOWN | BTAG_... ]\n" << 
        "   run with a systematic that runs on the normal ntuples\n" <<
        "-c channel (ee, emu, mumu)\n" <<
        "   provide this for all MC samples, automatically known for data" <<
        "   if no channel is provided, run over all three channels" <<
        "\n";
    exit(1);
}

int main(int argc, char** argv) {
    //CLParameter<std::string> opt_f("f", "Restrict to filename pattern", false, 1, 1);
    
    //CLAnalyser::interpretGlobal(argc, (char**)argv);
    
    char opt;
    TString validFilenamePattern;
    TString syst;
    TString channel;
    while ((opt = getopt(argc, argv, "f:s:c:")) != -1) {
        if (opt == 'f') {
            validFilenamePattern = optarg;
        } else if (opt == 's') {
            syst = optarg;
        } else if (opt == 'c') {
            if (channel != "" && channel != "ee" && channel != "emu" && channel != "mumu") {
                std::cerr << "invalid channel!\n"; exit(1);
            }
            channel = optarg;
        } else {
            syntaxError();
        }
    }
    if (optind < argc) syntaxError();
//     TProof* p = TProof::Open(""); // not before ROOT 5.34
    load_Analysis(validFilenamePattern, channel, syst);
//     delete p;
}
