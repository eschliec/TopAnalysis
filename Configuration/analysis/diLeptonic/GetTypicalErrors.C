#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>

#include <TString.h>

using namespace std;


vector<TString> Channels(){
    
    vector<TString> channel {"ee", "emu", "mumu"};
    
    return channel;
}

vector<TString> Variables(){
    
    vector<TString> variables {"ToppTLead", "ToppTNLead",
                               "TopEtaLead", "TopEtaNLead",
                               "BJetpTLead", "BJetpTNLead",
                               "BJetEtaLead", "BJetEtaNLead",
                               "LeptonpTLead", "LeptonpTNLead", "LeptonEtaLead", "LeptonEtaNLead",
                               "TTBarpT", "TTBarMass", "LLBarpT", "LLBarMass", 
                               "LeptonBjetMass"
                               };
    return variables;
}

vector<TString> Systematics (){

    vector<TString> systematics {"BTAG_", "BTAG_LJET_", "BTAG_BEFF_", "BTAG_CEFF_", "BTAG_LEFF_",
                                 "BTAG_PT_", "BTAG_ETA_", "BTAG_LJET_PT_", "BTAG_LJET_ETA_", 
                                 "MASS_", "SCALE_", "MATCH_",
                                 "JES_", "JER_", "PU_", "TRIG_", "DY_"
                                };

    return systematics;
}


vector<TString> Files(TString channel = "", TString variable = ""){

    vector<TString> WhichVariable;
    vector<TString> FileVector;

    if ( variable != "" ){WhichVariable.push_back(variable);}
    else{WhichVariable = Variables();}

    for (int j=0; j<(int)WhichVariable.size(); j++){
        FileVector.push_back(TString("Plots/").Append(channel).Append("/").Append(WhichVariable.at(j)).Append("SystematicsLaTeX.txt"));
    }

    return FileVector;
}


vector<string> SplitLine(string Line ){

    //Returns a vector which its elements will be the content of 'Line' separated by blank space ' '

    vector<string> output_vector;
    istringstream iss(Line);
    copy(istream_iterator<string>(iss),
             istream_iterator<string>(),
             back_inserter<vector<string> >(output_vector));
    return output_vector;
}




double ReadLineFromFile (TString Filename, TString Systematic){

    //Returns typical error for systematic 'Systematic' in file 'Filename'

    if (Filename == "" || Systematic == ""){
        cout<<"\n\n******** ERROR ******** ERROR ******** ERROR ******** ERROR ********"<<endl;
        cout<<"You didn't provide a file path neither a systematic. Exiting!"<<endl;
        cout<<"\n\n******** ERROR ******** ERROR ******** ERROR ******** ERROR ********"<<endl;
        exit(9);
    }

    string LineToSplit;
    ifstream infile;
    vector<string> SplittedLine;

    infile.open(Filename, ios_base::in);
    if (infile.fail()) {
        cout<<"\n\n******** WARNING ******** WARNING ******** WARNING ******** WARNING ********"<<endl;
        cout<<"The file "<<Filename<<" you requested doesn't exist."<<endl;
        cout<<"Tis file will be skiped in the calculation of 'typical error'"<<endl;
        cout<<"******** WARNING ******** WARNING ******** WARNING ******** WARNING ********\n\n"<<endl;
        return -1;
    }
    while ( !infile.eof() ) {
        LineToSplit.clear();
        getline(infile, LineToSplit);
        SplittedLine.clear();
        SplittedLine = SplitLine(LineToSplit);
        if(SplittedLine.size()<=0){return -1;}
        if (SplittedLine.at(0) == Systematic){
            for(int i=0; i< (int) SplittedLine.size(); i++){
                if(SplittedLine.at(i) == "Lin.Avg.(%)="){
                    double return_value = atof(SplittedLine.at(i+1).c_str());
                    return return_value;
                }
            }
        }
    }
    infile.close();
    return -1.;
}




void TypicalError( TString channel = "", TString systematic = "", TString variable = ""){

    vector<TString> Channel, Systematic;
    
    if ( channel != ""){Channel.push_back(channel);}
    else { Channel = Channels(); }
    
    if ( systematic != "") {Systematic.push_back(systematic);}
    else { Systematic = Systematics();}
    
    for (int l=0; l<(int)Channel.size(); l++){
        vector<TString> FileList;
        FileList = Files(Channel.at(l), variable);

        vector<double> error;
        double total_error =0.0;

        for (int j=0; j<(int)Systematic.size(); j++){
            for (int i=0; i<(int)FileList.size(); i++){
                double Typ_error = ReadLineFromFile(FileList.at(i), Systematic.at(j));
                if ( Typ_error >= 0. ) { error.push_back(Typ_error);}
                cout<<"In file "<<FileList.at(i)<<" typical error "<<Typ_error<<endl;
            }

            for (int k=0; k<(int)error.size(); k++){
                total_error = total_error + error.at(k);
            }
            total_error = total_error/error.size();
            cout<<"\n\nTotal typical error for systematic "<<Systematic.at(j)<<" in channel "<<Channel.at(l)<<" is: "<<total_error<<endl;
        }
    }
}


int main(int argc, const char * const argv[]) {
    //COMPILE USING
    //g++ -o TypicalError GetTypicalErrors.C `root-config --cflags` -Wall -Wextra -pedantic -std=c++0x `root-config --ldflags --libs`
    
    //Use it via: ./TypicalError channel systematic variable
    // this will return the typical error for systematic 'systematic' in the chanel 'channel'
    // in case you want the value in one specific variable please add a 3rd option: variable
    
    
    TString channel = argc > 1 ? argv[1] : "";
    TString systematic = argc > 2 ? argv[2] : "";
    TString variable = argc > 3 ? argv[3] : "";
    
    vector<TString> ValidChannels = Channels();
    vector<TString> ValidSystematics = Systematics();
    
    if (find(ValidChannels.begin(), ValidChannels.end(), channel) != ValidChannels.end() || channel == "") {}
    else{
        cout<<"\n\nThe proposed channel '"<<channel<<"' is not valid. Exiting!\n"<<endl;
        exit(2);
    }
    
    if (find(ValidSystematics.begin(), ValidSystematics.end(), systematic) != ValidSystematics.end() || systematic == "") {}
    else{
        cout<<"\n\nThe proposed systematic '"<<systematic<<"' is not valid (or is not implemented yet). Exiting!\n"<<endl;
        exit(22);
    }
    
    TypicalError(channel, systematic, variable);

    return 0;
}


