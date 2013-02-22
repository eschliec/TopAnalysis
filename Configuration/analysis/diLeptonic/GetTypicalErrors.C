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

#include "CommandLineParameters.hh"

using namespace std;


vector<TString> Channels(){
    
  vector<TString> channel {"ee", "emu", "mumu", "combined"};
    
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

  vector<TString> systematics {"BTAG_", "BTAG_LJET_", "BTAG_BEFF_",
                                 "BTAG_PT_", "BTAG_ETA_", "BTAG_LJET_PT_", "BTAG_LJET_ETA_", 
      "MASS_", "SCALE_", "MATCH_", "HAD_", "KIN_", "LEPT_",
      "JES_", "JER_", "PU_", "TRIG_", "DY_", "BG_", "PDF_"
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



void SanityCheck( TString channel = "", TString systematic = "", TString variable = ""){

    vector<TString> ValidChannel = Channels(), ValidSystematic = Systematics(), ValidVariable = Variables();
    
    if (find(ValidChannel.begin(), ValidChannel.end(), channel) != ValidChannel.end() || channel == "") {}
    else{
        cout<<"\n\nThe proposed channel '"<<channel<<"' is not valid. Exiting!\n"<<endl;
        exit(2);
    }
    
    if (find(ValidSystematic.begin(), ValidSystematic.end(), systematic) != ValidSystematic.end() || systematic == "") {}
    else{
        cout<<"\n\nThe proposed systematic '"<<systematic<<"' is not valid (or is not implemented yet). Exiting!\n"<<endl;
        exit(22);
    }
    
    if (find(ValidVariable.begin(), ValidVariable.end(), variable) != ValidVariable.end() || variable == "") {}
    else{
        cout<<"\n\nThe proposed variable '"<<variable<<"' is not valid (or is not implemented yet). Exiting!\n"<<endl;
        exit(222);
    }

}

void TypicalError( TString channel = "", TString systematic = "", TString variable = ""){

    SanityCheck(channel, systematic, variable);
    
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


int main(int argc, char** argv) {
    CLParameter<std::string> opt_v("v", "Return the typical error for certain variable, e.g. 'ToppTLead', 'LLBarMass', ...", false, 1, 1);
    CLParameter<std::string> opt_s("s", "Return the typical systematic uncertainty for a certain systematic variation, e.g. 'PU_', 'TRIG_', 'BTAG_LJET_ETA_', 'BTAG_PT_', ...", false, 1, 1);
    CLParameter<std::string> opt_c("c", "Return the typical systematic uncertainty for an specific channel (ee, emu, mumu). No channel specified = run on all channels", false, 1, 1,
            [](const std::string &ch){return ch == "" || ch == "ee" || ch == "emu" || ch == "mumu" || ch == "combined";});
    CLAnalyser::interpretGlobal(argc, argv);
    
    TString ValidSystematics = opt_s.isSet() ? opt_s[0] : "";
    TString ValidVariable    = opt_v.isSet() ? opt_v[0] : "";
    TString ValidChannel     = opt_c.isSet() ? opt_c[0] : "";
        
    TypicalError(ValidChannel, ValidSystematics, ValidVariable);
    
    return 0;
}
