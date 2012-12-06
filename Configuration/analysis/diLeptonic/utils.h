#ifndef util_h
#define util_h

#include "classes.h"
#include <map>
#include <TString.h>
#include <TFile.h>
#include <TStyle.h>

//convert LorentzVector to an array[E, px, py, pz]
void LVtod4(const LV lv, double *d);

// convert double to string (smart number of digits)
std::string d2s(double d);

void setHHStyle(TStyle& HHStyle);
void DrawDecayChLabel(TString decaychannel="", double textSize=0.04);
void DrawCMSLabels(int cmsprelim=true, double luminosity=0.0, double energy=8, double textSize=0.04);

class RootFileReader {
    std::map<TString, TFile*> fileMap;
    std::vector<TString> fileOrder;
    std::map<TString, int> accessed, opened;
    
    TObject* GetObj(const char* filename, const char* histoname, bool allowNonexisting);
    RootFileReader() {};
    ~RootFileReader();
public:
    //returns the singleton instnce
    static RootFileReader* getInstance();
    
    //get a histogram from the file
    template <typename T>
    void Get(const char* filename, const char* histoname, T*& result, bool allowNonexisting = false);
    
    //get a histogram from the file, you need to pass the type here
    template <typename T>
    T Get(const char* filename, const char* histoname, bool allowNonexisting = false);
};

#endif