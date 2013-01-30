#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h>
#include "utils.h"

TH1 *getRatio(TH1 *bkr, TH1 *akr) {
    bkr->Sumw2();
    akr->Sumw2();
    akr->Divide(akr, bkr, 1, 1, "B");
    return akr;    
}


int main() {
    TFile out("Plots/kinRecoPlots.root", "RECREATE");
    auto reader = RootFileReader::getInstance();
    
    for (TString ch : {"ee", "emu", "mumu", "combined"}) 
    for (TString plot : {"LeptonpT", "LeptonEta"})
    {
        TString source("Plots/");
        source.Append(ch).Append("/Nominal/").Append(plot);
        
        TH1 *dataRatio = getRatio(reader->GetClone<TH1>(source + "bkr_data.root", plot + "bkr"),
                                  reader->GetClone<TH1>(source + "akr_data.root", plot + "akr"));

        TH1 *mcsignalRatio = getRatio(reader->GetClone<TH1>(source + "bkr_signalmc.root", plot + "bkr"),
                                      reader->GetClone<TH1>(source + "akr_signalmc.root", plot + "akr"));
        
        TH1 *mcallRatio = getRatio(reader->GetClone<TH1>(source + "bkr_summc.root", plot + "bkr"),
                                   reader->GetClone<TH1>(source + "akr_summc.root", plot + "akr"));
        
        auto sf = static_cast<TH1*>(dataRatio->Clone());
        sf->Divide(mcsignalRatio);

        auto sfall = static_cast<TH1*>(dataRatio->Clone());
        sfall->Divide(mcallRatio);
        
        out.cd();
        dataRatio->Write(ch + "_" + plot + "_dataRatio");
        mcsignalRatio->Write(ch + "_" + plot + "_mcSignalRatio");
        mcallRatio->Write(ch + "_" + plot + "_mcAllRatio");
        sf->Write(ch + "_" + plot + "_sf_signal");
        sfall->Write(ch + "_" + plot + "_sf_allmc");
    }
    out.Write();
}
