#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TSystem.h>
#include <Rtypes.h>
#include <TAxis.h>
#include <TLegend.h>
#include "utils.h"

TH1 *getRatio(TH1 *bkr, TH1 *akr) {
    bkr->Sumw2();
    akr->Sumw2();
    akr->Divide(akr, bkr, 1, 1, "B");
    return akr;    
}

void saveRootAndEps(TH1 *h, TString name) {
    h->Write(name);
    TCanvas c;
    h->Draw();
    c.SaveAs("Plots/kinReco/" + name + ".eps");    
}

int main() {
    gSystem->mkdir("Plots/kinReco", true);
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
        saveRootAndEps(dataRatio, ch + "_" + plot + "_dataRatio");
        saveRootAndEps(mcsignalRatio, ch + "_" + plot + "_mcSignalRatio");
        saveRootAndEps(mcallRatio, ch + "_" + plot + "_mcAllRatio");
        saveRootAndEps(sf, ch + "_" + plot + "_sf_signal");
        saveRootAndEps(sfall, ch + "_" + plot + "_sf_allmc");
        
        TCanvas c;
        dataRatio->Draw();
        dataRatio->SetTitle(ch + " channel: " + plot + " - Kin. Reco. behaviour");
        dataRatio->GetYaxis()->SetRangeUser(0.6, 1.1);
        dataRatio->GetYaxis()->SetTitleOffset(1.1);
        mcallRatio->SetMarkerColor(kRed);
        mcallRatio->SetLineColor(kRed);
        mcallRatio->Draw("same");
        sfall->SetMarkerColor(kBlue);
        sfall->Draw("same");
        //sfall->Fit("pol0", "", "same");
        TLegend l(0.78, 0.95, 0.99, 0.7);
        l.AddEntry(dataRatio, "eff. in data");
        l.AddEntry(mcallRatio, "eff. in MC");
        l.AddEntry(sfall, "SF");
        l.Draw("same");
        c.SaveAs("Plots/kinReco/" + ch + "_" + plot + "_3in1.eps");
    }
    out.Write();
}
