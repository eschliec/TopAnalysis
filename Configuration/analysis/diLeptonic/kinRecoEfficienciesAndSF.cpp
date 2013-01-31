#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TSystem.h>
#include <Rtypes.h>
#include <TAxis.h>
#include <TLegend.h>
#include <string>
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
        source.Append(ch).Append("/Nominal/");
        
        TH1 *dataRatio = getRatio(reader->GetClone<TH1>(source + plot + "bkr_source.root", plot + "bkr" + "_data"),
                                  reader->GetClone<TH1>(source + plot + "akr_source.root", plot + "akr" + "_data"));

        TH1 *mcsignalRatio = getRatio(reader->GetClone<TH1>(source + plot + "bkr_source.root", plot + "bkr" + "_signalmc"),
                                      reader->GetClone<TH1>(source + plot + "akr_source.root", plot + "akr" + "_signalmc"));
        
        TH1 *mcallRatio = getRatio(reader->GetClone<TH1>(source + plot + "bkr_source.root", plot + "bkr" + "_allmc"),
                                   reader->GetClone<TH1>(source + plot + "akr_source.root", plot + "akr" + "_allmc"));
        
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
        dataRatio->GetYaxis()->SetRangeUser(0.5, 1.2);
        dataRatio->GetYaxis()->SetTitleOffset(1.1);
        mcallRatio->SetMarkerColor(kRed);
        mcallRatio->SetLineColor(kRed);
        mcallRatio->Draw("same");
        sfall->SetMarkerColor(kBlue);
        sfall->Draw("same");
        //sfall->Fit("pol0", "", "same");
        
        double dataEff = reader->Get<TH1>(source + "step9_source.root", "step9_data")->GetBinContent(2) / 
                         reader->Get<TH1>(source + "step8_source.root", "step8_data")->GetBinContent(2);
        double allmcEff = reader->Get<TH1>(source + "step9_source.root", "step9_allmc")->GetBinContent(2) / 
                          reader->Get<TH1>(source + "step8_source.root", "step8_allmc")->GetBinContent(2);
        char dataEffString[100]; sprintf(dataEffString, "%.2f%%", 100*dataEff);
        char allmcEffString[100]; sprintf(allmcEffString, "%.2f%%", 100*allmcEff);
        char sfString[100]; sprintf(sfString, "%.2f%%", 100*dataEff/allmcEff);
        
        TLegend l(0.73, 0.95, 0.99, 0.7);
        l.AddEntry(dataRatio, TString("eff data: ") + dataEffString);
        l.AddEntry(mcallRatio, TString("eff MC: ") + allmcEffString);
        l.AddEntry(sfall, TString("SF: ") + sfString);
        l.Draw("same");
        c.SaveAs("Plots/kinReco/" + ch + "_" + plot + "_3in1.eps");
    }
    out.Write();
}
