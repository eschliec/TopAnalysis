#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TSystem.h>
#include <Rtypes.h>
#include <TAxis.h>
#include <TLegend.h>
#include <string>
#include <iostream>
#include "utils.h"

TH1 *getRatio(TH1 *bkr, TH1 *akr) {
    bkr->Sumw2();
    akr->Sumw2();
    akr->Divide(akr, bkr, 1, 1, "B");
    return akr;    
}

void saveRootAndEps(TH1 *h, TString name) {
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();
    h->GetXaxis()->SetDefaults();
    h->GetXaxis()->SetRangeUser(xmin, xmax);
    h->Write(name);
    
    TCanvas c;
    h->Draw();
    c.SaveAs("Plots/kinReco/" + name + ".eps");    
}

inline double sqr(double value) { return value*value; }

int main() {
    gSystem->mkdir("Plots/kinReco", true);
    TFile out("Plots/kinRecoPlots.root", "RECREATE");
    auto reader = RootFileReader::getInstance();
    
    for (TString ch : {"ee", "emu", "mumu", "combined"}) 
    for (TString plot : {"LeptonpT", "LeptonEta", "MET", "BjetEta"})
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
        
        double NdataAfter = reader->Get<TH1>(source + "step9_source.root", "step9_data")->GetBinContent(2);
        double NdataBefore = reader->Get<TH1>(source + "step8_source.root", "step8_data")->GetBinContent(2);
        double dataEff = NdataAfter / NdataBefore;
        double dataEffUnc = sqrt(dataEff * (1-dataEff) / NdataBefore);
        
        double NmcAfter = reader->Get<TH1>(source + "step9_source.root", "step9_allmc")->GetBinContent(2);
        double NmcBefore = reader->Get<TH1>(source + "step8_source.root", "step8_allmc")->GetBinContent(2); 
        double allmcEff =  NmcAfter / NmcBefore;
        double allmcEffUnc = sqrt(allmcEff * (1-allmcEff) / NmcBefore);
        
        double sfUnc = sqrt(sqr(dataEffUnc/dataEff) + sqr(allmcEffUnc/allmcEff));
                          
        char dataEffString[100]; sprintf(dataEffString, "%.2f%%", 100*dataEff);
        char allmcEffString[100]; sprintf(allmcEffString, "%.2f%%", 100*allmcEff);
        char sfString[100]; sprintf(sfString, "%.2f%% +- %.2f", 100*dataEff/allmcEff, 100*sfUnc);
        
        TLegend l(0.73, 0.95, 0.99, 0.7);
        l.AddEntry(dataRatio, TString("eff data: ") + dataEffString);
        l.AddEntry(mcallRatio, TString("eff MC: ") + allmcEffString);
        l.AddEntry(sfall, TString("SF: ") + sfString);
        l.Draw("same");
        c.SaveAs("Plots/kinReco/" + ch + "_" + plot + "_3in1.eps");
    }
    out.Write();
    std::cout << "Just a reminder: if you include the Kin Reco SF in the Analysis.C, the results\n"
        << "shown here should show a scale factor of 1.0 (test this, its a cross check)! To\n"
        << "get the scale factor, set weightKinFit=1 in the Analysis.C and rerun.\n";
}
