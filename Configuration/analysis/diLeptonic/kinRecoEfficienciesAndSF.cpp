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

void drawDoubleRatio(const TH1* histNumerator1, const TH1* histDenominator1, 
               const TH1* histNumerator2, const TH1* histDenominator2,
               const Double_t ratioMin, const Double_t ratioMax, 
               const TStyle& myStyle = *gStyle)
{
    // this function draws a pad with the ratio of 'histNumerator' and 'histDenominator'
    // the range of the ratio is 'ratioMin' to 'ratioMax'
    // to the systematic variation "sys" of the enumerator "systematicVariation"
    // per default only the gaussian error of the 'histNumerator' is considered:
    // (error(bin i) = sqrt(histNumerator->GetBinContent(i))/histDenominator->GetBinContent(i))
    // if 'err_' is present and its size equals the number of bins in the histos,
    // its valus are considered as error for the ratio
    // NOTE: x Axis is transferred from histDenominator to the bottom of the canvas
    // modified quantities: none
    // used functions: none
    // used enumerators: none

    // check that histos have the same binning
    if(histNumerator1->GetNbinsX()!=histDenominator1->GetNbinsX() ||
       histNumerator1->GetNbinsX()!=histDenominator2->GetNbinsX()){
        std::cout << "error when calling drawRatio - histos have different number of bins" << std::endl;
        return;
    }

    // create ratio
    TH1* ratio1 = (TH1*)histNumerator1->Clone();
    TH1* ratio2 = (TH1*)histNumerator2->Clone();
    
    ratio1->Divide(histDenominator1);
    ratio1->SetLineStyle(histDenominator1->GetLineStyle());
    ratio1->SetLineColor(histDenominator1->GetLineColor());
    ratio1->SetMarkerStyle(histDenominator1->GetMarkerStyle());
    ratio1->SetMarkerColor(histDenominator1->GetMarkerColor());
    
    ratio2->Divide(histDenominator2);
    ratio2->SetLineStyle(histDenominator2->GetLineStyle());
    ratio2->SetLineColor(histDenominator2->GetLineColor());
    ratio2->SetMarkerStyle(histDenominator2->GetMarkerStyle());
    ratio2->SetMarkerColor(histDenominator2->GetMarkerColor());
    
    // calculate error for ratio
    // b) default: only gaussian error of histNumerator
    for(int bin=1; bin<=histNumerator1->GetNbinsX(); bin++){
        ratio1->SetBinError(bin, sqrt(histNumerator1->GetBinContent(bin))/histDenominator1->GetBinContent(bin));
        ratio2->SetBinError(bin, sqrt(histNumerator2->GetBinContent(bin))/histDenominator2->GetBinContent(bin));
    }

    Int_t    logx  = myStyle.GetOptLogx();
    Double_t left  = myStyle.GetPadLeftMargin();
    Double_t right = myStyle.GetPadRightMargin();

    // y:x size ratio for canvas
    double canvAsym = 4./3.;
    // ratio size of pad with plot and pad with ratio
    double ratioSize = 0.36;
    // change old pad
    gPad->SetBottomMargin(ratioSize);
    gPad->SetRightMargin(right);
    gPad->SetLeftMargin(left);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(0);
    gPad->SetFillColor(10);
    // create new pad for ratio plot
    TPad *rPad;
    rPad = new TPad("rPad","",0,0,1,ratioSize+0.001);
    rPad->SetFillStyle(0);
    rPad->SetFillColor(0);
    rPad->SetBorderSize(0);
    rPad->SetBorderMode(0);
    rPad->Draw();
    rPad->cd();
    rPad->SetLogy(0);
    rPad->SetLogx(logx);
    rPad->SetTicky(1);
    // configure ratio plot
    double scaleFactor = 1./(canvAsym*ratioSize);
    ratio1->SetStats(kFALSE);
    ratio1->SetTitle("");
    ratio1->SetMaximum(ratioMax);
    ratio1->SetMinimum(ratioMin);
    ratio1->SetLineWidth(1);
    // configure axis of ratio1 plot
    ratio1->GetXaxis()->SetTitleSize(histNumerator1->GetXaxis()->GetTitleSize()*scaleFactor*1.3);
    ratio1->GetXaxis()->SetTitleOffset(histNumerator1->GetXaxis()->GetTitleOffset()*0.9);
    ratio1->GetXaxis()->SetLabelSize(histNumerator1->GetXaxis()->GetLabelSize()*scaleFactor*1.4);
    ratio1->GetXaxis()->SetTitle(histNumerator1->GetXaxis()->GetTitle());
    ratio1->GetXaxis()->SetNdivisions(histNumerator1->GetNdivisions());
    ratio1->GetYaxis()->CenterTitle();
    ratio1->GetYaxis()->SetTitle("Ratio free/fixed");
    ratio1->GetYaxis()->SetTitleSize(histNumerator1->GetYaxis()->GetTitleSize()*scaleFactor);
    ratio1->GetYaxis()->SetTitleOffset(histNumerator1->GetYaxis()->GetTitleOffset()/scaleFactor*0.7);
    ratio1->GetYaxis()->SetLabelSize(histNumerator1->GetYaxis()->GetLabelSize()*scaleFactor);
    ratio1->GetYaxis()->SetLabelOffset(histNumerator1->GetYaxis()->GetLabelOffset()*3.3);
    ratio1->GetYaxis()->SetTickLength(0.03);
    ratio1->GetYaxis()->SetNdivisions(505);
    ratio1->GetXaxis()->SetRange(histNumerator1->GetXaxis()->GetFirst(), histNumerator1->GetXaxis()->GetLast());

    // delete axis of initial plot
    histNumerator1->GetXaxis()->SetLabelSize(0);
    histNumerator1->GetXaxis()->SetTitleSize(0);
    
    // draw ratio plot
    ratio1->SetLineWidth(0);
    ratio1->SetLineColor(kBlack);
    ratio2->SetLineColor(kBlack);
    ratio1->SetMarkerColor(kBlack);
    //ratio1-
     ratio1->DrawClone("p ex0");
     ratio2->DrawClone("same");
    
    
    
    rPad->SetTopMargin(0.0);
    rPad->SetBottomMargin(0.15*scaleFactor);
    rPad->SetRightMargin(right);
    gPad->SetLeftMargin(left);
    gPad->RedrawAxis();
    
    rPad->SetGrid(1,1);
}



TH1 *getRatio(TH1 *bkr, TH1 *akr) {
    bkr->Sumw2();
    akr->Sumw2();
    akr->Divide(akr, bkr, 1, 1, "B");
    return akr;    
}

void saveRootAndEps(TString plotsDir, TH1 *h, TString name) {
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();
    h->GetXaxis()->SetDefaults();
    h->GetXaxis()->SetRangeUser(xmin, xmax);
    h->Write(name);
    
    TCanvas c;
    h->Draw();
    c.SaveAs(plotsDir + "/kinReco/" + name + ".eps");    
}

inline double sqr(double value) { return value*value; }


void makeEffiPlots(TString plotsDir) {
    gSystem->mkdir(plotsDir + "/kinReco", true);
    TFile out(plotsDir + "/kinRecoPlots.root", "RECREATE");
    auto reader = RootFileReader::getInstance();
    
    for (TString ch : {"ee", "emu", "mumu", "combined"}) 
    for (TString plot : {"LeptonpT", "LeptonEta", "MET", "BjetEta"})
    {
        TString source(plotsDir + "/");
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
        saveRootAndEps(plotsDir, dataRatio, ch + "_" + plot + "_dataRatio");
        saveRootAndEps(plotsDir, mcsignalRatio, ch + "_" + plot + "_mcSignalRatio");
        saveRootAndEps(plotsDir, mcallRatio, ch + "_" + plot + "_mcAllRatio");
        saveRootAndEps(plotsDir, sf, ch + "_" + plot + "_sf_signal");
        saveRootAndEps(plotsDir, sfall, ch + "_" + plot + "_sf_allmc");
        
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
        c.SaveAs(plotsDir + "/kinReco/" + ch + "_" + plot + "_3in1.eps");
    }
    out.Write();
    std::cout << "Just a reminder: if you include the Kin Reco SF in the Analysis.C, the results\n"
        << "shown here should show a scale factor of 1.0 (test this, its a cross check)! To\n"
        << "get the scale factor, set weightKinFit=1 in the Analysis.C and rerun.\n";
}

void createRelativeEfficiencies(TString dir1, TString dir2, TString quantity) {
    auto reader = RootFileReader::getInstance();
    
    TString filename1 = dir1 + "/combined/preunfolded_" + quantity + "_source.root";
    TString filename2 = dir2 + "/combined/preunfolded_" + quantity + "_source.root";
    
    TH1 *mcFree = reader->GetClone<TH1>(filename1, quantity + "_allmc");
    TH1 *dataFree = reader->GetClone<TH1>(filename1, quantity + "_data");
    TH1 *mcLimited = reader->GetClone<TH1>(filename2, quantity + "_allmc");
    TH1 *dataLimited = reader->GetClone<TH1>(filename2, quantity + "_data");
    
    gStyle->SetOptStat(0);
    TCanvas c;
    mcFree->SetLineColor(kRed); 
    mcFree->SetFillColor(kWhite);
    mcFree->SetMarkerSize(0);
    mcLimited->SetLineColor(kBlue);
    mcLimited->SetFillColor(kWhite);
    mcLimited->SetMarkerSize(0);
    
    dataFree->SetMarkerStyle(5);
    dataFree->SetMarkerColor(kRed);
    dataFree->SetMarkerSize(1.5);            

    dataLimited->SetMarkerStyle(5);
    dataLimited->SetMarkerColor(kBlue);
    dataLimited->SetMarkerSize(1.5);            
    
    mcFree->Draw();
    mcLimited->Draw("same");
    dataFree->Draw("P,same");
    dataLimited->Draw("P,same");
    
    TLegend leg(.64,.65,.94,.85);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry(mcFree, "MC, free mass");
    leg.AddEntry(mcLimited, "MC, fixed mass");
    leg.AddEntry(dataFree, "Data, free mass", "p");
    leg.AddEntry(dataLimited, "Data, fixed mass", "p");  
    leg.Draw("same");

    drawDoubleRatio(dataLimited, dataFree, mcLimited, mcFree, 0.6, 1.2, *gStyle);
    
    
    c.SaveAs("kinRecoRelativeEffi"+quantity+".eps");
    
}


int main() {
    makeEffiPlots("Plots_nokinSF_100300");
//     makeEffiPlots("Plots_nokinSF_163183");
    makeEffiPlots("Plots_nokinSF_173");
    
    createRelativeEfficiencies("Plots_nokinSF_100300", "Plots_nokinSF_173", "HypToppT");
    createRelativeEfficiencies("Plots_nokinSF_100300", "Plots_nokinSF_173", "HypTopRapidity");
    
    return 0;
}