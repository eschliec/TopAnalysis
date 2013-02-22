#include <TString.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include <TFrame.h>
#include <iomanip>

//run these commands to compare two results:
//rm compare_free_fixed.txt
//for i in `awk '{print $1}'<HistoList`;do root -b -q compare_free_fixed.C'("'$i'")' >>| compare_free_fixed.txt;done
//or for a single plot:
// root -b -q compare_free_fixed.C'("'HypToppT'")'

void compare_free_fixed(TString quantity) {

//     TFile ffixed("Plots_kinSF_173/combined/DiffXS_"+quantity+"_source.root");
//     TFile ffree("Plots_kinSF_100300/combined/DiffXS_"+quantity+"_source.root");
    TFile fvariation("Plots_fakedata_173/combined/DiffXS_"+quantity+"_source.root");
    TFile fnominal("Plots_fakedata_100300/combined/DiffXS_"+quantity+"_source.root");
//     TFile ffixed("Plots/PDF_1_UP/combined/DiffXS_"+quantity+"_source.root");
//     TFile ffree("Plots/combined/DiffXS_"+quantity+"_source.root");

    TCanvas canvas;
    canvas.SetMargin(0.155, 0.040, 0.13, 0.1);
    TLegend l(0.73, 0.77, 0.99, 0.96);
    
    gStyle->SetOptStat(0);
    
    fvariation.cd();
    int MaxBin = mc->GetMaximumBin();
    int MinBin = mc->GetMinimumBin();
    double Max = mc->GetBinContent(MaxBin);
    double Min = mc->GetBinContent(MinBin);
    mc->GetYaxis()->SetRangeUser(0.1*Min, 1.2*Max);
    mc->SetLineColor(kBlue);
    mc->SetMarkerSize(0);
    mc->Draw();
    data_staterror_only->SetMarkerColor(kBlue);
    data_staterror_only->SetMarkerStyle(32);
    data_staterror_only->SetLineColor(kWhite);
    data_staterror_only->SetFillColor(kWhite);
    TGraphAsymmErrors *fixed_data = data_staterror_only->Clone();
    
    data_staterror_only->Draw("p,same");
    l.AddEntry(data_staterror_only, "Pseudo data (fixed m_t)");
    l.AddEntry(mc, "MadGraph (fixed m_t)");

    fnominal.cd();
    mc->SetLineStyle(2);
    mc->SetMarkerSize(0);
    mc->Draw("same");
    double x,y,diff;
//     data_staterror_only->SetLineWidth(2);
    data_staterror_only->SetMarkerStyle(26);
    data_staterror_only->SetMarkerColor(kRed);
    data_staterror_only->SetLineColor(kWhite);
    data_staterror_only->SetFillColor(kWhite);
    data_staterror_only->GetPoint(1,diff,y);
    data_staterror_only->GetPoint(0,x,y);
    diff -= x;
    
    for (int i = 0; i < data_staterror_only->GetN(); ++i) {
        double fixed_y;
        fixed_data->GetPoint(i,x,fixed_y);
        data_staterror_only->GetPoint(i,x,y);
        data_staterror_only->SetPoint(i,x+diff/6,y);
        cout << quantity << " " << i << " " << std::setprecision(2) << 100*(1-y/fixed_y) << " %\n";
    }
    
    
    data_staterror_only->Draw("p,same");
    l.AddEntry(data_staterror_only, "Pseudo data (default)");
    l.AddEntry(mc, "MadGraph (default)");
    
    l.Draw();
    canvas.SaveAs("cmp_kinreco_" + quantity + ".eps");
    
}
