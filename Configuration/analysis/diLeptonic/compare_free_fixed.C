#include <TString.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <iomanip>

//run these commands to compare two results:
//rm compare_free_fixed.txt
//for i in `awk '{print $1}'<HistoList`;do root -b -q compare_free_fixed.C'("'$i'")' >>| compare_free_fixed.txt;done

void compare_free_fixed(TString quantity) {

//     TFile ffixed("Plots_kinSF_173/combined/DiffXS_"+quantity+"_source.root");
//     TFile ffree("Plots_kinSF_100300/combined/DiffXS_"+quantity+"_source.root");
    TFile ffixed("Plots_fakedata_173/combined/DiffXS_"+quantity+"_source.root");
    TFile ffree("Plots_fakedata_100300/combined/DiffXS_"+quantity+"_source.root");

    TCanvas canvas;
    TLegend l(0.67, 0.67, 0.93, 0.86);
    
    ffixed.cd();
    mc->SetLineColor(kBlue);
    mc->Draw();
    data_staterror_only->SetMarkerColor(kBlue);
    data_staterror_only->SetLineWidth(2);
    data_staterror_only->SetLineColor(kBlue);
    data_staterror_only->SetFillColor(kWhite);
    TGraphAsymmErrors *fixed_data = data_staterror_only->Clone();
    
    data_staterror_only->Draw("p,same");
    l.AddEntry(data_staterror_only, "m_t = 173 GeV");

    ffree.cd();
    mc->Draw("same");
    double x,y,diff;
    data_staterror_only->SetLineWidth(2);
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
    l.AddEntry(data_staterror_only, "m_t = 100 .. 300 GeV");
    l.AddEntry(mc, "MadGraph");
    
    l.Draw();
    canvas.SaveAs("cmp_kinreco_" + quantity + ".eps");
    
}
