#define SYS_ERROR_MACRO

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <sys/stat.h>
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPostScript.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include "TopAnalysis/Configuration/analysis/semiLeptonic/diffXSection/basicFunctions.h"
#include "TopAnalysis/Configuration/analysis/semiLeptonic/diffXSection/HHStyle.h"

using namespace std;

/// CONFIGURABLES

// path to the ingoing root histogram files
//const TString inpath("~markusm/cms/systematics/");
//const TString inpath("/afs/naf.desy.de/user/w/wbehrenh/cms/systematicsOct14/");
const TString inpath("/afs/naf.desy.de/user/w/wbehrenh/cms/systematicsNov4/");

// path where the output is written
//const TString outpath("Markus/DiffXS2011/Systematics/plots/");
const TString outpath("~wbehrenh/cms/systematicsNov4-output/");

// output format
// const TString outform(".png");
const TString outform(".eps");


///////////////////////////////////////////////////////////

Bool_t set_LogY = kTRUE;
//Bool_t set_LogY = kFALSE;

Bool_t isDiff   = kTRUE;
// Bool_t isDiff   = kFALSE;

// output file name
const TString outfileName(isDiff ? outpath+"Systematic_Errors_DIFF.root" : outpath+"Systematic_Errors_TOTAL.root");

Bool_t isTotal = !(isDiff);

//Bool_t isPSE = kTRUE;    // Make efficiency plots "_PSE"
Bool_t isPSE = kFALSE;

Bool_t isModel = kTRUE;
//Bool_t isModel = kFALSE;

Bool_t isExp = kTRUE;
//Bool_t isExp = kFALSE;

// Bool_t isStat  = kTRUE;
Bool_t isStat  = kFALSE;

///////////////////////////////////////////////////////////

TObjArray Hlist;

TFile    *files[100];
THStack  *sum_stat_errors[1000];
THStack  *sum_sys_errors[1000];
THStack  *sum_exp_errors[1000];
THStack  *sum_mod_errors[1000];

//TLegend* leg = new TLegend(0.70,0.58,0.95,0.87);
TLegend* leg = new TLegend(0.25,0.58,0.50,0.87);


void setupTotalErrors() {

  TList *ListOfKeys = files[0]->GetListOfKeys();

  for(Int_t i = 0; i < ListOfKeys->GetEntries(); ++i) {

    sum_stat_errors[i] = new THStack("","");
    sum_sys_errors[i]  = new THStack("","");
    sum_exp_errors[i]  = new THStack("","");
    sum_mod_errors[i]  = new THStack("","");

  }

}

void setHistogramStyle(TH1* hist, Color_t LineColor, Style_t LineStyle, Width_t LineWidth, Color_t HistColor) {

  hist->SetLineColor(LineColor);
  hist->SetLineStyle(LineStyle);
  hist->SetLineWidth(LineWidth);
  hist->SetFillColor(HistColor);

  return;

}

TH1* SummedStackHisto(const THStack *stack) {

  TList* l = stack->GetHists();

  if ( l->GetEntries() == 0 )  return 0;

  TH1* result = (TH1*)l->At(0)->Clone();

  for(int i = 1; i < l->GetEntries(); ++i) {

    result->Add((TH1*)l->At(i));

  }

  return result;

}

void SymmetricAroundZero(TH1* h_ref, TH1* h_var_up, TH1* h_var_down, TH1* h_sys, TH1* h_sys2, double &Sum_Errors ) {

  Int_t N_bins = h_ref->GetNbinsX();

  Double_t Sys_Error      = 0.;
  Double_t Sys_Error2     = 0.;
  Double_t Sys_Error_Up   = 0.;
  Double_t Sys_Error_Down = 0.;

  for(Int_t bin = 1; bin <= h_ref->GetNbinsX()/2; ++bin) {

    if( h_ref->GetBinContent(bin) == 0 )  continue;

    Double_t PlusMinus_Average_Ref  = 0.;
    Double_t PlusMinus_Average_Up   = 0.;
    Double_t PlusMinus_Average_Down = 0.;

    PlusMinus_Average_Ref  = (h_ref->GetBinContent(bin)      + h_ref->GetBinContent(N_bins+1-bin)     ) / 2.;
    PlusMinus_Average_Up   = (h_var_up->GetBinContent(bin)   + h_var_up->GetBinContent(N_bins+1-bin)  ) / 2.;
    PlusMinus_Average_Down = (h_var_down->GetBinContent(bin) + h_var_down->GetBinContent(N_bins+1-bin)) / 2.;

    Sys_Error_Up   = abs(PlusMinus_Average_Ref - PlusMinus_Average_Up  ) / PlusMinus_Average_Ref;
    Sys_Error_Down = abs(PlusMinus_Average_Ref - PlusMinus_Average_Down) / PlusMinus_Average_Ref;

    Sys_Error  = (Sys_Error_Up+Sys_Error_Down)/2.;
    Sys_Error2 = Sys_Error * Sys_Error;

    h_sys->SetBinContent(           bin, Sys_Error  );
    h_sys->SetBinContent(  N_bins+1-bin, Sys_Error  );
    h_sys2->SetBinContent(          bin, Sys_Error2 );
    h_sys2->SetBinContent( N_bins+1-bin, Sys_Error2 );

    cout << "Sys. Error: " << Sys_Error << endl;
    Sum_Errors += (2.*Sys_Error);

  }

}


void SystematicErrors() {

  TH1 *dysys[1000], *dysys2[1000];
  
  leg->Clear();
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  files[0] = new TFile(inpath.Copy().Append("standard.root"));

  TList *ListOfKeys = files[0]->GetListOfKeys();

  setupTotalErrors();

  TH1 *h_ref;
  TH1 *h_style;
  TH1 *h_var_up;
  TH1 *h_var_down;

  Double_t Sys_Error      = 0.;
  Double_t Sys_Error2     = 0.;
  Double_t Sys_Error_Up   = 0.;
  Double_t Sys_Error_Down = 0.;

  Color_t LineColor = kBlack;
  Style_t LineStyle = 1;
  Width_t LineWidth = 1;
  Color_t HistColor = kRed+1;

  //  Statistical Error --------------------------------------------------

  LineColor = kBlack;
  LineStyle = 1;
  LineWidth = 1;
  HistColor = kGray+3;

  if( files[0] && !files[0]->IsZombie()  && isStat ) {

    cout << "   >>>   " << "Statistical Error" << endl;
    cout << "--------------------------------" << endl << endl;

    for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

      TString title = ListOfKeys->At(j)->GetName();
      cout << "   >>>   " << title << endl;
      if( title.Contains("PSE_") )  continue;

      Double_t Stat_Error  = 0.;
      Double_t Stat_Error2 = 0.;
      Double_t Sum_Errors  = 0.;

      files[0]->GetObject( title, h_ref);
      files[0]->GetObject( title, h_style);

      TH1* h_stat  = (TH1*) h_ref->Clone();   h_stat->Reset();
      TH1* h_stat2 = (TH1*) h_ref->Clone();   h_stat2->Reset();

      Int_t N_bins = h_ref->GetNbinsX();

      for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {

	if( h_ref->GetBinContent(bin) == 0 )  continue;

	Stat_Error  = h_ref->GetBinError(bin)/h_ref->GetBinContent(bin);
	Stat_Error2 = Stat_Error * Stat_Error;

	h_stat->SetBinContent(  bin, Stat_Error  );
	h_stat2->SetBinContent( bin, Stat_Error2 );
	cout << "Stat. Error: " << Stat_Error << endl;
	Sum_Errors = Sum_Errors+Stat_Error;
      }

      cout << "-------------" << endl;
      cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
      cout << "-------------" << endl << endl;

      setHistogramStyle(h_stat,  LineColor, LineStyle, LineWidth, HistColor);
      setHistogramStyle(h_stat2, LineColor, LineStyle, LineWidth, HistColor);
      h_style = h_stat;

      title.Append("_STAT_error");
      h_stat->SetName(title);
      Hlist.Add(h_stat);

      sum_sys_errors[j]->Add(h_stat2);
      sum_stat_errors[j]->Add(h_stat2);

    }

    leg->AddEntry(h_style,"Statistical","F");
    cout << "-------------" << endl << endl;

  }

  //  JER --------------------------------------------------

  LineColor = kBlack;
  LineStyle = 1;
  LineWidth = 1;
  HistColor = kRed+1;

  files[1] = new TFile(inpath.Copy().Append("jer_up.root"));
  files[2] = new TFile(inpath.Copy().Append("jer_down.root"));

  if( files[1] && !files[1]->IsZombie() && files[2] && !files[2]->IsZombie() && isExp ) {

    cout << "   >>>   " << "Jet Energy RESOLUTION" << endl;
    cout << "------------------------------------" << endl << endl;

    for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

      TString title = ListOfKeys->At(j)->GetName();
      cout << "   >>>   " << title << endl;
      if( title.Contains("PSE_") )  continue;

      Sys_Error  = 0.;
      Sys_Error2 = 0.;

      files[0]->GetObject( title, h_ref);
      files[0]->GetObject( title, h_style);
      files[1]->GetObject( title, h_var_up);
      files[2]->GetObject( title, h_var_down);

      TH1* h_sys       = (TH1*) h_var_up->Clone();   h_sys->Reset();
      TH1* h_sys2      = (TH1*) h_var_down->Clone(); h_sys2->Reset();
      TH1* h_sys_up    = (TH1*) h_var_up->Clone();   h_sys_up->Reset();
      TH1* h_sys_down  = (TH1*) h_var_down->Clone(); h_sys_down->Reset();
      TH1* h_sys2_up   = (TH1*) h_var_up->Clone();   h_sys2_up->Reset();
      TH1* h_sys2_down = (TH1*) h_var_down->Clone(); h_sys2_down->Reset();

      TCanvas* Canvas = new TCanvas("plot", "plot", 800, 800);

      setHistogramStyle(h_ref,      kBlack, 1, 2, HistColor);
      setHistogramStyle(h_var_up,   kBlue,  1, 1, HistColor);
      setHistogramStyle(h_var_down, kRed,   1, 1, HistColor);

      Canvas->Clear();
      h_ref->SetMarkerColor(kBlack);
      h_ref->Draw();
      h_var_up->SetMarkerColor(kBlue);
      h_var_up->Draw("SAME");
      h_var_down->SetMarkerColor(kRed);
      h_var_down->Draw("SAME");
      if( set_LogY ) gPad->SetLogy();
      Canvas->Print(outpath.Copy().Append(title).Append("_JER").Append(outform));

      Int_t N_bins = h_ref->GetNbinsX();
      Double_t Sum_Errors = 0.;

      for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {

	if( h_ref->GetBinContent(bin) == 0 )  continue;

 	Sys_Error_Up   = abs(h_ref->GetBinContent(bin) - h_var_up->GetBinContent(bin))/h_ref->GetBinContent(bin);
	Sys_Error_Down = abs(h_ref->GetBinContent(bin) - h_var_down->GetBinContent(bin))/h_ref->GetBinContent(bin);
  
	Sys_Error  = (Sys_Error_Up+Sys_Error_Down)/2.;
	Sys_Error2 = Sys_Error * Sys_Error;

	h_sys->SetBinContent(  bin, Sys_Error  );
	h_sys2->SetBinContent( bin, Sys_Error2 );
	cout << "Sys. Error: " << Sys_Error << endl;
	Sum_Errors = Sum_Errors+Sys_Error;
      }

      cout << "-------------" << endl;
      cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
      cout << "-------------" << endl << endl;

      setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
      setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
      h_style = h_sys;

      title.Append("_sys_error_JER");
      h_sys->SetName(title);
      Hlist.Add(h_sys);

      sum_sys_errors[j]->Add(h_sys2);
      sum_exp_errors[j]->Add(h_sys2);

    }

    leg->AddEntry(h_style,"JER","F");
    cout << "-------------" << endl << endl;

  }

  //  JES --------------------------------------------------

  LineColor = kBlack;
  LineStyle = 1;
  LineWidth = 1;
  HistColor = kRed-7;

  files[3] = new TFile(inpath.Copy().Append("jes_up.root"));
  files[4] = new TFile(inpath.Copy().Append("jes_down.root"));

  if( files[3] && !files[3]->IsZombie() && files[4] && !files[4]->IsZombie() && isExp ) {

    cout << "   >>>   " << "Jet Energy SCALE" << endl;
    cout << "-------------------------------" << endl << endl;

    for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

      TString title = ListOfKeys->At(j)->GetName();
      cout << "   >>>   " << title << endl;
      if( title.Contains("PSE_") )  continue;

      Sys_Error      = 0.;
      Sys_Error2     = 0.;
      Sys_Error_Up   = 0.;
      Sys_Error_Down = 0.;

      files[0]->GetObject( title, h_ref);
      files[0]->GetObject( title, h_style);
      files[3]->GetObject( title, h_var_up);
      files[4]->GetObject( title, h_var_down);

      TH1* h_sys       = (TH1*) h_var_up->Clone();   h_sys->Reset();
      TH1* h_sys2      = (TH1*) h_var_down->Clone(); h_sys2->Reset();
      TH1* h_sys_up    = (TH1*) h_var_up->Clone();   h_sys_up->Reset();
      TH1* h_sys_down  = (TH1*) h_var_down->Clone(); h_sys_down->Reset();
      TH1* h_sys2_up   = (TH1*) h_var_up->Clone();   h_sys2_up->Reset();
      TH1* h_sys2_down = (TH1*) h_var_down->Clone(); h_sys2_down->Reset();

      TCanvas* Canvas = new TCanvas("plot", "plot", 800, 800);

      setHistogramStyle(h_ref,      kBlack, 1, 2, HistColor);
      setHistogramStyle(h_var_up,   kBlue,  1, 1, HistColor);
      setHistogramStyle(h_var_down, kRed,   1, 1, HistColor);

      Canvas->Clear();
      h_ref->SetMarkerColor(kBlack);
      h_ref->Draw();
      h_var_up->SetMarkerColor(kBlue);
      h_var_up->Draw("SAME");
      h_var_down->SetMarkerColor(kRed);
      h_var_down->Draw("SAME");
      if( set_LogY ) gPad->SetLogy();
      Canvas->Print(outpath.Copy().Append(title).Append("_JES").Append(outform));

      Int_t N_bins = h_ref->GetNbinsX();
      Double_t Sum_Errors = 0.;

      for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {

	if( h_ref->GetBinContent(bin) == 0 )  continue;

	Sys_Error_Up   = abs(h_ref->GetBinContent(bin) - h_var_up->GetBinContent(bin))/h_ref->GetBinContent(bin);
	Sys_Error_Down = abs(h_ref->GetBinContent(bin) - h_var_down->GetBinContent(bin))/h_ref->GetBinContent(bin);
  
	Sys_Error  = (Sys_Error_Up+Sys_Error_Down)/2.;
	Sys_Error2 = Sys_Error * Sys_Error;

	h_sys->SetBinContent(  bin, Sys_Error  );
	h_sys2->SetBinContent( bin, Sys_Error2 );
	cout << "Sys. Error: " << Sys_Error << endl;
	Sum_Errors = Sum_Errors+Sys_Error;
      }

      cout << "-------------" << endl;
      cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
      cout << "-------------" << endl << endl;

      setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
      setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
      h_style = h_sys;

      title.Append("_sys_error_JES");
      h_sys->SetName(title);
      Hlist.Add(h_sys);

      sum_sys_errors[j]->Add(h_sys2);
      sum_exp_errors[j]->Add(h_sys2);

    }

    leg->AddEntry(h_style,"JES","F");
    cout << "-------------" << endl << endl;

  }


  //  BG -------------------------------------------------

  LineColor = kBlack;
  LineStyle = 1;
  LineWidth = 1;
  //  HistColor = kMagenta;
  HistColor = kPink;

  files[5] = new TFile(inpath.Copy().Append("bg_up.root"));
  files[6] = new TFile(inpath.Copy().Append("bg_down.root"));

  if( files[5] && !files[5]->IsZombie() && files[6] && !files[6]->IsZombie() && isExp ) {

    cout << "   >>>   " << "BACKGROUNDS" << endl;
    cout << "--------------------------" << endl << endl;

    for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

      TString title = ListOfKeys->At(j)->GetName();
      cout << "   >>>   " << title << endl;
      if( title.Contains("PSE_") )  continue;

      Sys_Error      = 0.;
      Sys_Error2     = 0.;
      Sys_Error_Up   = 0.;
      Sys_Error_Down = 0.;

      files[0]->GetObject( title, h_ref);
      files[0]->GetObject( title, h_style);
      files[5]->GetObject( title, h_var_up);
      files[6]->GetObject( title, h_var_down);

      TH1* h_sys       = (TH1*) h_var_up->Clone();   h_sys->Reset();
      TH1* h_sys2      = (TH1*) h_var_down->Clone(); h_sys2->Reset();
      TH1* h_sys_up    = (TH1*) h_var_up->Clone();   h_sys_up->Reset();
      TH1* h_sys_down  = (TH1*) h_var_down->Clone(); h_sys_down->Reset();
      TH1* h_sys2_up   = (TH1*) h_var_up->Clone();   h_sys2_up->Reset();
      TH1* h_sys2_down = (TH1*) h_var_down->Clone(); h_sys2_down->Reset();

      TCanvas* Canvas = new TCanvas("plot", "plot", 800, 800);

      setHistogramStyle(h_ref,      kBlack, 1, 2, HistColor);
      setHistogramStyle(h_var_up,   kBlue,  1, 1, HistColor);
      setHistogramStyle(h_var_down, kRed,   1, 1, HistColor);

      Canvas->Clear();
      h_ref->SetMarkerColor(kBlack);
      h_ref->Draw();
      h_var_up->SetMarkerColor(kBlue);
      h_var_up->Draw("SAME");
      h_var_down->SetMarkerColor(kRed);
      h_var_down->Draw("SAME");
      if( set_LogY ) gPad->SetLogy();
      Canvas->Print(outpath.Copy().Append(title).Append("_BG").Append(outform));

      Int_t N_bins = h_ref->GetNbinsX();
      Double_t Sum_Errors = 0.;

      for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {

	if( h_ref->GetBinContent(bin) == 0 )  continue;

	Sys_Error_Up   = abs(h_ref->GetBinContent(bin) - h_var_up->GetBinContent(bin))/h_ref->GetBinContent(bin);
	Sys_Error_Down = abs(h_ref->GetBinContent(bin) - h_var_down->GetBinContent(bin))/h_ref->GetBinContent(bin);
  
	Sys_Error  = (Sys_Error_Up+Sys_Error_Down)/2.;
	Sys_Error2 = Sys_Error * Sys_Error;

	h_sys->SetBinContent(  bin, Sys_Error  );
	h_sys2->SetBinContent( bin, Sys_Error2 );
	cout << "Sys. Error: " << Sys_Error << endl;
	Sum_Errors = Sum_Errors+Sys_Error;
      }

      cout << "-------------" << endl;
      cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
      cout << "-------------" << endl << endl;

      setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
      setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
      h_style = h_sys;

      title.Append("_sys_error_BG");
      h_sys->SetName(title);
      Hlist.Add(h_sys);

      sum_sys_errors[j]->Add(h_sys2);
      sum_exp_errors[j]->Add(h_sys2);

    }

    leg->AddEntry(h_style,"BG","F");
    cout << "-------------" << endl << endl;

  }

  //  DY ----------------------------------------------

  LineColor = kBlack;
  LineStyle = 1;
  LineWidth = 1;
  HistColor = kAzure+8;

  files[7] = new TFile(inpath.Copy().Append("dy_up.root"));
  files[8] = new TFile(inpath.Copy().Append("dy_down.root"));

  if( files[7] && !files[7]->IsZombie() && files[8] && !files[8]->IsZombie() && isExp ) {

    cout << "   >>>   " << "DRELL-YAN" << endl;
    cout << "------------------------" << endl << endl;

    for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

      TString title = ListOfKeys->At(j)->GetName();
      cout << "   >>>   " << title << endl;
      if( title.Contains("PSE_") )  continue;

      Sys_Error      = 0.;
      Sys_Error2     = 0.;
      Sys_Error_Up   = 0.;
      Sys_Error_Down = 0.;

      files[0]->GetObject( title, h_ref);
      files[0]->GetObject( title, h_style);
      files[7]->GetObject( title, h_var_up);
      files[8]->GetObject( title, h_var_down);

      TH1* h_sys       = (TH1*) h_var_up->Clone();   h_sys->Reset();
      TH1* h_sys2      = (TH1*) h_var_down->Clone(); h_sys2->Reset();
      TH1* h_sys_up    = (TH1*) h_var_up->Clone();   h_sys_up->Reset();
      TH1* h_sys_down  = (TH1*) h_var_down->Clone(); h_sys_down->Reset();
      TH1* h_sys2_up   = (TH1*) h_var_up->Clone();   h_sys2_up->Reset();
      TH1* h_sys2_down = (TH1*) h_var_down->Clone(); h_sys2_down->Reset();

      TCanvas* Canvas = new TCanvas("plot", "plot", 800, 800);

      setHistogramStyle(h_ref,      kBlack, 1, 2, HistColor);
      setHistogramStyle(h_var_up,   kBlue,  1, 1, HistColor);
      setHistogramStyle(h_var_down, kRed,   1, 1, HistColor);

      Canvas->Clear();
      h_ref->SetMarkerColor(kBlack);
      h_ref->Draw();
      h_var_up->SetMarkerColor(kBlue);
      h_var_up->Draw("SAME");
      h_var_down->SetMarkerColor(kRed);
      h_var_down->Draw("SAME");
      if( set_LogY ) gPad->SetLogy();
      Canvas->Print(outpath.Copy().Append(title).Append("_DY").Append(outform));

      Int_t N_bins = h_ref->GetNbinsX();
      Double_t Sum_Errors = 0.;

      for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {

	if( h_ref->GetBinContent(bin) == 0 )  continue;

	Sys_Error_Up   = abs(h_ref->GetBinContent(bin) - h_var_up->GetBinContent(bin))/h_ref->GetBinContent(bin);
	Sys_Error_Down = abs(h_ref->GetBinContent(bin) - h_var_down->GetBinContent(bin))/h_ref->GetBinContent(bin);
  
	Sys_Error  = (Sys_Error_Up+Sys_Error_Down)/2.;
	Sys_Error2 = Sys_Error * Sys_Error;

	h_sys->SetBinContent(  bin, Sys_Error  );
	h_sys2->SetBinContent( bin, Sys_Error2 );
	cout << "Sys. Error: " << Sys_Error << endl;
	Sum_Errors = Sum_Errors+Sys_Error;
      }

      cout << "-------------" << endl;
      cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
      cout << "-------------" << endl << endl;

      setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
      setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
      h_style = h_sys;

      title.Append("_sys_error_DY");
      h_sys->SetName(title);
      Hlist.Add(h_sys);

      sum_sys_errors[j]->Add(h_sys2);
      sum_exp_errors[j]->Add(h_sys2);
      
      dysys[j] = h_sys;
      dysys2[j] = h_sys2;

    }

    leg->AddEntry(h_style,"DY","F");
    cout << "-------------" << endl << endl;

  }


  // BTAG ----------------------------------------------

  LineColor = kBlack;
  LineStyle = 1;
  LineWidth = 1;
  HistColor = kAzure-2;

  if(isExp) {

    cout << "   >>>   " << "B-TAGGING" << endl;
    cout << "------------------------" << endl << endl;

    for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

      TString title = ListOfKeys->At(j)->GetName();
      cout << "   >>>   " << title << endl;
      if( title.Contains("PSE_") )  continue;

      Sys_Error  = 0.;
      Sys_Error2 = 0.;

      files[0]->GetObject( title, h_ref);
      files[0]->GetObject( title, h_style);

      TH1* h_sys  = (TH1*) h_ref->Clone(); h_sys->Reset();
      TH1* h_sys2 = (TH1*) h_ref->Clone(); h_sys2->Reset();

      Int_t N_bins = h_ref->GetNbinsX();
      Double_t Sum_Errors = 0.;

      for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {

	//  TOTAL XS
	if(isTotal) {
	  if( title.Contains("mumu_") )      Sys_Error  = 0.06;
	  if( title.Contains("ee_") )        Sys_Error  = 0.06;
	  if( title.Contains("emu_") )       Sys_Error  = 0.06;
	  if( title.Contains("combined_") )  Sys_Error  = 0.035;
	}

	//  DIFFERENTIAL XS
	if(isDiff) {
	  if( title.Contains("mumu_") )      Sys_Error  = 0.03;
	  if( title.Contains("ee_") )        Sys_Error  = 0.03;
	  if( title.Contains("emu_") )       Sys_Error  = 0.03;
	  if( title.Contains("combined_") )  Sys_Error  = 0.017;
	}

	Sys_Error2 = Sys_Error * Sys_Error;
	h_sys->SetBinContent(  bin, Sys_Error  );
	h_sys2->SetBinContent( bin, Sys_Error2 );
	cout << "Sys. Error: " << Sys_Error << endl;
	Sum_Errors = Sum_Errors+Sys_Error;
      }

      cout << "-------------" << endl;
      cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
      cout << "-------------" << endl << endl;

      setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
      setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
      h_style = h_sys;

      title.Append("_sys_error_BTAG");
      h_sys->SetName(title);
      Hlist.Add(h_sys);

      sum_sys_errors[j]->Add(h_sys2);
      sum_exp_errors[j]->Add(h_sys2);

    }

    leg->AddEntry(h_style,"BTAG","F");
    cout << "-------------" << endl << endl;

  }

  // TRIGGER EFFICIENCIES -----------------------------------------

  LineColor = kBlack;
  LineStyle = 1;
  LineWidth = 1;
  HistColor = kGreen-3;

  if(isExp) {

    cout << "   >>>   " << "TRIGGER Efficiencies" << endl;
    cout << "-----------------------------------" << endl << endl;

    for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

      TString title = ListOfKeys->At(j)->GetName();
      cout << "   >>>   " << title << endl;
      if( title.Contains("PSE_") )  continue;

      Sys_Error  = 0.;
      Sys_Error2 = 0.;

      files[0]->GetObject( title, h_ref);
      files[0]->GetObject( title, h_style);

      TH1* h_sys  = (TH1*) h_ref->Clone(); h_sys->Reset();
      TH1* h_sys2 = (TH1*) h_ref->Clone(); h_sys2->Reset();

      Int_t N_bins = h_ref->GetNbinsX();
      Double_t Sum_Errors = 0.;

      for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {

	//  TOTAL XS
	if(isTotal) {
	  if( title.Contains("mumu_") )      Sys_Error  = 0.015;
	  if( title.Contains("ee_") )        Sys_Error  = 0.025;
	  if( title.Contains("emu_") )       Sys_Error  = 0.020;
	  if( title.Contains("combined_") )  Sys_Error  = 0.011;
	}

	//  DIFFERENTIAL XS
	if(isDiff) {
	  if( title.Contains("mumu_") )      Sys_Error  = 0.015;
	  if( title.Contains("ee_") )        Sys_Error  = 0.025;
	  if( title.Contains("emu_") )       Sys_Error  = 0.020;
	  if( title.Contains("combined_") )  Sys_Error  = 0.011;
	}

	Sys_Error2 = Sys_Error * Sys_Error;
	h_sys->SetBinContent(  bin, Sys_Error  );
	h_sys2->SetBinContent( bin, Sys_Error2 );
	cout << "Sys. Error: " << Sys_Error << endl;
	Sum_Errors = Sum_Errors+Sys_Error;
      }

      cout << "-------------" << endl;
      cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
      cout << "-------------" << endl << endl;

      setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
      setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
      h_style = h_sys;

      title.Append("_sys_error_TRIGGER_EFF");
      h_sys->SetName(title);
      Hlist.Add(h_sys);

      sum_sys_errors[j]->Add(h_sys2);
      sum_exp_errors[j]->Add(h_sys2);

    }

    leg->AddEntry(h_style,"Trigger Eff.","F");
    cout << "-------------" << endl << endl;

  }

  // LEPTON SELECTION MODEL ------------------------------------

  LineColor = kBlack;
  LineStyle = 1;
  LineWidth = 1;
  HistColor = kGreen+3;

  if(isExp) {

    cout << "   >>>   " << "Lepton MODEL" << endl;
    cout << "---------------------------" << endl << endl;

    for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

      TString title = ListOfKeys->At(j)->GetName();
      cout << "   >>>   " << title << endl;
      if( title.Contains("PSE_") )  continue;

      Sys_Error  = 0.;
      Sys_Error2 = 0.;

      files[0]->GetObject( title, h_ref);
      files[0]->GetObject( title, h_style);

      TH1* h_sys  = (TH1*) h_ref->Clone(); h_sys->Reset();
      TH1* h_sys2 = (TH1*) h_ref->Clone(); h_sys2->Reset();

      //      Int_t N_bins = h_ref->GetNbinsX();
      Double_t Sum_Errors = 0.;

      for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {

	//  TOTAL XS
	if(isTotal) {
	  if( title.Contains("mumu_") )      Sys_Error  = 0.04;
	  if( title.Contains("ee_") )        Sys_Error  = 0.04;
	  if( title.Contains("emu_") )       Sys_Error  = 0.04;
	  if( title.Contains("combined_") )  Sys_Error  = 0.023;
	}

	//  DIFFERENTIAL XS
	if(isDiff) {
	  if( title.Contains("mumu_") )      Sys_Error  = 0.02;
	  if( title.Contains("ee_") )        Sys_Error  = 0.02;
	  if( title.Contains("emu_") )       Sys_Error  = 0.02;
	  if( title.Contains("combined_") )  Sys_Error  = 0.012;
	}

	Sys_Error2 = Sys_Error * Sys_Error;
	h_sys->SetBinContent(  bin, Sys_Error  );
	h_sys2->SetBinContent( bin, Sys_Error2 );
	cout << "Sys. Error: " << Sys_Error << endl;
	Sum_Errors = Sum_Errors+Sys_Error;
      }

      setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
      setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
      h_style = h_sys;

      title.Append("_sys_error_LEPTON_MODEL");
      h_sys->SetName(title);
      Hlist.Add(h_sys);

      sum_sys_errors[j]->Add(h_sys2);
      sum_exp_errors[j]->Add(h_sys2);

    }

    leg->AddEntry(h_style,"Lepton Model","F");
    cout << "-------------" << endl << endl;

  }

  // PILE-UP ----------------------------------------------

  LineColor = kBlack;
  LineStyle = 1;
  LineWidth = 1;
  HistColor = kYellow-6;

  files[9]  = new TFile(inpath.Copy().Append("pu_up.root"));
  files[10] = new TFile(inpath.Copy().Append("pu_down.root"));

  if( files[9] && !files[9]->IsZombie() && files[10] && !files[10]->IsZombie() && isExp ) {

    cout << "   >>>   " << "PILE-UP" << endl;
    cout << "----------------------" << endl << endl;

    for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

      TString title = ListOfKeys->At(j)->GetName();
      cout << "   >>>   " << title << endl;
      if( title.Contains("PSE_") )  continue;

      Sys_Error      = 0.;
      Sys_Error2     = 0.;
      Sys_Error_Up   = 0.;
      Sys_Error_Down = 0.;

      files[0]->GetObject(  title, h_ref);
      files[0]->GetObject(  title, h_style);
      files[9]->GetObject(  title, h_var_up);
      files[10]->GetObject( title, h_var_down);

      TH1* h_sys       = (TH1*) h_var_up->Clone();   h_sys->Reset();
      TH1* h_sys2      = (TH1*) h_var_down->Clone(); h_sys2->Reset();
      TH1* h_sys_up    = (TH1*) h_var_up->Clone();   h_sys_up->Reset();
      TH1* h_sys_down  = (TH1*) h_var_down->Clone(); h_sys_down->Reset();
      TH1* h_sys2_up   = (TH1*) h_var_up->Clone();   h_sys2_up->Reset();
      TH1* h_sys2_down = (TH1*) h_var_down->Clone(); h_sys2_down->Reset();

      TCanvas* Canvas = new TCanvas("plot", "plot", 800, 800);

      setHistogramStyle(h_ref,      kBlack, 1, 2, HistColor);
      setHistogramStyle(h_var_up,   kBlue,  1, 1, HistColor);
      setHistogramStyle(h_var_down, kRed,   1, 1, HistColor);

      Canvas->Clear();
      h_ref->SetMarkerColor(kBlack);
      h_ref->Draw();
      h_var_up->SetMarkerColor(kBlue);
      h_var_up->Draw("SAME");
      h_var_down->SetMarkerColor(kRed);
      h_var_down->Draw("SAME");
      if( set_LogY ) gPad->SetLogy();
      Canvas->Print(outpath.Copy().Append(title).Append("_PILE_UP").Append(outform));

      Int_t N_bins = h_ref->GetNbinsX();
      Double_t Sum_Errors = 0.;

      for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {

	if( h_ref->GetBinContent(bin) == 0 )  continue;

	Sys_Error_Up   = abs(h_ref->GetBinContent(bin) - h_var_up->GetBinContent(bin))/h_ref->GetBinContent(bin);
	Sys_Error_Down = abs(h_ref->GetBinContent(bin) - h_var_down->GetBinContent(bin))/h_ref->GetBinContent(bin);
  
	Sys_Error  = (Sys_Error_Up+Sys_Error_Down)/2.;
	Sys_Error2 = Sys_Error * Sys_Error;

	h_sys->SetBinContent(  bin, Sys_Error  );
	h_sys2->SetBinContent( bin, Sys_Error2 );
	cout << "Sys. Error: " << Sys_Error << endl;
	Sum_Errors = Sum_Errors+Sys_Error;
      }

      cout << "-------------" << endl;
      cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
      cout << "-------------" << endl << endl;

      setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
      setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
      h_style = h_sys;

      title.Append("_sys_error_PILE_UP");
      h_sys->SetName(title);
      Hlist.Add(h_sys);

      sum_sys_errors[j]->Add(h_sys2);
      sum_exp_errors[j]->Add(h_sys2);

    }

    leg->AddEntry(h_style,"Pile-Up","F");
    cout << "-------------" << endl << endl;

  }

  // KIN FIT !!!(OLD)!!!  ----------------------------------------------

//   LineColor = kBlack;
//   LineStyle = 1;
//   LineWidth = 1;
//   //  HistColor = kMagenta-6;
//   HistColor = kYellow;

//   if(isExp) {

//     cout << "   >>>   " << "KIN FIT" << endl;
//     cout << "----------------------" << endl << endl;

//     for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

//       TString title = ListOfKeys->At(j)->GetName();
//       cout << "   >>>   " << title << endl;
//       if( title.Contains("PSE_") )  continue;

//       Sys_Error  = 0.;
//       Sys_Error2 = 0.;

//       files[0]->GetObject( title, h_ref);
//       files[0]->GetObject( title, h_style);

//       TH1* h_sys  = (TH1*) h_ref->Clone(); h_sys->Reset();
//       TH1* h_sys2 = (TH1*) h_ref->Clone(); h_sys2->Reset();

//       Int_t N_bins = h_ref->GetNbinsX();
//       Double_t Sum_Errors = 0.;

//       for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {

// 	//  TOTAL XS
// 	if(isTotal) {
// 	  if( title.Contains("mumu_") )      Sys_Error  = 0.03;
// 	  if( title.Contains("ee_") )        Sys_Error  = 0.03;
// 	  if( title.Contains("emu_") )       Sys_Error  = 0.03;
// 	  if( title.Contains("combined_") )  Sys_Error  = 0.017;
// 	}

// 	//  DIFFERENTIAL XS
// 	if(isDiff) {
// 	  if( title.Contains("mumu_") )      Sys_Error  = 0.03;
// 	  if( title.Contains("ee_") )        Sys_Error  = 0.03;
// 	  if( title.Contains("emu_") )       Sys_Error  = 0.03;
// 	  if( title.Contains("combined_") )  Sys_Error  = 0.017;
// 	}

// 	Sys_Error2 = Sys_Error * Sys_Error;
// 	h_sys->SetBinContent(  bin, Sys_Error  );
// 	h_sys2->SetBinContent( bin, Sys_Error2 );
// 	cout << "Sys. Error: " << Sys_Error << endl;
// 	Sum_Errors = Sum_Errors+Sys_Error;
//       }

//       cout << "-------------" << endl;
//       cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
//       cout << "-------------" << endl << endl;

//       setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
//       setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
//       h_style = h_sys;

//       title.Append("_sys_error_KIN_FIT");
//       h_sys->SetName(title);
//       Hlist.Add(h_sys);

//       sum_sys_errors[j]->Add(h_sys2);
//       sum_exp_errors[j]->Add(h_sys2);

//     }

//     leg->AddEntry(h_style,"KIN FIT","F");
//     cout << "-------------" << endl << endl;

//   }

  //  KIN FIT  !!!(NEW)!!! --------------------------------------------------

  LineColor = kBlack;
  LineStyle = 1;
  LineWidth = 1;
  //  HistColor = kMagenta-6;
  HistColor = kYellow;

  files[11] = new TFile(inpath.Copy().Append("no_kin_scale.root"));
  files[12] = new TFile(inpath.Copy().Append("standard.root"));

  //  files[11] = new TFile(inpath.Copy().Append("standard.root"));
  //  files[12] = new TFile(inpath.Copy().Append("pythiaBug.root"));
  //  files[12] = new TFile(inpath.Copy().Append("mll30.root"));

  //  files[11] = new TFile(inpath.Copy().Append("standard.root"));
  //  files[12] = new TFile(inpath.Copy().Append("kin_scale2.root"));

  if( files[11] && !files[11]->IsZombie() &&  files[12] && !files[12]->IsZombie() && isExp ) {

    cout << "   >>>   " << "KIN FIT" << endl;
    cout << "----------------------" << endl << endl;

    for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

      TString title = ListOfKeys->At(j)->GetName();
      cout << "   >>>   " << title << endl;
      if( title.Contains("PSE_") )    continue;

      Sys_Error  = 0.;
      Sys_Error2 = 0.;

      files[11]->GetObject( title, h_ref);
      files[11]->GetObject( title, h_style);
      files[11]->GetObject( title, h_var_up);
      files[12]->GetObject( title, h_var_down);

      TH1* h_sys       = (TH1*) h_var_up->Clone();   h_sys->Reset();
      TH1* h_sys2      = (TH1*) h_var_down->Clone(); h_sys2->Reset();
      TH1* h_sys_up    = (TH1*) h_var_up->Clone();   h_sys_up->Reset();
      TH1* h_sys_down  = (TH1*) h_var_down->Clone(); h_sys_down->Reset();
      TH1* h_sys2_up   = (TH1*) h_var_up->Clone();   h_sys2_up->Reset();
      TH1* h_sys2_down = (TH1*) h_var_down->Clone(); h_sys2_down->Reset();

      TCanvas* Canvas = new TCanvas("plot", "plot", 800, 800);

      setHistogramStyle(h_ref,      kBlack, 1, 2, HistColor);
      setHistogramStyle(h_var_up,   kBlue,  1, 1, HistColor);
      setHistogramStyle(h_var_down, kRed,   1, 1, HistColor);

      Canvas->Clear();
      //      h_ref->SetMarkerColor(kBlack);
      //      h_ref->Draw();
      h_var_up->SetMarkerColor(kBlue);
      h_var_up->Draw();
      h_var_down->SetMarkerColor(kRed);
      h_var_down->Draw("SAME");
      if( set_LogY ) gPad->SetLogy();
      Canvas->Print(outpath.Copy().Append(title).Append("_KIN_FIT").Append(outform));

      Int_t N_bins = h_ref->GetNbinsX();
      Double_t Sum_Errors = 0.;

      for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {

	if( h_var_up->GetBinContent(bin) == 0 || h_var_down->GetBinContent(bin) == 0 )  continue;

	Sys_Error  = abs(h_var_up->GetBinContent(bin) - h_var_down->GetBinContent(bin))/h_var_up->GetBinContent(bin);
	Sys_Error2 = Sys_Error * Sys_Error;

	h_sys->SetBinContent(  bin, Sys_Error  );
	h_sys2->SetBinContent( bin, Sys_Error2 );
	cout << "Sys. Error: " << Sys_Error << endl;
	Sum_Errors = Sum_Errors+Sys_Error;
      }

      cout << "-------------" << endl;
      cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
      cout << "-------------" << endl << endl;

      setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
      setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
      h_style = h_sys;

      title.Append("_sys_error_KIN_FIT");
      h_sys->SetName(title);
      Hlist.Add(h_sys);

      sum_sys_errors[j]->Add(h_sys2);
      sum_mod_errors[j]->Add(h_sys2);

    }

    leg->AddEntry(h_style,"KIN FIT","F");
    cout << "-------------" << endl << endl;

  }





  LineColor = kGray;
  LineStyle = 1;
  LineWidth = 1;
  //  HistColor = kMagenta-6;
  HistColor = kGray;

  files[11] = new TFile(inpath.Copy().Append("mll30_new.root"));
  files[12] = new TFile(inpath.Copy().Append("standard_new.root"));

  //  files[11] = new TFile(inpath.Copy().Append("standard.root"));
  //  files[12] = new TFile(inpath.Copy().Append("pythiaBug.root"));
  //  files[12] = new TFile(inpath.Copy().Append("mll30.root"));

  //  files[11] = new TFile(inpath.Copy().Append("standard.root"));
  //  files[12] = new TFile(inpath.Copy().Append("kin_scale2.root"));

  if( files[11] && !files[11]->IsZombie() &&  files[12] && !files[12]->IsZombie() && isExp ) {

    cout << "   >>>   " << "MLL 30 " << endl;
    cout << "----------------------" << endl << endl;

    for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

      TString title = ListOfKeys->At(j)->GetName();
      cout << "   >>>   " << title << endl;
      if( title.Contains("PSE_") )    continue;

      Sys_Error  = 0.;
      Sys_Error2 = 0.;

      files[11]->GetObject( title, h_ref);
      files[11]->GetObject( title, h_style);
      files[11]->GetObject( title, h_var_up);
      files[12]->GetObject( title, h_var_down);

      TH1* h_sys       = (TH1*) h_var_up->Clone();   h_sys->Reset();
      TH1* h_sys2      = (TH1*) h_var_down->Clone(); h_sys2->Reset();
      TH1* h_sys_up    = (TH1*) h_var_up->Clone();   h_sys_up->Reset();
      TH1* h_sys_down  = (TH1*) h_var_down->Clone(); h_sys_down->Reset();
      TH1* h_sys2_up   = (TH1*) h_var_up->Clone();   h_sys2_up->Reset();
      TH1* h_sys2_down = (TH1*) h_var_down->Clone(); h_sys2_down->Reset();

      TCanvas* Canvas = new TCanvas("plot", "plot", 800, 800);

      setHistogramStyle(h_ref,      kBlack, 1, 2, HistColor);
      setHistogramStyle(h_var_up,   kBlue,  1, 1, HistColor);
      setHistogramStyle(h_var_down, kRed,   1, 1, HistColor);

      Canvas->Clear();
      //      h_ref->SetMarkerColor(kBlack);
      //      h_ref->Draw();
      h_var_up->SetMarkerColor(kBlue);
      h_var_up->Draw();
      h_var_down->SetMarkerColor(kRed);
      h_var_down->Draw("SAME");
      if( set_LogY ) gPad->SetLogy();
      Canvas->Print(outpath.Copy().Append(title).Append("_MLL_30").Append(outform));

      Int_t N_bins = h_ref->GetNbinsX();
      Double_t Sum_Errors = 0.;

      if(title.Contains("Leptons_Eta") || title.Contains("TopQuarks_Rapidity") ) {

        SymmetricAroundZero(h_ref, h_var_down, h_var_down, h_sys, h_sys2, Sum_Errors);

      } else {

        for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {

            if( h_var_up->GetBinContent(bin) == 0 || h_var_down->GetBinContent(bin) == 0 )  continue;

            Sys_Error  = fabs(h_var_up->GetBinContent(bin) - h_var_down->GetBinContent(bin))/h_var_up->GetBinContent(bin);
            Sys_Error2 = Sys_Error * Sys_Error;

            h_sys->SetBinContent(  bin, Sys_Error  );
            h_sys2->SetBinContent( bin, Sys_Error2 );
            cout << "Sys. Error: " << Sys_Error << endl;
            Sum_Errors = Sum_Errors+Sys_Error;
        }
      }
      
      cout << "-------------" << endl;
      cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
      cout << "-------------" << endl << endl;

      setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
      setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
      h_style = h_sys;

      title.Append("_sys_error_MLL_30");
      h_sys->SetName(title);

      //dont add this as own error, but rather add to DY
//       Hlist.Add(h_sys);
//       sum_sys_errors[j]->Add(h_sys2);
//       sum_exp_errors[j]->Add(h_sys2);

      //loop over all bins
      dysys2[j]->Add(h_sys2);
      for (int i = 1; i <= dysys[j]->GetNbinsX(); ++i) {
          dysys[j]->SetBinContent(i, TMath::Sqrt(dysys2[j]->GetBinContent(i)));
      }
      

    }

    //leg->AddEntry(h_style,"MLL 30","F");
    cout << "-------------" << endl << endl;

  }










  // MODEL ----------------------------------------------

//   LineColor = kBlack;
//   LineStyle = 1;
//   LineWidth = 1;
//   HistColor = kMagenta+2;

//   cout << "   >>>   " << "MODEL" << endl;
//   cout << "--------------------" << endl << endl;

//   if (isModel) {

//     for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

//       TString title = ListOfKeys->At(j)->GetName();
//       cout << "   >>>   " << title << endl;
//       if( !(isPSE) &&   title.Contains("PSE_")  )  continue;
//       if(   isPSE  && !(title.Contains("PSE_")) )  continue;

//       Sys_Error  = 0.;
//       Sys_Error2 = 0.;

//       files[0]->GetObject( title, h_ref);
//       files[0]->GetObject( title, h_style);

//       TH1* h_sys  = (TH1*) h_ref->Clone(); h_sys->Reset();
//       TH1* h_sys2 = (TH1*) h_ref->Clone(); h_sys2->Reset();

//       Int_t N_bins = h_ref->GetNbinsX();
//       Double_t Sum_Errors = 0.;

//       for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {

// 	//  TOTAL XS
// 	if(isTotal) {
// 	  if( title.Contains("mumu_") )      Sys_Error  = 0.052;
// 	  if( title.Contains("ee_") )        Sys_Error  = 0.052;
// 	  if( title.Contains("emu_") )       Sys_Error  = 0.031;
// 	  if( title.Contains("combined_") )  Sys_Error  = 0.024;
// 	}

// 	//  DIFFERENTIAL XS
// 	if(isDiff) {
// 	  if( title.Contains("mumu_") )      Sys_Error  = 0.022;
// 	  if( title.Contains("ee_") )        Sys_Error  = 0.022;
// 	  if( title.Contains("emu_") )       Sys_Error  = 0.014;
// 	  if( title.Contains("combined_") )  Sys_Error  = 0.010;
// 	}

// 	Sys_Error2 = Sys_Error * Sys_Error;
// 	h_sys->SetBinContent(  bin, Sys_Error  );
// 	h_sys2->SetBinContent( bin, Sys_Error2 );
// 	cout << "Sys. Error: " << Sys_Error << endl;
// 	Sum_Errors = Sum_Errors+Sys_Error;
//       }

//       cout << "-------------" << endl;
//       cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
//       cout << "-------------" << endl << endl;

//       setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
//       setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
//       h_style = h_sys;

//       title.Append("_sys_error_MODEL");
//       h_sys->SetName(title);
//       Hlist.Add(h_sys);

//       sum_sys_errors[j]->Add(h_sys2);
//       sum_mod_errors[j]->Add(h_sys2);

//     }

//     leg->AddEntry(h_style,"MODEL","F");
//     cout << "-------------" << endl << endl;

//   }

  // ISR/FSR ----------------------------------------------

//   LineColor = kBlack;
//   LineStyle = 1;
//   LineWidth = 1;
//   HistColor = kMagenta+4;

//   files[0]  = new TFile(inpath.Copy().Append("madgraph_Fall10.root"));
//   files[13] = new TFile(inpath.Copy().Append("ISRFSR_up_Fall10.root"));
//   files[14] = new TFile(inpath.Copy().Append("ISRFSR_down_Fall10.root"));

//   if( files[13] && !files[13]->IsZombie() && files[14] && !files[14]->IsZombie() && isModel ) {

//     cout << "   >>>   " << "ISR/FSR" << endl;
//     cout << "----------------------" << endl << endl;

//     for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

//       TString title = ListOfKeys->At(j)->GetName();
//       cout << "   >>>   " << title << endl;
//       if( !(isPSE) &&   title.Contains("PSE_")  )  continue;
//       if(   isPSE  && !(title.Contains("PSE_")) )  continue;

//       Sys_Error      = 0.;
//       Sys_Error2     = 0.;
//       Sys_Error_Up   = 0.;
//       Sys_Error_Down = 0.;

//       files[0]->GetObject(  title, h_ref);
//       files[0]->GetObject(  title, h_style);
//       files[13]->GetObject( title, h_var_up);
//       files[14]->GetObject( title, h_var_down);

//       TH1* h_sys       = (TH1*) h_var_up->Clone();   h_sys->Reset();
//       TH1* h_sys2      = (TH1*) h_var_down->Clone(); h_sys2->Reset();
//       TH1* h_sys_up    = (TH1*) h_var_up->Clone();   h_sys_up->Reset();
//       TH1* h_sys_down  = (TH1*) h_var_down->Clone(); h_sys_down->Reset();
//       TH1* h_sys2_up   = (TH1*) h_var_up->Clone();   h_sys2_up->Reset();
//       TH1* h_sys2_down = (TH1*) h_var_down->Clone(); h_sys2_down->Reset();

//       TCanvas* Canvas = new TCanvas("plot", "plot", 800, 800);

//       setHistogramStyle(h_ref,      kBlack, 1, 2, HistColor);
//       setHistogramStyle(h_var_up,   kBlue,  1, 1, HistColor);
//       setHistogramStyle(h_var_down, kRed,   1, 1, HistColor);

//       Canvas->Clear();
//       h_ref->SetMarkerColor(kBlack);
//       h_ref->Draw();
//       h_var_up->SetMarkerColor(kBlue);
//       h_var_up->Draw("SAME");
//       h_var_down->SetMarkerColor(kRed);
//       h_var_down->Draw("SAME");
//       if( !(isPSE) && set_LogY ) gPad->SetLogy();
//       Canvas->Print(outpath.Copy().Append(title).Append("_ISR_FSR").Append(outform));

//       Int_t N_bins = h_ref->GetNbinsX();
//       Double_t Sum_Errors = 0.;

//       //  Leptons Eta  ---------------------------------------

//       if(title.Contains("Leptons_Eta") || title.Contains("TopQuarks_Rapidity")){
//           SymmetricAroundZero(h_ref, h_var_up, h_var_down, h_sys, h_sys2, Sum_Errors);
//       }

//       else {

// 	for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {

// 	  if( h_ref->GetBinContent(bin) == 0 )  continue;

// 	  Sys_Error_Up   = abs(h_ref->GetBinContent(bin) - h_var_up->GetBinContent(bin))/h_ref->GetBinContent(bin);
// 	  Sys_Error_Down = abs(h_ref->GetBinContent(bin) - h_var_down->GetBinContent(bin))/h_ref->GetBinContent(bin);
  
// 	  //	  Sys_Error  = TMath::Max ( Sys_Error_Up, Sys_Error_Down);
// 	  Sys_Error  = (Sys_Error_Up+Sys_Error_Down)/2.;
// 	  Sys_Error2 = Sys_Error * Sys_Error;

// 	  h_sys->SetBinContent(  bin, Sys_Error  );
// 	  h_sys2->SetBinContent( bin, Sys_Error2 );
// 	  cout << "Sys. Error: " << Sys_Error << endl;
// 	  Sum_Errors = Sum_Errors+Sys_Error;

// 	}

//       }

//       cout << "-------------" << endl;
//       cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
//       cout << "-------------" << endl << endl;

//       setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
//       setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
//       h_style = h_sys;

//       title.Append("_sys_error_ISR_FSR");
//       h_sys->SetName(title);
//       Hlist.Add(h_sys);

//       sum_sys_errors[j]->Add(h_sys2);
//       sum_mod_errors[j]->Add(h_sys2);

//     }

//     leg->AddEntry(h_style,"ISR/FSR","F");
//     cout << "-------------" << endl << endl;

//   }

  // Q^2 Scale ----------------------------------------------

  LineColor = kBlack;
  LineStyle = 1;
  LineWidth = 1;
  HistColor = kMagenta-4;

  files[0]  = new TFile(inpath.Copy().Append("madgraph_Fall10.root"));
  files[15] = new TFile(inpath.Copy().Append("scaleup_Fall10.root"));
  files[16] = new TFile(inpath.Copy().Append("scaledown_Fall10.root"));

  if( files[15] && !files[15]->IsZombie() && files[16] && !files[16]->IsZombie() && isModel ) {

    cout << "   >>>   " << "Q scale" << endl;
    cout << "----------------------" << endl << endl;

    for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

      TString title = ListOfKeys->At(j)->GetName();
      cout << "   >>>   " << title << endl;
      if( !(isPSE) &&   title.Contains("PSE_")  )  continue;
      if(   isPSE  && !(title.Contains("PSE_")) )  continue;

      Sys_Error      = 0.;
      Sys_Error2     = 0.;
      Sys_Error_Up   = 0.;
      Sys_Error_Down = 0.;

      files[0]->GetObject(  title, h_ref);
      files[0]->GetObject(  title, h_style);
      files[15]->GetObject( title, h_var_up);
      files[16]->GetObject( title, h_var_down);

      TH1* h_sys       = (TH1*) h_var_up->Clone();   h_sys->Reset();
      TH1* h_sys2      = (TH1*) h_var_down->Clone(); h_sys2->Reset();
      TH1* h_sys_up    = (TH1*) h_var_up->Clone();   h_sys_up->Reset();
      TH1* h_sys_down  = (TH1*) h_var_down->Clone(); h_sys_down->Reset();
      TH1* h_sys2_up   = (TH1*) h_var_up->Clone();   h_sys2_up->Reset();
      TH1* h_sys2_down = (TH1*) h_var_down->Clone(); h_sys2_down->Reset();

      TCanvas* Canvas = new TCanvas("plot", "plot", 800, 800);

      setHistogramStyle(h_ref,      kBlack, 1, 2, HistColor);
      setHistogramStyle(h_var_up,   kBlue,  1, 1, HistColor);
      setHistogramStyle(h_var_down, kRed,   1, 1, HistColor);

      Canvas->Clear();
      h_ref->SetMarkerColor(kBlack);
      h_ref->Draw();
      h_var_up->SetMarkerColor(kBlue);
      h_var_up->Draw("SAME");
      h_var_down->SetMarkerColor(kRed);
      h_var_down->Draw("SAME");
      if( !(isPSE) && set_LogY ) gPad->SetLogy();
      Canvas->Print(outpath.Copy().Append(title).Append("_SCALE").Append(outform));

      Int_t N_bins = h_ref->GetNbinsX();
      Double_t Sum_Errors = 0.;

      // Symmetrise Eta / Rapidity  ---------------------------------------

      if( title.Contains("Leptons_Eta") || title.Contains("TopQuarks_Rapidity") ) {

          SymmetricAroundZero(h_ref, h_var_up, h_var_down, h_sys, h_sys2, Sum_Errors);

      }

      else {

	for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {

	  if( h_ref->GetBinContent(bin) == 0 )  continue;

	  Sys_Error_Up   = abs(h_ref->GetBinContent(bin) - h_var_up->GetBinContent(bin))/h_ref->GetBinContent(bin);
	  Sys_Error_Down = abs(h_ref->GetBinContent(bin) - h_var_down->GetBinContent(bin))/h_ref->GetBinContent(bin);
  
	  //	  Sys_Error  = TMath::Max ( Sys_Error_Up, Sys_Error_Down);
	  Sys_Error  = (Sys_Error_Up+Sys_Error_Down)/2.;
	  Sys_Error2 = Sys_Error * Sys_Error;

	  h_sys->SetBinContent(  bin, Sys_Error  );
	  h_sys2->SetBinContent( bin, Sys_Error2 );
	  cout << "Sys. Error: " << Sys_Error << endl;
	  Sum_Errors = Sum_Errors+Sys_Error;

	}

      }

      cout << "-------------" << endl;
      cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
      cout << "-------------" << endl << endl;

      setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
      setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
      h_style = h_sys;

      title.Append("_sys_error_SCALE");
      h_sys->SetName(title);
      Hlist.Add(h_sys);

      sum_sys_errors[j]->Add(h_sys2);
      sum_mod_errors[j]->Add(h_sys2);

    }

    leg->AddEntry(h_style,"Q2 scale","F");
    cout << "-------------" << endl << endl;

  }


  // Top mass ----------------------------------------------

  LineColor = kBlack;
  LineStyle = 1;
  LineWidth = 1;
  HistColor = kMagenta-2;

  files[0]  = new TFile(inpath.Copy().Append("madgraph_Fall10.root"));
  files[17] = new TFile(inpath.Copy().Append("mass178_5_Fall10.root"));
  files[18] = new TFile(inpath.Copy().Append("mass166_5_Fall10.root"));

  if( files[17] && !files[17]->IsZombie() && files[18] && !files[18]->IsZombie() && isModel ) {

    cout << "   >>>   " << "Top mass" << endl;
    cout << "-----------------------" << endl << endl;

    for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

      TString title = ListOfKeys->At(j)->GetName();
      cout << "   >>>   " << title << endl;
      if( !(isPSE) &&   title.Contains("PSE_")  )  continue;
      if(   isPSE  && !(title.Contains("PSE_")) )  continue;

      Sys_Error      = 0.;
      Sys_Error2     = 0.;
      Sys_Error_Up   = 0.;
      Sys_Error_Down = 0.;

      files[0]->GetObject(  title, h_ref);
      files[0]->GetObject(  title, h_style);
      files[17]->GetObject( title, h_var_up);
      files[18]->GetObject( title, h_var_down);

      TH1* h_sys       = (TH1*) h_var_up->Clone();   h_sys->Reset();
      TH1* h_sys2      = (TH1*) h_var_down->Clone(); h_sys2->Reset();
      TH1* h_sys_up    = (TH1*) h_var_up->Clone();   h_sys_up->Reset();
      TH1* h_sys_down  = (TH1*) h_var_down->Clone(); h_sys_down->Reset();
      TH1* h_sys2_up   = (TH1*) h_var_up->Clone();   h_sys2_up->Reset();
      TH1* h_sys2_down = (TH1*) h_var_down->Clone(); h_sys2_down->Reset();

      TCanvas* Canvas = new TCanvas("plot", "plot", 800, 800);

      setHistogramStyle(h_ref,      kBlack, 1, 2, HistColor);
      setHistogramStyle(h_var_up,   kBlue,  1, 1, HistColor);
      setHistogramStyle(h_var_down, kRed,   1, 1, HistColor);

      Canvas->Clear();
      h_ref->SetMarkerColor(kBlack);
      h_ref->Draw();
      h_var_up->SetMarkerColor(kBlue);
      h_var_up->Draw("SAME");
      h_var_down->SetMarkerColor(kRed);
      h_var_down->Draw("SAME");
      if( !(isPSE) && set_LogY ) gPad->SetLogy();
      Canvas->Print(outpath.Copy().Append(title).Append("_MASS").Append(outform));

      Int_t N_bins = h_ref->GetNbinsX();
      Double_t Sum_Errors = 0.;

      // Symmetrise Eta / Rapidity  ---------------------------------------

      if( title.Contains("Leptons_Eta") || title.Contains("TopQuarks_Rapidity") ) {

        for(Int_t bin = 1; bin <= h_ref->GetNbinsX()/2.; ++bin) {

	  if( h_ref->GetBinContent(bin) == 0 )  continue;

	  Double_t PlusMinus_Average_Ref  = 0.;
	  Double_t PlusMinus_Average_Up   = 0.;
	  Double_t PlusMinus_Average_Down = 0.;

	  PlusMinus_Average_Ref  = (h_ref->GetBinContent(bin)      + h_ref->GetBinContent(N_bins+1-bin)     ) / 2.;
	  PlusMinus_Average_Up   = (h_var_up->GetBinContent(bin)   + h_var_up->GetBinContent(N_bins+1-bin)  ) / 2.;
	  PlusMinus_Average_Down = (h_var_down->GetBinContent(bin) + h_var_down->GetBinContent(N_bins+1-bin)) / 2.;

	  Sys_Error_Up   = abs(PlusMinus_Average_Ref - PlusMinus_Average_Up  ) / PlusMinus_Average_Ref;
	  Sys_Error_Down = abs(PlusMinus_Average_Ref - PlusMinus_Average_Down) / PlusMinus_Average_Ref;

	  // ++++++++ WARNING ++++++++ Divide by 5
	  Sys_Error  = (Sys_Error_Up+Sys_Error_Down)/(2.*5.);
	  // ++++++++ WARNING ++++++++ Divide by 5
	  Sys_Error2 = Sys_Error * Sys_Error;

	  h_sys->SetBinContent(           bin, Sys_Error  );
	  h_sys->SetBinContent(  N_bins+1-bin, Sys_Error  );
	  h_sys2->SetBinContent(          bin, Sys_Error2 );
	  h_sys2->SetBinContent( N_bins+1-bin, Sys_Error2 );

	  cout << "Sys. Error: " << Sys_Error << endl;
	  Sum_Errors += (2.*Sys_Error);

	}

      }

      else {

	for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {

	  if( h_ref->GetBinContent(bin) == 0 )  continue;

	  Sys_Error_Up   = abs(h_ref->GetBinContent(bin) - h_var_up->GetBinContent(bin))/h_ref->GetBinContent(bin);
	  Sys_Error_Down = abs(h_ref->GetBinContent(bin) - h_var_down->GetBinContent(bin))/h_ref->GetBinContent(bin);

	  // ++++++++ WARNING ++++++++ Divide by 5
	  Sys_Error  = (Sys_Error_Up+Sys_Error_Down)/(2.*5.);
	  // ++++++++ WARNING ++++++++ Divide by 5
	  Sys_Error2 = Sys_Error * Sys_Error;

	  h_sys->SetBinContent(  bin, Sys_Error  );
	  h_sys2->SetBinContent( bin, Sys_Error2 );
	  cout << "Sys. Error: " << Sys_Error << endl;
	  Sum_Errors = Sum_Errors+Sys_Error;

	}

      }

      cout << "-------------" << endl;
      cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
      cout << "-------------" << endl << endl;

      setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
      setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
      h_style = h_sys;

      title.Append("_sys_error_MASS");
      h_sys->SetName(title);
      Hlist.Add(h_sys);

      sum_sys_errors[j]->Add(h_sys2);
      sum_mod_errors[j]->Add(h_sys2);

    }

    leg->AddEntry(h_style,"Top mass","F");
    cout << "-------------" << endl << endl;

  }

  // ME/PS matching ----------------------------------------------

  LineColor = kBlack;
  LineStyle = 1;
  LineWidth = 1;
  HistColor = kMagenta+4;

  files[0]  = new TFile(inpath.Copy().Append("madgraph_Fall10.root"));
  files[19] = new TFile(inpath.Copy().Append("matchingup_Fall10.root"));
  files[20] = new TFile(inpath.Copy().Append("matchingdown_Fall10.root"));

  if( files[19] && !files[19]->IsZombie() && files[20] && !files[20]->IsZombie() && isModel ) {

    cout << "   >>>   " << "ME/PS Matching" << endl;
    cout << "-----------------------------" << endl << endl;

    for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

      TString title = ListOfKeys->At(j)->GetName();
      cout << "   >>>   " << title << endl;
      if( !(isPSE) &&   title.Contains("PSE_")  )  continue;
      if(   isPSE  && !(title.Contains("PSE_")) )  continue;

      Sys_Error      = 0.;
      Sys_Error2     = 0.;
      Sys_Error_Up   = 0.;
      Sys_Error_Down = 0.;

      files[0]->GetObject(  title, h_ref);
      files[0]->GetObject(  title, h_style);
      files[19]->GetObject( title, h_var_up);
      files[20]->GetObject( title, h_var_down);

      TH1* h_sys       = (TH1*) h_var_up->Clone();   h_sys->Reset();
      TH1* h_sys2      = (TH1*) h_var_down->Clone(); h_sys2->Reset();
      TH1* h_sys_up    = (TH1*) h_var_up->Clone();   h_sys_up->Reset();
      TH1* h_sys_down  = (TH1*) h_var_down->Clone(); h_sys_down->Reset();
      TH1* h_sys2_up   = (TH1*) h_var_up->Clone();   h_sys2_up->Reset();
      TH1* h_sys2_down = (TH1*) h_var_down->Clone(); h_sys2_down->Reset();

      TCanvas* Canvas = new TCanvas("plot", "plot", 800, 800);

      setHistogramStyle(h_ref,      kBlack, 1, 2, HistColor);
      setHistogramStyle(h_var_up,   kBlue,  1, 1, HistColor);
      setHistogramStyle(h_var_down, kRed,   1, 1, HistColor);

      Canvas->Clear();
      h_ref->SetMarkerColor(kBlack);
      h_ref->Draw();
      h_var_up->SetMarkerColor(kBlue);
      h_var_up->Draw("SAME");
      h_var_down->SetMarkerColor(kRed);
      h_var_down->Draw("SAME");
      if( !(isPSE) && set_LogY ) gPad->SetLogy();
      Canvas->Print(outpath.Copy().Append(title).Append("_MATCH").Append(outform));

      Int_t N_bins = h_ref->GetNbinsX();
      Double_t Sum_Errors = 0.;

      // Symmetrise Eta / Rapidity  ---------------------------------------

      if( title.Contains("Leptons_Eta") || title.Contains("TopQuarks_Rapidity") ) {

        SymmetricAroundZero(h_ref, h_var_up, h_var_down, h_sys, h_sys2, Sum_Errors);

      }

      else {

	for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {

	  if( h_ref->GetBinContent(bin) == 0 )  continue;

	  Sys_Error_Up   = abs(h_ref->GetBinContent(bin) - h_var_up->GetBinContent(bin))/h_ref->GetBinContent(bin);
	  Sys_Error_Down = abs(h_ref->GetBinContent(bin) - h_var_down->GetBinContent(bin))/h_ref->GetBinContent(bin);
  
	  //	  Sys_Error  = TMath::Max ( Sys_Error_Up, Sys_Error_Down);
	  Sys_Error  = (Sys_Error_Up+Sys_Error_Down)/2.;
	  Sys_Error2 = Sys_Error * Sys_Error;

	  h_sys->SetBinContent(  bin, Sys_Error  );
	  h_sys2->SetBinContent( bin, Sys_Error2 );
	  cout << "Sys. Error: " << Sys_Error << endl;
	  Sum_Errors = Sum_Errors+Sys_Error;

	}

      }

      cout << "-------------" << endl;
      cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
      cout << "-------------" << endl << endl;

      setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
      setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
      h_style = h_sys;

      title.Append("_sys_error_MATCH");
      h_sys->SetName(title);
      Hlist.Add(h_sys);

      sum_sys_errors[j]->Add(h_sys2);
      sum_mod_errors[j]->Add(h_sys2);

    }

    leg->AddEntry(h_style,"ME/PS matching","F");
    cout << "-------------" << endl << endl;

  }

  //  HADRONISATION --------------------------------------------------

  LineColor = kBlack;
  LineStyle = 1;
  LineWidth = 1;
  HistColor = kMagenta-6;

  files[21] = new TFile(inpath.Copy().Append("powheg_pythia6_Fall10.root"));
  files[22] = new TFile(inpath.Copy().Append("powheg_herwig6_Fall10.root"));

  if( files[21] && !files[21]->IsZombie() && files[22] && !files[22]->IsZombie() && isModel ) {

    cout << "   >>>   " << "HADRONISATION" << endl;
    cout << "----------------------------" << endl << endl;

    for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

      TString title = ListOfKeys->At(j)->GetName();
      cout << "   >>>   " << title << endl;
      if( !(isPSE) &&   title.Contains("PSE_")  )  continue;
      if(   isPSE  && !(title.Contains("PSE_")) )  continue;

      Sys_Error  = 0.;
      Sys_Error2 = 0.;

      files[0]->GetObject(  title, h_ref);
      files[0]->GetObject(  title, h_style);
      files[21]->GetObject( title, h_var_up);
      files[22]->GetObject( title, h_var_down);

      //      cout << "Hist.Name: " << h_var_up->GetName() << endl;
      //      cout << "Hist.Name: " << h_var_down->GetName() << endl;

      TH1* h_sys       = (TH1*) h_var_up->Clone();   h_sys->Reset();
      TH1* h_sys2      = (TH1*) h_var_down->Clone(); h_sys2->Reset();
      TH1* h_sys_up    = (TH1*) h_var_up->Clone();   h_sys_up->Reset();
      TH1* h_sys_down  = (TH1*) h_var_down->Clone(); h_sys_down->Reset();
      TH1* h_sys2_up   = (TH1*) h_var_up->Clone();   h_sys2_up->Reset();
      TH1* h_sys2_down = (TH1*) h_var_down->Clone(); h_sys2_down->Reset();

      TCanvas* Canvas = new TCanvas("plot", "plot", 800, 800);

      setHistogramStyle(h_ref,      kBlack, 1, 2, HistColor);
      setHistogramStyle(h_var_up,   kBlue,  1, 1, HistColor);
      setHistogramStyle(h_var_down, kRed,   1, 1, HistColor);

      Canvas->Clear();

      if( title.Contains("PSE_") ) {
	h_var_up->SetMaximum(1.04);
	h_var_up->SetMinimum(0.);
      }

      Canvas->Clear();
      //      h_ref->SetMarkerColor(kBlack);
      //      h_ref->Draw();
      h_var_up->SetMarkerColor(kBlue);
      h_var_up->Draw();
      h_var_down->SetMarkerColor(kRed);
      h_var_down->Draw("SAME");
      if( !(isPSE) && set_LogY ) gPad->SetLogy();
      Canvas->Print(outpath.Copy().Append(title).Append("_HAD").Append(outform));

      Int_t N_bins = h_ref->GetNbinsX();
      Double_t Sum_Errors = 0.;

      // Symmetrise Eta / Rapidity  ---------------------------------------

      if( title.Contains("Leptons_Eta") || title.Contains("TopQuarks_Rapidity") ) {

	for(Int_t bin = 1; bin <= h_ref->GetNbinsX()/2.; ++bin) {

	  if( h_var_up->GetBinContent(bin) == 0 || h_var_down->GetBinContent(bin) == 0 )  continue;

	  //	  Double_t PlusMinus_Average_Ref  = 0.;
	  Double_t PlusMinus_Average_Up   = 0.;
	  Double_t PlusMinus_Average_Down = 0.;

	  //	  PlusMinus_Average_Ref  = (h_ref->GetBinContent(bin)      + h_ref->GetBinContent(N_bins+1-bin)     ) / 2.;
	  PlusMinus_Average_Up   = (h_var_up->GetBinContent(bin)   + h_var_up->GetBinContent(N_bins+1-bin)  ) / 2.;
	  PlusMinus_Average_Down = (h_var_down->GetBinContent(bin) + h_var_down->GetBinContent(N_bins+1-bin)) / 2.;

	  Sys_Error  = abs(PlusMinus_Average_Up - PlusMinus_Average_Down) / PlusMinus_Average_Up;
	  Sys_Error2 = Sys_Error * Sys_Error;

	  h_sys->SetBinContent(           bin, Sys_Error  );
	  h_sys->SetBinContent(  N_bins+1-bin, Sys_Error  );
	  h_sys2->SetBinContent(          bin, Sys_Error2 );
	  h_sys2->SetBinContent( N_bins+1-bin, Sys_Error2 );

	  cout << "Sys. Error: " << Sys_Error << endl;
	  Sum_Errors += (2.*Sys_Error);

	}

      }

      else {

	for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {

	  if( h_var_up->GetBinContent(bin) == 0 || h_var_down->GetBinContent(bin) == 0 )  continue;

	  Sys_Error  = abs(h_var_up->GetBinContent(bin) - h_var_down->GetBinContent(bin))/h_var_up->GetBinContent(bin);
	  Sys_Error2 = Sys_Error * Sys_Error;

	  h_sys->SetBinContent(  bin, Sys_Error  );
	  h_sys2->SetBinContent( bin, Sys_Error2 );
	  cout << "Sys. Error: " << Sys_Error << endl;
	  Sum_Errors = Sum_Errors+Sys_Error;

	}

      }

      cout << "-------------" << endl;
      cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
      cout << "-------------" << endl << endl;

      setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
      setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
      h_style = h_sys;

      title.Append("_sys_error_HAD");
      h_sys->SetName(title);
      Hlist.Add(h_sys);

      sum_sys_errors[j]->Add(h_sys2);
      sum_mod_errors[j]->Add(h_sys2);

    }

    leg->AddEntry(h_style,"Hadronisation","F");
    cout << "-------------" << endl << endl;

  }

  // BR  ----------------------------------------------

  LineColor = kBlack;
  LineStyle = 1;
  LineWidth = 1;
  HistColor = kBlue+1;

  if(isTotal) {

    cout << "   >>>   " << "Branching Ratio" << endl;
    cout << "------------------------------" << endl << endl;

    for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

      TString title = ListOfKeys->At(j)->GetName();
      cout << "   >>>   " << title << endl;
      if( title.Contains("PSE_") )  continue;

      Sys_Error  = 0.;
      Sys_Error2 = 0.;

      files[0]->GetObject( title, h_ref);
      files[0]->GetObject( title, h_style);

      TH1* h_sys  = (TH1*) h_ref->Clone(); h_sys->Reset();
      TH1* h_sys2 = (TH1*) h_ref->Clone(); h_sys2->Reset();

      Int_t N_bins = h_ref->GetNbinsX();
      Double_t Sum_Errors = 0.;

      for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {
	Sys_Error  = 0.015; // TOTAL XS
	Sys_Error2 = Sys_Error * Sys_Error;
	h_sys->SetBinContent(  bin, Sys_Error  );
	h_sys2->SetBinContent( bin, Sys_Error2 );
	cout << "Sys. Error: " << Sys_Error << endl;
	Sum_Errors = Sum_Errors+Sys_Error;
      }

      cout << "-------------" << endl;
      cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
      cout << "-------------" << endl << endl;

      setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
      setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
      h_style = h_sys;

      title.Append("_sys_error_BR");
      h_sys->SetName(title);
      Hlist.Add(h_sys);
      sum_sys_errors[j]->Add(h_sys2);
    }

    leg->AddEntry(h_style,"BR (t #bar{t} #rightarrow l l X)","F");
    cout << "-------------" << endl << endl;

  }

  // LUMINOSITY  ----------------------------------------------

  LineColor = kBlack;
  LineStyle = 1;
  LineWidth = 1;
  HistColor = kAzure+5;

  if(isTotal) {

    cout << "   >>>   " << "LUMINOSITY" << endl;
    cout << "-------------------------" << endl << endl;

    for(Int_t j = 0; j < ListOfKeys->GetEntries(); ++j) {

      TString title = ListOfKeys->At(j)->GetName();
      cout << "   >>>   " << title << endl;
      if( title.Contains("PSE_") )  continue;

      Sys_Error  = 0.;
      Sys_Error2 = 0.;

      files[0]->GetObject( title, h_ref);
      files[0]->GetObject( title, h_style);

      TH1* h_sys  = (TH1*) h_ref->Clone(); h_sys->Reset();
      TH1* h_sys2 = (TH1*) h_ref->Clone(); h_sys2->Reset();

      Int_t N_bins = h_ref->GetNbinsX();
      Double_t Sum_Errors = 0.;

      for(Int_t bin = 1; bin <= h_ref->GetNbinsX(); ++bin) {
	Sys_Error  = 0.045;
	Sys_Error2 = Sys_Error * Sys_Error;
	h_sys->SetBinContent(  bin, Sys_Error  );
	h_sys2->SetBinContent( bin, Sys_Error2 );
	cout << "Sys. Error: " << Sys_Error << endl;
	Sum_Errors = Sum_Errors+Sys_Error;
      }

      cout << "-------------" << endl;
      cout << " Ave. Error: " << Sum_Errors/N_bins << endl;
      cout << "-------------" << endl << endl;

      setHistogramStyle(h_sys,  LineColor, LineStyle, LineWidth, HistColor);
      setHistogramStyle(h_sys2, LineColor, LineStyle, LineWidth, HistColor);
      h_style = h_sys;

      title.Append("_sys_error_LUMI");
      h_sys->SetName(title);
      Hlist.Add(h_sys);
      sum_sys_errors[j]->Add(h_sys2);
    }

    leg->AddEntry(h_style,"Luminosity","F");
    cout << "-------------------------" << endl << endl;

  }

}

void CompactSystematicError() {

  setHHStyle(*gStyle);
  SystematicErrors();

  Double_t hist_max = 1.;
  Double_t hist_min = 0.;

  TList *ListOfKeys = files[0]->GetListOfKeys();

  for(Int_t i = 0; i < ListOfKeys->GetEntries(); ++i) {

    //    cout << "Name: " << ListOfKeys->At(i)->GetName() << endl;
    TString title = ListOfKeys->At(i)->GetName();

    if( !(isModel) && !(isExp) && !(isStat) ) {
      cout << "No uncertainty selected!" << endl << endl;
      return;
    }
    if( (!(isModel) || isExp || isStat) && title.Contains("PSE_") )  continue;

    if( isModel && !(isExp) && !(isStat) && !(isPSE) &&   title.Contains("PSE_")  )  continue;
    if( isModel && !(isExp) && !(isStat) &&   isPSE  && !(title.Contains("PSE_")) )  continue;

    TH1 *h;
    files[0]->GetObject(title, h);
    //    TH1 *h = dynamic_cast<TH1*> (files[0]->FindObject(ListOfKeys->At(i)->GetName()));
    assert(h);

    TAxis *xaxis = h->GetXaxis();

    // Total Systematic Error

    title.Append("_total_sys_error2");
    sum_sys_errors[i]->SetName(title);

    TCanvas* Canvas = new TCanvas("plot", "plot", 800, 800);
    Canvas->Clear();
    sum_sys_errors[i]->Draw("HIST");
    leg->Draw("SAME");

    if( sum_sys_errors[i]->GetMaximum() > 0. )  hist_max = sum_sys_errors[i]->GetMaximum();
    sum_sys_errors[i]->SetMaximum(1.8 * hist_max);
    //    sum_sys_errors[i]->SetMaximum(0.01);
    sum_sys_errors[i]->SetMinimum(hist_min);
    sum_sys_errors[i]->GetXaxis()->SetTitle(xaxis->GetTitle());
    //    sum_sys_errors[i]->GetXaxis()->SetTitle(ListOfKeys->At(i)->GetName());
    sum_sys_errors[i]->GetYaxis()->SetTitle("#sum #left( #frac{#Delta #sigma}{#sigma} #right)^{2}");
    sum_sys_errors[i]->GetYaxis()->SetTitleOffset(2.4);
    sum_sys_errors[i]->GetYaxis()->SetTitleSize(0.035);

    Canvas->Print(outpath.Copy().Append(title).Append(outform));

    Hlist.Add(sum_sys_errors[i]);

    TH1* total_sys_error = SummedStackHisto(sum_sys_errors[i]);

    for(Int_t bin = 1; bin <= total_sys_error->GetNbinsX(); ++bin) {
      Double_t Total_Sys_Error = sqrt(total_sys_error->GetBinContent(bin));
      total_sys_error->SetBinContent( bin, Total_Sys_Error );
    }

    // reset title
    title = ListOfKeys->At(i)->GetName();
    title.Append("_total_sys_error");
    total_sys_error->SetName(title);

    Hlist.Add(total_sys_error);

    // Summary of Systematic and Statistical Errors

    if ( isStat && isExp && isModel && !(isPSE) ) {

      TH1* stat_error = SummedStackHisto(sum_stat_errors[i]);
      // reset title
      title = ListOfKeys->At(i)->GetName();
      title.Append("_stat_error");
      stat_error->SetName(title);

      TH1*  exp_error = SummedStackHisto( sum_exp_errors[i]);
      // reset title
      title = ListOfKeys->At(i)->GetName();
      title.Append("_exp_error");
      exp_error->SetName(title);

      TH1*  mod_error = SummedStackHisto( sum_mod_errors[i]);
      // reset title
      title = ListOfKeys->At(i)->GetName();
      title.Append("_mod_error");
      mod_error->SetName(title);

      for(Int_t bin = 1; bin <= stat_error->GetNbinsX(); ++bin) {

	Double_t Stat_Error = TMath::Sqrt(stat_error->GetBinContent(bin));
	stat_error->SetBinContent( bin, Stat_Error );

	Double_t  Exp_Error = TMath::Sqrt(exp_error->GetBinContent(bin));
	exp_error->SetBinContent( bin, Exp_Error );

	Double_t Model_Error = TMath::Sqrt(mod_error->GetBinContent(bin));
	mod_error->SetBinContent( bin, Model_Error );

      }

      TCanvas* Canvas2 = new TCanvas("plot2", "plot2", 800, 800);

      setHistogramStyle(stat_error, kBlue,  1., 3., 0.);
      stat_error->SetFillStyle(0);
      setHistogramStyle(exp_error,  kRed,   1., 3., 0.);
      exp_error->SetFillStyle(0);
      setHistogramStyle(mod_error,  kGreen, 1., 3., 0.);
      mod_error->SetFillStyle(0);
      setHistogramStyle(total_sys_error, kBlack, 1., 4., 0.);
      total_sys_error->SetFillStyle(0);

      TLegend* leg2 = new TLegend(0.70,0.58,0.95,0.87);
      leg2->Clear();
      leg2->SetFillStyle(0);
      leg2->SetBorderSize(0);

      leg2->AddEntry(     stat_error,"Statistical","L");
      leg2->AddEntry(      exp_error,"Experimental","L");
      leg2->AddEntry(      mod_error,"Model","L");
      leg2->AddEntry(total_sys_error,"Total error","L");

      if( total_sys_error->GetMaximum() > 0. )  hist_max = total_sys_error->GetMaximum();
      total_sys_error->SetMaximum(2. * hist_max);
      total_sys_error->SetMinimum(0.);
      total_sys_error->GetXaxis()->SetTitle(xaxis->GetTitle());
      total_sys_error->GetYaxis()->SetTitle("#frac{#Delta #sigma}{#sigma}");
      total_sys_error->GetYaxis()->SetTitleOffset(2.4);
      total_sys_error->GetYaxis()->SetTitleSize(0.035);

      Hlist.Add(stat_error);
      Hlist.Add(exp_error);
      Hlist.Add(mod_error);

      Canvas2->Clear();
      total_sys_error->Draw("F");
      stat_error->Draw("F,SAME");
      exp_error->Draw("F,SAME");
      mod_error->Draw("F,SAME");
      leg2->Draw("SAME");

      // reset title
      title = ListOfKeys->At(i)->GetName();
      title.Append("_all_errors");

      Canvas2->Print(outpath.Copy().Append(title).Append(outform));

    }

  }

  files[0]->Close();

  TFile outfile(outfileName,"recreate");
  Hlist.Write();
  outfile.Close();

  cout << "-------------------------" << endl << endl;

}

int main() {

  CompactSystematicError();

}
