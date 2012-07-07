//
// 1. initialize an up-to-date ROOT version, e.g. with: ini ROOT532
//    (there are issues when drawing graphs with errors in old ROOT versions)
// 2. compile like: g++ -Wall -Wextra -O2 `root-config --cflags --libs` -lRooFit -o foldedLikelihoods foldedLikelihoods.C
// 3. execute with: foldedLikelihoods [options]
//    where options could be: 'pole'/'MSbar' and 'MSTW'/'HERA' or 'all'
//

#include "Riostream.h"

#include "TAxis.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"

#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooPlot.h"
#include "RooPolyVar.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"

#include <iostream>

#include "tdrstyle.C"
#include "foldedLikelihoods_helpers.h"
#include "foldedLikelihoods_classes.h"

double alpha_s(const double Q)
{
  //
  // c.f. arXiv:0908.1135 (Bethke: The 2009 World Average of alphaS)
  //
  const unsigned Nf = 5; // number of active flavors (assumed range here: from b mass threshold up to top mass scale!)
  const double lambda2 = 0.213*0.213; // from lambda_MSbar with alpha_s(mZ) = 0.1184 and Nf = 5
  const double beta0 = (33-2*Nf)/(12*TMath::Pi());
  const double beta1 = (153-19*Nf)/(24*TMath::Power(TMath::Pi(),2));
  const double beta2 = (77139-15099*Nf+325*Nf*Nf)/(3456*TMath::Power(TMath::Pi(),3));
  const double L = TMath::Log(Q*Q/lambda2);
  const double lnL = TMath::Log(L);
  return 1./(beta0*L) - beta1*lnL/(beta0*beta0*beta0*L*L)
    + 1./(beta0*beta0*beta0*L*L*L)*(TMath::Power(beta1/beta0,2)*(lnL*lnL-lnL-1)+beta2/beta0);
  //  return alpha_s_mZ / (1+alpha_s_mZ*beta0*TMath::Log(Q*Q/mZ/mZ)); // just as a cross check: the 1-loop approximation
}

double mMSbar(const double mPole, const double alpha_s)
{
  //
  // c.f. arXiv:0805.1333 (Bernreuther: Top Quark Physics at the LHC)
  // note that the value of alpha_s has to be given in the MSbar scheme!
  //
  const double aOverPi = alpha_s / TMath::Pi();
  return mPole / (1. + 4./3*aOverPi + 8.2364*TMath::Power(aOverPi,2) + 73.638*TMath::Power(aOverPi,3));
}

double mMSbar(const double mPole)
{
  //
  // including iterative calculation of alpha_s(mMSbar)
  //
  double mMSbarPre = mPole;
  double mMSbarNew = mPole;
  do {
    const double alpha = alpha_s(mMSbarNew);
    mMSbarPre = mMSbarNew;
    mMSbarNew = mMSbar(mPole, alpha);
  } while(TMath::Abs(mMSbarNew-mMSbarPre) > 0.01); // convergence criterion
  return mMSbarNew;
}

TGraph* getRunningAlphaGraph()
{
  int n = 111;
  double mass [n];
  double alpha[n];
  double m = 80;
  for(int i=0; i<n; i++) {
    mass [i] = m;
    alpha[i] = alpha_s(m);
    m += 1;
  }
  TGraph* graph = new TGraph(n, mass, alpha);
  graph->GetXaxis()->SetTitle("Q (GeV)");
  graph->GetYaxis()->SetTitle("#alpha_{S}");
  graph->SetTitle("At three-loop level and with N(flavors)=5, using 2009 world average for #alpha_{S}(M_{Z})");
  graph->GetYaxis()->CenterTitle();
  return graph;
}

TGraph* getMassShiftGraph()
{
  int n = 81;
  double mass [n];
  double shift[n];
  double m = 140;
  for(int i=0; i<n; i++) {
    mass [i] = mMSbar(m);
    shift[i] = m-mMSbar(m);
    m += 1;
  }
  TGraph* graph = new TGraph(n, mass, shift);
  graph->GetXaxis()->SetTitle("m_{t}^{#bar{MS}} (GeV)");
  graph->GetYaxis()->SetTitle("m_{t}^{pole}-m_{t}^{#bar{MS}} (GeV)");
  return graph;
}

std::vector<TGraphAsymmErrors*> readTheory(const TString name, const bool pole, const bool useAlphaUnc, const bool heraPDF)
{
  const unsigned nPoints = (name=="ahrens" ? 51 : 91);

  TString fileName = "theories_7TeV/" + name;
  if(pole || name=="kidonakis")
    fileName += "_pole";
  else
    fileName += "_msbar";
  if(heraPDF && name!="kidonakis")
    fileName += "_hera.tab";
  else
    fileName += "_mstw.tab";

  ifstream in;
  in.open(fileName);

  double mass     [nPoints];
  double sigma    [nPoints];
  double dummy_err[nPoints];
  double err_yu[4][nPoints];
  double err_yd[4][nPoints];

  unsigned i=0;
  while( !in.eof() ) {
    double mt, xsec, err_u[3], err_d[3];
    in >> mt >> xsec >> err_u[0] >> err_d[0] >> err_u[1] >> err_d[1] >> err_u[2] >> err_d[2];
    if( in.eof() ) break;
    mass [i] = mt;
    sigma[i] = xsec;
    for(unsigned j=0; j<3; j++) {
      err_yu[j][i] = err_u[j];
      err_yd[j][i] = err_d[j];
    }
    if(useAlphaUnc) {
      err_yu[3][i] = TMath::Sqrt(err_u[0]*err_u[0] + err_u[1]*err_u[1] + err_u[2]*err_u[2]);
      err_yd[3][i] = TMath::Sqrt(err_d[0]*err_d[0] + err_d[1]*err_d[1] + err_d[2]*err_d[2]);
    }
    else {
      err_yu[3][i] = TMath::Sqrt(err_u[0]*err_u[0] + err_u[1]*err_u[1]);
      err_yd[3][i] = TMath::Sqrt(err_d[0]*err_d[0] + err_d[1]*err_d[1]);
    }
    dummy_err[i] = 0;
    if(sigma[i]<1) {
      std::cout << "Spurious cross section of " << sigma[i] << " in line " << i+1 << " of " << fileName << std::endl;
      abort();
    }
    i++;      
  }

  in.close();
  if(i!=nPoints){
    std::cout << "Spurious number of lines in " << fileName << std::endl;
    std::cout << "Found / expected: " << i << " / " << nPoints << std::endl;
    abort();
  }

  std::vector<TGraphAsymmErrors*> vec;
  for(unsigned j=0; j<4; j++)
    vec.push_back(new TGraphAsymmErrors(nPoints, mass, sigma, dummy_err, dummy_err, err_yd[j], err_yu[j]));
  return vec;
}

TString epsString(const TString& label, const bool pole, const bool heraPDF)
{
  TString name = "figs/"+label;
  name += (pole ? "_pole" : "_MSbar");
  name += (heraPDF ? "_hera" : "_mstw");
  return (name + ".eps");
}

void drawTheoryGraph(TGraphAsymmErrors* graph, TCanvas* canvas, const bool pole,
		     const TString label, const TString printNameBase)
{
  if(pole)
    graph->GetXaxis()->SetTitle("m_{t}^{pole} (GeV)");
  else
    graph->GetXaxis()->SetTitle("m_{t}^{#bar{MS}} (GeV)");
  graph->GetYaxis()->SetTitle("#sigma_{t#bar{t}} (pb)");
  graph->GetYaxis()->CenterTitle();

  graph->Draw("AP");
  graph->Draw("2");
  graph->GetFunction("f1")->Draw("same");

  canvas->Print(printNameBase+".ps");
  canvas->Print("figs/"+label+".eps");
}

void removeEmptyPoints(TGraph* graph)
{
  unsigned n = graph->GetN();
  for(unsigned i=0; i<n; i++) {
    if((graph->GetY())[i]<0.0001) {
      graph->RemovePoint(i);
      i--;
      n--;
    }
  }
}

std::vector<TF1*> getAndDrawRelativeUncertainty(const TGraphAsymmErrors* graph, TCanvas* canvas, const bool pole,
						const TString title, const TString epsName, const TString printNameBase,
						const bool symmetrize)
{
  const int n = graph->GetN();
  double uncUp  [n];
  double uncDown[n];
  double uncMean[n];
  double* y = graph->GetY();
  for(int i=0; i<n; i++) {
    uncUp  [i] = graph->GetErrorYhigh(i)/y[i];
    uncDown[i] = graph->GetErrorYlow (i)/y[i];
    uncMean[i] = (uncUp[i] + uncDown[i])/2;
  }
  TGraph* relUncUp   = new TGraph(n, graph->GetX(), uncUp  );
  TGraph* relUncDown = new TGraph(n, graph->GetX(), uncDown);
  TGraph* relUncMean = new TGraph(n, graph->GetX(), uncMean);
  removeEmptyPoints(relUncUp  );
  removeEmptyPoints(relUncDown);
  removeEmptyPoints(relUncMean);
  TString titleApp = ": " + title + " uncertainty";
  relUncUp->SetTitle(graph->GetTitle() + titleApp);
  if(pole)
    relUncUp->GetXaxis()->SetTitle("m_{t}^{pole} (GeV)");
  else
    relUncUp->GetXaxis()->SetTitle("m_{t}^{#bar{MS}} (GeV)");
  relUncUp->GetYaxis()->SetTitle("#delta#sigma_{t#bar{t}} / #sigma_{t#bar{t}}");
  relUncUp->GetYaxis()->CenterTitle();
  relUncUp->GetYaxis()->SetRangeUser(.0, .115);
  relUncUp  ->SetMarkerStyle(26);
  relUncDown->SetMarkerStyle(32);
  relUncUp->Draw("AP");
  relUncDown->Draw("P");
  std::vector<TF1*> vec;
  if(symmetrize) {
    relUncMean->Fit("pol1");
    relUncMean->GetFunction("pol1")->Draw("same");
    vec.push_back((TF1*) relUncMean->GetFunction("pol1")->Clone());
  }
  else {
    relUncUp->Fit("pol1");
    relUncUp->GetFunction("pol1")->Draw("same");
    relUncDown->Fit("pol1");
    relUncDown->GetFunction("pol1")->Draw("same");
    vec.push_back((TF1*) relUncUp  ->GetFunction("pol1")->Clone());
    vec.push_back((TF1*) relUncDown->GetFunction("pol1")->Clone());
  }
  TLegend leg = TLegend(0.8, 0.2, 0.9, 0.35);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  leg.AddEntry(relUncUp  , "#delta_{+}", "P");
  leg.AddEntry(relUncDown, "#delta_{-}", "P");
  leg.Draw();
  gStyle->SetOptFit(0000);
  if(!(epsName.Contains("kidonakis") && (!epsName.Contains("pole") || epsName.Contains("hera")))) {
    canvas->Print(printNameBase+".ps");
    canvas->Print(epsName);
  }
  delete relUncUp;
  delete relUncDown;
  delete relUncMean;
  gStyle->SetOptFit(0011);
  return vec;
}

void drawConvolution(const PredXSec& predXSec, RooRealVar& xsec, RooRealVar& mass, const TF1* xsecFunc, const TString title,
		     TCanvas* canvas, const TString printNameBase, const TString epsName)
{
  if(mass.getTitle().Contains("pole"))
    mass.setVal(172.5);
  else
    mass.setVal(165.);
  RooPlot* convFrame = xsec.frame(0.85*xsecFunc->Eval(mass.getVal()), 1.15*xsecFunc->Eval(mass.getVal()));
  predXSec.gaussianProb.plotOn(convFrame, RooFit::LineColor(kGreen),
			       RooFit::NormRange("xsec_fullRange"),
			       RooFit::FillStyle(1001), RooFit::FillColor(kGreen), RooFit::DrawOption("F"));
  predXSec.gaussianProb.plotOn(convFrame, RooFit::Range(xsecFunc->Eval(mass.getVal())-0.01,
							xsecFunc->Eval(mass.getVal())+0.01),
			       RooFit::NormRange("xsec_fullRange"),
			       RooFit::LineColor(kBlack), RooFit::LineWidth(1), RooFit::VLines());
  predXSec.rectangularProb.plotOn(convFrame, RooFit::Normalization(1.), RooFit::LineStyle(2), RooFit::FillColor(0));
  predXSec.prob.plotOn(convFrame, RooFit::LineColor(kRed), RooFit::FillColor(0));
  convFrame->GetYaxis()->SetTitle("Probability density");
  convFrame->Draw();
  gPad->RedrawAxis();
  TLegend convLeg = TLegend(0.66, 0.6, 0.92, 0.9);
  convLeg.SetFillStyle(0);
  convLeg.SetBorderSize(0);
  convLeg.AddEntry(convFrame->findObject(predXSec.name+"_gaussianProb_Norm[xsec]"   ), "Gauss(PDF#oplus#alpha_{S})"  , "F");
  //  convLeg.AddEntry(convFrame->findObject(predXSec.name+"_rectangularProb_Norm[xsec]"), "Rect(#mu_{r/f}) (#times0.59)", "F");
  convLeg.AddEntry(convFrame->findObject(predXSec.name+"_rectangularProb_Norm[xsec]"), "Rect(#mu_{r/f})", "F");
  convLeg.AddEntry(convFrame->findObject(predXSec.name+"_prob_Norm[xsec]"           ), "Gauss#otimesRect"            , "F");
  convLeg.Draw();
  char tmpTxt[99];
  sprintf(tmpTxt, "%s, %s = %.1f GeV", title.Data(), mass.getTitle().Data(), mass.getVal());
  TLatex text(0.,0.,tmpTxt);
  text.SetNDC();
  text.SetTextAlign(13);
  text.SetX(0.16);
  text.SetY(0.99);
  text.SetTextFont(43);
  text.SetTextSizePixels(16);
  text.Draw();
  canvas->Print(printNameBase+".ps");
  canvas->Print(epsName);
  delete convFrame;
}

void plot2DimProb(const PredXSec& predXSec, TCanvas* canvas, const TString title,
		  const TString printNameBase, const TString epsName)
{
  TH1* h2 = predXSec.prob.createHistogram("mass,alpha,xsec",50,30,90);
  h2->Draw("box");
  TLatex text(0.,0.,title);
  text.SetNDC();
  text.SetTextAlign(13);
  text.SetX(0.16);
  text.SetY(0.99);
  text.SetTextFont(43);
  text.SetTextSizePixels(16);
  text.Draw();
  canvas->Print(printNameBase+".ps");
  canvas->Print(epsName);
  delete h2;
}

void plotProjectedPDF(RooAbsPdf* pdf, RooPlot* frame, const int color, const int fillStyle,
		      const double max, const double errLow, const double errHig)
{
  pdf->plotOn(frame, RooFit::LineColor(color), RooFit::NormRange("mass_fullRange"));
  pdf->plotOn(frame, RooFit::Range(max-errLow, max+errHig), RooFit::NormRange("mass_fullRange"),
	      //	      RooFit::LineColor(color), RooFit::LineWidth(1),
	      RooFit::FillStyle(fillStyle), RooFit::FillColor(color), RooFit::DrawOption("F"), RooFit::VLines());
  pdf->plotOn(frame, RooFit::Range(max-.01, max+.01), RooFit::NormRange("mass_fullRange"),
	      RooFit::LineColor(color), RooFit::LineWidth(3), RooFit::VLines());
}

TLatex* cmsTxt(const bool full2011dilep)
{
  const TString txt = (full2011dilep ?
		       "CMS, 2.3 fb^{-1} at  #sqrt{s} = 7 TeV" :
		       "CMS Preliminary, 1.14 fb^{-1} at  #sqrt{s} = 7 TeV");
  TLatex* text = new TLatex(3.570061,23.08044,txt);
  text->SetNDC();
  text->SetTextAlign(13);
  text->SetX(0.16);
  text->SetY(1.003);
  text->SetTextFont(42);
  text->SetTextSizePixels(24);
  return text;
}

int foldedLikelihoods(const bool pole, const bool heraPDF)
{
  const bool useAlphaUnc = false;
  const bool full2011dilep = true;

  setTDRStyle();
  gStyle->SetTitleBorderSize(1);
  gStyle->SetOptFit(0011);

  RooAbsReal::defaultIntegratorConfig()->method2D().setLabel("RooSegmentedIntegrator2D");
  RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setCatLabel("sumRule", "Midpoint");
  //  RooAbsReal::defaultIntegratorConfig()->methodND().setLabel("RooBinIntegrator");
  RooAbsReal::defaultIntegratorConfig()->methodND().setLabel("RooMCIntegrator");
  //  RooAbsReal::defaultIntegratorConfig()->printMultiline(std::cout, 1, true);

  gSystem->mkdir("figs");

  TString printNameBase = "figs/foldedLikelihoods";
  printNameBase += (pole ? "_pole" : "_MSbar");
  printNameBase += (heraPDF ? "_hera" : "_mstw");

  TCanvas* canvas = new TCanvas("canvas", "canvas", 10, 10, 900, 600);
  canvas->cd()->SetRightMargin(0.04);
  canvas->Print(printNameBase+".ps[");

  TGraph* runningAlpha = 0;
  if(!pole) {
    runningAlpha = getRunningAlphaGraph();
    gStyle->SetOptTitle(1);
    runningAlpha->Draw("AP");
    canvas->Print(printNameBase+".ps");
    canvas->Print("figs/alpha_s.eps");
    gStyle->SetOptTitle(0);
  }

  TGraph* mShift = 0;
  if(!pole) {
    mShift = getMassShiftGraph();
    mShift->Draw("AP");
    mShift->Fit("pol1");
    canvas->Print(printNameBase+".ps");
    canvas->Print("figs/MSbar_vs_pole_mass.eps");
  }

  std::vector<TGraphAsymmErrors*> ahr_vec = readTheory("ahrens"   , pole, useAlphaUnc, heraPDF);
  std::vector<TGraphAsymmErrors*> kid_vec = readTheory("kidonakis", pole, useAlphaUnc, heraPDF);
  std::vector<TGraphAsymmErrors*> moc_vec = readTheory("moch"     , pole, useAlphaUnc, heraPDF);

  TGraphAsymmErrors* ahr = ahr_vec.at(3);
  TGraphAsymmErrors* kid = kid_vec.at(3);
  TGraphAsymmErrors* moc = moc_vec.at(3);

  TF1* f1 = new TF1("f1", "([0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3))/TMath::Power(x,4)", 0, 1000);
  kid->Fit(f1, "0");
  ahr->Fit(f1, "0");
  moc->Fit(f1, "0");

  const TString pdfName = (heraPDF ? "(with HERAPDF)" : "(with MSTW)");
  const TString title[3] = {"Ahrens et al. "     + pdfName,
			    "Kidonakis "         + pdfName,
			    "Langenfeld et al. " + pdfName};
  for(unsigned j=0; j<4; j++) {
    ahr_vec.at(j)->SetTitle(title[0]);
    kid_vec.at(j)->SetTitle(title[1]);
    moc_vec.at(j)->SetTitle(title[2]);
  }
    
  gStyle->SetOptTitle(1);

  const TString errName [4] = {"Scale", "Experimental PDF", "#alpha_{S}", "Total"};
  const TString errLabel[4] = {"scale", "expPDF", "alphaS", "total"};

  std::vector<TF1*> kid_funcs[4];
  std::vector<TF1*> moc_funcs[4];
  std::vector<TF1*> ahr_funcs[4];

  for(unsigned j=0; j<4; j++)
    kid_funcs[j] = getAndDrawRelativeUncertainty(kid_vec.at(j), canvas, pole, errName[j],
						 epsString("relUnc_kidonakis_"+errLabel[j], pole, heraPDF),
						 printNameBase, (j==1 || j==2));
  if(pole && !heraPDF)
    drawTheoryGraph(kid, canvas, pole, "kidonakis", printNameBase);
  for(unsigned j=0; j<4; j++)
    moc_funcs[j] = getAndDrawRelativeUncertainty(moc_vec.at(j), canvas, pole, errName[j],
						 epsString("relUnc_moch_"+errLabel[j], pole, heraPDF),
						 printNameBase, (j==1 || j==2));
  drawTheoryGraph(moc, canvas, pole, "moch"  , printNameBase);
  for(unsigned j=0; j<4; j++)
    ahr_funcs[j] = getAndDrawRelativeUncertainty(ahr_vec.at(j), canvas, pole, errName[j],
						 epsString("relUnc_ahrens_"+errLabel[j], pole, heraPDF),
						 printNameBase, (j==1 || j==2));
  drawTheoryGraph(ahr, canvas, pole, "ahrens", printNameBase);

  gStyle->SetOptTitle(0);

  RooRealVar mass("mass", "m_{t}^{pole}", 140., 190., "GeV");
  if(!pole)
    mass.SetTitle("m_{t}^{#bar{MS}}");

  RooRealVar alpha("alpha", "#alpha_{S}(M_{Z})", 0.10, 0.13);
  alpha.setVal(0.117);
  //  alpha.setConstant(kTRUE);

  RooRealVar xsec("xsec" , "#sigma_{t #bar{t}}", 0., 900., "pb");

  mass.setRange("mass_fullRange", 140., 190.);
  xsec.setRange("xsec_fullRange", 10., 900.);

  RooRealVar shift_p0("shift_p0", "shift_p0", 0.);
  RooRealVar shift_p1("shift_p1", "shift_p1", 0.);
  if(!pole) {
    shift_p0.setVal(mShift->GetFunction("pol1")->GetParameter(0));
    shift_p1.setVal(mShift->GetFunction("pol1")->GetParameter(1));
  }
  RooPolyVar shift("shift", "shift", mass, RooArgSet(shift_p0, shift_p1));

  TString deltaM_formula = "@0-172.5+@1";
  if(full2011dilep)
    deltaM_formula = "@0+@1";
  RooFormulaVar deltaM("deltaM", "deltaM", deltaM_formula, RooArgSet(mass, shift));
  // e-mu
  //  RooRealVar p0("p0", "p0",  1.01    ); // e-mu
  //  RooRealVar p1("p1", "p1", -0.93E-02); // e-mu
  //  RooRealVar p2("p2", "p2",  1.2E-04 ); // e-mu
  // mu-mu
  //  RooRealVar p0("p0", "p0",  1.00    ); // mu-mu
  //  RooRealVar p1("p1", "p1", -1.01E-02); // mu-mu
  //  RooRealVar p2("p2", "p2",  1.1E-04 ); // mu-mu
  // e-e
  //  RooRealVar p0("p0", "p0",  1.02    ); // e-e
  //  RooRealVar p1("p1", "p1", -1.02E-02); // e-e
  //  RooRealVar p2("p2", "p2",  1.0E-04 ); // e-e
  // weighted average
  RooRealVar p0("p0", "p0",  1.01    ); // e-mu
  RooRealVar p1("p1", "p1", -0.97E-02); // e-mu
  RooRealVar p2("p2", "p2",  1.125E-04 ); // e-mu

  if(full2011dilep) {
    p0.setVal(-1.467);
    p1.setVal(0.03715);
    p2.setVal(-0.0001324);
  }

  RooPolyVar measXSecMassDepRel("measXSecMassDepRel", "measXSecMassDepRel", deltaM, RooArgSet(p0, p1, p2));
  RooRealVar measXSec("measXSec", "measXSec", 169.9, "pb");
  // total exp uncert of 18.23 from adding quadr. 3.9 (stat), 16.3 (syst) and 7.6 (lumi)
  // and subtracting 1.5 % = 2.5485 pb for the top mass (in e-mu, it is 2.6 % for ee and mumu)
  RooRealVar measXSecErr("measXSecErr", "measXSecErr", 18.23, "pb");

  if(full2011dilep) {
    measXSec.setVal(161.9);
    measXSecErr.setVal(6.725);
  }

  RooFormulaVar measXSecMassDep("measXSecMassDep", "measXSecMassDep", "@0*@1", RooArgSet(measXSecMassDepRel,measXSec));
  RooGaussian measXSecPDF("measXSecPDF", "measXSecPDF", xsec, measXSecMassDep, measXSecErr);

  TString alpha_funcFileName = "theories_7TeV/xSecVsMass_";
  if(pole)
    alpha_funcFileName += "pole_";
  else
    alpha_funcFileName += "msbar_";
  if(heraPDF)
    alpha_funcFileName += "hera_alphaScan.root";
  else
    alpha_funcFileName += "mstw_alphaScan.root";
  PredXSec kidPredXSec("kidPredXSec", xsec, mass, alpha, kid->GetFunction("f1"), kid_funcs, alpha_funcFileName);
  PredXSec mocPredXSec("mocPredXSec", xsec, mass, alpha, moc->GetFunction("f1"), moc_funcs, alpha_funcFileName);
  PredXSec ahrPredXSec("ahrPredXSec", xsec, mass, alpha, ahr->GetFunction("f1"), ahr_funcs, alpha_funcFileName);

  if(pole && !heraPDF) {
    drawConvolution(kidPredXSec, xsec, mass, kid->GetFunction("f1"), title[1], canvas, printNameBase,
		    epsString("convolution_kidonakis", pole, heraPDF));
    plot2DimProb(kidPredXSec, canvas, title[1], printNameBase, epsString("2d_kidonakis", pole, heraPDF));
  }
  drawConvolution(mocPredXSec, xsec, mass, moc->GetFunction("f1"), title[2], canvas, printNameBase,
		  epsString("convolution_moch", pole, heraPDF));
  plot2DimProb(mocPredXSec, canvas, title[2], printNameBase, epsString("2d_moch", pole, heraPDF));
  drawConvolution(ahrPredXSec, xsec, mass, ahr->GetFunction("f1"), title[0], canvas, printNameBase,
		  epsString("convolution_ahrens", pole, heraPDF));
  plot2DimProb(ahrPredXSec, canvas, title[0], printNameBase, epsString("2d_ahrens", pole, heraPDF));

  const int colorKid = kRed+1;
  const int colorAhr = kMagenta;
  const int colorMoc = kGreen-3;
  const int colorDil = kBlue;

  RooPlot* frame = mass.frame();
  if(pole && !heraPDF)
    kidPredXSec.xsec.plotOn(frame, RooFit::LineColor(colorKid));
  ahrPredXSec.xsec.plotOn(frame, RooFit::LineColor(colorAhr));
  mocPredXSec.xsec.plotOn(frame, RooFit::LineColor(colorMoc));
  measXSecMassDep.plotOn(frame, RooFit::LineColor(colorDil));
  frame->GetYaxis()->SetTitle("#sigma_{t#bar{t}} (pb)");
  frame->Draw();
  canvas->Print(printNameBase+".ps");
  canvas->Print(epsString("xsec_vs_mass", pole, heraPDF));

  RooPlot* frame_alpha = alpha.frame();
  mass.setVal(172.5);
  mocPredXSec.xsec.plotOn(frame_alpha, RooFit::LineColor(colorMoc));
  measXSecMassDep.plotOn(frame_alpha, RooFit::LineColor(colorDil));
  frame_alpha->GetYaxis()->SetTitle("#sigma_{t#bar{t}} (pb)");
  frame_alpha->Draw();
  canvas->Print(printNameBase+".ps");
  canvas->Print(epsString("xsec_vs_alpha", pole, heraPDF));

  RooRealVar alpha2009_mean("alpha2009_mean", "alpha2009_mean", 0.1184);
  RooRealVar alpha2009_unc ("alpha2009_unc" , "alpha2009_unc" , 0.0007);

  RooRealVar alphaMSTW_mean("alphaMSTW_mean", "alphaMSTW_mean", 0.1171); //Eur.Phys.J.C64:653-680,2009
  RooRealVar alphaMSTW_unc ("alphaMSTW_unc" , "alphaMSTW_unc" , 0.0014); //Eur.Phys.J.C64:653-680,2009

  RooRealVar alphaHERA_mean("alphaHERA_mean", "alphaHERA_mean", 0.1176);
  RooRealVar alphaHERA_unc ("alphaHERA_unc" , "alphaHERA_unc" , 0.0020);

  RooGaussian alphaMSTW_prob("alphaMSTW_prob", "alphaMSTW_prob", alpha, alphaMSTW_mean, alphaMSTW_unc);
  RooGaussian alphaHERA_prob("alphaHERA_prob", "alphaHERA_prob", alpha, alphaHERA_mean, alphaHERA_unc);

  RooRealVar mass_mean("mass_mean", "mass_mean", 173.2); // Tevatron July 2012
  RooRealVar mass_unc ("mass_unc" , "mass_unc" , 1.4); // Tevatron July 2012 plus 1 GeV for m(pole)/m(MC)?
  RooGaussian massProb("massProb", "massProb", mass, mass_mean, mass_unc);

  RooProdPdf kidProdPDF("kidProdPDF", "kidProdPDF", RooArgList(measXSecPDF, kidPredXSec.prob));
  RooProdPdf mocProdPDF("mocProdPDF", "mocProdPDF", RooArgList(measXSecPDF, mocPredXSec.prob, massProb));
  RooProdPdf ahrProdPDF("ahrProdPDF", "ahrProdPDF", RooArgList(measXSecPDF, ahrPredXSec.prob));

  RooAbsPdf* kidProjectedPDF = kidProdPDF.createProjection(xsec);
  RooAbsPdf* mocProjectedPDF = mocProdPDF.createProjection(xsec);
  RooAbsPdf* ahrProjectedPDF = ahrProdPDF.createProjection(xsec);

  FinalLikeliResults mocResultAlpha("mocResultAlpha",
				    xsec, mass, alpha, RooArgList(measXSecPDF, mocPredXSec.prob, massProb));

  FinalLikeliResults mocResultMass("mocResultMass",
				   xsec, mass, alpha, RooArgList(measXSecPDF, mocPredXSec.prob, alphaMSTW_prob));

  TF1* kidTF = 0;//kidProjectedPDF->asTF(RooArgList(mass, alpha), RooArgSet(), RooArgSet(mass, alpha));
  TF2* mocTF = (TF2*) mocProjectedPDF->asTF(RooArgList(alpha, mass), RooArgList(), RooArgSet(alpha, mass));
  TF1* ahrTF = 0;//ahrProjectedPDF->asTF(RooArgList(mass, alpha), RooArgSet(), RooArgSet(mass, alpha));

//  frame = mass.frame(RooFit::Range(165., 180.));
//  mocProjectedPDF->plotOn(frame);
//  frame->Draw();
//  mocTF_mass->Draw("L same");
//  canvas->Print(printNameBase+".ps");
//
//  frame = alpha.frame(RooFit::Range(.11, .12));
//  mocProjectedPDF->plotOn(frame);
//  frame->Draw();
//  mocTF_mass->Draw("L same");
//  canvas->Print(printNameBase+".ps");
//
//  gStyle->SetPalette(1);
//  mocTF->SetLineColor(kGray);
//  mocTF->Draw("surf7");
//  mocTF->Draw("cont1 same");
//  mocTF->GetZaxis()->SetRangeUser(0.1, 100.);
//  canvas->Print(printNameBase+".ps");
//
//  mocTF->GetYaxis()->SetRangeUser(170., 180.);
//  mocTF->GetXaxis()->SetRangeUser(0.114, 0.12);
//  mocTF->Draw("cont1");

  gStyle->SetEndErrorSize(5*gStyle->GetEndErrorSize());

  const double alphaAxisMin = 0.11;
  const double alphaAxisMax = 0.13;

  const double massAxisMin = 170.;
  const double massAxisMax = 180.;

  mocResultAlpha.point.GetXaxis()->SetLimits(alphaAxisMin, alphaAxisMax);
  mocResultAlpha.point.GetYaxis()->SetRangeUser(massAxisMin , massAxisMax);
  mocResultAlpha.point.SetPointError(0, mocResultAlpha.point.GetErrorX(0), 0);

  mocResultMass.point.GetXaxis()->SetLimits(alphaAxisMin, alphaAxisMax);
  mocResultMass.point.GetYaxis()->SetRangeUser(massAxisMin, massAxisMax);
  mocResultMass.point.SetPointError(0, 0, mocResultMass.point.GetErrorY(0));

  TBox boxAlpha2009(alpha2009_mean.getVal()-alpha2009_unc.getVal(), massAxisMin,
		    alpha2009_mean.getVal()+alpha2009_unc.getVal(), massAxisMax);
  boxAlpha2009.SetFillColor(kYellow-9);
  TLine lineAlpha2009(alpha2009_mean.getVal(), massAxisMin,
		      alpha2009_mean.getVal(), massAxisMax);
  lineAlpha2009.SetLineColor(kGray);

  TBox boxAlphaMSTW(alphaMSTW_mean.getVal()-alphaMSTW_unc.getVal(), massAxisMin,
		    alphaMSTW_mean.getVal()+alphaMSTW_unc.getVal(), massAxisMax);
  boxAlphaMSTW.SetFillColor(kGreen-9);
  TLine lineAlphaMSTW(alphaMSTW_mean.getVal(), massAxisMin,
		      alphaMSTW_mean.getVal(), massAxisMax);
  lineAlphaMSTW.SetLineColor(kGray);

  TBox boxAlphaHERA(alphaHERA_mean.getVal()-alphaHERA_unc.getVal(), massAxisMin,
		    alphaHERA_mean.getVal()+alphaHERA_unc.getVal(), massAxisMax);
  boxAlphaHERA.SetFillColor(kCyan-9);
  TLine lineAlphaHERA(alphaHERA_mean.getVal(), massAxisMin,
		      alphaHERA_mean.getVal(), massAxisMax);
  lineAlphaHERA.SetLineColor(kGray);

  TBox boxMass(alphaAxisMin, mass_mean.getVal()-mass_unc.getVal(),
	       alphaAxisMax, mass_mean.getVal()+mass_unc.getVal());
  boxMass.SetFillColor(kMagenta-9);
  TLine lineMass(alphaAxisMin, mass_mean.getVal(),
		 alphaAxisMax, mass_mean.getVal());
  lineMass.SetLineColor(kGray);

  mocResultAlpha.point.Draw("AP");
  boxMass.Draw();
  lineMass.Draw();
  if(heraPDF) {
    boxAlphaHERA.Draw();
    lineAlphaHERA.Draw();
  }
  else {
    boxAlphaMSTW.Draw();
    lineAlphaMSTW.Draw();
  }
  boxAlpha2009.Draw();
  lineAlpha2009.Draw();
  gPad->RedrawAxis();
  mocResultAlpha.point.Draw("P");
  mocResultAlpha.ellipse.Draw("L");
  canvas->Print(printNameBase+".ps");

  mocResultMass.point.Draw("AP");
  boxMass.Draw();
  lineMass.Draw();
  if(heraPDF) {
    boxAlphaHERA.Draw();
    lineAlphaHERA.Draw();
  }
  else {
    boxAlphaMSTW.Draw();
    lineAlphaMSTW.Draw();
  }
  boxAlpha2009.Draw();
  lineAlpha2009.Draw();
  gPad->RedrawAxis();
  mocResultMass.point.Draw("P");
  mocResultMass.ellipse.Draw("L");
  canvas->Print(printNameBase+".ps");

  mocResultAlpha.point.Draw("AP");
  boxMass.Draw();
  lineMass.Draw();
  if(heraPDF) {
    boxAlphaHERA.Draw();
    lineAlphaHERA.Draw();
  }
  else {
    boxAlphaMSTW.Draw();
    lineAlphaMSTW.Draw();
  }
  gPad->RedrawAxis();
  mocResultAlpha.point.Draw("P");
  mocResultAlpha.ellipse.Draw("L");
  mocResultMass.point.Draw("P");
  mocResultMass.ellipse.Draw("L");
  canvas->Print(printNameBase+".ps");

  //  delete mocTF;

  canvas->Print(printNameBase+".ps]");
  return 0;

  const double kidMax = kidTF->GetMaximumX();
  const double mocMax = mocTF->GetMaximumX();
  const double ahrMax = ahrTF->GetMaximumX();

  if(pole) {
    std::cout << "============================================================" << std::endl;
    std::cout << "Langenfeld et al.: Extracted pole mass = " << mocMax
	      << " GeV (converted to MSbar mass: " << mMSbar(mocMax) << " GeV)" << std::endl;
    if(!heraPDF)
      std::cout << "Kidonakis        : Extracted pole mass = " << kidMax
		<< " GeV (converted to MSbar mass: " << mMSbar(kidMax) << " GeV)" << std::endl;
    std::cout << "Ahrens et al.    : Extracted pole mass = " << ahrMax
	      << " GeV (converted to MSbar mass: " << mMSbar(ahrMax) << " GeV)" << std::endl;
    std::cout << "============================================================" << std::endl;
  }

  double kidLowErr, kidHigErr;
  double mocLowErr, mocHigErr;
  double ahrLowErr, ahrHigErr;

  if(pole && !heraPDF)
    getUncertaintiesFromIntegral(kidTF, kidLowErr, kidHigErr, 140., 200., 1., 0.01);
  getUncertaintiesFromIntegral(ahrTF, ahrLowErr, ahrHigErr  , 140., 200., 1., 0.01);
  getUncertaintiesFromIntegral(mocTF, mocLowErr, mocHigErr  , 140., 200., 1., 0.01);

  if(full2011dilep)
    frame = mass.frame(RooFit::Range(160., 190.));
  else
    frame = mass.frame(RooFit::Range(150., 190.));

  plotProjectedPDF(mocProjectedPDF, frame, colorMoc, 3505, mocMax, mocLowErr, mocHigErr);
  if(pole && !heraPDF)
    plotProjectedPDF(kidProjectedPDF, frame, colorKid, 3545, kidMax, kidLowErr, kidHigErr);
  plotProjectedPDF(ahrProjectedPDF, frame, colorAhr, 3554, ahrMax, ahrLowErr, ahrHigErr);

  frame->GetYaxis()->SetTitle("Probability density");
  if(full2011dilep)
    frame->GetYaxis()->SetRangeUser(0., 0.056);
  else
    frame->GetYaxis()->SetRangeUser(0., 0.034);
  frame->Draw();

  TLatex* cmsLabel = cmsTxt(full2011dilep);
  cmsLabel->Draw();

  kidTF->SetLineColor(kBlack);
  kidTF->SetFillColor(colorKid);
  kidTF->SetFillStyle(3545);

  mocTF->SetLineColor(kBlack);
  mocTF->SetFillColor(colorMoc);
  mocTF->SetFillStyle(3305);

  ahrTF->SetLineColor(kBlack);
  ahrTF->SetFillColor(colorAhr);
  ahrTF->SetFillStyle(3554);

  double yMin = 0.70;
  if(!pole || heraPDF)
    yMin = 0.77;
  TLegend leg = TLegend(0.2, yMin, 0.9, 0.9);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);

  char legTmpTxt[99];
  sprintf(legTmpTxt, "Langenfeld et al. (%.1f +%.1f -%.1f GeV)", mocMax, mocHigErr, mocLowErr);
  leg.AddEntry(mocTF, legTmpTxt, "F");
  if(pole && !heraPDF) {
    sprintf(legTmpTxt, "Kidonakis (%.1f +%.1f -%.1f GeV)", kidMax, kidHigErr, kidLowErr);
    leg.AddEntry(kidTF, legTmpTxt, "F");
  }
  sprintf(legTmpTxt, "Ahrens et al. (%.1f +%.1f -%.1f GeV)", ahrMax, ahrHigErr, ahrLowErr);
  leg.AddEntry(ahrTF, legTmpTxt, "F");
  leg.Draw();

  canvas->Print(printNameBase+".ps");
  canvas->Print(epsString("densities", pole, heraPDF));

  // producing final overview plots

  //  gStyle->SetOptFit(0000); ???
  ahr->GetXaxis()->SetRangeUser(140., 190.);
  ahr->SetFillColor(colorAhr);
  if(pole && !heraPDF)
    kid->SetFillColor(colorKid);
  moc->SetFillColor(colorMoc);
  //  moc->SetFillStyle(3002);
  ahr->Draw("A3");
  if(pole && !heraPDF)
    kid->Draw("3 same");
  moc->Draw("3 same");
  moc->Draw("L same");
  ahr->GetYaxis()->CenterTitle(false);
  gPad->RedrawAxis();
  canvas->Print(printNameBase+".ps");
  
  // cleaning up

  canvas->Print(printNameBase+".ps]");

  if(!pole) {
    delete runningAlpha;
    delete mShift;
  }
  for(unsigned j=0; j<4; j++) {
    delete ahr_vec.at(j);
    delete kid_vec.at(j);
    delete moc_vec.at(j);
  }
  delete canvas;
  delete f1;
  delete cmsLabel;

  return 0;
}

//////////////////////////////////////////////////
// main function
//////////////////////////////////////////////////

int main(const int argc, const char** argv) {
  TString arguments = "";
  for(int i=1; i<argc; i++)
    arguments += argv[i];
  arguments.ToLower();
  if(arguments=="all") {
    return
      foldedLikelihoods(0, 0) +
      foldedLikelihoods(0, 1) +
      foldedLikelihoods(1, 0) +
      foldedLikelihoods(1, 1);
  }
  else if(arguments.Contains("all")) {
    std::cout << "Option 'all' should not be combined with other arguments!" << std::endl;
    return 999;
  }
  if(arguments.Contains("pole") && arguments.Contains("msbar")) {
    std::cout << "Should choose either 'pole' or 'MSbar'...!" << std::endl;
    return 999;
  }
  if(arguments.Contains("hera") && arguments.Contains("mstw")) {
    std::cout << "Should choose either 'MSTW' or 'HERA'...!" << std::endl;
    return 999;
  }
  const bool pole    = !arguments.Contains("msbar");
  const bool heraPDF =  arguments.Contains("hera");
  return foldedLikelihoods(pole, heraPDF);
}
