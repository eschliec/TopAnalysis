#ifndef foldedLikelihoods_classes_h
#define foldedLikelihoods_classes_h

#include "TFile.h"
#include "TMatrixD.h"
#include "TVirtualFitter.h"

#include "RooEllipse.h"
#include "RooNumIntConfig.h"

#include "foldedLikelihoods_helpers.h"

//////////////////////////////////////////////////
/// Class for linear dependencies
//////////////////////////////////////////////////

class LinearDependence {

public:

  LinearDependence(const TString& label, const TF1* fittedFunction, RooRealVar& x);

private:

  const RooRealVar p0;
  const RooRealVar p1;

public:

  const RooPolyVar polyVar;

};

LinearDependence::LinearDependence(const TString& label, const TF1* fittedFunction, RooRealVar& x):
  p0(label+"_p0", label+"_p0", fittedFunction->GetParameter(0)),
  p1(label+"_p1", label+"_p1", fittedFunction->GetParameter(1)),
  polyVar(label, label, x, RooArgSet(p0, p1))
{}

//////////////////////////////////////////////////
/// Class for quadratic dependencies
//////////////////////////////////////////////////

class QuadraticDependence {

public:

  QuadraticDependence(const TString& label, const TF1* fittedFunction, RooRealVar& x);

private:

  const RooRealVar p0;
  const RooRealVar p1;
  const RooRealVar p2;

public:

  const RooPolyVar polyVar;

};

QuadraticDependence::QuadraticDependence(const TString& label, const TF1* fittedFunction, RooRealVar& x):
  p0(label+"_p0", label+"_p0", fittedFunction->GetParameter(0)),
  p1(label+"_p1", label+"_p1", fittedFunction->GetParameter(1)),
  p2(label+"_p2", label+"_p2", fittedFunction->GetParameter(2)),
  polyVar(label, label, x, RooArgSet(p0, p1, p2))
{}

//////////////////////////////////////////////////
/// Class for predicted cross sections with
/// mass dependence, alphaS dependence and uncertainties
//////////////////////////////////////////////////

class PredXSec {

public:

  PredXSec(const TString& label, RooRealVar& xsec_var, RooRealVar& mass_var, RooRealVar& alpha_var,
	   const TF1* xsec_func, const std::vector<TF1*>* unc_funcs, const TString& alpha_funcFileName);

public:

  const TString name;

private:

  const RooRealVar p0;
  const RooRealVar p1;
  const RooRealVar p2;
  const RooRealVar p3;

  const LinearDependence relUncPdf;
  const LinearDependence relUncScaleUp;
  const LinearDependence relUncScaleDown;

  TFile alpha_funcFile;

  const QuadraticDependence alpha_p0;
  const QuadraticDependence alpha_p1;
  const QuadraticDependence alpha_p2;

  const RooPolyVar alphaDep;

public:

  RooFormulaVar xsec;

private:

  RooFormulaVar xsecScaleUp;
  RooFormulaVar xsecScaleDown;
  RooFormulaVar gaussianUnc;

  RooNumIntConfig mcIntegratorCfg;

public:

  RooGaussian gaussianProb;
  RooGenericPdf rectangularProb;
  RooGenericPdf prob;

};

PredXSec::PredXSec(const TString& label, RooRealVar& xsec_var, RooRealVar& mass_var, RooRealVar& alpha_var,
		   const TF1* xsec_func, const std::vector<TF1*>* unc_funcs, const TString& alpha_funcFileName):
  name(label),
  p0(label+"_p0", label+"_p0", xsec_func->GetParameter(0)),
  p1(label+"_p1", label+"_p1", xsec_func->GetParameter(1)),
  p2(label+"_p2", label+"_p2", xsec_func->GetParameter(2)),
  p3(label+"_p3", label+"_p3", xsec_func->GetParameter(3)),
  relUncPdf      (label+"_relUncPdf"      , unc_funcs[1].at(0), mass_var),
  relUncScaleUp  (label+"_relUncScaleUp"  , unc_funcs[0].at(0), mass_var),
  relUncScaleDown(label+"_relUncScaleDown", unc_funcs[0].at(1), mass_var),
  alpha_funcFile(alpha_funcFileName, "READ"),
  alpha_p0(label+"_alpha_p0", ((TGraph*)alpha_funcFile.Get("graph_p0"))->GetFunction("pol2"), mass_var),
  alpha_p1(label+"_alpha_p1", ((TGraph*)alpha_funcFile.Get("graph_p1"))->GetFunction("pol2"), mass_var),
  alpha_p2(label+"_alpha_p2", ((TGraph*)alpha_funcFile.Get("graph_p2"))->GetFunction("pol2"), mass_var),
  alphaDep(label+"_alphaDep", label+"_alphaDep", alpha_var, RooArgSet(alpha_p0.polyVar,
								      alpha_p1.polyVar,
								      alpha_p2.polyVar)),
  xsec(label+"_xsec", label+"_xsec", "(@1+@2*@0+@3*@0*@0+@4*@0*@0*@0)*@5/(@0*@0*@0*@0)",
       RooArgSet(mass_var, p0, p1, p2, p3, alphaDep)),
  xsecScaleUp(label+"_xsecScaleUp", label+"_xsecScaleUp", "@0+@0*TMath::Abs(@1)",
	      RooArgSet(xsec, relUncScaleUp.polyVar)),
  xsecScaleDown(label+"_xsecScaleDown", label+"_xsecScaleDown", "@0-@0*TMath::Abs(@1)",
		RooArgSet(xsec, relUncScaleDown.polyVar)),
  gaussianUnc(label+"_gaussianUnc", label+"_gaussianUnc", "@0*TMath::Abs(@1)",
	      RooArgSet(xsec, relUncPdf.polyVar)),
  mcIntegratorCfg(*RooAbsReal::defaultIntegratorConfig()),
  gaussianProb(label+"_gaussianProb", label+"_gaussianProb", xsec_var, xsec, gaussianUnc),
  rectangularProb(label+"_rectangularProb", label+"_rectangularProb",
		  "(@0 >= @1) && (@0 < @2)", RooArgList(xsec_var, xsecScaleDown, xsecScaleUp)),
  prob(label+"_prob", label+"_prob",
       "1/(2*(@3-@2))*(TMath::Erf((@3-@0)/(@1*TMath::Sqrt(2)))-TMath::Erf((@2-@0)/(@1*TMath::Sqrt(2))))",
       RooArgList(xsec_var, gaussianUnc, xsecScaleDown, xsecScaleUp))
{
  alpha_funcFile.Close();
  mcIntegratorCfg.method1D().setLabel("RooMCIntegrator");
  rectangularProb.setIntegratorConfig(mcIntegratorCfg);
}

//////////////////////////////////////////////////
/// Class for vertical reflections of TF1 functions
//////////////////////////////////////////////////

class VerticalReflector {

public:
  
  VerticalReflector(TF1* func): f(func) {}
  double operator() (double *x, double *p) const { return -f->EvalPar(x, p); };

private:

  TF1* f;

};

//////////////////////////////////////////////////
/// Final likelihood and extracted results
/// for m_top and alpha_S
//////////////////////////////////////////////////

class FinalLikeliResults {

public:
  
  FinalLikeliResults(const TString& label, RooRealVar& xsec_var, RooRealVar& mass_var, RooRealVar& alpha_var,
		     const RooArgList& prodPdfList);

private:

  RooProdPdf prodPdf;

public:

  RooAbsPdf* projPdf;
  TF2* f2;
  TF1* f1_alpha;
  TF1* f1_mass;
  TGraphErrors point;
  TMatrixD covM;
  RooEllipse ellipse;

};

FinalLikeliResults::FinalLikeliResults(const TString& label,
				       RooRealVar& xsec_var, RooRealVar& mass_var, RooRealVar& alpha_var,
				       const RooArgList& prodPdfList):
  prodPdf(label+"_prodPdf", label+"_prodPdf", prodPdfList)
{
  projPdf = prodPdf.createProjection(xsec_var);
  f2 = (TF2*) projPdf->asTF(RooArgList(alpha_var, mass_var), RooArgList(), RooArgSet(alpha_var, mass_var));
  VerticalReflector reflector(f2);
  TF2 f2Reflected((TString)f2->GetName()+"_reflected", reflector,
		  f2->GetXmin(),
		  f2->GetXmax(),
		  f2->GetYmin(),
		  f2->GetYmax(),
		  f2->GetNpar(), "VerticalReflector");
  double bestX, bestY;
  f2Reflected.GetMinimumXY(bestX, bestY);
  TVirtualFitter::GetFitter()->PrintResults(2, 0.);
  covM.ResizeTo(2,2);
  covM = TMatrixD(2, 2, TVirtualFitter::GetFitter()->GetCovarianceMatrix());
  covM.Print();
  const double correl = covM(0,1) / (TMath::Sqrt(covM(0,0)) * TMath::Sqrt(covM(1,1)));
  std::cout << "Correlation factor: " << correl << std::endl;
  alpha_var.setVal(bestX);
  mass_var .setVal(bestY);
  f1_alpha = projPdf->asTF(RooArgList(alpha_var), RooArgList(), RooArgSet(alpha_var));
  f1_mass  = projPdf->asTF(RooArgList(mass_var ), RooArgList(), RooArgSet(mass_var ));
  double alphaLowErr, alphaHighErr;
  double massLowErr, massHighErr;
  getUncertaintiesFromIntegral(f1_alpha, alphaLowErr, alphaHighErr, .10 , .13 , .01, .00001);
  getUncertaintiesFromIntegral(f1_mass , massLowErr , massHighErr , 140., 200., .1 , .01  );
  const double xBest[1] = {bestX};
  const double yBest[1] = {bestY};
  const double xBestErr[1] = {(alphaLowErr+alphaHighErr)/2};
  const double yBestErr[1] = {(massLowErr +massHighErr )/2};
  point = TGraphErrors(1, xBest, yBest, xBestErr, yBestErr);
  point.GetXaxis()->SetTitle(alpha_var.getTitle());
  point.GetYaxis()->SetTitle(mass_var .getTitle(true));
  point.SetMarkerStyle(29);
  point.SetMarkerSize(3.);
  ellipse = RooEllipse(label+"_ellipse", bestX, bestY, xBestErr[0], yBestErr[0], correl);
}

//  std::cout << mocTF->Eval(xMinTest-alphaLowErr, yMinTest) << " "
//	    << mocTF->Eval(xMinTest+alphaLowErr, yMinTest) << " "
//	    << mocTF->Eval(xMinTest, yMinTest-massLowErr) << " "
//	    << mocTF->Eval(xMinTest, yMinTest+massLowErr) << std::endl;
//
//  const int nLevels = 4;
//  const double levels[nLevels] = {mocTF->Eval(xMinTest+alphaLowErr, yMinTest),
//				  mocTF->Eval(xMinTest, yMinTest-massLowErr),
//				  mocTF->Eval(xMinTest-alphaLowErr, yMinTest),
//				  mocTF->Eval(xMinTest, yMinTest+massLowErr)};
//  mocTF->SetContour(nLevels, levels);

//  std::cout << "Latest Integral: " << mocTF->Integral(0.1, 0.13, 140., 190.) << std::endl;
//  std::cout << "Latest Integral: " << mocTF->Integral(xMinTest-alphaLowErr,
//						      xMinTest+alphaLowErr,
//						      yMinTest-massLowErr,
//						      yMinTest+massLowErr) << std::endl;


#endif

