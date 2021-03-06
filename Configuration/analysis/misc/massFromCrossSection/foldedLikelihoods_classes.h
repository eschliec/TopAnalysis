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
	   const TF1* xsec_func, const std::vector<TF1*>* unc_funcs, TFile* alpha_funcFile, RooRealVar& alpha_def);

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

  const QuadraticDependence alpha_p0;
  const QuadraticDependence alpha_p1;
  const QuadraticDependence alpha_p2;

  const RooPolyVar alphaDep;
  const RooPolyVar alphaDepCorr;

public:

  RooFormulaVar xsec;
  RooFormulaVar xsecScaleUp;
  RooFormulaVar xsecScaleDown;
  RooFormulaVar gaussianUnc;

private:

  RooNumIntConfig mcIntegratorCfg;

public:

  RooGaussian gaussianProb;
  RooGenericPdf rectangularProb;
  RooGenericPdf prob;
  //  RooGaussian prob;

};

PredXSec::PredXSec(const TString& label, RooRealVar& xsec_var, RooRealVar& mass_var, RooRealVar& alpha_var,
		   const TF1* xsec_func, const std::vector<TF1*>* unc_funcs, TFile* alpha_funcFile,
		   RooRealVar& alpha_def):
  name(label),
  p0(label+"_p0", label+"_p0", xsec_func->GetParameter(0)),
  p1(label+"_p1", label+"_p1", xsec_func->GetParameter(1)),
  p2(label+"_p2", label+"_p2", xsec_func->GetParameter(2)),
  p3(label+"_p3", label+"_p3", xsec_func->GetParameter(3)),
  relUncPdf      (label+"_relUncPdf"      , unc_funcs[1].at(0), mass_var),
  relUncScaleUp  (label+"_relUncScaleUp"  , unc_funcs[0].at(0), mass_var),
  relUncScaleDown(label+"_relUncScaleDown", unc_funcs[0].at(1), mass_var),
  alpha_p0(label+"_alpha_p0", ((TGraph*)alpha_funcFile->Get("graph_p0"))->GetFunction("pol2"), mass_var),
  alpha_p1(label+"_alpha_p1", ((TGraph*)alpha_funcFile->Get("graph_p1"))->GetFunction("pol2"), mass_var),
  alpha_p2(label+"_alpha_p2", ((TGraph*)alpha_funcFile->Get("graph_p2"))->GetFunction("pol2"), mass_var),
  alphaDep(label+"_alphaDep", label+"_alphaDep", alpha_var, RooArgSet(alpha_p0.polyVar,
								      alpha_p1.polyVar,
								      alpha_p2.polyVar)),
  alphaDepCorr(label+"_alphaDepCorr", label+"_alphaDepCorr", alpha_def, RooArgSet(alpha_p0.polyVar,
										  alpha_p1.polyVar,
										  alpha_p2.polyVar)),
  xsec(label+"_xsec", label+"_xsec", "(@1+@2*@0+@3*@0*@0+@4*@0*@0*@0)*(@5/@6)/(@0*@0*@0*@0)",
       RooArgSet(mass_var, p0, p1, p2, p3, alphaDep, alphaDepCorr)),
  //  xsec(label+"_xsec", label+"_xsec", "(@1+@2*@0+@3*@0*@0+@4*@0*@0*@0)*(@5/@6)/(@0*@0*@0*@0)*(1+TMath::Abs(@7))",
  //       RooArgSet(mass_var, p0, p1, p2, p3, alphaDep, alphaDepCorr, relUncScaleUp.polyVar)),
  //  xsec(label+"_xsec", label+"_xsec", "(@1+@2*@0+@3*@0*@0+@4*@0*@0*@0)*(@5/@6)/(@0*@0*@0*@0)*(1-TMath::Abs(@7))",
  //       RooArgSet(mass_var, p0, p1, p2, p3, alphaDep, alphaDepCorr, relUncScaleDn.polyVar)),
  xsecScaleUp(label+"_xsecScaleUp", label+"_xsecScaleUp", "@0+@0*TMath::Abs(@1)",
	      RooArgSet(xsec, relUncScaleUp.polyVar)),
  xsecScaleDown(label+"_xsecScaleDown", label+"_xsecScaleDown", "@0-@0*TMath::Abs(@1)",
		RooArgSet(xsec, relUncScaleDown.polyVar)),
  gaussianUnc(label+"_gaussianUnc", label+"_gaussianUnc", "@0*TMath::Abs(@1)",
	      RooArgSet(xsec, relUncPdf.polyVar)),
  //  gaussianUnc(label+"_gaussianUnc", label+"_gaussianUnc",
  //  "@0*TMath::Sqrt(TMath::Abs(@1)*TMath::Abs(@1)+TMath::Power(TMath::Max(TMath::Abs(@2),TMath::Abs(@3)),2))", //max
  //  "@0*TMath::Sqrt(TMath::Abs(@1)*TMath::Abs(@1)+TMath::Power((TMath::Abs(@2)+TMath::Abs(@3))/2.,2))",        //sym
  //	      RooArgSet(xsec, relUncPdf.polyVar, relUncScaleUp.polyVar, relUncScaleDown.polyVar)),
	      //	      RooArgSet(xsec, relUncPdf.polyVar)),
  mcIntegratorCfg(*RooAbsReal::defaultIntegratorConfig()),
  gaussianProb(label+"_gaussianProb", label+"_gaussianProb", xsec_var, xsec, gaussianUnc),
  rectangularProb(label+"_rectangularProb", label+"_rectangularProb",
		  "(@0 >= @1) && (@0 < @2)", RooArgList(xsec_var, xsecScaleDown, xsecScaleUp)),
  prob(label+"_prob", label+"_prob",
       "1/(2*(@3-@2))*(TMath::Erf((@3-@0)/(@1*TMath::Sqrt(2)))-TMath::Erf((@2-@0)/(@1*TMath::Sqrt(2))))",
       RooArgList(xsec_var, gaussianUnc, xsecScaleDown, xsecScaleUp))
  //  prob(label+"_prob", label+"_prob", xsec_var, xsec, gaussianUnc)
{
  mcIntegratorCfg.method1D().setLabel("RooMCIntegrator");
  rectangularProb.setIntegratorConfig(mcIntegratorCfg);
}

//////////////////////////////////////////////////
/// Final likelihood and extracted results
/// for either m_top or alpha_S
//////////////////////////////////////////////////

class FinalLikeliResults1D {

public:
  
  FinalLikeliResults1D(const TString& label, RooRealVar& meas_xsec, RooRealVar& xsec_var, RooRealVar& target_var,
		       const RooArgList& prodPdfList, RooRealVar& constrained_var,
		       const RooRealVar& constrained_var_mean, const RooRealVar& constrained_var_unc);

  void addPointToGraphs(TGraphAsymmErrors& graphInnerError, TGraphAsymmErrors& graphTotalError,
			const unsigned iPoint, const double y) const;

private:

  RooProdPdf prodPdf;

public:

  RooAbsPdf* projPdf;
  TF1* f1;
  double bestX;
  double lowErrFromIntegral;
  double lowErrFromConstraintUncertainty;
  double lowErrFromBeamUncertainty;
  double lowErrTotal;
  double highErrFromIntegral;
  double highErrFromConstraintUncertainty;
  double highErrFromBeamUncertainty;
  double highErrTotal;

};

FinalLikeliResults1D::FinalLikeliResults1D(const TString& label, RooRealVar& meas_xsec,
					   RooRealVar& xsec_var, RooRealVar& target_var, const RooArgList& prodPdfList,
					   RooRealVar& constrained_var,
					   const RooRealVar& constrained_var_mean,
					   const RooRealVar& constrained_var_unc):
  prodPdf(label+"_prodPdf", label+"_prodPdf", prodPdfList)
{
  projPdf = prodPdf.createProjection(xsec_var);
  constrained_var.setVal(constrained_var_mean.getVal());
  f1 = projPdf->asTF(RooArgList(target_var), RooArgList(), RooArgSet(target_var));
  const bool quickDebug = false;
  if(quickDebug) {
    bestX = f1->GetMaximumX();
    lowErrFromIntegral  = (!strcmp(target_var.GetName(),"mass") ? 4. : 0.002);
    highErrFromIntegral = (!strcmp(target_var.GetName(),"mass") ? 4. : 0.002);
  }
  else {
    if(!strcmp(target_var.GetName(),"mass"))
      bestX = getMaxAndUncertaintiesFromIntegral(f1, lowErrFromIntegral, highErrFromIntegral, 140., 190., .1 , .01  );
    else if(!strcmp(target_var.GetName(),"alpha"))
      bestX = getMaxAndUncertaintiesFromIntegral(f1, lowErrFromIntegral, highErrFromIntegral, .10 , .13 , .01, .00001);
    else {
      std::cout << "Target variable " << target_var.GetName() << " not supported in " << label << "!" << std::endl;
      abort();
    }
  }
  constrained_var.setVal(constrained_var_mean.getVal()-constrained_var_unc.getVal());
  double variationA = f1->GetMaximumX();
  constrained_var.setVal(constrained_var_mean.getVal()+constrained_var_unc.getVal());
  double variationB = f1->GetMaximumX();
  lowErrFromConstraintUncertainty  = bestX - TMath::Min(TMath::Min(variationA, variationB), bestX);
  highErrFromConstraintUncertainty = TMath::Max(TMath::Max(variationA, variationB), bestX) - bestX;
  constrained_var.setVal(constrained_var_mean.getVal());
  const double measXSecSave = meas_xsec.getVal();
  meas_xsec.setVal(1.018*measXSecSave);
  variationA = f1->GetMaximumX();
  meas_xsec.setVal(0.982*measXSecSave);
  variationB = f1->GetMaximumX();
  lowErrFromBeamUncertainty  = bestX - TMath::Min(TMath::Min(variationA, variationB), bestX);
  highErrFromBeamUncertainty = TMath::Max(TMath::Max(variationA, variationB), bestX) - bestX;
  meas_xsec.setVal(measXSecSave);
  lowErrTotal  = TMath::Sqrt(TMath::Power(lowErrFromIntegral , 2) + TMath::Power(lowErrFromConstraintUncertainty , 2)
			     + TMath::Power(lowErrFromBeamUncertainty , 2));
  highErrTotal = TMath::Sqrt(TMath::Power(highErrFromIntegral, 2) + TMath::Power(highErrFromConstraintUncertainty, 2)
			     + TMath::Power(highErrFromBeamUncertainty , 2));
}

void
FinalLikeliResults1D::addPointToGraphs(TGraphAsymmErrors& graphInnerError, TGraphAsymmErrors& graphTotalError,
				       const unsigned iPoint, const double y) const
{
  graphInnerError.SetPoint(iPoint, bestX, y);
  graphTotalError.SetPoint(iPoint, bestX, y);
  graphInnerError.SetPointError(iPoint, lowErrFromIntegral, highErrFromIntegral, 0., 0.);
  graphTotalError.SetPointError(iPoint, lowErrTotal       , highErrTotal       , 0., 0.);
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

class FinalLikeliResults2D {

public:
  
  FinalLikeliResults2D(const TString& label, RooRealVar& xsec_var, RooRealVar& mass_var, RooRealVar& alpha_var,
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

FinalLikeliResults2D::FinalLikeliResults2D(const TString& label,
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
  getMaxAndUncertaintiesFromIntegral(f1_alpha, alphaLowErr, alphaHighErr, .10 , .13 , .01, .00001);
  getMaxAndUncertaintiesFromIntegral(f1_mass , massLowErr , massHighErr , 140., 200., .1 , .01  );
  const double xBest[1] = {bestX};
  const double yBest[1] = {bestY};
  const double xBestErr[1] = {(alphaLowErr+alphaHighErr)/2};
  const double yBestErr[1] = {(massLowErr +massHighErr )/2};
  point = TGraphErrors(1, xBest, yBest, xBestErr, yBestErr);
  point.GetXaxis()->SetTitle(alpha_var.getTitle());
  point.GetYaxis()->SetTitle(mass_var .getTitle(true));
  point.SetMarkerSize(3.);
  ellipse = RooEllipse(label+"_ellipse", bestX, bestY, xBestErr[0], yBestErr[0], correl);
}

#endif

