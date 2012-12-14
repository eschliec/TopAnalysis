#include "basicFunctions.h"

void createNNLOplot(TString outputRootFile="AhrensNNLL.root")
{
  gStyle->SetOptStat(0);
  // list of variables
  bool usequad=true;
  std::vector<TString> xSecVariables_, xSecLabel_;
  TString xSecVariablesUsed[] ={"ttbarMass"};
  xSecVariables_ .insert( xSecVariables_.begin(), xSecVariablesUsed, xSecVariablesUsed + sizeof(xSecVariablesUsed)/sizeof(TString) );
  // binning
  std::map<TString, std::vector<double> > bins_=makeVariableBinning();
  double low =bins_["ttbarMass"][0];
  double high=bins_["ttbarMass"][bins_["ttbarMass"].size()-1];

  // gotta get em all
  for(unsigned var=0; var<xSecVariables_.size(); ++var){
    // get raw nnlo points
    TString variable=xSecVariables_[var];
    TGraph * rawHist = new TGraph("AhrensTheory_Mtt_7000_172.5_Mtt_norm.dat");
    // convert to TH1F*
    //double *Xvalues = rawHist->GetX();
    int Nbins=rawHist->GetN();
    double *Yvalues = rawHist->GetY();
    double xmax=2720.;
    double xmin=345.;
    double binwidth=25.;
    //
    TH1F* hist = new TH1F ( variable, variable, (int)((xmax-xmin)/binwidth), xmin, xmax); 
    std::cout << "load values from .dat file: " << std::endl;
    for(int bin=1; bin<=Nbins; ++bin){
      std::cout << "bin " << bin;
      std::cout << "(<x>=" << hist->GetBinCenter(bin) << ", ";
      std::cout << hist->GetBinLowEdge(bin) << ".." << hist->GetBinLowEdge(bin+1);
      std::cout << "): " << Yvalues[bin-1] << std::endl;
      hist->SetBinContent(bin, Yvalues[bin-1]);
    }
    // check theory curve
    std::cout << std::endl << "analyze theory curve:" << std::endl;
    double integralRawTheory=hist->Integral(0.,hist->GetNbinsX()+1);
    std::cout << "Integral: " << integralRawTheory << std::endl;
    std::cout << "binwidth: " << binwidth << std::endl;
    std::cout << "Integral*binwidth: " << integralRawTheory*binwidth << std::endl;
    std::cout << "binRange: " << xmin << ".." << xmax << std::endl;
    std::cout << "-> reco range: " << low << ".." << high << std::endl;
    // styling
    hist->GetXaxis()->SetRange(low,high);
    hist->SetMarkerColor(kMagenta);
    hist->SetMarkerStyle(24);
    // create canvas
    std::cout << std::endl << "create canvas " << std::endl;
    TCanvas *canv = new TCanvas("canv","canv",800,600);
    canv->cd();
    std::cout << "draw original theory curve " << std::endl;
    hist->Draw("p");
    hist->Draw("hist same");
    // create rebinned plot
    std::cout << std::endl << "create rebinned curve:" << std::endl;
    TH1F* binnedPlot=new TH1F("nnll", "nnll", bins_["ttbarMass"].size()-1, &bins_["ttbarMass"][0]);
    for(int bin=1; bin<=hist->GetNbinsX(); ++bin){
      double y=hist->GetBinContent(bin)*hist->GetBinWidth(bin);
      double xlow=hist->GetBinLowEdge(bin);
      double xhigh=hist->GetBinLowEdge(bin+1);
      // search corresponding bin in rebinned curve
      bool found=false;
      for(int binRebinned=1; binRebinned<=binnedPlot->GetNbinsX(); ++binRebinned){
	if(binnedPlot->GetBinLowEdge(binRebinned)<=xlow&&binnedPlot->GetBinLowEdge(binRebinned+1)>=xhigh){
	  found=true;
	  binnedPlot->SetBinContent(binRebinned, binnedPlot->GetBinContent(binRebinned)+y);
	  break;
	}
      }
      // use linear interpolation for edge bins
      if(hist->GetBinCenter(bin)<high&&!found){
	std::cout << "need interpolation for bin" << bin << "(<x>=" << hist->GetBinCenter(bin) << ")!"<< std::endl;
	// search for the two bins involved
	double binLow=0;
	double binHigh=0;
	for(int binRebinned=1; binRebinned<=binnedPlot->GetNbinsX(); ++binRebinned){
	  // search for bin in binned histo where upper border of is close to lower border of unbinned histo
	  if(std::abs(binnedPlot->GetBinLowEdge(binRebinned+1)-xlow)<binwidth){
	    binLow =binRebinned; 
	    binHigh=binRebinned+1;
	    break;
	  }
	}
	std::cout << "theory bin " << xlow << ".." << xhigh << "-> reco bins ";
	std::cout << binnedPlot->GetBinLowEdge(binLow) << ".." << binnedPlot->GetBinLowEdge(binLow+1) << " & ";
	std::cout << binnedPlot->GetBinLowEdge(binHigh) << ".." << binnedPlot->GetBinLowEdge(binHigh+1) << std::endl;
	// get the two points (this bin and the previous one)
	double x2=hist->GetBinCenter (bin  );
	double x1=hist->GetBinCenter (bin-1);
	double x3=hist->GetBinCenter (bin+1);
	double y2=hist->GetBinContent(bin  );
	double y1=hist->GetBinContent(bin-1);
	double y3=hist->GetBinContent(bin+1);	
	// calculate linear funtion
	double a=(y2-y1)/(x2-x1);
	double b=y1-a*x1;
	TF1* linInterpol=new TF1("linInterpol"+getTStringFromInt(bin), "[0]*x+[1]", x1, x2);
	linInterpol->SetParameter(0,a);
	linInterpol->SetParameter(1,b);
	// calculate the corresponding area of linear function to binned curve
	double contributionLowerBin=linInterpol->Integral(xlow,binnedPlot->GetBinLowEdge(binHigh));
	double contributionUpperBin=linInterpol->Integral(binnedPlot->GetBinLowEdge(binHigh),xhigh);
	// draw interpolation function for checking
	linInterpol->SetRange(xlow,xhigh);
	linInterpol->SetLineColor(kMagenta);
	linInterpol->DrawClone("same");
	if(usequad&&x2<450){
	  // calculate quadratic funtion
	  double a2=(y2-((y3-y1))*(x2-x1)/(x3-x1))/(x2*x2-x1*x1-(x2-x1)*(x3*x3+x1*x1)/(x3-x1));
	  double b2=((y3-y1)-a2*(x3*x3-x1*x1))/(x3-x1);
	  double c2=y1-a2*x1*x1-b*x1;
	  TF1* quadInterpol=new TF1("quadInterpol"+getTStringFromInt(bin), "[0]*x*x+[1]*x+[2]", x1, x2);
	  quadInterpol->SetParameter(0,a2);
	  quadInterpol->SetParameter(1,b2);
	  quadInterpol->SetParameter(2,c2);
	  // draw quad interpolation function for checking
	  quadInterpol->SetRange(xlow,xhigh);
	  quadInterpol->SetLineColor(kGreen);
	  quadInterpol->SetLineStyle(2);
	  hist->Fit(quadInterpol, "", "same", x1, x3);
	  // calculate the corresponding area of linear function to binned curve
	  quadInterpol->DrawClone("same");
	  double areaLow =quadInterpol->Integral(xlow,binnedPlot->GetBinLowEdge(binHigh) );
	  double areaHigh=quadInterpol->Integral(binnedPlot->GetBinLowEdge(binHigh),xhigh);
	  std::cout << "ratio(high/low) linear/quadratic: " << contributionUpperBin/contributionLowerBin << "/"<< areaHigh/areaLow << std::endl;
	  contributionLowerBin=y2*binwidth*areaLow /(areaLow+areaHigh);
	  contributionUpperBin=y2*binwidth*areaHigh/(areaLow+areaHigh);

	}
	// add fitted result
	binnedPlot->SetBinContent(binLow , binnedPlot->GetBinContent(binLow )+contributionLowerBin);
	binnedPlot->SetBinContent(binHigh, binnedPlot->GetBinContent(binHigh)+contributionUpperBin);
      }
    }
    // ensure normalization
    binnedPlot->Scale(1./(binnedPlot->Integral(0.,binnedPlot->GetNbinsX()+1)));
    // divide by binwidth
    for(int bin=1; bin<=binnedPlot->GetNbinsX(); ++bin){
      binnedPlot->SetBinContent(bin, binnedPlot->GetBinContent(bin)/binnedPlot->GetBinWidth(bin));
    }
    // styling
    binnedPlot->SetLineColor(kBlue);
    binnedPlot->SetLineWidth(2);
    std::cout << "draw rebinned theory curve " << std::endl;
    binnedPlot->Draw("hist same");
    // draw bin boundaries
    std::cout << "draw bin boundaries " << std::endl;
    for(int bin=0; bin<(int)bins_.size(); ++bin){
      drawLine(bins_["ttbarMass"][bin], 0, bins_["ttbarMass"][bin], hist->GetMaximum(), kRed, 2, 6);
    }
    std::cout << "done" << std::endl;
    // save in png and rootfile
    std::cout << std::endl << "do saving..." << std::endl;
    canv->SaveAs(variable+"Norm_ahrens.png");
    TH1F* out=(TH1F*)binnedPlot->Clone();
    out->SetTitle("ttbarMass");
    out->SetName ("ttbarMass");
    out->GetXaxis()->SetTitle(xSecLabelName("ttbarMass"));
    out->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dm^{t#bar{t}}} [GeV^{-1}]");
    out->SetLineColor(kOrange-3);
    std::cout << std::endl << "draw final result " << std::endl;
    TCanvas *canv2 = new TCanvas("canv2","canv2",800,600);
    canv2->cd();
    out->Draw();
    saveToRootFile(outputRootFile, out, true, 0,"");
    std::cout << "done!" << std::endl;
  }
}
