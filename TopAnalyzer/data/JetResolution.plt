$ -- configBlockIO -- $
$=====================================================================
$
$ input/output files
$
$=====================================================================

  histInput    = analyzeJets.hist                                     $ hist file, which contains the names of the 
                                                                      $ histograms to be plottet 
  rootInput    = analyzeJets.root                                     $ root input file(s) corresponding to samples
$                                                                     $ (separated by blanks);

  inputDirs    = analyzeAllJets                                       $ directory in which the histograms are kept
$                analyzeSelJets                                       $ in within the root file; for edm::Analyzers 
                                                                      $ this  corresponds to the module name of the 
                                                                      $ Analyzer in the cfg file (separated by blanks)
  histFilter   = fit cal res                                          $ list of filter strings; only histograms con-
                                                                      $ taining this stings as substrings are plotted
                                                                      $ (separated by blanks)
  plotFilter   = res cal                                              $ list of filter strings; only histograms con-
                                                                      $ taining these stings as substrings are plotted
                                                                      $ (separated by blanks); empty means everything
                                                                      $ is plotted
                                                                      $
  filterOption = begins                                               $ histogram filter option ('begins', 'ends' & 
                                                                      $ 'contains' are supported
  rootOutput   = fittedHists.root                                     $ root output file to write histograms to
                                                                      $
  outputDir    = analyzeAllJets                                       $ output directory in root output file 
$                analyzeSelJets                                       $
								      $
  outputLabels = allJets                                              $ output labels for multiple ps/eps files 
$                selJets                                              $ (separated by blanks)
                                                                      $
                                                                      $
  writePlotsTo = .                                                    $ directory to save the plottet histograms to;
                                                                      $ '.' is the working directory
  writePlotsAs = ps                                                   $ decide wether to write histograms to [ps] or 
                                                                      $ [eps] 

$ -- configBlockHist -- $
$=====================================================================
$
$ canvas and histogram steering
$
$=====================================================================

  xLog        =                                                       $ logs can be declared for each histogram 
                                                                      $ individually; per default they are  
  yLog        =                                                       $ switched off
                                                                      $
  xGrid       =                                                       $ grids can be declared for each histogram
                                                                      $ individually; per default they are 
  yGrid       =                                                       $ switched off

  histScale   =                                                       $ histogram scale; can be steered for each 
                                                                      $ histogram individually
  histMaximum =    						      $ histogram maximum; can be steered for 
                                                                      $ each histogram individually
  histMinimum =   			                              $ histogram minimum; can be steered for 
                                                                      $ each histogram individually

  histErrors  =    1                                                  $ draw histogram errors for given sample?

  histType    =    1                                                  $ defines wether histogram should be plotted
                                                                      $ as line[0], with poly markers[1] or filled 
                                                                      $ [2]; default is line [0]; can be steered 
                                                                      $ for each sample individually
  histStyle   =    1                                                  $ defines line or fill style for each sample

  histColor   =    2                                                  $ defines line/marker/fill color for each 
                                                                      $ sample
  lineWidth   =    5                                                  $ defines line width for each sample

  markerStyle =   20                                                  $ defines marker style for each sample

  markerSize  =  2.3                                                  $ defines marker size for each sample

                                                                      $ set axes titles of histograms;
  xAxes      = "E_{t,jet_{1}}^{gen}"; 				      $ has to be givin in " 
               "E_{t,jet_{1}}^{gen}"; 				      $ and has to end with;
               "E_{t,jet_{2}}^{gen}"; 				      $
               "E_{t,jet_{2}}^{gen}"; 				      $
               "E_{t,jet_{3}}^{gen}"; 				      $
               "E_{t,jet_{3}}^{gen}"; 				      $
               "E_{t,jet_{4}}^{gen}"; 				      $
               "E_{t,jet_{4}}^{gen}"; 				      $
               "E_{t,alljets}^{gen}"; 				      $
               "E_{t,alljets}^{gen}"; 				      $
               "eta_{jet_{1}}^{gen}"; 				      $
               "eta_{jet_{1}}^{gen}"; 				      $
               "eta_{jet_{2}}^{gen}"; 				      $
               "eta_{jet_{2}}^{gen}"; 				      $
               "eta_{jet_{3}}^{gen}"; 				      $
               "eta_{jet_{3}}^{gen}"; 				      $
               "eta_{jet_{4}}^{gen}"; 				      $
               "eta_{jet_{4}}^{gen}"; 				      $
               "eta_{alljets}^{gen}"; 				      $
               "eta_{alljets}^{gen}"; 				      $
               "phi_{jet_{1}}^{gen}"; 				      $
               "phi_{jet_{1}}^{gen}"; 				      $
               "phi_{jet_{2}}^{gen}"; 				      $
               "phi_{jet_{2}}^{gen}"; 				      $
               "phi_{jet_{3}}^{gen}"; 				      $
               "phi_{jet_{3}}^{gen}"; 				      $
               "phi_{jet_{4}}^{gen}"; 				      $
               "phi_{jet_{4}}^{gen}"; 				      $
               "phi_{alljets}^{gen}"; 				      $
               "phi_{alljets}^{gen}"; 				      $
								      $
  yAxes      = "mean of StabilizedGauss";			      $ has to be adapted to evaluation type
	       "sigma of StabilizedGauss";
	       "mean of StabilizedGauss";
	       "sigma of StabilizedGauss";
	       "mean of StabilizedGauss";
	       "sigma of StabilizedGauss";
	       "mean of StabilizedGauss";
	       "sigma of StabilizedGauss";
	       "mean of StabilizedGauss";
	       "sigma of StabilizedGauss";
	       "mean of StabilizedGauss";
	       "sigma of StabilizedGauss";
	       "mean of StabilizedGauss";
	       "sigma of StabilizedGauss";
	       "mean of StabilizedGauss";
	       "sigma of StabilizedGauss";
	       "mean of StabilizedGauss";
	       "sigma of StabilizedGauss";
	       "mean of StabilizedGauss";
	       "sigma of StabilizedGauss";
	       "mean of StabilizedGauss";
	       "sigma of StabilizedGauss";
	       "mean of StabilizedGauss";
	       "sigma of StabilizedGauss";
	       "mean of StabilizedGauss";
	       "sigma of StabilizedGauss";
	       "mean of StabilizedGauss";
	       "sigma of StabilizedGauss";
	       "mean of StabilizedGauss";
	       "sigma of StabilizedGauss";

  legEntries = "all jets";                                            $ should contain a legend entry for each 
               "after jet kin";                                       $ sample in " and separated by ;
                                                                      $
                                                                      $
  legXLeft   = 0.35                                                   $ so far a steering of the legend coord's 
                                                                      $ is only supported globally
  legXRight  = 0.85                                                   $
                                                                      $
  legYLower  = 0.70

  legYUpper  = 0.95


$ -- configBlockFit -- $
$=====================================================================
$
$ histogram fitting
$
$=====================================================================

  evalType         = 0                                                $ mean/sigma evaluation type;
                                                                      $ 0 is StabilizedGauss
                                                                      $ 1 is HistogramMean
                                                                      $ 2 is MaximalValue
                                                                      $ 3 is Quantile
  targetLabel      = cal res                                          $ contains a list of prefixes of target
                                                                      $ histograms to write the fit results to
  fitFunctionName  = func                                             $ fit function name
                                                                      $
  fitFunctionTitle = Fit(Gauss)                                       $ title of the fit function to be shown in
                                                                      $ the legend e.g.
  fitFunctionType  = 0                                                $ fit function type; 0 is stable fitted sngl
                                                                      $ gauss
  fitLowerBound    = -2                                               $ lower boundary of fit range; to be adjusted
                                                                      $ for ratio and difference
  fitUpperBound    = 2.                                               $ upper boundary of fit range; to be adjusted
                                                                      $ for ratio and difference
								      $
  xAxesFit         = "(E_{t,jet_{1}}^{rec}-E_{t,jet_{1}}^{gen})/E_{t,jet_{1}}^{gen}";$ axes labels for fit histograms;
                     "(E_{t,jet_{2}}^{rec}-E_{t,jet_{2}}^{gen})/E_{t,jet_{2}}^{gen}";$ keep different labels here 
                     "(E_{t,jet_{3}}^{rec}-E_{t,jet_{3}}^{gen})/E_{t,jet_{3}}^{gen}";$ to adress them in the
                     "(E_{t,jet_{4}}^{rec}-E_{t,jet_{4}}^{gen})/E_{t,jet_{4}}^{gen}";$ following variable (titleIndex),
                     "(E_{t,alljets}^{rec}-E_{t,alljets}^{gen})/E_{t,alljets}^{gen}";$ should be in " and end with ;
                     "(eta_{jet_{1}}^{rec}-eta_{jet_{1}}^{gen})/eta_{jet_{1}}^{gen}";$
                     "(eta_{jet_{2}}^{rec}-eta_{jet_{2}}^{gen})/eta_{jet_{2}}^{gen}";$
                     "(eta_{jet_{3}}^{rec}-eta_{jet_{3}}^{gen})/eta_{jet_{3}}^{gen}";$
                     "(eta_{jet_{4}}^{rec}-eta_{jet_{4}}^{gen})/eta_{jet_{4}}^{gen}";$
                     "(eta_{alljets}^{rec}-eta_{alljets}^{gen})/eta_{alljets}^{gen}";$
                     "(phi_{jet_{1}}^{rec}-phi_{jet_{1}}^{gen})/phi_{jet_{1}}^{gen}";$
                     "(phi_{jet_{2}}^{rec}-phi_{jet_{2}}^{gen})/phi_{jet_{2}}^{gen}";$
                     "(phi_{jet_{3}}^{rec}-phi_{jet_{3}}^{gen})/phi_{jet_{3}}^{gen}";$
                     "(phi_{jet_{4}}^{rec}-phi_{jet_{4}}^{gen})/phi_{jet_{4}}^{gen}";$
                     "(phi_{alljets}^{rec}-phi_{alljets}^{gen})/phi_{alljets}^{gen}";$
									               $
                     						      $
  yAxesFit         =                                                  $ following variable (titleIndex), should be
                                                                      $ in " and end with ;
  titleIndex       = 0 0 0 0 0 0 0                                    $ titleIndex of x/yAxesFit for each if the
                     1 1 1 1 1 1 1                                    $ the fit histograms
                     2 2 2 2 2 2 2                                    $
                     3 3 3 3 3 3 3
		     4 4 4 4 4 4 4
		     5 5 5 5 5 5
		     6 6 6 6 6 6
		     7 7 7 7 7 7
		     8 8 8 8 8 8
		     9 9 9 9 9 9
		     10 10 10 10
		     11 11 11 11
		     12 12 12 12
		     13 13 13 13
		     14 14 14 14