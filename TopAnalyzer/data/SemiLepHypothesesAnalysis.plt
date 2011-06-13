$ -- configBlockIO -- $
$=====================================================================
$
$ input/output files
$
$=====================================================================

  histInput    = analyzeSemiLepHypotheses.hist                        $ hist file, which contains the names of the 
                                                                      $ histograms to be plottet 
                                                                      $
  rootInput    = analyzeSemiLepHypotheses.root                        $ root input file(s) corresponding to samples
                                                                      $ (separated by blanks)
                                                                      $
  inputDirs    = analyzeGeom                                          $ directory in which the histograms are kept
                 analyzeWMassMaxSumPt                                 $ within the root file; for edm::Analyzers 
                 analyzeMaxSumPtWMass                                 $ this  corresponds to the module name of the 
                 analyzeMVADisc                                       $ Analyzer in the cfg file (separated by blanks)
                 analyzeKinFit                                        $
                 analyzeGenMatch                                      $
                                                                      $
  histFilter   = kin qual                                             $ list of filter strings; only histograms con-
                                                                      $ taining these stings as substrings are loaded
                                                                      $ (separated by blanks); empty means everything
                                                                      $ is plotted 
                                                                      $
  plotFilter   =                                                      $ list of filter strings; only histograms con-
                                                                      $ taining these stings as substrings are plotted
                                                                      $ (separated by blanks); empty means everything
                                                                      $ is plotted
                                                                      $
  filterOption = begins                                               $ histogram filter option ('begins', 'ends' & 
                                                                      $ 'contains' are supported
                                                                      $
  rootOutput   =                                                      $ root output file to write histograms to
                                                                      $
  outputDir    =                                                      $ output directory in root output file 
                                                                      $
  outputLabels =                                                      $ output labels for multiple ps/eps files 
                                                                      $ (separated by blanks)
                                                                      $
  writePlotsTo = ./inspect_SemiLepHypotheses                          $ directory to save the plottet histograms to;
                                                                      $ '.' is the working directory
                                                                      $
  writePlotsAs = ps                                                   $ decide wether to write histograms to [ps] or 
                                                                      $ [eps] files

$ -- configBlockHist -- $
$=====================================================================
$
$ canvas and histogram steering
$
$=====================================================================

  xLog         =                                                      $ logs can be declared for each histogram 
                                                                      $ individually; per default they are  
  yLog         =                                                      $ switched off
                                                                      $
                                                                      $ grids can be declared for each histogram
                                                                      $ individually; per default they are 
                                                                      $ switched off
  xGrid        = 

  yGrid        = 
                                                                      $
  histScale    =                                                      $ histogram scale; can be steered for each 
                                                                      $ histogram individually; empty means scale
                                                                      $ is _1_
                                                                      $
  histMinimum  =                                                      $ histogram minimum; can be steered for 
                                                                      $ each histogram individually; empty means 
                                                                      $ set to _0_
                                                                      $
  histMaximum  =                                                      $ histogram maximum; can be steered for 
                                                                      $ each histogram individually; empty means
                                                                      $ set automatically
                                                                      $
  histErrors   =                                                      $ draw histogram errors; can be steered for 
                                                                      $ each sample; empty means no errors
                                                                      $
  histType     = 0 0 0 0 0 0                                          $ defines wether histogram should be plotted
                                                                      $ as line[0], with poly markers[1] or filled 
                                                                      $ [2]; can be steered for each sample individually
                                                                      $
  histStyle    = 1 1 1 1 1 1                                          $ defines line or fill style for each sample
                                                                      $
  histColor    = 2 3 4 6 7 28                                         $ defines line/marker/fill color for each 
                                                                      $ sample
                                                                      $
  lineWidth    = 5 5 5 5 5 5                                          $ defines line width for each sample
                                                                      $
  markerStyle  =                                                      $ defines marker style for each sample
                                                                      $
  markerSize   =                                                      $ defines marker size for each sample
                                                                      $
                                                                      $ set axes titles for each histograms;
                                                                      $ has to be given in '"' and end with ';'
  xAxes        = 

  yAxes        = 


                                                                      $
                                                                      $ should contain a legend entry for each 
                                                                      $ sample in '"' and separated by ';'
  legEntries   = "Geom";
                 "WMassMaxSumPt"; 
                 "MaxSumPtWMass";
                 "MVADisc";
                 "KinFit";
                 "GenMatch";
                                                                      $
  legXLeft     = 0.25                                                 $ so far a steering of the legend coord's 
                                                                      $ is only supported globally
  legXRight    = 0.95                                                 $
                                                                      $
  legYLower    = 0.70                                                 $
                                                                      $
  legYUpper    = 0.95                                                 $


