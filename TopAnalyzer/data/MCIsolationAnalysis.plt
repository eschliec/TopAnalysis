$ -- configBlockIO -- $
$=====================================================================
$
$ input/output files
$
$=====================================================================

  histInput    = ttbarMC_Isolation.hist                                    $ hist file, which contains the names of the 
                                                                      $ histograms to be plottet 
                                                                      $
  rootInput    = rerun_240708_Matrix_ttbar_0jet.root
  rerun_240708_Matrix_ttbar_1jet.root
  rerun_240708_Matrix_ttbar_2jet.root
  rerun_240708_Matrix_ttbar_3jet.root
  rerun_240708_Matrix_ttbar_4jet.root
  $analyzeQCDBackground_ttbar_240708.root                 $ root input file(s) corresponding to samples
                                                              $ (separated by blanks)
                                                                      $
  inputDirs    = analyzeisolationMET                                       $ directory in which the histograms are kept
                 			                                          $ within the root file; for edm::Analyzers 
                                                                      $ this  corresponds to the module name of the 
                                                                      $ Analyzer in the cfg file (separated by blanks)
                                                                      $
  histFilter   =                                            $ list of filter strings; only histograms con-
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
  writePlotsTo = .                                                    $ directory to save the plottet histograms to;
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
  xGrid        = 1 1 1 1 1 1 1 1 1 1 
                 1 1 1 1 1 1 1 1 1 1   
                 1 1 1 1 1 1 1 1 1 1   
                 1 1 1 1 1 1 1 1 1 1

  yGrid        = 1 1 1 1 1 1 1 1 1 1
                 1 1 1 1 1 1 1 1 1 1   
                 1 1 1 1 1 1 1 1 1 1   
                 1 1 1 1 1 1 1 1 1 1 
                                                                      $
  histScale    =                                                   $ histogram scale; can be steered for each 
                                                                      $ histogram individually; empty means scale
                                                                      $ is _1_
                                                                      $
  histMinimum  =				                                  $ histogram minimum; can be steered for 
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
  histType     =    2   2   2   2   2                                  $ defines wether histogram should be plotted
                                                                      $ as line[0], with poly markers[1] or filled 
                                                                      $ [2]; default is line [0]; can be steered 
                                                                      $ for each sample individually
                                                                      $
  histStyle    = 3004  3005                                           $ defines line or fill style for each sample
                                                                      $
  histColor    =    2     4   6    8       10                             $ defines line/marker/fill color for each 
                                                                      $ sample
                                                                      $
  lineWidth    =    2     2    2      2     2                            $ defines line width for each sample
                                                                      $
  markerStyle  =   20    21    22    23    24                               $ defines marker style for each sample
                                                                      $
  markerSize   =  2.   2.  2. 2.         2.                                $ defines marker size for each sample
                                                                      $
                                                                      $ set axes titles for each histograms;
                                                                      $ has to be givin in '"' and end with ';'
  xAxes        = "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "missing ET";
                 "had. top pt";
                 "lep. top pt";
                 "had. top pt bins";
                 "lep. top pt bins";
                 "had. top pt bins";
                 "lep. top pt bins";
                 

  yAxes        = "track isolation(GeV)";
                 "calo isolation(GeV)";
                 "leading muon pt";
                 "track isolation(GeV)";
                 "calo isolation(GeV)";
                 "leading muon pt";
                 "track isolation(GeV)";
                 "calo isolation(GeV)";
                 "leading muon pt";
                 "track isolation(GeV)";
                 "calo isolation(GeV)";
                 "leading muon pt";
                 "track isolation(GeV)";
                 "calo isolation(GeV)";
                 "leading muon pt";
                 "track isolation(GeV)";
                 "calo isolation(GeV)";
                 "leading muon pt";
                 "track isolation(GeV)";
                 "calo isolation(GeV)";
                 "leading muon pt";
                 "track isolation(GeV)";
                 "calo isolation(GeV)";
                 "leading muon pt";
                 "events";
                 "events";
                 "calo correlation";
                 "calo correlation";
                 "track correlation";
                 "track correlation";
                 



                                                                      $
                                                                      $ should contain a legend entry for each 
                                                                      $ sample in '"' and separated by ';'
  legEntries   = "ttbar 0j";
  				"ttbar 1j";
  				"ttbar 2j";
  				"ttbar 3j";
  				"ttbar 4j";      
  				                               
                 
                                                                      $
  legXLeft     = 0.25                                                 $ so far a steering of the legend coord's 
                                                                      $ is only supported globally
  legXRight    = 0.95                                                 $
                                                                      $
  legYLower    = 0.70                                                 $
                                                                      $
  legYUpper    = 0.95                                                 $
