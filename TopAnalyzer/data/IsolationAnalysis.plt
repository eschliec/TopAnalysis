$ -- configBlockIO -- $
$=====================================================================
$
$ input/output files
$
$=====================================================================

  histInput    = Isolation.hist                                     $ hist file, which contains the names of the 
                                                                      $ histograms to be plottet 
                                                                      $
  rootInput    = rerun_240708_Matrix_qcd_0-15.root
rerun_240708_Matrix_qcd_1000-1400.root
rerun_240708_Matrix_qcd_120-170.root
rerun_240708_Matrix_qcd_15-20.root
rerun_240708_Matrix_qcd_170-230.root
rerun_240708_Matrix_qcd_20-30.root
rerun_240708_Matrix_qcd_230-300.root
rerun_240708_Matrix_qcd_30-50.root
rerun_240708_Matrix_qcd_300-380.root
rerun_240708_Matrix_qcd_380-470.root
rerun_240708_Matrix_qcd_470-600.root
rerun_240708_Matrix_qcd_50-80.root
 rerun_240708_Matrix_qcd_600-800.root
 rerun_240708_Matrix_qcd_80-120.root
 rerun_240708_Matrix_qcd_800-1000.root
$rerun_240708_Matrix_ttbar_0jet.root
 $ rerun_240708_Matrix_ttbar_1jet.root
  $rerun_240708_Matrix_ttbar_2jet.root
  $rerun_240708_Matrix_ttbar_3jet.root
  $rerun_240708_Matrix_ttbar_4jet.root
  $analyzeQCDBackground_wjet_240708.root                 $ root input file(s) corresponding to samples
                 $analyzeQCDBackground_zjet_240708.root
  $               analyzeQCDBackground_ttbar_240708.root
                 $analyzeQCDBackground_qcd_240708.root                                              $ (separated by blanks)
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
  histType     =    2     2   2   2 2 2 2 2 2  2 2 2 2 2 2 2 2 2 2                                  $ defines wether histogram should be plotted
                                                                      $ as line[0], with poly markers[1] or filled 
                                                                      $ [2]; default is line [0]; can be steered 
                                                                      $ for each sample individually
                                                                      $
  histStyle    = 3004  3005                                           $ defines line or fill style for each sample
                                                                      $
  histColor    =    2     4   6    8   3                                 $ defines line/marker/fill color for each 
                                                                      $ sample
                                                                      $
  lineWidth    =    2     2    2      2   2 2                              $ defines line width for each sample
                                                                      $
  markerStyle  =   20    21    22    23   23                                 $ defines marker style for each sample
                                                                      $
  markerSize   =  2.   2.  2. 2. 2.                                        $ defines marker size for each sample
                                                                      $
                                                                      $ set axes titles for each histograms;
                                                                      $ has to be givin in '"' and end with ';'
  xAxes        = $"missing ET";
                 $"missing ET";
                 $"missing ET";
                 "min #Delta #phi (MET,4Jets)";
                 "#Delta #phi (MET, 1st jet)";
                 "#Delta #phi (MET, 2nd jet)";
                 "#Delta #phi (MET, 3rd jet)";
                 "#Delta #phi (MET, 4th jet)";
                 "#Delta #phi (MET, leading muon)";
                 

  yAxes        = $"track isolation(GeV)";
                 $"calo isolation(GeV)";
                 $"leading muon pt";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 



                                                                      $
                                                                      $ should contain a legend entry for each 
                                                                      $ sample in '"' and separated by ';'
  legEntries   = "qcd 0-15";
  "qcd 1000-1400";
  "qcd 120-170";
  "qcd 15-20";
  "qcd 170-230";
  "qcd 20-30";
  "qcd 230-300";
  "qcd 30-50";
  "qcd 300-380";
  "qcd 380-470";
  "qcd 470-600";
  "qcd 50-80";
  "qcd 600-800";
  "qcd 80-120";
  "qcd 800-1000";
  $"ttbar 0j";
  				$"ttbar 1j";
  				$"ttbar 2j";
  				$"ttbar 3j";
  				$"ttbar 4j";  
  $"wjet";      
  				 $"zjet";
  				 $"ttbar";
  				 $"qcd";                               
                 
                                                                      $
  legXLeft     = 0.25                                                 $ so far a steering of the legend coord's 
                                                                      $ is only supported globally
  legXRight    = 0.95                                                 $
                                                                      $
  legYLower    = 0.70                                                 $
                                                                      $
  legYUpper    = 0.95                                                 $
  
  $ -- configBlockSum -- $
$=====================================================================
$
$ histogram addition
$
$=====================================================================

  histWeights = 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.                                                  $ define weights which should be applied
                                                                      $ to each sample before adding them up 
