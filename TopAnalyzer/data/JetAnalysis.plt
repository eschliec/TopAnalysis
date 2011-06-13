$ -- configBlockIO -- $
$=====================================================================
$
$ input/output files
$
$=====================================================================

  histInput    = analyzeJets.hist                                     $ hist file, which contains the names of the 
                                                                      $ histograms to be plottet 
                                                                      $
  rootInput    = analyzeJets.root                                     $ root input file(s) corresponding to samples
                                                                      $ (separated by blanks)
                                                                      $
  inputDirs    = analyzeAllJets                                       $ directory in which the histograms are kept
                 analyzeSelJets                                       $ within the root file; for edm::Analyzers 
                                                                      $ this  corresponds to the module name of the 
                                                                      $ Analyzer in the cfg file (separated by blanks)
                                                                      $
  histFilter   = id pro kin                                           $ list of filter strings; only histograms con-
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
  outputFileName = analyzeJet								  		  $ gives the filename of the ps-output
  								  									  $ the filename is 'inspect' in the case of left blank


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
                 1 1 1 1 1 1 1 1 1 1
                 1 1 1 1 1 1 1 1 1 1

  yGrid        = 1 1 1 1 1 1 1 1 1 1
                 1 1 1 1 1 1 1 1 1 1   
                 1 1 1 1 1 1 1 1 1 1   
                 1 1 1 1 1 1 1 1 1 1 
                 1 1 1 1 1 1 1 1 1 1
                 1 1 1 1 1 1 1 1 1 1
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
  histType     =    2     2                                           $ defines wether histogram should be plotted
                                                                      $ as line[0], with poly markers[1] or filled 
                                                                      $ [2]; default is line [0]; can be steered 
                                                                      $ for each sample individually
                                                                      $
  histStyle    = 3004  3005                                           $ defines line or fill style for each sample
                                                                      $
  histColor    =    2     4                                           $ defines line/marker/fill color for each 
                                                                      $ sample
                                                                      $
  lineWidth    =    5     5                                           $ defines line width for each sample
                                                                      $
  markerStyle  =   20    21                                           $ defines marker style for each sample
                                                                      $
  markerSize   =  2.3   2.3                                           $ defines marker size for each sample
                                                                      $
                                                                      $ set axes titles for each histograms;
                                                                      $ has to be givin in '"' and end with ';'
  xAxes        = "em Fraction (all)";
                 "had Fraction (all)";
                 "E_{HO} (all) [GeV]";
                 "em Fraction (1.Jet)";
                 "had Fraction (1.Jet)";
                 "E_{HO} (1.Jet) [GeV]";
                 "em Fraction (2.Jet)";
                 "had Fraction (2.Jet)";
                 "E_{HO} (2.Jet) [GeV]";
                 "em Fraction (3.Jet)";
                 "had Fraction (3.Jet)";
                 "E_{HO} (3.Jet) [GeV]";
                 "em Fraction (4.Jet)";
                 "had Fraction (4.Jet)";
                 "E_{HO} (4.Jet) [GeV]";
                 "N_{Tower}/Jet (all)";
                 "(#eta_{Jet}-#eta_{Tower})(all)";
                 "(#phi_{Jet}-#phi_{Tower})(all)";
                 "(#eta_{Jet}-#eta_{Tower})(1.Jet)";
                 "(#phi_{Jet}-#phi_{Tower})(1.Jet)";
                 "(#eta_{Jet}-#eta_{Tower})(2.Jet)";
                 "(#phi_{Jet}-#phi_{Tower})(2.Jet)";
                 "(#eta_{Jet}-#eta_{Tower})(3.Jet)";
                 "(#phi_{Jet}-#phi_{Tower})(3.Jet)";
                 "(#eta_{Jet}-#eta_{Tower})(4.Jet)";
                 "(#phi_{Jet}-#phi_{Tower})(4.Jet)";
                 "E(all) [GeV]";
                 "E_{t}(all) [GeV]";
                 "#eta(all)";
                 "#phi(all)";
                 "#Sigma_{1.Jet-4.Jet}(E_{t}) [GeV]";
                 "#Sigma_{1.Jet-6.Jet}(E_{t}) [GeV]";
                 "E(1.Jet) [GeV]";
                 "E_{t}(1.Jet) [GeV]";
                 "#eta(1.Jet)";
                 "#phi(1.Jet)";
                 "E(2.Jet) [GeV]";
                 "E_{t}(2.Jet) [GeV]";
                 "#eta(2.Jet)";
                 "#phi(2.Jet)";
                 "E(3.Jet) [GeV]";
                 "E_{t}(3.Jet) [GeV]";
                 "#eta(3.Jet)";
                 "#phi(3.Jet)";
                 "E(4.Jet) [GeV]";
                 "E_{t}(4.Jet) [GeV]";
                 "#eta(4.Jet)";
                 "#phi(4.Jet)";
                 "N_{Jets}(all)";
                 "N_{Jets}(E_{t}>10 GeV)";
                 "N_{Jets}(E_{t}>20 GeV)";
                 "N_{Jets}(E_{t}>30 GeV)";
                 "N_{Jets}(E_{t}>40 GeV)";

  yAxes        = "jets";
                 "jets";
                 "jets";
                 "jets";
                 "jets";
                 "jets";
                 "jets";
                 "jets";
                 "jets";
                 "jets";
                 "jets";
                 "jets";
                 "jets";
                 "jets";
                 "jets";
                 "jets";
                 "1/N_{Jet} #Sigma(E_{t}^{Twr}/E_{t}^{Jet})(all) [GeV]";
                 "1/N_{Jet} #Sigma(E_{t}^{Twr})/E_{t}^{Jet})(all) [GeV]";
                 "1/N_{Evt} #Sigma(E_{t}^{Twr})/E_{t}^{Jet})(1.Jet) [GeV]";
                 "1/N_{Evt} #Sigma(E_{t}^{Twr})/E_{t}^{Jet})(1.Jet) [GeV]";
                 "1/N_{Evt} #Sigma(E_{t}^{Twr})/E_{t}^{Jet})(2.Jet) [GeV]";
                 "1/N_{Evt} #Sigma(E_{t}^{Twr})/E_{t}^{Jet})(2.Jet) [GeV]";
                 "1/N_{Evt} #Sigma(E_{t}^{Twr})/E_{t}^{Jet})(3.Jet) [GeV]";
                 "1/N_{Evt} #Sigma(E_{t}^{Twr})/E_{t}^{Jet})(3.Jet) [GeV]";
                 "1/N_{Evt} #Sigma(E_{t}^{Twr})/E_{t}^{Jet})(4.Jet) [GeV]";
                 "1/N_{Evt} #Sigma(E_{t}^{Twr})/E_{t}^{Jet})(4.Jet) [GeV]";
                 "jets";
                 "jets";
                 "jets";
                 "jets";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";
                 "events";



                                                                      $
                                                                      $ should contain a legend entry for each 
                                                                      $ sample in '"' and separated by ';'
  legEntries   = "All";                                     
                 "Selected";
                                                                      $
  legXLeft     = 0.25                                                 $ so far a steering of the legend coord's 
                                                                      $ is only supported globally
  legXRight    = 0.95                                                 $
                                                                      $
  legYLower    = 0.70                                                 $
                                                                      $
  legYUpper    = 0.95                                                 $
