$ -- configBlockIO -- $
$=====================================================================
$
$ input/output files
$
$=====================================================================

  histInput    = analyzeElec.hist                                     $ hist file, which contains the names of the 
                                                                      $ histograms to be plottet 
                                                                      $
  rootInput    = analyzeElec.root                                     $ root input file(s) corresponding to samples
                                                                      $ (separated by blanks)
                                                                      $
  inputDirs    = analyzeAllElec                                       $ directory in which the histograms are kept
                 analyzeSelElec                                       $ within the root file; for edm::Analyzers 
                                                                      $ this  corresponds to the module name of the 
                                                                      $ Analyzer in the cfg file (separated by blanks)
                                                                      $
  histFilter   = id kin iso                                           $ list of filter strings; only histograms con-
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
  xAxes        = "(E/p)_{SC}(electron)";
                 "(E/P)_{BC}(electron)";
                 "Number of clusters";
                 "#delta #eta (track, cluster)";
                 "#delta #phi (track, cluster)";
                 "log_{10}(E_{had}/E_{em})";
                 "(p_{t}^{inner}-p_{t}^{outer})/p_{t}^{inner}(track)";
                 "(#eta^{inner} - #eta^{outer})(track)";
                 "(#phi^{inner} - #phi^{outer})(track)";
                 "N_{Hits}(track)";
                 "N_{Hits}^{used}(track)";
                 "#chi^{2} Track Fit";
                 "dca(track, cluster)";
                 "S1OverS9";
                 "E(electron) [GeV]";
                 "p_{t}(electron) [GeV]";
                 "#eta(electron)";
                 "#phi(electron)";
                 "min(#delta R (electron, jet))";
                 "#Sigma_{#delta R<0.3}(p_{t}(track))";
                 "#Sigma_{#delta R<0.3}(E_{t}(tower))";
                 "#Sigma_{#delta R<0.3}(E_{t}^{em}(tower))";
                 "#Sigma_{#delta R<0.3}(E_{t}^{had}(tower))";
                 "dR(electron,track)";
                 "dR(electron,track)";
                 "dR(electron,tower)";
                 "dR(electron,tower)";
                 "N_{track}^{dR<0.3}";
                 "N_{tower}^{dR<0.3}";
                 "#delta(cft,gsf)";

  yAxes        = "events";
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
                 "1/N_{electron} #Sigma (p_{t}(track)) [GeV]";
                 "1/N_{electron} #Sigma (N_{track})";
                 "1/N_{electron} #Sigma (E_{t}(tower)) [GeV]";
                 "1/N_{electron} #Sigma (N_{tower})";
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
