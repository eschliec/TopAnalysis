$ -- configBlockIO -- $
$=====================================================================
$
$ input/output files
$
$=====================================================================

  histInput    = analyzeMuon.hist                                     $ hist file, which contains the names of the 
                                                                      $ histograms to be plottet 
                                                                      $
  rootInput    = analyzeMuon.root                                     $ root input file(s) corresponding to samples
                                                                      $ (separated by blanks)
                                                                      $
  inputDirs    = analyzeAllMuon                                       $ directory in which the histograms are kept
                 analyzeSelMuon                                       $ within the root file; for edm::Analyzers 
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
                                                                      
  outputFileName = analyzeMuon										$ gives the filename of the ps-output
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
  xAxes        = "#mu compatibility";
                 "E_{em} [GeV]";
                 "E_{em}^{3X3} [GeV]";
                 "E_{had} [GeV]";
                 "E_{had}^{3x3} [GeV]";
                 "E_{HO} [GeV]";
                 "E_{HO}^{3x3} [GeV]";
                 "E(muon) [GeV]";
                 "p_{t}(muon) [GeV]";
                 "#eta(muon)";
                 "#phi(muon)";
                 "min(#delta R (muon, jet))_{all}";
                 "min(#delta R (muon, jet))_{p_{t}>5 GeV}";
                 "min(#delta R (muon, jet))_{p_{t}>10 GeV}";
                 "min(#delta R (muon, jet))_{p_{t}>15 GeV}";
                 "min(#delta R (muon, jet))_{p_{t}>20 GeV}";
                 "#Sigma_{#delta R<0.3}(p_{t}(track))";
                 "#Sigma_{#delta R<0.3}(E_{t}(tower))";
                 "#Sigma_{#delta R<0.3}(E_{t}^{em}(tower))";
                 "#Sigma_{#delta R<0.3}(E_{t}^{had}(tower))";
                 "dR(muon,track)";
                 "dR(muon,track)";
                 "dR(muon,tower)";
                 "dR(muon,tower)";
                 "N_{track}^{dR<0.3}";
                 "N_{tower}^{dR<0.3}";
                 "log_{10}(#delta(cft,combined))";

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
                 "1/N_{muon} #Sigma (p_{t}(track)) [GeV]";
                 "1/N_{muon} #Sigma (N_{track})";
                 "1/N_{muon} #Sigma (E_{t}(tower)) [GeV]";
                 "1/N_{muon} #Sigma (N_{tower})";
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
