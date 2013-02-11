#!/bin/sh

## this shell script runs all macro parts of the 2012 differential cross section analysis 

## HOW TO USE THIS SHELL SCRIPT
## note: has to be called from .../TopAnalysis/Configuration/analysis/semiLeptonic/diffXSection
## a) set up the folder structure:
## mkdir -p diffXSecFromSignal/plots/muon/2012/monitoring/withoutRatioPlots
## mkdir -p diffXSecFromSignal/plots/muon/2012/partonLevel
## mkdir -p diffXSecFromSignal/plots/muon/2012/recoYield
## mkdir -p diffXSecFromSignal/plots/muon/2012/uncertainties
## mkdir -p diffXSecFromSignal/plots/muon/2012/uncertaintyDistributions
## mkdir -p diffXSecFromSignal/plots/muon/2012/uncertaintyDistributionsOverview
## mkdir -p diffXSecFromSignal/plots/muon/2012/xSec
## mkdir -p diffXSecFromSignal/plots/muon/2012/binning
## mkdir -p diffXSecFromSignal/plots/muon/2012/effAndAcc
## mkdir -p diffXSecFromSignal/plots/muon/2012/genRecoCorrPlots
## mkdir -p diffXSecFromSignal/plots/muon/2012/kinFitPerformance
## mkdir -p diffXSecFromSignal/plots/muon/2012/shapeReweighting
## mkdir -p diffXSecFromSignal/plots/muon/2012/unfolding
## mkdir -p diffXSecFromSignal/plots/electron/2012/monitoring/withoutRatioPlots
## mkdir -p diffXSecFromSignal/plots/electron/2012/partonLevel
## mkdir -p diffXSecFromSignal/plots/electron/2012/recoYield
## mkdir -p diffXSecFromSignal/plots/electron/2012/uncertainties
## mkdir -p diffXSecFromSignal/plots/electron/2012/uncertaintyDistributions
## mkdir -p diffXSecFromSignal/plots/electron/2012/uncertaintyDistributionsOverview
## mkdir -p diffXSecFromSignal/plots/electron/2012/xSec
## mkdir -p diffXSecFromSignal/plots/electron/2012/binning
## mkdir -p diffXSecFromSignal/plots/electron/2012/effAndAcc
## mkdir -p diffXSecFromSignal/plots/electron/2012/genRecoCorrPlots
## mkdir -p diffXSecFromSignal/plots/electron/2012/kinFitPerformance
## mkdir -p diffXSecFromSignal/plots/electron/2012/shapeReweighting 
## mkdir -p diffXSecFromSignal/plots/electron/2012/unfolding
## mkdir -p diffXSecFromSignal/plots/combined/2012/xSec
## mkdir -p diffXSecFromSignal/plots/combined/2012/binning
## mkdir -p diffXSecFromSignal/plots/combined/2012/partonLevel
## mkdir -p diffXSecFromSignal/plots/combined/2012/recoYield
## mkdir -p diffXSecFromSignal/plots/combined/2012/binning
## mkdir -p diffXSecFromSignal/plots/combined/2012/effAndAcc
## mkdir -p diffXSecFromSignal/plots/combined/2012/genRecoCorrPlots
## mkdir -p diffXSecFromSignal/plots/combined/2012/kinFitPerformance
## mkdir -p diffXSecFromSignal/plots/combined/2012/uncertainties
## mkdir -p diffXSecFromSignal/plots/combined/2012/uncertaintyDistributions
## mkdir -p diffXSecFromSignal/plots/combined/2012/uncertaintyDistributionsOverview
## mkdir -p diffXSecFromSignal/plots/combined/2012/monitoring/withoutRatioPlots

## b) you don't need to copy root files needed for the Analysis 
##    the are loaded automatically from /afs/naf.desy.de/group/cms/scratch/tophh/
## c) when using the shell script for the very first time, do "chmod a+x topDiffXSecAnalysis2012.sh
## find final plots in ./diffXSecFromSignal/plots/ after running the analysis via ./topDiffXSecAnalysis2012.sh
## d) The decay channel, PS and level can be specified via 3 arguments, e.g.
##    ./topDiffXSecAnalysis2012.sh muon (or electron or combined) extrapolate (or visible) parton (or hadron)
##    if 0 arguments are given -> default values as specified below
##    if not 0 or 3 arguments -> Abort!!!

########################
## configure settings ##
########################

## default values

# lepton flavour in semi leptonic decay
# choose \"muon\" or \"electron\" or \"combined\"
decayChannel=\"combined\" 

## extrapolate xSec to full PS?
## extrapolate = true / false (default: true)
extrapolate=true

## use hadron level PS instead of parton level PS?
## hadron = true / false (default: false)
hadron=false

## Unfolding closure test -> use pseudo data
## closureTestSpecifier = \"\" / \"Up\" or \"Down\" for ttbar shape distortions
##                           / \"500\" or \"750\" for corresponding Zprime pseudo data  (default: "")
closureTestSpecifier=\"\"
  
## run combination based on event yield combination instead of 
## combinedEventYields = true / false (default: false)
## attention: this affects only bothDecayChannelsCombination.C 
##            and is automatically adjusted with input string "combined2"
combinedEventYields=false

## take arguments
clear
echo "-------------------------------------------------------------------------------------"
if [ $# -eq 0 ]; then
   echo "The default settings for decay channel, phase space extrapolation and level are used:"
elif [ $# -eq 3 -o $# -eq 4 ]; then
   ## set first argument
   if   [ $1 == "muon"      ];   then decayChannel=\"muon\"
   elif [ $1 == "electron"  ];   then decayChannel=\"electron\"
   elif [ $1 == "combined"  ];   then decayChannel=\"combined\"
   elif [ $1 == "combined2" ];   then 
       decayChannel=\"combined\"
       combinedEventYields=true
   else                          echo "1st argument ( $1 ) is not valid! Choose \"muon\", \"electron\" or \"combined\" as decay channel! Abort!"; exit
   fi
   ## set second argument
   if   [ $2 == "extrapolate" ]; then extrapolate=true
   elif [ $2 == "visible"     ]; then extrapolate=false
   else                          echo "2nd argument ( $2 ) is not valid! Choose \"extrapolate\" or \"visible\" as phase space! Abort!"; exit
   fi
   ## set third argument
   if   [ $3 == "parton" ];      then hadron=false
   elif [ $3 == "hadron" ];      then hadron=true
   else                          echo "3rd argument ( $3 ) is not valid! Choose \"parton\" or \"hadron\" as level! Abort!"; exit
   fi
## in case of fourth argument:
   if [ $# -eq 4 ]; then
       closureTestSpecifier=\"$4\"
   fi
   echo "Decay channel, phase space extrapolation and level were specified by user:"
else
   echo "Wrong number of arguments! Choose 0 arguments for default values and 3 arguments to specify" 
   echo "decay channel, phase space extrapolation and level. Abort!"
   exit
fi

echo "-------------------------------------------------------------------------------------"
echo "Decay channel:                              $decayChannel    "
echo "extrapolate:                                $extrapolate     "
echo "hadron:                                     $hadron          "
echo "closureTestSpecifier                        $closureTestSpecifier          "
echo "-------------------------------------------------------------------------------------"
echo

## folder on /afs/naf.desy.de/group/cms/scratch/tophh where MC and data files are stored
## inputFolderName=\"RecentAnalysisRun8TeV\" (default)
inputFolderName=\"RecentAnalysisRun8TeV\"

## Dataset and luminosity [/pb]
## has to fit to current dataset

mudataSample=\"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV/analyzeDiffXData2012ABCAllMuon.root\"
eldataSample=\"/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV/analyzeDiffXData2012ABCAllElec.root\"

if [ $decayChannel == \"electron\" ]; then
    dataLuminosity=12148
    dataSample=$eldataSample
else
    if [ $decayChannel == \"muon\" ]; then
	dataLuminosity=12148
	dataSample=$mudataSample
    else
	dataLuminosity=12148 # mean value
	dataSample=$eldataSample\":\"$mudataSample
    fi
fi

## Data label, required for filename
## dataLabel=2012 (default)
dataLabel=2012

## save all plots? (.eps and .root)
save=true

## detail level of output 
## 0: no output (default)
## 1: std output
## 2: output for debugging
verbose=0

## Re-create monitoring plots
## redoControlPlots = true / false (default: true)
redoControlPlots=true

## Re-create systematic plots
## redoSystematics = true / false (default: true)
redoSystematics=true

## Make pt plots logarithmic
## makeLogPlots = true / false (default: false)
makeLogPlots=false

## last systematic to proceed (0: only std analysis without variation)
## has to be consistent with the enumerator "systematicVariation" in "basicFunctions.h"
## maxSys>0 needs a lot of time (must be <= 48 (default), see list of systematics below)
maxSys=48

## Include cross-check variables to get additional differential cross-sections for
## a) pT(top) and pT(antitop)
## b) y(top)  and y(antitop)
## c) eta(l+) and eta(l-)
## Attention: The exectution mainly of analyzeHypothesisKinFit.C lasts longer if this parameter is set to true
##
## inclCCVars = true / false (default: false)
inclCCVars=false

## Shape variations:
## a) Calculate them at all
##    shapeVar = true / false (default: true) 
## b) Exclude them from total systematic uncertainty 
##    exclShapeVar = true / false (default: true)
shapeVar=true
exclShapeVar=true

## disable waiting time to read output
## fast = true / false (default: true)
fast=true

## delete all (old) existing .eps, .png and .pdf plots?
## clean = true / false (default: false)
clean=false

## use SVD unfolding?
## SVD = true / false (default: true)
SVD=true

## redetermine regularisation parameter tau?
## redetTau = true / false (default: false)
redetTau=false

## Re-create purity/stability/resolution plots
## redoPurStab = true / false (default: true)
redoPurStab=true

## Use bin-centre corrections (BCC)
## useBCC = true / false (default: true)
useBCC=true

#### =====================
####  Prepare running
#### =====================

BoolArray=( true false )

## start the timer to stop the time needed for running the whole analysis
START=$(date +%s)

# get proper name of rootfile
PS=""
LV="Parton"
if [ $extrapolate == false ] 
    then 
    PS="PhaseSpace"
    echo $PS
fi
if [ $extrapolate == false -a $hadron == true ] 
    then
    LV="Hadron"
fi

# switches for unfolding closure test
if [ $closureTestSpecifier != \"\" ] 
    then
    echo
    echo "CLOSURE TEST FOR UNFOLDING!"
    echo "closure test type: " $closureTestSpecifier
    redoControlPlots=false
    redoSystematics=false
    maxSys=0
    redoPurStab=false
    useBCC=false
fi

muonFile=./diffXSecTopSemiMu$dataLabel$LV$PS.root
elecFile=./diffXSecTopSemiElec$dataLabel$LV$PS.root
combFile=./diffXSecTopSemiLep$dataLabel$LV$PS.root

## print out configuration
if [ $decayChannel == \"combined\" ]
    then
    echo
    echo "combining the electron and muon channel"
    if [ -f $muonFile -a -f $elecFile -a $1 != "combined2" ]; then
	echo
	echo "Doing the full differential top xSec analysis. "
	echo 
	echo "Data Label:                                 $dataLabel       "
	echo "Used data:                                  $dataSample      "
	echo "Decay channel:                              $decayChannel    "
	echo "Luminosity:                                 $dataLuminosity  " 
	echo "Re-do control plots:                        $redoControlPlots" 
	echo "Re-do systematic uncertainties:             $redoSystematics "
	echo "Number of considered systematics:           $maxSys          "
	echo "Consider shape variation:                   $shapeVar        " 
	echo "Exclude shape unc. from total uncertainty:  $exclShapeVar    "
	echo "Save plots:                                 $save            " 
	echo
    elif [  $1 == "combined2" ]; then
	echo
	echo "Doing the full differential top xSec analysis. "
	echo "COMBINING THE CHANNELS AT EVENT YIELD LEVEL"
	echo 
	echo "Data Label:                                 $dataLabel       "
	echo "Used data:                                  $dataSample      "
	echo "Decay channel:                              $decayChannel    "
	echo "Luminosity:                                 $dataLuminosity  " 
	echo "Re-do control plots:                        $redoControlPlots" 
	echo "Re-do systematic uncertainties:             $redoSystematics "
	echo "Number of considered systematics:           $maxSys          "
	echo "Consider shape variation:                   $shapeVar        " 
	echo "Exclude shape unc. from total uncertainty:  $exclShapeVar    "
	echo "Save plots:                                 $save            " 
	echo
    else
	echo
	echo "NOTE: The combination requires two files"
	echo "a) "$muonFile
        echo "b) "$elecFile
	echo
	echo "Please get them by running the e/mu channel first"
	echo
	exit
    fi
fi

if [ $fast = false ]
    then
    sleep 5
fi

#### ===================================
####  Delete existing root file/ plots
#### ===================================
echo
echo "part A: Delete existing files and plots (if applicable)"
if [ $fast = false ]; then
    sleep 3
fi

if [ $clean = true ]; then

    if [ $redoSystematics = false ]; then
	echo
	echo "Flag 'redoSystematics' is set to $redoSystematics "
	echo "Flag 'clean' set to 'false' to avoid deleting files/plots which are not recreated "
	echo
	clean=false
    else
	
        ## ============================
        ##  Delete existing root files
        ## ============================
   
	echo
	echo "Part A1: delete existing root file"
	if [ $decayChannel == \"electron\" ]; then
	    echo $elecFile
	    rm $elecFile
	else
	    if [ $decayChannel == \"muon\" ]; then
		echo $muonFile
		rm $muonFile
	    else
		if [ $decayChannel == \"combined\" ]; then
		    echo $muonFile
		    rm $combFile		    
		fi
	    fi
	fi 
	
        ## ============================    
        ##  Delete existing plots
	## ============================
    
	echo
	echo "Part A2: delete existing plots within diffXSecFromSignal/plots/$decayChannel/$dataLabel/*/*.*"
	if [ $decayChannel == \"electron\" ]; then
	    rm ./diffXSecFromSignal/plots/electron/$dataLabel/*/*.*
	else
	    if [ $decayChannel == \"muon\" ]; then
		rm ./diffXSecFromSignal/plots/muon/$dataLabel/*/*.*
	    else
		if [ $decayChannel == \"combined\" ]; then
		    rm ./diffXSecFromSignal/plots/combined/$dataLabel/*/*.*
		fi
	    fi
	fi
    fi
fi

#### =====================
####  Run cut monitoring
#### =====================

BEFOREB=$(date +%s)
echo
echo "Part B: process cut monitoring macro"
if [ $fast = false ]
    then
    sleep 3
fi

if [ $redoControlPlots = true ]; then
    
    ## Compile library
    
    if [ -f commandsMonPrepare.cint ]; then    
	rm commandsMonPrepare.cint
	rm analyzeTopDiffXSecMonitoring_C.so
	rm analyzeTopDiffXSecMonitoring_C.d
    fi
    
    cat >> commandsMonPrepare.cint << EOF
.L analyzeTopDiffXSecMonitoring.C++
EOF
    
    root -l -b < commandsMonPrepare.cint
    
    ## Execute macro

    if [ $redoControlPlots = true ]; then
	
	for label in "${BoolArray[@]}"; do
	    
        ## label = 1: Control plots with ratio plots
        ## label = 0: Control plots without ratio plots
	    
	    if [ -f commandsMonRun.cint ]; then    
		rm commandsMonRun.cint       
	    fi    
	    
	    cat >> commandsMonRun.cint << EOF
.L analyzeTopDiffXSecMonitoring_C.so
analyzeTopDiffXSecMonitoring($dataLuminosity, $save, $verbose, $inputFolderName, $dataSample, $decayChannel, $label, $extrapolate, $hadron) 
EOF
	    echo ""
	    echo " Processing .... analyzeTopDiffXSecMonitoring.C++($dataLuminosity, $save, $verbose, $inputFolderName, $dataSample, $decayChannel, true, $extrapolate, $hadron)"
	    root -l -b < commandsMonRun.cint
	done
    fi
fi

### ===================================
###  Run migration macro for binning 
### ===================================

BEFOREC=$(date +%s)
echo
echo "Part C: process migration macro to validate binning"
if [ $fast = false ]
    then
    sleep 3
fi

if [ $redoPurStab = true -a $redoControlPlots = true ]
    then
    # Array of differential variables
    listVar_=( \"topPt\" \"topY\" \"ttbarPt\" \"ttbarY\" \"ttbarMass\" \"lepPt\" \"lepEta\" \"bqPt\" \"bqEta\")
    plotAcceptance=true
    if [ $hadron = true ]; 
	then 
	listVar_=( \"lepPt\" \"lepEta\" \"bqPt\" \"bqEta\")
	plotAcceptance=false
    fi
    	
    echo "purity and stability will be calculated for the following variables: "
    echo
    echo "${listVar_[@]}"
  
    # loop over all systematic variations
    for (( iVar=0; iVar<${#listVar_[@]}; iVar++ )); do
	root -l -q -b './purityStabilityEfficiency.C++('${listVar_[$iVar]}','$save', '$decayChannel', '$inputFolderName', '$plotAcceptance', true, false, 99999, 0, '$hadron')'
    done
fi



#### ============================
####  Prepare PDF uncertainties 
#### ============================
BEFOREPDF=$(date +%s)
echo
echo "Part PDF: Prepare files for pdf uncertainties"

if [ $decayChannel != \"combined\" -a $redoSystematics = true ]; then
    echo
    root -l -q -b './analyzeTopDiffXSecMCdependency.C++('$dataLuminosity','$decayChannel', '$save', '$verbose', '$inputFolderName', '$dataSample', 'true', '$inclCCVars')' 
elif [ $1 == "combined2" -a $redoSystematics = true ]; then
    root -l -q -b './analyzeTopDiffXSecMCdependency.C++('$dataLuminosity', '\"muon\"',     '$save', '$verbose', '$inputFolderName', '$mudataSample', 'true', '$inclCCVars')' 
    root -l -q -b './analyzeTopDiffXSecMCdependency.C++('$dataLuminosity', '\"electron\"', '$save', '$verbose', '$inputFolderName', '$eldataSample', 'true', '$inclCCVars')' 
else
    echo "Done for 2012 analysis (in e/mu channel separate or when combining event yields) and if systematics are requested to be re-done (redoSystematics set to $redoSystematics)."
fi

#### ======================================================================
####  Run shape distortion macro to get ROOT files for MC dependency 
#### ======================================================================
BEFORED=$(date +%s)
echo
echo "Part D: Create rootfiles with shape variations"
echo

if [ $shapeVar = true -a $redoSystematics = true ]; then
    
    if [ $decayChannel != \"combined\" ]; then
	
	echo "will be done"
	root -l -q -b './analyzeTopDiffXSecMCdependency.C++('$dataLuminosity','$decayChannel', '$save', '$verbose', '$inputFolderName', '$dataSample', 'false', '$inclCCVars')'
    elif [ $1 == "combined2" ]; then
	root -l -q -b './analyzeTopDiffXSecMCdependency.C++('$dataLuminosity', '\"muon\"',     '$save', '$verbose', '$inputFolderName', '$mudataSample', 'false', '$inclCCVars')'
	root -l -q -b './analyzeTopDiffXSecMCdependency.C++('$dataLuminosity', '\"electron\"', '$save', '$verbose', '$inputFolderName', '$eldataSample', 'false', '$inclCCVars')'
    else
	echo "only done for 2012 analysis in e/mu channel separate"
    fi
else
    echo "choose shapeVar = true and redoSystematics = true!"
fi

#### ==========================================
####  Run efficiency & cross section macro 
####  ==========================================
BEFOREE=$(date +%s)
echo
echo "part E1: process cross section calculation macro for all systematics"
echo "INFO: missing files must not be problematic"
echo "      either all WZ, WW and ZZ or the combined VV sample are necessary"
echo "      same is true for the single top samples (s, t, tW)"
echo
if [ $fast = false ]
    then
    sleep 5
fi

## print key for systematic variations
## has to be consistent with the enumerator "systematicVariation" in "basicFunctions.h"

echo

echo "  0: sysNo                                                      "
echo "  1: sysLumiUp                   2: sysLumiDown                 "
echo "  3: sysPUUp                     4: sysPUDown                   "
echo "  5: sysJESUp                    6: sysJESDown                  "
echo "  7: sysJERUp                    8: sysJERDown                  "
echo "  9: sysLepEffSFNormUp          10: sysLepEffSFNormDown         "
echo " 11: sysLepEffSFShapeEtaUp      12: sysLepEffSFShapeEtaDown     "
echo " 13: sysLepEffSFShapePtUp       14: sysLepEffSFShapePtDown      "
echo " 15: sysBtagSFUp                16: sysBtagSFDown               "
echo " 17: sysBtagSFShapePt65Up       18: sysBtagSFShapePt65Down      "
echo " 19: sysBtagSFShapeEta0p7Up     20: sysBtagSFShapeEta0p7Down    "
echo " 21: sysMisTagSFUp              22: sysMisTagSFDown             "
echo " 23: sysTopScaleUp              24: sysTopScaleDown             "
echo " 25: sysVBosonScaleUp           26: sysVBosonScaleDown          "
echo " 27: sysSingleTopScaleUp        28: sysSingleTopScaleDown       "
echo " 29: sysTopMatchUp              30: sysTopMatchDown             "
echo " 31: sysVBosonMatchUp           32: sysVBosonMatchDown          "
echo " 33: sysTopMassUp               34: sysTopMassDown              "
echo " 35: sysQCDUp                   36: sysQCDDown                  "
echo " 37: sysSTopUp                  38: sysSTopDown                 "
echo " 39: sysDiBosUp                 40: sysDiBosDown                "
echo " 41: sysPDFUp                   42: sysPDFDown                  "
echo " 43: sysHadUp                   44: sysHadDown                  "
echo " 45: sysGenMCatNLO              46: sysGenPowheg                "
echo " 47: sysShapeUp                 48: sysShapeDown                "
echo " 49: ENDOFSYSENUM                                               "

echo

if [ $fast = false ]; then
    sleep 5
fi

#### ============================
####  Systematic Uncertainties
#### ============================

BEFORESYS=$(date +%s)

if [ $decayChannel != \"combined\" -o $1 == "combined2" ]; then

    ## ====================================================================================
    ##  Compile library, required only once before processing systematic uncertainties
    ## ====================================================================================
    
    if [ -f commandsSysPrepare.cint ]; then    
	rm commandsSysPrepare.cint
	rm analyzeHypothesisKinFit_C.so
	rm analyzeHypothesisKinFit_C.d
    fi

    cat >> commandsSysPrepare.cint << EOF
.L analyzeHypothesisKinFit.C++g
EOF

    root -l -b < commandsSysPrepare.cint

    ## ================================================================================================================
    ##  Processing reference data (noSys), required always thus excluded from looping over systematic uncertainties 
    ## ================================================================================================================

    if [ -f commandsNoSysRun.cint ]; then    
	rm commandsNoSysRun.cint
    fi
    
    cat >> commandsNoSysRun.cint << EOF
.L analyzeHypothesisKinFit_C.so
analyzeHypothesisKinFit($dataLuminosity, $save, 0, $verbose, $inputFolderName, $dataSample, $decayChannel, $SVD, $extrapolate, $hadron, $inclCCVars, $redetTau, $closureTestSpecifier)
EOF

    echo ""
    echo " Processing .... analyzeHypothesisKinFit($dataLuminosity, $save, 0, $verbose, $inputFolderName, $dataSample, $decayChannel, $SVD, $extrapolate, $hadron ,$inclCCVars, $redetTau, $closureTestSpecifier)"
    root -l -b < commandsNoSysRun.cint


    ## ==========================================
    ##  Processing systematic uncertainties 
    ## ==========================================

    if [ $redoSystematics = true ]; then
    
        ## loop all systematic variations (excluding shape variations)
	
	for (( systematicVariation = 1; systematicVariation <= $maxSys;  systematicVariation++ )); do

            ## exclude shape variation


	    if [ $systematicVariation == 47 -o $systematicVariation == 48 ]; then
		echo " Shape variations are executed separately."
	    else
	    ## run macro for 2012 analysis
		
		if [ -f commandsSysRun.cint ]; then    
		    rm commandsSysRun.cint
		fi
		
		cat >> commandsSysRun.cint << EOF
.L analyzeHypothesisKinFit_C.so
analyzeHypothesisKinFit($dataLuminosity, $save, $systematicVariation, $verbose, $inputFolderName, $dataSample, $decayChannel, $SVD, $extrapolate, $hadron, $inclCCVars, $redetTau, $closureTestSpecifier)
EOF
		echo ""
		echo " Processing .... analyzeHypothesisKinFit($dataLuminosity, $save, $systematicVariation, $verbose, $inputFolderName, $dataSample, $decayChannel, $SVD, $extrapolate, $hadron, $inclCCVars, $redetTau, $closureTestSpecifier)"
		root -l -b < commandsSysRun.cint
	    fi  
	done
	
        ##  Processing shape variations
        
	if [ $shapeVar = true ]; then
	    
	    echo ""
	    echo " All regular systematic uncertainties processed .... Now running shape variations."
	    echo ""
	    
	    for (( systematicVariation = 47; systematicVariation <= 48;  systematicVariation++ )); do
		
		if [ -f commandsSysShapeVarRun.cint ]; then    
		    rm commandsSysShapeVarRun.cint
		fi
		
		cat >> commandsSysShapeVarRun.cint << EOF
.L analyzeHypothesisKinFit_C.so
analyzeHypothesisKinFit($dataLuminosity, $save, $systematicVariation, $verbose, $inputFolderName, $dataSample, $decayChannel, $SVD, $extrapolate, $hadron, $inclCCVars, $redetTau, $closureTestSpecifier)
EOF
		echo ""
		echo " Processing .... analyzeHypothesisKinFit($dataLuminosity, $save, $systematicVariation, $verbose, $inputFolderName, $dataSample, $decayChannel, $SVD, $extrapolate, $hadron, $inclCCVars, $redetTau, $closureTestSpecifier)"
		root -l -b < commandsSysShapeVarRun.cint
	    done
	fi
    fi
else
    echo "will be ignored, only done for decayChannel=muon/electron"
fi

#### ===================================
####  Combine electron and muon channel 
#### ===================================
echo
echo "Part E2: Combine electron and muon channel"
if [ $fast = false ]; then
    sleep 2
fi

if [ $decayChannel == \"combined\" ]; then
    
    echo "Cross sections for all systematic variations and decay channels"
    echo
    
    ## Compile library
    
    if [ -f commandsCombineChannelsPrepare.cint ]; then    
	rm commandsCombineChannelsPrepare.cint
	rm bothDecayChannelsCombination_C.so
	rm bothDecayChannelsCombination_C.d
    fi
    
    cat >> commandsCombineChannelsPrepare.cint << EOF
.L bothDecayChannelsCombination.C++g
EOF
    
    root -l -b < commandsCombineChannelsPrepare.cint
    
    ## Execute macro

    if [ -f commandsCombineChannelsRun.cint ]; then    
	rm commandsCombineChannelsRun.cint       
    fi    
    
    cat >> commandsCombineChannelsRun.cint << EOF
.L bothDecayChannelsCombination_C.so
bothDecayChannelsCombination($dataLuminosity, $save, $verbose, $inputFolderName, $makeLogPlots, $extrapolate, $hadron, $inclCCVars, $combinedEventYields, $closureTestSpecifier)
EOF
    echo ""
    echo " Processing .... bothDecayChannelsCombination($dataLuminosity, $save, $verbose, $inputFolderName, $makeLogPlots, $extrapolate, $hadron, $inclCCVars, $combinedEventYields, $closureTestSpecifier)"
    root -l -b < commandsCombineChannelsRun.cint

else
    echo "will be ignored, only done for decayChannel=combined"
fi
echo

AFTERSYS=$(date +%s)

#### ==========================================
####  Combine uncertainties for final xSecs 
#### ==========================================

BEFOREF=$(date +%s)

echo
echo "Part F: Calculate systematic errors and draw final cross section"
if [ $fast = false ]; then
    sleep 3
fi

## FIXME: need to restore compatibility with root 532, skip it 
exit 0
## and do it by hand using the older rootversion in the meantime via
#ini root526

## Compile library

if [ -f commandsCombineUncPrepare.cint ]; then    
  
    rm commandsCombineUncPrepare.cint
    rm combineTopDiffXSecUncertainties_C.so
    rm combineTopDiffXSecUncertainties_C.d
fi

cat >> commandsCombineUncPrepare.cint << EOF
.L BCC.C++
.L combineTopDiffXSecUncertainties.C++
EOF

root -l -b < commandsCombineUncPrepare.cint

## Execute macro

if [ -f commandsCombineUncRun.cint ]; then    
    rm commandsCombineUncRun.cint
fi
    
cat >> commandsCombineUncRun.cint << EOF
.L BCC_C.so
.L combineTopDiffXSecUncertainties_C.so
combineTopDiffXSecUncertainties($dataLuminosity, $save, $verbose, $inputFolderName, $decayChannel, $exclShapeVar, $extrapolate, $hadron, $inclCCVars, $closureTestSpecifier, $useBCC)
EOF
    
echo ""
echo " Processing .... combineTopDiffXSecUncertainties($dataLuminosity, $save, $verbose, $inputFolderName, $decayChannel, $exclShapeVar, $extrapolate, $hadron, $inclCCVars, $closureTestSpecifier, $useBCC)"
root -l -b < commandsCombineUncRun.cint

## FIXME: need to restore compatibility with root 532, use older version in the meantime
#ini -d root526

#### ==========================================
####  Create ratio plots for final xSecs 
#### ==========================================

if [ $decayChannel == \"combined\" -a $closureTestSpecifier == \"\" ]; then
    echo ""
    echo " Processing .... createTheoryDataRatios($extrapolate, $hadron, $verbose)"
    root -l -q -b './createTheoryDataRatios.C++('$extrapolate', '$hadron', '$verbose')'
fi


#### ===================================================
####  Create latex code result tables for final xSecs 
#### ===================================================
if [ $decayChannel == \"combined\" -a $closureTestSpecifier == \"\" ]; then
    echo ""
    echo " Processing .... makeResultTables($decayChannel, $extrapolate, $hadron, $inclCCVars)"
    root -l -q -b './makeResultTables.C++('$decayChannel', '$extrapolate', '$hadron', '$inclCCVars')'
fi

#### ==========================================

echo ""
echo " All analysis steps finished!"
echo ""

#### ==========================================

#### ============================
####  After running the analysis
#### ============================
## stop the timer and echo time
END=$(date +%s)
TIME=$(( $END - $START ))
echo "time needed: $TIME seconds"
if [ $maxSys -ge 1 ]
    then
    SYS=$(( $AFTERSYS - $BEFORESYS ))
    echo "($SYS seconds due to systematic variations)"
fi
echo "part A: $(( $BEFOREB   - $START    )) seconds (clean up  )"
echo "part B: $(( $BEFOREC   - $BEFOREB  )) seconds (monitoring)"
echo "part C: $(( $BEFOREPDF - $BEFOREC  )) seconds (migration)"
echo "part D: $(( $BEFORED   - $BEFOREPDF)) seconds (prepare PDF uncertainty run)"
echo "part E: $(( $BEFOREE   - $BEFORED  )) seconds (shape variations)"
echo "part F: $(( $BEFOREF   - $BEFOREE  )) seconds (xSec, $maxSys systematic variations considered)"
echo "part G: $(( $END       - $BEFOREF  )) seconds (errors and final xSec)"
