#!/bin/sh

##################################################
# look for argument on command line
##################################################
option="ps"
if [ $1 ]
    then
    if [ $1 = "ps" -o $1 = "eps" ]
	then
	option=$1
    else
	if [ $1 = "pdf" ]
	    then
	    option="eps"
	else
	    echo "You can choose one of the following options"
	    echo "and use them as argument to this script:"
	    echo "+ ps  (default; a single output file)"
	    echo "+ eps (one file per histogram)"
	    echo "+ pdf (one file per histogram)"
	    echo "Please try again..."
	    exit 0
	fi
    fi
    echo "Option chosen for output format: $1"
else
    option="ps"
    echo "Option chosen for output format: ps"
fi



##################################################
# choose steering file to use
##################################################
cfgFile="TopAnalysis/TopFilter/test/makeSemiLepMuonSelectionAnalysisPlot_cutflow.cfg"

if [ ! -e $cfgFile ]
    then
    echo "ERROR! The following file is missing:"
    echo "... $cfgFile"
    exit 0
fi

echo "The following steering file will be used:"
echo "... $cfgFile"

##################################################
# choose hist files to use
##################################################
muonHist="analyzeMuon.hist"
jetsHist="analyzeJets.hist"

if [ ! -e $muonHist ]
    then
    echo "ERROR! The following file is missing:"
    echo "... $muonHist"
    exit 0
fi
if [ ! -e $jetsHist ]
    then
    echo "ERROR! The following file is missing:"
    echo "... $jetsHist"
    exit 0
fi

echo "The following hist files will be used:"
echo "... $muonHist"
echo "... $jetsHist"

##################################################
# choose output directory
##################################################
outDir="inspect_SemiLepCutflow"

echo "Output file(s) will be written to the following directory"
echo "(directory will be created if it is not yet existing):"
echo "... $outDir"

##################################################
# ask for confirmation before starting to work
##################################################
status=0
while [ $status = 0 ]
  do
  echo -n "Start job now? [y/n] "
  read answer
  case $answer in
      y|Y) status=1;;
      n|N) exit 0;;
      *) echo "Illegal answer, please try again...";;
  esac
done
echo "-------------------------------------------------------"

##################################################
# make backup of the original steering file
##################################################
cp $cfgFile $cfgFile.tmp

##################################################
# set output directory
##################################################
if [ ! -e $outDir ]
    then
    echo "Creating output directory..."
    mkdir $outDir
fi
 
TEXT1="writePlotsTo.*\$"
TEXT2="writePlotsTo = $outDir"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile

##################################################
# set output format option
##################################################
TEXT1="writePlotsAs.*\$"
TEXT2="writePlotsAs = $option"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile

##################################################
# 1: inspect before muon eta cut
##################################################
echo "Going to inspect before muon eta cut..."
# choose hist file
TEXT1="histInput.*\$"
TEXT2="histInput = $muonHist"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
# choose input directory
TEXT1="inputDirs.*\$"
TEXT2="inputDirs = semiMuonNon"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
# choose histogram
TEXT1="histFilter.*\$"
TEXT2="histFilter = kin_eta"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
Inspect_top $cfgFile
if [ $option = "ps" -a -e $outDir/inspect.ps ]
    then
    mv $outDir/inspect.ps $outDir/inspect1.tmp.ps
fi

##################################################
# 2: inspect before muon pt cut
##################################################
echo "Going to inspect before muon pt cut..."
# choose hist file
TEXT1="histInput.*\$"
TEXT2="histInput = $muonHist"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
# choose input directory
TEXT1="inputDirs.*\$"
TEXT2="inputDirs = semiMuonMuEta"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
# choose histogram
TEXT1="histFilter.*\$"
TEXT2="histFilter = kin_pt"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
Inspect_top $cfgFile
if [ $option = "ps" -a -e $outDir/inspect.ps ]
    then
    mv $outDir/inspect.ps $outDir/inspect2.tmp.ps
fi

##################################################
# 3: inspect before muon trk iso
##################################################
echo "Going to inspect before muon trk iso cut..."
# choose hist file
TEXT1="histInput.*\$"
TEXT2="histInput = $muonHist"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
# choose input directory
TEXT1="inputDirs.*\$"
TEXT2="inputDirs = semiMuonMuPt"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
# choose histogram
TEXT1="histFilter.*\$"
TEXT2="histFilter = iso_isoTrk"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
Inspect_top $cfgFile
if [ $option = "ps" -a -e $outDir/inspect.ps ]
    then
    mv $outDir/inspect.ps $outDir/inspect3.tmp.ps
fi

##################################################
# 4: inspect before muon cal iso
##################################################
echo "Going to inspect before muon cal iso cut..."
# choose hist file
TEXT1="histInput.*\$"
TEXT2="histInput = $muonHist"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
# choose input directory
TEXT1="inputDirs.*\$"
TEXT2="inputDirs = semiMuonTrkIso"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
# choose histogram
TEXT1="histFilter.*\$"
TEXT2="histFilter = iso_isoCal"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
Inspect_top $cfgFile
if [ $option = "ps" -a -e $outDir/inspect.ps ]
    then
    mv $outDir/inspect.ps $outDir/inspect4.tmp.ps
fi

##################################################
# 5: inspect before muon jet iso
##################################################
echo "Going to inspect before muon jet iso cut..."
# choose hist file
TEXT1="histInput.*\$"
TEXT2="histInput = $muonHist"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
# choose input directory
TEXT1="inputDirs.*\$"
TEXT2="inputDirs = semiMuonCalIso"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
# choose histogram
TEXT1="histFilter.*\$"
TEXT2="histFilter = iso_isoJet"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
Inspect_top $cfgFile
if [ $option = "ps" -a -e $outDir/inspect.ps ]
    then
    mv $outDir/inspect.ps $outDir/inspect5.tmp.ps
fi

##################################################
# 6: inspect before jet eta cuts
##################################################
echo "Going to inspect before jet eta cuts..."
# choose hist file
TEXT1="histInput.*\$"
TEXT2="histInput = $jetsHist"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
# choose input directory
TEXT1="inputDirs.*\$"
TEXT2="inputDirs = semiJetsJetIso"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
# choose histogram
TEXT1="histFilter.*\$"
TEXT2="histFilter = ta"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
TEXT1="filterOption.*\$"
TEXT2="filterOption = contains"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
Inspect_top $cfgFile
if [ $option = "ps" -a -e $outDir/inspect.ps ]
    then
    mv $outDir/inspect.ps $outDir/inspect6.tmp.ps
fi

##################################################
# 7: inspect before jet pt cuts
##################################################
echo "Going to inspect before jet pt cuts..."
# choose hist file
TEXT1="histInput.*\$"
TEXT2="histInput = $jetsHist"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
# choose input directory
TEXT1="inputDirs.*\$"
TEXT2="inputDirs = semiJetsJetEta"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
# choose histogram
TEXT1="histFilter.*\$"
TEXT2="histFilter = kin_n"
sed -i "s%$TEXT1%$TEXT2%" $cfgFile
Inspect_top $cfgFile
if [ $option = "ps" -a -e $outDir/inspect.ps ]
    then
    mv $outDir/inspect.ps $outDir/inspect7.tmp.ps
fi

echo "-------------------------------------------------------"
##################################################
# concatenate ps files if option = "ps"
##################################################
if [ $option = "ps" ]
    then
    echo "Concatenating all ps files from inspect into a single inspect.ps ..."
    gs -q -dNOPAUSE -dSAFER -sOutputFile=$outDir/inspect.ps -sDEVICE=pswrite \
	$outDir/inspect1.tmp.ps $outDir/inspect2.tmp.ps $outDir/inspect3.tmp.ps $outDir/inspect4.tmp.ps \
	$outDir/inspect5.tmp.ps $outDir/inspect6.tmp.ps $outDir/inspect7.tmp.ps \
	quit.ps
    rm -f $outDir/inspect*.tmp.ps
    echo "-------------------------------------------------------"
else
##################################################
# convert eps files to pdf files if argument was "pdf"
##################################################
    if [ $1 = "pdf" ]
	then
	echo "Converting eps files into pdf files..."
	for FILE in `ls $outDir/*.eps`
	  do
	  epstopdf $FILE
	done
	status=0
	while [ $status = 0 ]
	  do
	  echo -n "Remove all eps files in $outDir ? [y/n] "
	  read answer
	  case $answer in
	      y|Y) status=1; echo "Removing eps files..."; rm -f $outDir/*.eps;;
	      n|N) status=1;;
	      *) echo "Illegal answer, please try again...";;
	  esac
	done
	echo "-------------------------------------------------------"
    fi
fi

##################################################
# restore the original steering file
##################################################
mv $cfgFile.tmp $cfgFile

echo "Done."
