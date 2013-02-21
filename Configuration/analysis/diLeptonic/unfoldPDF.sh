#!/bin/sh


runSpecificVariation() {
    variation="$1"
    echo "*********************************************************"
    echo "* starting to unfold variation $variation"
    echo "*********************************************************"
    for ch in ee emu mumu; do
        rm FileLists/HistoFileList_Nominal_$ch.txt
        cp HistoFileList_Nominal_$ch.txt FileLists/
        echo "selectionRoot/$variation/$ch/${ch}_ttbarsignalplustau.root" >> FileLists/HistoFileList_Nominal_$ch.txt
    done
    mv -f Plots Plots_temp
    # calculate inclusive xsection
    if [ -d "Plots_temp/$variation" ] ; then mv "Plots_temp/${variation}" Plots ; fi
    mkdir -p Plots/combined
    ./Histo -t cp -p XSec -s Nominal
    # now calculate differential distributions
    #for plot in `awk '{print $1}' HistoList | grep Hyp`; do
    for plot in HypTTBarMass; do
        ./Histo -t unfold -p +$plot -s Nominal &
    done
    wait
    mv -f Plots "Plots_temp/${variation}"
    mv -f Plots_temp Plots
}


./mk_HistoFileList.sh

for i in ee emu mumu combined; do
    grep -v ttbarsignalplustau.root < FileLists/HistoFileList_Nominal_$i.txt >| HistoFileList_Nominal_$i.txt 
done

runSpecificVariation PDF_CENTRAL
for no in `seq 1 22`; do
#for no in 1; do
    for var in UP DOWN; do
        runSpecificVariation "PDF_${no}_${var}"
    done
done

./mk_HistoFileList.sh

echo "Done"

