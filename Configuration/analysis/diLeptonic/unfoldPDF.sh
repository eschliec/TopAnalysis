#!/bin/sh


runSpecificVariation() {
    variation="$1"
    for ch in ee emu mumu; do
        rm FileLists/HistoFileList_Nominal_$ch.txt
        cp HistoFileList_Nominal_$ch.txt FileLists/
        echo "selectionRoot/$variation/$ch/${ch}_ttbarsignalplustau.root" >> FileLists/HistoFileList_Nominal_$ch.txt
    done
    mv -f Plots Plots_temp
    # calculate inclusive xsection
    mkdir -p Plots/combined
    ./Histo -t cp -p XSec -s Nominal
    # now calculate differential distributions
    #for plot in HypToppT; do
        #./Histo -t unfold -p +$plot -s Nominal
    #done
    if [ -d "Plots_temp/$variation" ] ; then rm -rf "Plots_temp/${variation}" ; fi
    mv -f Plots "Plots_temp/${variation}"
    mv -f Plots_temp Plots
}


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

