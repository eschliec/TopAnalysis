#!/bin/sh

for i in ee emu mumu combined; do
    grep -v ttbarsignalplustau.root < FileLists/HistoFileList_Nominal_$i.txt >| HistoFileList_Nominal_$i.txt 
done

for no in `seq 1 22`; do
#for no in 1; do
    for var in UP DOWN; do
        for ch in ee emu mumu; do
            rm FileLists/HistoFileList_Nominal_$ch.txt
            cp HistoFileList_Nominal_$ch.txt FileLists/
            echo "selectionRoot/PDF_${no}_${var}/$ch/${ch}_ttbarsignalplustau.root" >> FileLists/HistoFileList_Nominal_$ch.txt
        done
        mv -f Plots Plots_temp
        # calculate inclusive xsection
        mkdir -p Plots/combined
        ./Histo -t cp -p XSec -s Nominal
        # now calculate differential distributions
        #for plot in HypToppT; do
            #./Histo -t unfold -p +$plot -s Nominal
        #done
        if [ -d "Plots_temp/PDF_${no}_${var}" ] ; then rm -rf "Plots_temp/PDF_${no}_${var}" ; fi
        mv -f Plots Plots_temp/PDF_${no}_${var}
        mv -f Plots_temp Plots
    done
done

