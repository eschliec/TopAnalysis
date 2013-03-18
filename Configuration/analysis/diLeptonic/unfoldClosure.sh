#!/bin/sh

set +e
rm -r Plots
rm -r PlotsPt*
set -e

for slope in 0.002 0.004 0.006 0.008 -0.002 -0.004 -0.006 -0.008; do

   for ch in ee emu mumu; do
       rm FileLists/HistoFileList_Nominal_${ch}.txt
       echo "selectionRoot/closure/$ch/${ch}_ttbarsignalplustau_fakerun_pttop${slope}.root" > FileLists/HistoFileList_Nominal_${ch}.txt
       cat closure_${ch}mc.txt >> FileLists/HistoFileList_Nominal_${ch}.txt
   done

   ./Histo -t cp -p +hyptoppt
   ./Histo -t unfold -p +hyptoppt
   grep -v Bin < Plots/combined/ToppTResultsLaTeXAfter.txt | awk '{print $7}' > closure_result_pt_${slope}.txt

   mv Plots PlotsPt$slope
done

echo "Pt resultsPt"
for s in 0.002 0.004 0.006 0.008 -0.002 -0.004 -0.006 -0.008; do echo -n "1.
2.
3.
4.
5." | paste - trueMadPt.txt trueMadPt.txt closure_result_pt_$s.txt ; done





set +e
rm -r Plots
rm -r PlotsRap*
set -e

for slope in 0.080 0.160 0.240 0.320 -0.080 -0.160 -0.240 -0.320; do

    for ch in ee emu mumu; do
        rm FileLists/HistoFileList_Nominal_${ch}.txt
        echo "selectionRoot/closure/$ch/${ch}_ttbarsignalplustau_fakerun_ytop${slope}.root" > FileLists/HistoFileList_Nominal_${ch}.txt
        cat closure_${ch}mc.txt >> FileLists/HistoFileList_Nominal_${ch}.txt
    done

    ./Histo -t cp -p +hyptoprapidity
    ./Histo -t unfold -p +hyptoprapidity
    grep -v Bin < Plots/combined/TopRapidityResultsLaTeXAfter.txt | awk '{print $7}' > closure_result_y_${slope}.txt
    
    mv Plots PlotsRap$slope
done


#paste closure_result_y_-0.008.txt closure_result_y_-0.006.txt closure_result_y_-0.004.txt closure_result_y_-0.002.txt closure_result_y_0.002.txt closure_result_y_0.004.txt closure_result_y_0.006.txt closure_result_y_0.008.txt  >| closuretesty.csv
echo "Y resultsY"
for s in 0.080 0.160 0.240 0.320 -0.080 -0.160 -0.240 -0.320; do echo -n "1.
2.
3.
4.
5.
6.
7.
8." | paste - trueMadY.txt trueMadY.txt closure_result_y_$s.txt ; done



