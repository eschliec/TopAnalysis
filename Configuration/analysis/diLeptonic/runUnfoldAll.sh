#!/bin/sh

for i in `awk '{print $1}' < HistoList`;
do 
    ./Histo -t unfold -p "+$i" & 
done
wait

echo "Unfolded everything!"

