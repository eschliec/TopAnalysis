#!/bin/sh

for c in ee emu mumu; do
    ./load_Analysis -f dy -d 11 -c $c &
    ./load_Analysis -f dy -d 13 -c $c &
    ./load_Analysis -f dy -d 15 -c $c &
    ./load_Analysis -f ttbarsignal -c $c &
done

wait

for c in ee emu mumu; do
    ./load_Analysis -f ${c}_run2012A -c $c &
    ./load_Analysis -f ${c}_run2012B -c $c &
    ./load_Analysis -f ${c}_run2012C -c $c &
done 

for i in qcd single ttbarbg wtol ww wz zz; do
    ./load_Analysis -f $i &
done

wait

echo "Processing all nominal samples finished!"
