#!/bin/sh

#
# ./load_Analysis -f dy -d 11 -c ee
# ./load_Analysis -f dy -d 13 -c mumu
# ./load_Analysis -f ee_run2012B -c ee
# ./load_Analysis -f ee_run2012C -c ee
# ./load_Analysis -f mumu_run2012B -c mumu
# ./load_Analysis -f mumu_run2012C -c mumu
#

w() { while [ `ps ax | grep load_Analysis | wc -l` -gt 10 ]; do sleep 1; done }

for c in ee emu mumu; do
    ./load_Analysis -f dy -d 11 -c $c &
    ./load_Analysis -f dy -d 13 -c $c &
    ./load_Analysis -f dy -d 15 -c $c &
    ./load_Analysis -f ttbarsignalplustau.root -c $c &
done

for c in ee emu mumu; do
    ./load_Analysis -f ${c}_run2012A -c $c &
    ./load_Analysis -f ${c}_run2012B -c $c &
    ./load_Analysis -f ${c}_run2012C -c $c &
done 

for i in zz qcd single ttbarbg wtol ww wz ; do
    w
    ./load_Analysis -f $i -c ee&
    ./load_Analysis -f $i -c emu&
    ./load_Analysis -f $i -c mumu&
done

wait

echo "Processing all nominal samples finished!"

