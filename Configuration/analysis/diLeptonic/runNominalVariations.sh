#!/bin/sh

w() { while [ `ps ax | grep load_Analysis | wc -l` -gt 10 ]; do sleep 1; done }

# be careful when using this script, you will submit a lot of jobs at the same. You may slow down the work group server performance for a while

for sys in PU_UP PU_DOWN \
           TRIG_UP TRIG_DOWN \
           LEPT_UP LEPT_DOWN \
           BTAG_UP BTAG_DOWN BTAG_LJET_UP BTAG_LJET_DOWN \
           BTAG_PT_UP BTAG_PT_DOWN BTAG_ETA_UP BTAG_ETA_DOWN \
           BTAG_LJET_PT_UP BTAG_LJET_PT_DOWN BTAG_LJET_ETA_UP BTAG_LJET_ETA_DOWN \
           BTAG_BEFF_UP BTAG_BEFF_DOWN BTAG_CEFF_UP BTAG_CEFF_DOWN BTAG_LEFF_UP BTAG_LEFF_DOWN; do

    for c in ee emu mumu; do
        w
        ./load_Analysis -f dy -d 11 -c $c -s $sys&
        ./load_Analysis -f dy -d 13 -c $c -s $sys&
        ./load_Analysis -f dy -d 15 -c $c -s $sys&
        ./load_Analysis -f ttbarsignalplustau.root -c $c -s $sys&
    done

    for i in qcd single ttbarbg wtol ww wz zz; do
        w
        ./load_Analysis -f $i -s $sys&
    done

done
wait

