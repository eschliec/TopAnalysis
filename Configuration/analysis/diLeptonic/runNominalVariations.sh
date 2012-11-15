#!/bin/sh

for sys in PU_UP PU_DOWN TRIG_UP TRIG_DOWN BTAG_UP BTAG_DOWN BTAG_PT_UP BTAG_PT_DOWN BTAG_ETA_UP BTAG_ETA_DOWN; do
    ./load_Analysis -s $sys -f Nominal &
done
wait

