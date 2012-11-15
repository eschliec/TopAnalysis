#!/bin/sh

for sys in PU_UP PU_DOWN TRIG_UP TRIG_DOWN BTAG_UP BTAG_DOWN; do
    ./load_Analysis -s $sys -f Nominal &
done
wait

