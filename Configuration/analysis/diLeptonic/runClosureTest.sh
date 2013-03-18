#!/bin/sh

for i in ee emu mumu; do 
    for s in 0.08 0.16 0.24 0.32 -0.08 -0.16 -0.24 -0.32; do 
        ./load_Analysis -f tau.root -s closure --closure ytop --slope $s -c $i &
    done
done
wait

for i in ee emu mumu; do 
    for s in 0.002 0.004 0.006 0.008 -0.002 -0.004 -0.006 -0.008; do 
        ./load_Analysis -f tau.root -s closure --closure pttop --slope $s -c $i &
    done
done
wait

echo "Done!"

echo "Now run unfoldClosure.sh"

