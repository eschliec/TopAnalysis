#!/bin/sh

w() { while [ `ps ax | grep load_Analysis | wc -l` -gt 12 ]; do sleep 1; done }


for ch in emu mumu ee; do
    for no in `seq 0 44`; do # don't care if maximum > number of pdf sets! It will just do nothing.
        w
        ./load_Analysis --pdf $no -c $ch -f ttbarsignalplustau.root &
    done
done

wait
echo "Done!"

