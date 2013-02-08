#!/bin/sh

w() { while [ `ps ax | grep load_Analysis | wc -l` -gt 10 ]; do sleep 1; done }

for Syst in match mass scale \
            powheg mcatnlo SpinCorrelation; do
    w
    ./load_Analysis -f $Syst &
done
wait
echo "Done!"
