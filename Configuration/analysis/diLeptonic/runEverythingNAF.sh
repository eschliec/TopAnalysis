#!/bin/sh

LA="qsub -@ optionsLA.txt ./load_Analysis"

#for c in ee emu mumu; do
#    $LA -f dy -d 11 -c $c &
#    $LA -f dy -d 13 -c $c &
#    $LA -f dy -d 15 -c $c &
#    $LA -f ttbarsignalplustau.root -c $c &
#done
#exit
for c in ee emu mumu; do
    $LA -f ${c}_run2012A -c $c &
    $LA -f ${c}_run2012B -c $c &
    $LA -f ${c}_run2012C -c $c &
done 

for i in zz qcd single ttbarbg wtol ww wz ; do
    $LA -f $i -c ee&
    $LA -f $i -c emu&
    $LA -f $i -c mumu&
done

wait

for sys in JES_UP JES_DOWN JER_UP JER_DOWN \
           PU_UP PU_DOWN \
           TRIG_UP TRIG_DOWN \
           LEPT_UP LEPT_DOWN \
           KIN_UP KIN_DOWN \
           BTAG_UP BTAG_DOWN BTAG_LJET_UP BTAG_LJET_DOWN \
           BTAG_PT_UP BTAG_PT_DOWN BTAG_ETA_UP BTAG_ETA_DOWN \
           BTAG_LJET_PT_UP BTAG_LJET_PT_DOWN BTAG_LJET_ETA_UP BTAG_LJET_ETA_DOWN; do # \
#            BTAG_BEFF_UP BTAG_BEFF_DOWN BTAG_CEFF_UP BTAG_CEFF_DOWN BTAG_LEFF_UP BTAG_LEFF_DOWN; do

    for c in ee emu mumu; do
        $LA -f dy -d 11 -c $c -s $sys&
        $LA -f dy -d 13 -c $c -s $sys&
        $LA -f dy -d 15 -c $c -s $sys&
        $LA -f ttbarsignalplustau.root -c $c -s $sys&
    done
    wait
    for i in qcd single ttbarbg wtol ww wz zz; do
        $LA -f $i -s $sys&
    done
    wait
done

for ch in emu mumu ee; do
    for no in `seq 0 44`; do # don't care if maximum > number of pdf sets! It will just do nothing.
        $LA --pdf $no -c $ch -f ttbarsignalplustau.root &
    done
    wait
done

for ch in emu mumu ee; do
    for Syst in match mass scale \
                powheg mcatnlo SpinCorrelation; do
        $LA -f $Syst -c $ch &
        wait
    done
done
