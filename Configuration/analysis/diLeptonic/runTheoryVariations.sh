#!/bin/sh

source parallelTools.sh

#start the longest job first
for syst in FullLeptMadgraphWithSpinCorrelation SemiLeptMadgraphWithSpinCorrelation HadronicMadgraphWithSpinCorrelation; do
    for ch in ee emu mumu; do
        $LA -f $syst -c $ch &
        w
    done
done

#these have rather low statistics, don't need to split them into ee/emu/mumu
for Syst in match mass scale powheg mcatnlo; do
    $LA -f $Syst &
    w
done

wait

echo "Done!"

