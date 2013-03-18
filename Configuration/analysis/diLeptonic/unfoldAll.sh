#!/bin/sh

source parallelTools.sh

# in the excludeList, put all distributions that you dont want. Separate them with a |, i.e. HypLLBarDPhi|HypNeutrinopT
excludeList='HypNeutrinopT'
unfoldList=`awk '{print $1}' HistoList | grep Hyp | grep -Ev $excludeList`

echo "Please press any key to start unfolding the following distributions in parallel or press Ctrl-C to cancel:"
echo "$unfoldList" | perl -l40 -pe ''
read -n 1 -s
echo ""

for i in $unfoldList; do 
    echo -n "Unfolding $i - "
    $HISTO -t unfold -s all -p +$i &
    w
done
