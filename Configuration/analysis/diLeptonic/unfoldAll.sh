#!/bin/sh

source parallelTools.sh

excludeList='HypNeutrinopT'
unfoldList=`awk '{print $1}' HistoList | grep Hyp | grep -Ev $excludeList`

# echo ""
# echo "*************** Information ******************"
# echo "Please unfold HypTTBarMass manually, is uses too much memory to be run in parallel or on the batch system"
# echo "i.e. run   ./Histo -t unfold -s all -p +HypTTBarMass"
# echo "**********************************************"
# echo ""
echo "Please press any key to start unfolding the following distributions in parallel or press Ctrl-C to cancel:"
echo "$unfoldList" | perl -l40 -pe ''
read -n 1 -s
echo ""

for i in $unfoldList; do 
    echo -n "Unfolding $i - "
    $HISTO -t unfold -s all -p +$i &
    w
done
