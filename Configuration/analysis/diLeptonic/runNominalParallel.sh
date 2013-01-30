#!/bin/sh

#
#wbehrenh 18804 81.3  0.1 189564 126852 pts/7   R+   18:10  34:19 ./load_Analysis -f dy -d 11 -c ee
#wbehrenh 18813 77.9  0.1 189564 126852 pts/7   R+   18:10  32:52 ./load_Analysis -f dy -d 13 -c mumu
#wbehrenh 18817 77.8  0.1 164752 101888 pts/7   R+   18:10  32:52 ./load_Analysis -f ee_run2012B -c ee
#wbehrenh 18818 80.8  0.1 166664 103900 pts/7   R+   18:10  34:07 ./load_Analysis -f ee_run2012C -c ee
#wbehrenh 18823 76.6  0.1 166600 103752 pts/7   R+   18:10  32:21 ./load_Analysis -f mumu_run2012B -c mumu
#wbehrenh 18824 77.8  0.1 169824 107136 pts/7   R+   18:10  32:50 ./load_Analysis -f mumu_run2012C -c mumu
#

w() { while [ `ps ax | grep load_Analysis | wc -l` -gt 10 ]; do sleep 1; done }

for c in ee emu mumu; do
    ./load_Analysis -f dy -d 11 -c $c &
    ./load_Analysis -f dy -d 13 -c $c &
    ./load_Analysis -f dy -d 15 -c $c &
    ./load_Analysis -f ttbarsignalplustau.root -c $c &
done

#wait


for c in ee emu mumu; do
    w
    ./load_Analysis -f ${c}_run2012A -c $c &
    ./load_Analysis -f ${c}_run2012B -c $c &
    ./load_Analysis -f ${c}_run2012C -c $c &
done 

for i in qcd single ttbarbg wtol ww wz zz; do
    w
    ./load_Analysis -f $i &
done

wait

echo "Processing all nominal samples finished!"

