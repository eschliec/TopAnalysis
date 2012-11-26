#!/bin/sh

# A relative or absolute directory can be given as argument.
# In this directory (or the current directory if non is given),
# the trigger reports and root files are added using
# sumTriggerReports2.pl and fasthadd.pl, respectively.
# In addition, the root files (including *.root.NNN) are deleted.
# The number of cores (nCores) that is used by fasthadd.pl
# can be changed below. Please, only use Integers!

nCores=16

if [ ${1} ]; then
    if [ -d ${1} ]; then
	cd ${1}
    else
	echo Directory \"$1\" does not exist! Aborting script!
	exit 1
    fi
fi

for fileName in `ls *.root.*`; do
    
    newFileName=`echo ${fileName} | sed -r 's/root.[0-9]+//'`

    echo "sumTriggerReports2.pl out*.txt >| ${newFileName}txt"
    eval `echo "sumTriggerReports2.pl out*.txt >| ${newFileName}txt"`

    echo "fasthadd.pl -j ${nCores} ${newFileName}root *.root"
    eval `echo "fasthadd.pl -j ${nCores} ${newFileName}root *.root"`

    if [ -e ${newFileName}root ]; then

	echo "rm ${fileName}"
	rm ${fileName}

	echo "rm analyzeMuon*.root"
	for rootFile in `ls analyzeMuon*.root`; do

	    rm ${rootFile}

	done

    fi

done
