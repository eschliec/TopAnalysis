#!/bin/sh

merged=''

if [ -n "$1" ]; then
    echo "Files are copied to $1"
    for FOLDER in `ls -d naf_*`; do
	
	if [[ "${FOLDER}" == *Elqcd* ]] || [[ "${FOLDER}" == *WW* ]] || [[ "${FOLDER}" == *WZ* ]] || [[ "${FOLDER}" == *ZZ* ]] || [[ "${FOLDER}" == *Top* ]]; then
	    merged='MergedFiles/'
	else
	    merged=''
	fi

	if [[ "${FOLDER}" == *JERDn* ]]; then
	    cp -vui ${FOLDER}/????DiffXSec*JERDown*.root     /afs/naf.desy.de/group/cms/scratch/tophh/$1/JERDown/${merged}
	    cp -vui ${FOLDER}/????DiffXSec*JERDown*.txt      /afs/naf.desy.de/group/cms/scratch/tophh/$1/JERDown/TriggerReports/
	elif [[ "${FOLDER}" == *JERUp* ]]; then
	    cp -vui ${FOLDER}/????DiffXSec*JERUp*.root	     /afs/naf.desy.de/group/cms/scratch/tophh/$1/JERUp/${merged}
	    cp -vui ${FOLDER}/????DiffXSec*JERUp*.txt	     /afs/naf.desy.de/group/cms/scratch/tophh/$1/JERUp/TriggerReports/
	elif [[ "${FOLDER}" == *JESDn* ]]; then
	    cp -vui ${FOLDER}/????DiffXSec*JESDown*.root     /afs/naf.desy.de/group/cms/scratch/tophh/$1/JESDown/${merged}
	    cp -vui ${FOLDER}/????DiffXSec*JESDown*.txt     /afs/naf.desy.de/group/cms/scratch/tophh/$1/JESDown/TriggerReports/
	elif [[ "${FOLDER}" == *JESUp* ]]; then
	    cp -vui ${FOLDER}/????DiffXSec*JESUp*.root	     /afs/naf.desy.de/group/cms/scratch/tophh/$1/JESUp/${merged}
	    cp -vui ${FOLDER}/????DiffXSec*JESUp*.txt	     /afs/naf.desy.de/group/cms/scratch/tophh/$1/JESUp/TriggerReports/
	elif [[ "${FOLDER}" == *MassDn* ]]; then
	    cp -vui ${FOLDER}/????DiffXSec*TopMassDown*.root /afs/naf.desy.de/group/cms/scratch/tophh/$1/TopMassDown/${merged}
	    cp -vui ${FOLDER}/????DiffXSec*TopMassDown*.txt /afs/naf.desy.de/group/cms/scratch/tophh/$1/TopMassDown/TriggerReports/
	elif [[ "${FOLDER}" == *MassUp* ]]; then
	    cp -vui ${FOLDER}/????DiffXSec*TopMassUp*.root   /afs/naf.desy.de/group/cms/scratch/tophh/$1/TopMassUp/${merged}
	    cp -vui ${FOLDER}/????DiffXSec*TopMassUp*.txt   /afs/naf.desy.de/group/cms/scratch/tophh/$1/TopMassUp/TriggerReports/
	elif [[ "${FOLDER}" == *MatchDn* ]]; then
	    cp -vui ${FOLDER}/????DiffXSec*MatchDown*.root   /afs/naf.desy.de/group/cms/scratch/tophh/$1/MatchDown/${merged}
	    cp -vui ${FOLDER}/????DiffXSec*MatchDown*.txt   /afs/naf.desy.de/group/cms/scratch/tophh/$1/MatchDown/TriggerReports/
	elif [[ "${FOLDER}" == *MatchUp* ]]; then
	    cp -vui ${FOLDER}/????DiffXSec*MatchUp*.root     /afs/naf.desy.de/group/cms/scratch/tophh/$1/MatchUp/${merged}
	    cp -vui ${FOLDER}/????DiffXSec*MatchUp*.txt     /afs/naf.desy.de/group/cms/scratch/tophh/$1/MatchUp/TriggerReports/
	elif [[ "${FOLDER}" == *ScaleDn* ]]; then
	    cp -vui ${FOLDER}/????DiffXSec*ScaleDown*.root   /afs/naf.desy.de/group/cms/scratch/tophh/$1/ScaleDown/${merged}
	    cp -vui ${FOLDER}/????DiffXSec*ScaleDown*.txt   /afs/naf.desy.de/group/cms/scratch/tophh/$1/ScaleDown/TriggerReports/
	elif [[ "${FOLDER}" == *ScaleUp* ]]; then
	    cp -vui ${FOLDER}/????DiffXSec*ScaleUp*.root     /afs/naf.desy.de/group/cms/scratch/tophh/$1/ScaleUp/${merged}
	    cp -vui ${FOLDER}/????DiffXSec*ScaleUp*.txt     /afs/naf.desy.de/group/cms/scratch/tophh/$1/ScaleUp/TriggerReports/
	else
	    cp -vui ${FOLDER}/????DiffXSec*.root             /afs/naf.desy.de/group/cms/scratch/tophh/$1/${merged}
	    cp -vui ${FOLDER}/????DiffXSec*.txt             /afs/naf.desy.de/group/cms/scratch/tophh/$1/TriggerReports/
	fi

    done

else
    echo "no destination folder given!"
fi