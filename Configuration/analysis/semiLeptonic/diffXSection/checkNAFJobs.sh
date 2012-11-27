# CONFIGURE time between checks [seconds]
wait=30

# do some printout and parameter definition
echo "check status every $wait seconds" 
start=$(($(date +%s)-$wait))
time=$(date +%s)

# check 
while [ $wait -ge 0 ]
  do 
  # check waiting time
  now=$(date +%s)
  waited=$(( $now - $start ))
  remaining=$(( $wait - $waited ))
  # if waiting time is over:
  if [ $waited -ge $wait ]
      then
      # perform check
      # LIST all folders here
      echo "$(date)"
      #nafJobSplitter.pl check ./naf_analyzeMuonDiffXSecCorrectedMtop_cfg_MuttbarSGMtop
      #rm ./naf_analyzeMuonDiffXSecCorrectedMtop_cfg_MuttbarSGMtop/err*.txt
      #nafJobSplitter.pl check ./naf_analyzeMuonDiffXSecCorrectedMtop_cfg_MuttbarSGMtop
      nafJobSplitter.pl check ./naf_analyzeMuonDiffXSecCorrectedFourJets_cfg_MuttbarSGFourJets
      #rm ./naf_analyzeMuonDiffXSecCorrectedFourJets_cfg_MuttbarSGFourJets/err*.txt
      #nafJobSplitter.pl check ./naf_analyzeMuonDiffXSecCorrectedFourJets_cfg_MuttbarSGFourJets
      #nafJobSplitter.pl check ./naf_analyzeElecDiffXData2011BPromptRecoV1_cfg
      #nafJobSplitter.pl check ./naf_analyzeMuonDiffXData2011AAug05ReReco_cfg 
      #nafJobSplitter.pl check ./naf_analyzeMuonDiffXData2011APromptRecoV6_cfg
      #nafJobSplitter.pl check ./naf_analyzeMuonDiffXData2011BPromptRecoV1_cfg

      # check all naf folders
      #for FOLDER in `ls -d naf_*`; do
      #    nafJobSplitter.pl check ${FOLDER}
      #done
      # automatic removal of error files
      # rm naf_*/err*.txt

      # reset waiting time
      start=$(date +%s)
  fi
  # else: print remaining waiting time
  # but only once every 15 seconds
  if [ $time -ne $now ]
      then
      if [ $remaining -eq 45 ]
	  then
	  time=$now
	  echo "next check in $remaining seconds"
      fi
      if [ $remaining -eq 30 ]
	  then
	  echo "next check in $remaining seconds"
	  time=$now
      fi
      if [ $remaining -eq 15 ]
	  then
	  echo "next check in $remaining seconds"
	  time=$now
      fi
  fi
done
