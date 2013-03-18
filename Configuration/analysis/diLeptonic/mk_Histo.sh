#!/bin/zsh

./lnData.sh
./mk_HistoFileList.sh
root -l -b -q ControlPlots.C++

