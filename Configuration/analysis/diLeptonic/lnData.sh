#!/bin/zsh

BASEDIR=`pwd`

foreach channel (ee emu mumu)
   
   foreach syst (JESDOWN JESUP JERDOWN JERUP MATCHUP MATCHDOWN MASSUP MASSDOWN PU_UP PU_DOWN TRIG_UP TRIG_DOWN SCALEUP SCALEDOWN MCATNLO POWHEG SPINCORR \
                 BTAG_UP BTAG_DOWN BTAG_LJET_UP BTAG_LJET_DOWN \
                 BTAG_PT_UP BTAG_PT_DOWN BTAG_ETA_UP BTAG_ETA_DOWN \
                 BTAG_LJET_PT_UP BTAG_LJET_PT_DOWN BTAG_LJET_ETA_UP BTAG_LJET_ETA_DOWN \
                 BTAG_BEFF_UP BTAG_BEFF_DOWN BTAG_CEFF_UP BTAG_CEFF_DOWN BTAG_LEFF_UP BTAG_LEFF_DOWN)
   #foreach syst (JESDOWN JESUP RESDOWN RESUP MATCHUP MATCHDOWN MASSUP MASSDOWN PU_UP PU_DOWN SCALEUP SCALEDOWN MCNLOUP MCNLODOWN POWHEGUP POWHEGDOWN)
   #foreach syst (MCNLOUP MCNLODOWN POWHEGUP POWHEGDOWN)
        if [ -d selectionRoot/$syst/$channel ] ; then
            cd selectionRoot/$syst/$channel/
            echo
            echo "Creating Links in ... " 
            pwd
            echo 
            ln -s ../../Nominal/$channel/*.root .
            cd $BASEDIR
        fi
   end

end

foreach SignalSyst (MATCHUP MATCHDOWN MASSUP MASSDOWN SCALEUP SCALEDOWN MCATNLO POWHEG SPINCORR)

    if [ -d selectionRoot/$SignalSyst ] ; then
        rm selectionRoot/$SignalSyst/*/*_ttbarbg.root
        rm selectionRoot/$SignalSyst/*/*_ttbarbgviatau.root
        rm selectionRoot/$SignalSyst/*/*_ttbarsignalplustau.root
    fi
end

# mkdir selectionRoot/HADUP selectionRoot/HADDOWN
# 
# foreach channel (ee emu mumu)
# 
#    mkdir selectionRoot/HADUP/$channel
#    mkdir selectionRoot/HADDOWN/$channel
#    
#    cp selectionRoot/MCATNLO/$channel/* selectionRoot/HADUP/$channel
#    cp selectionRoot/POWHEG/$channel/* selectionRoot/HADDOWN/$channel
# 
#    cd selectionRoot/HADUP/$channel/
#    echo
#    echo "Creating Links in ... " 
#    pwd
#    echo 
#    ln -s ../../Nominal/$channel/* .
#    cd $BASEDIR
# 
#    cd selectionRoot/HADDOWN/$channel/
#    echo
#    echo "Creating Links in ... " 
#    pwd
#    echo 
#    ln -s ../../Nominal/$channel/* .
#    cd $BASEDIR
# 
#    rm selectionRoot/HADUP/$channel/ttbarbg.root
#    rm selectionRoot/HADUP/$channel/ttbarsignalplustau.root
#    rm selectionRoot/HADDOWN/$channel/ttbarbg.root
#    rm selectionRoot/HADDOWN/$channel/ttbarsignalplustau.root
# 
# end

