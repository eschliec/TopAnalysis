## ---
##    this configfile does the same like analyzeMuonDiffXSec_cfg.py
##    but shifting JES up and without filtering on gen level
## ---

additionalEventWeights  = False
## get the mother file
execfile("analyzeMuonDiffXSecAll_cfg.py")

# JES up
process.scaledJetEnergy.scaleType   = "jes:up"
#process.scaledJetEnergy.scaleFactor = 1.015
if(jetType=="particleFlow"):
    process.scaledJetEnergy.payload = "AK5PF"
elif(jetType=="Calo"):
    process.scaledJetEnergy.payload = "AK5Calo"
else:
    print "unknown jetType"

## change output name 
#process.TFileService.fileName = 'analyzeDiffXSecJESUp_testAll.root'
process.TFileService.fileName = outputFileName+"JESupPF.root"
print "output file name = ", process.TFileService.fileName