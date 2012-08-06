import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/0A4F7A26-718E-E111-9E4B-00261894392B.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/0CB23B2C-6D8E-E111-BCA1-002618943875.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/0E5A73CF-6A8E-E111-9A75-003048D15DB6.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/10196E13-708E-E111-87F3-0018F34D0D62.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/1AE4F5F6-6B8E-E111-B6D6-003048678FB4.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/1EE22DDA-6C8E-E111-9985-003048678BAE.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/044ACEF0-6C8E-E111-A7A1-00261894390B.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/34A033ED-848E-E111-8A4D-003048FFD728.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/38C13299-6E8E-E111-A049-002618943980.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/12A06FB0-6B8E-E111-9C35-0026189437E8.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/14C7B2F8-6B8E-E111-8E21-00261894390C.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/282EFE40-6C8E-E111-AD14-003048678B38.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/4451563B-5D8E-E111-9E87-001A92811736.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/286AEC83-6D8E-E111-B28C-002618943980.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/3C102F45-A58E-E111-B4DE-001BFCDBD1BC.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/4C8B66C8-878E-E111-BB38-001A9281172A.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/5405D1CD-6D8E-E111-98FC-0026189438DF.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/566BE3D2-4D8E-E111-9344-001A92810AA6.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/587D3E25-698E-E111-AF3B-0018F3D0965C.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/5C029AB8-698E-E111-A098-001A92971B90.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/5E215741-6C8E-E111-8C77-002618943964.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/600CEB61-6B8E-E111-86E8-003048678ADA.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/3C1B9667-6B8E-E111-8A5A-002618943876.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/3C3F511D-6D8E-E111-8D19-003048679084.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/6266F69D-868E-E111-868C-0018F3D0966C.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/629746CC-6D8E-E111-BFD9-002618943866.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/3CFAB686-6A8E-E111-838D-003048678F6C.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/3E25BB30-858E-E111-92FE-003048679076.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/78086B1C-6D8E-E111-94FD-003048679236.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/78BB0113-6E8E-E111-B960-0026189438DA.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/7AE21AB5-6B8E-E111-A70E-0018F3D096E6.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/84C92B06-6E8E-E111-9AE5-002618943930.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/86376E30-878E-E111-8EEC-001A92971B48.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/3EFDD245-6C8E-E111-BC4D-003048679214.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/8C758CF8-6B8E-E111-BC22-003048679296.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/60D01994-6E8E-E111-8D98-002618943901.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/643CF899-888E-E111-B217-001A92810AEA.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/A22EFF77-878E-E111-8C16-0018F3D096EE.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/A8B39549-8A8E-E111-A6AB-002354EF3BDB.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/76A5DD41-6C8E-E111-92F6-002618943978.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/868C6CEB-9B8E-E111-BC9F-0018F3D09692.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/8CA0D4A5-918E-E111-BBAC-0018F3D096E0.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/B8937C8B-6C8E-E111-A54F-0026189438D3.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/9C2F1266-6D8E-E111-8A60-003048678C26.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/CA7BC062-6B8E-E111-9D3A-003048679214.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/AA7F2086-6D8E-E111-B342-0026189438CB.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/D4A2A26D-898E-E111-8D5C-002354EF3BDB.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/D604D2C6-838E-E111-9E15-00261894388B.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/B0240F1C-6B8E-E111-BF8E-003048678ADA.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/DC977B08-6E8E-E111-9439-0026189438E3.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/DEF4FBB2-6B8E-E111-93A0-001A92971AD0.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/B426FF69-6D8E-E111-9F7E-00304867929E.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/E83D268A-6C8E-E111-A265-002354EF3BDE.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/C6C23C83-6D8E-E111-9EBE-00261894397A.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/E8DC9C42-6C8E-E111-B7EF-003048678B06.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/F406D1F6-6B8E-E111-A813-00261894388F.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/F4513F28-6D8E-E111-8E6C-00261894398B.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/CE827DD2-6A8E-E111-97F2-003048678FEA.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/D8A201D7-6C8E-E111-86AF-002618943943.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/FAD2E60D-868E-E111-AEFE-001BFCDBD11E.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/E49EFABB-8B8E-E111-84F8-002354EF3BDB.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/E85EA89E-5A8E-E111-BC83-003048678BC6.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/F6CCD873-6F8E-E111-A55F-002618943876.root',
       '/store/mc/Summer12/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S7_START52_V9-v1/0000/F6FC3DF0-6C8E-E111-A2A1-002618943800.root' ] );


secFiles.extend( [
               ] )

