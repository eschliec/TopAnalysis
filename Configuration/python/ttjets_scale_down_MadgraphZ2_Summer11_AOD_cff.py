import FWCore.ParameterSet.Config as cms

#########################################################
# Dataset: TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola #
# Events : 967055 (70 files)                            # 
# eff    : 1.0                                          #
#########################################################

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [

"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/048A88FA-2BDB-E011-AC35-001A64789D80.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/0E351C63-48DB-E011-881E-001A64789E50.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/12E958D8-38DB-E011-9922-001A6478ABA0.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/18152B0C-51DB-E011-8CBA-002590200B68.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/1846CEAB-55DB-E011-80E4-0025902008F0.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/1898F152-66DB-E011-8932-0025902008C4.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/1A37DDCE-67DB-E011-BC9A-002590200B14.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/1E2BF28A-30DB-E011-9B64-001A64789D18.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/209E2570-34DB-E011-898C-002590200844.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/221C7861-61DB-E011-AF27-002590200A28.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/24B2FF53-BDDB-E011-9B47-003048673E9A.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/2633E0FC-2EDB-E011-9924-001A6478949C.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/264ECB70-2DDB-E011-B120-001A64789508.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/2EF6A75B-2EDB-E011-A3E7-001A64787064.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/34644B45-69DB-E011-B27A-002590200B14.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/36E24C1B-3ADB-E011-93E2-0025901248FA.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/38954110-69DB-E011-9701-001A6478ABB4.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/3A2CF088-5DDB-E011-8641-002590200998.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/3ABF4DD9-62DB-E011-A0B8-002590200A28.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/485A45A6-D1DB-E011-82FC-002590200A6C.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/4863F22C-31DB-E011-B1D4-001A64789508.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/4C909D73-3BDB-E011-AE95-001A6478945C.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/4CB22729-57DB-E011-9912-002590200A40.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/4EF555A7-5BDB-E011-B87C-002590200858.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/50016FCC-32DB-E011-AF0D-001A64789E58.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/5EA2631B-4FDB-E011-9606-002590200A80.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/626EB1CF-31DB-E011-A094-001A64789E20.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/661755C1-3FDB-E011-A8A4-0025902008AC.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/6629B976-37DB-E011-AB3D-001A64789D54.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/6A0E554C-46DB-E011-9EDB-001A647894F0.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/7837AACD-33DB-E011-8D2F-001A64789DEC.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/7E1625F7-4ADB-E011-95F8-0025901248FA.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/86207F80-50DB-E011-96A8-001A64789D94.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/966C26C1-39DB-E011-9770-001A64789508.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/9A731354-2CDB-E011-BB77-001A64789D18.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/9E12056D-4CDB-E011-85D5-002590200B70.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/9E2CBFAE-59DB-E011-BE95-002590200998.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/9E93CB4F-5EDB-E011-84D2-002590200998.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/A20F683A-52DB-E011-A2BB-002590200B74.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/A2D8B7D3-36DB-E011-B609-002590200844.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/A85FBBEF-2FDB-E011-B760-0030486709B4.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/AEC98B76-45DB-E011-BE10-002590200A1C.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/B06600FB-2ADB-E011-BC64-003048673E98.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/BA782B39-5ADB-E011-A054-0025901248FA.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/C0253DF7-5EDB-E011-A467-002590200998.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/C05C31CE-5CDB-E011-B7CB-001A6478AB00.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/CC7CED48-54DB-E011-97D0-001A6478AB00.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/CCBADF6D-68DB-E011-AAB9-003048674048.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/D092B160-4ADB-E011-A34B-002590200A80.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/D23E3BA7-53DB-E011-B711-0025902008B8.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/D656AD3E-56DB-E011-A61E-003048673F9E.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/D65FA18A-58DB-E011-B8AE-002590200964.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/D853E68A-54DB-E011-9B19-00259020080C.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/D8B2B6C8-AFDB-E011-A76D-002590200978.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/DA5C1720-6ADB-E011-ABA2-002590200B14.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/DA8BB0F1-2FDB-E011-B57E-001A64789358.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/DC0B7DA8-66DB-E011-8687-001A6478ABB4.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/DC53612F-36DB-E011-8F6A-001A647894DC.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/DE7E48E8-4FDB-E011-8579-0025902008F0.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/E2101B33-33DB-E011-9340-001A64789D1C.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/E4225A39-C5DB-E011-9F8E-003048673E9A.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/E820EB22-4EDB-E011-A29C-0025902008AC.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/E85E0925-3BDB-E011-8DF4-001A64789D14.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/EEC15720-3DDB-E011-AC9E-001A64789E24.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/F4484252-5CDB-E011-9F01-0025902009BC.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/F489CBF0-43DB-E011-B4A3-001A647894C4.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/F89BEDB9-4CDB-E011-A186-001A64789D94.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/FCADB51D-5BDB-E011-A680-002590200A28.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/FEDB7E2F-32DB-E011-BC50-001A64787064.root",
"/store/mc/Summer11/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/FEEA2016-59DB-E011-8BBB-001A64789D94.root"

     ] );

secFiles.extend( [
               ] )

