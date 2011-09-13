import FWCore.ParameterSet.Config as cms

############################################################
# Dataset: TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola #
# Events : 1065323 (83 files)                              # 
# eff    : 1.0                                             #
############################################################

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [

"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/00F082EB-D9D5-E011-AF2A-002590200B48.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/0234417B-C8D5-E011-BBDE-001A64789354.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/0C3C5967-B7D5-E011-9BDE-001A64789DD0.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/120BE00A-C2D5-E011-9CD7-001A6478ABB4.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/163B7F42-DFD5-E011-BC56-0025902008C8.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/16C9DC0C-D9D5-E011-B467-001A64789DA0.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/1A42A97E-C3D5-E011-84DC-001A64789E20.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/223C7DC0-B5D5-E011-A117-001A64789E40.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/224438E0-B4D5-E011-9B11-00304866C368.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/228EB8A4-CBD5-E011-A77A-002590200868.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/289415D5-CDD5-E011-980F-001A64789D84.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/2A362CD2-B1D5-E011-867F-00304867098C.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/342EC188-E6D5-E011-B928-0025902009A4.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/34F891A0-D9D5-E011-B9B6-002590200B44.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/38074261-ECD5-E011-9030-002590200B34.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/38967981-BFD5-E011-A202-001A64789D10.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/44F781B8-CCD5-E011-8DA5-001A64789454.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/4631C901-DDD5-E011-B47D-002590200898.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/486080CB-C6D5-E011-A151-001A64789464.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/4CA68D5D-E2D5-E011-A84F-0025902008C8.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/503CB43F-D6D5-E011-A11E-002590200878.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/5A2B0518-E6D5-E011-8F9B-002590200998.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/5A904513-D2D5-E011-BED9-0025902008F4.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/5E0DC89A-C9D5-E011-92FD-001A64789D84.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/605121F7-BED5-E011-A6FE-001A64789D10.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/66F19DA9-B7D5-E011-9FC5-001A64789458.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/68EEA36D-D0D5-E011-96EE-0025902008C8.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/6EA32B81-CAD5-E011-94AD-002590200974.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/76D0730E-D4D5-E011-9C7B-0025902008CC.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/78AD9DC4-BAD5-E011-B345-001A6478AA34.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/7C3094ED-4AD6-E011-A37B-002590200A68.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/7C50C058-CBD5-E011-A0C1-002590200838.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/7EC08853-D2D5-E011-9D1C-001A64789464.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/8277F131-C1D5-E011-9A89-001A64789D84.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/84825260-51D6-E011-82E8-0025901248FA.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/84C9B82E-CCD5-E011-B970-002590200AB8.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/8C348AD6-DFD5-E011-B3F6-0025902008C8.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/8EBABF1E-D5D5-E011-9773-0025902008FC.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/909FD6BD-D3D5-E011-A19E-001A64789354.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/988D7988-CFD5-E011-9DDA-002590200B0C.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/9C18CF3B-BAD5-E011-883D-001A64789D7C.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/9C697243-CDD5-E011-9561-001A64789354.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/A04855EE-BBD5-E011-B859-001A64789458.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/A25CB3C4-DAD5-E011-A49D-001A64789DA0.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/A2740B97-5BD6-E011-83B4-002590200B74.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/A2E01775-D3D5-E011-AC89-002590200934.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/A44781E9-D7D5-E011-9089-003048673F9E.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/A6D8FBD4-EDD5-E011-BACB-0025902009BC.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/A6FE66CC-DDD5-E011-BD3D-001A64789D84.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/AC2EA04E-C4D5-E011-B5E4-001A64789E20.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/AE2A8B5A-BBD5-E011-8732-001A64789E04.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/B4B6A8EF-EED5-E011-A036-0025902009BC.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/B641AC5E-30D6-E011-8637-002590200B6C.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/B80632ED-D2D5-E011-9507-002590200934.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/BC43DE12-C0D5-E011-A0C2-003048670BCC.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/BC9A5EA9-D4D5-E011-A8BC-002590200878.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/BE54DCA6-C0D5-E011-8BA5-001A64789464.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/C418DFC2-BED5-E011-9900-001A64789474.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/C4564A58-E9D5-E011-9F0B-002590200B6C.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/C81044AC-CED5-E011-9E6A-00304866C368.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/CA93DED4-B8D5-E011-AA15-001A64789508.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/CE0F0BF8-C5D5-E011-AF4D-003048670BF8.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/CE2F4A66-D5D5-E011-8128-002590200B34.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/D00F0EEB-B9D5-E011-B08D-001A64789490.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/D641FA3B-9FD5-E011-BA68-0025902008B8.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/D6787A98-BDD5-E011-B501-001A64789DA8.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/DCAB9117-D0D5-E011-B21D-002590200934.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/E03F309A-E4D5-E011-A416-002590200AB8.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/E05476CA-D6D5-E011-ACA4-001A64789DF4.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/E065BE79-D8D5-E011-BFDF-001A64789DF4.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/E0BB368F-B8D5-E011-ADBE-001A64789DF0.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/E8FBC9C8-BCD5-E011-AA56-001A6478ABB4.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/F008C2B0-E0D5-E011-9CC4-0025902009A4.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/F20FEA6A-D7D5-E011-B767-001A64789DF4.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/F22C7062-C9D5-E011-9B11-0025902008C4.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/F8956DEE-D0D5-E011-B3EB-0025901248FA.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/F8D3E999-B6D5-E011-9713-001A64789D28.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/FA4D7800-B4D5-E011-872C-00304867098C.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/FAC6A845-42D6-E011-83FB-0025901248FA.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/FE637CBF-3FD6-E011-8F66-0030486709FE.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/FE7FD210-B9D5-E011-AE9D-003048673F24.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/FE9CCED7-C4D5-E011-8EED-0025902009A4.root",
"/store/mc/Summer11/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola//AODSIM/PU_S4_START42_V11-v1/0000/FEFC7AEB-DBD5-E011-9CE1-002590200B7C.root"

     ] );

secFiles.extend( [
               ] )

