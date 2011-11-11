import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/312/98EC336A-7959-E011-B557-0030487CD6B4.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/312/4C2695AA-2C58-E011-A585-0030487CD6DA.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/312/2CC668D2-FC57-E011-B9E8-003048F1C58C.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/312/04850A2C-F757-E011-9A74-003048F024DC.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/311/D20D1C6F-CE57-E011-9DFD-001617E30F48.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/311/3A157852-CA57-E011-B4FB-001D09F253C0.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/311/34081632-5A58-E011-86DB-003048F11CF0.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/311/0E032AAC-DE57-E011-9CA7-0030487C6090.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/310/A6018673-7859-E011-8E0B-0030487C60AE.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/305/26629610-1F57-E011-91B6-0030487A18A4.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/303/90EE0F27-2057-E011-862F-003048F024FA.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/301/B4F34CAA-1D57-E011-9C16-0030487C6A66.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/233/5A323803-6857-E011-8893-001D09F2924F.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/224/8EA23F59-B756-E011-9137-001617C3B6E2.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/223/D2EFD897-095A-E011-94A8-0030487CBD0A.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/222/56910CF0-0557-E011-96BB-001617C3B70E.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/217/FCFC8D4A-FF56-E011-849F-0030487CD718.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/217/E8CD8838-0B5A-E011-A2F6-003048F1BF66.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/217/20E71ED7-0557-E011-B839-0030487CAF5E.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/217/02454E67-FC56-E011-9962-0030487A18A4.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/216/5EE00567-7057-E011-9249-001D09F2546F.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/213/5CDBFD3A-5C56-E011-A277-003048D2C174.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/212/30822193-5656-E011-9366-001D09F276CF.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/210/58731CC3-5856-E011-BC9C-001D09F23A3E.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/205/92FF8D75-5556-E011-84E1-003048F11942.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/176/4E93656C-B456-E011-BB43-001D09F24D67.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/165/141ECB0B-3B57-E011-A1E3-001617C3B76E.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/156/5052B3E9-9B56-E011-8B6F-0030487C6A66.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/119/88B9725B-E956-E011-84FA-003048F117EC.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/117/C8AC5F14-BB56-E011-9399-003048F118AA.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/116/80D44F46-9656-E011-9E8D-001D09F29619.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/113/74343010-C056-E011-9A19-0030487CAEAC.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/107/CEB2E2B2-C956-E011-A052-003048F024C2.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/106/D4BBE1D7-B456-E011-88D2-001617C3B76E.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/103/504F11B1-8257-E011-84BC-001D09F2516D.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/078/C2C1FEE2-0756-E011-8290-0030487C7392.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/076/96391EC6-0456-E011-9C47-000423D9A212.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/020/22EC0FB7-0656-E011-A677-003048F118C4.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/016/F0DECCCB-4B56-E011-818F-0030487CD840.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/161/008/6A776BFE-2656-E011-942F-001D09F2A465.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/998/84FDC215-FA55-E011-B32B-0030487C2B86.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/994/70D87067-7055-E011-9189-0030486780B4.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/957/F2F3DA5E-7155-E011-9699-0019B9F581C9.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/957/D09DF586-5556-E011-8ADF-003048D2C020.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/957/AE52B9F9-C855-E011-8C71-001D09F24448.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/956/24785975-7155-E011-B22E-001D09F290BF.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/955/5E2C0EBC-DC55-E011-8A5B-003048F11114.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/954/1C6B8AE3-2355-E011-AC54-0030487C5CE2.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/943/744486EB-FD54-E011-B8CA-003048F01E88.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/942/E0E5CFBF-F054-E011-BC25-003048F1182E.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/940/FC22868F-0C55-E011-B8FF-001617DBD556.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/939/A861E0C8-1555-E011-8D30-003048F11114.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/938/BCDE6CFD-1355-E011-BECC-001617E30CD4.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/937/6CE913A4-1A55-E011-8C8D-0030487D1BCC.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/936/1AB1F424-EE54-E011-ADA0-003048D2C0F0.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/935/985E9E5C-0955-E011-B0CE-001617E30E2C.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/916/2AA5CB7F-A954-E011-8ABE-003048F118DE.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/915/B64A7F40-AE54-E011-A24A-0030487CBD0A.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/914/E6027900-9454-E011-B7DC-0030487A1FEC.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/913/88115DF8-9854-E011-8E7C-0030487CD710.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/911/AE937F62-C354-E011-8309-003048F1C424.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/907/B8413A26-8D54-E011-9A5B-0030487A3232.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/898/1ACFA1BE-F853-E011-A66C-0030486780A8.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/894/30A2F5D4-5E54-E011-8625-001D09F253D4.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/890/E027ADA1-6854-E011-B954-001617C3B76A.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/888/5EB85F41-6754-E011-9C9F-003048D2C01A.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/877/98C6D13C-E453-E011-84AB-0030487CD6D8.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/876/B8E40FC6-D453-E011-BD4A-003048F024F6.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/875/9E801041-1854-E011-8750-0030487C7828.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/874/F8FC674C-0C54-E011-AE45-0030487CD812.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/873/52BD2EA0-0454-E011-8952-0030486733B4.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/872/72531733-EC53-E011-9B68-003048F110BE.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/871/98BC716B-1554-E011-9D47-0030487A322E.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/853/561A4011-5453-E011-A9AF-001D09F24399.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/835/7EB90ADF-CB53-E011-9376-001D09F25393.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/827/A0EF567C-8A53-E011-B2C3-001617DBD556.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/819/1E76E852-8753-E011-A12B-0030487C8E00.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/815/4CEC100D-EC53-E011-8C9E-001617C3B70E.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/808/7E401B28-8654-E011-A12B-000423D33970.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/744/E8350996-AE52-E011-A994-0030487CD76A.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/737/DAC5CDD7-9452-E011-B97E-003048F024DC.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/579/C6D0574C-8751-E011-AA56-00304879EDEA.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/578/BC0A9E66-8D51-E011-8DCD-000423D9997E.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/577/4827D98B-8951-E011-9445-0030487C5CE2.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/547/6EB39189-B850-E011-8168-0030487C7392.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/502/2E5DA2C5-6B50-E011-9FD8-0030487CF41E.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/500/B220F1E3-6B50-E011-AF3A-0030487CD718.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/499/0096DCAC-9450-E011-B240-0030487CD6D8.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/498/B811B646-8450-E011-BA2D-0030487CD718.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/497/0C7C0EDA-4C51-E011-876F-001617C3B706.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/495/A81FD061-3850-E011-8CE4-0030487CD7B4.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/494/90A3F476-3350-E011-8943-000423D94908.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/488/B48A4365-1850-E011-B59B-001D09F24F65.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/486/10AC6E80-1850-E011-B0D2-0030487CD7E0.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/484/F8331BDE-1650-E011-BB11-003048F0258C.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/472/90BED387-1450-E011-91A4-001617C3B66C.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/471/54BE31A2-1450-E011-AFA4-000423D987E0.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/469/EAE15C8B-1950-E011-B53E-001D09F24F65.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/467/BC94A0DE-1B50-E011-98BE-003048F1C420.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/466/BCDE6B9F-1D50-E011-A25E-000423D987E0.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/463/DCB6FF44-2350-E011-B411-001D09F2432B.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/462/7890956E-2550-E011-A52C-0030487A18D8.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/456/FA31211F-4D50-E011-8563-001617E30D12.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/455/D0882DF8-4950-E011-97C1-001D09F24489.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/454/18D69B16-DC50-E011-8B72-0030487CD178.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/450/FEC7255F-5050-E011-AFC5-0030487CD906.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/449/1444354E-4C50-E011-9FD9-001617E30D40.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/447/322C0F37-5650-E011-A075-001617C3B654.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/446/B00F2D7F-5C50-E011-8244-003048F1BF68.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/445/4052C901-5750-E011-A404-003048D2BDD8.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/444/14E61949-AC50-E011-86D8-0030487C6A66.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/443/C0D9B3E6-5750-E011-9B56-0030487CAEAC.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/442/FC1C6553-5B50-E011-8F17-0019B9F705A3.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/439/DE2FEA95-B34F-E011-8AA6-003048F024E0.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/433/FABB93AB-B54F-E011-A173-001617C3B77C.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/432/060DB184-1150-E011-BB7D-000423D987E0.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/431/10318C91-2D50-E011-A4E1-0030487C7392.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/428/1CF8AC81-B34F-E011-B0E5-003048CFB40C.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/427/50DA7269-814F-E011-8983-00304879FA4A.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/425/143AE138-844F-E011-9685-0030487C6A66.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/423/821466E1-7C4F-E011-B2BA-0030487C6090.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/421/540FAECC-774F-E011-BAC5-003048D2BD66.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/413/0E0C912D-EB4F-E011-A18F-0030487CD718.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/410/965527ED-5F4F-E011-A75C-003048F11CF0.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/406/64A61E4F-604F-E011-A82E-0030487CD840.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/405/F6F4C668-D64F-E011-964C-003048F11C28.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/404/9E3B63A5-374F-E011-A1F3-0030487CD76A.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/403/003DB456-364F-E011-B2D5-001D09F28D4A.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/386/3217F5B7-0A4F-E011-BE53-003048F117B6.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/384/0AA74BB0-084F-E011-915C-0030487CD718.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/383/987A33B3-0C4F-E011-A9DE-0030487CD906.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/379/1E36073F-0C4F-E011-B252-0030487CD6B4.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/377/A2B3291B-0C4F-E011-B868-0030487CD6DA.root',
       '/store/data/Run2011A/MuEG/AOD/PromptReco-v1/000/160/329/A6B1B1DD-2E4E-E011-82DE-0030487CBD0A.root' ] );


secFiles.extend( [
               ] )

