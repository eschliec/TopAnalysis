import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/379/FC52D477-0B4F-E011-B845-003048F1C58C.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/383/DC8108ED-0C4F-E011-AECD-003048F118AA.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/386/D2D5C498-0A4F-E011-8A3E-0030487CD77E.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/403/1A5D2F51-364F-E011-AF2D-003048F11C58.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/405/58F4BD7E-D64F-E011-8747-0030487CF41E.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/406/D68A88D0-5F4F-E011-BFAC-003048F118C6.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/410/AEE7EEBC-5F4F-E011-A1DF-003048F118C4.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/329/86ABC157-2F4E-E011-9004-0030487C778E.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/377/F8B59222-0C4F-E011-B15D-0030487CD14E.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/384/408DE5C7-084F-E011-A63D-0030487CD710.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/404/2E055ABE-374F-E011-A659-001D09F291D2.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/413/BC0268DA-EA4F-E011-B139-0030487CF41E.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/421/B65508AC-774F-E011-90F4-000423D94E70.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/431/9089836C-2D50-E011-A661-0030487CD812.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/423/4028D013-7D4F-E011-96C4-0030487C8CB6.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/425/C846587E-844F-E011-8D2C-0030487A195C.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/439/7C677417-B44F-E011-9B26-001D09F24D67.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/427/8EF3129F-814F-E011-BC8D-003048F117EC.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/428/A89255C6-B34F-E011-855A-000423D9890C.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/444/F8E31566-AC50-E011-86CF-0030487CD76A.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/432/3EFDAAD0-1150-E011-A786-001D09F290CE.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/433/5C7370AB-B54F-E011-950F-0019B9F72CE5.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/447/9CDBA3F7-5550-E011-A0B6-0030487CD716.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/449/F49AF06A-4C50-E011-B35F-003048F024E0.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/442/C817DE13-5B50-E011-B811-001617E30CC2.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/443/A298111A-5850-E011-AAA2-0030487CAEAC.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/445/B08D573B-5750-E011-B267-001D09F24498.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/456/0E955726-4D50-E011-AC95-003048D375AA.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/446/6638D85E-5C50-E011-8089-000423D94494.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/450/427D557C-5050-E011-97FC-003048F118AA.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/466/1CF33D1C-1E50-E011-9E47-0030487A17B8.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/467/105C5AC8-1B50-E011-B67F-0030487CD17C.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/454/7CAEAA07-DC50-E011-9C82-0030487CD17C.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/471/FCAE0FDB-1450-E011-A940-0030487A18D8.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/472/9E1502E5-1450-E011-A36D-000423D996C8.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/455/DC710FF9-4950-E011-AE57-001D09F251BD.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/462/8C0DB6A3-2550-E011-A775-001D09F28F1B.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/488/32D7C707-1850-E011-A381-001D09F282F5.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/494/400F249A-3350-E011-B009-003048F1110E.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/495/88266F95-3850-E011-9374-003048F118D4.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/463/FA5FB76B-2350-E011-AD20-0030487C635A.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/469/7015E1A2-1950-E011-9006-00304879EE3E.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/498/C26D6392-8450-E011-87AE-003048F11942.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/484/B43147D5-1650-E011-8367-001D09F2527B.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/500/2EB6D038-6C50-E011-AEAF-0030487C7392.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/486/7C1D2C04-1850-E011-9B47-001617E30D12.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/497/3AF418A7-4C51-E011-BAAB-001617C3B654.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/577/DC040FD8-5451-E011-9385-000423D9997E.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/578/1E802607-8D51-E011-AF95-0030487CD716.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/497/9A11C408-8750-E011-B733-0030487CD6DA.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/499/2E5865AC-9450-E011-A90D-0030487C2B86.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/502/D0D688EB-6B50-E011-A71C-0030487C7E18.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/547/3CB3AC28-B850-E011-904C-0030487CD6B4.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/578/D21DCE6D-5F51-E011-A17D-0030487C7828.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/579/D236C873-8751-E011-8FF5-0030487A18D8.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/819/5696FE52-8753-E011-A826-0030487CD6E8.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/737/5CB739D0-9452-E011-BEF7-001D09F2924F.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/744/10BEBABB-AE52-E011-BC0E-0030487C7392.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/808/82CAF258-8654-E011-890A-003048F1BF68.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/835/B6CE0D8F-CB53-E011-A605-003048F118AC.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/815/0C4C6809-EC53-E011-B373-003048F118AA.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/871/2EFC0B8D-1554-E011-AD06-0030487CBD0A.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/827/C80C856D-8A53-E011-9C96-001617C3B6CC.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/835/9E483420-8A53-E011-8AF7-001D09F24FBA.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/835/A47B0830-9053-E011-ADAE-0030487CD13A.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/853/14717021-5453-E011-BFA3-001D09F29146.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/875/AA4AFBD2-F753-E011-8640-003048F024E0.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/872/5ADB6543-EB53-E011-AD21-003048F110BE.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/873/F28C7DA0-0454-E011-98A9-003048D2C0F2.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/874/40A4ED22-0C54-E011-98C2-001617DBD5AC.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/875/40D5C869-1854-E011-AEF7-0030487C8E00.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/876/367129C0-D453-E011-9BE9-0030487CD6F2.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/877/F430BC90-E453-E011-86FB-001D09F23944.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/890/44DDA8D2-4554-E011-8524-000423D33970.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/888/6056C9A5-3554-E011-B033-001D09F253C0.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/888/C45F1631-4054-E011-9DFF-001D09F291D7.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/888/F8799441-3B54-E011-A426-001D09F2A690.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/907/1C5B9D31-8D54-E011-85EE-0030487CD7B4.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/890/34346E5D-6854-E011-9353-000423D987E0.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/911/5A2F688C-C354-E011-AF51-0030487C608C.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/890/68438EE0-4054-E011-AB93-003048CFB40C.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/894/963910BD-5E54-E011-97C7-000423D9890C.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/915/64FA48EE-9754-E011-B022-000423D996C8.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/898/AC7C51CA-F853-E011-B24F-001617E30D4A.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/911/2C4072F5-8B54-E011-9481-0030487CD7E0.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/913/2AE3BC0B-9954-E011-86CB-003048F11C5C.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/914/1E542106-9454-E011-940B-0030487CD14E.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/915/B289BC36-AE54-E011-B242-001D09F24498.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/916/9C908216-A954-E011-A230-003048CFB40C.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/937/8E785192-1A55-E011-8D03-001617C3B70E.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/935/E2F8B39E-0955-E011-8157-001D09F23944.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/935/F88D52FF-F054-E011-AA2E-001D09F24DDA.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/936/C487B619-EE54-E011-BD7C-0030487CD178.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/937/6CDD1E85-E554-E011-9B51-001617C3B6CC.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/938/82D09EAF-1355-E011-8683-001D09F2841C.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/939/B0613FDA-1555-E011-9778-003048F118D2.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/955/2874ECC9-4D55-E011-BACF-001617C3B70E.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/940/54623440-0C55-E011-894C-001D09F2906A.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/942/C8B7E91D-F154-E011-B7BA-001617E30F50.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/943/C2B9928F-FD54-E011-B0CC-001D09F23A84.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/954/CE87F535-2455-E011-B345-001617C3B6E2.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/955/8878DAAB-DC55-E011-AD74-003048F110BE.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/957/52E7DFEC-6C55-E011-BE20-0030487C7392.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/956/D2A536C1-7155-E011-BB0F-003048D375AA.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/957/5692014D-6955-E011-AEB8-003048F024E0.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/957/04859906-6855-E011-BD7D-000423D98950.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/957/346A08A2-6B55-E011-BD84-003048D2BED6.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/957/4009983C-6455-E011-A4B2-003048F11DE2.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/957/B8FDF505-7655-E011-8C57-0030487CD17C.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/957/56556E7D-7C55-E011-99CA-003048D374F2.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/957/68D643E6-7355-E011-AB8D-0030487CF41E.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/957/9081D154-A455-E011-8D98-001D09F24DA8.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/957/B0C4CF58-6255-E011-8C79-003048F117B4.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/957/BA2D7A4E-5955-E011-AA08-000423D98950.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/957/D02E34A2-6B55-E011-9FB6-003048D2BC42.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/957/E61348FB-6255-E011-91A0-003048F024FA.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/957/EA746616-5556-E011-9575-003048D2C16E.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/994/726CAB23-7055-E011-A7A5-0030487CD76A.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/998/78A60216-FA55-E011-8EEB-0030487C7392.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/998/A677EA47-C855-E011-B55F-001D09F25460.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/016/AEAB2A55-1B56-E011-B4FE-003048D2BE08.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/160/998/EC30B831-C655-E011-A540-003048F024DE.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/008/00C8884A-2756-E011-9664-0030487C90C2.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/008/8468EAFF-C855-E011-BB51-003048F117B4.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/016/F2ED8F55-E455-E011-9D3A-0016177CA7A0.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/020/8CF47390-0656-E011-A9C1-001617DBD556.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/008/E0E3BA24-D455-E011-9A66-003048F11942.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/016/BC1553FC-DD55-E011-ADD1-000423D996C8.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/016/C68F03B4-4B56-E011-B66C-0016177CA778.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/016/ECC00879-DA55-E011-AE3B-003048F1C424.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/076/D0833818-0556-E011-8092-0030487C778E.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/107/0A8FCE1C-CA56-E011-8260-0030487C60AE.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/078/CAB9DBEC-0756-E011-B743-0030487CD812.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/103/C0CFC57F-8257-E011-B06E-0030487C90EE.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/103/CC6ADB47-9A56-E011-949A-00304879FC6C.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/106/E8123EB1-B456-E011-A961-001D09F24024.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/113/801B64FD-BF56-E011-97C9-001617C3B65A.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/116/C4F81437-9656-E011-AE54-001617C3B6C6.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/117/FABC76AB-BB56-E011-9107-001D09F24DA8.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/165/708B19CC-3A57-E011-A917-003048F118DE.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/119/B4B6D1B9-E956-E011-B842-003048D2BD66.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/119/FA85F4DE-9156-E011-B081-0030486780B4.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/156/2240DB6A-9756-E011-8E59-001617C3B654.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/156/92461CFF-9B56-E011-A86C-0030487CAEAC.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/176/70610D6B-B456-E011-96BE-001D09F27003.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/205/904ABFA5-5556-E011-88EA-003048D2BC42.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/217/1C5F5E2A-0B5A-E011-972F-0030487C6A66.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/210/28C2688F-5856-E011-867A-001617C3B710.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/212/E46F24C1-5656-E011-B625-001617E30F48.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/213/6405AA25-5C56-E011-AE18-001617E30D12.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/216/70EE0E4D-7057-E011-BD5D-001D09F24498.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/217/1CBCE6F0-F356-E011-B213-0030487CAF5E.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/217/22304155-F756-E011-906F-001617E30E2C.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/217/A49D8DA3-EF56-E011-864D-003048F118C6.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/217/24AC4841-0857-E011-905C-003048F024E0.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/217/547D46A5-DC56-E011-93C4-0030486733D8.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/217/EADFA0B9-0B57-E011-AE06-003048F118E0.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/217/824E6956-F056-E011-90C9-000423D9890C.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/217/9CF5950C-1557-E011-941B-001D09F24691.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/217/E2F8D4EA-0857-E011-8831-003048F024DE.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/223/70333401-FD56-E011-9DEF-003048F118D4.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/217/E85E7124-0D57-E011-9B43-003048F024F6.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/222/3C7B8CF0-0557-E011-A5A8-003048F117EC.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/223/082EF30C-0457-E011-B316-003048F117B6.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/223/50599BD9-095A-E011-81DA-003048F024E0.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/223/8E11753A-FF56-E011-8624-0030487CD6E6.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/224/A03C1672-B756-E011-A2C4-001D09F2516D.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/233/D80940E2-6757-E011-860F-001617C3B66C.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/310/DA87DC85-7859-E011-BFC4-003048F11114.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/301/AED0A75A-1D57-E011-893B-003048F024DE.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/303/AE1FB190-2057-E011-9D6F-003048F117B6.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/311/20212DA7-BF57-E011-94AC-000423D9A2AE.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/305/56667E4E-1F57-E011-8F21-0030487C2B86.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/310/3E6A99E5-B757-E011-8EB4-001617E30D0A.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/311/0AFAE2E6-CA57-E011-AF97-001D09F29114.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/311/1A43F738-CA57-E011-89EE-0030487A3232.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/311/48E70A0E-BA57-E011-AD8E-0030487CD906.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/311/7EE0BA64-C057-E011-882A-001D09F29114.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/311/9CD4CC5A-5A58-E011-8F2D-0030487CD76A.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/311/AE6E3D99-D757-E011-A66C-003048F117B4.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/311/C244A2CB-BA57-E011-AAB9-001617C3B6E2.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/311/EE26DF9F-D257-E011-83F9-003048F117B4.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/311/FE75C6F7-C557-E011-A2B5-0030487C6062.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/312/0E8D9A70-E857-E011-B3DC-001D09F29321.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/312/2C8E8D20-0358-E011-B7E8-00304879FBB2.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/312/32AE26CD-E257-E011-A0C8-003048F024F6.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/312/40D0CE52-7959-E011-90DB-003048F118C2.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/312/8CC8FA1B-E257-E011-A644-001D09F26C5C.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/312/70491EF9-F957-E011-9024-003048F11942.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/312/D81D11DF-E957-E011-B822-001617DBD472.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/312/7284A8EE-0558-E011-9AC0-003048F1C424.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/312/7609160B-EE57-E011-9149-001617E30D12.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/312/ECDD3154-DF57-E011-9088-001D09F24D4E.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/312/822C8270-E857-E011-AB78-001D09F25460.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/312/902B9A5F-ED57-E011-9C03-003048F117B6.root',
       '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/312/E82604AF-E057-E011-A89E-001D09F2915A.root' ] );


secFiles.extend( [
               ] )

