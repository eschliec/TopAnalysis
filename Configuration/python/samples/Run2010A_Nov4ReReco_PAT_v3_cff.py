import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_100_1_j4M.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_101_1_cU1.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_10_1_J54.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_102_1_YyA.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_103_1_RRE.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_104_1_Q7A.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_104_3_Dey.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_105_1_zS3.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_106_1_x9k.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_107_1_bdt.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_108_1_FBV.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_109_1_Knb.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_110_1_s3J.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_111_1_wei.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_11_1_SIQ.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_112_1_tlq.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_113_1_iZX.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_114_1_xTg.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_115_1_Wce.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_116_1_C3o.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_117_1_1Sz.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_118_3_pC5.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_119_2_5om.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_1_1_D4v.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_120_1_xEv.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_121_1_mud.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_12_1_ruS.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_122_1_HjR.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_123_1_6Zi.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_124_1_mCl.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_125_3_l8C.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_126_4_N7f.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_127_1_jUQ.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_128_2_b1W.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_129_1_a2u.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_130_1_Ds0.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_131_1_P2P.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_132_1_mea.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_132_2_q1j.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_13_2_ERk.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_133_1_Swj.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_134_1_GUr.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_135_1_3ox.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_136_1_hoL.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_137_1_6Wd.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_138_3_lSt.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_139_3_1wy.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_140_3_rs5.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_141_1_WUV.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_141_2_OIl.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_142_2_6w6.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_143_4_wpL.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_144_3_Dc5.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_14_4_OHm.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_145_1_xRk.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_145_2_uTz.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_146_1_Mwj.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_147_1_1lw.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_148_1_GkU.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_149_1_pr2.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_150_1_Ubn.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_151_3_5YH.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_15_1_7GD.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_152_4_zHg.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_153_4_P1R.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_154_2_uYE.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_155_1_X42.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_156_1_zG8.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_157_1_zqa.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_16_1_96O.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_17_1_MFv.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_18_1_U5t.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_19_3_22X.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_20_3_SuN.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_21_1_Wbb.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_2_1_T4l.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_22_1_hRz.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_23_1_Syd.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_24_4_FP4.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_25_3_6Fk.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_26_1_dBg.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_27_2_7Zy.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_28_2_9MT.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_28_2_W4k.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_29_1_kh9.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_30_1_KWm.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_31_1_LB7.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_3_1_q8r.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_32_3_lbY.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_33_1_Rq6.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_34_3_UE3.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_35_1_nrB.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_36_1_04W.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_37_1_iJ7.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_37_2_iKh.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_38_1_FdH.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_39_1_bWt.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_40_1_Hjf.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_41_1_yBV.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_4_1_zBx.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_42_1_Pf2.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_43_1_M2d.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_44_1_UK6.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_45_1_5GL.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_46_1_O5r.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_47_1_mqj.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_48_1_LFx.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_49_1_CqX.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_50_1_jPr.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_51_1_Ipu.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_5_1_yeK.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_52_1_Zmo.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_53_1_cyD.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_54_1_vn4.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_55_3_VDW.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_56_1_3P0.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_57_1_TWb.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_58_1_K4N.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_59_1_DnQ.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_60_1_gph.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_60_2_gSr.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_60_3_LsT.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_61_1_7aL.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_6_1_LY1.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_62_1_rBO.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_63_1_9WD.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_64_1_j5k.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_65_1_ktY.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_66_1_MjO.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_67_1_kaK.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_68_1_KLE.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_69_1_SRI.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_70_1_24d.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_71_1_zLz.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_7_1_ziZ.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_72_1_v4r.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_73_1_Jmb.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_74_1_wDO.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_75_1_CVP.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_76_1_mWq.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_77_1_6lI.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_78_1_m60.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_79_1_UIB.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_80_1_kBK.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_81_1_Ph9.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_8_1_tXf.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_82_1_ugq.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_83_1_AsG.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_84_1_Mfl.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_85_1_Aq3.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_86_1_AB6.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_87_1_hjT.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_88_1_bEa.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_89_1_4UI.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_90_1_RsO.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_91_1_NqM.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_92_1_ZlO.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_93_1_ccC.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_9_3_8ck.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_94_1_WjO.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_95_1_Diq.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_96_1_pLy.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_97_1_OAZ.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_98_1_g8X.root',
        '/store/user/wbehrenh/Mu/Run2010A-Nov4ReReco-PAT-v3/267a23dfb7d045a9a0ecfb5f145848ff/rereco_99_1_VDa.root'] );

secFiles.extend( [
               ] )
