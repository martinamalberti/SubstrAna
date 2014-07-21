import FWCore.ParameterSet.Config as cms

process = cms.Process("PlotROCansSignificance")

process.Options = cms.PSet(

   ## Name of the tree on which perform th training                                                                                                                                
   TreeName           = cms.string("chs"),

   ## Name for ggH signal --> match what is written in the sample list file                                                                                                             
   SignalggHName      = cms.string("RSGWW1000"),

   ## Name for ggH signal --> match what is written in the sample list file                                                                                                             
   SignalqqHName      = cms.string(""),

   ## string which use tree branches to re-weight the events                                                                                                                              
   EventWeight        = cms.string(""),

   ## 0 use both ggH and qqH as signal, 1 use only ggH as signal, 2 use only qqH as signal                                                                                                
   useTypeOfSignal    = cms.uint32(0),

   ## string which is used in the TMVATraining class to define a cut to be applied on the events                                                                                         
   PreselectionCutType = cms.string("basicJetsCutCSA14"),

   ## luminosity in order to  compute signal and bk expectations
   Lumi   = cms.double(19297),

   ## Lepton Type: Muon, Electron ,EleMu and Jets (fully hadronic)                                                                                                                        
   LeptonType         = cms.string("Jets"),

   ## output directory for root and weight file                                                                                                                                           
   outputPlotDirectory  = cms.string("output/TMVATrainingPlots_CHS_SingleVaribales/"),

   ##input file and variables   
   InputInformationParam =  cms.VPSet(

#      cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_mraw_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{raw}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("mass_raw")),

      #cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_mclean_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{clean}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("mass_clean")),


      ### safe mass

     #cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_mtrimsafe_Rtrim_010_Ptfrac_003_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{trim} safe R = 0.1 p_{T}^{frac} = 0.03"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("mass_trim_R01_Pt_003")),

#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_mtrimsafe_Rtrim_020_Ptfrac_005_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{trim} safe R = 0.2 p_{T}^{frac} = 0.05"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("mtrim_R02_Pt005")),

     #cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_msoftdropsafe_beta00_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft drop} safe #beta = 0"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("mass_soft00")),
     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_msoftdropsafe_beta10_0__Gen_Training_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft drop} safe #beta = 1.0"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("mass_soft10")),
     #cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_msoftdropsafe_beta20_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft drop} safe #beta = 2.0"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("mass_soft20")),

     #cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_mprunedsafe_zcut_010_R_cut_050_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{pruned} safe z = 0.1 R_{cut} = 0.5"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("mass_prun_z01_R05")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_mprunedsafe_zcut_010_R_cut_075_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{pruned} safe z = 0.1 R_{cut} = 0.75"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("mass_prun_z01_R075")),

      ##QGlikelihood

#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_QGLikelihood_pr_zcut_010_R_cut_050_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("QGLikelihood Pruned z = 0.1 R_{cut} = 0.5"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("QGLike")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_QGLikelihood_pr_sub1_zcut_010_R_cut_050_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("QGLikelihood Pruned subjet 1 z = 0.1 R_{cut} = 0.5"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("QGLike_pr1")),
     #cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_QGLikelihood_pr_sub2_zcut_010_R_cut_050_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("QGLikelihood Pruned subjet 2 z = 0.1 R_{cut} = 0.75"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("QGLike_pr2")),

     #cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_charge_k05_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("jet charge #kappa = 0.5"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("charge05")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_charge_k07_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("jet charge #kappa = 0.7"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("charge07")),
     #cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_charge_k10_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("jet charge #kappa = 1.0"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("charge10")),

     #cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_tau1_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{1}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("tau1")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_tau1_pr_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{1} pruned"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("tau1_pr")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_tau1_softdrop_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{1} soft drop"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("tau1_soft")),

#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_Qjets_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#Gamma_{qjet}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("tau1_qjets")),

     #cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_ecf_beta_05_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("C_{2}(#beta = 0.5)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("tau1_ECFbeta05")),
     #cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_ecf_beta_10_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("C_{2}(#beta = 1.0)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("tau1_ECFbeta10")),
     #cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_ecf_beta_15_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("C_{2}(#beta = 1.5)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("tau1_ECFbeta15")),
     #cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_ecf_beta_20_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("C_{2}(#beta = 2.0)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("tau1_ECFbeta20")),
  

     #cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_tau2_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("tau2")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_tau2_pr_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2} pruned"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("tau2_pr")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_tau2_softdrop_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2} soft drop"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("tau2_soft")),

     #cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_tau2_0__tau1_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2}/#tau1_{1}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("tau2tau1")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_tau2_pr_0__tau1_pr_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2}/#tau_{1} pruned"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("tau2tau1_pr")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS/TMVATrainingResult_tau2_softdrop_0__tau1_softdrop_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2}/#tau_{1} soft drop"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("tau2tau1_soft")),

#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS_MLP/TMVATrainingResult_msoft10_mraw_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("MLP: m_{Trim} + m_{Raw}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("msoft10_mraw")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS_MLP/TMVATrainingResult_msoft10_mprun_PTBin_475_600_PU_0_100.root"), inputVariableOrMethodName = cms.vstring("MLP: m_{Trim} + m_{prun}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("msoft10_mprun")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS_MLP/TMVATrainingResult_msoft10_mclean_PTBin_475_600_PU_0_100.root"),inputVariableOrMethodName = cms.vstring("MLP: m_{Trim} + m_{clean}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("msoft10_mclean")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS_MLP/TMVATrainingResult_msoft10_tau1_PTBin_475_600_PU_0_100.root"),inputVariableOrMethodName = cms.vstring("MLP: m_{Trim} + #tau_{1}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("msoft10_tau1")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS_MLP/TMVATrainingResult_msoft10_tau2_PTBin_475_600_PU_0_100.root"),inputVariableOrMethodName = cms.vstring("MLP: m_{Trim} + #tau_{2}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("msoft10_tau2")),
     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS_MLP/TMVATrainingResult_msoft10_tau2tau1_PTBin_475_600_PU_0_100.root"),inputVariableOrMethodName = cms.vstring("MLP: m_{Trim} + #tau_{2}/#tau_{1}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("msoft10_tau2tau1")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS_MLP/TMVATrainingResult_msoft10_QGLikelihood_PTBin_475_600_PU_0_100.root"),inputVariableOrMethodName = cms.vstring("MLP: m_{Trim} + QGtag"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("msoft10_QGLike")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS_MLP/TMVATrainingResult_msoft10_QGLikelihood_pr2_PTBin_475_600_PU_0_100.root"),inputVariableOrMethodName = cms.vstring("MLP: m_{Trim} + QGtag sub2"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("msoft10_QGLike_sub2")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS_MLP/TMVATrainingResult_msoft10_ECF20_PTBin_475_600_PU_0_100.root"),inputVariableOrMethodName = cms.vstring("MLP: m_{Trim} + C_{2}(#beta=2)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("msoft10_ECF20")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS_MLP/TMVATrainingResult_msoft10_ECF15_PTBin_475_600_PU_0_100.root"),inputVariableOrMethodName = cms.vstring("MLP: m_{Trim} + C_{2}(#beta=1.5)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("msoft10_ECF15")),

     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS_MLP/TMVATrainingResult_mtrim_tau2tau1_msoft_PTBin_475_600_PU_0_100.root"),inputVariableOrMethodName = cms.vstring("MLP: m_{Trim} + #tau_{2}/#tau[1] + m_{soft}(#beta=1.0)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("mtrim_tau2tau1_msoft")),
     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS_MLP/TMVATrainingResult_mtrim_tau2tau1_mclean_PTBin_475_600_PU_0_100.root"),inputVariableOrMethodName = cms.vstring("MLP: m_{Trim} + #tau_{2}/#tau[1] + m_{clean}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("mtrim_tau2tau1_mclean")),
     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS_MLP/TMVATrainingResult_mtrim_tau2tau1_QGLikelihood_sub2_PTBin_475_600_PU_0_100.root"),inputVariableOrMethodName = cms.vstring("MLP: m_{Trim} + #tau_{2}/#tau[1] + QGtag sub2"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("mtrim_tau2tau1_QGLike")),
     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS_MLP/TMVATrainingResult_mtrim_tau2tau1_ECFbeta15_PTBin_475_600_PU_0_100.root"),inputVariableOrMethodName = cms.vstring("MLP: m_{Trim} + #tau_{2}/#tau[1] + C_{2}(#beta=1.5)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("mtrim_tau2tau1_ECF15")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS_MLP/TMVATrainingResult_mtrim_tau2tau1_ECFbeta20_PTBin_475_600_PU_0_100.root"),inputVariableOrMethodName = cms.vstring("MLP: m_{Trim} + #tau_{2}/#tau[1] + C_{2}(#beta=2.0)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("mtrim_tau2tau1_ECF20")),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS_MLP/TMVATrainingResult_mtrim_combo_PTBin_475_600_PU_0_100.root"),inputVariableOrMethodName = cms.vstring("MLP: All variables"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("mtrim_tau2tau1_combo")),

     cms.PSet( fileName = cms.string("output/outputTMVATraining_CHS_MLP//TMVATrainingResult_AllVariables_PTBin_475_600_PU_0_100.root"),inputVariableOrMethodName = cms.vstring("MLP: All variables"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100),ReducedName = cms.string("all variables")),


   )
)





#  LocalWords:  vdouble
