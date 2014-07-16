import FWCore.ParameterSet.Config as cms

process = cms.Process("PlotROCansSignificance")

process.Options = cms.PSet(

   ## Name of the tree on which perform th training                                                                                                                                
   TreeName           = cms.string("gen"),

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
   outputPlotDirectory  = cms.string("output/TMVATrainingPlots_puppi_SingleVaribales/"),

   ##input file and variables   
   InputInformationParam =  cms.VPSet(

#      cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_mraw_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{raw}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),

#      cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_mclean_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{cleansing}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),

      ### Not safe mass

#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_mtrim_Rtrim_010_Ptfrac_003_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{trim} R = 0.1 p_{T}^{frac} = 0.03"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),

#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_mtrim_Rtrim_020_Ptfrac_005_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{trim} R = 0.2 p_{T}^{frac} = 0.05"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),

#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_msoftdrop_beta00_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft drop} #beta = 0"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_msoftdrop_beta20_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft drop} #beta = 2.0"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),

#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_mpruned_zcut_010_R_cut_050_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{pruned} z = 0.1 R_{cut} = 0.5"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_mpruned_zcut_010_R_cut_075_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{pruned} z = 0.1 R_{cut} = 0.75"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),

      ### safe mass

#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_mtrimsafe_Rtrim_010_Ptfrac_003_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{trim} safe R = 0.1 p_{T}^{frac} = 0.03"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_mtrimsafe_Rtrim_020_Ptfrac_005_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{trim} safe R = 0.2 p_{T}^{frac} = 0.05"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),

#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_msoftdropsafe_beta00_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft drop} safe #beta = 0"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_msoftdropsafe_beta20_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft drop} safe #beta = 2.0"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),

#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_mprunedsafe_zcut_010_R_cut_050_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{pruned} safe z = 0.1 R_{cut} = 0.5"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_mprunedsafe_zcut_010_R_cut_075_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{pruned} safe z = 0.1 R_{cut} = 0.75"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),

      ##QGlikelihood

#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_QGLikelihood_pr_zcut_010_R_cut_050_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("QGLikelihood Pruned z = 0.1 R_{cut} = 0.5"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_QGLikelihood_pr_sub1_zcut_010_R_cut_050_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("QGLikelihood Pruned subjet 1 z = 0.1 R_{cut} = 0.5"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_QGLikelihood_pr_sub2_zcut_010_R_cut_050_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("QGLikelihood Pruned subjet 2 z = 0.1 R_{cut} = 0.75"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),

#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_charge_k05_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("jet charge #kappa = 0.5"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_charge_k07_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("jet charge #kappa = 0.7"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_charge_k10_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("jet charge #kappa = 1.0"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),

#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_tau1_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{1}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_tau1_pr_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{1} pruned"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_tau1_softdrop_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{1} soft drop"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),

#     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_Qjets_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#Gamma_{qjet}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),

     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_ecf_beta_05_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("C_{2}(#beta = 0.5)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_ecf_beta_10_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("C_{2}(#beta = 1.0)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_ecf_beta_15_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("C_{2}(#beta = 1.5)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_ecf_beta_20_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("C_{2}(#beta = 2.0)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
  


     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_tau2_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_tau2_pr_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2} pruned"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_tau2_softdrop_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2} soft drop"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),

     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_tau2_0__tau1_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2}/#tau1_{1}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_tau2_pr_0__tau1_pr_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2}/#tau_{1} pruned"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
     cms.PSet( fileName = cms.string("output/outputTMVATraining_puppi/TMVATrainingResult_tau2_softdrop_0__tau1_softdrop_0__GenTraining_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2}/#tau_{1} soft drop"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100)),
   
   )
)




