import FWCore.ParameterSet.Config as cms

process = cms.Process("PlotCorrealtionBDT")

process.Options = cms.PSet(

  InputLowPUFiles = cms.VPSet(
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_ECFbeta10_PTBin_475_600_PU_0_39.root"), variableName=cms.string("C2(#beta=1)")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_ECFbeta15_PTBin_475_600_PU_0_39.root"), variableName=cms.string("C2(#beta=1.5)")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_ECFbeta20_PTBin_475_600_PU_0_39.root"), variableName=cms.string("C2(#beta=2)")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_QGLikelihood_PTBin_475_600_PU_0_39.root"), variableName=cms.string("QG Likelihood")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_QGLikelihood_sub1_PTBin_475_600_PU_0_39.root"), variableName=cms.string("QG Likelihood sub1")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_QGLikelihood_sub2_PTBin_475_600_PU_0_39.root"), variableName=cms.string("QG Likelihood sub2")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_QGLikelihood_comb_PTBin_475_600_PU_0_39.root"), variableName=cms.string("QG Likelihood comb")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_QGLikelihood_comb_PTBin_475_600_PU_0_39.root"), variableName=cms.string("QG Likelihood comb")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_Qjets_PTBin_475_600_PU_0_39.root"), variableName =  cms.string("#Gamma_{Qjets}")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_mclean_PTBin_475_600_PU_0_39.root"), variableName =  cms.string("M_{clean}")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_mconst_PTBin_475_600_PU_0_39.root"), variableName =  cms.string("M_{const}")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_mpruned_PTBin_475_600_PU_0_39.root"), variableName =  cms.string("M_{pruned}")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_msoftbeta00_PTBin_475_600_PU_0_39.root"), variableName =  cms.string("M_{soft}(#beta=0)")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_msoftbeta10_PTBin_475_600_PU_0_39.root"), variableName =  cms.string("M_{soft}(#beta=1)")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_msoftbeta20_PTBin_475_600_PU_0_39.root"), variableName =  cms.string("M_{soft}(#beta=2)")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_mtrim_PTBin_475_600_PU_0_39.root"), variableName =  cms.string("M_{trim}")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_pullangle_PTBin_475_600_PU_0_39.root"), variableName =  cms.string("Pull Angle")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_tau1_PTBin_475_600_PU_0_39.root"), variableName =  cms.string("#tau_{1}")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_tau2_PTBin_475_600_PU_0_39.root"), variableName =  cms.string("#tau_{2}")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_tau2_tau1_PTBin_475_600_PU_0_39.root"), variableName =  cms.string("#tau_{2}/#tau_{1}")),
     cms.PSet(fileName = cms.string("AllVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU_BDTG/TMVATrainingResult_allvariables_PTBin_475_600_PU_0_39.root"), variableName =  cms.string("all")),
  ),

  InputHighPUFiles = cms.VPSet(
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_ECFbeta10_PTBin_475_600_PU_39_100.root"), variableName=cms.string("C2(#beta=1)")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_ECFbeta15_PTBin_475_600_PU_39_100.root"), variableName=cms.string("C2(#beta=1.5)")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_ECFbeta20_PTBin_475_600_PU_39_100.root"), variableName=cms.string("C2(#beta=2)")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_QGLikelihood_PTBin_475_600_PU_39_100.root"), variableName=cms.string("QG Likelihood")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_QGLikelihood_sub1_PTBin_475_600_PU_39_100.root"), variableName=cms.string("QG Likelihood sub1")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_QGLikelihood_sub2_PTBin_475_600_PU_39_100.root"), variableName=cms.string("QG Likelihood sub2")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_QGLikelihood_comb_PTBin_475_600_PU_39_100.root"), variableName=cms.string("QG Likelihood comb")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_Qjets_PTBin_475_600_PU_39_100.root"), variableName =  cms.string("#Gamma_{Qjets}")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_mclean_PTBin_475_600_PU_39_100.root"), variableName =  cms.string("M_{clean}")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_mconst_PTBin_475_600_PU_39_100.root"), variableName =  cms.string("M_{const}")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_mpruned_PTBin_475_600_PU_39_100.root"), variableName =  cms.string("M_{pruned}")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_msoftbeta00_PTBin_475_600_PU_39_100.root"), variableName =  cms.string("M_{soft}(#beta=0)")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_msoftbeta10_PTBin_475_600_PU_39_100.root"), variableName =  cms.string("M_{soft}(#beta=1)")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_msoftbeta20_PTBin_475_600_PU_39_100.root"), variableName =  cms.string("M_{soft}(#beta=2)")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_mtrim_PTBin_475_600_PU_39_100.root"), variableName =  cms.string("M_{trim}")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_pullangle_PTBin_475_600_PU_39_100.root"), variableName =  cms.string("Pull Angle")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_tau1_PTBin_475_600_PU_39_100.root"), variableName =  cms.string("#tau_{1}")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_tau2_PTBin_475_600_PU_39_100.root"), variableName =  cms.string("#tau_{2}")),
     cms.PSet(fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_tau2_tau1_PTBin_475_600_PU_39_100.root"), variableName =  cms.string("#tau_{2}/#tau_{1}")),
     cms.PSet(fileName = cms.string("AllVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU_BDTG/TMVATrainingResult_allvariables_PTBin_475_600_PU_39_100.root"), variableName =  cms.string("all")),
  )

)
