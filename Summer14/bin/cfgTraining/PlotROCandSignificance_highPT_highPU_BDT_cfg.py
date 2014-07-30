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
   outputPlotDirectory  = cms.string("TMVATrainingPlots_CHS_highPT_highPU/"),

   ##input file and variables   
   InputInformationParam =  cms.VPSet(
    cms.PSet( fileName = cms.string("PairVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVATrainingResult_mpruned_tau2tau1_PTBin_475_600_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("BDT: M_{pruned}+#tau_{2}/#tau_{1}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(39,100), ReducedName = cms.string("BDT_1")),
    cms.PSet( fileName = cms.string("PairVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVATrainingResult_msoftbeta10_tau2tau1_PTBin_475_600_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("BDT: M_{soft}(#beta=1)+#tau_{2}/#tau_{1}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(39,100), ReducedName = cms.string("BDT_2")),
    cms.PSet( fileName = cms.string("PairVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVATrainingResult_mpruned_QGLikelihood_PTBin_475_600_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("BDT: M_{pruned}+Q/G"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(39,100), ReducedName = cms.string("BDT_3")),
    cms.PSet( fileName = cms.string("TripletVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVATrainingResult_mpruned_msoftbeta10_tau2tau1_PTBin_475_600_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("BDT: M_{soft}(#beta=1)+M_{pruned}+#tau_{2}/#tau_{1}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(39,100), ReducedName = cms.string("BDT_4")),
    cms.PSet( fileName = cms.string("TripletVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVATrainingResult_mpruned_tau2tau1_QGLikelihood_PTBin_475_600_PU_39_100.root"), inputVariableOrMethodName = cms.vstring("BDT: M_{pruned}+Q/G+#tau_{2}/#tau_{1}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(39,100), ReducedName = cms.string("BDT_5")),
    cms.PSet( fileName = cms.string("AllVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVATrainingResult_allvariables_PTBin_475_600_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("all"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(39,100), ReducedName = cms.string("allvariables")),
  )
)
