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
   outputPlotDirectory  = cms.string("TMVATrainingPlots_CHS_lowPT_highPU/"),

   ##input file and variables   
   InputInformationParam =  cms.VPSet(
    cms.PSet( fileName = cms.string("PairVariablesTraining_lowPT_highPU_BDTG/outputTMVATraining_lowPT_highPU/TMVATrainingResult_mpruned_QGLikelihood_PTBin_300_450_PU_39_100.root"), inputVariableOrMethodName = cms.vstring("BDT: M_{pruned}+Q/G"), JetPtBinOfTraining = cms.vdouble(300,450), PileUpBinOfTraining = cms.vdouble(0,39), ReducedName = cms.string("BDT_1")),
    cms.PSet( fileName = cms.string("PairVariablesTraining_lowPT_highPU_BDTG/outputTMVATraining_lowPT_highPU/TMVATrainingResult_mpruned_tau2tau1_PTBin_300_450_PU_39_100.root"), inputVariableOrMethodName = cms.vstring("BDT: M_{pruned}+#tau_{2}/#tau_{1}"), JetPtBinOfTraining = cms.vdouble(300,450), PileUpBinOfTraining = cms.vdouble(0,39), ReducedName = cms.string("BDT_2")),
    cms.PSet( fileName = cms.string("PairVariablesTraining_lowPT_highPU_BDTG/outputTMVATraining_lowPT_highPU/TMVATrainingResult_mtrim_tau2tau1_PTBin_300_450_PU_39_100.root"), inputVariableOrMethodName = cms.vstring("BDT: M_{trim}+#tau_{2}/#tau_{1}"), JetPtBinOfTraining = cms.vdouble(300,450), PileUpBinOfTraining = cms.vdouble(0,39), ReducedName = cms.string("BDT_3")),
    cms.PSet( fileName = cms.string("PairVariablesTraining_lowPT_highPU_BDTG/outputTMVATraining_lowPT_highPU/TMVATrainingResult_mtrim_QGLikelihood_comb_PTBin_300_450_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("BDT: M_{trim}+Q/G combo"), JetPtBinOfTraining = cms.vdouble(300,450), PileUpBinOfTraining = cms.vdouble(0,39), ReducedName = cms.string("BDT_4")),
    cms.PSet( fileName = cms.string("TripletVariablesTraining_lowPT_highPU_BDTG/outputTMVATraining_lowPT_highPU/TMVATrainingResult_mpruned_mtrim_tau2tau1_PTBin_300_450_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("BDT: M_{pruned}+M_{trim}+#tau_{2}/#tau_{1}"), JetPtBinOfTraining = cms.vdouble(300,450), PileUpBinOfTraining = cms.vdouble(0,39), ReducedName = cms.string("BDT_5")),
    cms.PSet( fileName = cms.string("TripletVariablesTraining_lowPT_highPU_BDTG/outputTMVATraining_lowPT_highPU/TMVATrainingResult_mpruned_tau2tau1_ECFbeta10_PTBin_300_450_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("BDT: M_{pruned}+#tau_{2}/#tau_{1}+C_{2}(#beta=1)"), JetPtBinOfTraining = cms.vdouble(300,450), PileUpBinOfTraining = cms.vdouble(0,39), ReducedName = cms.string("BDT_6")),
    cms.PSet( fileName = cms.string("AllVariablesTraining_lowPT_highPU_BDTG/outputTMVATraining_lowPT_highPU/TMVATrainingResult_allvariables_PTBin_300_450_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("all"), JetPtBinOfTraining = cms.vdouble(300,450), PileUpBinOfTraining = cms.vdouble(0,39), ReducedName = cms.string("allvariables")),
  )
)
