import FWCore.ParameterSet.Config as cms

process = cms.Process("OptimizeSelections")

process.Options = cms.PSet(


#   InputVariableList = cms.vstring("mraw[0]","mclean[0]","mtrim_Rtrim_020_Ptfrac_005[0]","mtrim_Rtrim_010_Ptfrac_003[0]","mtrimsafe_Rtrim_020_Ptfrac_005[0]","mtrimsafe_Rtrim_010_Ptfrac_003[0]"),
# InputVariableList = cms.vstring("mpruned_zcut_010_R_cut_050[0]","mpruned_zcut_010_R_cut_075[0]","mprunedsafe_zcut_010_R_cut_050[0]","mprunedsafe_zcut_010_R_cut_075[0]","msoftdrop_beta20[0]","msoftdrop_beta00[0]","msoftdropsafe_beta20[0]","msoftdropsafe_beta00[0]"),
#InputVariableList = cms.vstring("tau1[0]","tau2[0]","tau2[0]/tau1[0]","QGLikelihood_pr_zcut_010_R_cut_050[0]","QGLikelihood_pr_zcut_010_R_cut_075[0]"),
#InputVariableList = cms.vstring("tau1_pr[0]","tau2_pr[0]","tau2_pr[0]/tau1_pr[0]","tau2_softdrop[0]","tau1_softdrop[0]","tau2_softdrop[0]/tau1_softdrop[0]"),
InputVariableList = cms.vstring("Qjets[0]","charge_k05[0]","charge_k07[0]","charge_k10[0]","ecf_beta_05[0]","ecf_beta_10[0]","ecf_beta_15[0]","ecf_beta_20[0]"),
#InputVariableList = cms.vstring("tau2[0]/tau1[0]","tau2_pr[0]/tau1_pr[0]","tau2_softdrop[0]/tau1_softdrop[0]"),


   ## Spectator variable to be used in the training
   InputSpectatorList = cms.vstring("pt","npu"), 

   ## Name of the tree on which perform th training
   TreeName           = cms.string("gen"), 

   ## Label to be used in the output file creation
   Label              = cms.string("GenTraining"),

   ## Lepton Type: Muon, Electron ,EleMu and Jets (fully hadronic)
   LeptonType         = cms.string("Jets"),
  
   ## if you want to run directly the TMVA macros
   isPrintResultwithTMVA  = cms.bool(False),
  
   ## Name for ggH signal --> match what is written in the sample list file
   SignalggHName      = cms.string("RSGWW1000"),

   ## Name for qqH signal
   SignalqqHName      = cms.string("qqHx600"),
  
   ## 0 use both ggH and qqH as signal, 1 use only ggH as signal, 2 use only qqH as signal
   useTypeOfSignal    = cms.uint32(0),

   ## string which use tree branches to re-weight the events
   EventWeight        = cms.string(""),

   ## string which is used in the TMVATraining class to define a cut to be applied on the events
   PreselectionCutType = cms.string("basicJetsCutCSA14"),

   ## List of MVA method to be used in the training
   UseMethodName       = cms.vstring("CutsSA"),

   ## W-jet pt bin for training
   JetPtBinOfTraining  = cms.vdouble(500,600),

   ## In time pile-up bin for the training
   PileUpBinOfTraining = cms.vdouble(0,100),

   ## Look at the variable list and train each of them with a rectangular cut skipping other non linear MVA possibilities
   TrainEachVariable   = cms.bool(True),
   
   ## output directory for root and weight file
   outputFileDirectory  = cms.string("output/outputTMVATraining/"),

   ## output file name
   outputFileName       = cms.string("TMVATrainingResult"),

   ## input directory where background trees are placed
   InputBackgroundParam = cms.VPSet(

#     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-300to470_Tune4C_13TeV_pythia8_Spring14/outtree_0.root"), ReducedName = cms.string("qcd300to470"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
#     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-300to470_Tune4C_13TeV_pythia8_Spring14/outtree_1.root"), ReducedName = cms.string("qcd300to470"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
#     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-300to470_Tune4C_13TeV_pythia8_Spring14/outtree_2.root"), ReducedName = cms.string("qcd300to470"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
#     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-300to470_Tune4C_13TeV_pythia8_Spring14/outtree_3.root"), ReducedName = cms.string("qcd300to470"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
#     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-300to470_Tune4C_13TeV_pythia8_Spring14/outtree_4.root"), ReducedName = cms.string("qcd300to470"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
#     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-300to470_Tune4C_13TeV_pythia8_Spring14/outtree_5.root"), ReducedName = cms.string("qcd300to470"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
#     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-300to470_Tune4C_13TeV_pythia8_Spring14/outtree_6.root"), ReducedName = cms.string("qcd300to470"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
#     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-300to470_Tune4C_13TeV_pythia8_Spring14/outtree_7.root"), ReducedName = cms.string("qcd300to470"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
#     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-300to470_Tune4C_13TeV_pythia8_Spring14/outtree_8.root"), ReducedName = cms.string("qcd300to470"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
#     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-300to470_Tune4C_13TeV_pythia8_Spring14/outtree_9.root"), ReducedName = cms.string("qcd300to470"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
#     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-300to470_Tune4C_13TeV_pythia8_Spring14/outtree_10.root"), ReducedName = cms.string("qcd300to470"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),

     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_0.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_1.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_2.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_3.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_4.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_5.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_6.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_7.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_8.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_9.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_11.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_12.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_13.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_14.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_16.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_17.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_18.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_19.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_20.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_21.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_22.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_23.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_24.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_25.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_26.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_27.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_28.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_29.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_30.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_31.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_32.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_33.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_34.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_35.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_36.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_37.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_38.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_39.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_40.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_41.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_42.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_43.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_44.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_45.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_46.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_47.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_48.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_49.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_50.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_51.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_52.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_53.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_54.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_55.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_56.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_57.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_58.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_59.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_60.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
   
   ),

   InputSignalParam = cms.VPSet(
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_0.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_1.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_2.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_3.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_4.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_5.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_6.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_7.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_8.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_9.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_10.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_11.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_11.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_12.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_13.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_14.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_15.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_16.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_17.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_18.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_19.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_20.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_21.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_22.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_23.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_24.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_25.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_26.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_27.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_28.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_29.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_30.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_31.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_32.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_33.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_34.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_35.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_36.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_37.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_38.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_39.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_40.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_41.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_42.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_43.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_44.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_45.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_46.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_47.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_48.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_49.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_50.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),   
     cms.PSet(inputFileName = cms.string("eos/cms/store/user/rgerosa/MiniNtuple_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_51.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),   
  )

)
