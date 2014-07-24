
import FWCore.ParameterSet.Config as cms

process = cms.Process("OptimizeSelections")

process.Options = cms.PSet(


   ## Single Variables secttion
   InputVariableList = cms.vstring(INPUFILELIST),

   ## Spectator variable to be used in the training
   InputSpectatorList = cms.vstring("pt","npu"), 

   ## Name of the tree on which perform th training
   TreeName           = cms.string("chs"), 

   ## Label to be used in the output file creation
   Label              = cms.string(LABELNAME),

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
   UseMethodName       = cms.vstring("METHOD"),

   ## W-jet pt bin for training
   JetPtBinOfTraining  = cms.vdouble(475,600),

   ## In time pile-up bin for the training
   PileUpBinOfTraining = cms.vdouble(0,39),

   ## Look at the variable list and train each of them with a rectangular cut skipping other non linear MVA possibilities
   TrainEachVariable   = cms.bool(False),
   
   ## output directory for root and weight file
   outputFileDirectory  = cms.string("OUTPUTDIR"),

   ## output file name
   outputFileName       = cms.string("TMVATrainingResult"),

   ## input directory where background trees are placed
   InputBackgroundParam = cms.VPSet(

     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_0.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_1.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_2.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_3.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_4.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_5.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_6.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_7.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_8.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_9.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_11.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_12.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_13.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_14.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_16.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_17.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_18.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_19.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_20.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_21.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_22.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_23.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_24.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_25.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_26.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_27.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_28.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_29.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_30.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_31.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_32.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_33.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_34.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_35.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_36.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_37.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_38.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_39.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_40.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_41.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_42.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_43.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_44.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_45.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_46.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_47.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_48.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_49.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_50.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_51.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_52.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_53.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_54.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_55.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_56.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_57.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_58.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_59.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_60.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_61.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_62.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_63.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_64.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_65.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_66.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_67.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_68.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_69.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_70.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_71.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_72.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_73.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_74.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_75.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_76.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_77.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_78.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_79.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_80.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/QCD_Pt-470to600_Tune4C_13TeV_pythia8_Spring14/outtree_81.root"), ReducedName = cms.string("qcd470to600"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
    
   ),

   InputSignalParam = cms.VPSet(
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_0.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_1.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_2.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_3.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_4.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_5.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_6.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_7.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_8.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_9.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_10.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_11.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_11.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_12.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_13.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_14.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_15.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_16.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_17.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_18.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_19.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_20.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_21.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_22.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_23.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_24.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_25.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_26.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_27.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_28.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_29.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_30.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_31.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_32.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_33.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_34.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_35.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_36.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_37.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_38.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_39.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_40.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_41.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_42.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_43.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_44.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_45.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_46.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_47.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_48.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_49.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_50.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),   
     cms.PSet(inputFileName = cms.string("eos/cms/store/caf/user/rgerosa/MiniNtuples_csa14/RSGravToWW_kMpl01_M-1000_Tune4C_13TeV-pythia8_Spring14/outtree_51.root"), ReducedName = cms.string("RSGWW1000"), CrossSection = cms.double(1.), NumberEntriesBefore = cms.int32(1)),   
  )

)
