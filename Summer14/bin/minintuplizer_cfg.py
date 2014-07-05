import FWCore.ParameterSet.Config as cms

process = cms.Process("MiniNtuplizer")

process.Options = cms.PSet(
    maxEvents       = cms.int32(10),
    jetR            = cms.double(0.8),
    jetPtCut        = cms.double(25.0),
   
    jetAlgo = cms.string('antikt_algorithm'), # ex: antikt_algorithm, ak, AK, cambridge_algorithm, ca, CA
    doCMSSWJets     = cms.bool(False),
    puppiConfig     = cms.string("Puppi_cff.py"),
    L1FastJetJEC    = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS162_V5_L1FastJet_AK7PF.txt"),
    L2RelativeJEC   = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS162_V5_L2Relative_AK7PF.txt"),
    L3AbsoluteJEC   = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS162_V5_L3Absolute_AK7PF.txt"),
    L2L3ResidualJEC = cms.string(""), 
    JECUncertainty  = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS162_V5_Uncertainty_AK7PF.txt"),
    dRMatiching     = cms.double(0.3),
    #mc matching
    DoMatchingToBoson = cms.bool(False), # use this to run on WW, ttbar...
    pdgIdBoson        = cms.int32(24),
    
    #softdrop
    beta = cms.double(2),
    symmetry_cut = cms.double(0.1),
    R0 = cms.double(1.),
    
    # trimming
    R_trimming = cms.double(0.2),
    PtFraction = cms.double(0.05),
    trimAlgo = cms.string('kt_algorithm'),
  
    #pruning
    z_cut = cms.double(0.1),
    R_Cut = cms.double(0.5),
    R_jet_def_pruning = cms.double(0.9),
    pruneAlgo = cms.string('cambridge_algorithm')

)
