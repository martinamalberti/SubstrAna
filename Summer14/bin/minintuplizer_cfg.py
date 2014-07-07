import FWCore.ParameterSet.Config as cms

process = cms.Process("MiniNtuplizer")

process.Options = cms.PSet(

    maxEvents       = cms.int32(10),    # maximum events to  run

    jetR            = cms.double(0.8),  # basic clustering cone size
    jetPtCut        = cms.double(25.0), # pt cut on pf and Gen jets  
    jetAlgo = cms.string('antikt_algorithm'), # ex: antikt_algorithm, ak, AK, cambridge_algorithm, ca, CA
    doCMSSWJets     = cms.bool(False),  # analyze also basic cmssw reconstructed jets 

    puppiConfig     = cms.string("Puppi_cff.py"), # puppi configuration to run

    L1FastJetJEC    = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS170_V6_L1FastJet_AK7PF.txt"),
    L2RelativeJEC   = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS170_V6_L2Relative_AK7PF.txt"),
    L3AbsoluteJEC   = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS170_V6_L3Absolute_AK7PF.txt"),
    JECUncertainty  = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS170_V6_Uncertainty_AK7PF.txt"),
    L2L3ResidualJEC = cms.string(""), 

    L1FastJetJEC_CHS    = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS170_V6_L1FastJet_AK7PFchs.txt"),
    L2RelativeJEC_CHS   = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS170_V6_L2Relative_AK7PFchs.txt"),
    L3AbsoluteJEC_CHS   = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS170_V6_L3Absolute_AK7PFchs.txt"),
    JECUncertainty_CHS  = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS170_V6_Uncertainty_AK7PFchs.txt"),
    L2L3ResidualJEC_CHS = cms.string(""), 

    #mc matching
    DoMatchingToBoson = cms.bool(False), # use this to run on WW, ttbar...
    pdgIdBoson        = cms.int32(24),
    dRMatiching       = cms.double(0.3),

    
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
