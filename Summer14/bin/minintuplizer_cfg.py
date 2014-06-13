import FWCore.ParameterSet.Config as cms

process = cms.Process("MiniNtuplizer2")

process.Options = cms.PSet(
    maxEvents      = cms.int32(10),
    jetR           = cms.double(0.8),
    doCMSSWJets    = cms.bool(False),
    puppiConfig    = cms.string("Puppi_cff.py")
)
