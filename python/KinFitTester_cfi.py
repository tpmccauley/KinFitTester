import FWCore.ParameterSet.Config as cms

KinFitTester = cms.EDAnalyzer('KinFitTester',
                              patMuonTag = cms.InputTag("slimmedMuons")
                              )
