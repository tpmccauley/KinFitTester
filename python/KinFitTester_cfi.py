import FWCore.ParameterSet.Config as cms

KinFitTester = cms.EDAnalyzer('KinFitTester',
                              patMuonTag = cms.InputTag("slimmedMuons"),
                              patGenParticleTag = cms.InputTag("packedGenParticles"),
                              ptMin = cms.double(1.0),
                              etaMax = cms.double(2.4)
                              )
