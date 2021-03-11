import FWCore.ParameterSet.Config as cms

KinFitTester = cms.EDAnalyzer('KinFitTester',
                              patMuonTag = cms.InputTag("slimmedMuons"),
                              packedCandidateTag = cms.InputTag("packedPFCandidates"),
                              prunedGenParticleTag = cms.InputTag("prunedGenParticles"),
                              packedGenParticleTag = cms.InputTag("packedGenParticles"),
                              isMC = cms.bool(False),
                              ptMin = cms.double(1.0),
                              etaMax = cms.double(2.4)
                              )
