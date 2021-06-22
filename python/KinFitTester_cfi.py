import FWCore.ParameterSet.Config as cms

KinFitTester = cms.EDAnalyzer('KinFitTester',
                              recoMuonTag = cms.InputTag("muons"),
                              patMuonTag = cms.InputTag("slimmedMuons"),
                              packedCandidateTag = cms.InputTag("packedPFCandidates"),
                              prunedGenParticleTag = cms.InputTag("prunedGenParticles"),
                              packedGenParticleTag = cms.InputTag("packedGenParticles"),
                              isMC = cms.bool(False),
                              isAOD = cms.bool(False),
                              isMiniAOD = cms.bool(False),
                              ptMin = cms.double(1.0),
                              etaMax = cms.double(2.4)
                              )
