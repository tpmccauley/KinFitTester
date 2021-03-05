import FWCore.ParameterSet.Config as cms

process = cms.Process("kinfit")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Services_cff')

process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source ("PoolSource",
    fileNames = cms.untracked.vstring(
    #'root://xrootd-cms.infn.it//store/data/Run2018A/Charmonium/MINIAOD/PromptReco-v2/000/316/239/00000/08BFAB4F-1359-E811-B01F-FA163E5285EC.root'
  
    'root://xrootd-cms.infn.it//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/B8DB27C5-57AC-E711-9112-FA163E7B5756.root'

    )

)

process.load('KinFitTester.KinFitTester.KinFitTester_cfi')
process.KinFitTester.patMuonTag = cms.InputTag('slimmedMuons')
process.KinFitTester.isMC = cms.bool(True)
process.KinFitTester.ptMin = cms.double(0.8)

process.p = cms.Path(
    process.KinFitTester
)

process.schedule = cms.Schedule(process.p)

