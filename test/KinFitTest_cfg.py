import FWCore.ParameterSet.Config as cms

process = cms.Process("kinfit")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )

process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')

'''
Use this for data:
'''
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

'''
Use this for MC:
Important!
CMSSW will not complain if you use the GlobalTag for data but
you will receive invalid kinematic fits.
This will be all down to consequently receiving magnetic field
values of 0 which will eventually produce nan values down the line.
Note to self: this was tracked down the hard way in
TrackKinematicStatePropagator.
'''
process.GlobalTag.globaltag = '106X_mc2017_realistic_v6'

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
  
    '/store/mc/RunIISummer20UL17MiniAOD/BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/106X_mc2017_realistic_v6-v1/00000/014122A9-ED46-6849-A0DC-6BC69F563EE5.root',
    #'/store/mc/RunIISummer20UL17MiniAOD/BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/106X_mc2017_realistic_v6-v1/00000/0261FDE9-2436-5F47-983A-947D28D6A7E3.root',
    #'/store/mc/RunIISummer20UL17MiniAOD/BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/106X_mc2017_realistic_v6-v1/00000/06D4D2CC-88BF-5C4A-A205-78497AB1B5B8.root',
    #'/store/mc/RunIISummer20UL17MiniAOD/BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/106X_mc2017_realistic_v6-v1/00000/084658A9-F3EA-804C-ABE4-7489CAD345D3.root',
    #'/store/mc/RunIISummer20UL17MiniAOD/BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/106X_mc2017_realistic_v6-v1/00000/0AD9B380-034C-AC42-8948-CD5C84F4CB57.root',
  
    #'/store/mc/RunIISummer20UL16RECO/BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/106X_mcRun2_asymptotic_v13-v2/280000/09E16B4E-95C1-934E-93BF-D0F54B999D46.root',
    #'/store/mc/RunIISummer20UL16RECO/BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/106X_mcRun2_asymptotic_v13-v2/280000/0A28A937-3856-4544-AB30-40CD14B0FCC1.root',
    #'/store/mc/RunIISummer20UL16RECO/BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/106X_mcRun2_asymptotic_v13-v2/280000/1F63348D-C859-3245-9DCA-1429A07B0C7E.root',
    #'/store/mc/RunIISummer20UL16RECO/BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/106X_mcRun2_asymptotic_v13-v2/280000/23D5DCCE-FA45-1844-A8DC-C8C840C971EE.root',
    #'/store/mc/RunIISummer20UL16RECO/BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/106X_mcRun2_asymptotic_v13-v2/280000/24F34997-8CFF-3F44-9900-48CCA2E457AA.root',

    #'/store/data/Run2016G/Charmonium/AOD/21Feb2020_UL2016-v1/20000/140F7D26-0ECD-4F4D-A597-C144F0D7112D.root',
    #'/store/data/Run2016G/Charmonium/AOD/21Feb2020_UL2016-v1/20000/74F1A473-848E-5342-BEC7-2160362ADB0A.root',
    #'/store/data/Run2016G/Charmonium/AOD/21Feb2020_UL2016-v1/20000/CB113922-6FA9-2C4F-8104-FE5AA8716394.root',
    #'/store/data/Run2016G/Charmonium/AOD/21Feb2020_UL2016-v1/210000/22215BCE-FF56-604E-AB9B-A74120F6DD74.root',
    #'/store/data/Run2016G/Charmonium/AOD/21Feb2020_UL2016-v1/210000/50ECD1B8-BA6A-414D-9774-E7E4362F67F4.root',
    ),

    # This is for AOD
    inputCommands=cms.untracked.vstring('keep *', 'drop *_muonReducedTrackExtras_*_*')


)

process.load('KinFitTester.KinFitTester.KinFitTester_cfi')
process.KinFitTester.patMuonTag = cms.InputTag('slimmedMuons')
process.KinFitTester.isMC = cms.bool(True)
process.KinFitTester.isMiniAOD = cms.bool(True)
process.KinFitTester.ptMin = cms.double(0.8)

process.p = cms.Path(
    process.KinFitTester
)

process.schedule = cms.Schedule(process.p)

