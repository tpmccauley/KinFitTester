import FWCore.ParameterSet.Config as cms

process = cms.Process("kinfit")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(400) )

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


    #'root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/016238DC-0238-7B47-BBAD-C1BC2E57E6A2.root',
    'root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/03D0CB2A-4355-674A-A94F-9587733F06BE.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/09745F94-76B9-FB47-8475-53C2FF4F9ADC.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/0BEEEC95-B25D-FF43-8C6A-01A8F3F045A0.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/0D293CCC-AD78-5641-8E27-8F24B3E8CE3C.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/0DA87150-D1AD-9B4E-B4CA-6CA4E3E496A0.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/0F73A8A1-AD94-8D4F-8DBA-7E4E26AEEACF.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/0F9F8225-BB31-E646-94A1-ED333500C979.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/1178049C-578C-7249-8341-296EBBFFAEE3.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/131F167D-1F8A-6E4B-BA35-2687F57B520D.root',
   
 
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/14422CC4-3BAD-E711-85D2-FA163E3BCC72.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/3A0A3362-36AC-E711-BBED-FA163ED8E79E.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/4469C355-B2AC-E711-BE8E-02163E012D69.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/4CD344EC-3BAD-E711-BD95-FA163E3201B4.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/549DB69D-3CAD-E711-849F-0CC47A7C3422.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/5AC2A21E-64AC-E711-B29B-FA163E7625E2.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/6E47188D-3BAD-E711-8B3F-0CC47A6C1874.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/7A22593E-4AAC-E711-84F7-02163E0176A5.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/80125A7C-7EAC-E711-8CDA-FA163ECC515A.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/8E5FA699-64AC-E711-9772-FA163EE771A0.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/949F3A31-6BAC-E711-B38C-FA163EF027B6.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/A0126530-5EAC-E711-98A1-02163E0164E3.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/A60BFA1C-3CAD-E711-AE31-008CFAC91B60.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/B8DB27C5-57AC-E711-9112-FA163E7B5756.root',
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

