#------------------------------------------------------------------------------------
# Imports
#------------------------------------------------------------------------------------
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
import FWCore.ParameterSet.VarParsing as VarParsing

#------------------------------------------------------------------------------------
# Declare the process and input variables
#------------------------------------------------------------------------------------
#process = cms.Process('NOISE',eras.Run2_50ns)#for 50ns 13 TeV data
#process = cms.Process('NOISE',eras.Run2_25ns)#for 25ns 13 TeV data
options = VarParsing.VarParsing ('analysis')
process = cms.Process("Trees",eras.Phase2) 

##
## Setup command line options
##
options.register ('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "no of skipped events")

##
## Default
##
options.maxEvents = -1 # means all events
options.skipEvents = 0 # default is 0.

##
## get and parse the command line arguments
##
options.parseArguments()

#
# Dataset e.g.
# dasgoclient --query 'dataset dataset=/RelValTTbar_13/CMSSW_10_2_0_pre3-*realistic*/GEN-SIM-RECO'                 
# dasgoclient --query 'file dataset=/RelValTTbar_13/CMSSW_10_2_0_pre3-101X_upgrade2018_realistic_v7-v1/GEN-SIM-RECO'
#
# TTbar sample
#
# RECO
options.inputFiles = [
     'file:/store/relval/CMSSW_11_0_0_pre7/RelValSinglePiFlatPt0p7To10/MINIAODSIM/110X_mcRun4_realistic_v1_2026D41noPU-v1/10000/9D938200-8C7C-6940-862D-067B7DCB4980.root'
#    'file:/eos/uscms/store/user/lpcjme/PFForwardTrack/SinglePiPt50_3p45_3p54_n1000/SinglePiPt50_3p45_3p54_n1000_step3.root',
#
# noTrackLink samples 
#
#    'file:/afs/cern.ch/work/g/gagarwal/public/SinglePiPt10Eta3p45_3p54_n1000_noTrackLink_step3.root',
#    'file:/afs/cern.ch/work/g/gagarwal/public/SinglePiMinusPt10Eta3p45_3p54_n1000_noTrackLink_step3.root',     
#    'file:/afs/cern.ch/work/g/gagarwal/public/SinglePiPt50Eta2p8_4p2_n500_noTrackLink_step3.root',
#    'file:/afs/cern.ch/work/g/gagarwal/public/SinglePiPt50Eta3p45_3p54_n1000_noTrackLink_step3.root',
#
# TrackLink samples 
#
#    'file:/afs/cern.ch/work/g/gagarwal/public/SinglePiPt50Eta3p45_3p54_n800_step3.root',
#
#    'file:/eos/uscms/store/user/lpcjme/PFForwardTrack/SinglePiPt5_3p45_3p54_n1000/SinglePiPt5_3p45_3p54_n1000_step3.root',
#
#
#
#   
######    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_1.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_2.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_3.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_4.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_5.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_6.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_7.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_8.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_9.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_10.root',
######    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_11.root',
######    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_12.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_13.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_14.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_15.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_16.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_17.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_18.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_19.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_20.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_21.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt5Eta3p45_3p54/SinglePiPt5Eta3p45_3p54_step3_22.root',
#
#
#
####    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_1.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_2.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_3.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_4.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_5.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_6.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_7.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_8.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_9.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_10.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_11.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_12.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_13.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_14.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_15.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_16.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_17.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_18.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_19.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_20.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_21.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_22.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_23.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_24.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_25.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_26.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_27.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_28.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_29.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_30.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiMinusPt5Eta3p45_3p54/SinglePiMinusPt5Eta3p45_3p54_step3_31.root',
#  FakeTracks
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt10Eta3p45_3p54_n1000_FakeTracks.root',
#
#   PU sample
#    'file:/afs/cern.ch/work/g/gagarwal/public/step3_SinglePiPt50Eta3p45_3p54_n10_PU10.root'
#
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_1.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_2.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_3.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_4.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_5.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_6.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_7.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_8.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_9.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_10.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_11.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_12.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_13.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_14.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_15.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_16.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_17.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_18.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_19.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_20.root',
#    'file:/afs/cern.ch/user/i/iashvili/public/PFLOW/SinglePiPt50Eta3p45_3p54_PU10/SinglePiPt50Eta3p45_3p54_PU10_21.root',
#
]

options.outputFile = 'hist_PFTrackHFAnalyzer.root'
#
#
#
print("maxEvents: ", options.maxEvents)
print("inputFiles: ", options.inputFiles)
print("outputFile: ", options.outputFile)

#------------------------------------------------------------------------------------
# Get and parse the command line arguments
#------------------------------------------------------------------------------------
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

readFiles = cms.untracked.vstring()
process.source = cms.Source("PoolSource",
    fileNames = readFiles,
    skipEvents = cms.untracked.uint32(options.skipEvents) # default is 0.
)
readFiles.extend(options.inputFiles);

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string(options.outputFile)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    Rethrow = cms.untracked.vstring("ProductNotFound"), # make this exception fatal
    fileMode  =  cms.untracked.string('NOMERGE') # no ordering needed, but calls endRun/beginRun etc. at file boundaries
)

#------------------------------------------------------------------------------------
# import of standard configurations
#------------------------------------------------------------------------------------
# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D41Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1

#------------------------------------------------------------------------------------
# Set up our analyzer
#------------------------------------------------------------------------------------

#process.load("PFAnalysis.PFAnalyzers.PFTrackHFAnalyzer_cfi")
eventN = [1,2,3,4,5,6,7,8,9,10]
process.pfTrackHFAnalyzer = cms.EDAnalyzer("PFTrackHFAnalyzer",
                                           source_genpars = cms.untracked.InputTag('genParticles', ''),
                                           source_calopars = cms.untracked.InputTag('mix', 'MergedCaloTruth'),
                                           source_vertices = cms.untracked.InputTag('offlinePrimaryVertices', ''),
                                           source_pfcands = cms.untracked.InputTag('particleFlow', ''),
                                           source_pfclustersHF = cms.untracked.InputTag('particleFlowClusterHF', ''),
                                           source_pfrechitsHF = cms.untracked.InputTag('particleFlowRecHitHF', ''),
                                           source_pftracks = cms.untracked.InputTag('pfTrack', ''),
                                           source_tracks = cms.untracked.InputTag('generalTracks', ''),
                                           source_hfrechits = cms.untracked.InputTag('hfreco', ''),
                                           source_pileup = cms.untracked.InputTag('addPileupInfo', ''),
                                           debug = cms.untracked.bool(True),
                                           debugRecHit = cms.untracked.bool(False),
                                           EventsToScan = cms.vint32(eventN),
                                           ptlow  = cms.double(0.),
                                           pthigh = cms.double(1000.),    
                                           etalow  = cms.double(0.),    
                                           etahigh = cms.double(10.),

)

#------------------------------------------------------------------------------------
# Specify Global Tag
#------------------------------------------------------------------------------------
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

#------------------------------------------------------------------------------------
# Sequence definition
#------------------------------------------------------------------------------------
process.ana_step = cms.Sequence(
    process.pfTrackHFAnalyzer
)

#-----------------------------------------------------------------------------------
# Path and EndPath definitions
#-----------------------------------------------------------------------------------
process.preparation = cms.Path(
    process.ana_step
)
