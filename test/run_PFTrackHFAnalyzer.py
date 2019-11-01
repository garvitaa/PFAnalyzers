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
options.maxEvents = -1 #-1 # means all events
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

      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_99.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_98.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_97.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_96.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_95.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_94.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_93.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_92.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_91.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_9.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_88.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_87.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_86.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_85.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_84.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_82.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_81.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_80.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_8.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_77.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_76.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_75.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_74.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_73.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_71.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_70.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_7.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_69.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_68.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_67.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_65.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_64.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_63.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_60.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_6.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_59.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_58.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_56.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_55.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_54.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_53.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_52.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_51.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_5.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_49.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_48.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_47.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_46.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_45.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_44.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_43.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_42.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_41.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_40.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_4.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_39.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_38.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_37.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_36.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_35.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_34.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_33.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_32.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_31.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_30.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_3.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_29.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_28.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_27.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_25.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_24.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_23.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_22.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_21.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_20.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_2.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_19.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_18.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_17.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_16.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_15.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_14.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_13.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_12.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_11.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_100.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_10.root',
      'root://cmseos.fnal.gov//store/group/lpcjme/PFForwardTrack/SinglePiPt5Eta3p45_3p54_2026D41/step3/191018_160332/0000/SinglePiPt5Eta3p45_3p54_step3_1.root'
]

options.outputFile = 'hist_PFTrackHFAnalyzer_SinglePiPt5Eta3p45_3p54.root'
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
