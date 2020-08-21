###https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile

from WMCore.Configuration import Configuration
config = Configuration() ###create a Configuration object

config.section_('General')###add a new section of type "General"
###General: In this section, the user specifies generic parameters about the request (e.g. request name).
config.General.workArea     = 'MC_hydjet_check3' ###fixed name for projects dir in my area
config.General.requestName  = 'MC_hydjet_check3'
config.General.transferLogs = True
config.General.transferOutputs = True

config.section_("Debug")
config.Debug.extraJDL = ["+CMS_ALLOW_OVERFLOW=False"]

################################

config.section_('JobType')###add a new section of type "JobType"
config.JobType.pluginName     = 'Analysis'
config.JobType.psetName       = 'ConfFile_cfg.py'#'run_PbPb_cfg.py'#'ConfFile_forCrab_DataHM_ntrackoffline80to104_cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.JobType.inputFiles = ['2018PbPb_Efficiency_GeneralTracks_MB_ChargePlus.root', '2018PbPb_Efficiency_GeneralTracks_MB_ChargeMinus.root', '2018PbPb_Efficiency_GeneralTracks_highPt.root', '2018PbPb_Efficiency_GeneralTracks_MB.root', '2018PbPb_Efficiency_PixelTracks.root']

config.section_('Data')###add a new section of type "Data"
###Data: This section contains all the parameters related to the data to be analyzed, 
config.Data.inputDataset      = '/MinBias_Hydjet_Drum5F_2018_5p02TeV/clindsey-RECODEBUG_20190625-5db5dfa073297cb96791f14c622e83e2/USER'#'/HydjetMB_10_3_1_patch2_UpdatedSmearing/abaty-crab_20181217_012503-5db5dfa073297cb96791f14c622e83e2/USER'#
config.Data.ignoreLocality = True

#'/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/davidlw-RecoSkim2016_pPb_V0_v1-2fc6918bc3c19ca88eae36cad5440243/USER'
#'/PAMinimumBias'+str(a)+'/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
config.Data.splitting         = 'FileBased'
#config.Data.unitsPerJob       = X ###files per job (but not impose)
#config.Data.totalUnits        = Y ###how many files to analyze
config.Data.unitsPerJob       = 10 ##High Multiplicity
config.Data.totalUnits        = 3958#263
config.Data.inputDBS          = 'phys03'
#config.Data.secondaryInputDataset = '/PAMinimumBias1/PARun2016C-PromptReco-v1/AOD'
#config.Data.useParent = True
config.Data.outLFNDirBase             = '/store/group/phys_heavyions/ddesouza/MC_hydjet_check/lemos/'

config.section_('Site')###add a new section of type "Site"
###Site: Grid site parameters are defined in this section, including the stage out information 
config.Site.whitelist = ['T2_CH_CERN','T2_US_MIT']
config.Site.storageSite       = 'T2_CH_CERN'
