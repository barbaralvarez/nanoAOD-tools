#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from framework.postprocessor import PostProcessor
import json
#from framework.crabhelper import inputFiles,runsAndLumis

isData = sys.argv[-1] == 'data'
dataverbose = 'Is data!!' if isData else 'is MC!!'
print dataverbose

jsonfile = ''
year = 17

### INPUT FILE
#filepath = ['/afs/cern.ch/work/j/jrgonzal/public/pruebaNanoAOD/8EAB6B64-9210-E811-B19D-FA163E759AE3.root']
#filepath = ['root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAOD/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/40000/2CE738F9-C212-E811-BD0E-EC0D9A8222CE.root']
#filepath = ['/nfs/fanae/user/juanr/nanoAOD/CMSSW_9_4_6/src/PhysicsTools/NanoAODTools/python/postprocessing/skimtree.root']
if isData:
  #filepath = ['root://cms-xrd-global.cern.ch//store/data/Run2017D/SingleMuon/NANOAOD/Nano14Dec2018-v1/80000/0726FD7F-D51D-4348-8B48-FA2B7D2B403D.root']
  if     year == 18: 
    jsonfile = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-322057_13TeV_PromptReco_Collisions18_JSON.txt'
    filepath = ['root://cms-xrd-global.cern.ch//store/data/Run2018A/SingleMuon/NANOAOD/Nano14Dec2018-v1/40000/FCA30DD8-1CDE-DE47-B1FB-F6192CE78748.root']
  elif   year == 17:
    #jsonfile = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
    filepath = ['root://cms-xrd-global.cern.ch//store/data/Run2017D/SingleMuon/NANOAOD/Nano14Dec2018-v1/80000/3E9394FA-EC63-114B-8AB1-18BD46BC5D43.root']
    #filepath = ['root://cms-xrd-global.cern.ch//store/data/Run2017D/SingleElectron/NANOAOD/Nano14Dec2018-v1/20000/123128F7-CB46-1043-935A-03C2A361E438.root']
  elif   year == 16: jsonfile = ''
  #jsn = open(jsonfile, 'r')
  #jsonInput = json.loads(jsn.read())
  #jsn.close()
  #jsn = open(jsonfile, 'r')
  #jsonInput = json.loads(jsn.read())
  #jsn.close()
  #print 'Using json: ', jsonfile

else: 
  if   year == 17:
    filepath = ['root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAODv4/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/20000/CCF670CF-51AD-FA46-8EEF-A10CE19C3474.root']
  elif year == 18:
    filepath = ['root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18NanoAODv4/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano14Dec2018_102X_upgrade2018_realistic_v16-v2/110000/FCD6697F-D52A-2D4F-8EE9-594503D09DCA.root']

### OUTPUT
outdir = '.'

### SKIM 
#cut = 'nMuon >= 2 && Muon_pt[0] > 25 && Muon_pt[1] >= 12 && Jet_pt > 200'#(nElectron + nMuon) >= 2 && Jet_pt > 200' # nGenDressedLepton >= 2
cut = 'nMuon >= 2 && Muon_pt[0] > 25 && Muon_pt[1] >= 12'#(nElectron + nMuon) >= 2 && Jet_pt > 200' # nGenDressedLepton >= 2
#cut = 'nElectron >= 2 && Electron_pt[0] > 25 && Electron_pt[1] >= 12 '

### SLIM FIL
slimfile = "SlimFileTnPmuon.txt"
#slimfile = "SlimFileTnPele.txt"

### MODULES
### Include modules to compute derivate quantities or calculate uncertainties
from modules.jme.jetmetUncertainties import *
from modules.common.puWeightProducer import *
from modules.skimNRecoLeps import *
from modules.addTnPvarMuon import *
from modules.addTnPvarEle import *
#from modules.addSUSYvar import *
#mod = [puAutoWeight(), skimRecoLeps(), addSUSYvarsMC()] # countHistogramsProducer(), jetmetUncertainties2017All()

mod = [puAutoWeight()] if not isData else []# countHistogramsProducer(), jetmetUncertainties2017All()
if isData:
  if   year == 17: mod.append(addTnPMuon17data())
  #if   year == 17: mod.append(addTnPEle17data())
  elif year == 18: mod.append(addTnPMuonForMoriond18data())
else:
  if   year == 17: mod.append(addTnPMuon17())
  #if   year == 17: mod.append(addTnPEle17())
  elif year == 18: mod.append(addTnPMuonForMoriond18())

p=  PostProcessor(outdir,filepath,cut,slimfile,mod,provenance=True,outputbranchsel=slimfile) if jsonfile == '' else PostProcessor(outdir,filepath,cut,slimfile,mod,provenance=True,outputbranchsel=slimfile,jsonInput = jsonInput)
p.run()
