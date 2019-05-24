#!/usr/bin/env python
import os, sys, json
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import * 

#this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles,runsAndLumis

### SKIM 
#cut = 'Jet_pt > 200 && (nElectron + nMuon) >= 2 && nGenDressedLepton >= 2'
cut = '(nElectron + nMuon) >= 2'

### SLIM FILE
slimfile = "SlimFile.txt"

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.skimNRecoLeps import *
from PhysicsTools.NanoAODTools.postprocessing.modules.addTnPvarMuon import *
from PhysicsTools.NanoAODTools.postprocessing.modules.addTnPvarEle import *

isData   = 'data' in sys.argv[-1] or 'Data' in sys.argv[-1]
doTnP    = 'TnP'  in sys.argv[-1]
doMuon   = 'muo'  in sys.argv[-1]
doEle    = 'ele'  in sys.argv[-1]
doJECunc = 'JEC'  in sys.argv[-1]
if         '18' in sys.argv[-1] : year = 18
elif       '16' in sys.argv[-1] : year = 16
else                            : year = 17

### Json file
jsonfile = runsAndLumis()

#if isData:
#  if     year == 18: jsonfile = 'Cert_314472-322057_13TeV_PromptReco_Collisions18_JSON.txt'
#  elif   year == 17: jsonfile = 'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt' 
#  elif   year == 16: jsonfile = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
  

mod = []
if not isData: 
  if   year == 16:  mod.append(puWeight2016())
  elif year == 17:  mod.append(puAutoWeight())
  elif year == 18:  mod.append(puAutoWeight())
  else           :  mod.append(puAutoWeight())

if doTnP:
  if doMuon:
    if isData:
      if   year == 17: mod.append(addTnPMuon17data())
      elif year == 16: mod.append(addTnPMuon16data())
      elif year == 18: mod.append(addTnPMuonForMoriond18data())
    else:
      if   year == 17: mod.append(addTnPMuon17())
      elif year == 16: mod.append(addTnPMuon16())
      elif year == 18: mod.append(addTnPMuonForMoriond18())
    cut = 'nMuon >= 2 && Muon_pt[0] > 25 && Muon_pt[1] >= 5'
    slimfile = "SlimFileTnPmuon.txt"
  elif doEle:
    if isData:
      if   year == 17: mod.append(addTnPEle17data())
      elif year == 16: mod.append(addTnPEle16data())
      elif year == 18: mod.append(addTnPEleForMoriond18data())
    else:
      if   year == 17: mod.append(addTnPEle17())
      elif year == 16: mod.append(addTnPEle16())
      elif year == 18: mod.append(addTnPEleForMoriond18())
    cut = 'nElectron >= 2 && Electron_pt[0] > 25 && Electron_pt[1] >= 5'
    slimfile = "SlimFileTnPele.txt"
  else:
    print 'Error No electron or muon selection'


  #cut = 'nMuon >= 2 && Muon_pt[0] > 25 && Muon_pt[1] >= 12'
else:
  mod.append(skimRecoLeps())
  if doJECunc and not isData: mod.append(jetmetUncertainties2017All())

print '>> Slim file: ', slimfile
print '>> cut: ', cut
print '>> ' + ('Is data!' if isData else 'Is MC!')
print '>> ' + ('Creating a TnP Tree' if doTnP else 'Creating a skimmed nanoAOD file')
if doJECunc: print '>> Adding JEC uncertainties'

#mod = [puAutoWeight(),jetmetUncertainties2017All(), skimRecoLeps()]
p=PostProcessor(".",inputFiles(),cut,slimfile,mod,provenance=True,fwkJobReport=True,jsonInput=jsonfile,outputbranchsel=slimfile) #jsonInput
p.run()

print "DONE"
os.system("ls -lR")
