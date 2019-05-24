import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import TLorentzVector

import os
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.JetReCalibrator import JetReCalibrator

doCalculateJecLepAwareFromNanoAOD = True

class addTnPvarEle(Module):
    def __init__(self, isdata = False, year = 17, recalibjets = '', era = ''):
      print 'Initializing addTnPvarEle...'
      #self.kProbePt = 12
      self.kProbePt = 7
      self.kTagPt   = 29
      self.kTagIso  = 0.20
      self.kMaxMass = 140
      self.kMinMass = 60
      self.isData = isdata
      self.year = year
      self.era = era
      self.i = 0
      self.filenameJECrecal = recalibjets
      self.filenameJEC = recalibjets
      if self.filenameJEC == '': self.filenameJEC = self.GetFileNameJEC(self.isData, self.filenameJEC, self.year, self.era)
      if not doCalculateJecLepAwareFromNanoAOD: self.jetReCalibrator = self.OpenJECcalibrator()

    def GetFileNameJEC(self, isdata, version = '', year = '', era = ''):
      f = version
      if f == '':
        if   year == 16: f = 'Summer16_23Sep2016V4'
        elif year == 17: f = 'Fall17_17Nov2017_V32'
        elif year == 18: f = 'Autumn18_V3'
      if isdata: f+= '_DATA'
      else:      f+= '_MC'
      return f

    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Tag_pt",   "F")
        self.out.branch("Tag_eta",  "F")
        self.out.branch("Tag_phi",  "F")
        self.out.branch("Tag_mass", "F")
        self.out.branch("Tag_iso",  "F")
        self.out.branch("Tag_dxy",  "F")
        self.out.branch("Tag_dz",   "F")
        self.out.branch("Tag_charge", "I")
        self.out.branch("Tag_isGenMatched", "I")
        self.out.branch("Probe_pt",   "F")
        self.out.branch("Probe_eta",  "F")
        self.out.branch("Probe_phi",  "F")
        self.out.branch("Probe_mass", "F")
        self.out.branch("Probe_charge", "I")
        self.out.branch("Probe_dxy",   "F")
        self.out.branch("Probe_dz",    "F")
        self.out.branch("Probe_SIP3D", "F")
        self.out.branch("Probe_iso",   "F")
        self.out.branch("Probe_miniiso",     "F")
        self.out.branch("Probe_passL",       "I")
        self.out.branch("Probe_passM",       "I")
        self.out.branch("Probe_passMP",      "I")
        self.out.branch("Probe_passT",       "I")
        self.out.branch("Probe_passRelIsoVL","I")
        self.out.branch("Probe_passRelIsoL", "I")
        self.out.branch("Probe_passRelIsoM", "I")
        self.out.branch("Probe_passRelIsoT", "I")
        self.out.branch("Probe_passMiniIsoL","I")
        self.out.branch("Probe_passMiniIsoM","I")
        self.out.branch("Probe_passMiniIsoT","I")
        self.out.branch("Probe_passMiniIsoVT", "I")
        self.out.branch("Probe_passMultiIsoL", "I")
        self.out.branch("Probe_passMultiIsoM", "I")
        self.out.branch("Probe_passMultiIsoM2017",   "I")
        self.out.branch("Probe_passMultiIsoM2017v2", "I")
        self.out.branch("Probe_passttH",   "I")
        self.out.branch("Probe_ptRatio",   "F")
        self.out.branch("Probe_ptRel",     "F")
        self.out.branch("Probe_jetRelIso", "F")
        self.out.branch("Probe_conept",    "F")
        self.out.branch("Probe_jetbtagdeepcsv", "F")
        self.out.branch("Probe_pdgId", "I")
        self.out.branch("Probe_mvaTTH",      "F")
        self.out.branch("Probe_mvaFall17V2noIso_WPL",  "F")
        self.out.branch("Probe_mvaFall17V2noIso",  "F")
        self.out.branch("Probe_lostHits",              "F")
        self.out.branch("Probe_tightCharge",           "I")
        self.out.branch("Probe_convVeto",              "I")
        self.out.branch("Probe_ttH_idEmu_cuts_E3",   "I")

        self.out.branch("Probe_passDptPt02",      "I")
        self.out.branch("Probe_passSIP4",         "I")
        self.out.branch("Probe_passSIP8",         "I")
        self.out.branch("Probe_passMVAL",         "I")
        self.out.branch("Probe_passMVAM",         "I")
        self.out.branch("Probe_passMVAT",         "I")
        self.out.branch("Probe_isGenMatched",     "I")
        self.out.branch("TnP_mass", "F")
        self.out.branch("TnP_ht",   "F")
        self.out.branch("TnP_met",  "F")
        self.out.branch("TnP_trigger", "I")
        self.out.branch("TnP_npairs",  "I")

    def jetLepAwareJEC(self,lep,jet,L1corr):
      p4l = lep.p4(); l = ROOT.TLorentzVector(p4l.Px(),p4l.Py(),p4l.Pz(),p4l.E())
      if not hasattr(jet,'rawFactor'): return l
      c = jet.rawFactor
      f = 1 - c # factor to go to rawpt
      p4j = jet.p4(); j = ROOT.TLorentzVector(p4j.Px(),p4j.Py(),p4j.Pz(),p4j.E())
      if ((j*c-l).Rho()<1e-4): return l
      print "origpt = %1.2f, mod pt = %1.2f"%(j.Pt(), (j*f).Pt())
      if L1corr == 0: L1corr = 0.1
      j = (j*f - l*(1.0/L1corr)) * (1/f) + l
      return j

    def ptRelv2(self,lep,jet,L1corr): # use only if jetAna.calculateSeparateCorrections==True
      m = self.jetLepAwareJEC(lep,jet,L1corr)
      p4l = lep.p4(); l = ROOT.TLorentzVector(p4l.Px(),p4l.Py(),p4l.Pz(),p4l.E())
      if ((m-l).Rho()<1e-4): return 0 # lep.jet==lep (no match) or jet containing only the lepton
      return l.Perp((m-l).Vect())
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def IsMatched(self, electronLorentzVector, trigObjCollection):
      dRmin = 0.1
      match = False
      for trigObj in trigObjCollection:
        dR = electronLorentzVector.DeltaR(trigObj)
        if dR < dRmin: match = True
      return match

    def OpenJECcalibrator(self, jetType = "AK4PF", doRes = True):
        # For jet re-calibrating
        fileNameJEC = self.filenameJEC
        jesInputFilePath = os.environ['CMSSW_BASE'] + "/src/PhysicsTools/NanoAODTools/data/jme/" # By default
        print 'Using the file: ', jesInputFilePath+fileNameJEC
        return JetReCalibrator(fileNameJEC, jetType , doRes, jesInputFilePath, upToLevel=1)

    def analyze(self, event):
        # Get Jet and Ele collections
        jet      = Collection(event, 'Jet')
        electron = Collection(event, 'Electron')
        trigObj  = Collection(event, 'TrigObj')

        #### Construct the trigger object collection containing electrons triggering a single iso electron trigger
        selTrigObj = []
        for tr in trigObj:
          if not abs(tr.id) == 11: continue
          #print 'Trig electron found! with filterBits: ', tr.filterBits
          if not (bool(tr.filterBits & 2) or bool(tr.filterBits & 8)): continue
          t = TLorentzVector()
          t.SetPtEtaPhiM(tr.pt, tr.eta, tr.phi, 0.000510)
          selTrigObj.append(t)
        #print 'len(selTrigObj) = ', len(selTrigObj)

        #### Selection of tag and probe electrons
        # Tag: pT, eta, tightId, iso
        # Probe: pT, eta
        tags = []; probes = []; pair = []
        index = 0
        for ele in electron:
          if not abs(ele.eta) < 2.5 or not ele.pt > self.kProbePt: 
            index += 1
            continue
          if ele.pt > self.kTagPt and ele.cutBased>= 4 and ele.jetRelIso < self.kTagIso:  
            tp = TLorentzVector()
            tp.SetPtEtaPhiM(ele.pt, ele.eta, ele.phi, ele.mass)
            if self.IsMatched(tp, selTrigObj): tags.append(index)
          else: probes.append(index)
          index+=1
        if len(tags) == 0:              return False
        if len(tags) + len(probes) < 2: return False
        #print 'len(tags) =', len(tags)
        #print 'len(probes) =', len(probes)

        ### Forming TnP pairs with leptons passing Tag selection
        ### If two tight electrons, you have two TnP pairs!
        ### If one tag has more than one probe candidate, the pairs are not selected
        for t1 in tags:
          nProbesFound = 0
          tp1 = TLorentzVector()
          tp1.SetPtEtaPhiM(electron[t1].pt, electron[t1].eta, electron[t1].phi, electron[t1].mass)
          for t2 in tags:
            if t2 == t1: continue
            tp2 = TLorentzVector()
            tp2.SetPtEtaPhiM(electron[t2].pt, electron[t2].eta, electron[t2].phi, electron[t2].mass)
            mass = (tp1+tp2).M()
            if mass > self.kMaxMass or mass < self.kMinMass: continue
            ### The tag must match a trigger object:
            if nProbesFound == 0: 
              nProbesFound += 1
              pair.append([t1, t2])
            else: # This tag has 2 or more probes!!
              pair = pair[:-1]

        ### Forming TnP pairs with a Tag lepton and Probe leptons
        for t in tags:
          nProbesFound = 0
          tp = TLorentzVector()
          tp.SetPtEtaPhiM(electron[t].pt, electron[t].eta, electron[t].phi, electron[t].mass)
          for p in probes:
            pp = TLorentzVector()
            pp.SetPtEtaPhiM(electron[p].pt, electron[p].eta, electron[p].phi, electron[p].mass)
            mass = (tp+pp).M()
            if mass > self.kMaxMass or mass < self.kMinMass: continue
            ### The tag must match a trigger object:
            if nProbesFound == 0: 
              nProbesFound += 1
              pair.append([t, p])
            else: # This tag has 2 or more probes!!
              pair = pair[:-1]

        # Check that we have at least one pair... calculate the mass of the pair
        if len(pair) == 0: return False # events with 1 or 2 pairs!

        rho = event.fixedGridRhoFastjetAll
        #print "[%i] rho = %1.2f" %(self.i, rho)
        self.i += 1

        # Set variables for tag, probe and event
        for thisPair in pair:
          ti, pi = thisPair
          ptag = TLorentzVector(); ppro = TLorentzVector()
          ptag.SetPtEtaPhiM(electron[ti].pt, electron[ti].eta, electron[ti].phi, electron[ti].mass)
          ppro.SetPtEtaPhiM(electron[pi].pt, electron[pi].eta, electron[pi].phi, electron[pi].mass)
          mass        = (ptag+ppro).M()

          # Trigger requirement... IsoEle
          if   self.year == 17:
            passTrigger = event.HLT_Ele35_WPTight_Gsf
          elif self.year == 16:
            passTrigger = event.XXX or event.XXX
          elif self.year == 18:
            passTrigger = event.XXX
        
          # Compute HT and MET
          ht = 0; met = event.METFixEE2017_pt if self.year == 17 else event.MET_pt
          for j in jet: 
            if self.filenameJECrecal != "":
              pass
            else:
              ht += j.pt if j.pt > 30 else 0

          isdata = 0 if hasattr(event, 'Elecrton_genPartFlav') else 1

          # Tag kinematics
          self.out.fillBranch("Tag_pt",    electron[ti].pt)
          self.out.fillBranch("Tag_eta",   electron[ti].eta)
          self.out.fillBranch("Tag_phi",   electron[ti].phi)
          self.out.fillBranch("Tag_mass",  electron[ti].mass)
          self.out.fillBranch("Tag_charge",electron[ti].charge)
          self.out.fillBranch("Tag_iso",   electron[ti].jetRelIso)
          self.out.fillBranch("Tag_dz",    electron[ti].dz)
          self.out.fillBranch("Tag_dxy",   electron[ti].dxy)

          tagMatch = 1 if isdata else (electron[ti].genPartFlav == 1 or electron[ti].genPartFlav == 15)
          self.out.fillBranch("Tag_isGenMatched", tagMatch)

          # Probe kinematics
          self.out.fillBranch("Probe_pt",                    electron[pi].pt)
          self.out.fillBranch("Probe_eta",                   electron[pi].eta)
          self.out.fillBranch("Probe_phi",                   electron[pi].phi)
          self.out.fillBranch("Probe_mass",                  electron[pi].mass)
          self.out.fillBranch("Probe_charge",                electron[pi].charge)
          self.out.fillBranch("Probe_dxy",                   electron[pi].dxy)
          self.out.fillBranch("Probe_dz",                    electron[pi].dz)
          self.out.fillBranch("Probe_SIP3D",                 electron[pi].sip3d)
          self.out.fillBranch("Probe_iso",                   electron[pi].jetRelIso)
          self.out.fillBranch("Probe_miniiso",               electron[pi].miniPFRelIso_all)
          self.out.fillBranch("Probe_mvaTTH",                electron[pi].mvaTTH)
          self.out.fillBranch("Probe_mvaFall17V2noIso_WPL",  electron[pi].mvaFall17V2noIso_WPL)
          self.out.fillBranch("Probe_lostHits",              electron[pi].lostHits)
          self.out.fillBranch("Probe_tightCharge",           electron[pi].tightCharge)

          # Probe ID and ISO flags (from nanoAOD tag IDs)
          self.out.fillBranch("Probe_passL",         1)
          self.out.fillBranch("Probe_passM",         electron[pi].cutBased==3)
          self.out.fillBranch("Probe_passMP",        electron[pi].cutBased==3)
          self.out.fillBranch("Probe_passT",         electron[pi].cutBased>= 4)
          self.out.fillBranch("Probe_passRelIsoVL",  electron[pi].jetRelIso <= 0.4)
          self.out.fillBranch("Probe_passRelIsoM",   electron[pi].jetRelIso <= 0.25)
          self.out.fillBranch("Probe_passRelIsoL",   electron[pi].jetRelIso <= 0.20)
          self.out.fillBranch("Probe_passRelIsoT",   electron[pi].jetRelIso <= 0.15)
          self.out.fillBranch("Probe_passMiniIsoL",  electron[pi].miniPFRelIso_all < 0.4)
          self.out.fillBranch("Probe_passMiniIsoT",  electron[pi].miniPFRelIso_all < 0.2)
          self.out.fillBranch("Probe_passMiniIsoVT", electron[pi].miniPFRelIso_all < 0.1)

          probeMatch = 1 if isdata else (electron[pi].genPartFlav == 1 or electron[pi].genPartFlav == 15)
          self.out.fillBranch("Probe_isGenMatched", probeMatch)

          # Other working points: SIP2D and dpt/pt
          self.out.fillBranch("Probe_passDptPt02",   electron[pi].pt/electron[pi].pt < 0.2) #pterr not defined for electrons
          self.out.fillBranch("Probe_passSIP4",      electron[pi].sip3d < 4)
          self.out.fillBranch("Probe_passSIP8",      electron[pi].sip3d < 8)

          # Probe Lepton MVA (from nanoAOD, for the moment)
          self.out.fillBranch("Probe_passMVAL",      1) #mvaID not defined for electrons
          self.out.fillBranch("Probe_passMVAM",      1)
          self.out.fillBranch("Probe_passMVAT",      1)

          # MultiIso... calculate jet-JecLepAware ptRel, ptRatio
          jetId = electron[pi].jetIdx
          MiniRelIso = electron[pi].miniPFRelIso_all
          jetRelIso  = electron[pi].jetRelIso # jetRelIso = 1/ptRatio - 1
          mvaTTH           = electron[pi].mvaTTH
          pdgId            = electron[pi].pdgId
          mvaFall17V2noIso = electron[pi].mvaFall17V2Iso
          lostHits         = electron[pi].lostHits
          tightCharge      = electron[pi].tightCharge
          convVeto         = electron[pi].convVeto
          if jetId == -1: # No jet matching....
            ptRel = 0
            ptRatio = 1/(jetRelIso + 1)
            jetbtagdeepcsv = 0
          else:
            # from Ele_jetRelIso in nanoAOD
            if doCalculateJecLepAwareFromNanoAOD:
              jetRelIso  = electron[pi].jetRelIso # jetRelIso = 1/ptRatio - 1
              ptRatio = 1/(jetRelIso + 1)
              jetbtagdeepcsv = jet[jetId].btagDeepB
              jetJECLepAwarePt = electron[pi].pt*(jetRelIso + 1)
              lp = TLorentzVector(); jt = TLorentzVector()
              lp.SetPtEtaPhiM(electron[pi].pt, electron[pi].eta, electron[pi].phi, electron[pi].mass)
              jt.SetPtEtaPhiM(jetJECLepAwarePt,  jet[jetId].eta,  jet[jetId].phi,  jet[jetId].mass)
              ptRel = lp.Perp((jt-lp).Vect())
            else:
              j = jet[jetId]
              corr = self.jetReCalibrator.getCorrection(j, rho)
              #print 'corr = %1.2f, jet Pt = %1.2f' %(corr, j.pt)
              ptRel = self.ptRelv2(electron[pi], j, corr)
              #print 'ptRel = ', ptRel
              ptRatio = electron[pi].pt/self.jetLepAwareJEC(electron[pi], j, corr).Pt()
              #print 'ptRatio = ', ptRatio


          # Definitions from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSLeptonSF
          passMultiIsoL =       MiniRelIso < 0.20 and (ptRatio > 0.69 or ptRel > 6.0) 
          passMultiIsoM =       MiniRelIso < 0.16 and (ptRatio > 0.76 or ptRel > 7.2) 
          passMultiIsoM2017   = MiniRelIso < 0.12 and (ptRatio > 0.80 or ptRel > 7.5) 
          passMultiIsoM2017v2 = MiniRelIso < 0.11 and (ptRatio > 0.74 or ptRel > 6.8) 

          #if (abs(pdgId)!=11): ttH_idEmu_cuts_E3 = True
          #if (lep.hoe>=(0.10-0.00*(abs(lep.deltaEtaSC+lep.eta)>1.479))): return False
          #if (lep.eInvMinusPInv<=-0.04): return False
          #if (lep.sieie>=(0.011+0.019*(abs(lep.deltaEtaSC+lep.eta)>1.479))): return False
          #return True

          if mvaTTH > 0.90: 
              conept =  electron[pi].pt 
          else:
              conept = 0.90 * electron[pi].pt * (1 + jetRelIso)
 
              
          ttH_idEmu_cuts_E3 = 1
          if (abs(electron[pi].pdgId)!=11)             : ttH_idEmu_cuts_E3 = 1
          if (electron[pi].eInvMinusPInv<=-0.04)       : ttH_idEmu_cuts_E3 = 0
          if (electron[pi].hoe>=(0.10-0.00*(abs(electron[pi].deltaEtaSC+electron[pi].eta)>1.479)))    : ttH_idEmu_cuts_E3 = 0
          if (electron[pi].sieie>=(0.011+0.019*(abs(electron[pi].deltaEtaSC+electron[pi].eta)>1.479))): ttH_idEmu_cuts_E3 = 0

          passttH = conept > 10 and jetbtagdeepcsv < 0.4941 and ttH_idEmu_cuts_E3 and convVeto and lostHits==0 and mvaTTH>0.90  
         
          #print 'passMultiIsoL = ', passMultiIsoL
          self.out.fillBranch("Probe_passMultiIsoL",       passMultiIsoL)
          self.out.fillBranch("Probe_passMultiIsoM",       passMultiIsoM)
          self.out.fillBranch("Probe_passMultiIsoM2017",   passMultiIsoM2017)
          self.out.fillBranch("Probe_passMultiIsoM2017v2", passMultiIsoM2017v2)
          self.out.fillBranch("Probe_passttH",             passttH)
          self.out.fillBranch("Probe_ptRatio",             ptRatio)
          self.out.fillBranch("Probe_ptRel",               ptRel)
          self.out.fillBranch("Probe_jetRelIso",           jetRelIso)
          self.out.fillBranch("Probe_conept",              conept)
          self.out.fillBranch("Probe_jetbtagdeepcsv",      jetbtagdeepcsv)
          self.out.fillBranch("Probe_mvaFall17V2noIso",    mvaFall17V2noIso)
          self.out.fillBranch("Probe_pdgId",               pdgId)
          self.out.fillBranch("Probe_convVeto",            convVeto)
          self.out.fillBranch("Probe_ttH_idEmu_cuts_E3",   ttH_idEmu_cuts_E3)

          # TnP variables
          self.out.fillBranch("TnP_mass",     mass);
          self.out.fillBranch("TnP_trigger",  passTrigger); 
          self.out.fillBranch("TnP_npairs",   len(pair)); 
          self.out.fillBranch("TnP_met",      met);
          self.out.fillBranch("TnP_ht",       ht);
          self.out.fill()
        return False

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
addTnPEle16 = lambda : addTnPvarEle(0,16)
#addTnPEle17 = lambda : addTnPvarEle(0,17,"Fall17_17Nov2017_V32_MC")
addTnPEle17 = lambda : addTnPvarEle(0,17)
addTnPEle18 = lambda : addTnPvarEle(0,18)

addTnPEle16data = lambda : addTnPvarEle(1,16)
#addTnPEle17data = lambda : addTnPvarEle(1,17,"Fall17_17Nov2017_V32_DATA")
addTnPEle17data = lambda : addTnPvarEle(1,17)
addTnPEle18data = lambda : addTnPvarEle(1,18)
addTnPEle   = lambda : addTnPvarEle(0,17)
#addTnPEleForMoriond18  = lambda : addTnPvarEle(0,18, "Autumn18_V3_MC")
#addTnPEleForMoriond18data  = lambda : addTnPvarEle(1,18, "Autumn18_V3_DATA")
addTnPEleForMoriond18  = lambda : addTnPvarEle(0,18)
addTnPEleForMoriond18data  = lambda : addTnPvarEle(1,18)
