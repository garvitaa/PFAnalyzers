import math
import copy
import array
import numpy as np

etabins = array.array("d",[0.0, 0.5, 1.3, 2.1, 2.5, 3.0, 5.2])
ptbins = np.arange(0.,21.,2)
#[10,24,32,43,56,74,97,133,174,245,300,362,430,507,592,686,846,1032,1248,1588,2000,2500,3000,4000,6000])
deltaRbins = 0.01*np.arange(1,31,1)

def deltaR2( e1, p1, e2=None, p2=None):
  """Take either 4 arguments (eta,phi, eta,phi) or two objects that have 'eta', 'phi' methods)"""
  if (e2 == None and p2 == None):
    return deltaR2(e1.eta(),e1.phi(), p1.eta(), p1.phi())
  de = e1 - e2
  dp = deltaPhi(p1, p2)
  return de*de + dp*dp

def deltaPhi( p1, p2):
  '''Computes delta phi, handling periodic limit conditions.'''
  res = p1 - p2
  while res > math.pi:
    res -= 2*math.pi
  while res < -math.pi:
    res += 2*math.pi
  return res

# import ROOT in batch mode
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

from ctypes import c_uint8

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.AutoLibraryLoader.enable()

# Create histograms, etc.
ROOT.gROOT.SetStyle('Plain') # white background
H_ParticleEta_gen = ROOT.TH1F ("ParticleEta_gen","ParticleEta_gen",int(len(etabins)-1),etabins)
H_ParticleEta_matched = ROOT.TH1F ("ParticleEta","ParticleEta",int(len(etabins)-1),etabins)
H_ParticleEta_fake = ROOT.TH1F ("ParticleEta_Fake","ParticleEta_Fake",int(len(etabins)-1),etabins)

H_ParticlePt_gen = ROOT.TH1F ("ParticlePt_gen","ParticlePt_gen",int(len(ptbins)-1),ptbins)
H_ParticlePt_matched = ROOT.TH1F ("ParticlePt_matched","ParticlePt_matched",int(len(ptbins)-1),ptbins)
H_ParticlePt_fake = ROOT.TH1F ("ParticlePt_fake","ParticlePt_fake",int(len(ptbins)-1),ptbins)

H_recoJet_eta = ROOT.TH1F ("recoJet_eta","recoJet_eta",int(len(etabins)-1),etabins)
#H_recoJet_eta_matched = ROOT.TH1F ("recoJet_eta_matched","recoJet_eta_matched",int(len(etabins)-1),etabins)
#H_recoJet_eta_fake = ROOT.TH1F ("recoJet_eta_fake","recoJet_eta_fake",int(len(etabins)-1),etabins)

H_recoJet_pt = ROOT.TH1F ("recojet_pt","recojet_pt",int(len(ptbins)-1),ptbins)
#H_recoJet_pt_matched = ROOT.TH1F ("recojet_pt_matched","recojet_pt_matched",int(len(ptbins)-1),ptbins)
#H_recoJet_pt_fake = ROOT.TH1F ("recojet_pt_fake","recojet_pt_fake",int(len(ptbins)-1),ptbins)

H_recoJetRaw_pt = ROOT.TH1F ("recojetRaw_pt","recojetRaw_pt",int(len(ptbins)-1),ptbins)
#H_recoJetRaw_pt_matched = ROOT.TH1F ("recojetRaw_pt_matched","recojetRaw_pt_matched",int(len(ptbins)-1),ptbins)
#H_recoJetRaw_pt_fake = ROOT.TH1F ("recojetRaw_pt_fake","recojetRaw_pt_fake",int(len(ptbins)-1),ptbins)

H_pfcand_eta = ROOT.TH1F ("pfcand_eta","pfcand_eta",int(len(etabins)-1),etabins)
H_pfcand_eta_matched = ROOT.TH1F ("pfcand_eta_matched","pfcand_eta_matched",int(len(etabins)-1),etabins)
H_pfcand_eta_fake = ROOT.TH1F ("pfcand_eta_fake","pfcand_eta_fake",int(len(etabins)-1),etabins)

H_pfcand_pt = ROOT.TH1F ("pfcand_pt","pfcand_pt",int(len(ptbins)-1),ptbins)
H_pfcand_pt_matched = ROOT.TH1F ("pfcand_pt_matched","pfcand_pt_matched",int(len(ptbins)-1),ptbins)
H_pfcand_pt_fake = ROOT.TH1F ("pfcand_pt_fake","pfcand_pt_fake",int(len(ptbins)-1),ptbins)

H_deltaR = ROOT.TH1F ("deltaR", "deltaR", 1000, 0., 1.)
H_trk_nPFcands = ROOT.TH1F ("trk_nPFcand", "trk_nPFcand", 1001, -0.5, 1000.5)
H_calo_nPFcands = ROOT.TH1F ("calo_nPFcand", "calo_nPFcand", 1001, -0.5, 1000.5)
H_deltaR_energy = ROOT.TProfile ("deltaR_energy", "deltaR_energy", 1000, 0., 1., 0., 2.)

# single Particle response
H_response_pt=[]
H_response_track_pt=[]
H_response_charge_pt = []
H_response_neutral_pt = []

# summed particle response
H_response_pt_dR=[]
H_response_track_pt_dR=[]
H_response_charge_pt_dR = []
H_response_neutral_pt_dR = []

H_ParticlePt_matched_eta=[]
H_ParticlePt_eta=[]

for i in etabins:
  sub=[]
  sub_t=[]
  sub_c=[]
  sub_n=[]
  hname = "particle_pT_eta" + (str(i).replace('.', 'p'))
  H_ParticlePt_eta.append(ROOT.TH1F(hname, hname, int(len(ptbins)-1),ptbins))
  hname = "particle_matched_pT_eta" + (str(i).replace('.', 'p'))
  H_ParticlePt_matched_eta.append(ROOT.TH1F(hname, hname, int(len(ptbins)-1),ptbins))
  for j in ptbins:
    hname = "response_pt" + (str(j).replace('.', 'p')) + "_eta" + (str(i).replace('.', 'p'))
    sub.append(ROOT.TH1F(hname, hname, 52, 0, 2))
    hname = "response_track_pt" + (str(j).replace('.', 'p')) + "_eta" + (str(i).replace('.', 'p'))
    sub_t.append(ROOT.TH1F(hname, hname, 52, 0, 2))
    hname = "response_charge_pt" + (str(j).replace('.', 'p')) + "_eta" + (str(i).replace('.', 'p'))
    sub_c.append(ROOT.TH1F(hname, hname, 52, 0, 2))
    hname = "response_neutral_pt" + (str(j).replace('.', 'p')) + "_eta" + (str(i).replace('.', 'p'))
    sub_n.append(ROOT.TH1F(hname, hname, 52, 0, 2))
  H_response_pt.append(sub)
  H_response_track_pt.append(sub_t)
  H_response_charge_pt.append(sub_c)
  H_response_neutral_pt.append(sub_n)

for k in deltaRbins:
  res=[]
  res_t=[]
  res_c=[]
  res_n=[]
  for i in etabins:
    sub=[]
    sub_t=[]
    sub_c=[]
    sub_n=[]
    for j in ptbins:
      hname = "response_pt" + (str(j).replace('.', 'p')) + "_eta" + (str(i).replace('.', 'p')) + "_deltaR" + (str(k).replace('.', 'p'))
      sub.append(ROOT.TH1F(hname, hname, 52, 0, 2))
      hname = "response_track_pt" + (str(j).replace('.', 'p')) + "_eta" + (str(i).replace('.', 'p')) + "_deltaR" + (str(k).replace('.', 'p'))
      sub_t.append(ROOT.TH1F(hname, hname, 52, 0, 2))
      hname = "response_charge_pt" + (str(j).replace('.', 'p')) + "_eta" + (str(i).replace('.', 'p')) + "_deltaR" + (str(k).replace('.', 'p'))
      sub_c.append(ROOT.TH1F(hname, hname, 52, 0, 2))
      hname = "response_neutral_pt" + (str(j).replace('.', 'p')) + "_eta" + (str(i).replace('.', 'p')) + "_deltaR" + (str(k).replace('.', 'p'))
      sub_n.append(ROOT.TH1F(hname, hname, 52, 0, 2))
    res.append(sub)
    res_t.append(sub_t)
    res_c.append(sub_c)
    res_n.append(sub_n)
  H_response_pt_dR.append(res)
  H_response_track_pt_dR.append(res_t)
  H_response_charge_pt_dR.append(res_c)
  H_response_neutral_pt_dR.append(res_n) 

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

pgenpars, pgenParLabel = Handle("std::vector<pat::PackedGenParticle>"), "packedGenParticles"
jets, jetLabel = Handle("std::vector<pat::Jet>"), "slimmedJets"
pjets, pjetLabel = Handle("std::vector<pat::Jet>"), "slimmedJetsPuppi"
pfcands, pfcandLabel = Handle("std::vector<pat::PackedCandidate>"), "packedPFCandidates"

verticesScore = Handle("edm::ValueMap<float>")

# open file (you can use 'edmFileUtil -d /store/whatever.root' to get the physical file name)
#events = Events('file:step4_inMINIAODSIM.root')

# test
# pfTICL
#listFiles=[
#   'file:8c20b216-dd17-43e2-878e-d6ceb8c735f5.root'
#]
# ref
listFiles=[
    'file:3161edd0-2375-4aea-ac8e-bc6c416ac9c4.root'
]
events = Events(listFiles)

for iev,event in enumerate(events):
  #if iev >= 1: break
  event.getByLabel(pgenParLabel, pgenpars)
  event.getByLabel(jetLabel, jets)
  event.getByLabel(pjetLabel, pjets)
  event.getByLabel(pfcandLabel,pfcands)
  print "\nEvent: run %6d, lumi %4d, event %12d" % (event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event())

  # PackedGenParticles
  for i,pgenp in enumerate(pgenpars.product()): 
    #if pgenp.pt() < 5 : continue
    print "pgenpar: run %6d, event %10d, pt %4.1f, eta %5.2f, phi %5.2f, pdgId %d." % (
      event.eventAuxiliary().run(), event.eventAuxiliary().event(), pgenp.pt(), pgenp.eta(), pgenp.phi(), pgenp.pdgId())
    i_eta = int(np.digitize([abs(pgenp.eta())],etabins))-1
    i_pt = int(np.digitize([pgenp.pt()],ptbins))-1
    H_ParticlePt_gen.Fill(pgenp.pt())
    H_ParticleEta_gen.Fill(abs(pgenp.eta()))
    if (i_eta <len(etabins)): H_ParticlePt_eta[i_eta].Fill(pgenp.pt())

    dr_min = 99999.0
    match = False # matching gen particles to pfcands
    for k,j in enumerate(pfcands.product()):
      if j.pt() < 0:continue
      H_deltaR.Fill(math.sqrt(deltaR2(pgenp,j)))
      H_deltaR_energy.Fill(math.sqrt(deltaR2(pgenp,j)), (j.energy()/pgenp.energy()), 1)
      if (abs(pgenp.pdgId())==211):
        if (deltaR2(pgenp,j) < dr_min):
          dr_min = deltaR2(pgenp,j)
          if (dr_min<0.0001):
            match=True
            recop = j

    if match:
      print "matched: run %6d, event %10d, pt %5.1f eta %5.2f phi %5.2f pdgId %5d %5.3f %5.3f " % ( event.eventAuxiliary().run(),event.eventAuxiliary().event(),recop.pt(), recop.eta(), recop.phi(), recop.pdgId(), recop.rawCaloFraction(), recop.rawHcalFraction())
      H_ParticlePt_matched.Fill(pgenp.pt())
      H_ParticleEta_matched.Fill(pgenp.eta())
      print "    RecopT: %5.3f GenpT: %5.3f Response: %5.3f i_pt: %3i i_eta: %3i" % (recop.pt(), pgenp.pt(), (recop.pt() / pgenp.pt()), i_pt, i_eta)
      if (i_eta <len(etabins)) & (i_pt < len(ptbins)):
        H_ParticlePt_matched_eta[i_eta].Fill(pgenp.pt())
        H_response_pt[i_eta][i_pt].Fill(recop.pt() / pgenp.pt())
        H_response_track_pt[i_eta][i_pt].Fill(recop.ptTrk() / pgenp.pt())
        if (abs(recop.charge()) > 0):
          H_response_charge_pt[i_eta][i_pt].Fill(recop.pt() / pgenp.pt())
        else:
          H_response_neutral_pt[i_eta][i_pt].Fill(recop.pt() / pgenp.pt())
        for d in range(len(deltaRbins)):
          pt_sum=0
          pt_track_sum=0
          for k,j in enumerate(pfcands.product()):
            if j.pt() < 0: continue
            if (deltaR2(pgenp,j) < (deltaRbins[d])**2):
              pt_sum+=j.pt()
              pt_track_sum+=j.ptTrk()
          print pt_sum
          H_response_pt_dR[d][i_eta][i_pt].Fill(pt_sum / pgenp.pt())
          H_response_track_pt_dR[d][i_eta][i_pt].Fill(pt_track_sum / pgenp.pt())
          if (abs(recop.charge()) > 0):
            H_response_charge_pt_dR[d][i_eta][i_pt].Fill(pt_sum / pgenp.pt())
          else:
            H_response_neutral_pt_dR[d][i_eta][i_pt].Fill(pt_sum / pgenp.pt())

    else:
      H_ParticlePt_fake.Fill(pgenp.pt())
      H_ParticleEta_fake.Fill(pgenp.eta())

    #count number of pf candidates in the vicinity of matched candidate
    trk_count = 0
    calo_count = 0
    for k,j in enumerate(pfcands.product()):
      if j.pt() < 0:continue
      if (deltaR2(recop,j) < 0.01):
        if ((recop.pt() / pgenp.pt()) > 0.8):
          trk_count+=1
        else:
          calo_count+=1
    H_trk_nPFcands.Fill(trk_count)
    H_calo_nPFcands.Fill(calo_count)

  # Jets (standard AK4)
  for i,j in enumerate(jets.product()):
    print "jet: run %6d, event %10d, pt %5.1f (raw pt %5.1f, matched-calojet pt %5.1f), eta %+4.2f, btag run1(CSV) ) %.3f, run2(pfCSVIVFV2) %.3f, pileup mva disc %+.2f" % (
      event.eventAuxiliary().run(), event.eventAuxiliary().event(), j.pt(), j.pt()*j.jecFactor('Uncorrected'), j.userFloat("caloJetMap:pt"), j.eta(), max(0,j.bDiscriminator("combinedSecondaryVertexBJetTags")), max(0,j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")), j.userFloat("pileupJetId:fullDiscriminant"))

    H_recoJet_pt.Fill(j.pt())
    H_recoJetRaw_pt.Fill(j.pt()*j.jecFactor('Uncorrected'))
    H_recoJet_eta.Fill(j.eta())

    # match = False
    # for k,pgenp in enumerate(pgenpars.product()):
    #   if abs(pgenp.pdgId())==211:
    #     if deltaR2(pgenp,j)<0.01:
    #       match=True

    # if match:
    #   H_recoJet_pt_matched.Fill(j.pt())
    #   H_recoJetRaw_pt_matched.Fill(j.pt()*j.jecFactor('Uncorrected'))
    #   H_recoJet_eta_matched.Fill(j.eta())
    # else:
    #   H_recoJet_pt_fake.Fill(j.pt())
    #   H_recoJetRaw_pt_fake.Fill(j.pt()*j.jecFactor('Uncorrected'))
    #   H_recoJet_eta_fake.Fill(j.eta())

  # pfcands
  for i,j in enumerate(pfcands.product()):
    if j.pt() < 0: continue
    print "pfcands: run %6d, event %10d, pt %5.1f eta %5.2f phi %5.2f pdgId %5d %5.3f %5.3f " % ( event.eventAuxiliary().run(), event.eventAuxiliary().event(), j.pt(), j.eta(), j.phi(), j.pdgId(), j.rawCaloFraction(), j.rawHcalFraction())
    
    H_pfcand_pt.Fill(j.pt())
    H_pfcand_eta.Fill(j.eta())

    match = False
    for k,pgenp in enumerate(pgenpars.product()):
      if abs(pgenp.pdgId())==211:
        if deltaR2(pgenp,j)<0.0001:
          match=True

    if match:
      H_pfcand_pt_matched.Fill(j.pt())
      H_pfcand_eta_matched.Fill(j.eta())
    else:
      H_pfcand_pt_fake.Fill(j.pt())
      H_pfcand_eta_fake.Fill(j.eta())
    

f = ROOT.TFile.Open("myfile_ref.root","RECREATE")
H_ParticleEta_gen.Write()
H_ParticleEta_matched.Write()

H_ParticlePt_gen.Write()
H_ParticlePt_matched.Write()

H_recoJet_eta.Write()
H_recoJet_pt.Write()
H_recoJetRaw_pt.Write()

H_pfcand_eta.Write()
H_pfcand_eta_matched.Write()
H_pfcand_eta_fake.Write()

H_pfcand_pt.Write()
H_pfcand_pt_matched.Write()
H_pfcand_pt_fake.Write()
H_deltaR.Write()
H_trk_nPFcands.Write()
H_calo_nPFcands.Write()
H_deltaR_energy.Write()

for i in range(len(etabins)): 
  H_ParticlePt_matched_eta[i].Write()
  H_ParticlePt_eta[i].Write()
  for j in range(len(ptbins)):
    H_response_pt[i][j].Write()
    H_response_track_pt[i][j].Write()
    H_response_charge_pt[i][j].Write()
    H_response_neutral_pt[i][j].Write()
    for k in range(len(deltaRbins)):
      H_response_pt_dR[k][i][j].Write()
      H_response_track_pt_dR[k][i][j].Write()
      H_response_charge_pt_dR[k][i][j].Write()
      H_response_neutral_pt_dR[k][i][j].Write()
    #H_responseRaw_pt[i][j].Write()
#ROOT.TFile.Close(f)
f.Write()
f.Close()
