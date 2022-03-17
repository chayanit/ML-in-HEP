#!/usr/bin/env python

import sys

import ROOT
from ROOT import TFile, TTree, TChain, TVector3, gROOT, addressof

# A C/C++ structure is required, to allow memory based access
gROOT.ProcessLine(
  "struct mssmhbb_str {\
    Int_t    njets;\
    Int_t    nbjets;\
    Float_t  mj1;\
    Float_t  ptj1;\
    Float_t  etaj1;\
    Float_t  phij1;\
    Int_t    flavorj1;\
    Float_t  mj2;\
    Float_t  ptj2;\
    Float_t  etaj2;\
    Float_t  phij2;\
    Int_t    flavorj2;\
    Float_t  mj3;\
    Float_t  ptj3;\
    Float_t  etaj3;\
    Float_t  phij3;\
    Int_t    flavorj3;\
    Float_t  pt12;\
    Float_t  eta12;\
    Float_t  phi12;\
    Float_t  m12;\
    Float_t  dR12;\
    Float_t  dEta12;\
    Float_t  dPhi12;\
    Float_t  pt13;\
    Float_t  eta13;\
    Float_t  phi13;\
    Float_t  m13;\
    Float_t  dR13;\
    Float_t  dEta13;\
    Float_t  dPhi13;\
    Float_t  pt23;\
    Float_t  eta23;\
    Float_t  phi23;\
    Float_t  m23;\
    Float_t  dR23;\
    Float_t  dEta23;\
    Float_t  dPhi23;\
  };");

try:
  input = raw_input
except:
  pass

if len(sys.argv) != 3:
  print(" Usage: %s input_file output_file "%(sys.argv[0]))
  sys.exit(1)

ROOT.gSystem.Load("libDelphes")

try:
  ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
  ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
except:
  pass

inputFile  = sys.argv[1]
outputFile = sys.argv[2]
print ("Reading from ", inputFile, "and writing to ", outputFile)

# Create chain of root trees
chain = ROOT.TChain("Delphes")
chain.Add(inputFile)

# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()

# Get pointers to branches used in this analysis
branchJet = treeReader.UseBranch("Jet")

# Output file
outfile = ROOT.TFile.Open(outputFile,"RECREATE")
outfile.cd()

# Book new Tree
mssmhbb = ROOT.TTree("mssmhbb", "mssmhbb")
tree_out = ROOT.mssmhbb_str()
# Create Tree branch
mssmhbb.Branch("njets",addressof( tree_out, "njets" ),"njets/I")
mssmhbb.Branch("nbjets",addressof( tree_out, "nbjets"),"nbjets/I")
# kinematic Jet1
mssmhbb.Branch("mj1",addressof( tree_out, "mj1"),"mj1/F")
mssmhbb.Branch("ptj1",addressof( tree_out, "ptj1"),"ptj1/F")
mssmhbb.Branch("etaj1",addressof( tree_out, "etaj1"),"etaj1/F")
mssmhbb.Branch("phij1",addressof( tree_out, "phij1"),"phij1/F")
mssmhbb.Branch("flavorj1",addressof( tree_out, "flavorj1"),"flavorj1/I")
# kinematic Jet2
mssmhbb.Branch("mj2",addressof( tree_out, "mj2"),"mj2/F")
mssmhbb.Branch("ptj2",addressof( tree_out, "ptj2"),"ptj2/F")
mssmhbb.Branch("etaj2",addressof( tree_out, "etaj2"),"etaj2/F")
mssmhbb.Branch("phij2",addressof( tree_out, "phij2"),"phij2/F")
mssmhbb.Branch("flavorj2",addressof( tree_out, "flavorj2"),"flavorj2/I")
# kinematic Jet3
mssmhbb.Branch("mj3",addressof( tree_out, "mj3"),"mj3/F")
mssmhbb.Branch("ptj3",addressof( tree_out, "ptj3"),"ptj3/F")
mssmhbb.Branch("etaj3",addressof( tree_out, "etaj3"),"etaj3/F")
mssmhbb.Branch("phij3",addressof( tree_out, "phij3"),"phij3/F")
mssmhbb.Branch("flavorj3",addressof( tree_out, "flavorj3"),"flavorj3/I")
# 4 vector Jet(1,2)
mssmhbb.Branch("pt12",addressof( tree_out, "pt12"),"pt12/F")
mssmhbb.Branch("eta12",addressof( tree_out, "eta12"),"eta12/F")
mssmhbb.Branch("phi12",addressof( tree_out, "phi12"),"phi12/F")
mssmhbb.Branch("m12",addressof( tree_out, "m12"),"m12/F")
mssmhbb.Branch("dR12",addressof( tree_out, "dR12"),"dR12/F")
mssmhbb.Branch("dEta12",addressof( tree_out, "dEta12"),"dEta12/F")
mssmhbb.Branch("dPhi12",addressof( tree_out, "dPhi12"),"dPhi12/F")
# 4 vector Jet(1,3)
mssmhbb.Branch("pt13",addressof( tree_out, "pt13"),"pt13/F")
mssmhbb.Branch("eta13",addressof( tree_out, "eta13"),"eta13/F")
mssmhbb.Branch("phi13",addressof( tree_out, "phi13"),"phi13/F")
mssmhbb.Branch("m13",addressof( tree_out, "m13"),"m13/F")
mssmhbb.Branch("dR13",addressof( tree_out, "dR13"),"dR13/F")
mssmhbb.Branch("dEta13",addressof( tree_out, "dEta13"),"dEta13/F")
mssmhbb.Branch("dPhi13",addressof( tree_out, "dPhi13"),"dPhi13/F")
# 4 vector Jet(2,3)
mssmhbb.Branch("pt23",addressof( tree_out, "pt23"),"pt23/F")
mssmhbb.Branch("eta23",addressof( tree_out, "eta23"),"eta23/F")
mssmhbb.Branch("phi23",addressof( tree_out, "phi23"),"phi23/F")
mssmhbb.Branch("m23",addressof( tree_out, "m23"),"m23/F")
mssmhbb.Branch("dR23",addressof( tree_out, "dR23"),"dR23/F")
mssmhbb.Branch("dEta23",addressof( tree_out, "dEta23"),"dEta23/F")
mssmhbb.Branch("dPhi23",addressof( tree_out, "dPhi23"),"dPhi23/F")

# Loop over all events
for entry in range(0, numberOfEntries):
  # Load selected branches with data from specified event
  treeReader.ReadEntry(entry)

  tree_out.njets  = int(branchJet.GetEntries())
  tree_out.nbjets = 0

  # Loop over Jet collection to count number of b-jets
  for i in range(0, branchJet.GetEntries()):
    jet = branchJet.At(i)
    if jet.Flavor == 5:
      tree_out.nbjets = tree_out.nbjets + 1 

  # Pre-selection (event contains at least 3 jet)
  if branchJet.GetEntries() > 3:

    # Take first 3 jets
    jet1 = branchJet.At(0)
    jet2 = branchJet.At(1)
    jet3 = branchJet.At(2)

    # momentum vector of first 3 jets
    #pjet1 = ROOT.TVector3()
    #pjet2 = ROOT.TVector3()
    #pjet3 = ROOT.TVector3()
 
    #pjet1.SetPtEtaPhi(jet1.PT,jet1.Eta,jet1.Phi)
    #pjet2.SetPtEtaPhi(jet2.PT,jet2.Eta,jet2.Phi)
    #pjet3.SetPtEtaPhi(jet3.PT,jet3.Eta,jet3.Phi)

    # 4-momentum vector of first 3 jets
    pjet1 = ROOT.TLorentzVector()
    pjet2 = ROOT.TLorentzVector()
    pjet3 = ROOT.TLorentzVector()

    pjet1.SetPtEtaPhiM(jet1.PT,jet1.Eta,jet1.Phi,jet1.Mass)
    pjet2.SetPtEtaPhiM(jet2.PT,jet2.Eta,jet2.Phi,jet2.Mass)
    pjet3.SetPtEtaPhiM(jet3.PT,jet3.Eta,jet3.Phi,jet3.Mass)
 
    # dijet system
    dijet12 = pjet1 + pjet2
    dijet13 = pjet1 + pjet3
    dijet23 = pjet2 + pjet3
    
    # Print jet transverse momentum
    #print(jet1.PT, jet2.PT, jet3.PT)

    tree_out.mj1      = float(jet1.Mass)
    tree_out.ptj1     = float(jet1.PT)
    tree_out.etaj1    = float(jet1.Eta)
    tree_out.phij1    = float(jet1.Phi)
    tree_out.flavorj1 = int(jet1.Flavor)

    tree_out.mj2      = float(jet2.Mass)
    tree_out.ptj2     = float(jet2.PT)
    tree_out.etaj2    = float(jet2.Eta)
    tree_out.phij2    = float(jet2.Phi)
    tree_out.flavorj2 = int(jet2.Flavor)

    tree_out.mj3      = float(jet3.Mass)
    tree_out.ptj3     = float(jet3.PT)
    tree_out.etaj3    = float(jet3.Eta)
    tree_out.phij3    = float(jet3.Phi)
    tree_out.flavorj3 = int(jet3.Flavor)

    tree_out.pt12     = float(dijet12.Pt())
    tree_out.eta12    = float(dijet12.Eta())
    tree_out.phi12    = float(dijet12.Phi())
    #tree_out.m12      = float(jet1.Mass + jet2.Mass)
    tree_out.m12      = float(abs(dijet12.M()))
    tree_out.dR12     = float(pjet1.DeltaR(pjet2))
    tree_out.dEta12   = float(jet1.Eta - jet2.Eta)
    tree_out.dPhi12   = float(pjet1.DeltaPhi(pjet2))

    tree_out.pt13     = float(dijet13.Pt())
    tree_out.eta13    = float(dijet13.Eta())
    tree_out.phi13    = float(dijet13.Phi())
    #tree_out.m13      = float(jet1.Mass + jet3.Mass)
    tree_out.m13      = float(abs(dijet13.M()))
    tree_out.dR13     = float(pjet1.DeltaR(pjet3))
    tree_out.dEta13   = float(jet1.Eta - jet3.Eta)
    tree_out.dPhi13   = float(pjet1.DeltaPhi(pjet3))

    tree_out.pt23     = float(dijet23.Pt())
    tree_out.eta23    = float(dijet23.Eta())
    tree_out.phi23    = float(dijet23.Phi())
    #tree_out.m23      = float(jet2.Mass + jet3.Mass)
    tree_out.m23      = float(abs(dijet23.M()))
    tree_out.dR23     = float(pjet2.DeltaR(pjet3))
    tree_out.dEta23   = float(jet2.Eta - jet3.Eta)
    tree_out.dPhi23   = float(pjet2.DeltaPhi(pjet3))

  # Fill TTree
  mssmhbb.Fill()

#Write and Close
mssmhbb.Write()
outfile.Close()  

