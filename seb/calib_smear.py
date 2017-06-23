#!/usr/bin/env python
import numpy as np
import ROOT
import plot_support as ps

# /home/vault/capm/mppi025h/EXOoffline/offline/utilities/calib/src/EXOEnergyResol.cc
ROOT.gSystem.Load('libEXOCalibUtilities')

preProcessDir = '/home/vault/capm/mppi025h/analysis/preprocess/'
# mcFile = ROOT.TFile(preProcessDir + 'SourceS5_Th228.root', 'READ')
# mcFile = ROOT.TFile(preProcessDir + 'pre_artDrift_c_m-1.593_c_m2-0_001_t_h1_177_c_A7_6714.root', 'READ')
# mcFile = ROOT.TFile(preProcessDir + 'pre_10m_Th_real_nodiff_pre.root', 'READ')
mcFile = ROOT.TFile(preProcessDir + 'pre_10m_Th_real_nodiff_pre.root', 'READ')

mcTree = mcFile.Get('mcTree')
cutData = ps.getCut(calibCut=True, energyCut=True, type='ms', MC=True, eMin=700, eMax=3500)
ROOT.gROOT.cd()
mcTreeCut = mcTree.CopyTree( cutData )

mcTreeCut.SetEstimate(mcTreeCut.GetEntries()+1)
mcTreeCut.Draw('energy_mc:standoff_distance','','para goff')

resol = ROOT.EXOEnergyResol.GetInstanceForFlavor('2016-v3-weekly', 'energy-mcbased-fit')

# TH1D SmearedMC(const std::string& channel,double* energies, double* weights, int length, int multiplicity, int nBins, double xLow, double xUp, long int seconds, int nano, int binMC = 1) const;
# print list( np.ndarray(mcTreeCut.GetEntries(), dtype=np.float32, buffer=mcTreeCut.GetV2() ) )
# w = mcTreeCut.GetV2()
w = np.array( [1.] * mcTreeCut.GetEntries() )
# print w
# print len(w)

h = resol.SmearedMC('Rotated', mcTreeCut.GetV1(), w, mcTreeCut.GetSelectedRows(), 1, 300, 0, 3000, 1608494400, 0, 1)
c = ROOT.TCanvas() 
c.SetLogy()
h.Draw('')
raw_input('')

'''
outFile = ROOT.TFile('hey.root','recreate')
h.Write()
outFile.Close()
'''
del h
mcFile.Close()

