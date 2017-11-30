#!/usr/bin/env python

import numpy as np
import sys
import cPickle
import os

import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')
ROOT.PyConfig.IgnoreCommandLineOptions = True

import plot_support as ps

def main():
    fName = '/home/vault/capm/mppi025h/analysis/preprocess/S5Th/Xe134_run_4628_tree.root'
    FV = [162, 5, 182]

    ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
    ROOT.gROOT.cd()

    cutData = ps.getCut(calibCut=True, energyCut=True, type='ms', MC=False, eMin=750, eMax=3500)
    f = ROOT.TFile.Open(fName, 'READ')
    t = f.Get('dataTree')
    treeCut = tree.CopyTree(cutData, 'fast')

    N = treeCut.GetEntries()
    for i in range( N ):
        treeCut.GetEntry(i)
        es = treeCut.EventSummary
        mul = es.multiplicity

        x, y, z = ps.getClusterPos(es, True)
        if np.isnan(x) or np.isnan(y) or np.isnan(z):
            continue

if __name__ == '__name__':
    main()

