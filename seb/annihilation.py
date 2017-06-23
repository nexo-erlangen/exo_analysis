#!/usr/bin/env python
import argparse
import os
import cPickle
import numpy as np
from collections import defaultdict
from itertools import compress
import plot_support as ps
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

import cluster as c

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

# Energy correction
energyCorr = 2614.5 / 2526.97

FV = (162., 5., 182.)

def main():
    fOutName, srcString = c.getPath('S5', False)
    dataChain = c.getDataChain(fOutName, srcString, False)

    getClusterAnnihilation(dataChain, FV, 'ss')

def getClusterTrack(tree, FV, art=None):
    ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
    ROOT.gROOT.cd()

    # Perform cuts
    if art:
        cut = ps.getCut(calibCut=True, energyCut=False, type=art, MC=False)
    else:
        cut = ps.getCutCombined(calibCut=True, energyCut=True, MC=False, eMin=750, eMax=3500)
    treeCut = tree # tree.CopyTree( cut )

    N = treeCut.GetEntries()
    N = 10000000
    NValid = 0

    thetaList = []
    energy0, energy90 = [], []
    for i in range( N ):
        tree.GetEntry( i )
        es = treeCut.EventSummary
        mul = es.multiplicity

        if mul < 3:
            continue

        x, y, z = np.array(es.cluster_x), np.array(es.cluster_y), np.array(es.cluster_z)

        if mul != len(x):
            continue

        pos = zip(x, y, z)
        # If one cluster position is missing, skip the event
        if (abs( x ) > 900).any() or (abs( y ) > 900).any() or (abs( z ) > 900).any():
            continue

        energy = np.array( es.cluster_energy )
        energyScint = np.array( es.e_scint )

        vecList, vecListNorm = [], []
        idxPairs = [(m, m+1) for m in range(mul-1)] 
        for idx in idxPairs:
            clusterFirst, clusterSecond = np.array(pos[idx[0]]), np.array(pos[idx[1]])
            vec = clusterSecond - clusterFirst

            vecList.append( vec )
            vecListNorm.append( vec / np.linalg.norm(vec) )

        firstVec = vecListNorm[0]
        areaList, zDiffList = [], []
        for m in range(1, len(vecListNorm)):
            areaList.append( np.linalg.norm(np.cross(firstVec, vecListNorm[m])) )
            zDiffList.append( abs(firstVec[2] - vecListNorm[m][2]) )
            # theta.append( vecAngle(firstVec, vecListNorm[m])*180./np.pi )

        if np.any([area > 0.01 for area in areaList]): # or np.any([zDiff < 1. for zDiff in zDiffList]):
            continue

        theta = [vecAngle(vec, [1., 0., 0.])*180./np.pi for vec in vecListNorm]
        # theta = [abs(t - 180) if t > 90 else t for t in theta]
        # theta = sum(theta) / len(theta)
        # thetaList.append( theta )
        thetaList += theta

        '''
        if theta < 5:
            energy0.append( energyScint )
        elif theta > 85:
            energy90.append( energyScint )
        '''

        print pos
        print energy, sum( energy )
        print vecListNorm, np.dot(vecListNorm[0], vecListNorm[1]), (np.linalg.norm(vecListNorm[0]), np.linalg.norm(vecListNorm[1]))
        print [np.linalg.norm(vec) for vec in vecList]
        # print thetaList
        print

        NValid += 1

    print thetaList
    plt.hist(np.nan_to_num(thetaList), 40)
    #plt.hist(energy0, label='energy0')
    #plt.hist(energy90, label='energy90', alpha=0.5)
    plt.show()
    print 'Valid events:', (float(NValid) / N)

def getClusterAnnihilation(tree, FV, art=None):
    ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
    ROOT.gROOT.cd()

    # Perform cuts
    if art:
        cut = ps.getCut(calibCut=True, energyCut=False, type=art, MC=False)
    else:
        cut = ps.getCutCombined(calibCut=True, energyCut=True, MC=False, eMin=750, eMax=3500)
    treeCut = tree # tree.CopyTree( cut )

    N = treeCut.GetEntries()
    NValid = 0
    for i in range( N ):
        tree.GetEntry( i )
        es = treeCut.EventSummary
        mul = es.multiplicity

        if mul != 3:
            continue

        x, y, z = np.array(es.cluster_x), np.array(es.cluster_y), np.array(es.cluster_z)
        pos = zip(x, y, z)
        # If one cluster position is missing, skip the event
        if (abs( x ) > 900).any() or (abs( y ) > 900).any() or (abs( z ) > 900).any():
            continue

        energy = np.array( es.cluster_energy ) * energyCorr
        if len(energy) != 3 or abs(sum(energy) - 2614.5) > 60:
            continue

        if abs(energy[0] - energy[1]) < 100 or abs(energy[2] - energy[1]) < 100 or abs(energy[0] - energy[2]) < 100:
            cnt = 0 
            idxList = []
            r = 0 

            # Get clusters with minimum energies
            for m, en in enumerate( energy ):
                if en > 511 - 100:
                    cnt += 1
                    idxList.append( m )
                else:
                    r = m
            idxList.append( r )

            if cnt < 2:
                continue

            # Decide which clusters are most similar
            fracList = [energy[0]/energy[1], energy[0]/energy[2], energy[1]/energy[2]]
            fracList = [abs(1 - 1./frac) if 1./frac < frac else 1-frac for frac in fracList]
            fracIdx = fracList.index(min(fracList))
            # print fracList, fracIdx

            if fracIdx == 0:
                idxList = [0, 1, 2]
            elif fracIdx == 1:
                idxList = [0, 2, 1]
            else:
                idxList = [1, 2, 0]

            # if abs(energy[idxList[0]] + energy[idxList[1]] - 2*511 - energy[idxList[2]]) > 100:
            #    continue

            # print fracList, fracIdx
            # print idxList
            # if np.linalg.norm(np.array(pos[idxList[0]]) - np.array(pos[idxList[1]])) > 50:
            #    continue

            print energy
            # print
            # print x, y, z
            # raw_input('')

            NValid += 1
            if not NValid % 100:
                print NValid
                print

    print NValid
    print float(NValid)/N

def vecAngle(v1, v2):
    return np.arccos(np.dot(v1, v2))

if __name__ == '__main__':
    main()

