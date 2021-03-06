#!/usr/bin/env python
import argparse
import os
import cPickle
import numpy as np
from collections import defaultdict
from itertools import compress
import plot_support as ps
from matplotlib.backends.backend_pdf import PdfPages

import energy_smear as ensm

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

# Phase 2 data?
PHASE2 = False

# Fiducial Volume
FV = (171, 10, 182)
# FV = (171, 0, 182)

# Multiplicity
MULTIPLICITY = 2

# Source position (S5)
if PHASE2:
    SRC_S5 = (258.22, 7.9511, 10.292)
else:
    SRC_S5 = (255., 3.9, -30.)
# Source position (S2)
SRC_S2 = (-25.4, 2.3, -292.)
# Source position (S8)
SRC_S8 = (35, 0., 294.8)

# Rotated energy resolution
if PHASE2:
    resolution = [6.44094447e-1, 1.66891442e+1, -3.09627265e-9]
else:
    resolution = [0.714076, 21.9255, 0.]

# Cluster energy resolution
# resolutionCluster = [1.81960744, 0.00306382, 0.03331489]
if PHASE2:
    resolutionCluster = [3.68667951e-5, 6.73349953e+1, 2.29624217e-2]
else:
    resolutionCluster = [-.000123185451, 74.7268349, 0.0270834510]

# MC scintillation energy resolution
energyScintRes = 130 # keV

# = BIN RANGES =
# Choose binning
zBins, thetaBins = 40, 1
rInit, drInit = 183.2356, 183.2356

# Energy correction
if PHASE2:
    energyCorr = 2614.5 / 2651.59 # Shift = -5.16
    energyCorrRot = 2614.5 / 2608.86 # Shift = -0.59
else:
    energyCorr = 1.034517 # 2614.5 / 2534.17 # 2614.5 / 2526.97
    energyCorrRot = 1.

energyCorrScint = 0.984497 # 0.965230
energyCorrMC = 2614.5 / 2608.31 # 2614.5 / 2601.87

# Get bin ranges
# Recursion formula for the bins in r
def dr(r, ri, dri):
    return -r + np.sqrt(r**2 + dri*(2*ri + dri))

def rGetBins(ri, dri):
    rList = [ri]
    r = ri
    while True:
        dr_ = dr(r, ri, dri)
        r -= dr_
        if r < 0:
            break

        rList.append( r )

    # Check if last entry is 0,
    # if not, append 0
    if rList[-1]:
        rList.append( 0 )
    return np.array( list(reversed(rList)) )

zRange = np.append(np.linspace(-FV[2], -FV[1], (zBins+1) // 2), np.linspace(FV[1], FV[2], (zBins+1) // 2))
thetaRange = np.linspace(0, 2*np.pi, thetaBins + 1)
rRange = rGetBins(rInit, drInit)
print rRange, thetaRange, zRange

def vectorDist():
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    vecGamma = 2614.5 * np.array( [-1, 0] )
    theta = np.linspace(0, np.pi, 20)
    for t in theta:
        Ee = comptonEnergy(2614.5, t)
        vecGammaSec = np.array([np.cos(np.pi - t), np.sin(np.pi - t)])
        vecGammaSec =  (2614.5 - Ee) * vecGammaSec/np.linalg.norm(vecGammaSec)

        vecElectron = vecGamma - vecGammaSec

        print vecGamma / np.linalg.norm(vecGamma)
        print vecGammaSec / np.linalg.norm(vecGammaSec)
        print vecElectron / np.linalg.norm(vecElectron)
        print vecAngle(vecElectron, np.array([0, -1]))*180/np.pi
        print np.sqrt(sum(np.square(vecElectron)))
        print

def testFunction():
    '''
    x = np.array( [np.array(range(100)), np.array(range(100))] )
    y = np.array( [500 + np.array(zip(314.2*np.sin(np.arange(100)), [0.001]*100)), 500 + np.array( zip(np.cos(np.arange(100)), [0.001]*100))] )
    print x[:,0]
    xlabel, ylabel = 'x', 'y'
    pp = PdfPages( 'cluster_test.pdf' )
    plotXY(x, y, xlabel, ylabel, pp=pp, smooth=False, title=None, show=True, light=True, polar=True)
    # plotHistList([np.sin(np.pi/100 * np.arange(100)), np.cos(np.pi/100 * np.arange(100))], title='', show=True, pp=pp, fit=False, light=False)
    pp.close()
    raw_input('')
    return
    '''

    from numpy.random import normal
    d = cPickle.load(open('cluster_test.p', 'rb'))
    d2 = cPickle.load(open('cluster_test2.p', 'rb'))
    x, y = np.array(d['x']), np.array(d['y'])
    x = (x - np.mean(x)) / abs(np.std(x))

    x2, y2 = np.array(d2['x']), np.array(d2['y'])
    x2 = (x2 - np.mean(x2)) / abs(np.std(x2))

    # x = normal(0, 20, 10000)
    # y = normal(0, 10, 10000)

    # pp = PdfPages( 'cluster_test.pdf' )
    pp = None
    kdeOverlay(x, y, x2, y2, '$\Delta E$ [$\sigma$]', r'Electron angle $\theta_\mathrm{e}$')
    # plotThree(x, y, xlabel='x', ylabel='y', title='default', pp=pp, iterpolate=True, show=False, project=True)
    # pp.close()

    return

def testPlot():
    import matplotlib.lines as mlines
    from scipy.stats import kde
    import seaborn as sns
    import matplotlib.pyplot as plt

    light = True

    # Change figure size
    # fig, ax = plt.subplots(figsize=(4, 8))
    fig, ax = plt.subplots(figsize=(.3 * 8.27, 11.69))
    # plt.subplots_adjust(top=.4)

    # Load data
    if light:
        d = np.array( cPickle.load(open('histPlotLight.p', 'rb')) )
    else:
        d = np.array( cPickle.load(open('histPlot.p', 'rb')) )

    # Electron angle
    elAngle = d[:,2]
    # Data and MC / Charge and Light
    data, mc = d[:,0], d[:,1]

    offset = 0
    offs = [offset]

    colors = sns.color_palette('deep').as_hex()

    peakList = []
    peakOffs = []

    # Loop over angles
    for i in range(0, len(elAngle), 1):
        d, m = np.array(data[i]), np.array(mc[i])

        bins = np.linspace(min([min(b) for b in [d, m]]), max([max(b) for b in [d, m]]), 200)
        fdata = kde.gaussian_kde(d)
        fmc = kde.gaussian_kde(m)
        dKde, mKde = fdata(bins), fmc(bins)

        peakList.append( [max(zip(dKde, bins))[-1], max(zip(mKde, bins))[-1]] )
        peakOffs.append([offset + max(dKde), offset + max(mKde)])

        if light:
            ax.fill_between(bins, mKde + offset, 0 * mKde + offset, zorder=-i, facecolor=colors[1], edgecolor='w', lw=1, alpha=1)
            ax.fill_between(bins, dKde + offset, 0 * dKde + offset, zorder=-i, facecolor=colors[0], edgecolor='w', lw=1, alpha=1)
            chargeLabel = mlines.Line2D([], [], color=colors[1], label='Charge')
            scintLabel = mlines.Line2D([], [], color=colors[0], label='Scintillation')
            plt.legend(handles=[chargeLabel, scintLabel], frameon=True, loc='upper left')

        else:
            ax.fill_between(bins, dKde + offset, 0 * dKde + offset, zorder=-i, facecolor=colors[1], edgecolor='w', lw=1, alpha=1)
            ax.fill_between(bins, mKde + offset, 0 * mKde + offset, zorder=-i, facecolor=colors[0], edgecolor='w', lw=1, alpha=1)
            dataLabel = mlines.Line2D([], [], color=colors[1], label='Data')
            mcLabel = mlines.Line2D([], [], color=colors[0], label='MC')
            plt.legend(handles=[dataLabel, mcLabel], frameon=True, loc='upper left')

        ax.axhline(y=offset, zorder=-i, linestyle='-', color='k', lw=2)

        if light:
            ax.set_xlim(2614.5 - 300, 2614.5 + 300)
        else:
            ax.set_xlim(1500, 2381.75)
            # ax.set_xlim(0, 2614.5)

        offset += .5 * max([dKde.max(), mKde.max()])
        offs.append(offset)

    print peakList, peakOffs
    peakColors = sns.color_palette('dark').as_hex()
    peakList, peakOffs = np.array(peakList), np.array(peakOffs)
    cMC, cData = '#dce1f0', '#d5f5d9'
    if light:
        cMC, cData = cData, cMC

    # Data
    plt.plot(peakList[:,0], peakOffs[:,0], c=cData, ls='--')
    # MC
    plt.plot(peakList[:,1], peakOffs[:,1], c=cMC, ls='--')

    sns.despine(left=True, bottom=False, offset=-22.5)
    # ax.set_yticks(offs, [r'$%.2f^\circ$' % elAngle[i] for i in range(0, elAngle, 2)])
    plt.yticks(offs[:-1], [r'$%.1f$' % elAngle[i] if not i % 2 else ''  for i in range(0, len(elAngle), 1)])
    plt.xlabel('Energy [keV]')
    plt.ylabel(r'Electron angle $\theta_\mathrm{e}$ [$^\circ$]')
    print offs

    plt.tight_layout()
    fig.show()
    if light:
        plt.savefig('histStackLight.pdf')
    else:
        plt.savefig('histStack.pdf')
    raw_input('')

# === MAIN ===
def main():
    generate, plot, zstudy = get_args()

    # if plot:
        # vectorDist()
        # return
        # testPlot()
        # return
    #      testFunction()
    #      return

    if zstudy:
        zStudy()
        return

    # clusterTest()
    # return

    if PHASE2:
        pdfOut = 'comptonScatterPhase2'
    else:
        pdfOut = 'comptonScatterPhase1'

    if MULTIPLICITY > 2:
        pdfOut += 'Multi.pdf'
    else:
        pdfOut += '.pdf'

    dS5 = getDict(src='S5', generate=generate, MC=False)
    resultListS5, resultDictS5 = comptonScattering(dS5, l=1, rLim=None, multiplicity=MULTIPLICITY)
    # PHASE2 = False
    # generate = True
    dS5MC = getDict(src='S5', generate=generate, MC=True)
    resultListS5MC, resultDictS5MC = comptonScattering(dS5MC, l=1, rLim=None, multiplicity=MULTIPLICITY)

    '''
    dS2 = getDict(src='S2', generate=generate, MC=False)
    resultListS2, resultDictS2 = comptonScattering(dS2, l=1, rLim=None, multiplicity=MULTIPLICITY)
    dS2MC = getDict(src='S2', generate=generate, MC=True)
    resultListS2MC, resultDictS2MC = comptonScattering(dS2MC, l=1, rLim=None, multiplicity=MULTIPLICITY)

    dS8 = getDict(src='S8', generate=generate, MC=False)
    resultListS8, resultDictS8 = comptonScattering(dS8, l=1, rLim=None, multiplicity=MULTIPLICITY)
    dS8MC = getDict(src='S8', generate=generate, MC=True)
    resultListS8MC, resultDictS8MC = comptonScattering(dS8MC, l=1, rLim=None, multiplicity=MULTIPLICITY)

    # Combine results
    resultList, resultListMC = resultListS5 + resultListS2 + resultListS8, resultListS5MC + resultListS2MC + resultListS8MC
    resultDict, resultDictMC = {}, {}
    for key in set(resultDictS5.keys() + resultDictS2.keys() + resultDictS5.keys()):
        resList, resListMC = [], []
        try:
            resS5, resS5MC = resultDictS5[key], resultDictS5MC[key]
            resList.append( resS5 )
            resListMC.append( resS5MC )
        except:
            pass

        try:
            resS2, resS2MC = resultDictS2[key], resultDictS2MC[key]
            resList.append( resS2 )
            resListMC.append( resS2MC )
        except:
            pass

        try:
            resS8, resS8MC = resultDictS8[key], resultDictS8MC[key]
            resList.append( resS8 )
            resListMC.append( resS8MC )
        except:
            pass

        subDict, subDictMC = {}, {}
        for subKey in resultDictS5[resultDictS5.keys()[0]].keys():
            if subKey == 'position':
                continue

            resAdd, resAddMC = [], []
            for item in resList:
                resAdd += list( item[subKey] )
            for item in resListMC:
                resAddMC += list( item[subKey] )

            subDict[subKey] = resAdd
            subDictMC[subKey] = resAddMC
        subDict['position'] = resultDictS5[key]['position']
        subDictMC['position'] = resultDictS5[key]['position']

        resultDict[key], resultDictMC[key] = subDict, subDictMC
    '''

    # Plot results
    if plot:
        pp = PdfPages( pdfOut )
        plotScattering(resultListS5, 'Th228 @S5', pp)
        plotScattering(resultListS5MC, 'Th228 @S5 (MC)', pp)
        '''
        plotScattering(resultListS2, 'Th228 @S2', pp)
        plotScattering(resultListS2MC, 'Th228 @S2 (MC)', pp)
        plotScattering(resultListS8, 'Th228 @S8', pp)
        plotScattering(resultListS8MC, 'Th228 @S8 (MC)', pp)
        plotScattering(resultList, 'Th228 @S5+S2+S8', pp)
        plotScattering(resultListMC, 'Th228 @S5+S2+S8 (MC)', pp)
        # raw_input('')
        '''
        pp.close()

    
    getEnergyInfo( resultDictS5, resultDictS5MC, smooth=True, src='ThS5' )
    getEnergyInfo( resultDictS5, resultDictS5MC, smooth=True, src='ThS5', light=True, MC=False )

    '''
    getEnergyInfo( resultDictS2, resultDictS2MC, smooth=True, src='ThS2' )
    getEnergyInfo( resultDictS2, resultDictS2MC, smooth=True, src='ThS2', light=True )

    getEnergyInfo( resultDictS8, resultDictS8MC, smooth=True, src='ThS8' )
    getEnergyInfo( resultDictS8, resultDictS8MC, smooth=True, src='ThS8', light=True )

    getEnergyInfo( resultDict, resultDictMC, smooth=True, src='' )
    getEnergyInfo( resultDict, resultDictMC, smooth=True, src='', light=True )
    '''

    # getDensity(d)
    return
    
def clusterTest():
    pdfOut = 'comptonTest.pdf'
    fOutName, srcString = getPath('S5', False)

    EIn = 2614.5

    # Use the same first cluster
    clusterFirst = (150., 3.9, -30.)
    energyFirst = 500.

    # Distance to second cluster
    dist = 10

    # Scattering angles
    thetaList = np.linspace(0., np.pi, 100)

    # Choose two vector components of second gamma vector
    y = 0
    
    d = {}
    dictList = []
    # for theta in thetaList:
    thetaList, thetaEList = [], []
    energyList = []
    vecList = []
    zList = np.linspace(-300., 240., 10000)
    scatterList = []
    for z in zList:
        subDict = {}

        vec1 = np.array( clusterFirst ) - SRC
        # vec2 = np.array([10., 10., (np.cos( theta ) * (dist * np.linalg.norm(vec1)) - 10*vec1[0] - 10*vec1[1]) / vec1[2]])
        # vec2 = np.array( (20., 20., z) )

        # print vec1
        # print (np.cos(theta) * dist * np.linalg.norm(vec1) - vec1[0] - vec1[1]) / vec1[2]
        # clusterSecond = tuple( np.array( clusterFirst ) + np.array(vec2) )
        clusterSecond = (100, 3.9, z)
        vec2 = np.array( clusterSecond ) - np.array( clusterFirst )

        theta = np.arccos(np.dot(vec1, vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))
        thetaList.append( theta )

        energyFirst = EIn * (1 - 1./(1 + (1 - np.cos( theta )) * EIn / 511))
        energySecond = EIn - energyFirst

        p1 = EIn * vec1/np.linalg.norm(vec1)
        p2 = energySecond * vec2/np.linalg.norm(vec2)
        pE = p1 - p2
        print 'Momentum'
        print np.linalg.norm(pE), energyFirst, np.sqrt(EIn**2 + energySecond**2 - 2*EIn*energySecond*np.cos(theta)), np.sqrt( np.linalg.norm(pE)**2 + 511**2 )
        energyList.append( energyFirst )

        thetaE = np.arccos(np.dot(pE, np.array([0, 0, 1.]))/np.linalg.norm(pE))
        print thetaE, np.arctan( (1+EIn)/511*np.tan(theta/2) ) - np.pi/2
        thetaEList.append( thetaE )

        vecList.append( np.arccos(1./energyFirst * (EIn*vec1[2]/np.linalg.norm(vec1) - (EIn-energyFirst)*vec2[2]/np.linalg.norm(vec2))) )

        vecE, vecE2 = np.abs(vecAngle(pE, p1)), np.abs(np.arctan((1+EIn/511)*np.tan(theta/2)) - np.pi/2.)
        if round(vecE,3) != round(vecE2,3):
            print vecE, vecE2
            raw_input('')

        subDict['position'] = [clusterFirst, clusterSecond]
        subDict['energy'] = (energyFirst, energySecond)
        subDict['energyScint'] = EIn

        dictList.append( subDict )

        scatterList.append( [p1, p2, pE] )

        print clusterFirst, clusterSecond, vec1, vec2, theta, np.arccos(np.dot(vec1, vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))), energyFirst, energySecond

    pList = [[(0., 0., 0.)]*3]*len(scatterList)
    plot3dVecs([pList[0]], [scatterList[0]])

    import matplotlib.pyplot as plt
    # plt.hist(thetaList, 100)
    # plt.plot(thetaList, vecList)
    plt.plot(np.array(thetaEList)*180./np.pi, np.array(thetaList)*180./np.pi)
    plt.show()

    for item in dictList:
        print item

    # di = {'position': [(150., 3.9, -30), (139.74, 3.9, -1.81)], 'energy': (2015.75, 598.76), 'energyScint': 2614.5}
    # Di = {}
    # Di[(1, 1, 1)] = [di] 
    Di = {}
    Di[(1, 1, 1)] = dictList
    resultListS5, resultDictS5 = comptonScattering(Di, l=0, rLim=None)

    pp = PdfPages( pdfOut )
    plotScattering(resultListS5, 'Th228 @S5', pp)
    pp.close()

    getEnergyInfo( resultDictS5, resultDictS5, smooth=True, src='ThS5' )
    print resultDictS5

    return

def plot3dVecs(pList, vecList):
    print pList
    print vecList
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    colorList = ['red', 'green', 'blue']
    for j, pElm in enumerate(pList):
        vecElm = vecList[j]
        print vecElm, pElm
        p, vec = np.array( pElm ), np.array( vecElm )
        X, Y, Z = p[:,0], p[:,1], p[:,2]
        U, V, W = vec[:,0], vec[:,1], vec[:,2]

        ax.quiver(X, Y, Z, U, V, W, color=colorList[j])

    fig.show()
    raw_input('')

# === GET DICT ===
# Fill dictionary with events
def getDict(src='S5', generate=False, MC=False):
    fOutName, srcString = getPath(src, MC)

    if not generate and os.path.isfile(fOutName):
        print 'Load data from file...',
        d = cPickle.load(open(fOutName, 'rb'))
        print 'done!'
        print
    else:
        dataChain = getDataChain(fOutName, srcString, MC)
        if MC:
            d = getClusterPos(dataChain, FV, art='mc', multiplicity=MULTIPLICITY)
        else:
            d = getClusterPos(dataChain, FV, art='ms', multiplicity=MULTIPLICITY)

        cPickle.dump(d, open(fOutName, 'wb'))
        print 'done!'
        print

    return d

def getDataChain(fOutName, srcString, MC=False):
    # Get data runs
    preProcessDir = '/home/vault/capm/mppi025h/analysis/preprocess/'
    print 'Generate data...',
    if not MC:
        dataList = []
        for fn in os.listdir(preProcessDir + srcString):
            # runNum = [int(s) for s in fn.split('_') if s.isdigit()][0] 
            # if runNum >= 4500 and runNum <= 4600:
            dataList.append( preProcessDir + srcString + '/%s' % fn )
        dataChain = ps.getChain('dataTree', dataList)

    else:
        dataList = [preProcessDir + srcString]
        dataChain = ps.getChain('mcTree', dataList)

    return dataChain

def getPath(src='S5', MC=False):
    global SRC
    if src == 'S5':
        SRC = SRC_S5
        if not MC:
            if PHASE2:
                srcString = 'S5ThPhase2'
                if MULTIPLICITY > 2:
                    fOutName = 'clusterPosS5MultiPhase2.p'
                else:
                    fOutName = 'clusterPosS5Phase2.p'
            else:
                srcString = 'S5Th'
                if MULTIPLICITY > 2:
                    fOutName = 'clusterPosS5Multi.p'
                else:
                    fOutName = 'clusterPosS5.p'
        else:
            srcString = 'SourceS5_Th228.root'
            if MULTIPLICITY > 2:
                fOutName = 'clusterPosS5MultiMC.p'
            else:
                fOutName = 'clusterPosS5MC.p'

    elif src == 'S2':
        SRC = SRC_S2
        if not MC:
            srcString = 'S2Th/p2_nz'
            fOutName = 'clusterPosS2.p'
        else:
            srcString = 'SourceS2_Th228.root'
            fOutName = 'clusterPosS2MCArt.p'

    elif src == 'S8':
        SRC = SRC_S8
        if not MC:
            srcString = 'S2Th/p2_pz'
            fOutName = 'clusterPosS8.p'
        else:
            srcString = 'SourceS8_Th228.root'
            fOutName = 'clusterPosS8MCArt.p'

    else:
        return False

    return fOutName, srcString

# = COMPTON SCATTERING =
# Discard bins with less than l events
def comptonScattering(d, l=10, rLim=None, multiplicity=2):
    import matplotlib.pyplot as plt

    EIn = 2614.5
    validPercList = []

    # For 2d histograms
    electronList = []

    # Dictionary for results
    dOut = {}

    # Loop over volume bins
    for idx, entry in d.iteritems():
        # Get bin position
        r, phi, z = getEdges(idx)
        print idx, (r, phi, z)

        # Use central position of edges as bin position
        r, phi, z = ( getCenter(r), getCenter(phi)*180./np.pi, getCenter(z) ) 
        # if abs(z + 60) > 5:
        #     continue
        # x, y, z = getCart(r, phi, z)

        # Skip bins with radius < rLim
        if rLim:
            if r > rLim:
                continue

        # Get total number of events and 
        # use only bins with high statistics
        NEvents = len(entry)
        if NEvents < l:
            continue
        
        # Store results in lists
        electronMassList = []
        thetaList = []
        thetaReList = []
        energyList = []
        energyReList = []
        energyTotalList = []
        energyScintList = []
        vecList = []

        # == FILTER EVENTS ==
        # Loop over events
        for val in entry:
            # Get energies of both clusters
            energy = val['energy']
            # Get positions of both clusters
            position = val['position']

            energyScint = val['energyScint'] * energyCorrScint

            if multiplicity > 2:
                energyListEntry, energyReListEntry, energyTotalListEntry, thetaListEntry, thetaReListEntry, vecListEntry = getComptonOrderMulti(position, energy)
            else:
                energyListEntry, energyReListEntry, energyTotalListEntry, thetaListEntry, thetaReListEntry, vecListEntry = getComptonOrder(position, energy, idx)

            if not energyListEntry:
                continue

            energyScintList += [energyScint] * len(energyListEntry)
            energyReList += energyReListEntry
            thetaReList += thetaReListEntry
            energyList += energyListEntry
            energyTotalList += energyTotalListEntry
            thetaList += thetaListEntry
            vecList += vecListEntry

        # Get rid of nan values
        thetaReList = np.array( thetaReList )
        # thetaReList = thetaReList[~np.isnan(np.array(thetaReList))]
        thetaReList = np.nan_to_num( thetaReList )

        NValid = len(energyList)      # Number of valid events

        '''
        # Make histogram
        nTheta, binsTheta, patchesTheta = plt.hist(thetaList, 20)
        # nTheta, binsTheta, patchesTheta = plt.hist(thetaReList, 20, alpha=0.5)
        # nEnergy, binsEnergy, patchesEnergy = plt.hist(energyList, 30)
        # nEnergy, binsEnergy, patchesEnergy = plt.hist(energyReList, 30, alpha=0.5)
        plt.legend(loc='best')
        plt.title('bin position (r=%.2f, phi=%.2f, z=%.2f)' % (r, theta, z))
        plt.show()
        '''

        # == EVALUATE ==
        EPhoto = 2614.5

        electronAngleList = []
        electronAngleReList = []
        electronGammaAngleList = []

        interactionEnergyList = []
        vecGamma1List = []
        vecGamma2AngleList = []
        xList, yList = [], []

        energyListFilt = []
        thetaListFilt = []
        thetaReListFilt = []
        for i in range( len(thetaList) ):
            energyListFilt.append( energyList[i] )
            thetaListFilt.append( thetaList[i] )
            thetaReListFilt.append( float(thetaReList[i]) )

            theta = thetaList[i] * np.pi/180
            gamma1Vec, gamma2Vec = vecList[i]

            # Reconstructed from angle
            energyTransfered = comptonEnergy(EPhoto, theta)
            energyGamma2Re = EPhoto - energyTransfered

            # Measured
            energyGamma2 = energyTotalList[i] - energyList[i]

            # Momentum vectors
            # Reconstructed
            p1VecRe = EPhoto * np.array( gamma1Vec ) / np.linalg.norm( gamma1Vec )
            p2VecRe = energyGamma2Re * np.array( gamma2Vec ) / np.linalg.norm( gamma2Vec )

            # Measured
            p1Vec = energyTotalList[i] * np.array( gamma1Vec ) / np.linalg.norm( gamma1Vec )
            p2Vec = energyGamma2 * np.array( gamma2Vec ) / np.linalg.norm( gamma2Vec )

            vecGamma1List.append( vecAngle(gamma1Vec, [0., 0., 1.]) * 180./np.pi )
            vecGamma2AngleList.append( vecAngle(gamma2Vec, [0., 0., 1.]) * 180./np.pi )
            xList.append( gamma2Vec[0] )
            yList.append( gamma2Vec[1] )

            pElectronRe = p1VecRe - p2VecRe
            pElectron = p1Vec - p2Vec

            # Angle relative to field lines
            fieldVec = np.array( [0., 0., 1.] )
        
            angleElectronRe = vecAngle(fieldVec, pElectronRe) 
            angleElectron = vecAngle(fieldVec, pElectron) 
            interactionEnergyList.append( energyTransfered )

            '''
            print 'theta:', theta
            print 'energyTrans:', energyTransfered
            print 'energyGamma2:', energyGamma2
            print 'energyReal:', energyList[i]
            print 'p1Vec:', p1Vec
            print 'p2Vec:', p2Vec
            print 'pElectron:', pElectron
            print 'angleElectron:', angleElectron * 180./np.pi
            print
            '''

            '''
            if angleElectron > np.pi/2.:
                angleElectron -= np.pi
                angleElectron = abs( angleElectron )
            if angleElectronRe > np.pi/2.:
                angleElectronRe -= np.pi
                angleElectronRe = abs( angleElectronRe )
            '''

            electronAngleList.append( angleElectron*180./np.pi )
            electronAngleReList.append( angleElectronRe*180./np.pi )

            electronGammaAngleList.append( vecAngle(gamma1Vec, pElectron)*180./np.pi )

        electronList += zip([z] * len(electronAngleList), [phi] * len(electronAngleList), thetaListFilt, energyListFilt, electronAngleList, interactionEnergyList, electronAngleReList, electronGammaAngleList, vecGamma1List, vecGamma2AngleList, energyTotalList, energyScintList, xList, yList)

        if NEvents < 1:
            continue

        dRes = {'position': (r, phi, z), 'theta': thetaListFilt, 'thetaRe': thetaReListFilt, 'energy': energyListFilt, 'electronAngle': electronAngleList, 'energyRe': interactionEnergyList, 'electronAngleRe': electronAngleReList, 'electronGammaAngle': electronGammaAngleList, 'vecGamma1Angle': vecGamma1List, 'energyTotal': energyTotalList, 'energyScint': energyScintList}

        dOut[idx] = dRes

        validPerc = float(NValid)/NEvents
        print 'Valid compton events:', validPerc      

        try:
            validPercList.append( ((r, theta, z), validPerc) )
        except:
            pass

    for v in validPercList:
        print v

    return electronList, dOut


def getComptonOrder(position, energy, idx):
    EIn = 2614.5
    energyRes = 100
    energyFirst, energySecond = energy
    if energyFirst == 0 or energySecond == 0:
        return [False] * 6

    # Get rid of events with small energy cluster
    # if energyFirst < 500 or energySecond < 500:
    #   return [False]*6

    # Problem: either use nominal value of photopeak
    # as the total energy or the sum of the energy
    # of the charge clusters.
    energyTotal = EIn # energyFirst + energySecond
    comptonEdge = energyTotal * (1 - 1./(1 + 2.*energyTotal / 511))

    posFirst, posSecond = position

    # Calculation of scattering angle
    vec1 = np.array(posFirst) - np.array(SRC)
    vec2 = np.array(posSecond) - np.array(posFirst)
    thetaClusterRel = vecAngle(vec1, vec2)

    # if not (thetaClusterRel >= 75.*np.pi/180 and thetaClusterRel <= 105.*np.pi/180):
    #    return [False]*6 

    # if abs(energyFirst - 511) < 50 or abs(energyFirst - (EIn - 511)) < 50 or abs(energyFirst - 250) < 50:
    #    return [False]*6 

    # Reset discard-flag
    discard = False

    # Check if cluster sit on top of each other
    posR = np.linalg.norm( abs( np.array(posFirst) - np.array(posSecond) )[:2] )
    posZ = np.abs( posFirst[2] - posSecond[2] )
    if posR < .5 or posZ < 0.01:
        return [False]*6

    # Calculate scattering angle from energy
    # if (energyFirst - comptonEdge > energyRes) or (energySecond - comptonEdge > energyRes):
    #    return [False]*6
    # else:
        # thetaRe = comptonAngle(2614., energyFirst)
    thetaRe = comptonAngle(energyTotal, energyFirst)

    thetaReErr = comptonAngleErr(energyTotal, energyFirst, energyRes)

    # If calculated and measured angle do not
    # match, try to switch cluster positions
    thetaReRot = 0.
    thetaClusterRelRot = 0.
    if (abs(thetaClusterRel - thetaRe) > thetaReErr) or np.isnan(thetaRe):
        vec1 = np.array(posSecond) - np.array(SRC)
        vec2 = np.array(posFirst) - np.array(posSecond)
        # thetaClusterRelRot = np.pi - thetaClusterRel
        thetaClusterRelRot = vecAngle(vec1, vec2) 
        thetaReRot = comptonAngle(energyTotal, energySecond)
        thetaReErrRot = comptonAngleErr(energyTotal, energySecond, energyRes)

        # If those angles also don't match, 
        # check distance of cluster
        if (abs(thetaClusterRelRot - thetaReRot) > thetaReErrRot) or np.isnan(thetaReRot):
            vec = np.array(posSecond) - np.array(posFirst)
            clusterDist = np.linalg.norm( vec ) 
            # Assumption: If cluster are close to each other,
            # there is a large error on the angle due to the error
            # on position reconstruction
            # print clusterDist
            if clusterDist > 15:
                # print 'Discarded       :',
                discard = True

            else:
                if energyFirst > energySecond:
                    # NEvents -= 1
                    # print 'Close           :',
                    discard = True
                else:
                    # NEvents -= 1
                    # print 'Close & switched:',
                    discard = True

        else:
            # print 'Switched        :',
            # Also discard switched events if they belong
            # to another volume bin
            # print idx, getIdx(*getCyl(*posSecond)),
            if idx != getIdx(*getCyl(*posSecond)):
                # print
                # NEvents -= 1
                discard = True
            else:
                # print 'True'
                thetaClusterRel, thetaRe = thetaClusterRelRot, thetaReRot
                # vec1, vec2 = vec2, vec1
                energyFirst, energySecond = energySecond, energyFirst
                # NValid += 1

    # else:
        # NValid += 1
        # print 'Normal          :',

    # discard = False

    # Calculate energy from measured angle
    energyRe = comptonEnergy(EIn, thetaClusterRel)

    # if abs(vecAngle(np.array(vec1), np.array([1.,0.,0.])) - 30.*np.pi/180) > 10.*np.pi/180:
    #   return [False]*6

    if discard:
        # print 'First:', posFirst, 'Second:', posSecond
        return [False]*6
    # if 180./np.pi*thetaRe >= 35 and 180./np.pi*thetaRe <= 38:

    '''
    if (thetaClusterRel * 180./np.pi - 90) < 5:
        print energyFirst, energySecond, 180./np.pi*thetaClusterRel, 180./np.pi * thetaRe, 180./np.pi*thetaClusterRelRot, 180./np.pi * thetaReRot
    '''
    
    # electronMass = comptonGetElectronMass(energyFirst, energySecond, thetaClusterRel)
    # electronMassList.append( electronMass )
    energyList = [ float(energyFirst) ]
    energyReList = [ energyRe ]
    energyTotalList = [ float(energyFirst) + float(energySecond) ]
    # energyScintList.append( energyScint )

    thetaList = [ thetaClusterRel * 180./np.pi ]
    thetaReList = [ thetaRe * 180./np.pi ]
    vecList = [ (vec1, vec2) ]

    return energyList, energyReList, energyTotalList, thetaList, thetaReList, vecList

def getComptonOrderMulti(position, energy):
    import itertools

    energyTotal = 2614.5 # sum(energy) # 2614.5

    # First entry has to be the source
    position = [SRC] + list(position)
    energy = [None] + list(energy)

    # Get indices of entries. Here, the first two
    # have to be the same while the rest consists
    # of permutations
    idxList = [[0, 1] + list(l) for l in list(itertools.permutations(range(2, len(position))))]

    # Loop over the pairs of indices
    for idx in idxList:
        # Create pairs for the entries to check
        # against each other
        triple = [(idx[i], idx[i+1], idx[i+2]) for i in range(len(idx) - 2)]

        # Results are stored in lists. Each scattering
        # in the pair list is interpreted as a single
        # scattering taking place 
        electronMassList = []
        thetaList = []
        thetaReList = []
        energyList = []
        energyReList = []
        energyTotalList = []
        vecList = []

        energyFirst = energyTotal
        for t in triple:
            comptonEdge = energyFirst * (1 - 1./(1 + 2.*energyFirst / 511))

            # First: incoming gamma
            # Second: compton scattering, ET is deposited
            # Third: outgoing gamma
            # print 'Triple:', 
            energySecond = energy[t[1]]
            energyThird = energyFirst - energySecond

            # Energy conservation
            if energySecond > energyFirst:
                break

            posFirst, posSecond, posThird = position[t[0]], position[t[1]], position[t[2]]
            # print 'Energies:', (energyFirst, energySecond, energyThird)
            # print 'Positions:', (posFirst, posSecond, posThird)

            vec1 = np.array(posSecond) - np.array(posFirst)
            vec2 = np.array(posThird) - np.array(posSecond)
            thetaClusterRel = vecAngle(vec1, vec2)

            # Check if cluster sit on top of each other
            # or have the same height. If so, discard
            posR = min( [np.linalg.norm( abs( np.array(posFirst) - np.array(posSecond) )[:2] ), np.linalg.norm( abs( np.array(posFirst) - np.array(posThird) )[:2] ), np.linalg.norm( abs( np.array(posSecond) - np.array(posThird) )[:2] )] )
            posZ = min( [np.abs( posFirst[2] - posSecond[2] ), np.abs( posFirst[2] - posThird[2] ), np.abs(posSecond[2] - posThird[2])] )
            if posR < 3 or posZ < 0.01:
                break

            # Check if the energy is larger than the 
            # Compton edge. Impossible -> discard
            if (energySecond > comptonEdge + 100) or (energyThird > comptonEdge + 100):
                break
            else:
                thetaRe = comptonAngle(energyFirst, energySecond)
            thetaReErr = comptonAngleErr(energyFirst, energySecond, 100)

            # print 'theta', thetaClusterRel
            # print 'thetaRe', (thetaRe, thetaReErr)

            # If geometrical angle and angle from 
            # energy mismatch -> discard
            if (abs(thetaClusterRel - thetaRe) > thetaReErr) or np.isnan(thetaRe):
                continue

            # Energy from geometrical angle.
            # Keep in mind that the energy used here is not
            # the theoretical value as it is not known for
            # multiple scatterings
            energyRe = comptonEnergy(energyFirst, thetaClusterRel)

            # Append to results
            energyList.append( float(energySecond) )
            energyReList.append( energyRe )
            energyTotalList.append( energyFirst )

            thetaList.append( thetaClusterRel * 180./np.pi )
            thetaReList.append( thetaRe * 180./np.pi )

            vecList.append( (vec1, vec2) )

            # For next iteration: Incoming gamma is outgoing
            # gamma of last scattering
            energyFirst = energyThird

            # print 'Result found'
            # print energyList
        # print

    # Check if possible configuration was found
    if len(energyList) == len(position)-2: 
        #print energyList, energyReList, energyTotalList, thetaList, thetaReList, vecList
        #print
        return energyList, energyReList, energyTotalList, thetaList, thetaReList, vecList
    else:
        return [False]*6

# = GET DENSITY =
# Get density of events per bin
# Only show bins with more or equal l events
def getDensity(d, l=60):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    posList = []
    valList = []
    for idx, val in d.iteritems():
        try:
            r, theta, z = getEdges(idx)
        except:
            pass

        # print idx, (r, theta, z)

        # Use central position of edges as bin position
        pos = ( getCenter(r), getCenter(theta), getCenter(z) ) 
        pos = getCart(*pos)

        posList.append( pos )
        valList.append( len(val) )
        # print pos, len(val)
        # print

    posList = np.array( posList )
    valList = [v if v > l else np.nan for v in valList]
    x, y, z = posList[:,0], posList[:,1], posList[:,2]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    s = ax.scatter(x, y, z, c=valList)
    cb = fig.colorbar(s)
    fig.show()

    raw_input('')

# Returns cluster positions and their energies in a dictionary
def getClusterPos(tree, FV, art=None, multiplicity=2):
    import re
    ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
    ROOT.gROOT.cd()

    # Perform cuts
    if art:
        cut = ps.getCut(calibCut=True, energyCut=False, type=art, MC=False, phase2=PHASE2)
    else:
        cut = ps.getCutCombined(calibCut=True, energyCut=True, MC=False, eMin=750, eMax=3500, phase2=PHASE2)

    '''
    if PHASE2:
        cut = re.sub('!EventSummary.isDiagonallyCut\(\) && ', '', cut)
    '''
    print cut
    treeCut = tree.CopyTree( cut )

    # Lists containing the results
    eventEnergyList = []
    eventChargeEnergyList = []
    clusterEnergyList = []
    posList = []

    # Create empty dictionary of lists
    d = {}
    for r in range(1, len(rRange)):
        for theta in range(1, len(thetaRange)):
            for z in range(1, len(zRange)):
                coord = (r, theta, z)
                d[coord] = []

    N = treeCut.GetEntries()

    # Loop over all events
    for i in range( N ):
        treeCut.GetEntry(i)
        es = treeCut.EventSummary

        # Use only events with amount of charge clusters
        # specified in multiplicity
        mul = es.multiplicity
        if mul != multiplicity:
            continue

        # Get cluster positions
        x, y, z = np.array(es.cluster_x), np.array(es.cluster_y), np.array(es.cluster_z)

        # If one cluster position is missing, skip the event
        if (abs( x ) > 900).any() or (abs( y ) > 900).any() or (abs( z ) > 900).any():
            continue

        # Combine positions to list
        pos = zip(x, y, z)
        if not pos:
            continue

        position = pos
        if art == 'ss' or art == 'ms':
            energy = tuple(np.array(es.cluster_energy) * energyCorr)
        else:
            energy = tuple(np.array(es.cluster_energy) * energyCorrMC)
            energyFrac = np.array(energy)/sum(energy)
            energy = ensm.getNewEnergy(resolutionCluster, sum(energy)) * energyFrac
            # energy = [ensm.getNewEnergy(resolutionCluster, np.array( en )) for en in energy]

        if art == 'mc':
            from numpy.random import normal
            energyScint = normal(2614.5, energyScintRes)

        else:
            energyScint = np.array(es.e_scint) * energyCorrScint

        if len(energy) != multiplicity:
            continue

        # Get rotated energy to cut on the photopeak
        if art == 'ss':
            energyRot = np.array( es.energy_ss )
            energySigma = 40
        elif art == 'ms':
            energyRot = np.array( es.energy_ms )
            energySigma = 40
        elif art == 'mc':
            energyRot = ensm.getNewEnergy(resolution, np.array( es.energy_mc ))
            energySigma = 40

            # energyRot = np.array( es.energy_mc )
            # energySigma = 1

        # TODO: Test if it works!
        # energyRot = sum( energy )
        # energySigma = 40

        # MC contains sharp lines, therefore
        # reduce sigma tremendously if MC not smeared
        energyMean = 2614.5
        if energyRot < (energyMean - 1.5*energySigma) or energyRot > (energyMean + 1.5*energySigma):
            continue

        # Loop over cluster tuple and its rotation
        idxOld = (0, 0, 0)  # Does not exist
        en = energy
        for i in range( len(pos) ):
            # Create empty dictionary to store information
            dp = {}

            # rotate for different entries
            if multiplicity > 2:
                position = rotate(pos, i)
                energy = rotate(en, i)

            else:
                # rotate for the second entry
                if i == 1:
                    position = (pos[1], pos[0])
                    energy = (energy[1], energy[0])

            dp['position'] = position
            dp['energy'] = energy
            dp['energyScint'] = energyScint

            # Get first entry of position tuple
            x, y, z = position[0]

            r, theta, z = getCyl(x, y, z)
            idx = getIdx(r, theta, z)

            # Make sure that a cluster tuple and its rotation
            # are not within the same volume element. But, it
            # is okay if they are in separate elements.
            if multiplicity == 2:
                if idx == idxOld:
                     break
                else:
                    d[idx].append( dp )
                    idxOld = idx
            else:
                d[idx].append( dp )

        '''
        wis = [isWithinSphere(np.array(p), np.array([50., 50., 50.]), 60.) for p in pos]
        if not np.array(wis).any():
            continue
        '''

        '''
        pos = list( compress(pos, wis) )
        clusterEnergy = list( compress(np.array( es.cluster_energy ), wis) )
        eventEnergy = np.array( es.energy )
        eventChargeEnergy = np.array( es.e_charge )

        eventEnergyList.append( eventEnergy )
        eventChargeEnergyList.append( eventChargeEnergy )
        clusterEnergyList.append( clusterEnergy )
        posList.append( pos )
        '''

        # h.Fill(scint, charge, 1)
        # en.append( (scint, charge) )

    '''
    for key, val in d.iteritems():
        if len( val ) > 10:
            print key, val
    '''

    '''
    for i in range(20):
        print 'Rotated energy:', eventEnergyList[i]
        print 'Total charge energy:', eventChargeEnergyList[i]
        print 'Cluster energy:', clusterEnergyList[i]
        print 'Cluster positions:', posList[i]
        print
    '''

    '''
    c = ROOT.TCanvas()
    h.Draw('colz')
    raw_input('')
    '''
    
    return d

def rotate(l, n):
    return list(l)[n:] + list(l)[:n]

# === GET ENERGY INFO ===
def logNormModified(x, A, s, mu, sigma, c):
    return -float(A)/(x - s) * np.exp(-(np.log( (s - x)/(s - mu) ) - sigma**2)**2 / (2*sigma**2)) + c

def getElectronAngleInfo(d, dMC):
    return

# Use intersecting boxes to smooth the curve
def getEnergyInfo(d, dMC, smooth=False, src='ThS5', light=False, MC=False):
    # Define areas to merge data in
    # Format: (thetaMin, zMin), (thetaMax, zMax)
    # boxList = [((70, -100), (110, -10)), ((0, 50), (50, 182)), ((130, -182), (180, -50))]
    # boxList = [((15, -182), (45, 182)), ((45, -182), (75, 182)), ((75, -182), (105, 182)), ((105, -182), (135, 182)), ((135, -182), (165, 182))]

    if smooth:
        boxList = []
        boxLow, boxHigh = -182, 182

        boxWidth = 15
        boxRange = (15, 165)
        boxStep = 5
        # boxRange = (0, 180-boxStep)

        for i in range( int(float(boxRange[1] - boxRange[0] - boxWidth) / boxStep) ):
            # Get box edges
            # boxCenter = boxRange[0] + boxStep * i
            # boxLeft = boxCenter - .5*boxStep
            # boxRight = boxCenter + .5*boxStep
            boxLeft = boxRange[0] + boxStep * i
            boxRight = boxLeft + boxWidth
            boxCenter = .5*(boxRight + boxLeft)

            box = ((boxLeft, boxLow), (boxRight, boxHigh))
            boxList.append( box )

    else:
        boxList = [((15, -182), (30, 182)), ((30, -182), (45, 182)), ((60, -182), (75, 182)), ((75, -182), (90, 182)), ((90, -182), (105, 182)), ((105, -182), (120, 182)), ((120, -182), (135, 182)), ((135, -182), (150, 182)), ((150, -182), (165, 182))]

    # Finally, fit and plot
    if light:
        pdfOut = 'electronAngleHists%sLight.pdf' % src
        pdfOutRes = 'electronAngleResults%sLight.pdf' %src
    else:
        pdfOut = 'electronAngleHists%s.pdf' % src
        pdfOutRes = 'electronAngleResults%s.pdf' %src

    pp = PdfPages( pdfOut )
    ppRes = PdfPages( pdfOutRes )

    angleBins = 10

    if light:
        mergeListListMC = mergeBoxData(d if not MC else dMC, boxList, angleBins, True, light)
    else:
        mergeListListMC = mergeBoxData(dMC, boxList, angleBins, True, light)

    mergeListList = mergeBoxData(d if not MC else dMC, boxList, angleBins, False, light)

    # print np.array(mergeListList)
    # print np.shape( mergeListList )

    scatteringAngles = np.linspace(0, 180*(1-1./angleBins), angleBins) + 180./(2*angleBins)
    titleAngle = [r'Scattering angle: $\theta = %.0f^\circ$' % angle for angle in scatteringAngles]

    # Contains results for different scattering angles
    mean2dNorm, mean2dNormErr = [], []
    mean2dNormMC, mean2dNormMCErr = [], []

    mean2d, mean2dMC = [], []
    mean2dErr, mean2dMCErr = [], []

    scatteringAnglesFilt = []
    angleList = [float(box[1][0] + box[0][0])/2 for box in boxList]

    for j in range( len(mergeListList) ):
        # Merge data points in specified boxes
        mergeList = mergeListList[j]
        mergeListMC = mergeListListMC[j]
        # print mergeList[0][:5]
        # print mergeListMC[0][:5]

        mergeList = [[] if not m else m for m in mergeList]
        mergeListMC = [[] if not m else m for m in mergeListMC]
        # if np.array([not bool(m) for m in mergeList]).any() or np.array([not bool(m) for m in mergeListMC]).any():
        #    continue
        scatteringAnglesFilt.append( scatteringAngles[j] ) 

        # Get maximum energy of each histogram using fit
        maxList = []
        if smooth:
            # Get only a few angles
            Nangle = 20
            # print len(mergeList)
            step = int(len(mergeList) // Nangle)
            if not step:
                step = 1
            step = 1
            loopRange = range( 0, len(mergeList), step )

        else:
            loopRange = range( len(mergeList) )
            print 'Break'

        for i in loopRange:
            if not mergeList[i] or not mergeListMC[i]:
                continue
            title = r'Electron angle $\theta_\mathrm{e} = %.2f^\circ$, Scattering angle $\theta = %.2f^\circ$' % (angleList[i], scatteringAngles[j])
            # print title
            # maxList.append( plotHistList( [mergeList[i], mergeListMC[i]], title, False, None, False, light ) )

        # print maxList

        # maxList = np.array(maxList)[:,0], np.array(maxList)[:,1]

        # Plot maximum energy vs electron angle
        # plotXY([angleListFilt] * 2, maxList, r'Electron angle $\theta_e$ [$^\circ$]', r'${E_T}_\textrm{max}$ [keV]', pp)

        # Get mean energy
        # meanList, meanListMC = [(float(np.sum(m))/len(m), float(np.sqrt(np.sum(m)))/len(m)) if m else (np.nan, np.nan) for m in mergeList], [(float(np.sum(m))/len(m), float(np.sqrt(np.sum(m)))/len(m)) if m else (np.nan, np.nan) for m in mergeListMC]
        meanList, meanListMC = [(np.mean(m), np.std(m)/np.sqrt(len(m))) if m else (np.nan, np.nan) for m in mergeList], [(np.mean(m), np.std(m)/np.sqrt(len(m))) if m else (np.nan, np.nan) for m in mergeListMC]
        if meanList and meanListMC:
            try:
                plotXY([angleList] * 2, [meanList, meanListMC], r'Electron angle $\theta_\mathrm{e}$ [$^\circ$]', r'$\overline{E_\mathrm{e}}$ [keV]', ppRes, smooth, titleAngle[j], False, light, False)
            except:
                pass

        mean2d.append( np.array(meanList)[:,0] )
        mean2dMC.append( np.array(meanListMC)[:,0] )
        mean2dErr.append( np.array(meanList)[:,1] )
        mean2dMCErr.append( np.array(meanListMC)[:,1] )

        mean2dNorm.append( np.array(meanList)[:,0]/sum(np.array(meanList)[:,0]) )
        mean2dNormErr.append( [(np.float64(np.sqrt(np.sum(m)))/len(m)) / sum(np.array(meanList)[:,0]) for m in mergeList] )
        mean2dNormMC.append( np.array(meanListMC)[:,0]/sum(np.array(meanListMC)[:,0]) )
        mean2dNormMCErr.append( [(np.float64(np.sqrt(np.sum(m)))/len(m)) / sum(np.array(meanListMC)[:,0]) for m in mergeListMC] )

        # Get mean energy relation
        print maxList
        meanEnList, meanEnListMC = getMean(mergeList), getMean(mergeListMC)
        if maxList:
            plotXY([angleList] * 2, [meanEnList, meanEnListMC], r'Electron angle $\theta_\mathrm{e}$ [$^\circ$]', r'$\frac{\overline{E_\mathrm{e}}}{\overline{E_\mathrm{e}}_\textrm{mean}}$', ppRes, smooth, titleAngle[j], False, light, False)

            # Divide histogram in half and compare areas 
            # Use mean of data as split point
            splitMean = .5*(max(np.array(meanList)[:,0]) + min(np.array(meanList)[:,0]))

            areaList, areaListMC = compareAreas(mergeList, splitMean), compareAreas(mergeListMC, splitMean)
            plotXY([angleList] * 2, [areaList, areaListMC], r'Electron angle $\theta_\mathrm{e}$ [$^\circ$]', r'Area fraction', ppRes, smooth, titleAngle[j], False, light, False)
        
    # Correction of angles
    angleListFilt = [2*angleList[0] - angleList[1]] + angleList
    scatteringAnglesFilt = [2*scatteringAnglesFilt[0] - scatteringAnglesFilt[1]] + scatteringAnglesFilt

    # == Sum over all scattering angles ==
    import itertools
    mergeListFlat = [list(itertools.chain(*np.array(mergeListList)[:,i])) for i in range(len(mergeListList[0]))]
    mergeListMCFlat = [list(itertools.chain(*np.array(mergeListListMC)[:,i])) for i in range(len(mergeListListMC[0]))]

    # Mean energy
    # meanListAll, meanListMCAll = [(float(np.sum(m))/len(m), float(np.sqrt(np.sum(m)))/len(m)) if m else (np.nan, np.nan) for m in mergeListFlat], [(float(np.sum(m))/len(m), float(np.sqrt(np.sum(m)))/len(m)) if m else (np.nan, np.nan) for m in mergeListMCFlat]
    meanListAll, meanListMCAll = [(np.mean(m), np.std(m)/np.sqrt(len(m))) if m else (np.nan, np.nan) for m in mergeListFlat], [(np.mean(m), np.std(m)/np.sqrt(len(m))) if m else (np.nan, np.nan) for m in mergeListMCFlat]
    plotXY([angleList] * 2, [meanListAll, meanListMCAll], r'Electron angle $\theta_\mathrm{e}$ [$^\circ$]', r'$\overline{E_\mathrm{e}}$ [keV]', ppRes, smooth, 'All scattering angles', show=False, light=light, polar=False)
    plotXY([angleList] * 2, [meanListAll, meanListMCAll], r'Electron angle $\theta_\mathrm{e}$ [$^\circ$]', r'$\overline{E_\mathrm{e}}$ [keV]', ppRes, smooth, 'All scattering angles', show=False, light=light, polar=True)

    # Normalized mean energy
    meanEnList, meanEnListMC = getMean(mergeListFlat), getMean(mergeListMCFlat)
    plotXY([angleList] * 2, [meanEnList, meanEnListMC], r'Electron angle $\theta_\mathrm{e}$ [$^\circ$]', r'$\overline{E_\mathrm{e}} / \overline{E_\mathrm{e}}_\textrm{mean}$', ppRes, smooth, 'All scattering angles', show=False, light=light)
    plotXY([angleList] * 2, [meanEnList, meanEnListMC], r'Electron angle $\theta_\mathrm{e}$ [$^\circ$]', r'$\overline{E_\mathrm{e}} / \overline{E_\mathrm{e}}_\textrm{mean}$', ppRes, smooth, 'All scattering angles', show=False, light=light, polar=True)

    # Area fraction
    splitMean = .5*(max(np.array(meanListAll)[:,0]) + min(np.array(meanListAll)[:,0]))

    areaList, areaListMC = compareAreas(mergeListFlat, splitMean), compareAreas(mergeListMCFlat, splitMean)
    plotXY([angleList] * 2, [areaList, areaListMC], r'Electron angle $\theta_\mathrm{e}$ [$^\circ$]', r'Area fraction', ppRes, smooth, 'All scattering angles', show=False, light=light)
    plotXY([angleList] * 2, [areaList, areaListMC], r'Electron angle $\theta_\mathrm{e}$ [$^\circ$]', r'Area fraction', ppRes, smooth, 'All scattering angles', show=False, light=light, polar=True)

    dOut = []

    # Plot energy spectra
    for k in loopRange:
    # for k, m in enumerate(mergeListFlat):
        if not mergeListFlat[k] or not mergeListMCFlat[k]:
            continue
        title = r'Electron angle $\theta_\mathrm{e} = %.2f^\circ$' % angleList[k]
        plotHistList( [mergeListFlat[k], mergeListMCFlat[k]], title, False, pp, False, light, False)
        dOut.append([mergeListFlat[k], mergeListMCFlat[k], angleList[k]])

    if light:
        cPickle.dump(dOut, open('histPlotLight.p', 'wb'))
    else:
        cPickle.dump(dOut, open('histPlot.p', 'wb'))

    # == Each scattering angle ==
    # Electron angle vs. scattering angle
    print mean2d
    mean2d, mean2dMC = np.array( mean2d ), np.array( mean2dMC )
    mean2dNorm, mean2dNormMC = np.array( mean2dNorm ), np.array( mean2dNormMC )

    # Data
    plotThreeOld(mean2dNorm, angleListFilt, scatteringAnglesFilt, xlabel=r'Electron angle $\theta_\mathrm{e}$ [$^\circ$]', ylabel=r'Scattering angle $\theta$ [$^\circ$]', title=r'Data', pp=ppRes, interpolate=False, show=False, project=False)
    plotList([angleList] * len(mean2dNorm), mean2dNorm, yerr=mean2dNormErr, xlabel=r'Electron angle $\theta_\mathrm{e}$ [$^\circ$]', ylabel=r'$\overline{E_\mathrm{e}} / \sum_i\overline{E_{T,i}}$', labelList=[r'$\theta = %.2f$' % angle for angle in scatteringAngles], pp=ppRes, smooth=smooth, show=False)

    # MC
    plotThreeOld(mean2dNormMC, angleListFilt, scatteringAnglesFilt, xlabel=r'Electron angle $\theta_\mathrm{e}$ [$^\circ$]', ylabel=r'Scattering angle $\theta$ [$^\circ$]', title=r'MC', pp=ppRes, interpolate=False, show=False, project=False)
    plotList([angleList] * len(mean2dNormMC), mean2dNormMC, yerr=mean2dNormMCErr, xlabel=r'Electron angle $\theta_\mathrm{e}$ [$^\circ$]', ylabel=r'$\overline{E_\mathrm{e}} / \sum_i\overline{E_{T,i}}$', labelList=[r'$\theta = %.2f$' % angle for angle in scatteringAngles], pp=ppRes, smooth=smooth, show=False)

    # Residuals
    resMean = (np.array(mean2d) - np.array(mean2dMC)) / np.array(mean2d)
    print resMean
    plotThreeOld(resMean, angleListFilt, scatteringAnglesFilt, xlabel=r'Electron angle $\theta_\mathrm{e}$ [$^\circ$]', ylabel=r'Scattering angle $\theta$ [$^\circ$]', title=r'(Data - MC)/Data', pp=ppRes, interpolate=False, show=False, project=False)
    plotList([angleList] * len(resMean), np.array(resMean), xlabel=r'Electron angle $\theta_\mathrm{e}$ [$^\circ$]', ylabel=r'(Data-MC)/Data', labelList=[r'$\theta = %.2f$' % angle for angle in scatteringAngles], pp=ppRes, smooth=smooth, show=False)

    # raw_input('')

    pp.close()
    ppRes.close()

# === MERGE BOX DATA ===
# Bin the data using the scattering angle
def mergeBoxData(d, boxList, binN=12, MC=False, light=False):
    # mergeList = [[] for i in range(len(boxList))]
    angleList = [[[] for j in range(len(boxList))] for i in range(binN)]

    # Loop over volume bins
    for key, val in d.iteritems():
        r, phi, z = val['position']
        # Scattering angle
        thetaList = val['thetaRe']
        # thetaList = val['vecGamma1Angle']

        # Energy of electron
        # energyList = val['energy']

        # When using the total energy, energyTotal in
        # ~ line 300 must be set to energyFirst + energySecond!
        # energyList = val['energyTotal']

        # Scintillation energy (only available for data)
        if light:
            if MC:
                energyList = val['energyTotal']
            else:
                energyList = val['energyScint']
        else:
            # TODO: away with it
            energyList = val['energy']

        # Angle of electron to the z-axis
        electronAngleList = val['electronAngle']

        print 'Key:', key
        print 'Energy:', np.array( energyList )
        print 'Position:', (r, phi, z)
        # Loop over boxes
        for i, box in enumerate(boxList):
            (thetaMin, zMin), (thetaMax, zMax) = box
            subBoxes = None
            if thetaMin < 0:
                subBoxes = [[(180 - thetaMin, zMin), (180, zMax)], [(0, zMin), (thetaMax, zMax)]]
            elif thetaMax > 180:
                subBoxes = [[(thetaMin, zMin), (180, zMax)], [(0, zMin), (thetaMax - 180, zMax)]]

            if subBoxes:
                bList = subBoxes
            else:
                bList = box

            for b in bList:
                (thetaMin, zMin), (thetaMax, zMax) = box
                print 'Box:', box

                if (z >= zMin and z <= zMax):
                    # Loop over electron angle
                    for j, electronAngle in enumerate(electronAngleList):
                        if (electronAngle >= thetaMin) and (electronAngle <= thetaMax):
                            # print i, energyList[j]
                            try:
                                en = energyList[j]
                                if light:
                                    if en > 500:
                                        angleList[np.digitize(thetaList[j], np.linspace(0, 180.*(1-1./binN), binN)) - 1][i].append( en )
                                else:
                                    angleList[np.digitize(thetaList[j], np.linspace(0, 180.*(1-1./binN), binN)) - 1][i].append( en )
                            except:
                                pass

                        # mergeList[i].append( energyList[j] )

    # mergeList now contains all energies within
    # the specified regions.
    return angleList

# === NERR ===
# N has to be a tuple where each element is poisson distributed.
# Function returns error of their fraction.
def Nerr(N0, N1):
    return np.sqrt( (1./N0*np.sqrt(N1))**2 + (float(N1)/(N0**2)*np.sqrt(N0))**2 )

# === COMPARE AREAS ===
def compareAreas(hList, energySplit=None):
    EIn = 2614.5
    if not energySplit:
        comptonEdge = EIn * (1 - 1./(1 + 2.*EIn / 511))
        energySplit = comptonEdge

    areaFracList = []
    for h in hList:
        N = [0, 0]
        for energy in h:
            if energy < energySplit:
                N[0] += 1
            else:
                N[1] += 1

        if N[0]:
            areaFracList.append( (float(N[1]) / N[0], Nerr(*N)) )
        else:
            areaFracList.append( (0., 0.) )

    return areaFracList

# === GET MEAN ===
def getMean(hList):
    meanList, meanListErr = [], []
    meanListTotal = []
    for h in hList:
        meanList.append( np.mean(h) )
        meanListErr.append( np.std(h)/np.sqrt(len(h)) )
        meanListTotal += h

    totalMean = np.float64( np.sum(meanListTotal) ) / len( meanListTotal )
    return [(np.float64(m)/totalMean, np.float64(meanListErr[i])/totalMean) for i, m in enumerate(meanList)]

# === SUPPORT ===
def getLastNonZeroIndex(d, default=None):
    rev = (len(d) - idx for idx, item in enumerate(reversed(d), 1) if item)
    return next(rev, default)

def isWithinSphere(pos, center, radius):
    if np.linalg.norm( np.array( pos ) - np.array( center ) ) > radius:
        return False
    else:
        return True

def getCyl(x, y, z):
    # Get cylindrical coordinates
    r = np.sqrt( x**2 + y**2 )
    theta = np.arctan2(y, x)
    if theta < 0:
        theta += 2*np.pi

    return r, theta, z

def getCart(r, theta, z):
    x = r * np.cos(theta)
    y = r * np.sin(theta)

    return x, y, z

# Function to get the index of the dictionary
def getIdx(r, theta, z):
    rIdx, thetaIdx, zIdx = int(np.digitize(r, rRange)), int(np.digitize(theta, thetaRange)), int(np.digitize(z, zRange))
    if thetaIdx >= thetaBins:
        thetaIdx = thetaBins
    if zIdx >= zBins:
        zIdx = zBins
    return (rIdx, thetaIdx, zIdx)
    
# Function to get the bin edges from an index
def getEdges(idx):
    ir, itheta, iz = idx
    r = (rRange[ir-1], rRange[ir])
    theta = (thetaRange[itheta-1], thetaRange[itheta])
    z = (zRange[iz-1], zRange[iz])

    return r, theta, z

# Returns the angle between v1 and v2
def vecAngle(v1, v2):
    v1, v2 = np.array( v1 ), np.array( v2 )
    return np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1)*np.linalg.norm(v2)))
    # return np.arctan2(np.linalg.norm( np.cross(v1, v2) ), np.sum(v1 * v2))

# Returns the angle between p1 and p2
def pointAngle(p1, p2):
    x, y, z = np.array(p2) - np.array(p1)
    theta = np.arctan2(z/np.sqrt(x**2 + y**2), 1)
    if theta < 0:
        theta += 2*np.pi

    return theta

def getCenter(x):
    return ( x[1] + x[0] ) / 2.

# === COMPTON ===
def comptonGetElectronMass(E, E_, theta):
    # return (1 - np.cos(theta)) * E*E_ / (E - E_)
    return E * (E - E_) / E_ * (1 - np.cos(theta))

def comptonEnergy(E, theta):
    return E * (1 - 1./(1 + E/511. * (1 - np.cos(theta))))

def comptonAngle(E, ET):
    m = 511.
    return np.arccos(1 - float(ET)/((E - ET)*E) * m)

def comptonAngleErr(E, ET, deltaET):
    err = np.sqrt( 511. ) * E / (np.sqrt(ET*(2*E**2 - ET*(2*E + 511))) * abs(ET - E)) * deltaET
    '''
    if err > 40.*np.pi/180:
        return 40.*np.pi/180
    else:
        return err
    '''
    return err

# === PLOT ===
# resultList has the format
# [z, phi, theta, E_T]
# where theta is the angle of the electron
def plotScattering(resultList, title='default', pp=None):
    zList = np.array( resultList )[:,0]
    phiList = np.array( resultList )[:,1]
    thetaScatterList = np.array( resultList )[:,2]
    energyList = np.array( resultList )[:,3]
    thetaElectronList = np.array( resultList )[:,4]
    energyReList = np.array( resultList )[:,5]
    thetaElectronReList = np.array( resultList )[:,6]
    thetaElectronGammaList = np.array( resultList )[:,7]
    vecGamma1List = np.array( resultList )[:,8]
    vecGamma2AngleList = np.array( resultList )[:,9]
    energyTotalList = np.array( resultList )[:,10]
    energyScintList = np.array( resultList )[:,11]
    xList = np.array( resultList )[:,12]
    yList = np.array( resultList )[:,13]

    BINS = (40, 20)
    bins = (40, 40)

    # z vs. electron angle
    # H, xedges, yedges = np.histogram2d(thetaElectronList, zList, bins=[160, 40])
    # H = H.T
    # plotThree(H, xedges, yedges, r"Electron angle [$^\circ$]", r"z [mm]", title, pp=pp, show=False, project=True)
    # d = {'x': thetaElectronList, 'y': zList}
    # cPickle.dump(d, open('cluster_test.p', 'wb'))
    plotThree(thetaElectronList, zList, BINS, r"Electron angle $\theta_\mathrm{e}$ [$^\circ$]", r"z [mm]", title, pp=pp, show=False, project=True)

    # phi vs. electron angle
    # H, xedges, yedges = np.histogram2d(thetaElectronList, phiList, bins=[160, 40])
    # H = H.T
    # plotThree(H, xedges, yedges, r"Electron angle [$^\circ$]", r"$\Phi$ [$^\circ$]", title, pp=pp, project=True)
    plotThree(thetaElectronList, phiList, BINS, r"Electron angle $\theta_\mathrm{e}$ [$^\circ$]", r"$\Phi$ [$^\circ$]", title, pp=pp, project=True)

    # energy vs. electron angle
    # H, xedges, yedges = np.histogram2d(thetaElectronList, energyList, bins=[160, 160])
    # H = H.T
    # plotThree(H, xedges, yedges, r"Electron angle [$^\circ$]", r"$E_T$ [keV]", title, pp=pp, project=True)
    plotThree(thetaElectronList, energyList, bins, r"Electron angle $\theta_\mathrm{e}$ [$^\circ$]", r"$E_\mathrm{e}$ [keV]", title, pp=pp, project=True)

    # energy vs. electron angle (reconstructed)
    # H, xedges, yedges = np.histogram2d(thetaElectronReList, energyReList, bins=[160, 160])
    # H = H.T
    # plotThree(H, xedges, yedges, r"Electron angle $\theta_{e,\, \text{rec}}$ [$^\circ$]", r"$E_{T,\, \text{rec}}$ [keV]", title, pp=pp, project=True)
    plotThree(thetaElectronReList, energyReList, bins, r"Electron angle $\theta_{e,\, \mathrm{rec}}$ [$^\circ$]", r"$E_{T,\, \mathrm{rec}}$ [keV]", title, pp=pp, project=True)

    # energy vs. scattering angle
    # H, xedges, yedges = np.histogram2d(thetaScatterList, energyList, bins=[160, 160])
    # H = H.T
    # plotThree(H, xedges, yedges, r"Scattering angle [$^\circ$]", r"$E_T$ [keV]", title, pp=pp, project=True)
    plotThree(thetaScatterList, energyList, bins, r"Scattering angle $\theta$ [$^\circ$]", r"$E_\mathrm{e}$ [keV]", title, pp=pp, project=True)

    # reconstructed energy vs. scattering angle
    # H, xedges, yedges = np.histogram2d(thetaScatterList, energyReList, bins=[160, 160])
    # H = H.T
    # plotThree(H, xedges, yedges, r"Scattering angle [$^\circ$]", r"$E_{T,\, \text{rec}}$ [keV]", title, pp=pp, project=True)
    plotThree(thetaScatterList, energyReList, bins, r"Scattering angle $\theta$ [$^\circ$]", r"$E_{T,\, \mathrm{rec}}$ [keV]", title, pp=pp, project=True)

    # scattering angle vs. electron angle
    # H, xedges, yedges = np.histogram2d(thetaElectronList, thetaScatterList, bins=[160, 160])
    # H = H.T
    # plotThree(H, xedges, yedges, r"Electron angle [$^\circ$]", "Scattering angle [$^\circ]", title, pp=pp, project=True)
    plotThree(thetaElectronList, thetaScatterList, bins, r"Electron angle $\theta_\mathrm{e}$ [$^\circ$]", r"Scattering angle $\theta$ [$^\circ$]", title, pp=pp, project=True)

    # scattering angle vs. electron angle (reconstructed)
    # H, xedges, yedges = np.histogram2d(thetaElectronReList, thetaScatterList, bins=[160, 160])
    # H = H.T
    # plotThree(H, xedges, yedges, r"Electron angle $\theta_{e,\, \text{rec}}$ [$^\circ$]", r"Scattering angle [$^\circ$]", title, pp=pp, project=True)
    plotThree(thetaElectronReList, thetaScatterList, bins, r"Electron angle $\theta_{e,\, \mathrm{rec}}$ [$^\circ$]", r"Scattering angle $\theta$ [$^\circ$]", title, pp=pp, project=True)

    # gamma1 angle vs. z
    # H, xedges, yedges = np.histogram2d(vecGamma1List, zList, bins=[160, 40])
    # H = H.T
    # plotThree(H, xedges, yedges, r"$\gamma_1$-angle [$^\circ$]", r"$z$ [mm]", title, pp=pp, project=True)
    plotThree(vecGamma1List, zList, BINS, r"$\gamma_1$-angle [$^\circ$]", r"$z$ [mm]", title, pp=pp, project=True)

    # gamma1 angle vs. scattering angle
    # H, xedges, yedges = np.histogram2d(vecGamma1List, thetaScatterList, bins=[160, 160])
    # H = H.T
    # plotThree(H, xedges, yedges, r"$\gamma_1$-angle [$^\circ$]", r"Scattering angle [$^\circ$]", title, pp=pp, project=True)
    plotThree(vecGamma1List, thetaScatterList, bins, r"$\gamma_1$-angle [$^\circ$]", r"Scattering angle $\theta$ [$^\circ$]", title, pp=pp, project=True)

    # gamma2 angle vs. z
    # H, xedges, yedges = np.histogram2d(vecGamma2AngleList, zList, bins=[160, 40])
    # H = H.T
    # plotThree(H, xedges, yedges, r"$\gamma_2$-angle [$^\circ$]", r"$z$ [mm]", title, pp=pp, project=True)
    plotThree(vecGamma2AngleList, zList, BINS, r"$\gamma_2$-angle [$^\circ$]", r"$z$ [mm]", title, pp=pp, project=True)

    # gamma2 angle vs. scattering angle
    # H, xedges, yedges = np.histogram2d(vecGamma2AngleList, thetaScatterList, bins=[160, 160])
    # H = H.T
    # plotThree(H, xedges, yedges, r"$\gamma_2$-angle [$^\circ$]", r"Scattering angle [$^\circ$]", title, pp=pp, project=True)
    plotThree(vecGamma2AngleList, thetaScatterList, bins, r"$\gamma_2$-angle [$^\circ$]", r"Scattering angle $\theta$ [$^\circ$]", title, pp=pp, project=True)

    # gamma2: y vs. x
    # H, xedges, yedges = np.histogram2d(xList, yList, bins=[160, 160])
    # H = H.T
    # plotThree(H, xedges, yedges, r"${\gamma_2}_x$ [mm]", r"${\gamma_2}_y$ [mm]", title, pp=pp, project=True)
    plotThree(xList, yList, bins, r"${\gamma_2}_x$ [mm]", r"${\gamma_2}_y$ [mm]", title, pp=pp, project=True)

    # ET vs. ETRe
    # H, xedges, yedges = np.histogram2d(energyList, energyReList, bins=[160, 160])
    # H = H.T
    # plotThree(H, xedges, yedges, r"$E_T$ [keV]", r"$E_{T,\ \text{rec}}$ [keV]", title, pp=pp, project=True)
    plotThree(energyList, energyReList, bins, r"$E_\mathrm{e}$ [keV]", r"$E_{T,\ \mathrm{rec}}$ [keV]", title, pp=pp, project=True)

    # Electron angle vs. Electron angle (Re)
    # H, xedges, yedges = np.histogram2d(thetaElectronList, thetaElectronReList, bins=[160, 160])
    # H = H.T
    # plotThree(H, xedges, yedges, r"Electron angle $\theta_{e}$ [$^\circ$]", r"Electron angle $\theta_{e,\, \text{rec}}$ [$^\circ$]", title, pp=pp, project=True)
    plotThree(thetaElectronList, thetaElectronReList, bins, r"Electron angle $\theta_{e}$ [$^\circ$]", r"Electron angle $\theta_{e,\, \mathrm{rec}}$ [$^\circ$]", title, pp=pp, project=True)

    # Electron-Gamma angle vs. z
    # H, xedges, yedges = np.histogram2d(thetaElectronGammaList, zList, bins=[160, 40])
    # H = H.T
    # plotThree(H, xedges, yedges, r"Electron angle $\theta_{e,\, \gamma}$ [$^\circ$]", r'z [mm]', title, pp=pp, project=True)
    plotThree(thetaElectronGammaList, zList, BINS, r"Electron angle $\theta_{e,\, \gamma}$ [$^\circ$]", r'z [mm]', title, pp=pp, project=True)

    # EnergyScint vs. z
    enScint, zScint = zip(*[item for item in zip(energyScintList, zList) if abs(item[0] - 2614.5) < 500])
    # H, xedges, yedges = np.histogram2d(zScint, enScint, bins=[40, 160])
    # H = H.T
    # plotThree(H, xedges, yedges, r'z [mm]', r'$E_{\mathrm{scint}}$ [keV]',  title, pp=pp, project=True)
    try:
        plotThree(zScint, enScint, (20, 40), r'z [mm]', r'$E_{\mathrm{scint}}$ [keV]',  title, pp=pp, project=True)
    except:
        pass

    # EnergyScint vs. Electron angle (Re)
    enScint, thetaEl = zip(*[item for item in zip(energyScintList, thetaElectronReList) if abs(item[0] - 2614.5) < 500])

    d = {'x': enScint, 'y': thetaEl}
    cPickle.dump(d, open('cluster_test.p', 'wb'))

    # H, xedges, yedges = np.histogram2d(enScint, thetaEl, bins=[160, 160])
    # H = H.T
    # plotThree(H, xedges, yedges, r"$E_{\text{scint}}$ [keV]", r"Electron angle $\theta_{e,\, \text{rec}}$ [$^\circ$]", title, pp=pp, project=True)
    try:
        plotThree(enScint, thetaEl, bins, r"$E_{\mathrm{scint}}$ [keV]", r"Electron angle $\theta_{\mathrm{e},\, \mathrm{rec}}$ [$^\circ$]", title, pp=pp, project=True)
    except:
        pass

    # EnergyScint vs. Electron angle
    enScint, thetaEl = zip(*[item for item in zip(energyScintList, thetaElectronList) if abs(item[0] - 2614.5) < 500])
    # H, xedges, yedges = np.histogram2d(enScint, thetaEl, bins=[160, 160])
    # H = H.T
    # plotThree(H, xedges, yedges, r"$E_{\text{scint}}$ [keV]", r"Electron angle $\theta_{e}$ [$^\circ$]", title, pp=pp, project=True)
    try:
        plotThree(enScint, thetaEl, bins, r"$E_{\mathrm{scint}}$ [keV]", r"Electron angle $\theta_\mathrm{e}$ [$^\circ$]", title, pp=pp, project=True)
    except:
        pass

    # EnergyTotal vs. Electron angle
    enTot, thetaEl = zip(*[item for item in zip(energyTotalList, thetaElectronList) if abs(item[0] - 2614.5) < 500])
    d2 = {'x': enTot, 'y': thetaEl}
    cPickle.dump(d2, open('cluster_test2.p', 'wb'))

    # H, xedges, yedges = np.histogram2d(enTot, thetaEl, bins=[160, 160])
    # H = H.T
    # plotThree(H, xedges, yedges, r"$E_{\text{charge}}$ [keV]", r"Electron angle $\theta_{e,\, \text{rec}}$ [$^\circ$]", title, pp=pp, project=True)
    try:
        plotThree(enTot, thetaEl, bins, r"$E_{\mathrm{charge}}$ [keV]", r"Electron angle $\theta_\mathrm{e}$ [$^\circ$]", title, pp=pp, project=True)
    except:
        pass

    # EnergyTotal vs. Electron angle
    enTot, thetaEl = zip(*[item for item in zip(energyTotalList, thetaElectronList) if abs(item[0] - 2614.5) < 500])
    d2 = {'x': enTot, 'y': thetaEl}
    cPickle.dump(d2, open('cluster_test2.p', 'wb'))

    # EnergyElectron vs. Electron angle
    try:
        plotThree(energyList, thetaElectronList, bins, r"$E_\mathrm{e}$ [keV]", r"Electron angle $\theta_\mathrm{e}$ [$^\circ$]", title, pp=pp, project=True)
    except:
        pass

    # EnergyTotal vs. EnergyScint
    enTot, enScint = zip(*[item for item in zip(energyTotalList, energyScintList) if abs(item[0] - 2614.5) < 500])
    # H, xedges, yedges = np.histogram2d(enTot, enScint, bins=[160, 160])
    # H = H.T
    # plotThree(H, xedges, yedges, r"$E_{\mathrm{charge}}$ [keV]", r"$E_{\mathrm{scint}}$ [keV]", title, pp=pp, project=True)
    try:
        plotThree(enTot, enScint, bins, r"$E_{\mathrm{charge}}$ [keV]", r"$E_{\mathrm{scint}}$ [keV]", title, pp=pp, project=True)
    except:
        pass

def plotThree(x, y, bins=(80, 20), xlabel='x', ylabel='y', title='default', pp=None, interpolate=True, show=False, project=True, kde=True):
    x, y = np.array(x), np.array(y)

    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib import cm
    from matplotlib import rc
    from copy import copy
    from scipy.stats import kendalltau
    import seaborn as sns
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    # rc('legend', fontsize=14)

    paper_rc = {'lines.linewidth': .5}
    sns.set_context("paper", rc = paper_rc)

    # Normalize
    # h = np.array(h) / np.sum( h )

    colorsData = sns.color_palette('viridis', 256) # sns.cubehelix_palette(20, start=.5, rot=-.75, dark=0.2, light=.7)
    colors = [(1., 1., 1.)] + colorsData # ['#FFFFFFFF'] + colorsData
    cmap = ListedColormap(sns.color_palette(colors)) #.as_hex())
    colorsLine = [[c*1.8 if c*1.8 <= 1 else 1. for c in color] for color in colorsData]
    cmapLines = ListedColormap(sns.color_palette(colorsLine)) # sns.cubehelix_palette(20, start=.5, rot=-.75, dark=0.15, light=.65, as_cmap=True)

    lm = sns.jointplot(x, y, kind='hex', stat_func=None, space=0, cmap=cmap, marginal_kws=dict(bins=40), joint_kws=dict(gridsize=bins)) # , color='b')
    try:
        if kde:
            lm.plot_joint(sns.kdeplot, n_levels=20, cmap=cmapLines)
        else:
            pass
    except:
        pass

    lm.set_axis_labels(xlabel, ylabel, fontsize=14)
    plt.xlim(min(x), max(x))
    plt.ylim(min(y), max(y))
    plt.subplots_adjust(left=0.15, bottom=0.1)
    if show:
        plt.show()

    if pp:
        pp.savefig( lm.fig )

def kdeOverlay(x1, y1, x2, y2, xlabel='x', ylabel='y'):
    x1, y1, x2, y2 = np.array(x1), np.array(y1), np.array(x2), np.array(y2)

    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib import cm
    from matplotlib import rc
    from copy import copy
    from scipy.stats import kendalltau
    import seaborn as sns
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)

    paper_rc = {'lines.linewidth': .3}
    sns.set_context("paper", rc = paper_rc)

    cmapBlue, cmapRed = cmapAlpha('Blues', alpha=1.), cmapAlpha('Reds', alpha=1.)

    cmapBlue = sns.cubehelix_palette(8, start=.5, rot=-.75, reverse=False, as_cmap=True)
    cmapRed = sns.cubehelix_palette(8, start=.5, rot=-.75, reverse=True, as_cmap=True)

    fig, ax = plt.subplots(figsize=(3, 6))
    sns.kdeplot(x1, y1, shade=True, shade_lowest=True, ax=ax, cmap=cmapBlue, n_levels=20, bw=.1) 
    sns.kdeplot(x2, y2, shade=False, shade_lowest=False, ax=ax, cmap='Purples_d', n_levels=20, bw=.1) 

    # ax.set_xlim(min([min(x1), min(x2)]), max([max(x1), max(x2)]))
    ax.set_xlim(-3, 3)
    ax.set_ylim(15, 160)
    # ax.set_ylim(min([min(y1), min(y2)]), max([max(y1), max(y2)]))

    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)

    plt.show()

def cmapAlpha(cmap='Blues', alpha=.2):
    import matplotlib.pylab as pl
    from matplotlib.colors import ListedColormap

    # Choose colormap
    cmap = pl.cm.get_cmap(cmap)

    # Get the colormap colors
    my_cmap = cmap(np.arange(cmap.N))

    # Set alpha
    my_cmap[:,-1] = np.array([alpha]*cmap.N) # np.linspace(0, 1, cmap.N)

    # Create new colormap
    my_cmap = ListedColormap(my_cmap)

    return my_cmap

def plotThreeOld(h, x, y, xlabel='x', ylabel='y', title='default', pp=None, interpolate=True, show=False, project=True):
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib import rc
    from matplotlib.ticker import NullFormatter
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    nullfmt = NullFormatter()

    def getBinWidth(x):
        return abs(x[1] - x[0])

    # Normalize
    h = np.array(h) / np.sum( h )

    # Axis definition
    if project:
        left, width = 0.1, 0.6
        bottom, height = 0.1, 0.6
        left_h = left + width + 0.05
        bottom_h= bottom + height + 0.05
    else:
        left, width = 0.1, 0.9
        bottom, height = 0.1, 0.8

    # Plot areas
    f = plt.figure(figsize=(8, 8))

    rect_main = [left, bottom, width, height]
    axMain = plt.axes(rect_main)
    if project:
        rect_histX = [left, bottom_h, width, 0.15]
        rect_histY = [left_h, bottom, 0.15, height]
        axHistX = plt.axes(rect_histX)
        axHistY = plt.axes(rect_histY)

        # No axis labels
        axHistX.xaxis.set_major_formatter(nullfmt)
        axHistY.yaxis.set_major_formatter(nullfmt)

    # Title
    f.suptitle(title)

    # Main plot
    if interpolate:
        im = axMain.imshow(h, interpolation='bicubic', origin='low', extent=[x[1], x[-1], y[0], y[-1]], aspect='auto', cmap='viridis')
    else:
        im = axMain.imshow(h, interpolation='nearest', origin='low', extent=[x[1], x[-1], y[0], y[-1]], aspect='auto', cmap='viridis')

    axMain.set_xlabel(xlabel, fontsize=16)
    axMain.set_ylabel(ylabel, fontsize=16)
    if project:
        cbar = f.colorbar(im, ax=axHistY, location='right')
    else:
        cbar = f.colorbar(im, ax=axMain, location='right')
    cbar.formatter.set_powerlimits((0, 0))
    cbar.ax.yaxis.set_offset_position('left')
    cbar.update_ticks()
    cbar.ax.set_ylabel('Normalized counts')

    if project:
        # HistX
        h = np.nan_to_num( h )
        x = np.array( x )
        y = np.array( y )

        print x, len(x)
        print h.sum(axis=0), len(h.sum(axis=0))
        print y, len(y)
        print h.sum(axis=1), len(h.sum(axis=1))
        axHistX.bar(x[:-1]+getBinWidth(x)*.5, h.sum(axis=0), width=getBinWidth(x)) 
        axHistX.set_xlim(axMain.get_xlim())

        # HistY
        axHistY.barh(y[:-1]+getBinWidth(y)*.5, h.sum(axis=1), height=getBinWidth(y)) 
        axHistY.set_ylim(axMain.get_ylim())

    if show:
        f.show()
    if pp:
        f.savefig(pp, format='pdf')

def plotHistList(hList, title, show=False, pp=None, fit=False, light=False, kde=False):
    from plot_functions import autoscale_y
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib import rc
    from matplotlib.ticker import NullFormatter
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    from scipy.optimize import curve_fit

    import seaborn as sns
    # sns.set_style('whitegrid', {'axes.grid' : False})
    # sns.set('ticks', {"ytick.major.size": 8})
    sns.set_style('ticks')
    sns.set(style='white')
    # paper_rc = {'lines.linewidth': 1.5, 'lines.markersize': 5, 'lines.markeredgecolor': 'auto', 'lines.markeredgewidth': 1.5}
    # sns.set_context("paper", rc = paper_rc)                                    

    # plt.clf()

    bins = np.linspace(min([min(b) for b in hList]), max([max(b) for b in hList]), 80)
    parInit = [0.6, 2500., 2200., 0.5, 0.0001]
    fitLim = 0

    colorList = ['#af8dc3', '#7fbf7b']
    fitResult = []

    fig, ax = plt.subplots(figsize=(7, 3))
    if light:
        fig.subplots_adjust(top=0.65, bottom=0.2)
    else:
        fig.subplots_adjust(bottom=0.2)

    sns.despine(fig=fig, ax=ax, left=True)
    ax.set_yticks([])

    if light:
        labelList = ['Scintillation', 'Charge']
    else:
        labelList = ['Data', 'MC']

    for i, h in enumerate(hList):
        print np.array( h )
        print colorList[i]
        # n, bins, patches = plt.hist(h, bins, alpha=0.5, label=labelList[i], normed=True) # , color=colorList[i])
        sns.distplot(h, bins=bins, hist=True, kde=True, rug=kde, ax=ax, kde_kws=dict(bw=.05))

        ax.set_xlabel('Energy [keV]', fontsize=16)
        # ax.set_ylabel('Probability')

        # Normalize
        # n = np.array([float(j) for j in n]) / sum(h)

        # Fit
        fit = False
        if fit:
            binsNew, nNew = filterList(bins, n, fitLim)
            try:
                popt, pcov = curve_fit(logNormModified, binsNew, nNew, p0=parInit)
                perr = np.sqrt(np.diag(pcov))
            except:
                popt = parInit
                perr = [0.] * len( popt )

            print popt[2], binsNew[-1]
            if popt[2] > binsNew[-1]:
                mean, meanErr = binsNew[-1], 0.
            else:
                mean, meanErr = popt[2], perr[2]

            binsNew = np.array( binsNew )
            print popt, perr

            # plt.plot(binsNew, logNormModified(binsNew, *popt), label='%s Fit' % labelList[i], color=colorList[i])

            fitResult.append( (mean, meanErr) )

    peakPositions = []
    for h in hList:
        # hNorm = np.array(h) / sum(h)
        hNorm, b = np.histogram(h, bins=bins, normed=True)
        print zip(h, b)
        hMax, binsMax = max(zip(hNorm, b)) # max(zip(bins, hNorm), key=lambda item:item[1])
        print binsMax
        peakPositions.append( (binsMax, hMax) )

    xSize = bins[-1] - bins[0]
    size = ax.get_ylim()
    ySize = size[-1] - size[0]
    print ySize

    print peakPositions
    yPeakData, yPeakMC = peakPositions[0][1], peakPositions[1][1]
    yShift = ySize * 0.05
    if yPeakData == yPeakMC or (abs(yPeakData - yPeakMC) / ySize < 0.05):
        yShiftAdd = ySize*0.03
        if yPeakData > yPeakMC:
            yShiftData = yShift + yShiftAdd
            yShiftMC = yShift - yShiftAdd
        else:
            yShiftData = yShift - yShiftAdd
            yShiftMC = yShift + yShiftAdd
    else:
        yShiftData, yShiftMC = yShift, yShift

    shiftList = [yShiftData, yShiftMC]

    for i, label in enumerate(labelList):
        peak = peakPositions[i]
        ax.text(peak[0]+0.005*xSize, peak[1] + shiftList[i], label, fontsize=12)

    if light:
        ax.set_xlim(2614.5 - 600, 2614.5 + 600)
        autoscale_y(ax, 0.)
    else:
        ax.set_xlim(0, 2500)
    ax.set_ylim(bottom=0.)

    # plt.legend(loc='best')
    plt.title(title)

    if show:
        fig.show()

    if pp:
        # plt.savefig(pp, format='pdf')
        # pp.savefig(fig)
        fig.savefig(pp, format='pdf')

    return fitResult

def plotList(x, y, yerr=None, xlabel='', ylabel='', labelList=None, pp=None, smooth=False, title=None, show=False):
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib import rc
    from matplotlib.ticker import NullFormatter
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)

    if not labelList:
        labelList = [None] * len(x)

    plt.clf()
    f, ax = plt.subplots()

    for i in range( len(x) ):
        x[i], y[i] = np.array(x[i]), np.array(y[i])
        if yerr:
            yerr[i] = np.array( yerr[i] )
            if smooth:
                ax.plot(x[i], y[i], label=labelList[i])
                ax.fill_between(x[i], y[i]-yerr[i], y[i]+yerr[i], alpha=.5)
            else:
                ax.errorbar(x[i], y[i], yerr=yerr[i], label=labelList[i])
        else:
            ax.plot(x[i], y[i], label=labelList[i])

    ax.legend(loc='best', ncol=2)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)

    if title:
        ax.suptitle( title )

    if show:
        f.show()

    if pp:
        f.savefig(pp, format='pdf')

def roundInteger(x, right=False, digit=1):
    flip = False
    if x == 0:
        return 0

    if x < 0:
        x = -x
        flip = True

    '''
    if x <= 1:
        if right:
            return 0
        else:
            if flip:
                return -1
            else:
                return 1
    '''

    y = np.floor(np.log10(np.abs( x )))
    if np.isnan( y ):
        return 0
    else:
        y = int(y)

    # e = y
    e = digit

    if right:
        res = (np.floor(x / (10.**(y-2)))) * 10.**(y-2)
    else:
        res = (int(x / (10.**(y-2))) + 1) * 10.**(y-2)

    if flip:
        return -res
    else:
        return res

def plotXY(x, y, xlabel, ylabel, pp=None, smooth=False, title=None, show=False, light=False, polar=False):
    from fractional_polar import fractional_polar_axes
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib import rc
    from matplotlib.ticker import NullFormatter
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)

    import seaborn as sns
    sns.set_style('whitegrid', {'axes.grid' : False})
    sns.set(style = 'ticks')
    paper_rc = {'lines.linewidth': 1.5, 'lines.markersize': 5, 'lines.markeredgecolor': 'auto', 'lines.markeredgewidth': 1.5}
    sns.set_context("paper", rc = paper_rc)                                    
    if polar:
        from matplotlib import rcParams
        rcParams['axes.titlepad'] = 10

    # Get data
    xData, xMC = np.nan_to_num(np.array( x[0] )), np.nan_to_num(np.array( x[1] ))
    yData, yMC = np.nan_to_num(np.array( y[0] )[:,0]), np.nan_to_num(np.array( y[1] )[:,0])
    res = np.array( (yData - yMC) / yData )
    print res

    xRange = (min([min(xData), min(xMC)]), max([max(xData), max(xMC)]))
    yRange = (min([min(yData), min(yMC)]), max([max(yData), max(yMC)]))
    xSize = xRange[1] - xRange[0]
    xBins = xSize / 5

    ySize = yRange[1] - yRange[0]
    yTop = roundInteger( yRange[1] + .1*ySize, right=False )
    if yRange[0] > 0:
        yBottom = roundInteger( yRange[0] - .1*ySize, right=True, digit=100 )
    else:
        yBottom = roundInteger( yRange[0] - .1*ySize, right=False, digit=100 )

    yBins = float(ySize) / 5
    print yRange, ySize
    print yRange[0] - ySize*.1, yRange[1] + ySize*.1
    print yBottom, yTop, yBins
    print 

    if polar:
        yMin, yMax = yBottom, yTop # yRange

        dataErr, mcErr = np.array( y[0] )[:,1], np.array( y[1] )[:,1]
        resErr = np.sqrt((yMC/yData**2 * dataErr)**2 + (1./yData * mcErr)**2)
        print 'yMin/Max', yMin, yMax
        yData = (yData - yMin) / (yMax - yMin)
        yMC = (yMC - yMin) / (yMax - yMin)
        # res = np.array( (yData - yMC) / yData )

        res[res == -np.inf] = 0
        resRange = (min([min(res), min(res)]), max([max(res), max(res)]))
        resSize = resRange[1] - resRange[0]
        resTop = roundInteger( resRange[1] + .1*resSize, right=False, digit=0.1 )
        if resRange[0] > 0:
            resBottom = roundInteger( resRange[0] - .1*resSize, right=True, digit=0.01 )
        else:
            resBottom = roundInteger( resRange[0] - .1*resSize, right=False, digit=0.01 )

        print resErr
        resMin, resMax = resBottom, resTop
        res = (res - resMin) / (resMax - resMin)
        print 'resMin/Max', resMin, resMax

        # plt.clf()
        # Create figure
        projection = 'polar'
        sharex=False
        #f, (axData, axRes) = plt.subplots(1, 2, figsize=(8, 3), sharex=False, sharey=False, subplot_kw=dict(projection='polar'))
        fData = plt.figure(figsize=(6, 3.5))
        axData = fractional_polar_axes(fData, thlim=(22.5, 152.5), rlim=(yMin, yMax), thlabel=r'', rlabel=ylabel, step=(15, 0.2))
        fData.suptitle(xlabel)
        print 

        fRes = plt.figure(figsize=(6, 3.5))
        if light:
            reslabel = r'(Scintillation - Charge)/Scintillation'
        else:
            reslabel = r'(Data - MC)/Data'
        
        axRes = fractional_polar_axes(fRes, thlim=(22.5, 152.5), thlabel=r'', rlabel=reslabel, rlim=(resMin, resMax), step=(15, 0.2))
        fRes.suptitle(xlabel)

    else:
        f, (axData, axRes) = plt.subplots(2, sharex=True, sharey=False, subplot_kw=dict(projection=None))
        sns.despine(fig=f)
        f.subplots_adjust(wspace=0, hspace=0.1)

    if not polar:
        axRes.axhline(y=0, color='white', linestyle='-', linewidth=1.5)
        axRes.axhline(y=0, color='k', linestyle='--', linewidth=.5)

    # Axis 1 - Data and MC
    if light:
        labelList = ['Scintillation', 'Charge']
    else:
        labelList = ['Data', 'MC']
    # axData.set_yscale('log', nonposy='clip')

    if len(y[0][0]) == 2:
        dataErr, mcErr = np.array( y[0] )[:,1], np.array( y[1] )[:,1]
        size = (xData[0], xData[-1]) # axRes.get_xlim()
        print size
        xSize = size[1] - size[0]
        size = axData.get_ylim()
        ySize = size[1] - size[0]

        yLabelData, yLabelMC = yData[-1], yMC[-1]
        if yLabelData == yLabelMC or (abs(yLabelData - yLabelMC) / ySize < 0.1):
            yShift = ySize*0.2
            if yData[-1] > yMC[-1]:
                yShiftData = yData[-1] + yShift
                yShiftMC = yMC[-1] - yShift
            else:
                yShiftData = yData[-1] - yShift
                yShiftMC = yMC[-1] + yShift
        else:
            yShiftData, yShiftMC = 0, 0

        print yShiftData, yShiftMC
        print xData[-1], xSize

        if smooth:
            axData.plot(xData, yData, label=labelList[0])
            if not polar:
                axData.text(xData[-1]+0.01*xSize, yData[-1], labelList[0], fontsize=12)
            axData.plot(xMC, yMC, label=labelList[1])
            if not polar:
                axData.text(xMC[-1]+0.01*xSize, yMC[-1], labelList[1], fontsize=12)
            axData.fill_between(xData, yData-dataErr, yData+dataErr, alpha=.5)
            axData.fill_between(xMC, yMC-mcErr, yMC+mcErr, alpha=.5)

        else:
            axData.errorbar(xData, yData, yerr=dataErr, label=labelList[0])
            if not polar:
                axData.text(xData[-1]+0.01*xSize, yShiftData, labelList[0], fontsize=12)
            axData.errorbar(xMC, yMC, yerr=mcErr, label=labelList[1])
            if not polar:
                axData.text(xMC[-1]+0.01*xSize, yShiftMC, labelList[1], fontsize=12)
    else:
        axData.plot(xData, yData, label=labelList[0])
        axData.plot(xMC, yMC, label=labelList[1])

    # Axis 2 - residuals
    if len(y[0][0]) == 2:
        if not polar:
            try:
                resErr = np.sqrt((yMC/yData**2 * dataErr)**2 + (1./yData * mcErr)**2)
                print resErr
            except:
                resErr = 0

        if smooth:
            axRes.plot(xData, res)
            axRes.fill_between(xData, res-resErr, res+resErr, alpha=.5)

        else:
            axRes.errorbar(xData, res, yerr=resErr)
    else:
        axRes.plot(xData, res)

    # axRes.set_ylim(-0.5, 0.5)

    # Labels
    if polar:
        axData.legend(loc='best', frameon=False)
    if polar:
        axData.set_xlabel(xlabel, fontsize=16)
        axData.set_title(ylabel, fontsize=16)

    else:
        axData.set_ylabel(ylabel)

    axRes.grid(linestyle='-', linewidth=.2)
    axRes.set_xlabel(xlabel, fontsize=16)
    if light:
        if polar:
            axRes.set_title(r'(Scintillation - Charge)/Scintillation')
        else:
            axRes.set_ylabel(r'(Scintillation - Charge)/Scintillation')
    else:
        if polar:
            axRes.set_title(r'(Data - MC)/Data')
        else:
            axRes.set_ylabel(r'(Data - MC)/Data')

    # plt.xlabel(xlabel)
    # plt.ylabel(ylabel)
    plt.legend()

    if title and not polar:
        f.suptitle( title )

    if pp:
        if polar:
            fData.savefig(pp, format='pdf')
            fRes.savefig(pp, format='pdf')
        else:
            f.savefig(pp, format='pdf')
            # pp.savefig( f )

    if show:
        if polar:
            fData.show()
            fRes.show()
        else:
            f.show()
        raw_input('')

def placeLabels(ax, posList, labelList, textSize=12):
    yBottom, yTop = ax.get_ylim()
    yRange = yTop - yBottom

    # Size of figure in pixels
    size = ax.get_size_inches()*fig.dpi

    # Do all labels fit in?
    # if len(labelList) > yRange / 

def filterList(x, h, val):
    xNew, hNew = zip(*((xEntry, h) for xEntry, h in zip(x, h) if xEntry > val))

    return xNew, hNew

# === Z STUDY ===
# Get MS events, calculate z-distance of cluster
# and show their distance in a histogram
def zStudy(MC=False):
    from itertools import combinations
    import matplotlib.pyplot as plt

    fOutName, srcString = getPath(src='S5', MC=MC)
    dataChain = getDataChain(fOutName, srcString, MC) 

    cut = ps.getCutCombined(calibCut=True, energyCut=True, MC=False, eMin=750, eMax=3500)

    dataChain = dataChain.CopyTree( cut )

    N = dataChain.GetEntries()

    distList = []
    for i in range( N ):
        dataChain.GetEntry(i)
        es = dataChain.EventSummary

        if es.multiplicity != 2:
            continue

        z = np.array(es.cluster_z)
        if len( z ) != 2:
            continue

        # Calculate distance between z-coordinates
        # for the cluster
        distList += [abs(j-i) for i, j in combinations(z, 2)]

    # Plot
    plt.clf()

    plt.hist(distList, 1000)
    plt.show()

# === ARGUMENT PARSER ===
def get_args():
    ap = argparse.ArgumentParser(description=' ')
    ap.add_argument('-g', '--generate', help='Generate files', required=False, action='store_true')
    ap.add_argument('-p', '--plot', help='Show plots', required=False, action='store_true')
    ap.add_argument('-z', '--zstudy', help='Study behaviour of z-coordinates in clusters', required=False, action='store_true')

    args = ap.parse_args()

    return args.generate, args.plot, args.zstudy

# === MAIN ===
if __name__ == '__main__':
    main()

