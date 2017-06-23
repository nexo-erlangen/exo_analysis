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

# Multiplicity
MULTIPLICITY = 3

# Fiducial Volume
FV = (300, 10, 182)

# Source position (S5)
SRC_S5 = (255., 3.9, -30.)
# Source position (S2)
SRC_S2 = (-25.4, 2.3, -292.)
# Source position (S8)
SRC_S8 = (35, 0., 294.8)

# Global source positions
SRC = SRC_S5

# Rotated energy resolution
resolution = [0.714076, 21.9255, 0.]
# Cluster energy resolution
resolutionCluster = [-.000123185451, 74.7268349, 0.0270834510]

# = BIN RANGES =
# Choose binning
zBins, thetaBins = 40, 1
rInit, drInit = 183.2356, 183.2356

# Energy correction
energyCorr = 2614.5 / 2526.97

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

# === MAIN ===
def main():
    generate, plot, zstudy = get_args()

    if zstudy:
        zStudy()
        return

    pdfOut = 'comptonTestMulti.pdf'

    dS5 = getDict(src='S5', generate=generate, MC=False)
    resultListS5, resultDictS5 = comptonScattering(dS5, l=1, rLim=None)
    dS5MC = getDict(src='S5', generate=generate, MC=True)
    resultListS5MC, resultDictS5MC = comptonScattering(dS5MC, l=1, rLim=None)

    '''
    dS2 = getDict(src='S2', generate=generate)
    resultListS2, resultDictS2 = comptonScattering(dS2, l=1, rLim=None)

    dS8 = getDict(src='S8', generate=generate)
    resultListS8, resultDictS8 = comptonScattering(dS8, l=1, rLim=None)

    # Combine results
    resultList = resultListS5 + resultListS2 + resultListS8
    resultDict = {}
    for key in set(resultDictS5.keys() + resultDictS2.keys() + resultDictS5.keys()):
        resList = []
        try:
            resS5 = resultDictS5[key]
            resList.append( resS5 )
        except:
            pass

        try:
            resS2 = resultDictS2[key]
            resList.append( resS2 )
        except:
            pass

        try:
            resS8 = resultDictS8[key]
            resList.append( resS8 )
        except:
            pass

        for subKey in resultDictS5[resultDictS5.keys()[0]].keys():
            subDict = {}
            if subKey == 'position':
                continue

            resAdd = []
            for item in resList:
                resAdd += list( item[subKey] )

            subDict[subKey] = resAdd

        resultDict[key] = subDict
    '''

    # Plot results
    if plot:
        pp = PdfPages( pdfOut )
        plotScattering(resultListS5, 'Th228 @S5', pp)
        plotScattering(resultListS5MC, 'Th228 @S5 (MC)', pp)
        # plotScattering(resultListS2, 'Th228 @S2', pp)
        # plotScattering(resultListS8, 'Th228 @S8', pp)
        # plotScattering(resultList, 'Th228 @S5+S2+S8', pp)
        raw_input('')
        pp.close()

    getEnergyInfo( resultDictS5, resultDictS5MC )

    # getDensity(d)
    return
    
    di = {'position': [(150., 3.9, -30), (139.74, 3.9, -1.81)], 'energy': (2015.75, 598.76)}
    Di = {}
    Di[(1, 1, 1)] = [di] 
    # comptonScattering(Di, l=0, rLim=None)

    return

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
            srcString = 'S5Th'
            fOutName = 'clusterPosS5Multi.p'
        else:
            srcString = 'SourceS5_Th228.root'
            fOutName = 'clusterPosS5MCMulti.p'

    elif src == 'S2':
        SRC = SRC_S2
        srcString = 'S2Th/p2_nz'
        fOutName = 'clusterPosS2.p'

    elif src == 'S8':
        SRC = SRC_S8
        srcString = 'S2Th/p2_pz'
        fOutName = 'clusterPosS8.p'

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

        NValid = 0      # Number of valid events

        # == FILTER EVENTS ==
        # Loop over events
        for val in entry:
            # Get energies of both clusters
            energy = val['energy']

            # Get positions of both clusters
            position = val['position']

            # print 'energy', energy
            # print 'position', position

            energyListEntry, energyReListEntry, energyTotalListEntry, thetaListEntry, thetaReListEntry, vecListEntry = getComptonOrder(position, energy)

            if not energyListEntry:
                continue

            energyReList += energyReListEntry
            thetaReList += thetaReListEntry
            energyList += energyListEntry
            energyTotalList += energyTotalListEntry
            thetaList += thetaListEntry
            vecList += vecListEntry

        # Get rid of nan values
        thetaReList = np.nan_to_num( np.array(thetaReList) )

        # == EVALUATE ==
        electronAngleList = []
        electronAngleReList = []
        interactionEnergyList = []
        vecGamma1List = []
        vecGamma2AngleList = []
        xList, yList = [], []

        for i in range( len(thetaList) ):
            EPhoto = energyTotalList[i]
            theta = thetaList[i] * np.pi/180
            gamma1Vec, gamma2Vec = vecList[i]

            # Assumption: cluster positions are more precise
            # than energy values
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

            electronAngleList.append( angleElectron*180./np.pi )
            electronAngleReList.append( angleElectronRe*180./np.pi )

        electronList += zip([z] * len(electronAngleList), [phi] * len(electronAngleList), thetaList, energyList, electronAngleList, interactionEnergyList, electronAngleReList, vecGamma1List, vecGamma2AngleList, xList, yList)

        if NEvents < 1:
            continue

        dRes = {'position': (r, phi, z), 'theta': thetaList, 'energy': energyList, 'electronAngle': electronAngleList, 'energyRe': interactionEnergyList, 'electronAngleRe': electronAngleReList}
        dOut[idx] = dRes

    return electronList, dOut

def getComptonOrder(position, energy):
    import itertools

    energyTotal = 2614.5

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
            if energySecond > comptonEdge or energyThird > comptonEdge:
                break
            else:
                thetaRe = comptonAngle(energyFirst, energySecond)
            thetaReErr = comptonAngleErr(energyFirst, energySecond, 100)

            # print 'theta', thetaClusterRel
            # print 'thetaRe', (thetaRe, thetaReErr)

            # If geometrical angle and angle from 
            # energy mismatch -> discard
            if (abs(thetaClusterRel - thetaRe) > thetaReErr) or np.isnan(thetaRe):
                break

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

        print idx, (r, theta, z)

        # Use central position of edges as bin position
        pos = ( getCenter(r), getCenter(theta), getCenter(z) ) 
        pos = getCart(*pos)

        posList.append( pos )
        valList.append( len(val) )
        print pos, len(val)
        print

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
    ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
    ROOT.gROOT.cd()

    # Perform cuts
    if art:
        cut = ps.getCut(calibCut=True, energyCut=False, type=art, MC=False)
    else:
        cut = ps.getCutCombined(calibCut=True, energyCut=True, MC=False, eMin=750, eMax=3500)
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
        if art == 'mc':
            energy = tuple(np.array(es.cluster_energy))
            energy = [ensm.getNewEnergy(resolutionCluster, np.array( en )) for en in energy]
        else:
            energy = tuple(np.array(es.cluster_energy) * energyCorr)

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
            energyRot = np.array( es.energy_mc )
            energySigma = 40

        energyMean = 2614.5
        if energyRot < (energyMean - 1.5*energySigma) or energyRot > (energyMean + 1.5*energySigma):
            continue

        # Loop over cluster tuple and its rotation
        for i in range( len(pos) ):
            # Create empty dictionary to store information
            dp = {}

            # rotate for different entries
            position = rotate(pos, i)
            energy = rotate(energy, i)

            dp['position'] = position
            dp['energy'] = energy

            # Get first entry of position tuple
            x, y, z = position[0]

            r, theta, z = getCyl(x, y, z)
            idx = getIdx(r, theta, z)
            d[idx].append( dp )

    return d

def rotate(l, n):
    return l[n:] + l[:n]

# === GET ENERGY INFO ===
def logNormModified(x, A, s, mu, sigma, c):
    return -float(A)/(x - s) * np.exp(-(np.log( (s - x)/(s - mu) ) - sigma**2)**2 / (2*sigma**2)) + c

def getEnergyInfo(d, dMC):
    # Define areas to merge data in
    # Format: (thetaMin, zMin), (thetaMax, zMax)
    # boxList = [((70, -100), (110, -10)), ((0, 50), (50, 182)), ((130, -182), (180, -50))]
    # boxList = [((15, -182), (45, 182)), ((45, -182), (75, 182)), ((75, -182), (105, 182)), ((105, -182), (135, 182)), ((135, -182), (165, 182))]
    boxList = [((15, -182), (30, 182)), ((30, -182), (45, 182)), ((60, -182), (75, 182)), ((75, -182), (90, 182)), ((90, -182), (105, 182)), ((105, -182), (120, 182)), ((120, -182), (135, 182)), ((135, -182), (150, 182)), ((150, -182), (165, 182))]

    # Merge data points in specified boxes
    mergeList = mergeBoxData(d, boxList)
    mergeListMC = mergeBoxData(dMC, boxList)

    # Finally, fit and plot
    pdfOut = 'electronAngleHistsMulti.pdf'
    pp = PdfPages( pdfOut )

    # Get maximum energy of each histogram using fit
    maxList = []
    angleList = [float(box[1][0] + box[0][0])/2 for box in boxList]
    for i in range(len( mergeList )):
        title = r'Electron angle $\Theta_e = %.2f$' % angleList[i]
        maxList.append( plotHistList( [mergeList[i], mergeListMC[i]], title, False, pp ) )

    maxList = np.array(maxList)[:,0], np.array(maxList)[:,1]

    # Plot maximum energy vs electron angle
    plotXY([angleList] * 2, maxList, r'Electron angle $\theta_e$ [$^\circ$]', r'${E_T}_\textrm{max}$ [keV]', pp)

    # Get mean energy
    meanList, meanListMC = [(float(np.sum(m))/len(m), float(np.sqrt(np.sum(m)))/len(m)) for m in mergeList], [(float(np.sum(m))/len(m), float(np.sqrt(np.sum(m)))/len(m)) for m in mergeListMC]
    plotXY([angleList] * 2, [meanList, meanListMC], r'Electron angle $\theta_e$ [$^\circ$]', r'$\overline{E_T}$ [keV]', pp)

    # Get mean energy relation
    print maxList
    meanEnList, meanEnListMC = getMean(mergeList), getMean(mergeListMC)
    plotXY([angleList] * 2, [meanEnList, meanEnListMC], r'Electron angle $\theta_e$ [$^\circ$]', r'$\frac{\overline{E_T}}{\overline{{E_T}_\textrm{mean}}}$', pp)

    # Divide histogram in half and compare areas 
    areaList, areaListMC = compareAreas(mergeList), compareAreas(mergeListMC)
    print areaList
    print areaListMC
    plotXY([angleList] * 2, [areaList, areaListMC], r'Electron angle $\theta_e$ [$^\circ$]', r'Area fraction', pp)
    
    pp.close()

# === MERGE BOX DATA ===
def mergeBoxData(d, boxList):
    mergeList = [[] for i in range(len(boxList))]

    # Loop over volume bins
    for key, val in d.iteritems():
        r, phi, z = val['position']
        thetaList = val['theta']
        energyList = val['energy']
        electronAngleList = val['electronAngle']

        print 'Key:', key
        print 'Energy:', np.array( energyList )
        print 'Position:', (r, phi, z)
        # Loop over boxes
        for i, box in enumerate(boxList):
            (thetaMin, zMin), (thetaMax, zMax) = box
            print 'Box:', box
            if (z >= zMin and z <= zMax):
                # Loop over electron angle
                for j, electronAngle in enumerate(electronAngleList):
                    if (electronAngle >= thetaMin) and (electronAngle <= thetaMax):
                        print i, energyList[j]
                        mergeList[i].append( energyList[j] )

    # mergeList now contains all energies within
    # the specified regions.
    return mergeList

# === NERR ===
# N has to be a tuple where each element is poisson distributed.
# Function returns error of their fraction.
def Nerr(N0, N1):
    return np.sqrt( (1./N0*np.sqrt(N1))**2 + (float(N1)/(N0**2)*np.sqrt(N0))**2 )

# === COMPARE AREAS ===
def compareAreas(hList):
    EIn = 2614.5
    comptonEdge = EIn * (1 - 1./(1 + 2.*EIn / 511))

    areaFracList = []
    for h in hList:
        N = [0, 0]
        for energy in h:
            if energy < comptonEdge/2:
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
    meanList = []
    meanListTotal = []
    for h in hList:
        meanList.append( np.sum(h) )
        meanListTotal += meanList

    totalMean = float( np.sum(meanListTotal) ) / len( meanListTotal )
    return [(float(m)/totalMean, Nerr(m, totalMean)) for m in meanList]

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
    if err > 40.*np.pi/180:
        return 40.*np.pi/180
    else:
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
    vecGamma1List = np.array( resultList )[:,7]
    vecGamma2AngleList = np.array( resultList )[:,8]
    xList = np.array( resultList )[:,9]
    yList = np.array( resultList )[:,10]

    # z vs. electron angle
    H, xedges, yedges = np.histogram2d(thetaElectronList, zList, bins=[160, 40])
    H = H.T
    plotThree(H, xedges, yedges, r'Electron angle [$^\circ$]', r'z [mm]', title, pp=pp)

    # phi vs. electron angle
    H, xedges, yedges = np.histogram2d(thetaElectronList, phiList, bins=[160, 40])
    H = H.T
    plotThree(H, xedges, yedges, r'Electron angle [$^\circ$]', r'$\Phi$ [$^\circ$]', title, pp=pp)

    # energy vs. electron angle
    H, xedges, yedges = np.histogram2d(thetaElectronList, energyList, bins=[160, 160])
    H = H.T
    plotThree(H, xedges, yedges, r'Electron angle [$^\circ$]', r'$E_T$ [keV]', title, pp=pp)

    # energy vs. electron angle (reconstructed)
    H, xedges, yedges = np.histogram2d(thetaElectronReList, energyReList, bins=[160, 160])
    H = H.T
    plotThree(H, xedges, yedges, r'Electron angle $\theta_{e,\, \text{rec}}$ [$^\circ$]', r'$E_{T,\, \text{rec}}$ [keV]', title, pp=pp)

    # energy vs. scattering angle
    H, xedges, yedges = np.histogram2d(thetaScatterList, energyList, bins=[160, 160])
    H = H.T
    plotThree(H, xedges, yedges, r'Scattering angle [$^\circ$]', r'$E_T$ [keV]', title, pp=pp)

    # reconstructed energy vs. scattering angle
    H, xedges, yedges = np.histogram2d(thetaScatterList, energyReList, bins=[160, 160])
    H = H.T
    plotThree(H, xedges, yedges, r'Scattering angle [$^\circ$]', r'$E_{T,\, \text{rec}}$ [keV]', title, pp=pp)

    # scattering angle vs. electron angle
    H, xedges, yedges = np.histogram2d(thetaElectronList, thetaScatterList, bins=[160, 160])
    H = H.T
    plotThree(H, xedges, yedges, r'Electron angle [$^\circ$]', r'Scattering angle [$^\circ]', title, pp=pp)

    # scattering angle vs. electron angle (reconstructed)
    H, xedges, yedges = np.histogram2d(thetaElectronReList, thetaScatterList, bins=[160, 160])
    H = H.T
    plotThree(H, xedges, yedges, r'Electron angle $\theta_{e,\, \text{rec}}$ [$^\circ$]', r'Scattering angle [$^\circ$]', title, pp=pp)

    # gamma1 angle vs. z
    H, xedges, yedges = np.histogram2d(vecGamma1List, zList, bins=[160, 40])
    H = H.T
    plotThree(H, xedges, yedges, r'$\gamma_1$-angle [$^\circ$]', r'$z$ [mm]', title, pp=pp)

    # gamma1 angle vs. scattering angle
    H, xedges, yedges = np.histogram2d(vecGamma1List, thetaScatterList, bins=[160, 160])
    H = H.T
    plotThree(H, xedges, yedges, r'$\gamma_1$-angle [$^\circ$]', r'Scattering angle [$^\circ$]', title, pp=pp)

    # gamma2 angle vs. z
    H, xedges, yedges = np.histogram2d(vecGamma2AngleList, zList, bins=[160, 40])
    H = H.T
    plotThree(H, xedges, yedges, r'$\gamma_2$-angle [$^\circ$]', r'$z$ [mm]', title, pp=pp)

    # gamma2 angle vs. scattering angle
    H, xedges, yedges = np.histogram2d(vecGamma2AngleList, thetaScatterList, bins=[160, 160])
    H = H.T
    plotThree(H, xedges, yedges, r'$\gamma_2$-angle [$^\circ$]', r'Scattering angle [$^\circ$]', title, pp=pp)

    # gamma2: y vs. x
    H, xedges, yedges = np.histogram2d(xList, yList, bins=[160, 160])
    H = H.T
    plotThree(H, xedges, yedges, r'${\gamma_2}_x$ [mm]', r'${\gamma_2}_y$ [mm]', title, pp=pp)

def plotThree(h, x, y, xlabel='x', ylabel='y', title='default', pp=None):
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib import rc
    from matplotlib.ticker import NullFormatter
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    nullfmt = NullFormatter()

    def getBinWidth(x):
        return abs(x[1] - x[0])

    # Axis definition
    left, width = 0.1, 0.6
    bottom, height = 0.1, 0.6
    left_h = left + width + 0.05
    bottom_h= bottom + height + 0.05

    # Plot areas
    rect_main = [left, bottom, width, height]
    rect_histX = [left, bottom_h, width, 0.2]
    rect_histY = [left_h, bottom, 0.15, height]

    f = plt.figure(figsize=(8, 8))
    axMain = plt.axes(rect_main)
    axHistX = plt.axes(rect_histX)
    axHistY = plt.axes(rect_histY)

    # No axis labels
    axHistX.xaxis.set_major_formatter(nullfmt)
    axHistY.yaxis.set_major_formatter(nullfmt)

    # Title
    f.suptitle(title)

    # Main plot
    im = axMain.imshow(h, interpolation='bicubic', origin='low', extent=[x[1], x[-1], y[0], y[-1]], aspect='auto', cmap='viridis')

    axMain.set_xlabel(xlabel)
    axMain.set_ylabel(ylabel)
    cbar = f.colorbar(im, ax=axHistY, location='right', format='%.0e')
    cbar.ax.set_ylabel('Counts')

    # HistX
    h = np.nan_to_num( h )
    axHistX.bar(x[:-1]+getBinWidth(x)*.5, h.sum(axis=0), width=getBinWidth(x)) 
    axHistX.set_xlim(axMain.get_xlim())

    # HistY
    axHistY.barh(y[:-1]+getBinWidth(y)*.5, h.sum(axis=1), height=getBinWidth(y)) 
    axHistY.set_ylim(axMain.get_ylim())

    f.show()
    if pp:
        f.savefig(pp, format='pdf')

def plotHistList(hList, title, show=False, pp=None):
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit

    plt.clf()

    bins = np.linspace(min([min(b) for b in hList]), max([max(b) for b in hList]), 80)
    parInit = [0.6, 2500., 2200., 0.5, 0.0001]
    fitLim = 0

    colorList = ['#af8dc3', '#7fbf7b']
    fitResult = []
    labelList = ['Data', 'MC']
    for i, h in enumerate(hList):
        print np.array( h )
        print colorList[i]
        n, bins, patches = plt.hist(h, bins, alpha=0.5, label=labelList[i], normed=True, color=colorList[i])

        # Normalize
        # n = np.array([float(j) for j in n]) / sum(h)

        # Fit
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
        plt.plot(binsNew, logNormModified(binsNew, *popt), label='%s Fit' % labelList[i], color=colorList[i])

        fitResult.append( (mean, meanErr) )

    plt.legend(loc='best')
    plt.title(title)

    if show:
        plt.show()

    if pp:
        plt.savefig(pp, format='pdf')

    return fitResult

def plotXY(x, y, xlabel, ylabel, pp=None):
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib import rc
    from matplotlib.ticker import NullFormatter
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)

    plt.clf()
    # Create figure
    f, (axData, axRes) = plt.subplots(2, sharex=True, sharey=False)
    f.subplots_adjust(wspace=0, hspace=0.1)

    # Axis 1 - Data and MC
    # labelList = ['Data', 'MC']
    # axData.set_yscale('log', nonposy='clip')

    xData, xMC = np.array( x[0] ), np.array( x[1] )
    yData, yMC = np.array( y[0] )[:,0], np.array( y[1] )[:,0]
    if len(y[0][0]) == 2:
        dataErr, mcErr = np.array( y[0] )[:,1], np.array( y[1] )[:,1]
        axData.errorbar(xData, yData, yerr=dataErr, label='Data')
        axData.errorbar(xMC, yMC, yerr=mcErr, label='MC')
    else:
        axData.plot(xData, yData, label='Data')
        axData.plot(xMC, yMC, label='MC')

    # Axis 2 - residuals
    res = (yData - yMC) / yData
    if len(y[0][0]) == 2:
        resErr = np.sqrt((yMC/yData**2 * dataErr)**2 + (1./yData * mcErr)**2)
        axRes.errorbar(xData, res, yerr=resErr)
    else:
        axRes.plot(xData, res)
    axRes.axhline(y=0, color='k', linestyle='--', linewidth=.5)
    # axRes.set_ylim(-0.5, 0.5)

    # Labels
    axData.legend(loc='best')
    axData.set_ylabel(ylabel)

    axRes.grid(color='k', linestyle='-', linewidth=.2)
    axRes.set_xlabel(xlabel)
    axRes.set_ylabel(r'(Data - MC)/Data')

    # plt.xlabel(xlabel)
    # plt.ylabel(ylabel)
    # plt.legend()

    if pp:
        f.savefig(pp, format='pdf')

    f.show()

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
    ap.add_argument('-g', '--generate', help='Settings to load', required=False, action='store_true')
    ap.add_argument('-p', '--plot', help='Show plots', required=False, action='store_true')
    ap.add_argument('-z', '--zstudy', help='Study behaviour of z-coordinates in clusters', required=False, action='store_true')

    args = ap.parse_args()

    return args.generate, args.plot, args.zstudy

# === MAIN ===
if __name__ == '__main__':
    main()

