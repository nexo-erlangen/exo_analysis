#!/usr/bin/env python
import argparse
import multiprocessing
import numpy as np
from numpy.random import normal, uniform
from scipy import spatial
from scipy.optimize import curve_fit
import sys
import os, os.path
import time
import cPickle
import random

import columnar_recombination as cr

# Doesn't work when run on cluster
try:
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
except:
    pass

def main():
    trackNumber, multi, combineresults, genscripts, single, randomtrajectories, parameter_b = get_args()

    global gifOut
    gifOut = 'cr/cr_gif_track'
    global pickleOut
    pickleOut = 'cr/cr_track.p'

    GAUSS = True
    INTERRUPT = 100
    SHOW = True
    SAVE = False
    PERPENDICULAR = True
    REPEAT = 3

    # Constants
    if parameter_b:
        b = parameter_b*1.e-9
    else:
        b = 50.e-9
        
    # b = 1.e-8
    A = 49.e-9
    vD = 1705
    Dt = 55.e-4
    Dl = 0.1 * Dt
    D = [Dt, Dt, Dl]

    # Time
    dt = 1.e-12
    tt = 1.e-9

    if genscripts:
        jobs = (0, 200)
        generateScripts('cr_job.sh', 'cr_exec_job.sh', jobs[0], jobs[1], lima=True, single=single, b=int(b*1.e9))
        return

    if combineresults:
        if single:
            bList = [50, 100, 150, 200, 250, 300]
            muList, sigmaList = [], []
            for bChoice in bList: 
                mu, sigma = combineResultsSingle('cr/track_results_single', bChoice)
                muList.append( mu ), sigmaList.append( sigma )

            plotSingleResults(bList, muList, sigmaList)

        else:
            combineResults('cr/track_results')
        return

    # trajectories = getTrajectories('cr/pe-trajectories.dat')
    # plotTrajectories(trajectories, [0, 1000], proj='xy')
    # plotTrajectories(trajectories, [0, 1000], proj='xz')
    # plotTrajectories(trajectories, [0, 1000], proj='yz')
    # raw_input('')

    if randomtrajectories:
        trajectories = getTrajectoriesRandom(None)
        cPickle.dump(trajectories, open('cr/random_trajectories.p', 'wb'))
        return

    # === SINGLE BATCH ===
    if single:
        trajectories = cPickle.load(open('cr/random_trajectories.p', 'rb'))
        # plotTrajectories(trajectories, [0, 1000], proj='xy')
        # plotTrajectories(trajectories, [0, 1000], proj='xz')
        # plotTrajectories(trajectories, [0, 1000], proj='yz')
        # raw_input('')
        # return

        trajectory = trajectories[trackNumber]
        getResultsSingle(b, A, D, vD, trajectory, dt, tt, gauss=GAUSS, interrupt=INTERRUPT, repeat=REPEAT, resultOut='cr/track_results_single/track_%d_b_%dnm.p' % (trackNumber, int(b*1.e9)))
        return 

    # === SINGLE ===
    # Need trajectories for the following functions
    trajectories = getTrajectories('cr/pe-trajectories.dat')
    if trackNumber <= len(trajectories):
        trajectory = trajectories[trackNumber]
    else:
        print 'Tried to access non-existing track!'
        return

    # Process trajectory once
    if not multi:
        generate(b, A, D, vD, trajectory, dt=dt, tt=tt, gauss=GAUSS, interrupt=INTERRUPT, perp=PERPENDICULAR, out_q=None, show=SHOW, save=SAVE)
        return

    # === MULTI BATCH ===
    # Process trajectory multiple times,
    # including the perpendicular and parallel case
    else:
        getResults(b, A, D, vD, trajectory, dt, tt, gauss=GAUSS, interrupt=INTERRUPT, repeat=REPEAT, resultOut='cr/track_results/track_%d.p' % trackNumber)
        return

def getResults(b, A, D, vD, trajectory, dt, tt, gauss, interrupt, repeat=1, resultOut='cr/track_results/std.p'):
    manager = multiprocessing.Manager()
    out_q = manager.Queue()

    procs = []
    # Start multiple processes
    for i in range(repeat):
        randIntPerp, randIntPar = np.random.randint(1, 256, 2)
        # Perpendicular
        pPerp = multiprocessing.Process(target=generate, args=(b, A, D, vD, trajectory, dt, tt, gauss, interrupt, False, out_q, False, False, randIntPerp))
        # Parallel
        pPar = multiprocessing.Process(target=generate, args=(b, A, D, vD, trajectory, dt, tt, gauss, interrupt, True, out_q, False, False, randIntPar))

        procs += [pPerp, pPar]
        pPerp.start(), pPar.start()

    resultDict = {}
    percListListPerp = []
    percListListPar = []
    # timeList is the same for all
    for i in range(repeat * 2):
        resultDict.update( out_q.get() )
        timeList, percList = resultDict['t'], resultDict['p']
        if resultDict['perp']:
            percListListPerp.append( percList )
        else:
            percListListPar.append( percList )

    # Wait for processes to finish
    for p in procs:
        p.join()

    # Get mean and std for the data
    if repeat > 1: 
        percListPerp = [(np.mean(p), np.std(p)) for p in zip(*percListListPerp)]
        percListPar = [(np.mean(p), np.std(p)) for p in zip(*percListListPar)]

    outDict = {'t': timeList, 'pPerp': percListPerp, 'pPar': percListPar}
    cPickle.dump(outDict, open(resultOut, 'wb'))

def getResultsSingle(b, A, D, vD, trajectory, dt, tt, gauss, interrupt, repeat=1, resultOut='cr/track_results_single/std.p'):
    manager = multiprocessing.Manager()
    out_q = manager.Queue()

    procs = []
    # Start multiple processes
    for i in range(repeat):
        randInt = np.random.randint(1, 256, 1)[0]
        p = multiprocessing.Process(target=generate, args=(b, A, D, vD, trajectory, dt, tt, gauss, interrupt, True, out_q, False, False, randInt))

        procs.append( p )
        p.start()

    resultDict = {}
    percListList = []
    # timeList is the same for all
    for i in range(repeat):
        resultDict.update( out_q.get() )
        timeList, percList = resultDict['t'], resultDict['p']
        percListList.append( percList )

    # Wait for processes to finish
    for p in procs:
        p.join()

    # Get mean and std for the data
    if repeat > 1: 
        percList = [(np.mean(p), np.std(p)) for p in zip(*percListList)]

    outDict = {'t': timeList, 'b': b, 'p': percList}
    cPickle.dump(outDict, open(resultOut, 'wb'))

def combineResults(resDir):
    # Figure and axis for recombination vs. time
    # Two axes are created - one for the plot, one for the colorbar
    figPerp = plt.figure()
    axPerp = figPerp.add_axes([0.1, 0.1, 0.7, 0.8])
    axPerpCBar = figPerp.add_axes([0.85, 0.1, 0.05, 0.8])

    figPar = plt.figure()
    axPar = figPar.add_axes([0.1, 0.1, 0.7, 0.8])
    axParCBar = figPar.add_axes([0.85, 0.1, 0.05, 0.8])

    # Number of files in the directory
    nameList = [name for name in os.listdir(resDir) if os.path.isfile('%s/%s' % (resDir, name))]
    N = len( [name for name in nameList if '.p' in name] )
    print N

    # Add color bar
    # cr.getColorBar(axPerpCBar, 1, N, 'Track #')
    # cr.getColorBar(axParCBar, 1, N, 'Track #')

    # Loop over all files in directory
    i = 0

    # Store recombination results in list
    resListPerp, resListPar = [], []
    for fn in os.listdir(resDir):

        fn = '%s/%s' % (resDir, fn)
        if os.path.isfile(fn) and fn[-2:] == '.p':
            f = cPickle.load(open(fn, 'rb'))
            trackNumber = [int(s) for s in fn.split() if s.isdigit()]
            timeList, percListPerp, percListPar = f['t'], f['pPerp'], f['pPar']
            if len(timeList) != len(percListPerp) or len(timeList) != len(percListPar):
                continue

            # Fit curve to percList
            poptPerp, pcovPerp = curve_fit(fitFunc, timeList, np.array(percListPerp)[:,0])
            perrPerp = np.sqrt( np.diag(pcovPerp) )
            poptPar, pcovPar = curve_fit(fitFunc, timeList, np.array(percListPar)[:,0])
            perrPar = np.sqrt( np.diag(pcovPar) )

            # Get plot data for fits
            timeListFit = timeList
            percListPerpFit, percListParFit = fitFunc(timeListFit, *poptPerp), fitFunc(timeListFit, *poptPar)

            # Get last entries of perc lists for the 
            # recombination results
            # resListPerp.append( percListPerp[-1] ), resListPar.append( percListPar[-1] )
            resListPerp.append( (poptPerp[0], perrPerp[0]) ), resListPar.append( (poptPar[0], perrPar[0]) )

            color = cr.getColor(N, i)
            # Perpendicular
            # cr.plotFigure(figPerp, axPerp, timeList, np.array(percListPerp)*100, xlabel=r'time [ns]', ylabel=r'Recombination [%]', label=None, title='Perpendicular', color=color)
            # cr.plotFigure(figPerp, axPerp, timeListFit, np.array(percListPerpFit)*100, xlabel=r'time [ns]', ylabel=r'Recombination [%]', label=None, title='Perpendicular', color=color)
            # Parallel
            # cr.plotFigure(figPar, axPar, timeList, np.array(percListPar)*100, xlabel=r'time [ns]', ylabel=r'Recombination [%]', label=None, title='Parallel', color=color)
            # cr.plotFigure(figPar, axPar, timeListFit, np.array(percListParFit)*100, xlabel=r'time [ns]', ylabel=r'Recombination [%]', label=None, title='Parallel', color=color)

            # Increment loop counter
            i += 1

    # Plot the results for all trajectories
    '''
    if len(resListPerp[0]) == 2:
        # True if (mean, std)
        resListPerp, resListPar = np.array(resListPerp)[:,0], np.array(resListPar)[:,0]
    '''

    # Transform to percent
    resListPerp = np.array(resListPerp) * 100
    resListPar = np.array(resListPar) * 100

    figResPerp, axResPerp = plt.subplots()
    KDE(figResPerp, axResPerp, resListPerp, smooth=5, xlabel=r'Recombination [%]', ylabel=r'Probability', label=None, title=r'Perpendicular', showDetail=True, fit=True)

    figResPar, axResPar = plt.subplots()
    KDE(figResPar, axResPar, resListPar, smooth=5, xlabel=r'Recombination [%]', ylabel=r'Probability', label=None, title=r'Parallel', showDetail=True, fit=True)
    raw_input('')

    figRes, axRes = plt.subplots()
    KDE(figRes, axRes, resListPerp, smooth=5, xlabel=r'Recombination [%]', ylabel=r'Probability', label='Perpendicular', title=None, showDetail=False)
    KDE(figRes, axRes, resListPar, smooth=5, xlabel=r'Recombination [%]', ylabel=r'Probability', label='Parallel', title=None, showDetail=False)

    figHist, axHist = plt.subplots()
    plotHistogram(figHist, axHist, np.array(resListPerp)[:,0], xlabel=r'Recombination [%]', ylabel='Probability', label='Perpendicular', title=None)
    plotHistogram(figHist, axHist, np.array(resListPar)[:,0], alpha=.7, xlabel=r'Recombination [%]', ylabel='Probability', label='Parallel', title=None)

    meanPerp = np.average(resListPerp[:,0], weights=1./np.square(resListPerp[:,1])) # np.mean( resListPerp[:,0] )
    stdPerp = np.sqrt(1./np.sum(1./np.square(resListPerp[:,1]))) # np.sum( np.square(resListPerp[:,1]) / len(resListPerp) )
    meanPar = np.average(resListPar[:,0], weights=1./np.square(resListPar[:,1])) # np.mean( resListPar[:,0] )
    stdPar = np.sqrt(1./np.sum(1./np.square(resListPar[:,1]))) # np.sum( np.square(resListPar[:,1]) / len(resListPar) )

    print 'Mean recombination: %f +/- %f (perpendicular), %f +/ - %f (parallel)' % (meanPerp, stdPerp, meanPar, stdPar)
    print np.mean(resListPerp[:,0]), np.mean(resListPar[:,0])
    print np.std(resListPerp[:,0]), np.std(resListPar[:,0])
    raw_input('')

def combineResultsSingle(resDir, bChoice):
    import re
    # Figure and axis for recombination vs. time
    # Two axes are created - one for the plot, one for the colorbar
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.7, 0.8])
    axCBar = fig.add_axes([0.85, 0.1, 0.05, 0.8])

    # Number of files in the directory
    nameList = [name for name in os.listdir(resDir) if os.path.isfile('%s/%s' % (resDir, name))]
    N = len( [name for name in nameList if '.p' in name] )
    print N

    # Add color bar
    # cr.getColorBar(axCBar, 1, N, 'Track #')

    # Loop over all files in directory
    i = 0

    # Store recombination results in list
    resList = []
    for fn in os.listdir(resDir):
        fn = '%s/%s' % (resDir, fn)
        if not os.path.isfile(fn):
            continue

        trackNumber, bFn = [int(s) for s in re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", fn)]

        if os.path.isfile(fn) and fn[-2:] == '.p' and bFn == bChoice:
            f = cPickle.load(open(fn, 'rb'))
            timeList, percList, b = f['t'], f['p'], f['b']
            if len(timeList) != len(percList):
                continue

            # Fit curve to percList
            popt, pcov = curve_fit(fitFunc, timeList, np.array(percList)[:,0])
            perr = np.sqrt( np.diag(pcov) )

            # Get plot data for fits
            timeListFit = timeList
            percListFit = fitFunc(timeListFit, *popt)

            # Get last entries of perc lists for the 
            # recombination results
            # resList.append( percList[-1] )
            resList.append( (popt[0], perr[0]) )

            color = cr.getColor(N, i)
            # cr.plotFigure(fig, ax, timeList, np.array(percList)*100, xlabel=r'time [ns]', ylabel=r'Recombination [%]', label=None, title='', color=color)
            # cr.plotFigure(fig, ax, timeListFit, np.array(percListFit)*100, xlabel=r'time [ns]', ylabel=r'Recombination [%]', label=None, title='', color=color)

            # Increment loop counter
            i += 1

    if i == 0:
        return (0, 0), (0, 0)

    # Plot the results for all trajectories
    # Transform to percent
    resList = np.array(resList) * 100

    figRes, axRes = plt.subplots()
    mu, sigma = KDE(figRes, axRes, resList, smooth=5, xlabel=r'Recombination [%]', ylabel=r'Probability', label=None, title=r'', showDetail=True, fit=True, cauchy=True)

    mean = np.average(resList[:,0], weights=1./np.square(resList[:,1]))
    std = np.sqrt(1./np.sum(1./np.square(resList[:,1])))

    print 'Weighted Mean recombination: %f +/- %f' % (mean, std)
    print 'Mean recombination: %f +/- %f' % (np.mean(resList[:,0]), np.std(resList[:,0]))
    print 'mu: %f +/- %f; sigma: %f +/- %f' % (mu[0], mu[1], sigma[0], sigma[1])
    raw_input('')
    return mu, sigma

# === GENERATE ===
def generate(b, A, D, vD, trajectory, dt=1.e-12, tt=10.e-9, gauss=True, interrupt=False, perp=False, out_q=None, show=False, save=False, seed=None):
    if seed:
        np.random.seed(seed)

    # Generate Xe+ positions and create KD tree from them
    print 'Generate Xe+ distribution...',
    XePosList = initialDistributionTrack(b, trajectory, gauss=gauss, perp=perp)
    xyz = np.c_[np.array(XePosList)[:,0], np.array(XePosList)[:,1], np.array(XePosList)[:,2]]
    XeTree = spatial.cKDTree(xyz)
    print 'Done!'

    # Generate e positions and check for intersection
    print 'Generate e distribution...',
    ePosList = initialDistributionTrack(b, trajectory, gauss=gauss, perp=perp, tree=XeTree, A=A)
    xyz = np.c_[np.array(ePosList)[:,0], np.array(ePosList)[:,1], np.array(ePosList)[:,2]]
    eTree = spatial.cKDTree(xyz)
    print 'Done!'
    
    # Store information for plots
    timeList = [0]
    NList = [len(ePosList)]
    N = NList[0]
    percList = [0.]

    if show:
        # Initialize plots
        fig, h, hBottom = plotInit(XePosList, timeList, NList, tt)

    if save:
        posDict = {}
        posDict['Xe'] = {0: XePosList}
        posDict['e'] = {0: ePosList}

    # == Time Loop ==
    print '=== TIME LOOP START ==='
    loopStart = time.clock()
    loop = list(enumerate(np.arange(dt, tt, dt)))
    for i, t in loop:
        # Show progress
        sys.stdout.write('\r')
        perc = int(float(i)/len(loop)*100)
        sys.stdout.write("[%-20s] %d%% (Ne = %d)" % ('='*int(0.2*perc), perc, len(ePosList)))
        sys.stdout.flush()

        # Loop over all electrons
        ePosListNew = []
        XeRemoveIdx = []
        for ePos in ePosList:
            # Drift electron
            ePosNew = cr.driftElectron(ePos, D, vD, dt, eTree, XeTree, ePosList, XePosList)

            # Check for intersection with Xe
            while True:
                intersect, idxList = cr.checkIntersect(ePosNew, XeTree, A, None, ePos, XePosList)
                if intersect:
                    for idx in idxList:
                        if idx not in XeRemoveIdx:
                            # Electron hit Xe+
                            XeRemoveIdx.append( idx )
                            break
                        else:
                            continue
                    else:
                        ePosListNew.append( ePosNew )
                        break
                    break
                else:
                    # Electron survives
                    ePosListNew.append( ePosNew )
                    break
            
        ePosList = ePosListNew
        # Renew the XeTree
        XePosList = list([XePos for u, XePos in enumerate(XePosList) if u not in XeRemoveIdx])

        if save:
            # Store Xe+ & e positions
            posDict['Xe'][t] = XePosList
            posDict['e'][t] = ePosList

        # Works faster than expected
        xyz = np.c_[np.array(XePosList)[:,0], np.array(XePosList)[:,1], np.array(XePosList)[:,2]]
        XeTree = spatial.cKDTree(xyz)
        xyz = np.c_[np.array(ePosList)[:,0], np.array(ePosList)[:,1], np.array(ePosList)[:,2]]
        eTree = spatial.cKDTree(xyz)

        # Number of electrons vs. time
        timeList.append(t*1.e9)
        NList.append(len(ePosList))
        p = 1 - float(len(ePosList)) / N
        percList.append( p )

        # Test on every 10th iteration if change in 
        # recombination is small. If so, interrupt
        if interrupt and (not i % 10):
            if len(NList) > interrupt:
                mean = .5*(NList[-interrupt] + NList[-2])
                if mean == NList[-1]:
                    print
                    print 'In region of small changes -> interrupting'
                    break

        if show:
            updatePlot(fig, h, hBottom, t, i, timeList, NList, ePosList, XePosList, show=show)

    # After the time loop is done, count remaining particles 
    print
    print 'Calculation time:', time.clock() - loopStart 
    print 'Electrons remaining:', len(ePosList)
    print 'Recombination percentage:', 1 - float(len(ePosList))/N

    if save:
        # Save Xe+ & e positions to disk
        cPickle.dump(posDict, open(pickleOut, 'wb'))

    if out_q:
        print
        print 'Reached end'
        outDict = {'perp': perp, 't': timeList, 'p': percList}
        out_q.put( outDict )

    return timeList, percList

def initialDistributionTrack(b, trajectory, gauss=True, perp=False, tree=None, A=None):
    # Energy needed to create an electron-ion pair
    W = 15.6    # eV

    if tree and not A:
        print 'Please pecify the recombination radius when using a tree!'
        return False

    ranPosList = []
    for i in range(len(trajectory) - 1):
        p1, p2 = trajectory[i], trajectory[i+1]
        if perp:
            # Switch x and z positions
            p1 = [-p1[2], p1[1], p1[0], p1[3]]
            p2 = [-p2[2], p2[1], p2[0], p2[3]]
            
        pos1, E1, pos2, E2 = np.array(p1[:3]), p1[3], np.array(p2[:3]), p2[3]

        deltaE = E1 - E2
        # There are processes in which no energy is transfered
        if not deltaE:
            continue

        # Number of particles created
        deltaN = int( deltaE / W )

        vec = pos2 - pos1
        vX, vY, vZ = vec
        # Vectors for the random variables
        rZ = vec / np.linalg.norm( vec )

        if vX == 0 and vZ == 0:
            rX = np.array([0, 0, 1])
            rY = np.array([0, 1, 0])
        elif vY == 0 and vZ == 0:
            rX = np.array([1, 0, 0])
            rY = np.array([0, 0, 1])
        else:
            rX = vZ / (vX**2 + vZ**2) * np.array([vZ, 0, vX])
            rY = vZ / (vY**2 + vZ**2) * np.array([0, vZ, vY])

        length = np.linalg.norm( vec )
        for n in range(deltaN): 
            ranX, ranY, ranZ = getParticleTrack(b, length, gauss)
            x, y, z = tuple(pos1 + ranX*rX + ranY*rY + ranZ*rZ)
            if np.isnan(x):
                print vX, vY, vZ
                print rX, rY, rZ
                continue

            '''
            if tree:
                intersect, idx = cr.checkIntersect((x, y, z), tree, A)
                while intersect:
                    ranX, ranY, ranZ = getParticleTrack(b, length, gauss)
                    x, y, z = tuple(pos1 + ranX*rX + ranY*rY + ranZ*rZ)
                    intersect, idx = cr.checkIntersect((x, y, z), tree, A)
                print n
            '''

            ranPosList.append( (x, y, z) )

    return ranPosList

def getParticleTrack(b, length, gauss=True):
    # Get random height
    ranZ = np.random.uniform(0, np.linalg.norm( length ))
    # Get random radius
    if gauss:
        ranX, ranY = np.random.normal(0, b, 2)
    else:
        ranX, ranY = np.random.uniform(0, b, 2)

    return ranX, ranY, ranZ

def getTrajectories(fn):
    trajectoryList = []

    with open(fn) as f:
        trajectory = []
        while True:
            line = f.readline()

            if line.lstrip().startswith('#'):
                continue
            elif line.startswith('00000'):
                # Use only tracks which deposited their whole energy
                if trajectory and trajectory[-1][3] == 0.:
                    trajectoryList.append( trajectory )
                line = f.readline()
                if line.strip()=='':
                    break

                # Create empty list to store trajectory data in
                while not line.startswith('11111'):
                    line = f.readline()
                    continue
                trajectory = []
                continue
            data = line.strip().split()
            if int(data[-2]) == 2:
                continue

            data = [float(entry) for entry in data[0:4]]
            # Transform cm to m
            data = [data[i]*1.e-2 for i in range(3)] + [data[3]]
            trajectory.append( data )

    return trajectoryList

# Load trajectories from multiple files,
# shuffle them and return N of them
def getTrajectoriesRandom(N=100, multi=False):
    if multi:
        # Load all trajectories in the folder 
        # trajectoryPath = 'cr/trajectory_files/'
        trajectoryList = []
        for fn in os.listdir(trajectoryPath):
            if os.path.isfile(trajectoryPath + fn):
                trajectoryList += list( getTrajectories(trajectoryPath + fn) )
    else:
        trajectoryList = getTrajectories('cr/pe-trajectories_Co60.dat')

    random.shuffle(trajectoryList)
    trajectoryOut = []

    if not N:
        N = len(trajectoryList)

    for i in range(N):
        t = trajectoryList[i]
        # R = getRandomRotation()
        R = rand_rotation_matrix()
        t = [list(R.dot(x[:3])) + [x[3]] for x in t]
        trajectoryOut.append( t )

    return trajectoryOut

# === PLOTS ===
def plotInit(XePosList, timeList, NList, tt):
    fig = plt.figure(figsize=(5, 6))

    # Plot area
    left, width = 0., 0.9
    bottom, height = 0.27, 0.7

    # Main
    rect_main = [left, bottom, width, height]
    ax = plt.axes(rect_main, projection='3d')

    # Bottom
    rect_bottom = [0.15, 0.1, 0.8, 0.13]
    axBottom = plt.axes(rect_bottom)

    # Draw figure
    fig.canvas.draw()
    
    # Scatter plot for Xe+
    scale = 1000    # m -> mm
    XePosList = np.array(XePosList) * scale
    h = ax.scatter(np.array(XePosList)[:,0], np.array(XePosList)[:,1], np.array(XePosList)[:,2], marker='.', c='r')

    # Number of electrons vs. time
    hBottom = axBottom.plot(timeList, NList)

    # Set axis labels
    ax.set_xlabel('x [mm]'), ax.set_ylabel('y [mm]'), ax.set_zlabel('z [mm]')
    axBottom.set_xlabel(r't [ns]')
    axBottom.set_ylabel(r'N')

    # Set axis limits
    axBottom.set_xlim(0., tt*1.e9) # in [ns]
    axBottom.set_yscale('symlog', nonposx='clip')
    axBottom.set_ylim(0, NList[0]*1.1) 

    return fig, h, hBottom

def updatePlot(fig, h, hBottom, t, i, timeList, NList, ePosList, XePosList, show=False):
    # Main
    posPlot = np.array([list(x) for xs in zip(XePosList, ePosList) for x in xs])
    scale = 1000 # m -> mm
    posPlot = posPlot*scale
    h._facecolor3d = [[1, 0, 0, 1], [0, 0, 1, 1]] * int(len(posPlot)*.5)
    h._edgecolor3d = h._facecolor3d
    h._offsets3d = (posPlot[:,0], posPlot[:,1], posPlot[:,2])

    # Bottom
    hBottom[0].set_data(timeList, NList)

    # Figure title
    fig.suptitle('t = %.3f ns' % (t*1.e9))

    # Draw figure
    fig.canvas.draw()
    fig.canvas.flush_events()
    fig.savefig('%s/step_%d.png' % (gifOut, i), format='png', dpi=100)
            
    plt.draw()
    if show:
        fig.show()

def plotHistogram(fig, ax, data, alpha=1, xlabel=None, ylabel=None, label=None, title=None):
    data = np.array( data )
    if len(data.shape) == 2 and x.shape[1] == 2:
        data = data[:,0]

    weights = np.ones_like(data)/float(len(data))
    ax.hist(data, 10, normed=0, weights=weights, alpha=alpha, label=label)

    # Set axes label
    ax.set_xlabel( xlabel )
    ax.set_ylabel( ylabel )

    ax.legend(loc='best')

    if title:
        fig.suptitle( title )

    fig.canvas.draw()
    fig.show()

def plotTrajectories(trajectories, plotRange, proj='xy'):
    fig, ax = plt.subplots()
    tList = trajectories[plotRange[0]:plotRange[1]+1]

    for i, t in enumerate(tList):
        x, y, z, E = zip(*t)
        x, y, z = np.array(x)*1.e3, np.array(y)*1.e3, np.array(z)*1.e3
        if proj == 'xy' or proj == 'yx':
            ax.plot(x, y, c=cr.getColor(len(tList), i))
            ax.set_xlabel('x [mm]'), ax.set_ylabel('y [mm]')
        elif proj == 'xz' or proj == 'zx':
            ax.plot(x, z, c=cr.getColor(len(tList), i))
            ax.set_xlabel('x [mm]'), ax.set_ylabel('z [mm]')
        elif proj == 'yz' or proj == 'zy':
            ax.plot(y, z, c=cr.getColor(len(tList), i))
            ax.set_xlabel('y [mm]'), ax.set_ylabel('z [mm]')
        else:
            print 'Wrong projection chosen!'
            return

    fig.show()
    return

def plotSingleResults(bList, muList, sigmaList):
    print bList
    print muList
    print sigmaList

    muList = np.array( muList )
    FWHMList = 2*np.sqrt(2*np.log(2)) * np.array( sigmaList )
    res = np.array( [FWHMList[i][0]/muList[i][0] for i in range(len(muList))] ) * 100
    res_err = np.array( [np.sqrt((FWHMList[i][1]/muList[i][0])**2 + (FWHMList[i][0]/muList[i][0]**2 * muList[i][1])**2) for i in range(len(muList))] ) * 100

    # Mu on first y-axis
    fig, axMu = plt.subplots()

    mu, mu_err = zip(*muList)
    axMu.errorbar(bList, mu, yerr=mu_err, capsize=3, markeredgewidth=1, label='Recombination', c='b')

    axMu.set_xlabel(r'b [nm]')
    axMu.set_ylabel(r'Recombination [%]', color='b')
    axMu.grid()

    # FWHM on second y-axis
    axFWHM = axMu.twinx()

    FWHM, FWHM_err = zip(*FWHMList)
    axFWHM.errorbar(bList, FWHM, yerr=FWHM_err, capsize=3, markeredgewidth=1, c='orange', label='FWHM')
    
    axFWHM.set_ylabel(r'FWHM [%]', color='orange')
    axFWHM.grid(ls='--', lw=.4)

    # Second figure for resolution
    figRes, axRes = plt.subplots()
    axRes.errorbar(bList, res, yerr=res_err, capsize=5, markeredgewidth=2)
    axRes.set_xlabel(r'b [nm]')
    axRes.set_ylabel(r'Recombination resolution [%]')

    # fig.tight_layout()
    fig.show(), figRes.show()
    raw_input('')

# === SUPPORT ===
def fitFunc(x, A, a, b):
    return A / (1 + np.exp(a*x**b))

def KDE(fig, ax, x, smooth=1, xlabel=None, ylabel=None, label=None, title=None, showDetail=False, fit=False, cauchy=False):
    x = np.array( x )
    N = len( x )
    if len(x.shape) != 2 and x.shape[1] != 2:
        print 'Need mean and std for each element in x!'
        return False

    # Each std is multiplied with smooth in order 
    # to get better visualization
    x = np.array( [(i[0], i[1]*smooth) for i in x])

    if showDetail:
        # Create normal distribution for each data point
        # and plot it
        for xi in x:
            m, s = xi

            plotRange = np.linspace(m-3*s, m+3*s, 1000)
            ax.plot(plotRange, 1./N * normal(plotRange, m, s), color='r', ls='--')

    # Create sum of all normal distributions and plot it
    rangeMin = min([i[0] - 3*i[1] for i in x])
    rangeMax = max([i[0] + 3*i[1] for i in x])

    if not showDetail:
        color = None
    else:
        color = 'b'
    plotRange = np.linspace(rangeMin, rangeMax, 10000)
    KDEfunction = functionSum(plotRange, normal, x, N)
    ax.plot(plotRange, KDEfunction, color=color, label=label)
    if fit:
        if cauchy:
            p0 = [80, 8, 1.]
            popt, pcov = curve_fit(cauchyDist, plotRange, KDEfunction, p0=p0)
            ax.plot(plotRange, cauchyDist(plotRange, *popt), color=color, label='')
        else:
            p0 = [80, 8]
            popt, pcov = curve_fit(normal, plotRange, KDEfunction, p0=p0)
            ax.plot(plotRange, normal(plotRange, *popt), color=color, label='')

        perr = np.sqrt( np.diag(pcov) )
        print popt, perr

        mu, sigma = (popt[0], perr[0]), (abs(popt[1]), abs(perr[1]))
        # Use: FWHM = 2*gamma = 2*sqrt(2*ln(2))*sigma
        if cauchy:
            sigma = [np.sqrt(2 * np.log(2))*float(s) for s in sigma]

    ax.legend(loc='best')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if title:
        fig.suptitle(title)

    fig.canvas.draw()
    fig.show()

    if fit:
        return mu, sigma

def functionSum(x, func, paramList, N):
    s = 0
    for i in range(N):
        s += 1./N * func(x, *paramList[i])
    return s

def cauchyDist(x, mu, gamma, A):
    return A/(np.pi*gamma*(1+((x - mu)/gamma)**2))

def normal(x, mu, sigma):
    return 1./np.sqrt(2*np.pi*sigma**2) * np.exp(-(x - mu)**2 / (2*sigma**2) )

# === ROTATION MATRICES ===
# http://blog.lostinmyterminal.com/python/2015/05/12/random-rotation-matrix.html
def rand_rotation_matrix(deflection=1.0, randnums=None):
    """
    Creates a random rotation matrix.

    deflection: the magnitude of the rotation. For 0, no rotation; for 1, competely random
    rotation. Small deflection => small perturbation.
    randnums: 3 random numbers in the range [0, 1]. If `None`, they will be auto-generated.
    """
    # from http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c

    if randnums is None:
        randnums = np.random.uniform(size=(3,))

    theta, phi, z = randnums
    theta = theta * 2.0*deflection*np.pi  # Rotation about the pole (Z).
    phi = phi * 2.0*np.pi  # For direction of pole deflection.
    z = z * 2.0*deflection  # For magnitude of pole deflection.

    # Compute a vector V used for distributing points over the sphere
    # via the reflection I - V Transpose(V).  This formulation of V
    # will guarantee that if x[1] and x[2] are uniformly distributed,
    # the reflected points will be uniform on the sphere.  Note that V
    # has length sqrt(2) to eliminate the 2 in the Householder matrix.

    r = np.sqrt(z)
    Vx, Vy, Vz = V = (
        np.sin(phi) * r,
        np.cos(phi) * r,
        np.sqrt(2.0 - z)
        )

    st = np.sin(theta)
    ct = np.cos(theta)

    R = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))

    # Construct the rotation matrix  ( V Transpose(V) - I ) R.

    M = (np.outer(V, V) - np.eye(3)).dot(R)
    return M

# Wrong: This will lead to a higher density at the poles!
def getRandomRotation():
    randomAngles = np.random.uniform(0, 2*np.pi, 3)
    return np.dot(np.dot(rotationMatrix(randomAngles[0], 'x'), rotationMatrix(randomAngles[1], 'y')), rotationMatrix(randomAngles[2], 'z'))

def rotationMatrix(angle, coordinate='x'):
    c, s = np.cos(angle), np.sin(angle)
    if coordinate == 'x':
        R = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

    elif coordinate == 'y':
        R = np.array([[c, 0, -s], [0, 1, 0], [s, 0, c]])

    elif coordinate == 'z':
        R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

    else: 
        return False

    return R

def generateScripts(fName, fNameExec, Nl, Nh, lima=False, single=False, b=None):
    if lima:
        ppn = 24
    else:
        ppn = 4

    # Job script
    with open(fName, 'w') as f:
        f.write('#!/bin/bash -l\n')
        f.write('#\n')
        f.write('#PBS -V -j oe\n')
        f.write('#PBS -l nodes=1:ppn=%d,walltime=05:00:00\n' % ppn)
        f.write('#PBS -N cr_track_${NUMBER}\n')
        f.write('#PBS -o cr/logs/cr_track_${NUMBER}.log\n')
        f.write('#\n')
        f.write('cd $PBS_O_WORKDIR\n')
        if lima:
            f.write('source $VAULT/.exorc_test\n')
        f.write('python columnar_recombination_track.py -t ${NUMBER} ')
        if b:
            f.write('-b %d ' % b)
        if single:
            f.write('-s\n\n')
        else:
            f.write('-m\n\n')

    # Executes the job script in a loop
    with open(fNameExec, 'w') as f:
        f.write('#!/usr/bin/env bash\n')
        f.write('for MY_OPT in {%d..%d}; do\n' % (Nl, Nh))
        f.write('\tqsub -v NUMBER=$MY_OPT %s\n' % fName)
        f.write('done\n\n')

def get_args():
    ap = argparse.ArgumentParser(description=' ')

    ap.add_argument('-t', '--track', type=int, help='Track to process, beginning with 0', required=False)
    ap.add_argument('-b', '--parameter_b', type=float, help='Parameter b to use (in nm)', required=False)
    ap.add_argument('-m', '--multi', help='Process a trajectory multiple times', action='store_true')
    ap.add_argument('-cr', '--combineresults', help='Combine the results of the files within the result-directory', action='store_true')
    ap.add_argument('-gs', '--genscripts', help='Generate job scripts', action='store_true')
    ap.add_argument('-s', '--single', help='Process parallel tracks only', action='store_true')
    ap.add_argument('-rt', '--randomtrajectories', help='Pick trajectories randomly and rotate them with random angle in 3D space. Store result in cPickle-file', action='store_true')

    args = ap.parse_args()
    return args.track, args.multi, args.combineresults, args.genscripts, args.single, args.randomtrajectories, args.parameter_b

if __name__ == '__main__':
    main()

