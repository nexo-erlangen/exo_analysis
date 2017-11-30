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
from uncertainties import ufloat, unumpy

try:
    import seaborn as sns
    sns.set_style('whitegrid', {'axes.grid' : False})
    sns.set(style='ticks')
except:
    pass

import columnar_recombination as cr

# Doesn't work when run on cluster
try:
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    # plt.style.use('seaborn-dark-palette')
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
    SHOW = False
    SAVE = True
    PERPENDICULAR = True
    REPEAT = 1

    # Constants
    if parameter_b:
        b = parameter_b*1.e-9
    else:
        # STD of gaussian distribution around columnar axis
        # OR starting distance of electron to Xe+ ion
        b = 10.e-9 # 50.e-9

    # Use 3*sigma interval
    if not GAUSS:
        b *= 3
        
    A = 5.e-9 # 49.e-9
    # b = 10.e-9
    vD = 1705
    Dt = 55.e-4 # 30.e-4 # 55.e-4
    Dl = 0.1 * Dt # 0.12 * Dt
    D = [Dt, Dt, Dl]

    # Time
    dt = 3.e-16 # 1.e-14 # 1./(2 * (2*55.e-4 + 5.5e-4)) * (4.e-9 / 3)**2 
    tt = 3333 * dt # 2.e-12 # 1000 * dt  # 2 * b / 1705

    if genscripts:
        jobs = (0, 50)
        generateScripts('cr_job.sh', 'cr_exec_job.sh', jobs[0], jobs[1], lima=True, single=single, b=int(b*1.e9))
        return

    if combineresults:
        if single:
            bList = [5, 7, 10, 12, 15] # , 300]
            bTitleList = [5, 7.5, 10, 12.5, 15] # , 300]
            muList, sigmaList = [], []
            figKDE, axKDE = plt.subplots(figsize=(7,3))
            figKDE.subplots_adjust(bottom=.2)
            for i, bChoice in enumerate(bList): 
                mu, sigma, res = combineResultsSingle('cr/track_results_single_uniform', bChoice=bChoice, show=False)
                muList.append( mu ), sigmaList.append( sigma )
                if len(res) > 0:
                    KDE(figKDE, axKDE, res, smooth=20, xlabel=r'Recombination fraction [\%]', ylabel=r'Probability', label=r'b = %.1f nm' % bTitleList[i], title=None, showDetail=False, color=cr.getColor(len(bList), i), fit=True, cauchy=False, show=False)
            axKDE.set_xlim(left=35, right=100)
            figKDE.show()
            raw_input('')

            plotSingleResults(bTitleList, muList, sigmaList)

        else:
            combineResults('cr/track_results_nodiff')
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
        # trajectories = getTrajectories('cr/pe-trajectories.dat')
        # plotTrajectories(trajectories, [0, 100], proj='xy')
        # plotTrajectories(trajectories, [0, 100], proj='xz')
        # plotTrajectories(trajectories, [0, 100], proj='yz')
        # raw_input('')
        # return

        trajectory = trajectories[trackNumber]
        # TODO: Remove workaround!
        print b
        getResultsSingle(b, A, D, vD, trajectory, dt, tt, gauss=GAUSS, interrupt=INTERRUPT, repeat=REPEAT, resultOut='cr/track_results_single_uniform/track_%d_b_%dnm.p' % (trackNumber, int(b*1.e9)))
        return 

    # === SINGLE ===
    # Need trajectories for the following functions
    trajectories = getTrajectories('cr/pe-trajectories.dat')
    # getMeanLET(trajectories, show=SHOW)
    # raw_input('')
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
        getResults(b, A, D, vD, trajectory, dt, tt, gauss=GAUSS, interrupt=INTERRUPT, repeat=REPEAT, resultOut='cr/track_results_nodiff/track_%d.p' % trackNumber)
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
    # if repeat > 1: 
    percListPerp = [(np.mean(p), np.std(p)) for p in zip(*percListListPerp)]
    percListPar = [(np.mean(p), np.std(p)) for p in zip(*percListListPar)]
    outDict = {'t': timeList, 'pPerp': percListPerp, 'pPar': percListPar}

    # else: 
    #    outDict = {'t': timeList, 'pPerp': percListListPerp[0], 'pPar': percListListPar[0]}

    print resultOut
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
            # cr.plotFigure(figPerp, axPerp, timeList, np.array(percListPerp)*100, xlabel=r'time [ps]', ylabel=r'Recombination [%]', label=None, title='Perpendicular', color=color)
            # cr.plotFigure(figPerp, axPerp, timeListFit, np.array(percListPerpFit)*100, xlabel=r'time [ps]', ylabel=r'Recombination [%]', label=None, title='Perpendicular', color=color)
            # Parallel
            # cr.plotFigure(figPar, axPar, timeList, np.array(percListPar)*100, xlabel=r'time [ps]', ylabel=r'Recombination [%]', label=None, title='Parallel', color=color)
            # cr.plotFigure(figPar, axPar, timeListFit, np.array(percListParFit)*100, xlabel=r'time [ps]', ylabel=r'Recombination [%]', label=None, title='Parallel', color=color)

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

    figResPerp, axResPerp = plt.subplots(figsize=(7,3))
    figResPerp.subplots_adjust(bottom=.2)
    KDE(figResPerp, axResPerp, resListPerp, smooth=20, xlabel=r'Recombination fraction [\%]', ylabel=r'Probability', label=None, title=r'Perpendicular', showDetail=True, fit=True)

    figResPar, axResPar = plt.subplots(figsize=(7,3))
    figResPar.subplots_adjust(bottom=.2)
    KDE(figResPar, axResPar, resListPar, smooth=20, xlabel=r'Recombination fraction [\%]', ylabel=r'Probability', label=None, title=r'Parallel', showDetail=True, fit=True)
    raw_input('')

    figRes, axRes = plt.subplots(figsize=(7, 3))
    figRes.subplots_adjust(bottom=.2)
    KDE(figRes, axRes, resListPerp, smooth=20, xlabel=r'Recombination fraction [\%]', ylabel=r'Probability', label='Perpendicular', title=None, showDetail=False)
    KDE(figRes, axRes, resListPar, smooth=20, xlabel=r'Recombination fraction [\%]', ylabel=r'Probability', label='Parallel', title=None, showDetail=False)

    figHist, axHist = plt.subplots()
    plotHistogram(figHist, axHist, np.array(resListPerp)[:,0], xlabel=r'Recombination fraction [\%]', ylabel='Probability', label='Perpendicular', title=None)
    plotHistogram(figHist, axHist, np.array(resListPar)[:,0], alpha=.7, xlabel=r'Recombination fraction [\%]', ylabel='Probability', label='Parallel', title=None)

    meanPerp = np.average(resListPerp[:,0], weights=1./np.square(resListPerp[:,1])) # np.mean( resListPerp[:,0] )
    stdPerp = np.sqrt(1./np.sum(1./np.square(resListPerp[:,1]))) # np.sum( np.square(resListPerp[:,1]) / len(resListPerp) )
    meanPar = np.average(resListPar[:,0], weights=1./np.square(resListPar[:,1])) # np.mean( resListPar[:,0] )
    stdPar = np.sqrt(1./np.sum(1./np.square(resListPar[:,1]))) # np.sum( np.square(resListPar[:,1]) / len(resListPar) )

    print 'Mean recombination: %f +/- %f (perpendicular), %f +/ - %f (parallel)' % (meanPerp, stdPerp, meanPar, stdPar)
    print np.mean(resListPerp[:,0]), np.mean(resListPar[:,0])
    print np.std(resListPerp[:,0]), np.std(resListPar[:,0])
    mPerp, mPar = ufloat(meanPerp, stdPerp), ufloat(meanPar, stdPar)
    mPerpFit, mParFit = ufloat(np.mean(resListPerp[:,0]), np.std(resListPerp[:,0])), ufloat(np.mean(resListPar[:,0]), np.std(resListPar[:,0]))
    mRes = (mPerp - mPar) / (mPerp + mPar) * 100
    mResFit = (mPerpFit - mParFit) / (mPerpFit + mParFit) * 100
    print mPerp, mPar, mRes
    print mPerpFit, mParFit, mResFit
    print
    raw_input('')

def combineResultsSingle(resDir, bChoice=50, show=True):
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
            # print percList
            try:
                p0 = [ 0.45578102, -6.3189622, 0.79577747]
                popt, pcov = curve_fit(fitFunc, timeList, np.array(percList), p0=p0)
            except:
                continue

            print popt, pcov
            perr = np.sqrt( np.diag(pcov) )

            # Get plot data for fits
            timeListFit = timeList
            percListFit = fitFunc(timeListFit, *popt)

            # Get last entries of perc lists for the 
            # recombination results
            # resList.append( percList[-1] )
            resList.append( (popt[0], perr[0]) )

            color = cr.getColor(N, i)
            # cr.plotFigure(fig, ax, timeList, np.array(percList)*100, xlabel=r'time [ps]', ylabel=r'Recombination [\%]', label=None, title='', color=color)
            # cr.plotFigure(fig, ax, timeListFit, np.array(percListFit)*100, xlabel=r'time [ps]', ylabel=r'Recombination [\%]', label=None, title='', color=color)

            # Increment loop counter
            i += 1

    if i == 0:
        return (0, 0), (0, 0), []

    # Plot the results for all trajectories
    # Transform to percent
    resList = np.array(resList) * 100

    figRes, axRes = plt.subplots()
    mu, sigma = KDE(figRes, axRes, resList, smooth=20, xlabel=r'Recombination fraction [%]', ylabel=r'Probability', label=None, title=r'', showDetail=True, fit=True, cauchy=True, show=show)

    mean = np.average(resList[:,0], weights=1./np.square(resList[:,1]))
    std = np.sqrt(1./np.sum(1./np.square(resList[:,1])))

    print 'Weighted Mean recombination: %f +/- %f' % (mean, std)
    print 'Mean recombination: %f +/- %f' % (np.mean(resList[:,0]), np.std(resList[:,0]))
    print 'mu: %f +/- %f; sigma: %f +/- %f' % (mu[0], mu[1], sigma[0], sigma[1])
    if show:
        raw_input('')
    return mu, sigma, resList

# === GENERATE ===
def generate(b, A, D, vD, trajectory, dt=1.e-12, tt=10.e-9, gauss=True, interrupt=False, perp=False, out_q=None, show=False, save=False, seed=None):
    if seed:
        np.random.seed(seed)

    # Generate Xe+ positions and create KD tree from them
    print 'Generate Xe+ distribution...',
    XePosList, ePosList = initialDistributionTrackStraight(trajectory, perp=perp, dist=b)
    # XePosList = initialDistributionTrack(b, trajectory, gauss=gauss, perp=perp)
    xyz = np.c_[np.array(XePosList)[:,0], np.array(XePosList)[:,1], np.array(XePosList)[:,2]]
    XeTree = spatial.cKDTree(xyz)
    print 'Done!'

    meanDist = cr.getMeanDistance(XePosList, XeTree)
    print meanDist
    print 'Mean distance of Xe+ ions: %f nm' % (meanDist*1.e9)
    # Compare mean distance to 5*sigma value of maximum diffusion
    if 1./max(D) * (5*meanDist)**2 < dt:
        print 'WARNING: chosen time step not small enough!'

    # Generate e positions and check for intersection
    print 'Generate e distribution...',
    # ePosList = initialDistributionTrack(b, trajectory, gauss=gauss, perp=perp, tree=XeTree, A=A)
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
    print 'Number of time steps: %d' % int(tt / dt)
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
        timeList.append(t*1.e12)
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
            updatePlot(fig, h, hBottom, t, i, timeList, np.array([float(n) for n in NList]) / N * 100, ePosList, XePosList, show=show)

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

def getMeanLET(trajectories, show=False):
    W = 15.6
    LETList = []
    lengthList = []
    EList, NList = [], []

    for trajectory in trajectories:
        for i in range(len(trajectory) - 1):
            p1, p2 = trajectory[i], trajectory[i+1]
            pos1, E1, pos2, E2 = np.array(p1[:3]), p1[3], np.array(p2[:3]), p2[3]
            vec = pos2 - pos1
            if np.linalg.norm(vec) == 0:
                continue

            deltaE = E1 - E2
            if not deltaE:
                continue
            deltaN = int( deltaE // W )
            if not deltaN:
                continue

            LET = deltaE * 1.e-6 / (np.linalg.norm(vec) * 1.e2) * 1./2.96
            if LET > 1000:
                continue

            LETList.append( LET )
            lengthList.append( np.linalg.norm(pos2 - pos1) )
            NList.append( deltaN )

    fig, ax = plt.subplots(figsize=(7,3))
    figL, axL = plt.subplots(figsize=(7,3))
    figN, axN = plt.subplots(figsize=(7,3))
    fig.subplots_adjust(bottom=.2)
    figL.subplots_adjust(bottom=.2)
    figN.subplots_adjust(bottom=.2)

    Nbins = 100
    Nmax = 200
    NList = np.array([float(N) for N in NList])

    # LET
    ax.hist(np.array(LETList), bins=np.linspace(0, 1000, Nbins), weights=None, normed=True)
    ax.set_xlabel(r'LET [MEV cm$^2$/g]')
    ax.set_ylabel(r'Normalised counts / (%.1f MEV cm$^2$/g)' % (float(Nmax)/Nbins))
    ax.set_xlim(0, 1000)
    ax.set_yscale('log', nonposx='clip')

    meanLET = np.mean(LETList) # (sum(np.array(LETList) * np.array(NList)) / sum(NList) )
    stdLET = np.std(LETList)
    ax.axvline(x=meanLET, c='k', ls='--')
    sns.despine(ax=ax, fig=fig)
    fig.show()

    # Length
    lengthList = np.array( lengthList ) * 1.e6
    axL.hist(lengthList, bins=np.linspace(0, 150, Nbins), normed=True, weights=None)
    axL.set_xlabel(r'Track length [$\mu$m]')
    axL.set_ylabel(r'Normalised counts / (%.2f $\mu$m)' % (1./Nbins))
    axL.set_xlim(0, 150)
    axL.set_yscale('log', nonposx='clip')
    meanLength = np.mean(lengthList) # ( sum(np.array(lengthList) * np.array(NList)) / sum(NList) )
    stdLength = np.std(lengthList)
    axL.axvline(x=meanLength, c='k', ls='--')
    sns.despine(ax=axL, fig=figL)
    figL.show()

    # Particles
    NList = np.array( NList )
    axN.hist(NList, normed=True, bins=np.linspace(0, 2000, Nbins))
    axN.set_xlabel(r'Electron-ion pairs')
    axN.set_ylabel(r'Normalised counts / (%.2f Pairs)' % (1./Nbins))
    axN.set_xlim(0, 2000)
    axN.set_yscale('log', nonposx='clip')
    
    meanN, stdN = np.mean( NList ), np.std( NList )
    axN.axvline(x=meanN, c='k', ls='--')
    sns.despine(ax=axN, fig=figN)
    figN.show()

    from uncertainties import ufloat
    print ufloat(meanLET, stdLET)
    print ufloat(meanLength, stdLength)
    print ufloat(meanN, stdN)
    print
    print 'Mean LET: (%f +/- %f) MeV cm^2 / g' % (meanLET, stdLET)
    print 'Mean length: (%f +/- %f) mum' % (meanLength, stdLength)
    print 'Mean N: (%f +/- %f)' % (meanN, stdN)

# Create Xe+ ions along a straight track for each
# trajectory segment. Additionally, create electron
# in specified distance and arbitrary angle to the ion. 
# If no distance is specified, only Xe+ ions are created.
def initialDistributionTrackStraight(trajectory, perp=False, dist=None):
    W = 15.6    # eV

    # Particle positions
    XePosList, ePosList = [], []

    for i in range(len(trajectory) - 1):
        p1, p2 = trajectory[i], trajectory[i+1]
        if perp:
            # Switch x and z positions
            p1 = [-p1[2], p1[1], p1[0], p1[3]]
            p2 = [-p2[2], p2[1], p2[0], p2[3]]
            
        pos1, E1, pos2, E2 = np.array(p1[:3]), p1[3], np.array(p2[:3]), p2[3]
        vec = pos2 - pos1

        deltaE = E1 - E2
        # There are processes in which no energy is transfered
        if not deltaE:
            continue

        # Number of electron ion pairs created
        deltaN = int( deltaE // W )

        # Normalized vector vec
        rZ = vec / np.linalg.norm( vec )

        for n in range( deltaN ):
            # Random length for rZ (uniform distribution)
            ranZ = np.random.uniform(0, np.linalg.norm( vec ))

            # Xe+ position
            XePos = tuple(pos1 + ranZ * rZ)
            XePosList.append( XePos )

            if dist:
                # https://math.stackexchange.com/questions/1585975/how-to-generate-random-points-on-a-sphere
                # http://mathworld.wolfram.com/SpherePointPicking.html
                z = np.random.uniform(-1., 1.)
                theta = np.random.uniform(0, 2*np.pi)
                x, y = np.sqrt(1 - z**2)*np.cos(theta), np.sqrt(1 - z**2)*np.sin(theta)
                pos = np.array([x, y, z]) * dist

                ePos = tuple( np.array(XePos) + pos )
                ePosList.append( ePos )

    return XePosList, ePosList
                
def initialDistributionTrack(b, trajectory, gauss=True, perp=False, tree=None, A=None):
    from operator import itemgetter

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

        # https://stackoverflow.com/questions/19337314/generate-random-point-on-a-2d-disk-in-3d-space-given-normal-vector
        # vec is the normal vector of the area to create points in
        # -> Need vectors which span this area
        vec = pos2 - pos1
        vX, vY, vZ = vec

        # Vectors for the random variables
        rZ = vec / np.linalg.norm( vec )

        # Find minimum value of vec and get index
        # idxMin = list(vec).index(min(abs(vec)))
        idxMin = min(enumerate(abs(vec)), key=itemgetter(1))[0] 
        A = vec[idxMin]

        # Remove index from idxList to get
        # remaining components
        idxList = range(3)
        idxList.remove(idxMin)
        B, C = [vec[j] for j in idxList]

        # Get orthogonal vectors 
        u = np.array([0, -C, B])
        v = np.array([B**2 + C**2, -A*B, -A*C])

        # Make them orthonormal
        u = u / np.linalg.norm( u )
        v = v / np.linalg.norm( v )

        # Rotate vectors by the minimum index
        u = rotate(u, idxMin)
        v = rotate(v, idxMin)

        # u & v are now perpendicular to the normal
        # vector vec!

        '''
        if vX == 0 and vZ == 0:
            rX = np.array([0, 0, 1])
            rY = np.array([0, 1, 0])
        elif vY == 0 and vZ == 0:
            rX = np.array([1, 0, 0])
            rY = np.array([0, 0, 1])
        else:
            rX = vZ / (vX**2 + vZ**2) * np.array([vZ, 0, vX])
            rY = vZ / (vY**2 + vZ**2) * np.array([0, vZ, vY])
        '''

        length = np.linalg.norm( vec )
        for n in range(deltaN): 
            ranX, ranY, ranZ = getParticleTrack(b, length, gauss)
            x, y, z = tuple(pos1 + ranX*u + ranY*v + ranZ*rZ)
            if np.isnan(x):
                print vX, vY, vZ
                print u, v, rZ
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

def rotate(l, n):
    return np.array( list(l)[-n:] + list(l)[:-n] )

def getParticleTrack(b, length, gauss=True):
    # Get random height
    ranZ = np.random.uniform(0, np.linalg.norm( length ))
    # Get random radius
    if gauss:
        ranX, ranY = np.random.normal(0, b, 2)
    else:
        r = b + 1
        while r > b:
            ranX, ranY = np.random.uniform(0, b, 2)
            r = np.sqrt( ranX**2 + ranY**2 )

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
    h = ax.scatter(np.array(XePosList)[:,0], np.array(XePosList)[:,1], np.array(XePosList)[:,2], marker='.', c='crimson')

    # Number of electrons vs. time
    hBottom = axBottom.plot(timeList, [100])

    # Set axis labels
    ax.set_xlabel('x [mm]'), ax.set_ylabel('y [mm]'), ax.set_zlabel('z [mm]')
    axBottom.set_xlabel(r't [ps]')
    axBottom.set_ylabel(r'N [%]')

    # Set axis limits
    axBottom.set_xlim(0., tt*1.e12) # in [ps]
    # axBottom.set_yscale('symlog', nonposx='clip')
    # axBottom.set_ylim(0.2*NList[0], NList[0]*1.1) 
    axBottom.set_ylim(0, 100) 

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
    fig.suptitle('t = %.3f ps' % (t*1.e12))

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
    sns.despine(fig=fig, ax=ax, left=True, bottom=True)
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
    from uncertainties import ufloat, unumpy

    print bList
    print muList
    print sigmaList

    muList = np.array( muList )
    FWHMList = np.array( sigmaList ) # 2*np.sqrt(2*np.log(2)) * np.array( sigmaList )
    res = np.array( [FWHMList[i][0]/muList[i][0] for i in range(len(muList))] ) * 100
    res_err = np.array( [np.sqrt((FWHMList[i][1]/muList[i][0])**2 + (FWHMList[i][0]/muList[i][0]**2 * muList[i][1])**2) for i in range(len(muList))] ) * 100
    
    muListText = [r'$%s \pm %s$' % tuple(l) for l in [str(m).split('+/-') for m in unumpy.uarray(*zip(*muList))]]
    FWHMListText = [r'$%s \pm %s$' % tuple(l) for l in [str(m).split('+/-') for m in unumpy.uarray(*zip(*FWHMList))]]
    resListText = [r'$%s \pm %s$' % tuple(l) for l in [str(m).split('+/-') for m in unumpy.uarray(res, res_err)]]

    latexExport('cr/single_means.tex', [r'$b$ [nm]', r'$\mu_R [\%]$', r'$\sigma_R$ [\%]', r'$\frac{\sigma_R}{\mu_R}$'], zip(bList, muListText, FWHMListText, resListText), 'Mean values', 'tab:means')

    # Mu on first y-axis
    fig, axMu = plt.subplots(figsize=(7, 3))
    fig.subplots_adjust(bottom=.2)
    pal = sns.color_palette().as_hex()

    def expFit(x, a, b, c):
        x = np.array( x )
        return a * np.exp(-b * x) + c

    def expFit2(x, a, b, c):
        x = np.array( x )
        return a / (1 + b*np.exp(-c*x))

    mu, mu_err = zip(*muList)
    poptMu, pcovMu = curve_fit(expFit, bList, mu)
    print poptMu,
    perrMu = np.sqrt(np.diag(pcovMu))
    print 'Mu-Fit:'
    print '======='
    print 'f(b) = (%.3f+/-%.4f)*b + (%.3f+/-%.3f)\n' % (poptMu[0], perrMu[0], poptMu[1], perrMu[1])
    for k in range( len(poptMu) ):
        print ufloat(poptMu[k], perrMu[k]),
    print

    axMu.errorbar(bList, mu, yerr=mu_err, capsize=3, markeredgewidth=1, label='Recombination fraction', c=pal[0]) # , c='dodgerblue')
    bListPlot = np.linspace(min(bList), max(bList), 1000)
    axMu.plot(bListPlot, expFit(bListPlot, *poptMu), c=pal[0]) #, c='dodgerblue')

    axMu.set_xlabel(r'b [nm]')
    axMu.set_ylabel(r'Recombination fraction $\mu_R$ [\%]', color=pal[0], fontsize=12, fontweight='bold') # , color='dodgerblue')
    axMu.grid()

    # FWHM on second y-axis
    axFWHM = axMu.twinx()

    FWHM, FWHM_err = zip(*FWHMList)
    p0 = [90, 20, .5]
    poptFWHM, pcovFWHM = curve_fit(expFit2, bList, FWHM, p0=p0)
    perrFWHM = np.sqrt(np.diag(pcovFWHM))
    print 'FWHM-Fit:'
    print '========='
    print 'f(b) = (%.3f+/-%.4f)*b + (%.3f+/-%.3f)\n' % (poptFWHM[0], perrFWHM[0], poptFWHM[1], perrFWHM[1])
    for k in range( len(poptFWHM) ):
        print ufloat(poptFWHM[k], perrFWHM[k]),
    print

    axFWHM.errorbar(bList, FWHM, yerr=FWHM_err, capsize=3, markeredgewidth=1, label='FWHM', color=pal[1])
    axFWHM.plot(bListPlot, expFit2(bListPlot, *poptFWHM), c=pal[1])# , c='orange')
    
    axFWHM.set_ylabel(r'Standard deviation $\sigma_R$ [\%]', color=pal[1], fontsize=12, fontweight='bold') #, color='orange')
    # axFWHM.grid(ls='--', lw=.4)

    # Second figure for resolution
    figRes, axRes = plt.subplots(figsize=(7, 3))
    figRes.subplots_adjust(bottom=.2)
    poptRes, pcovRes = curve_fit(linear, bList, res)
    perrRes = np.sqrt(np.diag(pcovRes))
    print 'Resolution-Fit:'
    print '==============='
    print 'f(b) = (%.3f+/-%.4f)*b + (%.3f+/-%.3f)\n' % (poptRes[0], perrRes[0], poptRes[1], perrRes[1])

    axRes.errorbar(bList, res, yerr=res_err, capsize=3, markeredgewidth=1)
    axRes.plot(bList, linear(bList, *poptRes), color=pal[0]) # , color='dodgerblue')
    axRes.grid()
    axRes.set_xlabel(r'b [nm]')
    axRes.set_ylabel(r'Recombination resolution [\%]')
    axRes.set_yticks(np.linspace(axRes.get_yticks()[0],axRes.get_yticks()[-1],len(axFWHM.get_yticks())))
    sns.despine(fig=figRes, ax=axRes)

    # fig.tight_layout()
    sns.despine(ax=axMu, left=False, right=True)
    sns.despine(ax=axFWHM, left=True, right=False)
    fig.show(), figRes.show()
    raw_input('')

# === SUPPORT ===
def fitFunc(x, A, a, b):
    return A / (1 + np.exp(a*x**b))
    # return A/(1 + a*x**(-b))

def KDE(fig, ax, x, smooth=1, xlabel=None, ylabel=None, label=None, title=None, showDetail=False, color=None, fit=False, cauchy=False, show=True):
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
            ax.plot(plotRange, 20./N * normal(plotRange, m, s), color='crimson', ls='--', lw=.5)

    # Create sum of all normal distributions and plot it
    rangeMin = min([i[0] - 3*i[1] for i in x])
    rangeMax = max([i[0] + 3*i[1] for i in x])

    if showDetail:
        color = 'dodgerblue'
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

    ax.legend(loc='upper left', ncol=3)
    ax.set_xlabel(xlabel)
    # ax.set_ylabel(ylabel)

    sns.despine(fig=fig, ax=ax, left=True, offset=-5)
    ax.set_yticks([])
    ax.set_xlim(right=100)

    if title:
        fig.suptitle(title)

    if show:
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

def linear(x, m, t):
    return m*np.array(x) + t

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
        f.write('#PBS -l nodes=1:ppn=%d,walltime=23:59:59\n' % ppn)
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

def latexExport(file, header, data, caption, label):
    f = open(str(file), 'w')
    f.write('\\begin{table}[!htp]\n\t\\centering\n\t\\begin{tabular}{')
    f.write('c')
    for i in range(1, len(header)):
        f.write('c')
    f.write('}' + '\n\t\t\\toprule ' + header[0] + ' ')
    for i in range(1, len(header)):
        f.write('& %s ' % str(header[i]))

    f.write('\\\\\n\t\t\\midrule ' + str(data[0][0]) + ' ')
    for i in range(1, len(data[0])):
        f.write('& %s ' % str(round(data[0][i], 5) if isinstance(data[0][i], float) else data[0][i]))

    for obj in data[1:len(data)+1]:
        f.write('\\\\\n\t\t')
        f.write(str(obj[0]) + ' ')
        for i in range(1, len(obj)):
            f.write('& %s ' % str(round(obj[i], 5) if isinstance(obj[i], float) else obj[i]))
    f.write('\\\\\\bottomrule\n\t\\end{tabular}\n')
    f.write('\t\\caption{' + str(caption) + '}\n\t\\label{tab:' + str(label) + '}\n')

    f.write('\\end{table}\n')
    f.close()

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

