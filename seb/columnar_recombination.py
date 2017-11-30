#!/usr/bin/env python
import argparse
import multiprocessing
import numpy as np
from numpy.random import normal, uniform
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

from scipy import spatial
from scipy.optimize import curve_fit
import sys
import os.path
import time
import cPickle

# Doesn't work when run on cluster
try:
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib import rc
    from matplotlib.ticker import NullFormatter
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    import seaborn as sns
    sns.set()
    sns.set_style('white')
    sns.set_style('ticks')

except:
    pass

# Constants
X = 3.8e4 # V/m, electric field
Dt = 55.e-4 # m^2/s
Dl = .1 * Dt # 12.9399369923e-4 # 0.1 * Dt
vD = 1705 # m/s
mobility = vD / X # mobility
# X = 0.

b = 1.e-9 # 1.6716e-8 # 1.0813e-6 # 2.6e-6 
d = 2.34e-3 # 6.e-6 # 5.e-4 # 6.e-7 # 1.e-5 # Length of the column
A = .1e-9 # .5e-9 # To be optimized!
angle = 0

# Number of particles
N = int( 6.4e4 ) # 200
print N
# N = 5000

# estimated number of survivors (for bottom plot)
N_rem = 1000

# Time step size
dt = 6.e-16 # 4.3e-17 # (10.e-9 / 3)**2 * 1./(2 * (2*Dt + Dl)) # 2.5e-18 # 1.6e-16 # s
print dt

# Status bar
STATUSBARWIDTH = 20

# How often to repeat a simulation in order to estimate
# a statistical error
REPEAT = 1

def main():
    # Load command line arguments
    uniform, parallel, show, save, MAKEPLOTS, TEST, PRINTRESULTS, SINGLE, VAL_A, VAL_B = get_args()
    global SHOW
    global SAVE
    global PARALLEL
    global GAUSS

    # Show plot during drift?
    SHOW = show
    # Save particle positions?
    SAVE = save

    # Parallel or perpendicular?
    PARALLEL = parallel
    # Gaussian or uniform distribution?
    GAUSS = not uniform

    # Declare paths
    global gifOut
    global pickleOut
    global resultOut
    global resOut
    global stdOut
    global stdOutVal
    global D

    global INTERRUPT
    if PARALLEL:
        # D = [Dt] * 3
        D = [Dt, Dt, Dl]
        INTERRUPT = 500      # Limit for interruption
        gifOut = 'cr/cr_gif_par'
        pickleOut = 'cr/cr_par'
        resultOut = 'cr/res_par_test'
        resOut = 'cr/res_out_par/'
        stdOut = 'cr/par.log'
        stdOutVal = resOut + 'log/'
    else:
        # D = [Dl] * 3
        D = [Dt, Dt, Dl]
        INTERRUPT = 500      # Limit for interruption
        gifOut = 'cr/cr_gif_perp'
        pickleOut = 'cr/cr_perp'
        resultOut = 'cr/res_perp_test'
        resOut = 'cr/res_out_perp/'
        stdOut = 'cr/perp.log'
        stdOutVal = resOut + 'log/'

    if not GAUSS:
        gifOut += '_uni'
        pickleOut += '_uni'
        resultOut += '_uni'

    pickleOut += '.p'
    resultOut += '.p'

    # For printresults
    resultOutList = ['cr/res_par.p', 'cr/res_par_uni.p', 'cr/res_perp.p', 'cr/res_perp_uni.p']
    # titleList = ['Parallel (Gauss)', 'Parallel (uniform)', 'Perpendicular (Gauss)', 'Perpendicular (uniform)']
    titleList = ['Parallel', 'Parallel', 'Perpendicular', 'Perpendicular']
    colorList = ['blue', 'blue', 'green', 'green']

    # === MAKE PLOTS ===
    if MAKEPLOTS:
        dic = cPickle.load(open(pickleOut, 'rb'))
        plotDict(dic, SHOW)
        return

    # === PROCESS SINGLE ===
    if SINGLE:
        global b
        global A
        inter = 100
        generate(N, d, b, A, angle=None, interrupt=INTERRUPT, out_q=None)
        return

    # === TEST ===
    if TEST:
        # driftElectronTest()
        dt = 1.e-6
        driftElectronTest()
        raw_input('')
        return

    # === PRINT RESULTS ===
    if PRINTRESULTS:
        printResults2D(resultOutList, titleList)
	return

        # Figure of recombination vs. b
        figRes, axRes = plt.subplots(figsize=(7, 4))
        figRes.subplots_adjust(bottom=0.15)

        resListList = []
        for m, resultOut in enumerate(resultOutList):
            if not os.path.isfile(resultOut):
                print resultOut
                print 'continue'
                continue

            # Figure of recombination vs. t
            # fig, ax = plt.subplots()
            fig = plt.figure(figsize=(7, 4))
            # fig.subplots_adjust(bottom=0.2, right=0.2)
            ax = fig.add_axes([0.1, 0.14, 0.7, 0.8])
            ax.set_xlim(xmin=0, xmax=1.e-3)
            # ax.set_ylim(ymin=0)
            sns.despine(fig=fig, ax=ax)
            # Add colorbar
            axCBar = fig.add_axes([0.85, 0.1, 0.05, 0.8])

            dic = cPickle.load(open(resultOut, 'rb'))
            bListOld, percListListOld, timeListListOld = dic['b'], dic['p'], dic['t']
	    bList, percListList, timeListList = [], [], []
            print bListOld
	    for k, Ab in enumerate( bListOld ): 
		A, b = Ab
		if A == 5.9e-9:
			bList.append( b )
			percListList.append( percListListOld[k] )
			timeListList.append( timeListListOld[k] )

            print bList
            bList, percListList, timeListList = zip(*sorted(zip(bList, percListList, timeListList)))

            # DEBUG!
            getColorBar(axCBar, bList[0], bList[-1], len(bList), 'b [nm]')

            resList = []
            for i, b in enumerate( bList ):
                timeList, percList = timeListList[i], percListList[i]
                color = getColor(len(bList), i)

                # Fit the curves in the parallel case
                if True: # m in [0, 1]:
                    # Length of the column: d 
		    # tt = 1./1000 * d / vD * 1.e9
                    dt = 3.e-16
                    tt = 1000 * dt * 1.e9

                    # Get only mean values
                    if not isinstance(percList[0], float):
                        pList = np.array(percList)[:,0]
                    else:
                        pList = percList

                    # p0 = [percList[-1], -5, 0.5]
                    popt, pcov = curve_fit(parallelFit, timeList, pList) # , p0=p0)
                    perr = np.sqrt( np.diag(pcov) )
                    print popt, perr
            
                    timeListFit = timeList # np.linspace(0, tt, float(len(timeList))/timeList[-1]*tt)
                    print timeListFit
                    percListFit = parallelFit(timeListFit, *popt)

                    resList.append( (1 - popt[0], perr[0]) ) # percListFit[-1] )
    
                    plotFigure(fig, ax, timeListFit, np.array(percListFit)*100, xlabel=r'time [ns]', ylabel=r'Recombination fraction [\%]', label=None, title=titleList[m], color=color, lw=1.5)
                
                # PERPENDICULAR
                else:
                    resList.append( percList[-1] )

                # No label, use colorbar
                plotFigure(fig, ax, timeList, np.array(percList)*100, xlabel=r'time [ns]', ylabel=r'Recombination fraction [\%]', label=None, title=titleList[m], color=color)
            ax.set_ylim(ymin=0)

            # Show diffusion loss
            if not isinstance(resList[0], float):
                resList = [(1-r[0], r[1]) for r in resList]
            else:
                resList = [1 - r for r in resList]

            print 'Results'
            print bList
            print resList
            print
            
            # Fit result curve
            try:
                p0 = [.08, .3]
                popt, pcov = curve_fit(lambda x, a, b: parallelFit(x, 1., a, b), bList, resList, p0=p0)
                perr = np.sqrt( np.diag(pcov) )
            except:
                popt = p0
                perr = [0, 0]
                # continue

            print 'Fit result', popt, perr
            a, b = popt
        
            bListFit = np.linspace(0, bList[-1], 1000) 
            resListFit = parallelFit(bListFit, 1., a, b)

            # plotFigure(figRes, axRes, bListFit, np.array(resListFit)*100, xlabel=r'b [m]', ylabel=r'Diffusion loss [\%]', color=colorList[m])
            sns.despine(fig=figRes, ax=axRes)
            plotFigure(figRes, axRes, np.array(bList)*1.e9, np.array(resList)*100, xlabel=r'b [nm]', ylabel=r'Recombination fraction [\%]', label=titleList[m], color=colorList[m])

        raw_input('')
        return

    if VAL_A or VAL_B:
        getResults('Ab', val_A=VAL_A, val_b=VAL_B)
    else:
        getResults('Ab')
    return

# === PRINT RESULTS 2D ===
def printResults2D(resultOutList, titleList):
    from scipy.interpolate import griddata

    # Loop over parallel and perpendicular
    results = []
    factorList = []
    for m, resultOut in enumerate(resultOutList):
        if not os.path.isfile(resultOut):
            continue

        dic = cPickle.load(open(resultOut, 'rb'))
        paramList, percListList, timeListList = dic['b'], dic['p'], dic['t']
        # Generate points on a grid
        resList = []
        AList, bList = [], []
        for i, param in enumerate( paramList ):
            timeList, percList = timeListList[i], percListList[i]
            A, b = param
            if A > b: 
                continue 
            AList.append( A ), bList.append( b )

            print A, b
            tt = d/vD * 1.e9
            if not isinstance(percList[0], float):
                pList = np.array(percList)[:,0]
            else:
                pList = percList

            p0 = [.5, 3.46205358e-04, 1.19802176e+00]
            try:
                popt, pcov = curve_fit(parallelFit, timeList, pList, p0=p0)
                perr = np.sqrt( np.diag(pcov) )
            except:
                popt = np.array(p0)
                popt[0] = pList[-1]
                perr = np.array([0] * len(popt))
            print popt, perr
            factorList.append( popt[0] / percList[-1] )

            timeListFit = np.linspace(0, 2.e-3, int(10e6)) # float(len(timeList))/timeList[-1]*tt)
            #  timeListFit = timeList
            percListFit = parallelFit(timeListFit, *popt)
            # plt.plot(timeListFit, percListFit)
            # plt.plot(timeList, percList)
            # plt.show()

            # Get value of fit for t -> oo
            resList.append( percList[-1] * 100 ) # popt[0] * 100 )

        # Convert to nm
        AList, bList, resList = np.array(AList)*1.e9, np.array(bList)*1.e9, np.array(resList)
        xi, yi = np.mgrid[AList.min():AList.max():20j, bList.min():bList.max():20j]
        a_rescale = griddata((AList, bList), resList, (xi, yi), method='cubic') 
        # a_rescale = rescaled_interp(AList, bList, resList, xi, yi)

        plot2D(AList, bList, resList, a_rescale, title=titleList[m])

        results.append( zip(zip(AList, bList), resList) )

    # Parallel first, then perpendicular
    parRes, perpRes = results
    parRes, perpRes = sorted(parRes), sorted(perpRes)
    parCoord, parR = zip(*parRes)
    perpCoord, perpR = zip(*perpRes)

    residual = (np.array(perpR) - np.array(parR)) / np.array(perpR)
    print residual
    AList, bList = zip(*perpCoord)
    AList, bList = np.array(AList), np.array(bList)

    xi, yi = np.mgrid[AList.min():AList.max():20j, bList.min():bList.max():20j]
    a_rescale = griddata((AList, bList), residual, (xi, yi), method='cubic') 

    # print a_rescale.shape
    # a_rescale = rescaled_interp(AList, bList, residual, xi, yi)
    plot2D(AList, bList, residual, a_rescale, cmap='RdBu_r', title='', residual=True)

    print np.mean(factorList), np.std(factorList)

    return

def plot2D(x, y, a, ai, cmap='Blues', title=None, residual=False):
    from scipy.interpolate import griddata
    fig, ax = plt.subplots()

    if residual:
        vmin, vmax = -1, 1
    else:
        vmin, vmax = 0, 100

    if residual:
        xi, yi = np.mgrid[x.min():x.max():1000j, y.min():y.max():1000j]
        a_rescale = griddata((x, y), a, (xi, yi), method='cubic') 
        scale = 0.7 # np.max(abs(ai))

        from matplotlib import colors
        im = ax.imshow(a_rescale.T, origin='lower', cmap=cmap, extent=[x.min(), x.max(), y.min(), y.max()], aspect='auto', vmin=-scale, vmax=scale)
        cs = ax.contour(xi, yi, a_rescale, 10, colors='k', vmin=-1., vmax=1, levels=[-0.12, -0.08, -0.04, 0., 0.16, 0.32, 0.48], linestyles='solid')
        ax.clabel(cs, inline=1, fontsize=10)
    else:
        im = ax.imshow(ai.T, origin='lower', cmap=cmap, extent=[x.min(), x.max(), y.min(), y.max()], aspect='auto', vmin=0, vmax=100)

    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())

    ax.scatter(x, y, c='k')
    ax.set(xlabel=r'$A$ [nm]', ylabel=r'$b$ [nm]', title=title)
    cb = fig.colorbar(im)

    if residual:
        cb.set_label(r'$(R_\mathrm{perp} - R_\mathrm{par}) / R_\mathrm{perp}$')
    else:
        cb.set_label(r'Recombination fraction [\%]')

    fig.show()
    raw_input('')

def normal_interp(x, y, a, xi, yi):
    import scipy.interpolate
    rbf = scipy.interpolate.Rbf(x, y, a)
    ai = rbf(xi, yi)
    return ai

def rescaled_interp(x, y, a, xi, yi):
    a_rescaled = np.nan_to_num( (a - a.min()) / a.ptp() )
    ai = normal_interp(x, y, a_rescaled, xi, yi)
    ai = a.ptp() * ai + a.min()
    return ai

# === GET RESULTS ===
# Loop over parameters and store results to dict
def getResults(param='b', show=False, val_A=None, val_b=None):
    print val_A, val_b
    if val_A:
        if val_b:
            stdOut = stdOutVal + 'A%d_b%d' % (val_A*1.e9, val_b*1.e9) + '.log'
        else:
            stdOut = stdOutVal + 'A%d' % (val_A*1.e9) + '.log'
    elif val_b:
        stdOut = stdOutVal + 'b%d' % (val_b*1.e9) + '.log'
    else:
        global stdOut

    stdFile = open(stdOut, 'w')
    sys.stdout = stdFile

    # Run over bList
    if show:
        # Create figure to use in plots
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.7, 0.8])
        # Add colorbar
        axCBar = fig.add_axes([0.85, 0.1, 0.05, 0.8])

    # bList = np.linspace(1.e-7, 30.e-7, 10)
    if val_A:
        if val_b:
            paramList = [(val_A, val_b)]
            param == 'Ab'
        else:
            paramList = [val_A]
            param == 'A'
    elif val_b:
        paramList = [val_b]
        param == 'b'

    elif param == 'b':
        paramList = np.linspace(.5e-9, 2.5e-9, 5) # np.linspace(1.e-8, 30.e-8, 10)
        colorBarLabel = 'b [m]'
    elif param == 'A':
        paramList = np.linspace(1.5e-9, 1.e-9, 5) # np.linspace(60.e-9, 100.e-9, 5)
        colorBarLabel = 'A [m]'
    elif param == 'Ab':
        colorBarLabel = 'Ab [m]'
        # paramList = zip(np.linspace(.1e-9, 100.e-9, 3), np.linspace(.1e-9, 100.e-9, 3)) 
        AList, bList = np.linspace(1.e-9, 25.5e-9, 6), np.linspace(1.e-9, 25.5e-9, 6)
        paramList = [(A, b) for A in AList for b in bList if A <= b]
    elif param == 'angle':
        paramList = np.linspace(0., np.pi, 5)
        colorBarLabel = 'Theta (rad)'

    #bList = [1.e-7, 2.e-7]
    if show:
        # getColorBar(axCBar, bList[0], bList[-1], 'b [m]')
        getColorBar(axCBar, paramList[0], paramList[-1], len(paramList), 'b [nm]')

    manager = multiprocessing.Manager()
    out_q = manager.Queue()

    procs = []
    # for b in bList:
    for para in paramList:
        if param == 'b':
            global angle
            global A
            b = para
        elif param == 'A':
            global b
            global angle
            A = para
        elif param == 'Ab':
            global angle
            A, b = para
        elif param == 'angle':
            global b
            global A
            angle = para

        if REPEAT > 1:
            p = multiprocessing.Process(target=generateRepeat, args=(N, d, b, A, angle, INTERRUPT, out_q, REPEAT, param))
        else:
            p = multiprocessing.Process(target=generate, args=(N, d, b, A, angle, INTERRUPT, out_q, param))

        procs.append(p)
        p.start()

    resultDict = {}
    resultList = []
    # bResList = []
    paramResList = []
    timeListList, percListList = [], []

    # for i in range(len(bList)):
    for i in range(len(paramList)):
        resultDict.update(out_q.get())

        # bRes = resultDict['b']
        # bResList.append( bRes )
        paramRes = resultDict['b']
        paramResList.append( paramRes )
        timeList, percList = resultDict['t'], resultDict['p']
        timeListList.append( timeList ), percListList.append( percList )

        if show:
            # color = getColor(len(bList), list(bList).index(bRes))
            color = getColor(len(paramList), list(paramList).index(paramRes))
            plotFigure(fig, ax, timeList, np.array(percList)*100, xlabel=r'time [ns]', ylabel=r'Recombination fraction [\%]', color=color)

    for p in procs:
        p.join()

    # Order lists
    # zipped = zip(bResList, timeListList, percListList)
    zipped = zip(paramResList, timeListList, percListList)
    zipped.sort()
    # bResList, timeListList, percListList = zip(*zipped)
    paramResList, timeListList, percListList = zip(*zipped)

    # Store lists to file
    # resultDict = {'b': bResList, 'p': percListList, 't': timeListList}
    resultDict = {'b': paramResList, 'p': percListList, 't': timeListList}
    if val_A:
        if val_b:
            resultOut = resOut + 'A%d_b%d' % (val_A*1.e9, val_b*1.e9) + '.p'
        else:
            resultOut = resOut + 'A%d' % (val_A*1.e9) + '.p'
    elif val_b:
            resultOut = resOut + 'b%d' % (val_b*1.e9) + '.p'
    else:
        global resultOut

    print resultOut
    cPickle.dump(resultDict, open(resultOut, 'wb'))

    if show:
        raw_input('')

    stdFile.close()

def generateRepeat(N, d, b, A, angle=None, interrupt=False, out_q=None, repeat=1, out='Ab'):
    percListList = []
    for k in range(repeat):
        timeList, percList = generate(N, d, b, A, angle, interrupt, None)
        percListList.append( percList )

    # Get list with minimum length
    minLen = min([len(p) for p in percListList])
    # Cut other lists to same length
    timeList = timeList[:minLen]
    percListList = [p[:minLen] for p in percListList]
    # Get mean and std of all points
    percList = [(np.mean(p), np.std(p)) for p in zip(*percListList)]

    if out_q:
        outDict = {'b': b, 't': timeList, 'p': percList}
        out_q.put( outDict )

    return timeList, percList

def generate(N, d, b, A, angle=None, interrupt=False, out_q=None, param='b'):
    if True: # PARALLEL:
        tt = d / vD
    else:
        tt = 2 * 5*b / vD

    tt = 2.e-12 # 6000 * dt # 1.e-9
    tt = 100 * dt

    # DEBUG
    # tt = 1.e-12

    # Generate Xe+ positions
    print 'Generate initial Xe+ positions...',
    # XePosList = initialDistribution(N, d, b, angle, PARALLEL, GAUSS)
    XePosList, ePosList = initialDistributionStraight(N, d, parallel=PARALLEL, dist=b)
    print 'Done!'

    # Create cKD tree
    print 'Fill Xe+ KD tree...',
    xyz = np.c_[np.array(XePosList)[:,0], np.array(XePosList)[:,1], np.array(XePosList)[:,2]]
    XeTree = spatial.cKDTree(xyz)
    print 'Done!'

    meanDist = getMeanDistance(XePosList, XeTree)
    print 'Mean distance of two Xe+ ions:', meanDist
    print 'Minimum time step length:', meanDist / vD
    # if dt > meanDist / vD:
    #    print 'dt is too small!'
    #    return

    # Generate electrons and check for intersection
    print 'Generate initial e positions...',
    # ePosList = initialDistribution(N, d, b, angle, parallel=PARALLEL, gauss=GAUSS) #, tree=XeTree, A=A)
    print 'Done!'
    print 'Fill e KD tree...',
    xyz = np.c_[np.array(ePosList)[:,0], np.array(ePosList)[:,1], np.array(ePosList)[:,2]]
    eTree = spatial.cKDTree(xyz)
    print 'Done!'
    print

    # Data for bottom 
    timeList = [0]
    NList = [N]
    percList = [0.]

    if SHOW:
        # === PLOT INIT ===
        fig, h, hBottom, hTop, hRight = plotInit(XePosList, timeList, NList, tt)

    if SAVE:
        posDict = {}
        # Store initial Xe+ positions
        posDict['Xe'] = {0: XePosList}
        # Store initial e positions
        posDict['e'] = {0: ePosList}
    
        NListSave = [N]
        tListSave = [0]

    # == Time loop == 
    # Initial distribution has no intersection, 
    # therefore start with moving the electrons
    print '=== TIME LOOP START ==='
    print 'Total drift duration t = %f ns' % (tt*1.e9)
    loopStart = time.clock()
    loop = list(enumerate(np.arange(dt, tt, dt)))
    for i, t in loop:
        # Show progress
        sys.stdout.write('\r')
        perc = int(float(i)/len(loop)*100)
        sys.stdout.write("[%-20s] %d%% (Ne = %d)" % ('='*int(0.2*perc), perc, len(ePosList)))
        sys.stdout.flush()

        '''
        # Diffusion of Xe-ions
        XePosNewList = []
        for XePos in XePosList:
            XePosNew = driftElectron(XePos, D, -vD, dt, eTree, XeTree, ePosList, XePosList)
            XePosNewList.append( XePosNew )
        XePosList = XePosNewList
        '''

        # Loop over all electrons
        ePosListNew = []
        XeRemoveIdx = []
        s = time.clock()
        for ePos in ePosList:
            # Drift the electron
            # [Dt, Dt, Dl]
            ePosNew = driftElectron(ePos, D, vD, dt, eTree, XeTree, ePosList, XePosList)

            # Check for intersection only if electron is close
            # enough to the column
            r = getRadius(ePosNew, -vD*t)
            # if r > 5*b or (PARALLEL and (ePosNew[2] < -.6*d)):
            if r > A or (PARALLEL and (ePosNew[2] < -(d + A))):
                ePosListNew.append( ePosNew )
                continue

            # Check for intersection with Xe
            while True:
                intersect, idxList = checkIntersect(ePosNew, XeTree, A, None, ePos, XePosList)
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

                    # else:
                    # Electron hit already recombined Xe+
                    #    k += 1
                else:
                    # Electron survives
                    ePosListNew.append( ePosNew )
                    break
                
        print time.clock() - s

        ePosList = ePosListNew
        # Renew the XeTree
        XePosList = list([XePos for u, XePos in enumerate(XePosList) if u not in XeRemoveIdx])

        # Save only every 100th frame
        # if not i % 10:
        if SAVE and not i % 40 and i != 0:
            print t
            print N
            # Store Xe+ & e positions
            posDict['Xe'][t] = XePosList
            posDict['e'][t] = ePosList

            NListSave.append(len(XePosList))
            tListSave.append(t)

        # Works faster than expected
        xyz = np.c_[np.array(XePosList)[:,0], np.array(XePosList)[:,1], np.array(XePosList)[:,2]]
        XeTree = spatial.cKDTree(xyz)
        xyz = np.c_[np.array(ePosList)[:,0], np.array(ePosList)[:,1], np.array(ePosList)[:,2]]
        eTree = spatial.cKDTree(xyz)
        print time.clock() - s 

        # Number of electrons vs. time
        timeList.append(t*1.e9)
        NList.append(len(ePosList))
        p = 1 - float(len(ePosList)) / N
        percList.append( p )

        # Test on every 10th iteration if change in 
        # recombination is small. If so, interrupt
        if interrupt and (not i % 10):
            if len(NList) > interrupt:
                # slope = np.mean([dt/(percList[n+1]-percList[n]) if percList[n] != percList[n+1] else 0 for n in range(len(percList)-100, len(percList), 2)]) / 100
                mean = .5*(NList[-interrupt] + NList[-2])
                # if abs(slope) < 1.e-13:
                if mean == NList[-1]:
                    print
                    print 'In region of small changes -> interrupting'
                    break

        if SHOW:
            if PARALLEL:
                xe, ye = np.array(ePosList)[:,0], np.array(ePosList)[:,1]
                xE, yE = np.array(XePosList)[:,0], np.array(XePosList)[:,1]
            else:
                xe, ye = np.array(ePosList)[:,1], np.array(ePosList)[:,2]
                xE, yE = np.array(XePosList)[:,1], np.array(XePosList)[:,2]
            rE = np.sqrt(xE**2 + yE**2)
            # print t, np.mean(rE), np.std(xE), np.std(yE), np.sqrt(np.std(xE)**2 + np.std(yE)**2)

            updatePlot(fig, h, hBottom, hTop, hRight, t, i, timeList, np.array([float(n) for n in NList]) / N * 100, ePosList, XePosList, show=SHOW)

    if SHOW:
        # plt.hist(xE, 50)
        plt.hist(ye, 50)
        plt.hist(yE, 50, alpha=.5)
        plt.hist(rE, 50, alpha=.5)
        plt.show()

    # After the time loop is done, count remaining particles 
    print
    print 'Calculation time:', time.clock() - loopStart 
    print 'Electrons remaining:', len(ePosList)
    print 'Recombination percentage:', 1 - float(len(ePosList))/N

    if SAVE:
        # Save Xe+ & e positions to disk
        cPickle.dump(posDict, open(pickleOut, 'wb'))

        plotFacet(posDict['Xe'], posDict['e'], tListSave, NListSave)

    if out_q:
        print
        print 'Reached end'
        # DEBUG: b -> A

        if param == 'b':
            outDict = {'b': b, 't': timeList, 'p': percList}
        elif param == 'A':
            outDict = {'b': A, 't': timeList, 'p': percList}
        elif param =='Ab': 
            outDict = {'b': (A, b), 't': timeList, 'p': percList}
            
        out_q.put( outDict )

    return timeList, percList

def driftElectronTest(tt=1.e-14, dt=1.e-18):
    import seaborn as sns
    paper_rc = {'lines.linewidth': 1., 'lines.markersize': 5, 'lines.markeredgecolor': 'auto', 'lines.markeredgewidth': 1.5}
    sns.set_context("paper", rc = paper_rc)                                    
    pos = (0, 0, 0)

    x, dx = [], []
    y, dy = [], []
    z, dz = [], []
    time = np.arange(0, tt, dt)
    for i in time:
        # [Dt, Dt, Dl]
        posNew = driftElectron(pos, D, vD, dt)
        x.append( posNew[0] ), dx.append( posNew[0] - pos[0] )
        y.append( posNew[1] ), dy.append( posNew[1] - pos[1] )
        z.append( posNew[2] ), dz.append( posNew[2] - pos[2] )
        pos = posNew

    print 'Mean:', np.mean(x)
    print 'Std:', np.std(x)
    print 'Std (theory):', np.sqrt(2*Dt*tt)

    f, ax = plt.subplots(figsize=(2, 12))
    # f.subplots_adjust(bottom=.08, left=.3, top=.95)
    f.subplots_adjust(bottom=.08, left=.35, top=.95)
    # plt.tight_layout()
    plt.gca().invert_yaxis()

    time = np.array( time ) * 1.e12
    x, y, z = np.array(x)*1.e9, np.array(y)*1.e9, np.array(z)*1.e9
    ax.plot(x, time, label=r'x')
    ax.plot(y, time, label=r'y')
    ax.plot(z, time, label=r'z')
    # plt.plot(time, dx)

    ax.set_ylim(ymax=0.)
    ax.set_xlabel(r'Position [nm]', fontsize=14)
    ax.set_ylabel(r'Time [ps]', fontsize=14)
    plt.legend(loc='best')
    sns.despine(fig=f, ax=ax, left=True, bottom=True, offset=-100)

    f.savefig('driftTest.pdf', format='pdf')
    f.show()

def getParticle(d, b, parallel, gauss, angle=None):
    # Rotation matrix
    if angle:
        c, s = np.cos(angle), np.sin(angle)
        R = np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

    if gauss:
        if parallel:
            # Get the x and y position
            x, y = np.random.normal(0, b, 2)

            # Get the z position
            z = np.random.uniform(-.5*d, .5*d)

        # perpendicular
        else:
            x = np.random.uniform(-.5*d, .5*d)
            y, z = np.random.normal(0, b, 2)
            # if angle:
            #    return tuple(np.dot(R, np.array([x, y, z])))

    # Uniform
    else:
        if parallel:
            x, y = getRandomRadius(3*b)
            z = np.random.uniform(-.5*d, .5*d)
        else:
            x = np.random.uniform(-.5*d, .5*d)
            y, z = getRandomRadius(3*b)

    return (x, y, z)

def getRandomRadius(r):
    R = r + 1
    while R > r:
        x, y = np.random.uniform(-r, r, 2)
        R = np.sqrt(x**2 + y**2)
    return x, y

# Create particles along an axis pointing in positive z-direction
# Check intersection with solid particles in tree
# N - number of particles
# d - length of the axis in cm
# b - standard deviation of the distance to the axis
def initialDistribution(N, d, b, angle=None, parallel=True, gauss=True, tree=None, A=None):
    if tree and not A:
        print 'Please pecify the recombination radius when using a tree!'
        return False

    posList = []
    for i in range(N):
        x, y, z = getParticle(d, b, parallel, gauss, angle)
        if tree:
            intersect, idx = checkIntersect((x, y, z), tree, A)
            while intersect:
                x, y, z = getParticle(d, b, parallel, gauss, angle)
                intersect, idx = checkIntersect((x, y, z), tree, A)

        # posList.append( getParticle(d, b, parallel, gauss, angle) )
        posList.append((x, y, z))

    return posList

def initialDistributionStraight(N, d, parallel=True, dist=None, gauss=False):
    XePosList, ePosList = [], []
    for i in range(N):
        if parallel:
            XePos = np.array( [0, 0, np.random.uniform(-.5*d, .5*d)] )
        else:
            XePos = np.array( [np.random.uniform(-.5*d, .5*d), 0, 0] )
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

def checkIntersect(ePos, XeTree, A, k=1, ePosOld=None, XePosList=None):
    ePos = np.array(ePos)

    # Point - sphere intersection
    if not ePosOld:
        if not k:
            k = 10

        searchDist = np.inf

        dist, idx = XeTree.query(np.array(ePos), k=k, n_jobs=0, distance_upper_bound=searchDist)
        dist, idx = dist[0], idx[0]
        if not dist:
            return False 

        if not isinstance(dist, float):
            dist, idx = float(dist[-1]), int(idx[-1])

        if dist <= A:
            return True, [idx]
        else:
            return False, [0]

    # Point - area intersection
    else:
        ePosOld = np.array(ePosOld)
        posDist = np.linalg.norm(ePos - ePosOld)
        # Central point
        center = ePos + .5*(ePosOld - ePos)
        
        # Search radius
        searchDist = .5*posDist + A

        if not k:
            k = 10

        # Get k nearest neighbours
        dist, idx = XeTree.query(np.array(center), k=k, n_jobs=0, distance_upper_bound=searchDist)

        # Filter distances and indices
        elemList = np.array( [item for item in zip(idx, dist) if item[0] != len(XePosList)] )
        if len(elemList) == 0:
            return False, [0]
        else:
            dist, idx = elemList[:,1], elemList[:,0].astype('int')
   
        # Check if electron is close to a Xe+
        if dist[0] <= A:
            idxIntersect = [idx[0]]
        else:
            idxIntersect = []

        # Check if electron track intersects with a Xe+
        a = posDist
        for j, i in enumerate(idx[1:]):
            b = np.linalg.norm(ePosOld - np.array(XePosList[i]))
            c = dist[j+1]

            h = getAltitude(a, b, c)
            if h and h <= A:
                idxIntersect.append(i)
        else:
            if idxIntersect:
                return True, idxIntersect
            else:
                return False, [0]

# https://en.wikipedia.org/wiki/Altitude_(triangle)
def getAltitude(a, b, c):
    s = .5 * (a + b + c)
    try:
        h = 2 * np.sqrt(s*(s-a)*(s-b)*(s-c)) / a
    except:
        return None

    return h

def driftElectron(pos, D, vD, dt, eTree=None, XeTree=None, ePosList=None, XePosList=None):
    if vD == 0:
        X = 0
        mobility = 0.3 # m^2/(V*s)
    else:
        global X
        global mobility

    newPos = []
    for i, p in enumerate(pos):
        sigma = np.sqrt( 2*D[i]*dt )
        newPos.append( p + np.random.normal(0, sigma, 1)[0] )

    # Calculate couloumb interactions
    if eTree and XeTree and ePosList and XePosList:
        from scipy import constants
        epsilon = 1.95
        coulombRange = 5.e-9 # 49.e-9 # old(1.e-8), Onsager Radius

        eDist, eIdx = eTree.query(np.array([pos]), n_jobs=0, k=10, distance_upper_bound=coulombRange)
        # print eDist, eIdx
        # Remove first element, as it is the electron itself
        eDist, eIdx = eDist[0][1:], eIdx[0][1:]
        eDist, eIdx = filterQuery(eDist, eIdx, len(ePosList))
        
        XeDist, XeIdx = XeTree.query(np.array([pos]), n_jobs=0, k=10, distance_upper_bound=coulombRange)
        XeDist, XeIdx = XeDist[0], XeIdx[0]
        # print XeDist, XeIdx
        XeDist, XeIdx = filterQuery(XeDist, XeIdx, len(XePosList))

        constantTerm = constants.e/(4*np.pi*constants.epsilon_0*epsilon)
        eTerm = np.array( [0, 0, 0] )
        for i in range(len(eDist)):
            vec = np.array(ePosList[eIdx[i]]) - np.array(pos)
            eTerm = vec / (np.linalg.norm( vec )**3)

        XeTerm = np.array( [0, 0, 0] )
        for i in range(len(XeDist)):
            vec = np.array(XePosList[XeIdx[i]]) - np.array(pos)
            XeTerm = vec / (np.linalg.norm( vec )**3)

        if np.isnan(eTerm).any() or np.isnan(XeTerm).any():
            newPos[2] -= vD * dt
        else:
            eField = constantTerm * (XeTerm - eTerm) + X*np.array([0, 0, -1])
            newPos = tuple(np.array(newPos) + mobility * eField * dt)
    
    # No interaction, consider only electric field
    else:
        # Electric field in z-direction
        newPos = list(np.array(newPos) + np.array([0, 0, -1]) * vD * dt)

    return newPos

def getRadius(pos, axisHeight):
    if PARALLEL:
        x, y = pos[0], pos[1]
        xA, yA = 0, 0
    else:
        x, y = pos[1], pos[2]
        xA, yA = 0, axisHeight

    return np.sqrt((x-xA)**2 + (y-yA)**2)

# Use it with a tree made from posList
def getMeanDistance(posList, tree):
    distList = []
    for pos in posList:
        dist, idx = tree.query(np.array([pos]), k=2, n_jobs=4)
        dist = dist[0][-1]
        distList.append( dist )

    return sum(distList)/len(distList)

def spatialProbability(ePos, XePos, D, A, vD, t):
    xProb = 1./(np.sqrt(4*np.pi*D[0]))*np.exp(-(XePos[0] - ePos[0])**2/(4*D[0]*t))
    yProb = 1./(np.sqrt(4*np.pi*D[1]))*np.exp(-(XePos[1] - ePos[1])**2/(4*D[1]*t))
    zProb = 1./(np.sqrt(4*np.pi*D[2]))*np.exp(-(XePos[2] - ePos[2] - vD*t)**2/(4*D[2]*t))

    return xProb * yProb * zProb

# === PLOT ===
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

def plotFacet(XeList, eList, timeList, NList):
    import seaborn as sns
    sns.set(style='ticks')
    paper_rc = {'lines.linewidth': 1.5, 'lines.markersize': 5, 'lines.markeredgecolor': 'auto', 'lines.markeredgewidth': 1.5}
    sns.set_context("paper", rc = paper_rc)                                    
    rows, cols = 4, 5

    # fig = plt.figure(figsize=(8.27, rows * .25 * 11.69))
    fig = plt.figure()
    fig_single = plt.figure()
    ax_single = plt.axes(projection='3d')

    fig.set_size_inches(.7 * 8.27, rows * .25 * 11.69)
    #fig.set_size_inches(5, 9)
    fig.subplots_adjust(wspace=.1, hspace=0., top=1., bottom=.05)
    # plt.tight_layout()

    axRow = []
    for r in range(rows):
        axTopEmpty = [plt.subplot2grid((6*rows, cols), (r*6, i), projection='3d', rowspan=3) for i in range(cols)]
        axTop3D = [plt.subplot2grid((6*rows, cols), (r*6+1, i), projection='3d', rowspan=3) for i in range(cols)]
        axTop = plt.subplot2grid((6*rows, cols), (r*6+4, 0), colspan=cols)
        axRow.append( [axTop3D, axTop] )

    # Get data from input lists
    timeList = []
    XePosListList, ePosListList = [], []
    for key in XeList.keys():
        timeList.append( key )
        XePosListList.append( XeList[key] ), ePosListList.append( eList[key] )

    # Sort lists
    zipped = zip(timeList, XePosListList, ePosListList)
    zipped.sort()
    timeList, XePosListList, ePosListList = zip(*zipped)

    Nlist = [len(ePosListList[0])]
    N = Nlist[0]

    # Draw scatter
    for r in range(rows):
        axTop3D, axTop = axRow[r]
        for i, ax in enumerate(axTop3D):
            XePosList, ePosList = XePosListList[r*cols + i], ePosListList[r*cols + i]
            
            # ax = plt.axes(projection='3d')
            h = ax.scatter(0, 0, 0, marker='.', c='r')
            posPlot = np.array([list(x) for xs in zip(XePosList, ePosList) for x in xs])
            h._facecolor3d = [[1, 0, 0, 1], [0, 0, 1, 1]] * int(len(posPlot)*.5)
            h._edgecolor3d = h._facecolor3d
            h._offsets3d = (posPlot[:,0], posPlot[:,1], posPlot[:,2])

            h_single = ax_single.scatter(0, 0, 0, marker='.', c='r')
            h_single._facecolor3d = [[1, 0, 0, 1], [0, 0, 1, 1]] * int(len(posPlot)*.5)
            h_single._edgecolor3d = h._facecolor3d
            h_single._offsets3d = (posPlot[:,0], posPlot[:,1], posPlot[:,2])

            # Set Scale
            # Main
            if b > d:
                scale = 3*b
            else:
                scale = .3*d
            ax.set_xlim(-scale, scale)
            ax.set_ylim(-scale, scale)
            ax.set_zlim(-scale, scale)
            ax.set_aspect('equal')

            # Set coordinate system
            ax.set_axis_off()

            ax_single.set_xlim(-scale, scale)
            ax_single.set_ylim(-scale, scale)
            ax_single.set_zlim(-scale, scale)
            ax_single.set_aspect('equal')

            # Set coordinate system
            ax_single.set_axis_off()

            '''
            if PARALLEL:
                a = Arrow3D([0, 0], [0, 0], [-.5*d, .5*d], mutation_scale=20, arrowstyle='-|>', color='k')
            else:
                a = Arrow3D([-.5*d, .5*d], [0, 0], [0, 0], mutation_scale=20, arrowstyle='-|>', color='k', lw=1)
            ax.add_artist(a)
            '''
            ang = np.linspace(0, 2*np.pi, 100)
            xcir = 3*b*np.cos(ang)
            ycir = 3*b*np.sin(ang)
            zcir = np.zeros((100,))
            if PARALLEL:
                ax.plot(xcir, ycir, zcir, 'gray')
                ax_single.plot(xcir, ycir, zcir, 'gray')
            else:
                ax.plot(zcir, xcir, ycir, 'gray')
                ax_single.plot(zcir, xcir, ycir, 'gray')

            fig_single.show()
            raw_input('')

        # Draw time line
        axTop.plot(timeList[(r*cols):(r*cols + cols)], np.array(NList[(r*cols):(r*cols + cols)]).astype(float) / N * 100, marker='x')
        axTop.set_ylim(0, 120)
        print N, np.array(NList[(r*cols):(r*cols + cols)]).astype(float) / N * 100
        # axTop.set_ylabel(r'N [\%]')
        if r >= 1:
            sns.despine(ax=axTop, left=True)
        else:
            sns.despine(ax=axTop)

    else:
        axTop.set_xlabel(r'time [s]', fontsize=14)
    plt.ylabel(r'N [\%]')

    fig.show()
    fig.savefig('cr_sim.pdf', format='pdf')
    import matplotlib
    # fig.savefig('cr_sim.pdf', bbox_inches=matplotlib.transforms.Bbox(np.array(((0, 0), (2, 2)))), format='pdf')
    raw_input('')

def plotInit(XePosList, timeList, NList, tt):
    fig = plt.figure(figsize=(10, 5))

    # Plot area 
    left, width = 0.18, 0.7
    bottom, height = 0.27, 0.7

    # Main
    rect_main = [left, bottom, width, height]
    # ax =  fig.add_subplot(111, projection='3d')
    ax = plt.axes(rect_main, projection='3d')

    # Bottom
    rect_bottom = [0.1, 0.1, width+0.1, 0.13]
    axBottom = plt.axes(rect_bottom)

    # Top
    rect_top = [0., bottom+0.05, width-0.25, .6]
    axTop = plt.axes(rect_top)

    # Right
    rect_right = [width+0.1, bottom+0.05, 0.1, .6]
    axRight = plt.axes(rect_right)

    # Draw figure
    fig.canvas.draw()
    # Scatter plot of electrons and Xe+
    h = ax.scatter(np.array(XePosList)[:,0], np.array(XePosList)[:,1], np.array(XePosList)[:,2], marker='.', c='r')

    # Number of electrons vs. percentage
    hBottom = axBottom.plot(timeList, [100])

    # Initial electron distribution
    # in (z)
    if PARALLEL:
        x, y, z = np.array(XePosList)[:,0], np.array(XePosList)[:,1], np.array(XePosList)[:,2]
    else:
        x, y, z = np.array(XePosList)[:,1], np.array(XePosList)[:,2], np.array(XePosList)[:,2]

    if PARALLEL:
        eDist, bins = np.histogram(z, 50, range=(-1.1*d, 1.1*d), normed=True)
    else:
        eDist, bins = np.histogram(z, 50, range=(-6*b, 6*b), normed=True)
    bins = np.insert(bins, 0, 2*bins[0] - bins[1])
    bins = np.append(bins, 0)
    eDist = np.insert(eDist, 0, 0)
    eDist = np.append(eDist, 0)
    hRight = axRight.step(eDist, bins[:-1], where='pre')

    # in (x, y)
    H, xedges, yedges = np.histogram2d(x, y, bins=(30, 30), range=[(-7*b, 7*b), (-8*b, 8*b)], normed=True)
    H = H.T

    hTop = axTop.imshow(H, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=1, cmap='viridis')

    # Define scale
    if b > d:
        scale = 3*b
        ticks = [-3*b, -2*b, -b, 0, b, 2*b, 3*b]
        ticklabels = ['-3b', '-2b', '-b', '0', 'b', '2b', '3b']
    else:
        scale = .6*d
        ticks = [-.6*d, -.3*d, 0, .3*d, .6*d]
        ticklabels = ['-0.6d', '-0.3d', '0', '0.3d', '0.6d']

    # Set axis labels
    ax.set_xlabel('x', fontsize=14), ax.set_ylabel('y', fontsize=14), ax.set_zlabel('z', fontsize=14)
    axBottom.set_xlabel(r't [ns]', fontsize=14)
    axBottom.set_ylabel(r'N [\%]', fontsize=14)
    axRight.set_ylabel(r'z', fontsize=14)

    if PARALLEL:
        axTop.set_xlabel('x', fontsize=14)
        axTop.set_ylabel('y', fontsize=14)
    else:
        axTop.set_xlabel('z', fontsize=14)
        axTop.set_ylabel('y', fontsize=14)

    # Set ticks
    # Main
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels)
    ax.set_yticks(ticks)
    ax.set_yticklabels(ticklabels)
    ax.set_zticks(ticks)
    ax.set_zticklabels(ticklabels)

    # use radial ticks
    ticksRad = [-5*b, -2.5*b, 0, 2.5*b, 5*b]
    tickRadLabels = ['-5b', '-2.5b', '0', '2.5b', '5b']

    # Top
    axTop.set_xticks(ticksRad)
    axTop.set_xticklabels(tickRadLabels)
    axTop.set_yticks(ticksRad)
    axTop.set_yticklabels(tickRadLabels)

    # Right
    if PARALLEL:
        axRight.set_yticks(ticks)
        axRight.set_yticklabels(ticklabels)
    else:
        axRight.set_yticks(ticksRad)
        axRight.set_yticklabels(tickRadLabels)

    # Set Scale
    # Main
    ax.set_xlim(-scale, scale)
    ax.set_ylim(-scale, scale)
    ax.set_zlim(-scale, scale)
    ax.set_aspect('equal')

    # Set coordinate system
    '''
    ax.set_axis_off()
    if PARALLEL:
        a = Arrow3D([0, 0], [0, 0], [-.5*d, .5*d], mutation_scale=20, arrowstyle='-|>', color='k')
    else:
        a = Arrow3D([-.5*d, .5*d], [0, 0], [0, 0], mutation_scale=20, arrowstyle='-|>', color='k', lw=1)
    ax.add_artist(a)
    ang = np.linspace(0, 2*np.pi, 100)
    xcir = 3*b*np.cos(ang)
    ycir = 3*b*np.sin(ang)
    zcir = np.zeros((100,))
    if PARALLEL:
        ax.plot(xcir, ycir, zcir, 'gray')
    else:
        ax.plot(zcir, xcir, ycir, 'gray')
    '''

    # Right
    if PARALLEL:
        axRight.set_ylim(-d, d)
    else:
        axRight.set_ylim(-5*b, 5*b)

    # Bottom
    axBottom.set_xlim(0., tt*1.e9) # in [ns]
    # axBottom.set_yscale('symlog', nonposx='clip')
    # Show percentage
    axBottom.set_ylim(0, 100) 

    return fig, h, hBottom, hTop, hRight

def updatePlot(fig, h, hBottom, hTop, hRight, t, i, timeList, NList, ePosList, XePosList, show=False):
    # Append for bottom plot
    if PARALLEL:
        x, y, z = np.array(ePosList)[:,0], np.array(ePosList)[:,1], np.array(ePosList)[:,2]
    else:
        x, y, z = np.array(ePosList)[:,1], np.array(ePosList)[:,2], np.array(ePosList)[:,2]

    # Histogram for right plot
    if PARALLEL:
        eDist, bins = np.histogram(z, 50, range=(-1.1*d, 1.1*d), normed=True)

    else:
        eDist, bins = np.histogram(z, 50, range=(-6*b, 6*b), normed=True)
    # Histogram for top plot
    H, xedges, yedges = np.histogram2d(x, y+vD*t, bins=(30, 30), range=[(-8*b, 8*b), (-8*b, 8*b)], normed=True)
    H = H.T

    # === PLOT ===
    # Main
    posPlot = np.array([list(x) for xs in zip(XePosList, ePosList) for x in xs])
    h._facecolor3d = [[1, 0, 0, 1], [0, 0, 1, 1]] * int(len(posPlot)*.5)
    h._edgecolor3d = h._facecolor3d
    h._offsets3d = (posPlot[:,0], posPlot[:,1], posPlot[:,2])

    # Top
    hTop.set_data(H)
    
    # Bottom
    hBottom[0].set_data(timeList, NList)

    # Right
    bins = np.insert(bins, 0, 2*bins[0] - bins[1])
    bins = np.append(bins, 0)
    eDist = np.insert(eDist, 0, 0)
    eDist = np.append(eDist, 0)
    hRight[0].set_data(eDist, bins[:-1])

    # Figure title
    fig.suptitle('t = %.3f ns' % (t*1.e9))

    fig.canvas.draw()
    fig.canvas.flush_events()
    fig.savefig('%s/step_%d.png' % (gifOut, i), format='png', dpi=100)
            
    # plt.pause(0.000000000001)
    plt.draw()
    if show:
        fig.show()

# Plot data stored in dictionary
def plotDict(d, show=False):
    timeList = []
    XePosListList, ePosListList = [], []
    for key in d['e'].keys():
        timeList.append( key )
        XePosListList.append( d['Xe'][key] ), ePosListList.append( d['e'][key] )

    # Sort lists
    zipped = zip(timeList, XePosListList, ePosListList)
    zipped.sort()
    timeList, XePosListList, ePosListList = zip(*zipped)
    # sortedList = np.array(sorted(zip(timeList, XePosListList, ePosListList)))
    # timeList, XePosListList, ePosListList = sortedList[:,0], sortedList[:,1], sortedList[:,2]

    Nlist = [len(ePosListList[0])]

    print 'Start plotting'
    fig, h, hBottom, hTop, hRight = plotInit(XePosListList[0], [0], Nlist, timeList[-1])
    for i, t in enumerate(timeList[1:]):
        ePosList, XePosList = ePosListList[i+1], XePosListList[i+1]

        # Status bar
        sys.stdout.write('\r')
        perc = int(float(i+1)/len(timeList)*100)
        sys.stdout.write("[%-20s] %d%% (Ne = %d)" % ('='*int(0.2*perc), perc, len(ePosList)))
        sys.stdout.flush()

        # Plot
        updatePlot(fig, h, hBottom, hTop, hRight, t, i+1, np.array(timeList[:i+1])*1.e9, np.array([float(n) for n in Nlist]) / N * 100, ePosList, XePosList, show)
        Nlist.append( len(ePosList) )

def plotFigure(fig, ax, timeList, percList, xlabel=None, ylabel=None, label=None, title=None, color=None, lw=2, errList=None):
    timeList, percList = np.array(timeList), np.array(percList)

    # Mean + std
    if not isinstance(percList[0], float):
    # if len(percList[0]) == 2:
        percListMean, percListStd = percList[:,0], percList[:,1]

        ax.plot(timeList, percListMean, label=label, c=color, lw=lw)
        ax.fill_between(timeList, percListMean-percListStd, percListMean + percListStd, alpha=.5, color=color)
        
    else:
        ax.plot(timeList, percList, label=label, c=color, lw=lw)

    # Set axis label
    ax.set_xlabel( xlabel , fontsize=14)
    ax.set_ylabel( ylabel , fontsize=14)

    # Use exponential format on x-axis
    ax.get_xaxis().get_major_formatter().set_powerlimits((0, 0))
    # ax.legend(loc='best')

    if title:
        fig.suptitle( title )

    fig.canvas.draw()
    fig.show()

# === COLORS ===
def getColorBar(ax, cbMin, cbMax, N=20, label=None):
    # Plot colorbar
    from matplotlib.colors import ListedColormap
    # cmap = mpl.cm.get_cmap('viridis', N)
    cmap = ListedColormap(sns.cubehelix_palette(N, start=.2, rot=-.75, dark=.2, light=.7).as_hex()) # mpl.cm.get_cmap('Dark20')
    norm = mpl.colors.Normalize(vmin=cbMin, vmax=cbMax)
    cBar = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')

    # cBar.ax.invert_yaxis()
    cBar.formatter.set_powerlimits((0, 0))
    cBar.ax.yaxis.set_offset_position('right')
    cBar.update_ticks()

    print cbMin, cbMax, float(cbMax - cbMin)/N, N
    labels = np.linspace(cbMin, cbMax, N + 1)
    locLabels = np.linspace(cbMin, cbMax, N)
    print labels
    loc = labels + abs(labels[1] - labels[0])*.5
    cBar.set_ticks(loc)
    cBar.ax.set_yticklabels(locLabels * 1.e9, rotation=90, verticalalignment='center')
    cBar.outline.set_visible(False)
    print loc
    cBar.set_label(label)

def getColor(N, idx):
    # import colorcet as cc
    from matplotlib.colors import ListedColormap

    cmap = ListedColormap(sns.cubehelix_palette(N, start=.2, rot=-.75, dark=.2, light=.7).as_hex()) # mpl.cm.get_cmap('Dark20')
    # cmap = ListedColormap(sns.color_palette("PuBuGn_d").as_hex())
    norm = mpl.colors.Normalize(vmin=0.0, vmax=N - 1)
    return cmap(norm(idx))

# === SUPPORT ===
def filterQuery(dist, idx, length):
    elemList = np.array( [item for item in zip(idx, dist) if item[0] != length] )
    if elemList.size:
        return elemList[:,1], elemList[:,0].astype('int')
    else:
        return [], []

# === FIT ===
def parallelFit(x, A, a, b): 
    # return A / (1 + np.exp(a*x**b))
    # return A/(1 + a*np.sqrt(1./x**b))
    return A/(1 + a*x**(-b))

# === STOPPING POWER ===
# E in [keV]
def getDeltaE(deltaX, E):
    import scipy.integrate

    # Stopping power
    density = 3.057 # g/cm^3
    x, y = np.loadtxt('stopping_power.dat', delimiter='\t', usecols=range(2), unpack=True)
    x = np.array( x ) * 1.e3                # Convert MeV to keV
    y = np.array( y ) * density * 1.e3      # Multiply with density

    # Filter by energy
    x, y = zip(*[item for item in zip(x, y) if item[0] <= E])

    dTrack = scipy.integrate.simps(1./y, x, even='avg')

def get_args():
    ap = argparse.ArgumentParser(description=' ')

    ap.add_argument('-u', '--uniform', help='Use uniform initial distribution', action='store_true')
    ap.add_argument('-par', '--parallel', help='Simulation of perpendicular track', action='store_true')
    ap.add_argument('-sh', '--show', help='Show plots during simulation', action='store_true')
    ap.add_argument('-s', '--save', help='Save position data', action='store_true')
    ap.add_argument('-mp', '--makeplots', help="Make plots. Use 'save' before using this command", action='store_true')
    ap.add_argument('-t', '--test', help="Test the simulation", action='store_true')
    ap.add_argument('-pr', '--printresults', help='Show results', action='store_true')
    ap.add_argument('-si', '--single', help='Process only single parameter set', action='store_true')
    ap.add_argument('-A', '--val_A', help='Fuck', type=float, required=False)
    ap.add_argument('-b', '--val_b', help='Fuck', type=float, required=False)

    args = ap.parse_args()
    return args.uniform, args.parallel, args.show, args.save, args.makeplots, args.test, args.printresults, args.single, args.val_A, args.val_b

if __name__ == '__main__':
    main()

