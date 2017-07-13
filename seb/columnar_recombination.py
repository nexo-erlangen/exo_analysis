#!/usr/bin/env python
import argparse
import multiprocessing
import numpy as np
from numpy.random import normal, uniform
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
except:
    pass

# Constants
X = 3.8e6 # V/m^2, electric field
Dt = 55.e-4 # m^2/s
Dl = 0.1 * Dt
vD = 1705 # m/s
mobility = vD / X # mobility
b = 1.e-8 # 2.5e-3 # 5.e-6
d = 1.e-4 # 2.34e-5 # 2.34e-4 # 5.e-5    # Length of the column
A = 49.e-9  # Onsager radius
angle = 0

# Number of particles
N = 5000

# estimated number of survivors (for bottom plot)
N_rem = 1000

# Time step size
dt = 1.e-12 # 5.e-12  # s

# Status bar
STATUSBARWIDTH = 20

# How often to repeat a simulation in order to estimate
# a statistical error
REPEAT = 3

def main():
    # Load command line arguments
    uniform, parallel, show, save, MAKEPLOTS, TEST, PRINTRESULTS, SINGLE = get_args()
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

    global INTERRUPT
    if PARALLEL:
        INTERRUPT = 500      # Limit for interruption
        gifOut = 'cr/cr_gif_par'
        pickleOut = 'cr/cr_par'
        resultOut = 'cr/res_par'
    else:
        INTERRUPT = 100      # Limit for interruption
        gifOut = 'cr/cr_gif_perp'
        pickleOut = 'cr/cr_perp'
        resultOut = 'cr/res_perp'

    if not GAUSS:
        gifOut += '_uni'
        pickleOut += '_uni'
        resultOut += '_uni'

    pickleOut += '.p'
    resultOut += '.p'

    # For printresults
    resultOutList = ['cr/res_par.p', 'cr/res_par_uni.p', 'cr/res_perp.p', 'cr/res_perp_uni.p']
    titleList = ['Parallel (Gauss)', 'Parallel (uniform)', 'Perpendicular (Gauss)', 'Perpendicular (uniform)']
    colorList = ['orange', 'red', 'blue', 'green']

    # === MAKE PLOTS ===
    if MAKEPLOTS:
        dic = cPickle.load(open(pickleOut, 'rb'))
        plotDict(dic, SHOW)
        return

    # === PROCESS SINGLE ===
    if SINGLE:
        global b
        inter = 100
        generate(N, d, b, A, angle=None, interrupt=inter, out_q=None)
        return

    # === TEST ===
    if TEST:
        driftElectronTest()
        dt = 1.e-6
        driftElectronTest()
        raw_input('')
        return

    # === PRINT RESULTS ===
    if PRINTRESULTS:
        # Figure of recombination vs. b
        figRes, axRes = plt.subplots()

        resListList = []
        for m, resultOut in enumerate(resultOutList):
            if not os.path.isfile(resultOut):
                continue

            # Figure of recombination vs. t
            # fig, ax = plt.subplots()
            fig = plt.figure()
            ax = fig.add_axes([0.1, 0.1, 0.7, 0.8])
            # Add colorbar
            axCBar = fig.add_axes([0.85, 0.1, 0.05, 0.8])

            dic = cPickle.load(open(resultOut, 'rb'))
            bList, percListList, timeListList = dic['b'], dic['p'], dic['t']
            getColorBar(axCBar, bList[0], bList[-1], 'b [m]')

            resList = []
            for i, b in enumerate( bList ):
                timeList, percList = timeListList[i], percListList[i]
                color = getColor(len(bList), i)

                # Fit the curves in the parallel case
                if m in [0, 1]:
                    # Length of the column: d 
                    tt = d / vD * 1.e9

                    # Get only mean values
                    if not isinstance(percList[0], float):
                        pList = np.array(percList)[:,0]
                    else:
                        pList = percList

                    popt, pcov = curve_fit(parallelFit, timeList, pList)
                    perr = np.sqrt( np.diag(pcov) )
                    print popt, perr
            
                    timeListFit = np.linspace(0, tt, float(len(timeList))/timeList[-1]*tt)
                    percListFit = parallelFit(timeListFit, *popt)

                    resList.append( percListFit[-1] )

                    plotFigure(fig, ax, timeListFit, np.array(percListFit)*100, xlabel=r'time [ns]', ylabel=r'Recombination [%]', label=None, title=titleList[m], color=color)
                
                # PERPENDICULAR
                else:
                    resList.append( percList[-1] )

                # No label, use colorbar
                plotFigure(fig, ax, timeList, np.array(percList)*100, xlabel=r'time [ns]', ylabel=r'Recombination [%]', label=None, title=titleList[m], color=color)

            # Show diffusion loss
            if not isinstance(resList[0], float):
                resList = [(1-r[0], r[1]) for r in resList]
            else:
                resList = [1 - r for r in resList]

            
            # Fit result curve
            print len(bList), len(resList)
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

            # plotFigure(figRes, axRes, bListFit, np.array(resListFit)*100, xlabel=r'b [m]', ylabel=r'Diffusion loss [%]', color=colorList[m])
            plotFigure(figRes, axRes, bList, np.array(resList)*100, xlabel=r'b [m]', ylabel=r'Diffusion loss [%]', label=titleList[m], color=colorList[m])

        raw_input('')
        return

    getResults('b')
    return

# === GET RESULTS ===
# Loop over paraneters and store results to dict
def getResults(param='b'):
    # Run over bList
    # Create figure to use in plots
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.7, 0.8])
    # Add colorbar
    axCBar = fig.add_axes([0.85, 0.1, 0.05, 0.8])

    # bList = np.linspace(1.e-7, 30.e-7, 10)
    if param == 'b':
        paramList = np.linspace(1.e-8, 30.e-8, 10)
        colorBarLabel = 'b [m]'
    elif param == 'angle':
        paramList = np.linspace(0., np.pi, 5)
        colorBarLabel = 'Theta (rad)'

    #bList = [1.e-7, 2.e-7]
    # getColorBar(axCBar, bList[0], bList[-1], 'b [m]')
    getColorBar(axCBar, paramList[0], paramList[-1], 'b [m]')

    manager = multiprocessing.Manager()
    out_q = manager.Queue()

    procs = []
    # for b in bList:
    for para in paramList:
        if param == 'b':
            global angle
            b = para
        else:
            global b
            angle = para

        if REPEAT > 1:
            p = multiprocessing.Process(target=generateRepeat, args=(N, d, b, A, angle, INTERRUPT, out_q, REPEAT))
        else:
            p = multiprocessing.Process(target=generate, args=(N, d, b, A, angle, INTERRUPT, out_q))

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
        paramRes = resultDict[param]
        paramResList.append( paramRes )
        timeList, percList = resultDict['t'], resultDict['p']
        timeListList.append( timeList ), percListList.append( percList )

        # color = getColor(len(bList), list(bList).index(bRes))
        color = getColor(len(paramList), list(paramList).index(paramRes))
        plotFigure(fig, ax, timeList, np.array(percList)*100, xlabel=r'time [ns]', ylabel=r'Recombination [%]', color=color)

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
    resultDict = {param: paramResList, 'p': percListList, 't': timeListList}
    cPickle.dump(resultDict, open(resultOut, 'wb'))
    raw_input('')

def generateRepeat(N, d, b, A, angle=None, interrupt=False, out_q=None, repeat=1):
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

def generate(N, d, b, A, angle=None, interrupt=False, out_q=None):
    if PARALLEL:
        tt = d / vD
    else:
        tt = 2 * 5*b / vD

    # Generate Xe+ positions
    print 'Generate initial Xe+ positions...',
    XePosList = initialDistribution(N, d, b, angle, PARALLEL, GAUSS)
    print 'Done!'

    # Create cKD tree
    print 'Fill Xe+ KD tree...',
    xyz = np.c_[np.array(XePosList)[:,0], np.array(XePosList)[:,1], np.array(XePosList)[:,2]]
    XeTree = spatial.cKDTree(xyz)
    print 'Done!'

    meanDist = getMeanDistance(XePosList, XeTree)
    print 'Mean distance of two Xe+ ions:', meanDist
    print 'Minimum time step length:', meanDist / vD
    if dt > meanDist / vD:
        print 'dt is too small!'
        return

    # Generate electrons and check for intersection
    print 'Generate initial e positions...',
    ePosList = initialDistribution(N, d, b, angle, parallel=PARALLEL, gauss=GAUSS, tree=XeTree, A=A)
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

        # Loop over all electrons
        ePosListNew = []
        XeRemoveIdx = []
        # s = time.clock()
        for ePos in ePosList:
            # Drift the electron
            ePosNew = driftElectron(ePos, [Dt, Dt, Dl], vD, dt, eTree, XeTree, ePosList, XePosList)

            # Check for intersection only if electron is close
            # enough to the column
            r = getRadius(ePosNew, -vD*t)
            if r > 3*b or (PARALLEL and (ePosNew[2] < -.6*d)):
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
                
        # print time.clock() - s

        ePosList = ePosListNew
        # Renew the XeTree
        XePosList = list([XePos for u, XePos in enumerate(XePosList) if u not in XeRemoveIdx])

        # Save only every 100th frame
        # if not i % 10:
        if SAVE:
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
                # slope = np.mean([dt/(percList[n+1]-percList[n]) if percList[n] != percList[n+1] else 0 for n in range(len(percList)-100, len(percList), 2)]) / 100
                mean = .5*(NList[-interrupt] + NList[-2])
                # if abs(slope) < 1.e-13:
                if mean == NList[-1]:
                    print
                    print 'In region of small changes -> interrupting'
                    break

        if SHOW:
            updatePlot(fig, h, hBottom, hTop, hRight, t, i, timeList, NList, ePosList, XePosList, show=SHOW)

    # After the time loop is done, count remaining particles 
    print
    print 'Calculation time:', time.clock() - loopStart 
    print 'Electrons remaining:', len(ePosList)
    print 'Recombination percentage:', 1 - float(len(ePosList))/N

    if SAVE:
        # Save Xe+ & e positions to disk
        cPickle.dump(posDict, open(pickleOut, 'wb'))

    if out_q:
        print
        print 'Reached end'
        outDict = {'b': b, 't': timeList, 'p': percList}
        out_q.put( outDict )

    return timeList, percList

def driftElectronTest():
    pos = (0, 0, 0)

    x, dx = [], []
    y, dy = [], []
    z, dz = [], []
    time = np.arange(0, tt, dt)
    for i in time:
        posNew = driftElectron(pos, [Dt, Dt, Dl], vD, dt)
        x.append( posNew[0] ), dx.append( posNew[0] - pos[0] )
        y.append( posNew[1] ), dy.append( posNew[1] - pos[1] )
        z.append( posNew[2] ), dz.append( posNew[2] - pos[2] )
        pos = posNew

    print 'Mean:', np.mean(x)
    print 'Std:', np.std(x)
    print 'Std (theory):', np.sqrt(2*Dt*tt)

    f, ax = plt.subplots()

    ax.plot(time, x)
    ax.plot(time, y)
    ax.plot(time, z)
    # plt.plot(time, dx)
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

def checkIntersect(ePos, XeTree, A, k=1, ePosOld=None, XePosList=None):
    ePos = np.array(ePos)

    # Point - sphere intersection
    if not ePosOld:
        searchDist = np.inf

        dist, idx = XeTree.query(np.array([ePos]), k=k, n_jobs=0)
        dist, idx = dist[0], idx[0]
        if not isinstance(dist, float):
            dist, idx = float(dist[-1]), int(idx[-1])

        if dist <= A:
            return True, idx
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
    newPos = []
    for i, p in enumerate(pos):
        sigma = np.sqrt( 2*D[i]*dt )
        newPos.append( p + np.random.normal(0, sigma, 1)[0] )

    # Calculate couloumb interactions
    if eTree and XeTree and ePosList and XePosList:
        from scipy import constants
        epsilon = 1.95
        coulombRange = 1.e-8

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
        newPos[2] -= vD * dt

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

    # Number of electrons vs. time
    hBottom = axBottom.plot(timeList, NList)

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
    ax.set_xlabel('x'), ax.set_ylabel('y'), ax.set_zlabel('z')
    axBottom.set_xlabel(r't [ns]')
    axBottom.set_ylabel(r'N')
    axRight.set_ylabel(r'z')

    if PARALLEL:
        axTop.set_xlabel('x')
        axTop.set_ylabel('y')
    else:
        axTop.set_xlabel('z')
        axTop.set_ylabel('y')

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

    # Right
    if PARALLEL:
        axRight.set_ylim(-d, d)
    else:
        axRight.set_ylim(-5*b, 5*b)

    # Bottom
    axBottom.set_xlim(0., tt*1.e9) # in [ns]
    axBottom.set_yscale('symlog', nonposx='clip')
    axBottom.set_ylim(N_rem, N*1.1) 

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
        updatePlot(fig, h, hBottom, hTop, hRight, t, i+1, np.array(timeList[:i+1])*1.e9, Nlist, ePosList, XePosList, show)
        Nlist.append( len(ePosList) )

def plotFigure(fig, ax, timeList, percList, xlabel=None, ylabel=None, label=None, title=None, color=None):
    timeList, percList = np.array(timeList), np.array(percList)

    # Mean + std
    if not isinstance(percList[0], float):
    # if len(percList[0]) == 2:
        percListMean, percListStd = percList[:,0], percList[:,1]

        ax.plot(timeList, percListMean, label=label, c=color)
        ax.fill_between(timeList, percListMean-percListStd, percListMean + percListStd, alpha=.5, color=color)
        
    else:
        ax.plot(timeList, percList, label=label, c=color)

    # Set axis label
    ax.set_xlabel( xlabel )
    ax.set_ylabel( ylabel )

    # Use exponential format on x-axis
    ax.get_xaxis().get_major_formatter().set_powerlimits((0, 0))
    ax.legend(loc='best')

    if title:
        fig.suptitle( title )

    fig.canvas.draw()
    fig.show()

# === COLORS ===
def getColorBar(ax, cbMin, cbMax, label=None):
    # Plot colorbar
    cmap = mpl.cm.get_cmap('viridis')
    norm = mpl.colors.Normalize(vmin=cbMin, vmax=cbMax)
    cBar = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')

    cBar.ax.invert_yaxis()
    cBar.formatter.set_powerlimits((0, 0))
    cBar.ax.yaxis.set_offset_position('right')
    cBar.update_ticks()

    cBar.set_label(label)

def getColor(N, idx):
    cmap = mpl.cm.get_cmap('viridis')
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
    return A / (1 + np.exp(a*x**b))
    # return A/(1 + a*np.sqrt(1./x**b))

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

    args = ap.parse_args()
    return args.uniform, args.parallel, args.show, args.save, args.makeplots, args.test, args.printresults, args.single

if __name__ == '__main__':
    main()

