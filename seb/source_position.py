#!/usr/bin/env python
import numpy as np
import os
import ROOT
from matplotlib.backends.backend_pdf import PdfPages

import plot_support as ps

# Fiducial volume
# FV = [171, 5, 182]
# Apothem = 171
FV = [171, 5, 182]
Apothem = 171
REFLECTORINNERRAD = 183.2356 # mm

# Paths
preProcessDir = '/home/vault/capm/mppi025h/analysis/preprocess/'

def main():
    # folderList = ['S17Th', 'S17ThShift', 'S11Th', 'S11ThShift']
    topList = [False, False, True, True]
    colorList = ['red', 'blue', 'blue', 'red']
    lsList = ['-', '--', '-', '--']
    muList = [-25., 100., 50., -90.]

    # MC only
    folderList = [None, None, None, None]
    fileList = ['SourceS17_Th228.root', 'SourceS17Shift_Th228.root', 'SourceS11_Th228.root', 'SourceS11Shift_Th228.root']
    # compareSourcePositions(folderList, fileList, topList, colorList, lsList, muList)

    # Data only
    folderList = ['S17Th', 'S17ThShift', 'S11Th', 'S11ThShift']
    fileList = [None]*4
    # compareSourcePositions(folderList, fileList, topList, colorList, lsList, muList)

    # MC Shift vs MC
    folderList = [None]*8
    fileList = ['SourceS17_Th228.root', 'SourceS17Shift_Th228.root', 'SourceS11_Th228.root', 'SourceS11Shift_Th228.root', 'SourceS17Art_Th228.root', 'SourceS17ShiftArt_Th228.root', 'SourceS11Art_Th228.root', 'SourceS11ShiftArt_Th228.root']

    topList = [False, False, True, True, False, False, True, True]
    colorList = ['red', 'blue', 'blue', 'red', 'orange', 'darkviolet', 'darkviolet', 'orange']
    lsList = ['-', '--', '-', '--', '-', '--', '-', '--']
    muList = [-25., 100., 50., -90.]*2

    # compareSourcePositions(folderList, fileList, topList, colorList, lsList, muList)

    # MC & Data
    folderList = ['S17Th', 'S17ThShift', 'S11Th', 'S11ThShift'] + [None]*4
    fileList = [None]*4 + ['SourceS17_Th228.root', 'SourceS17Shift_Th228.root', 'SourceS11_Th228.root', 'SourceS11Shift_Th228.root']
    topList = [False, False, True, True]*2
    colorList = ['red', 'blue', 'blue', 'red', 'orange', 'darkviolet', 'darkviolet', 'orange']
    lsList = ['-', '--', '-', '--']*2
    muList = [-25., 100., 50., -90.]*2

    compareSourcePositions(folderList, fileList, topList, colorList, lsList, muList, 'source_position_MCvsData.pdf')

    # MCArt & Data
    folderList = ['S17Th', 'S17ThShift', 'S11Th', 'S11ThShift'] + [None]*4
    fileList = [None]*4 + ['SourceS17Art_Th228.root', 'SourceS17ShiftArt_Th228.root', 'SourceS11Art_Th228.root', 'SourceS11ShiftArt_Th228.root']
    topList = [False, False, True, True]*2
    colorList = ['red', 'blue', 'blue', 'red', 'orange', 'darkviolet', 'darkviolet', 'orange']
    lsList = ['-', '--', '-', '--']*2
    muList = [-25., 100., 50., -90.]*2

    compareSourcePositions(folderList, fileList, topList, colorList, lsList, muList, 'source_position_MCArtvsData.pdf')
    raw_input('')

def compareSourcePositions(folderList, fileList, topList, colorList, lsList, muList=None, pdfOut=None):
    xHistList, xHistTheoList = [], []
    zHistList, zHistTheoList = [], []
    poptList, perrList = [], []

    rHistList, rHistErrList = [], []
    rCumSumList, rCumSumErrList = [], []

    # Number of events in data runs to get the same number in MC runs
    dataEvents = []
    k = 0
    for i, folder in enumerate( folderList ):
        if folder:
            dataList = getDataList( folder )
            dataChain = ps.getChain('dataTree', dataList)

            # Get z positions
            rList, xList, yList, zList, NEvents = getZDistribution(dataChain, art='ss', thetaCut=False, top=topList[i], MC=False, NEvents=None)
            dataEvents.append( NEvents )

        else:
            dataChain = ROOT.TFile.Open(preProcessDir + fileList[i]).Get('mcTree')
            if dataEvents:
                NEvents = dataEvents[k]
            else:
                NEvents = None

            rList, xList, yList, zList, NEvents = getZDistribution(dataChain, art='ss', thetaCut=False, top=topList[i], MC=True, NEvents=NEvents)
            k += 1

        pInit = [0., 50., 0.01, 0.001]
        x, xHist, xTheo, xHistTheo, poptX, perrX = getFit(xList, pInit)

        if muList:
            pInit = [muList[i]]
        else:
            if topList[i]:
                pInit = [100.]
            else: 
                pInit = [-100.]

        pInit += [50., 0.01, 0.001]
        z, zHist, zTheo, zHistTheo, popt, perr = getFit(zList, pInit)
        
        xHistList.append( xHist ), zHistList.append( zHist )
        xHistTheoList.append( xHistTheo ), zHistTheoList.append( zHistTheo )
        poptList.append( popt )
        perrList.append( perrList )

        # Get radii around mean
        print 'Source%s: (%.3f, y, %.3f)' % (folderList[i], poptX[0],  popt[0])
        mu = popt[0]
        if mu > 0:
            mu += 0
            # mu = 100
        else: 
            mu -= 0
            # mu = -100
            
        r, rHist, rHistErr, rCumSum, rCumSumErr = getRadiusDist(zList, rList, mu, 5.)

        rCum = list(reversed(r))
        rCum = list( np.array(rCum) + abs(rCum[1] - rCum[0]) )
        print len(rCumSum), len(rHist), len(rCum), len(r)

        rHistErrList.append( rHistErr )
        rHistList.append( rHist )

        # Get cumulative sum
        # As the standoff distance is shown, begin the
        # sum from high to low values
        # rCumSum = [0] + list( np.cumsum( list(reversed(rHist[:-1])) )  )
        # rCumSumErr = [0] + list(np.sqrt( np.cumsum( np.square(list(reversed(rHistErr[:-1]))) ) ))
        # print rCumSumErr

        rCumSumList.append( rCumSum )
        rCumSumErrList.append( rCumSumErr )

    # Residuals for standoff distance
    # Compare the source positions to each other
    if len( folderList ) == 4:
        resIdxList = [0]*3
        idxRange = range( len(folderList) )
        idxRange.remove( 0 )

    # Compare two data sets
    elif len( folderList ) == 8:
        resIdxList = range(4)
        idxRange = range(4, 8)

    resList, resListErr = [], []
    resCumList, resCumErrList = [], []
    for j, i in enumerate( idxRange ):
        resIdx = resIdxList[j]

        normEntry = np.array( rHistList[resIdx] )
        normEntryErr = np.array( rHistErrList[resIdx] )

        res = (normEntry - np.array(rHistList[i]))/normEntry
        resErr = np.sqrt( (1./normEntry * np.array( rHistErrList[i] ))**2 + (np.array( rHistList[i] )/(normEntry)**2 * normEntryErr)**2 )
        resErr = [0 if np.isneginf(item) or np.isinf(item) or np.isnan(item) else item for item in resErr]

        normEntry = np.array( rCumSumList[resIdx] )
        normEntryErr = np.array( rCumSumErrList[resIdx] )

        resCum = (normEntry - np.array(rCumSumList[i]))/normEntry
        resCumErr = np.sqrt( (1./normEntry * np.array( rCumSumErrList[i] ))**2 + (np.array( rCumSumList[i] )/(normEntry)**2 * normEntryErr)**2 )
        resCumErr = [0 if np.isneginf(item) or np.isinf(item) or np.isnan(item) else item for item in resCumErr]
    
        print 'Residuals'
        print (resIdx, i)
        print (rHistList[resIdx], rHistList[i])
        print (rHistErrList[resIdx], rHistErrList[i])
        print (res, resErr), (resCum, resCumErr)

        resList.append( res ), resListErr.append( resErr )
        resCumList.append( resCum ), resCumErrList.append( resCumErr )

        # resList = [(normEntry - np.array(h))/np.array(normEntry) for h in (rHistList[:resIdx] + rHistList[resIdx+1:])]
        # resListErr = [abs(np.array(h)/normEntry) * np.sqrt( (rHistErrList[i+1]/normEntry)**2 + (normEntryErr/np.array(h))**2 ) for i, h in enumerate(rHistList[:resIdx] + rHistList[resIdx+1:])]

    # Residuals for cumulative sum
    '''
    normEntry = np.array( rCumSumList[resIdx] )
    normEntryErr = np.array( rCumSumErrList[resIdx] )
    resCumList = [(normEntry - np.array(h))/np.array(normEntry) for h in (rCumSumList[:resIdx] + rCumSumList[resIdx+1:])]
    resCumErrList = [abs(np.array(h)/normEntry) * np.sqrt( (rCumSumErrList[i+1]/normEntry)**2 + (normEntryErr/np.array(h))**2 ) for i, h in enumerate(rCumSumList[:resIdx] + rCumSumList[resIdx+1:])]
    '''
    print poptList

    # == Plot ==
    import re
    labelList = [l if l else re.sub(r'.*S1', 'S1', fileList[m].split('_')[0]) for m, l in enumerate(folderList)]
    if len(resList) == 3:
        labelListDiff = ['(S17-S17S)/S17', '(S17-S11)/S17', '(S17-S11S)/S17']
    elif len(resList) == 4:
        labelListDiff = ['(S17-S17M)/S17', '(S17S-S17SM)/S17S', '(S11-S11M)/S11', '(S11S - S11SM)/S11S']
    elif len(resList) == 7:
        labelListDiff = ['(S17-S17S)/S17', '(S17-S11)/S17', '(S17-S11S)/S17', '(S17-S17M)/S17', '(S17-S17SM)/S17', '(S17-S11M)/S17', '(S17-S11SM)/S17']

    # Colors for fits
    # lists are doubled to give fits the same color as data
    if len(xHistList) == 4:
        # fitColor = ['#6dac5a', '#895eb2', '#c18a3b', '#be485c']*2
        fitColor = ['#64b8e8', '#e89064', '#e8d264', '#be64e8']*2
    elif len(xHistList) == 8:
        # fitColor = ['#6dac5a', '#895eb2', '#c18a3b', '#be485c', '#3A7927', '#562B7F', '#8E5708', '#8B1529']*2
        fitColor = ['#64b8e8', '#e89064', '#e8d264', '#be64e8', '#3185B5', '#B55D31', '#B59F31', '#8B31B5']*2
    
    # Colors for residuals
    if len(resList) == 3:
        customColor = ['#a26bbd', '#7ea24e', '#cb6450']
    elif len(resList) == 4:
        customColor = ['#a26bbd', '#7ea24e', '#cb6450', '#da680b']
    elif len(resList) == 7:
        # 3 colors + 3 darkened colors + 1 extra color
        customColor = ['#a26bbd', '#7ea24e', '#cb6450', '#da680b', '#6F388A', '#4B6F1B', '#98311D']

        # 7 different colors
        # customColor = ['#69a742', '#c05db0', '#4eaa86', '#ca673a', '#ad943f', '#ce546a', '#737cce']

    meanLabel = ['%s: %.2f mm' % (folderList[i] if folderList[i] else fileList[i].split('_')[0], np.array(poptList)[:,0][i]) for i in range(len(folderList))]

    # = PLOT =
    if pdfOut:
        pp = PdfPages( pdfOut )
    else:
        pp = None

    # Gauss fits to find the mean x-positions
    plotLists([x]*len(xHistList) + [xTheo]*len(xHistTheoList), xHistList + xHistTheoList, labelList+[None]*len(labelList), xlabel='x [mm]', title='x-Distribution', colorList=fitColor, pp=pp)

    # Gauss fits to find the mean z-positions
    plotLists([z]*len(zHistList) + [zTheo]*len(zHistTheoList), zHistList + zHistTheoList, labelList+[None]*len(labelList), title='z-Distribution', colorList=fitColor, pp=pp)

    # Standoff distributions
    plotLists([r]*len(rHistList), rHistList, labelList, log=False, xlabel='Standoff distance [mm]', title='Standoff distribution', colorList=fitColor, pp=pp)

    # Standoff distributions w/ different color scheme
    # and error bars
    plotLists([r]*len(rHistList), rHistList, meanLabel, xlabel='Standoff distance [mm]', errList=rHistErrList, lsList=lsList, colorList=colorList, title='Standoff distribution', pp=pp)

    # Residuals of standoff distributions
    plotLists([r]*len(resList), resList, labelListDiff, log=False, xlabel='Standoff distance [mm]', errList=resListErr, colorList=customColor, title='Standoff residuals', pp=pp)

    # Cumulative sum
    plotLists([r]*len(rCumSumList), rCumSumList, meanLabel, log=False, xlabel='Standoff distance [mm]', errList=rCumSumErrList, lsList=lsList, colorList=fitColor, title='Cumulative sum', pp=pp)

    # Cumulative sum w/ different color scheme
    plotLists([r]*len(rCumSumList), rCumSumList, meanLabel, log=True, xlabel='Standoff distance [mm]', errList=rCumSumErrList, lsList=lsList, colorList=colorList, title='Cumulative sum', pp=pp)

    # Residuals of cumulative sum
    plotLists([r]*len(resCumList), resCumList, labelListDiff, log=False, xlabel='Standoff distance [mm]', errList=resCumErrList, colorList=customColor, title='Cumulative sum residuals', pp=pp)

    if pp:
        pp.close()

def getDataList(folder):
    dataList = []
    for fn in os.listdir(preProcessDir + folder + '/'):
            dataList.append( preProcessDir + folder + '/%s' % fn )

    return dataList

def getFit(zList, pInit):
        # Histogramize
        zList = [z for z in zList if (abs(z) > 20)]
        zHist, z = np.histogram(zList, 80, (-160, 160), density=True)
        zHist = np.nan_to_num( zHist )
        zS, zHist = removeRange(z[:-1], zHist, [-20, 20])
        z = zS

        # Fit
        zTheo = np.linspace(-160, 160, 1000)

        popt, perr = fitGauss(z, zHist, pInit)
        zHistTheo = gauss(zTheo, *popt)

        return z, zHist, zTheo, zHistTheo, popt, perr

def getZDistribution(tree, art='ss', thetaCut=False, top=True, MC=False, NEvents=None):
	ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
	ROOT.gROOT.cd()
    
	cut = ps.getCut(calibCut=True, energyCut=True, type=art, MC=MC, eMin=750, eMax=3500)
        cut = cut.replace('&& !isSolicitedTrigger', '')

        # cut = 'multiplicity > 0 && !isMissingPosition && nsc == 1 && (Min$(abs(cluster_z)) > 5 && Max$(abs(cluster_z)) < 182 && Max$(abs(cluster_x)*(abs(cluster_x)<900)) < %d && Max$(abs( 0.5*( ((cluster_z>0) - (cluster_z<=0))*cluster_x*(abs(cluster_x)<900) + sqrt(3.0)*cluster_y*(abs(cluster_y)<900))  )) < %d && Max$(abs( 0.5*( ((cluster_z<=0) - (cluster_z>0))*cluster_x*(abs(cluster_x)<900) + sqrt(3.0)*cluster_y*(abs(cluster_y)<900))  )) < %d && Max$(cluster_x*cluster_x*(abs(cluster_x)<900) + cluster_y*cluster_y*(abs(cluster_y)<900)) < 33575.3) && !isNoise && isWithinDriftTime && !EventSummary.isDiagonallyCut() && !isVetoed && (energy > 750. && energy < 3500.) &&' % (FV[0], FV[0], FV[0])
        
        # isSolicitedTrigger makes trouble
        if art == 'ss':
            cut += ' && abs(multiplicity - 1) < 0.1'
        else:
            cut += ' && multiplicity > 1'
        print cut

	ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
	treeCut = tree.CopyTree( cut )

        if not NEvents:
            NEvents = treeCut.GetEntries()

        xList, yList, zList = [], [], []
        rList = []
        for j in range( NEvents ):
		treeCut.GetEntry(j)
		es = treeCut.EventSummary

		x, y, z = np.array(es.cluster_x), np.array(es.cluster_y), np.array(es.cluster_z)
                if thetaCut:
                    if top:
                        theta = np.arctan2(x, y)
                    else:
                        theta = np.arctan2(x, -y)
                    thetaCutVal = 30.*np.pi/180 - np.arccos(float(Apothem)/REFLECTORINNERRAD)
                    thetaOffset = 0 # -5.*np.pi/180
                    if not (theta <= (thetaCutVal+thetaOffset) and theta >= (-thetaCutVal+thetaOffset)):
                        continue

                rList += list( REFLECTORINNERRAD - np.sqrt( x**2 + y**2 ) )
                zList += list( z )
                xList += list( x )
                yList += list( y )

        return rList, xList, yList, zList, NEvents

def removeRange(x, y, ran):
    xFilt, yFilt = [], []
    for i in range( len(x) ):
        if not (x[i] <= ran[1] and x[i] >= ran[0]):
            xFilt.append( x[i] )
            yFilt.append( y[i] )

    return xFilt, yFilt

def gauss(x, mu, sigma, A, C):
    return A*np.exp(-(x - mu)**2/(2*sigma**2)) + C

def fitGauss(x, y, pInit):
    import scipy.optimize
    import scipy.stats

    try:
        popt, pcov = scipy.optimize.curve_fit(gauss, x, y, p0=pInit)
        perr = np.sqrt(np.diag(pcov))
    except:
        popt, perr = pInit, [0.]*len(pInit)

    return popt, perr

def plotLists(xList, yList, labelList, log=True, xlabel='z [mm]', title='', errList=None, lsList=None, colorList=None, pp=None):
    from matplotlib import pyplot as plt
    import matplotlib.ticker as ticker

    if not colorList:
        colorList = len( yList ) * [None]
    if not lsList:
        lsList = len( yList ) * ['-']

    '''
    if len( labelList ) != len( yList ):
        if 2*len( labelList ) == len( yList ):
            labelList = labelList + [label + 'MC' for label in labelList]
            # labelList = labelList + (len(lsList) - len(labelList))*[None]
        else:
            labelList = labelList + (len(yList) - len(labelList))*[None]
    if len( colorList ) != len( yList ):
        colorList = colorList + (len(yList)-len(colorList))*[None]
    if len( lsList ) != len( yList ):
        lsList = lsList + (len(yList)-len(lsList))*[None]

    '''
    print labelList
    print lsList
    print colorList

    f, ax = plt.subplots()
    for i, y in enumerate( yList ):
        bw = xList[i][1] - xList[i][0]
        if not errList:
            ax.plot(np.array( xList[i] ) + .5*bw, y, label=labelList[i], color=colorList[i], ls=lsList[i])
        else:
            ax.errorbar(np.array( xList[i] ) + .5*bw, y, yerr=errList[i], fmt='x', color=colorList[i], label=labelList[i], ls=lsList[i])
        # ax.step(xList[i], y, where='post', label=labelList[i], c=colorList[i], ls=lsList[i])

    ax.set_xlabel( xlabel )
    ax.set_ylabel( 'Normalized counts' )
    f.suptitle( title )

    leg = ax.legend()
    # set the linewidth of each legend object
    for legobj in leg.legendHandles:
            legobj.set_linewidth(2.0)

    if log:
        ax.set_yscale("log", nonposy='clip')
    else:
        ax.axhline(y=0., linewidth=.5, ls='--')

    if pp:
        f.savefig(pp, format='pdf')
    else:
        f.show()

def getRadiusDist(z, r, mu, ran):
    z, r = removeRange(z, r, [-200, mu-ran])
    z, r = removeRange(z, r, [mu+ran, 200])
    
    rHist, r = np.histogram(r, 200, (0, 200), density=False)
    rHistErr = np.sqrt( np.array(rHist).astype(np.float32) ) / sum( rHist )
    rHist = np.array(rHist).astype(np.float32) / sum( rHist )

    # Get only specific range
    rS, rHist = removeRange(r[:-1], rHist, [20, 200])
    rS, rHistErr = removeRange(r[:-1], rHistErr, [20, 200])
    # r = r[:-1]

    # Get cumulative sum
    # As the standoff distance is shown, begin the
    # sum from high to low values
    rCumSum = list( reversed( list( np.cumsum( list(reversed(rHist)) )  ) ) )
    rCumSumErr = list( reversed( list(np.sqrt( np.cumsum( np.square(list(reversed(rHistErr))) ) )) ) )

    rS, rCumSum = removeRange(r[:-1], rCumSum, [20, 200])
    rS, rCumSumErr = removeRange(r[:-1], rCumSumErr, [20, 200])
    r = rS
    return r, rHist, rHistErr, rCumSum, rCumSumErr

if __name__ == '__main__':
	main()

