import numpy as np
import plot_support as ps

import ROOT
from matplotlib.backends.backend_pdf import PdfPages
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')

pdfOut = 'energy_smear_nodiff2.pdf'
stdResolution = [1.81960744, 0.00306382, 0.03331489]

PHASE2 = False

def energySmear(tData, dataList, tMC, FV, art='ss'):
    import re
    # Return fraction of SS and MS events and compare MC and data
    # SSfrac(dataList, tData, tMC, FV)
    # raw_input('')
    # return

    bins = [100, 20]
    hRange = [[2000, 3000], [0, 200]]

    ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
    ROOT.gROOT.cd()

    pp = PdfPages( pdfOut )

    # === DATA ===
    # Cut trees and get histograms
    print 'Processing data...'
    cutData = ps.getCut(calibCut=True, energyCut=True, type=art, MC=False, eMin=hRange[0][0] - 100, eMax=hRange[0][1] + 100)
    # Gets rid of pileup in double escape peak
    cutData += ' && e_scint > 0'
    if PHASE2:
        cutData = re.sub('!EventSummary.isDiagonallyCut\(\) && ', '', cutData)
    # cutData += ' && ( (energy_%s > 1000 && energy_%s < 1200) || (energy_%s > 1561 && energy_%s < 1673 ) || (energy_%s > 1800 && energy_%s < 2000) )' % (art, art, art, art, art, art)
    
    print cutData
    tDataCut = tData.CopyTree( cutData )
    hData, hDataE, hDataSO, xData, yData = fillEnergy(tDataCut, None, bins, hRange, art, False)
    
    # Get resolution of data by fits
    # resolution = getResolution(xData, hDataE, src='charge', art=art)
    resolution = [-.000123185451, 74.7268349, 0.0270834510]

    # Standoff vs. energy
    # plotEnergy(hData, xData, yData)
    plotEnergyThree(hData, xData, yData, False, pp)

    # === MC ===
    print 'Processing MC...'
    cutMC = ps.getCut(calibCut=True, energyCut=True, type=art, MC=True, eMin=hRange[0][0] - 100, eMax=hRange[0][1] + 100)
    # cutMC += ' && ( (energy_%s > 1000 && energy_%s < 1200) || (energy_%s > 1561 && energy_%s < 1673 ) || (energy_%s > 1800 && energy_%s < 2000) )' % ('mc', 'mc', 'mc', 'mc', 'mc', 'mc')
    cutMC = re.sub('!EventSummary.isDiagonallyCut\(\) && ', '', cutMC)
    print cutMC
    tMCCut = tMC.CopyTree( cutMC )

    # Method using one resolution
    hMC, hMCE, hMCSO, xMC, yMC = fillEnergy(tMCCut, resolution, bins, hRange, art, True)
    # Method using multiple resolutions
    # hMC, hMCE, hMCSO, xMC, yMC = smearMCRandom(dataList, tMCCut, bins, hRange, art)

    # === PLOTS ===
    plotEnergyStandoff(hData, hMC, xData, bins, hRange, pp)
    raw_input('')

    # Standoff distance vs. energy
    # plotEnergy(hMC, xMC, yMC)
    plotEnergyThree(hMC, xMC, yMC, False, pp)

    # === DATA vs. MC ===
    # Residual energy
    plotEnergyRes(hDataE, hMCE, xData, None, pp)

    # Normalization
    # normRange = [int(np.digitize(1750, xData)), int(np.digitize(1950, xData))]
    normRange = [int(np.digitize(1000, xData)), int(np.digitize(2000, xData))]

    hDataNorm = getNorm( hData, hDataE, normRange )
    hMCNorm = getNorm( hMC, hMCE, normRange )

    # Residual standoff vs. energy
    plotEnergy((hDataNorm - hMCNorm)/hDataNorm, xData, yData, True, pp)

    pp.close()
    raw_input('')

# === SS FRAC ===
def SSfrac(dataList, tData, tMC, FV, art='energy'):
    resolution = {}
    resolution['ss'] = [0.727855, 17.7408, 0.]
    resolution['ms'] = [0.714076, 21.9255, 0.]

    # Histogram settings
    if art == 'standoff':
        bins = 20
        hRange = [0, 200]
    elif art == 'energy':
        bins = 50
        hRange = [750, 3500]
    else:
        return False

    # Set ROOT to operate in memory
    ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
    ROOT.gROOT.cd()

    # Data
    cutData = ps.getCutCombined(calibCut=True, energyCut=True, MC=False, eMin=700, eMax=3500)
    # Gets rid of pileup in double escape peak
    cutData += ' && e_scint > 0'
    print cutData
    tDataCut = tData.CopyTree( cutData )
    print

    # MC
    cutMC = ps.getCutCombined(calibCut=True, energyCut=True, MC=True, eMin=700, eMax=3500)
    print cutMC
    tMCCut = tMC.CopyTree( cutMC )
    print

    x, DataFrac, DataFracErr = fillSSFrac(tDataCut, bins, hRange, art, False)
    x, MCFrac, MCFracErr = fillSSFracMC(dataList, tMCCut, bins, hRange, art)

    plotSSfrac(x, DataFrac, DataFracErr, MCFrac, MCFracErr, art)
    # raw_input('')

# == PLOT FUNCTIONS ==
# Use not normalized histograms!
def plotEnergyRes(hData, hMC, x, title=None, pp=None):
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib import rc

    # = Error and Normalization =
    x = x[:-1]
    # Data
    hData, hDataErr = getErrorNorm( hData )

    # MC
    hMC, hMCErr = getErrorNorm( hMC )

    # Init figure and axes
    f, (axData, axRes) = plt.subplots(2, sharex=True, sharey=False, gridspec_kw = {'height_ratios':[2.5, 1]})
    f.subplots_adjust(wspace=0, hspace=0.1)

    # Axis 1 - Data and MC
    axData.set_yscale('log', nonposy='clip')

    binW = getBinWidth(x)
    axData.step(x, hMC, where='post', color='blue', label='MC')
    axData.errorbar(x+binW*0.5, hMC, yerr=hMCErr, fmt='none', color='blue')
    axData.errorbar(x+binW*0.5, hData, yerr=hDataErr, fmt='.', color='black', label='Data')

    # Axis 2 - residuals
    res = (hData - hMC) / hData
    resErr = abs(hMC/hData) * np.sqrt( (hMCErr/hMC)**2 + (hDataErr/hData)**2 )
    axRes.errorbar(x+binW*0.5, res, xerr=binW*.5, yerr=resErr, fmt='.', color='black')
    axRes.axhline(y=0, color='k', linestyle='--')
    axRes.set_ylim(-0.5, 0.5)

    # Labels
    # ax1.set_xlim(0, 180)
    # ax1.set_ylim(0, 0.25)
    # ax1.xaxis.set_ticks(np.arange(0, 200+20, 20))
    # ax1.yaxis.set_ticks(np.arange(0, 0.25+0.05, 0.05))
    axData.legend(loc='best')
    axData.set_ylabel(r'Normalized counts / (%.1f mm)' % binW)
    axRes.set_xlabel(r'Energy [keV]')
    axRes.set_ylabel(r'$\frac{n_\text{Data} - n_\text{MC}}{n_\text{Data}}$')

    if title:
        f.suptitle(title)

    if pp:
        f.savefig(pp, format='pdf')

    f.show()

def plotEnergy(h, x, y, diff=False, pp=None):
    import matplotlib.pyplot as plt
    from matplotlib import cm
    # binWidth = float(hRange[1][1]-hRange[1][0])/bins[1]
    # plt.bar(yData[:-1], hDataSO, width=binWidth)

    # Create figure
    f, ax = plt.subplots(1)

    if diff:
        im = ax.imshow(h, interpolation='nearest', origin='low', extent=[x[1], x[-1], y[0], y[-1]], aspect='auto', cmap='RdBu', vmin=-0.4, vmax=0.4)
    else:
        im = ax.imshow(h, interpolation='nearest', origin='low', extent=[x[1], x[-1], y[0], y[-1]], aspect='auto', cmap='inferno')

    plt.xlabel('Energy [keV]')
    plt.ylabel('Standoff Distance [mm]')
    cbar = plt.colorbar(im)
    cbar.ax.set_ylabel('Normalized Counts')

    if pp:
        f.savefig(pp, format='pdf')

    f.show()
    # raw_input('')

def plotEnergyThree(h, x, y, diff=False, pp=None):
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib import rc
    from matplotlib.ticker import NullFormatter
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    nullfmt = NullFormatter()

    # Normalization
    N = float( np.sum( h ) )
    h = h / N

    # Axis definition
    left, width = 0.1, 0.6
    bottom, height = 0.1, 0.35
    left_h = left + width + 0.05
    bottom_h= bottom + height + 0.05

    # Plot areas
    rect_main = [left, bottom, width, height]
    rect_histE = [left, bottom_h, width, 0.2]
    rect_histSo = [left_h, bottom, 0.15, height]

    f = plt.figure(figsize=(8, 8))
    axMain = plt.axes(rect_main)
    axHistE = plt.axes(rect_histE)
    axHistSo = plt.axes(rect_histSo)

    # No axis labels
    axHistE.xaxis.set_major_formatter(nullfmt)
    axHistSo.yaxis.set_major_formatter(nullfmt)
    # Exponential formatting
    # axHistE.yaxis.get_major_formatter().set_powerlimits((0, 0))
    # axHistSo.xaxis.get_major_formatter().set_powerlimits((0, 0))

    # Main plot
    if diff:
        im = axMain.imshow(h, interpolation='nearest', origin='low', extent=[x[1], x[-1], y[0], y[-1]], aspect='auto', cmap='RdBu', vmin=-0.4, vmax=0.4)
    else:
        im = axMain.imshow(h, interpolation='nearest', origin='low', extent=[x[1], x[-1], y[0], y[-1]], aspect='auto', cmap='inferno')

    axMain.set_xlabel('Energy [keV]')
    axMain.set_ylabel('Standoff Distance [mm]')
    cbar = f.colorbar(im, ax=axHistSo, location='right', format='%.0e')
    cbar.ax.set_ylabel('Normalized Counts')

    # HistE
    h = np.nan_to_num( h )
    print h
    print h.sum(axis=0)
    axHistE.bar(x[:-1]+getBinWidth(x)*.5, h.sum(axis=0), width=getBinWidth(x)) 
    axHistE.set_xlim(axMain.get_xlim())

    # HistSo
    axHistSo.barh(y[:-1]+getBinWidth(y)*.5, h.sum(axis=1), height=getBinWidth(y)) 
    axHistSo.set_ylim(axMain.get_ylim())

    if pp:
        f.savefig(pp, format='pdf')

    f.show()

def plotSSfrac(x, dataFrac, dataFracErr, MCFrac, MCFracErr, art='standoff', pp=None):
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib import rc
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)

    # Create figure
    f, (axData, axRes) = plt.subplots(2, sharex=True, sharey=False)
    f.subplots_adjust(wspace=0, hspace=0.1)

    # Axis 1 - Data and MC
    axData.set_yscale('log', nonposy='clip')
    axData.errorbar(x, dataFrac, yerr=dataFracErr, fmt='x-', label='Data')
    axData.errorbar(x, MCFrac, yerr=MCFracErr, fmt='x-', label='MC')

    # Axis 2 - residuals
    res = (dataFrac - MCFrac) / dataFrac
    resErr = np.sqrt((MCFrac/dataFrac**2 * dataFracErr)**2 + (1/dataFrac * MCFracErr)**2) # abs(MCFrac/dataFrac) * np.sqrt( (MCFracErr/MCFrac)**2 + (dataFracErr/dataFracErr)**2 )
    axRes.errorbar(x, res, yerr=resErr)
    axRes.set_ylim(-0.5, 0.5)

    # Labels
    axData.legend(loc='best')
    axData.set_ylabel(r'SS/(SS + MS)')
    if art == 'standoff':
        axRes.set_xlabel(r'Standoff Distance [mm]')
    elif art == 'energy':
        axRes.set_xlabel(r'Energy [keV]')
    axRes.set_ylabel(r'(Data - MC)/Data')

    if pp:
        f.savefig(pp, format='pdf')

    f.show()

# == FILL ENERGY ==
# Get histograms
def fillEnergy(t, res, bins, hRange, art='ss', MC=True):
    Elist = []
    SOlist = []
    N = t.GetEntries()
    for i in range( N ):
        t.GetEntry(i)
        es = t.EventSummary
        if not MC:
            if art == 'ss':
                # E = es.energy_ss
                E = np.sum( np.array(es.cluster_energy) )
            else:
                # E = es.energy_ms
                E = np.sum( np.array(es.cluster_energy) ) * 2614.5 / 2526.97 # 2614.5 / 2651.59 
        else:
            # E = getNewEnergy(res, es.energy_mc)[0]
            energyCorrMC = 2614.5 / 2601.87
            # E = getNewEnergy(res, np.sum(np.array(es.cluster_energy)))[0] * energyCorrMC
            energy = np.array( es.cluster_energy ) * energyCorrMC
            # energyFrac = energy/sum(energy)
            E = getNewEnergy(res, sum(energy))[0] # * energyFrac

            if not (E >= hRange[0][0] and E <= hRange[0][1]):
                continue

        # x, y, z = ps.getClusterPos(es, cut=True)
        so = es.standoff_distance

        Elist.append( E )
        SOlist.append( so )

    hE, edges = np.histogram(np.array(Elist), bins[0], hRange[0])
    hSO, edges = np.histogram(np.array(SOlist), bins[1], hRange[1])
    h, xedges, yedges = np.histogram2d(np.array(Elist), np.array(SOlist), bins=bins, range=hRange, normed=False)
    # Normalize by number of events
    # h /= N

    h = np.rot90(h)
    h = np.flipud(h)
    return h, hE, hSO, xedges, yedges

# === SMEAR MC RANDOM ===
def smearMCRandom(dataList, MCtree, bins, hRange, art='ss'):
    import re

    runNumList = [int(re.findall(r'\d+', data.split('/')[-1])[-1]) for data in dataList]
    print runNumList
    resList = [getResolutionDB(runNum, art) for runNum in runNumList]
    print resList

    # Store entries
    Elist = []
    SOlist = []

    N = MCtree.GetEntries()
    entryList, nList = getEntryList(dataList, N, MCtree)
    MCtree.SetEntryList( entryList )

    for k in range(entryList.GetN()):
        MCtree.GetEntry( k )
        es = MCtree.EventSummary

        res = resList[np.digitize(k, nList)]
        E = getNewEnergy(res, es.energy_mc)[0]

        if not (E >= hRange[0][0] and E <= hRange[0][1]):
            continue

        Elist.append( E )
        SOlist.append( es.standoff_distance )

    hE, edges = np.histogram(np.array(Elist), bins[0], hRange[0])
    hSO, edges = np.histogram(np.array(SOlist), bins[1], hRange[1])
    h, xedges, yedges = np.histogram2d(np.array(Elist), np.array(SOlist), bins=bins, range=hRange, normed=False)
    # Normalize by number of events
    # h /= N

    h = np.rot90(h)
    h = np.flipud(h)
    return h, hE, hSO, xedges, yedges

# === GET ENTRY LIST ===
# dataList: Contains the names of the data runs
#           to extract their run numbers
# N:        Number of events in the MC run
def getEntryList(dataList, N, t):
    from random import randint

    NList = getEntriesFromList(dataList)
    # As there might be more data than MC, 
    # scale MC events accordingly
    NList = [float(n)/sum(NList) * N for n in NList]
    M = range( N )

    # Get cumulative values of entries
    nList = [NList[0]]
    for i in range(1, len(NList)):
        nList.append( sum(NList[0:i]) + NList[i] )

    entryList = ROOT.TEntryList()
    while len( M ) > 0:
        k = M.pop(randint(0, len(M)-1))
        entryList.Enter(k, t)

    return entryList, nList

# === ENTRIES FROM LIST ===
def getEntriesFromList(dataList):
    NList = []
    for data in dataList:
        f = ROOT.TFile.Open(data, 'READ')
        t = f.Get('dataTree')

        NList.append( t.GetEntries() )
        f.Close()

    return NList

# === FILL SS FRAC ===
def fillSSFrac(t, bins, hRange, art='standoff', MC=False, resolution=None):
    N = t.GetEntries()
    print N
    typeList = []
    SOlist = []
    bins = np.linspace(hRange[0], hRange[1], bins + 1)[:-1]
    SSlist = [0]*bins
    MSlist = [0]*bins

    for i in range( N ):
        t.GetEntry(i)
        es = t.EventSummary
        if art == 'standoff':
            soIdx = np.digitize(es.standoff_distance, bins) - 1
        elif art == 'energy':
            if MC:
                if es.multiplicity > 1.1:
                    en = getNewEnergy(resolution['ms'], es.energy_mc)[0]
                else:
                    en = getNewEnergy(resolution['ss'], es.energy_mc)[0]

                soIdx = np.digitize(en, bins) - 1
            else:
                soIdx = np.digitize(es.energy, bins) - 1

        if es.multiplicity > 1.1:
            MSlist[soIdx] += 1
        else:
            typeList.append( 0 )
            SSlist[soIdx] += 1

    SSlist, MSlist = np.array( SSlist ), np.array( MSlist )
    SSlist = SSlist.astype(np.float)
    SSfrac = SSlist / (SSlist + MSlist)
    SSfracErr = np.sqrt( (MSlist/(SSlist + MSlist)**2 * np.sqrt(SSlist))**2 + (SSlist/(SSlist + MSlist)**2 * np.sqrt(MSlist))**2 )
    print SSlist
    print MSlist
    print
    print SSfrac
    print SSfracErr

    return bins, np.nan_to_num( SSfrac ), np.nan_to_num( SSfracErr )

# === FILL SS FRAC MC ===
def fillSSFracMC(dataList, t, bins, hRange, art='standoff'):
    import re
    typeList = []
    SOlist = []
    bins = np.linspace(hRange[0], hRange[1], bins + 1)[:-1]
    SSlist = [0]*bins
    MSlist = [0]*bins

    # Get resolutions
    runNumList = [int(re.findall(r'\d+', data.split('/')[-1])[-1]) for data in dataList]
    resListSS = [getResolutionDB(runNum, 'ss') for runNum in runNumList]
    resListMS = [getResolutionDB(runNum, 'ms') for runNum in runNumList]

    N = t.GetEntries()
    if art == 'energy':
        entryList, nList = getEntryList(dataList, N, t)
        t.SetEntryList( entryList )
        N = entryList.GetN()

    for i in range( N ):
        t.GetEntry(i)
        es = t.EventSummary

        if art == 'standoff':
            soIdx = np.digitize(es.standoff_distance, bins) - 1
        elif art == 'energy':
            if es.multiplicity > 1.1:
                res = resListMS[np.digitize(i, nList)]
            else:
                res = resListSS[np.digitize(i, nList)]
            en = getNewEnergy(res, es.energy_mc)[0]
            soIdx = np.digitize(en, bins) - 1

        if es.multiplicity > 1.1:
            MSlist[soIdx] += 1
        else:
            typeList.append( 0 )
            SSlist[soIdx] += 1

    SSlist, MSlist = np.array( SSlist ), np.array( MSlist )
    SSlist = SSlist.astype(np.float)
    SSfrac = SSlist / (SSlist + MSlist)
    SSfracErr = np.sqrt( (MSlist/(SSlist + MSlist)**2 * np.sqrt(SSlist))**2 + (SSlist/(SSlist + MSlist)**2 * np.sqrt(MSlist))**2 )
    print SSlist
    print MSlist
    print
    print SSfrac
    print SSfracErr

    return bins, np.nan_to_num( SSfrac ), np.nan_to_num( SSfracErr )

# === PLOT ENERGY STANDOFF ===
def plotEnergyStandoff(hData, hMC, x, bins, hRange, pp=None):
    hData, hMC = np.nan_to_num( hData ), np.nan_to_num( hMC )
    binWidth = (hRange[1][1] - hRange[1][0]) / bins[1]

    # Loop over standoff bins
    for i in range(bins[1]):
        binCenter = (i + .5) * binWidth
        plotEnergyRes(hData[i], hMC[i], x, '@ Standoff distance %.2f' % binCenter, pp)

# == SUPPORT ==
def getNorm(h, hE, normRange):
    N = float( np.sum(hE[normRange[0]:normRange[1]]) )
    return np.array( h ) / N

# For 1d histograms
def getErrorNorm(h):
    hSum = float( sum( h ) )
    hErr = np.sqrt(h) / hSum
    h = h / hSum
    return h, hErr

def getBinWidth(x):
    return abs(x[1] - x[0])

def getRes(res, E):
    return np.sqrt( res[0]**2*E + res[1]**2 + (res[2]*E)**2 )

def getNewEnergy(res, E):
    sigma = getRes(res, E)
    mu = E

    return np.random.normal(mu, sigma, 1)

# === GETRESOLUTION ===
# Determine resolution of data by fitting gauss distributions
def getResolution(energyBins, energyValues, src='Th', art='ss'):
    import matplotlib.pyplot as plt
    from matplotlib import rc
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)

    import scipy.special
    import scipy.optimize
    print 'Determing resolution from data...'

    if src == 'Th':
        # Peak positions
        # energies = [1078, 1170, 1620, 2103, 2614]
        energies = [2614, 2103, 1620, 1078, 725]
        # Fit ranges
        # fitRange = [[2450, 2800], [1990, 2200], [1510, 1700], [1015, 1200], [800, 950]]
        if art == 'ss':
            fitRange = [[2465, 2800], [1990, 2200], [1480, 1760], [810, 920], [650, 800]]
            initPars = [[2614, 40, 1.e3, 5.e2, 1.e2], [2103, 40, 5.e2, 2.5e2, 1.e2], [1620, 40, 5.e2, 2.5e2, 1.e2], [855, 30, 10.e3, 5.e3, 20.e3], [725, 40, 100, -10, 120]]
        elif art == 'ms':
            fitRange = [[2400, 2800], [1990, 2200], [1505, 1760], [810, 915], [655, 800]]
            initPars = [[2610, 40, 1.e3, 5.e2, 1.e2], [2103, 40, 5.e2, 2.5e2, 1.e2], [1620, 40, 5.e2, 2.5e2, 1.e2], [855, 30, 10.e3, 5.e3, 20.e3], [721, 30, 10.e3, 0.1, 20.e3]]
        else:
            return False

        # Initial parameters
        # initPars = [[2614, 40, 1.e3, 5.e2, 1.e2], [2103, 40, 5.e2, 2.5e2, 1.e2], [1620, 40, 5.e2, 2.5e2, 1.e2], [1170, 30, 4000, 500., 8000], [850, 40, 100, -10, 120]]
    elif src =='charge':
        energies = [1510, 2014, 2532]
        # Fit ranges
        if art == 'ss':
            fitRange = [[1250, 1750], [1750, 2180], [2180, 2930]]
            initPars = [[1510, 80, 1.e3, 5.e2, 1.e2], [2014, 80, 5.e2, 2.5e2, 1.e2], [2532, 40, 5.e2, 2.5e2, 1.e2]]
        elif art == 'ms':
            if PHASE2:
                fitRange = [[1250, 1800], [1800, 2300], [2300, 3000]]
            else:
                fitRange = [[1250, 1800], [1800, 2250], [2260, 2990]]
            initPars = [[1560, 80, 1.e3, 5.e2, 1.e2], [2100, 80, 5.e2, 2.5e2, 1.e2], [2615, 60, 175.e2, -1.e3, 20.e3]]
        else:
            return False

    else:
        return False

    # energyBins contains the edges of the histogram.
    # Therefore, shift it by the half bin width to center the energy values
    binWidth = energyBins[1] - energyBins[0]
    energyBins = np.array( energyBins ) + 0.5 * binWidth

    # Fit functions
    def gauss(x, mu, sigma, A, off):
        return A * np.exp( -(x-mu)**2 / (2.*sigma**2) ) + off

    def shift(a, b, mu, sigma):
        return np.sqrt(2./np.pi)*float(b)/a*sigma

    def gaussErf(x, mu, sigma, a, b, c):
        return gauss(x+shift(a, b, mu, sigma), mu, sigma, a, c) + b * scipy.special.erf((x+shift(a, b, mu, sigma) - mu) / (np.sqrt(2) * sigma)) + abs(b)

    # To present results
    symbolList = ['mu', 'sigma', 'a', 'b', 'c']
    symbolListSigma = ['sigma', 'a', 'c']
    poptList = []
    perrList = []
    poptSigmaList = []

    # Perform fits
    for i in range( len(energies) ):
        # Get data in range
        lim = fitRange[i]
        energy = [en for en in zip(energyBins, energyValues) if (en[0] >= lim[0] and en[0] <= lim[1])]
        en, val = np.array(energy)[:,0], np.array(energy)[:,1]

        # Get mean values of peaks
        p0 = initPars[i]
        try:
            popt, pcov = scipy.optimize.curve_fit(gaussErf, en, val, p0)
            perr = np.sqrt(np.diag(pcov))
        except:
            popt = p0
            perr = np.zeros(len(popt))

        poptList.append( popt )
        perrList.append( perr )

        # Get sigma of peaks
        mean, sigma = np.abs( popt[0:2] )              # Fix the mean in fit
        print mean, sigma
        # Select data range
        energy = [en for en in zip(energyBins, energyValues) if (en[0] >= (mean - 2* sigma) and en[0] <= (mean + 2 * sigma))]
        print energy
        en, val = np.array(energy)[:,0], np.array(energy)[:,1]

        p0Sigma = [popt[1], 200, 200]        # sigma, amplitude, offset
        poptSigma = p0Sigma
        try:
            poptSigma, pcovSigma = scipy.optimize.curve_fit(lambda x, sigma, A, off: gauss(x, mean, sigma, A, off), en, val, p0Sigma)
        except:
            poptSigma = p0Sigma
            pcovSigma = np.zeros(len(poptSigma))
        poptSigmaList.append(poptSigma)
        # poptSigma = poptSigma
        perrSigma = np.sqrt(np.diag(pcovSigma))

        print 'Results for peak @%d keV' % p0[0]
        print '---------'
        for j, p in enumerate(popt):
            print '%s = %.2f +/- %.2f' % (symbolList[j], p, perr[j])
        print 'Shift = %.2f' % shift(popt[2], popt[3], popt[0], popt[1])
        print '---------'
        for j, p in enumerate(poptSigma):
            print '%s = %.2f +/- %.2f' % (symbolListSigma[j], p, perr[j])
        print

    # Plot fit results
    plt.plot(energyBins[:-1], energyValues)
    for j, popt in enumerate(poptList):
        mu, sig, a, b, c = popt

        lim = fitRange[j]
        x = np.linspace(lim[0], lim[1], 1000)
        y = gaussErf(x, *popt)
        yErf = b * (scipy.special.erf((x+shift(a, b, mu, sig) - mu) / (np.sqrt(2) * sig))) + abs(b)
        yGauss = gauss(x+shift(a, b, mu, sigma), mu, sig, a, c)

        # ySigma = gauss(x, popt[0], *poptSigmaList[j])
        # plt.plot(en, val, label='Data')
        plt.plot(x, y, label='%s keV Peak' % initPars[j][0])
        plt.plot(x, yErf, ls='--')
        plt.plot(x, yGauss, ls='--')
        # plt.plot(x, ySigma)
        plt.legend()
    plt.xlabel('Energy [keV]')
    plt.ylabel('Counts')
    plt.show()

    # Get resolution by fitting function to sigma
    mean = np.array( poptList )[:,0]
    meanErr = np.array( perrList )[:,0]
    # sigma = poptSigma[0] 
    sigma = np.abs(np.array( poptList ))[:,1]
    # sigmaErr = perrSigma[0]
    sigmaErr = np.abs(np.array( perrList ))[:,1]

    # Resolution function
    def resFunc(x, p0, p1, p2):
        return np.sqrt(p0**2*x + p1**2 + p2**2*x**2)

    p0 = [0.7, 28, 0.001]
    popt, pcov = scipy.optimize.curve_fit(resFunc, mean, sigma, p0)
    perr = np.sqrt(np.diag( pcov ))
    print 'Results for resolution:'
    print '---------'
    print popt, perr

    # Plot result
    plt.clf()

    # Data
    plt.errorbar(mean, sigma, xerr=meanErr, yerr=sigmaErr, fmt='x', label='Data')

    # Fit
    x = np.linspace(mean[-1] - 100, mean[0] + 100, 1000)
    y = resFunc(x, *popt)
    plt.plot(x, y, label='Fit')
    plt.legend()

    plt.xlabel(r'Energy')
    plt.ylabel(r'$\sigma$ [keV]')
    plt.show()

    # Return result
    return popt

# === GET RESOLUTION DB ===
# Get resolution from database
def getResolutionDB(runNum, art='ss'):
    import time
    from datetime import datetime
    import re

    resol = ROOT.EXOEnergyResol.GetInstanceForFlavor('2016-v1-weekly', 'energy-mcbased-fit')

    # Get datetime of run
    try:
        t = getXML(runNum, 'startTime')
    except:
        return stdResolution

    t = int( datetime.strptime(t.split('.')[0], '%Y-%m-%dT%H:%M:%S').strftime('%s') )

    if art == 'ss':
        resString = resol.ResolutionString('Rotated', 'x', 1, t, 0)
    else:
        resString = resol.ResolutionString('Rotated', 'x', 2, t, 0)

    res = [float(x) for x in re.findall("\d+\.\d+", resString)]
    return [np.sqrt(r) for r in res]

def getXML(runNum, art='startTime'):
    import gzip, urllib2, xmltodict, ast
    from StringIO import StringIO

    url = 'http://exo-data.slac.stanford.edu/ExoDatacat/rest/runs/%d' % runNum
    request = urllib2.Request( url )
    request.add_header('Accept', 'application/xml')
    request.add_header('Accept-encoding', 'gzip')

    response = urllib2.urlopen( request )
    if response.info().get('Content-Encoding') == 'gzip':
        uzdata = response.read()
        buf = StringIO( uzdata )
        f = gzip.GzipFile(fileobj=buf)
        data = f.read()
    else:
        return

    data = xmltodict.parse( data )
    return data['run'][art]

