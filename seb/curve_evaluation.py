#!/usr/bin/env python

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import signal
from scipy.optimize import curve_fit
from scipy.stats import ks_2samp
from scipy.spatial import distance
from uncertainties import ufloat
from uncertainties import umath
import logging
import numpy as np
import numpy.fft as fft
import cPickle
import os

dataDir = 'result/data/'
plotDir = 'result/plot/'
RING_DIST = 16.8656

def main():
	fnReal = 'real/realSide_60.p'
        evalPara(fnReal, fnReal)
        return

	evalResDict = {}
        for fn in os.listdir( dataDir ):
		if os.path.isfile( dataDir + fn ):
			evalResDict[fn] = evalPara(fn, fnReal)

        evalResDictSort = sorted(evalResDict.items(), key=lambda e:e[1][2])
	for key, val in evalResDictSort:
		print key, val

def evalPara(fn, fnReal, cut=True):
	print 'Processing %s' % fn
	print '============='
	print
	# fnReal = 'real/real_60.p'
	pdfOut = plotDir + fn.split('.')[0] + '.pdf'
	pp = PdfPages( pdfOut )

	# === DIFF ===
	# Bin 0
    	dataRange = (-170, -20)
	plotRange = (0.1, -0.6)

	# resultRealBin0 = evaluateFile(fnReal, 0, dataRange, plotRange, None, 'diff')
        try:
            resultBin0 = evaluateFile(fn, 0, dataRange, plotRange, pp, 'diff')
            resultBin0.update( compareFile(fnReal, fn, 0, dataRange, plotRange, pp, 'diff') )
            evalResBin0 = eval(resultBin0, cut, True)
        except:
            evalResBin0 = (0., 0., 0.)
	
	print
	print '===================='
	print 'Result of evaluation (Bin0):', evalResBin0
        print '\n'

	# Other bins
	plotRange = (0.1, -0.1)

	# Most of the time garbage
	# for nBin in range(1, 2):
	# evaluateFile(fnReal, nBin, left, right)

        nBin = 1
        try:
            resultBin1 = evaluateFile(fn, nBin, dataRange, plotRange, pp, 'diff')
            resultBin1.update( compareFile(fnReal, fn, nBin, dataRange, plotRange, pp, 'diff') )
            evalResBin1 = eval(resultBin1, cut, True)
        except:
            evalResBin1 = (0., 0., 0.)

	print
	print '===================='
	print 'Result of evaluation (Bin1):', evalResBin1
        print '\n'

	# === DATA ===
	# Bin 0
	plotRange = (0.12, 0.035)

        try:
            resultBin0Data = evaluateFile(fn, 0, dataRange, plotRange, pp, 'data')
            resultBin0Data.update( compareFile(fnReal, fn, 0, dataRange, plotRange, pp, 'data') )
            evalResData = eval(resultBin0Data, cut, False)
        except:
            evalResData = (0., 0., 0.)

	# Other bins
	plotRange = (0.1, -0.1)
	for nBin in range(1, 1):
            try:
                evaluateFile(fn, nBin, dataRange, plotRange, pp, 'data')
                compareFile(fnReal, fn, nBin, dataRange, plotRange, pp, 'data')
            except:
                pass

	pp.close()
	return 1./7 * (2*np.array(evalResBin0) + np.array(evalResBin1) + 4*np.array(evalResData))

def evaluateFile(fn, nBin, dataRange=(-170, -20), plotRange=(0.1, -0.6), pp=None, key='diff', show=False):
	d = cPickle.load( open(dataDir + fn, 'rb') )

	x = np.array( d['x'] )
	data = np.array( d[key] )[:,nBin]
	if key == 'diff':
	    data_err = np.array( d['diff_err'] )[:, nBin]
	else:
	    data_err = np.array([])

	return evaluateData(x, data, data_err, key,dataRange, plotRange, pp, show)

def compareFile(fn1, fn2, nBin, dataRange=(-170, -20), plotRange=(0.1, -0.6), pp=None, key='diff', show=False):
	d1 = cPickle.load( open(dataDir + fn1, 'rb') )
	d2 = cPickle.load( open(dataDir + fn2, 'rb') )

	x = np.array( d1['x'] )
	data1 = np.array( d1[key] )[:,nBin]
	data2 = np.array( d2[key] )[:,nBin]

	return compareData(x, data1, data2, key, dataRange, plotRange, pp, show)

def evaluateData(x, data, data_err=np.array([]), key='diff', dataRange=(-170, -20), plotRange=(0.1, -0.6), pp=None, show=False):
	resultReturnDict = {}

	# = Data within range =
	rangeMin, rangeMax = dataRange

	'''
	# Test for FFT
	x = np.linspace(-1000, 1000, 10000)
	data = np.sin(2*np.pi * x/16.8)
	data_err = np.array([])
	'''

	xFit, dataFit = getDataRange(x, data, rangeMin, rangeMax)

	# = Fit =
        if key == 'diff':
            popt, perr, dataLin = fitEval(xFit, dataFit, linear)
            dataLin = linear(x, *popt)
            print 'Fit result:'
            print 'f(z) = (%f +/- %f) * z + (%f +/- %f)\n' % (popt[0], perr[0], popt[1], perr[1])
        else:
            pParabola = [-100, -0.02, 0.08]
            popt, perr, dataLin = fitEval(xFit, dataFit, parabola, pParabola)
            dataLin = parabola(x, *popt)

	# = Deviation =
        if key == 'diff':
            yFit = linear(xFit, *popt)
        else:
            yFit = parabola(xFit, *popt)

        # resultReturnDict['meanAbsDev'], resultReturnDict['stdDev'] = getDeviations(dataFit, yFit)

	# = FFT =
	dataFlat = dataFit - yFit
	xFFT, dataFFT = fftEval(xFit, dataFlat)

	wavelength = xFFT[list(dataFFT).index( max( dataFFT ) )]
	resultReturnDict['waveFFT'] = wavelength
	print 'Wavelength of FFT: %f' % wavelength

	# = Auto correlation =
	dataAC = autoCorr( dataFlat )
	t = xFit[1] - xFit[0]
	xAC = np.linspace(0, xFit[-1] - xFit[0], len(dataAC))

	# AC Fit
	pAC = [0., 0.5, 16.8]
	poptAC, perrAC, dataACFit = fitEval(xAC, dataAC, sine, pAC)
	xACFit = np.linspace(xAC[0], xAC[-1], 10000)
	dataACFit = sine(xACFit, *poptAC)

	# AC Fit result
	aAC = ufloat(poptAC[1], perrAC[1])
	bAC = ufloat(poptAC[2], perrAC[2])
	z0AC = ufloat(poptAC[0] + bAC.n/4, perrAC[0])
	resultReturnDict['ACb'] = bAC
	resultReturnDict['ACz0'] = z0AC
	print 'f(z) = {:2.3f} * sin(2*pi*(z - {:2.2f}) / {:2.2f})'.format(aAC, z0AC, bAC)

	# AC FFT
	xAC_FFT, dataAC_FFT = fftEval(xFit, dataAC)

	autoWavelength = xFFT[list(dataAC_FFT).index( max( dataAC_FFT ) )]
	resultReturnDict['waveAC'] = autoWavelength
	print 'Wavelength of auto correlation: %f\n' % autoWavelength

	# = Plots =
	fig1 = plt.figure()
	ax1 = fig1.add_subplot(111)

	ax1.set_xlabel('z [mm]')
	ax1.set_ylabel('data')

	plotLeft, plotRight = plotRange
	ax1.set_xlim(left=-190, right=-10)
        ax1.set_ylim(top=plotLeft, bottom=plotRight)

	# diff
	if data_err.any():
		ax1.errorbar(x, data, yerr=data_err, fmt='-x', ecolor='g')
	else:
		ax1.plot(x, data)

	# linear fit
	ax1.plot(x, dataLin)

	# FFT
	# Limits and labels
	fig2, axFFT = plt.subplots(2, sharex=False)
	axFFT[0].set_ylabel('FFT')
	axFFT_ax1 = axFFT[0].twinx() 	# Second y-axis
        axFFT[0].set_xlabel('Wavelength [mm]')
	axFFT_ax1.set_ylabel('FFT of autocorrelation')

	axFFT[1].set_xlabel('Shift [mm]')
	axFFT[1].set_ylabel('Flat')
	axFFT_ax2 = axFFT[1].twinx()
	axFFT_ax2.set_ylabel('AC')

	# Plot
	# Auto correlation
	ln1_1 = axFFT[0].plot(xFFT, dataFFT, '-x', label='FFT', color='b')
	ln1_2 = axFFT_ax1.plot(xAC_FFT, dataAC_FFT, '.-', label='FFT of auto correlation', color='g')

	axFFT[1].plot(xAC, dataAC, label='AC', color='b')
	axFFT[1].plot(xACFit, dataACFit, label='Fit', color='g')

	axFFT[0].set_xlim(left=8, right=64)
        axFFT[1].set_xlim(left=0, right=140)

	# Combine labels for legend
	ln1 = ln1_1 + ln1_2
	labs1 = [l.get_label() for l in ln1]

	axFFT[0].legend(ln1, labs1)
	axFFT[1].legend(loc='lower center')

	if show:
		fig1.show()
		fig2.show()
		raw_input('')

	if pp:
		fig1.savefig(pp, format='pdf')
		fig2.savefig(pp, format='pdf')

	plt.cla()

	return resultReturnDict

def compareData(x, data1, data2, key, dataRange, plotRange, pp=None, show=False):
	resultReturnDict = {}

	# = Data within range =
	rangeMin, rangeMax = dataRange

	xRange, data1Range = getDataRange(x, data1, rangeMin, rangeMax)
	xRange, data2Range = getDataRange(x, data2, rangeMin, rangeMax)

	# = Fit =
        if key == 'diff':
            fitFunc = linear
        else:
            fitFunc = parabola

        popt1, perr1, data1Fit = fitEval(xRange, data1Range, fitFunc)
        popt2, perr2, data2Fit = fitEval(xRange, data2Range, fitFunc)

        if key == 'diff':
            m1 = ufloat(popt1[0], perr1[0])
            m2 = ufloat(popt2[0], perr2[0])
            t1 = ufloat(popt1[1], perr1[1])
            t2 = ufloat(popt2[1], perr2[1])

            # The closer to 1 the better
            mFrac = (m1 / m2).n
            tFrac = (t1 / t2).n
            
            sumFrac = 0
            if mFrac <= 1:
                sumFrac += np.exp(mFrac - 1)
            else:
                sumFrac += np.exp(-(mFrac - 1))

            if tFrac <= 1:
                sumFrac += np.exp(tFrac - 1)
            else:
                sumFrac += np.exp(-(tFrac - 1))

            resultReturnDict['m'] = 0.5 * sumFrac

        else:
            a1 = ufloat(popt1[1], perr1[1])
            a2 = ufloat(popt2[1], perr2[1])
            b1 = ufloat(popt1[2], perr1[2])
            b2 = ufloat(popt2[2], perr2[2])

            aFrac = (a1 / a2).n
            bFrac = (b1 / b2).n

            sumFrac = 0
            if aFrac <= 1:
                sumFrac += np.exp(aFrac - 1)
            else:
                sumFrac += np.exp(-(aFrac - 1))

            if bFrac <= 1:
                sumFrac += np.exp(bFrac - 1)
            else:
                sumFrac += np.exp(-(bFrac - 1))

            resultReturnDict['m'] = 0.5 * sumFrac

        # = Deviation =
        if key == 'diff':
            yFit1 = linear(xRange, *popt1)
            yFit2 = linear(xRange, *popt2)
        else:
            yFit1 = parabola(xRange, *popt1)
            yFit2 = parabola(xRange, *popt2)

        print 'Deviation of data1: '
        meanAbsDev1, stdDev1 = getDeviations(data1Range, yFit1)
        print 'Deviation of data2: '
        meanAbsDev2, stdDev2 = getDeviations(data2Range, yFit2)

        # Store fractions of deviations
        resultReturnDict['meanAbsDev'] = meanAbsDev1 / meanAbsDev2
        resultReturnDict['stdDev'] = stdDev1 / stdDev2

	# = Cumulative sum =
	data1Cum = np.cumsum( data1 )
	data2Cum = np.cumsum( data2 )

	# = Kolmogorov Smirnov =
	D_ks_min, p_ks = ks_2samp(data1Range, data2Range)
	D_ks_max = abs(np.cumsum(data1Range) - np.cumsum(data2Range))[-1]
	resultReturnDict['D_min'] = D_ks_min
	resultReturnDict['D_max'] = D_ks_max
	print 'Kolmogorov-Smirnov (min, max) = (%.2f, %.2f)\n' % (D_ks_min, D_ks_max)

	# = Cross correlation =
	xCC, dataCC = crossCorr(xRange, data1Range - data1Fit, data2Range - data2Fit)

	pCC = [0., 0.01, 16.8]
	poptCC, perrCC, dataCCFit = fitEval(xCC, dataCC, sine, pCC)
	xCCFit = np.linspace(xCC[0], xCC[-1], 10000)
	dataCCFit = sine(xCCFit, *poptCC)

	aCC = ufloat(poptCC[1], perrCC[1])
	bCC = ufloat(poptCC[2], perrCC[2])
	z0CC = ufloat(poptCC[0] + bCC.n/4, perrCC[0])
	resultReturnDict['CCb'] = bCC
	resultReturnDict['CCz0'] = z0CC

        print 'Cross Correlation:'
	print 'f(z) = {:2.3f} * sin(2*pi*(z - {:2.2f}) / {:2.2f})\n'.format(aCC, z0CC, bCC)

	# = Distance =
	# distDList = [euclideanDist, jaccardDist, pearsonChiSqDist]
	distDList = [euclideanDist, pearsonChiSqDist]
	distSList = [innerProductDist, harmonicMeanDist, correlationDist]

	distDRes = []
	distSRes = []
	for func in distDList:
		distDRes.append( func(data1Range, data2Range) )
	for func in distSList:
		distSRes.append( func(data1Range, data2Range) )
	resultReturnDict['distS'] = distSRes
	resultReturnDict['distD'] = distDRes

        print 'Distances (dissimilarity):', distDRes
        print 'Distances (similarity):', distSRes

	# = Plots =
	# Cross correlation
	# Limits and labels
	fig1, axarr1 = plt.subplots(2, sharex=False)
	axarr1[0].set_xlim(left=-190, right=0)

	plotLeft, plotRight = plotRange
	axarr1[0].set_ylim(top=plotLeft, bottom=plotRight)
	axarr1[0].set_ylabel('Data')

	axarr1[1].set_xlim(left=-100, right=100)
	axarr1[1].set_xlabel('z [mm]')
	axarr1[1].set_ylabel('Cross correlation')

	# Plot
	axarr1[0].plot(x, data1, label='Data1')
	axarr1[0].plot(x, data2, label='Data2')
	axarr1[0].legend()
	axarr1[1].plot(xCC, dataCC)
	axarr1[1].plot(xCCFit, dataCCFit)

	# Cumulative sum
	# Limits and labels
	fig2, ax2_1 = plt.subplots()
	ax2_1.set_xlabel('z [mm]')
	ax2_1.set_ylabel('Cumulative sum')

	ax2_1.set_xlim(left=-190, right=0)
	ax2_1.set_ylim(top=max(max(data1Cum), max(data2Cum)) + 1, bottom=-20)

	# Plot
	ax2_1.plot(x, data1Cum, label='Data1')
	ax2_1.plot(x, data2Cum, label='Data2')

	ax2_1.legend()

	# = Show Plots =
	if show:
		fig1.show()
		fig2.show()
		raw_input('')

	if pp:
		fig1.savefig(pp, format='pdf')
		fig2.savefig(pp, format='pdf')

	plt.cla()

	return resultReturnDict

# === DATA RANGE ===
def getDataRange(x, data, left=-170, right=-20):
	xRange = []
	dataRange = []
	for i in range( len(x) ):
		if x[i] >= left and x[i] <= right:
			xRange.append( x[i] )
			dataRange.append( data[i] )

	return np.array( xRange ), np.array( dataRange )

# === FFT ===
def fftEval(x, data):
	n = len( data )
	k = np.arange( n )
	T = x[-1] - x[0]
	frq = k / T
	frq = frq[range(n/2)]

	dataFFT = np.abs( fft.rfft( data ) ) / n
	dataFFT = dataFFT[range(n/2)]

	xFFT = 1 / frq
	return xFFT, dataFFT

# === FIT EVAL ===
def fitEval(x, data, func, p=None):
        data = np.nan_to_num(data)
	popt, pcov = curve_fit(func, x, data, p)
	perr = np.sqrt(np.diag(pcov))
	dataFit = func(x, *popt)

	return popt, perr, dataFit

# === FIT FUNCTIONS ===
def linear(x, m, t):
	return m * x + t

def parabola(x, x0, a, b):
        return a * (x - x0)**2 + b

def sine(x, x0, a, b):
	return a * np.sin(2*np.pi*(x-x0) / b)

# === DEVIATIONS ===
def meanAbsoluteDev(x, x_):
	n = len( x )
	return 1./n * sum( abs(x-x_) )

def standardDev(x, x_):
	n = len( x )
	return np.sqrt( 1./n * sum( (x-x_)**2 ) )

def getDeviations(y1, y2):
	meanAbsDev = meanAbsoluteDev(y1, y2)
	stdDev = standardDev(y1, y2)
	print 'Absolute deviation: %f' % meanAbsDev
	print 'Standard deviation: %f\n' % stdDev

        return meanAbsDev, stdDev

# === CORRELATION ===
# Estimated of auto correlation, see:
# http://en.wikipedia.org/wiki/Autocorrelation
def autoCorr(x):
	n = len(x)
	variance = x.var()
	x = x-x.mean()
	r = np.correlate(x, x, mode='full')[-n:]
	result = r/(variance*(np.arange(n, 0, -1)))
	# result = np.correlate(x, x, mode='full')
	return result

def crossCorr(x, data1, data2):
        n = len(x)
        t = x[-1] - x[0]
        k = np.arange(-n+1, n)
        frq = k / (n / t)

        d = np.correlate(data1, data2, mode='full')
	return frq, d

# === WAVELET ===
def wavelet(x, data, show=False):
    left, right = x[0], x[-1]

    widths = np.arange(1, 2)
    cwtmatr = signal.cwt(data, signal.ricker, widths)
    cwtproj = cwtmatr.sum(axis=0)
    cwtproj /= sum( cwtproj )

    popt, perr, cwtprojFit = fitEval(x, data, sine)
    xFit = np.linspace(left, right, 10000)
    cwtprojFit = sine(xFit, *popt)

    if show:
        fig, axes = plt.subplots(2, sharex=True)

        axes[0].plot(x, cwtproj, '-x')
        axes[0].plot(xFit, cwtprojFit)
        axes[1].imshow(cwtmatr, extent=[left, right, 1, 10], cmap='jet', aspect='auto', vmax=abs(cwtmatr).max(), vmin=-abs(cwtmatr).max())

        fig.show()

# === DISTANCE ===
# dissimilarity
def euclideanDist(data1, data2):
	return np.sqrt( sum( np.abs(data1 - data2)**2 ) )

# dissimilarity
def jaccardDist(data1, data2):
        return distance.jaccard(data1, data2)

# similarity
def innerProductDist(data1, data2):
	return sum( data1 * data2 )

# similarity
def harmonicMeanDist(data1, data2):
	return 2 * sum( (data1 * data2)/(data1 + data2) )

# dissimilarity
def diceDist(data1, data2):
	return distance.dice(data1, data2)

# similarity
def correlationDist(data1, data2):
	return distance.correlation(data1, data2)

# dissimilarity
def pearsonChiSqDist(data1, data2):
	return sum( (data1 - data2)**2 / data2 )

# === EVALUATION ===
def eval(resDict, cut=False, diff=True):
        # if not diff => data
        # NOTE: Currently overridden!
        diff = True

	# Divide in five components:
	# Kolmogorov-Smirnov, distance, slope, wavelength, deviation

	# = KS =
        if diff:
            D_min, D_max = resDict['D_min'], resDict['D_max']
            d_KS = 1./2 * ( np.exp(-D_min) + np.exp(-D_max) )

	# = Distance =
            distS, distD = resDict['distS'], resDict['distD']
            distS_sum = sum( [ np.exp( -1./abs(dS) ) for dS in distS ] ) / len( distS )
            distD_sum = sum( [ np.exp( -abs(dD) ) for dD in distD ] ) / len( distD )
            # d_dist = 1./2 * ( distS_sum + distD_sum )
            d_dist = distD_sum

	# = Slope =
        # Normalisation to 1 happens in compareData
	d_m = resDict['m']

	# = Wavelength =
        if not cut:
            bCC = resDict['CCb'].n
            if bCC <= RING_DIST:
                    d_w = np.exp(bCC - RING_DIST)
            else:
                    d_w = np.exp(-(bCC - RING_DIST))

        # = Deviation =
        meanAbsDev = resDict['meanAbsDev']
        stdDev = resDict['stdDev']

        d_dev = 0
        if meanAbsDev <= 1:
            d_dev += np.exp(meanAbsDev - 1)
        else:
            d_dev += np.exp(-(meanAbsDev - 1))

        if stdDev <= 1:
            d_dev += np.exp(stdDev - 1)
        else:
            d_dev += np.exp(-(stdDev - 1))

        d_dev *= 0.5

	# = Total =
        if not cut:
            dShape = 1./6 * (3*d_m + d_w + 2*d_dev)
        else:
            dShape = 1./5 * (3*d_m + 2*d_dev)

        if diff:
            dDist = 1./3 * (d_KS + 2*d_dist)
            d_g = 1./3 * (dDist + 2*dShape)
        else:
            # Set dDist = 1 for correct end result display
            # Has no effect on d_g
            dDist = 1.
            d_g = dShape

	return dShape, dDist, d_g

if __name__ == '__main__':
	main()

