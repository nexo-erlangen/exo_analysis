import sys
import numpy as np
from scipy import interpolate
from scipy import signal

import plot_support as ps
import plot_functions as pf

import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')

def findPeaks(xLs, dataLs, dataErrLs, doBox=False, title='default', outFn='standoff_field_shape.pdf', show=True):
        if not show:
            ps.setShow(False)

	ROOT.gStyle.SetOptStat(0)

	print 'Finding peaks...'
	for i in range(len(xLs)):
		print xLs[i], dataLs[i], dataErrLs[i]

	lHist = ROOT.TH1F('lHist', title, len(xLs), min(xLs), max(xLs))
	for i, item in enumerate( dataLs ):
		lHist.SetBinContent(i+1, item)
		lHist.SetBinError(i+1, dataErrLs[i])
	
	ROOT.TGraph(len(xLs), np.array(xLs), np.array(dataLs))

	# = SPLINE INTERPOLATION =
	ipl = interpolate.CubicSpline(np.array(xLs), np.array(dataLs))

	xNew = []
	yNew = []
	for x in np.linspace(min(xLs), max(xLs), 10000):
		xNew.append( x ) 
		yNew.append( ipl(x) )

	l = ROOT.TGraph(len(xNew), np.array(xNew), np.array(yNew))

	# = PEAK FINDING =
	'''
	# Bad results
	peaks = signal.find_peaks_cwt(-np.array(yNew), np.arange(2, 10))
	lPeak = ROOT.TGraph()
		for i in range(len(peaks)):
		xP, yP = xNew[peaks[i]], yNew[peaks[i]]
		lPeak.SetPoint(i, xP, yP)
	'''

	xRange = (-185, -5)
	yRange = (-0.4, 0.4)
	
	# Good results
	lPeak = ROOT.TGraph()
	peakMin = []
	for item in peakdet(yNew, 0.000001, xNew):
		peakMin += list(item)
	peakMin = np.array(peakMin)

	for i in range(len(peakMin)):
		xP, yP = peakMin[i]

		# Quick and dirty fix
		if xP > xRange[0] and xP < xRange[1] and yP > yRange[0] and yP < yRange[1]:
			lPeak.SetPoint(i, xP, yP)
		else:
			lPeak.SetPoint(i, 999, 999)

	latex = [ROOT.TLatex] * len(peakMin)
	for i, item in enumerate(latex):
		# item = ROOT.TLatex(lPeak.GetX()[i] + 2.5, lPeak.GetY()[i] - 0.002, str(round(peakMin[i][0], 2)))
		item = ROOT.TLatex(lPeak.GetX()[i], lPeak.GetY()[i], str(round(peakMin[i][0], 2)))
		item.SetTextSize(0.03)
		lPeak.GetListOfFunctions().Add(item)

	# = PLOT =
	c = ROOT.TCanvas()

	l.SetTitle(title)
	l.GetXaxis().SetTitle('z [mm]')
	l.Draw('AC')
	lHist.Draw('HIST E SAME')

	if doBox:
		l.GetXaxis().SetRangeUser(xRange[0], xRange[1])
		l.GetYaxis().SetRangeUser(yRange[0], yRange[1])
		lPeak.GetXaxis().SetRangeUser(xRange[0], xRange[1])
		lPeak.GetYaxis().SetRangeUser(yRange[0], yRange[1])
		lHist.GetXaxis().SetRangeUser(xRange[0], xRange[1])
		lHist.GetYaxis().SetRangeUser(yRange[0], yRange[1])
		boxes = [ROOT.TBox()] * 10
		for n in range(10):
			top, bottom = fieldShape(n)
			boxes[n].SetFillColorAlpha(15, 0.5)
			boxes[n].DrawBox(bottom, yRange[0], top, yRange[1]) 

	lPeak.SetMarkerStyle(3)
	lPeak.Draw('same P')

	c.SaveAs(outFn)
        if show:
            raw_input('')
	del c

        ps.setShow(True)

def fieldShape(n):
	FIELDRING_WIDTH = 9.652
	top = -18.2753 - n*16.8656
	bottom = top - FIELDRING_WIDTH
	return (top, bottom)

def peakdet(v, delta, x = None):
	"""
	Converted from MATLAB script at http://billauer.co.il/peakdet.html
	
	Returns two arrays
	
	function [maxtab, mintab]=peakdet(v, delta, x)
	%PEAKDET Detect peaks in a vector
	%        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
	%        maxima and minima ("peaks") in the vector V.
	%        MAXTAB and MINTAB consists of two columns. Column 1
	%        contains indices in V, and column 2 the found values.
	%      
	%        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
	%        in MAXTAB and MINTAB are replaced with the corresponding
	%        X-values.
	%
	%        A point is considered a maximum peak if it has the maximal
	%        value, and was preceded (to the left) by a value lower by
	%        DELTA.
	
	% Eli Billauer, 3.4.05 (Explicitly not copyrighted).
	% This function is released to the public domain; Any use is allowed.
	
	"""
	maxtab = []
	mintab = []
	   
	if x is None:
		x = np.arange(len(v))
	
	v = np.asarray(v)
	
	if len(v) != len(x):
		sys.exit('Input vectors v and x must have same length')
	
	if not np.isscalar(delta):
		sys.exit('Input argument delta must be a scalar')
	
	if delta <= 0:
		sys.exit('Input argument delta must be positive')
	
	mn, mx = np.Inf, -np.Inf
	mnpos, mxpos = np.NaN, np.NaN
	
	lookformax = True
	
	for i in np.arange(len(v)):
		this = v[i]
		if this > mx:
			mx = this
			mxpos = x[i]
		if this < mn:
			mn = this
			mnpos = x[i]
		
		if lookformax:
			if this < mx-delta:
				maxtab.append((mxpos, mx))
				mn = this
				mnpos = x[i]
				lookformax = False
		else:
			if this > mn+delta:
				mintab.append((mnpos, mn))
				mx = this
				mxpos = x[i]
				lookformax = True

	return np.array(maxtab), np.array(mintab)

