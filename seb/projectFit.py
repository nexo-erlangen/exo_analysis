#!/usr/bin/env python
import numpy as np
import importlib
import argparse

import plot_support as ps
import plot_functions as pf

import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')

REFLECTORINNERRAD = 183.2356 # mm
CATHODE_ANODE_y_DISTANCE = 192.23657 # mm
GAUSS_EXCLUDE = (-10, 10)
# GAUSS_EXCLUDE = None

def main():
	settings = get_args()
	loadSettings(settings)

	f = ROOT.TFile(SSoutDir + SSoutFile)
	hList = ps.openAllHistos(f)
	mcHist = hList['mcHistSS']

	# parList = getParList(mcHist)
	parList = (255.0, 64.96189303640774, 3.9, 18554.095519082686, -30.0, 104.85530141525176)
	posList, h = generateRandom(parList, 15.e3)

	h.Scale(1./h.Integral())

	fOut = ROOT.TFile('mcArt.root', 'RECREATE')
	h.Write()
	fOut.Close()

	pf.plotHistoSliceMulti(h, H_BINS, proj=True, cyl=False, show=True, out='')

	pf.plotHistoSliceSingle(h, H_BINS, option='x', value=0, proj2d=False, proj1d=True, show=True, out='')
	pf.plotHistoSliceSingle(h, H_BINS, option='y', value=0, proj2d=False, proj1d=True, show=True, out='')
	pf.plotHistoSliceSingle(h, H_BINS, option='z', value=0, proj2d=False, proj1d=True, show=True, out='')

	f.Close()
	# projectFit(H_BINS, SSoutDir, SSoutFile)
	return

def getParList(h):
	center = [255., 3.9, -30.]
	parList = gaus3dFit(h, center)[1:]
	mu1, sig1, mu2, sig2, mu3, sig3 = parList
	# sig1, sig2, sig3 = abs(sig1)*np.sqrt(2), abs(sig2)*np.sqrt(2), abs(sig3)*np.sqrt(2)

	return mu1, sig1, mu2, sig2, mu3, sig3

def generateRandom(parList, N):
	mu1, sig1, mu2, sig2, mu3, sig3 = parList

        H_BINS_ALL = 60
        H_BINS = [H_BINS_ALL, -200, 200, H_BINS_ALL, -200, 200, H_BINS_ALL, -200, 200]
	h = ROOT.TH3F('hist', 'hist', *H_BINS)
	num = 0
	
	ratio = 3172
	posList = []
	for i in range(int(N * ratio)):
		if not i % 100000:
			print '%d Events processed' % i
		x, y, z = np.random.normal(mu1, sig1, 1)[0], np.random.normal(mu2, sig2, 1)[0], np.random.normal(mu3, sig3, 1)[0]
		r = np.sqrt( x**2 + y**2 )

		# if not ps.isFiducial(x, y, z, 162, 0, 200) or r > REFLECTORINNERRAD or abs( z ) > 184:
                if r > REFLECTORINNERRAD or abs(z) > CATHODE_ANODE_y_DISTANCE:
			continue
		else:
			num += 1

		h.Fill(x, y, z, 1)
		posList.append( list(np.array([x, y, z])/1000.) )

	print 'Number of Events: %d' % num
	return posList, h

def gaus3d(x, par):
	xx = (x[0] - par[1]) / par[2]
	yy = (x[1] - par[3]) / par[4]
	zz = (x[2] - par[5]) / par[6]

	rx = ROOT.TMath.Exp(-0.5*xx**2)
	ry = ROOT.TMath.Exp(-0.5*yy**2)
	rz = ROOT.TMath.Exp(-0.5*zz**2)

	return par[0]*rx*ry*rz + par[7]

def gaus3dFit(h, center):
	c = center
	# func = ROOT.TF3('func', gaus3d, c[0]-200, c[0]+200, c[1]-200, c[1]+200, c[2]+200, c[2]-200, 7)
	func = ROOT.TF3('func', gaus3d, -160, 160, -70, 70, -100, 100, 8)
	func.SetParameters(h.GetMaximum(), c[0], 70, c[1], 70, c[2], 70, 0.01)
	func.FixParameter(1, c[0])
	func.FixParameter(3, c[1])
	func.FixParameter(5, c[2])

	ret = h.Fit(func, 'BSR+')
	parList = []
	for i in range(7):
		parList.append( ret.Parameter(i) )

	'''
	func.Draw('goff')
	h3 = ROOT.TH3F( func.GetHistogram() )
	h3.Draw()

	x, y, z = np.array([0.]), np.array([0.]), np.array([0.])
	for i in range(20):
		h3.GetRandom3(x, y, z)
		print x, y, z

	'''
	return parList

def gausExclude(x, par):
	if GAUSS_EXCLUDE:
		if x[0] > GAUSS_EXCLUDE[0] and x[0] < GAUSS_EXCLUDE[1]:
			ROOT.TF1.RejectPoint()
			return -999

	return par[0] * np.exp(-(x[0]-par[1])**2/(2*par[2]**2)) + par[3] * (x[0] - par[5]) + par[4]

def projectFit(H_BINS, inDir, inFile):
	# Read histograms from file
	f = ROOT.TFile(inDir + inFile)
	hList = ps.openAllHistos(f)

	mcHist = hList['mcHistSS'].Project3D('zy')
	dataHist = hList['dataHistSS'].Project3D('zy')

	mcHistY = mcHist.ProjectionY()
	dataHistY = dataHist.ProjectionY()

	# Set correct titles and axis labels
	mcHist.SetTitle('MC (zy-projection)')
	mcHist.GetXaxis().SetTitle('y / [mm]')
	mcHist.GetYaxis().SetTitle('z / [mm]')

	dataHist.SetTitle('Data (zy-projection)')
	dataHist.GetXaxis().SetTitle('y / [mm]')
	dataHist.GetYaxis().SetTitle('z / [mm]')

	# Fits
	gausExcludeFunc = ROOT.TF1('gausExclude', gausExclude, -180, 180, 6)

	gausExcludeFunc.SetParName(0, 'amplitude')
	gausExcludeFunc.SetParName(1, 'mean')
	gausExcludeFunc.SetParName(2, 'sigma')

	gausExcludeFunc.SetParName(3, 'm')
	gausExcludeFunc.SetParName(4, 't')
	gausExcludeFunc.SetParName(5, 's')

	# Estimate fit parameters
	m = ( mcHistY.GetBinContent(38) - mcHistY.GetBinContent(3) ) / ( (38 - 3) * 10 )

	print mcHistY.GetBinContent(38)
	print mcHistY.GetBinContent(3)
	t = mcHistY.GetBinContent(3) - 0.005

	gausExcludeFunc.SetParameter(0, 0.03)
	gausExcludeFunc.SetParameter(1, -30.)
	gausExcludeFunc.SetParameter(2, 70.)
	gausExcludeFunc.SetParameter(3, m)
	gausExcludeFunc.SetParameter(4, 0.)
	gausExcludeFunc.FixParameter(5, 0.015)

	# Plot & Fit
	c = ROOT.TCanvas()
	mcHistY.Draw('HIST')
	dataHistY.Draw('HIST same')

	print '== MC FIT ==\n============'
	mcHistY.Fit('gausExclude', 'S same', '', -170, 170)
	gausExcludeMC = gausExcludeFunc.Clone()
	gausExcludeMC.Draw('lsame')
	print

	print '== DATA FIT ==\n=============='
	dataHistY.Fit('gausExclude', 'S same', '', -170, 170)
	gausExcludeFunc.Draw('lsame')
	print 

	c.Update()
	raw_input('end')

# === SUPPORT ===
def loadSettings(settings):
	settingsLib = importlib.import_module(settings.split('.')[0])
	settingsMod = settingsLib.__dict__
	try:
		to_import = settingsLib.__all__
	except AttributeError:
		to_import = [name for name in settingsMod if not name.startswith('_')]
	globals().update({name: settingsMod[name] for name in to_import})

def get_args():
	ap = argparse.ArgumentParser(description=' ')

	# Create histogram
	ap.add_argument('-s', '--settings', type=str, help='Settings to load', required=True)
	args = ap.parse_args()

	return args.settings

if __name__ == '__main__':
	main()

