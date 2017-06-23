import argparse
import os
import importlib

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gStyle.SetPalette(55)

from saveHistograms import *
from drawHistograms import *
from plot_support import *
from drawHistoCyl import *
from generate_random import *
import standoff_comparison as sc

# from settingsQxu import *

def main():
	settings, histogram, ms, ss, plot, show, cart, nodiff, drawcyl, random = get_args()
	loadSettings(settings)
	cyl = not cart

	# ==== SAVE HISTOGRAM ====
	if(histogram):
		print 'Saving histograms:'
		print '=================='
		if(ms):
			saveHistograms(H_BINS, FV, EN, preProcessDir, dataList, mcList, MSoutDir, MSoutFile, True, nodiff, cyl)
			print
		if(ss):
			saveHistograms(H_BINS, FV, EN, preProcessDir, dataList, mcList, SSoutDir, SSoutFile, False, nodiff, cyl)
			print

	# ==== PLOT HISTOGRAM ====
	if(plot):
		print 'Plotting histograms...'
		if(ms):
			makeDir(plotDir + 'MS/')
			drawHistograms(H_BINS, plotDir + 'MS/', MSoutDir, MSoutFile, True, cyl, nodiff, show)
		if(ss):
			makeDir(plotDir + 'SS/')
			drawHistograms(H_BINS, plotDir + 'SS/', SSoutDir, SSoutFile, False, cyl, nodiff, show)

	# ==== TEST RANGE ====
	if ms:
		f = ROOT.TFile(MSoutDir + MSoutFile)
	else:
		f = ROOT.TFile(SSoutDir + SSoutFile)

	mcHist, dataHist, diffHist, diffRelHist = getHist(f)

	dataChain = getChain('dataTree', dataList)
	mcChain = getChain('mcTree', mcList)

	sc.compStandoff2d(dataChain, mcChain, bins=40)
	#sc.compStandoff(dataChain, mcChain, 'Theta', bins=40)

	if(drawcyl):
		drawHistoCyl(diffRelHist)

	if(random):
		# FV is the fiducial volume set for the 
		# histogram. If random is executed, the apothem
		# should be set to a large value, so no radial
		# cut is performed. Instead, the cut is done
		# afterwards according to the specified apothem
		# ApoRand. Keep in mind that no z-cut is 
		# performed here, so it has to be set via FV.
		ApoRand = 162./1000
		# randomProcess(mcHist, H_BINS_RAND, H_BINS, ApoRand, int(10000))

		'''
		ApoRand = 1
		h = ROOT.TH2D('h', 'h', 1000, -1.5, 1.5, 1000, -1.5, 1.5)
		import random
		for i in range(1000000):
			x = random.uniform(-2, 2)
			y = random.uniform(-2, 2)
			if isFiducial(x, y, ApoRand):
				h.Fill(x, y)

		h.Draw('colz')
		raw_input('end')
		'''

# ==== FUNCTIONS ====
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
	ap.add_argument('-hist', '--histogram', help='Generate histograms. Needs to be called before plotting', required=False, action='store_true', default=False)
	ap.add_argument('-ms', '--multisite', help='Process MS events', required=False, action='store_true', default=False)
	ap.add_argument('-ss', '--singlesite', help='Process SS events', required=False, action='store_true', default=False)
	ap.add_argument('-p', '--plot', help='Plot histograms', required=False, action='store_true', default=False)
	ap.add_argument('-sh', '--show', help='Show plots', required=False, action='store_true', default=False)
	ap.add_argument('-c', '--cart', help='Use cartesian coordinates', action='store_true', default=False)
	ap.add_argument('-nd', '--nodifference', help='Make no difference-plots', action='store_true', default=False)
	ap.add_argument('-dc', '--drawcyl', help='Do plots in cylindrical coordinates', action='store_true', default=False)
	ap.add_argument('-r', '--random', help='Generate random data from MC', action='store_true')

	args = ap.parse_args()

	if not args.histogram and not args.plot and not args.drawcyl and not args.random:
		ap.error('Either --histogram, --drawcyl or --plot has to be executed.')
	if not args.multisite and not args.singlesite:
		ap.error('Either --multisite or --singlesite or both have to be selected.')

	return args.settings, args.histogram, args.multisite, args.singlesite, args.plot, args.show, args.cart, args.nodifference, args.drawcyl, args.random

if __name__ == '__main__':
	main()
