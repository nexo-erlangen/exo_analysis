#!/usr/bin/env python

import argparse
import os
import importlib
import subprocess

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gStyle.SetPalette(55)

import time
from saveHistograms import *
from drawHistograms import *
from plot_support import *
from drawHistoCyl import *
from generate_random import *
import sortData as sd
import standoff_comparison as sc
import vtk_to_vec as vtv
import vtk_to_potential as vtp
import peak_finder as peak
import projectFit as fit
import plot_functions as pf
import gen_ran_qsub as grq

# from settingsQxu import *

def main():
        start_time = time.time()

	settings, histogram, ms, ss, plot, show, cart, nodiff, drawcyl, random, standoff, test, energy = get_args()
	loadSettings('settings/', settings)
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
        if standoff:
		nBins = 120
		# showBin = 5 shows sum over first five bins
                showBin = 0

                sortDataOutFn = settings.split('.')[0]
                sortDataOutDir = '/home/vault/capm/mppi025h/plot_scripts/standoffFiles/'

		# == GENERATE PLOT DATA ==
		if not show:
			# Currently wo cut!
			# But using this data is the better choice due to much
			# higher statistics in bin #0
			
                        if not random:
				dataChain = getChain('dataTree', dataList)
				mcChain = getChain('mcTree', mcList)

				# Fill standoff histograms and save them to
				# root files.

                                if not os.path.isdir(sortDataOut + sortDataOutFn):
                                    os.makedirs(sortDataOut + sortDataOutFn)

				if ms:
					sd.sortDataSave(dataChain, mcChain, FV, nBins, sortDataOutFn + 'MS', 'ms', sortDataOutDir)
				if ss:
					sd.sortDataSave(dataChain, mcChain, FV, nBins, sortDataOutFn + 'SS', 'ss', sortDataOutDir)
			
                        else:
				# Random data can only be used with SS evenst!
				# Positions wo cut
				print 'Processing event positions...'
				# mergeLists(getFileList('drift_data', 'randPos_'), 'randPos.p', True)
				# mergeLists(getFileList('drift_data', 'randDriftPos_'), 'randDriftPos.p', True)

				sd.sortDataSaveRand('randPos.p', 'randDriftPos.p', nBins, 'Standoff')
				# Positions w cut
				# sd.sortDataSaveRand('randPosCut.p', 'randDriftPosCut.p', nBins)

			return

		# == PLOT ==
		else:
			if not random:
				# outDir = '/home/vault/capm/mppi025h/plot_scripts/NMExoAnalysis/'
				# outDir = './'
				if ms:
					# sd.getPlot(False, nBins, showBin, 'StandoffMS', outDir)
					sd.getPlot(False, nBins, showBin, sortDataOutFn + 'MS', sortDataOutDir)
				if ss: 
					# sd.getPlot(False, nBins, showBin, 'StandoffSS', outDir)
					sd.getPlot(False, nBins, showBin, sortDataOutFn + 'SS', sortDataOutDir)

			else:
				sd.getPlot(True, nBins, showBin, 'Standoff')

			return

	if drawcyl:
                # fName = 'standoffS5fid.root'
                # fName = 'standoffHex_artDrift_neg.root'
                # fName = 'standoffArtDriftNM.root'
                fName = sortDataOutFn + '.root'

		# ==== GET TRUE STANDOFF ====
		# Store standoff histograms of investigated
		# data to root file. 
		# ATTENTION: Keep the fiducial cut in mind!
		dataChain = getChain('dataTree', dataList)
		mcChain = getChain('mcTree', mcList)

                # scatterEnergy(dataChain, FV, 'ss')
                # return

                # Compare characteristics of MC and data root-files
                # getCompHisto(fName, dataChain, mcChain, FV, type='ss', selectList=['standoff_distance', 'multiplicity', 'energy', 'num_coll_wires', 'num_ind_wires', 'u_mst_metric', 'v_mst_metric'], name='')
                # getCompHisto(fName, dataChain, mcChain, FV, type='ms', selectList=['standoff_distance', 'multiplicity', 'energy', 'num_coll_wires', 'num_ind_wires', 'u_mst_metric', 'v_mst_metric'], name='')

                # getCompHisto(fName, dataChain, mcChain, FV, type='ss', selectList=['standoff_distance', 'multiplicity', 'energy', 'u_mst_metric', 'v_mst_metric'], name='')
                # getCompHisto(fName, dataChain, mcChain, FV, type='ms', selectList=['standoff_distance', 'multiplicity', 'energy', 'u_mst_metric', 'v_mst_metric'], name='')

                # getCompHisto(fName, dataChain, mcChain, FV, type='ss', selectList=['energy'], name='')
                # return

		if ss:
			getStandoffHisto('standoffRoot/' + fName, dataChain, FV, 'ss', False, 'DataSS')
			# NEvents = getApothemThetaHisto('standoffRoot/' + fName, dataChain, FV, 'ss', False, 0, 'DataSS', 3)
			# getApothemHisto('standoffRoot/' + fName, dataChain, FV, 'ss', False, 'DataSS')
                        # print 'NEvents:', NEvents

                        if random:
                            # Rate agreement
                            # NEvents = int( 1.e6 )
                            generate2nbb('standoffRoot/' + fName, NEvents, FV, name='McSS', N=3)
                        else:
                            getStandoffHisto('standoffRoot/' + fName, mcChain, FV, 'ss', True, 'McSS')
                            # getApothemThetaHisto('standoffRoot/' + fName, mcChain, FV, 'ss', True, NEvents, 'McSS', 3)
                            # getApothemHisto('standoffRoot/' + fName, mcChain, FV, 'ss', True, 'McSS')
		if ms:
			getStandoffHisto('standoffRoot/' + fName, dataChain, FV, 'ms', False, 'DataMS')
			# getApothemThetaHisto('standoffRoot/' + fName, dataChain, FV, 'ms', False, 'DataMS')

		    	getStandoffHisto('standoffRoot/' + fName, mcChain, FV, 'ms', True, 'McMS')
		    	# getApothemThetaHisto('standoffRoot/' + fName, mcChain, FV, 'ms', True, 'McMS')

                # plotLbkg('standoffRoot/' + fName, fit=True, art='ss', output='hexSplit6_lbkgCut.pdf', N=3)

		# ==== READ & COMPARE STANDOFF ====
                if not os.path.isdir('standoffPlots/%s' % sortDataOutFn):
                    os.makedirs('standoffPlots/%s' % sortDataOutFn)

                if ss:
			plotStandoffHisto('standoffRoot/' + fName, 'ss', False, 'standoffPlots/%s' % sortDataOutFn + fName.split('.')[0] + 'SS.pdf')
                if ms:
			plotStandoffHisto('standoffRoot/' + fName, 'ms', False, 'standoffPlots/' + fName.split('.')[0] + 'MS.pdf')

		'''
		# ==== COMPARE DRIFTED RND DATA ====
		spacing = [np.linspace(-0.0000001, -0.000001, 10), np.linspace(-0.2, 0.2, 10), np.linspace(-0.05, 0.05, 10)]

		fName = 'standoffDrift.root'
		fNameTrue = 'trueStandoff.root'
		rootToHist(fName, fNameTrue, spacing, 5)

		# drawHistoCyl(diffRelHist)
		'''

	if(random):
		if ss:
			f = ROOT.TFile(SSoutDir + SSoutFile)
		else:
			f = ROOT.TFile(MSoutDir + MSoutFile)

		mcHist, dataHist, diffHist, diffRelHist = getHist(f)
		
		# Clear root file if it exists.
		# Later used to store the resulting standoff
		# distance histograms from drift
		#if os.path.exists('standoffDrift.root'):
		#	os.remove('standoffDrift.root')

		# FV is the fiducial volume set for the 
		# histogram. If random is executed, the apothem
		# should be set to a large value, so no radial
		# cut is performed. Instead, the cut is done
		# afterwards according to the specified apothem
		# ApoRand. Keep in mind that no z-cut is 
		# performed here, so it has to be set via FV.
		ApoRand = [162./1000, 5./1000, 182./1000]
		# ATTENTION: Currently z-value is unused!

		nParticles = 1000
		outName = 'standoffDrift.root'

		ROOT.gROOT.cd()
		randTree = ROOT.TTree('random', 'random')
		dumpTree = ROOT.TTree('dump', 'dump')

		# Generate random positions:
		# randPos used to generate normal distribution
		# randTree used to get drifted distribution

		print 'Random positions:'
		randPos = generateRandomArt(mcHist, dumpTree, nParticles)
		print 'Random drifted positions:'
		generateRandomArt(mcHist, randTree, nParticles)
		del dumpTree

		sfepyDir = '/home/hpc/capm/mppi025h/sfepy'
		#sigma_0 = -0.001 # C/m**2
		sigma_0 = -0.0000001
		theta = 0.
		z = -0.037
		ra = 0.085
		rb = 0.075

		for sigma_0 in [0]: #np.linspace(-1.e-6, -1.e-5, 10):
			# for theta in np.linspace(-0.2, 0.2, 9):
				# for z in np.linspace(-0.13, -0.05, 9):

					sigma_0 = 0.
					# sigma_0 = -1.e-5
					# theta = 0.
					# z = -0.024

					name = 'sigma0: %0.2e, theta: %0.2f, z: %0.3f' % (sigma_0, theta, z)
					print 'Processing %s...' % name
					
					# write settings for sfepy
					#writeSfepySettings(sigma_0, theta, z, ra, rb, sfepyDir)
					
					# writes exo200_wo_edges.vtk in dir the
					# process is executed
					
					#subprocess.call([ 'python', sfepyDir + '/simple.py', 'exo_field_shape_fem.py' ])
					# os.system('python %s/simple.py %s' % (sfepyDir, 'exo200_fem_local.py'))

					inFile = 'exo200_wo_edges.vtk'
					# Writes efield_hist.root in sfepyDir
					#vtv.storeEField('.', inFile, H_BINS_RAND, False)

					os.system('free -m')

					# Reads the efield file and drifts randomly
					# generated particles. In conclusion, plots
					# the standoff distance of the initial and
					# drifted particles.
					#try:
					efieldName = '/home/vault/capm/mppi025h/g/exo_data6/groups/3DFieldMaps/3Dmaxwell/3D_Efield.root'
					H_BINS_RAND = None

					randomProcess(randPos, randTree, H_BINS_RAND, H_BINS, ApoRand, name, nParticles, outName, True, efieldName, True, inFile, True, False)
					#except:
					#	print
					#	print 'Failed'

		f.Close()

        if energy:
                import energy_smear as es
                dataChain = getChain('dataTree', dataList)
                mcChain = getChain('mcTree', mcList)
                if ss:
                    art = 'ss'
                elif ms:
                    art = 'ms'

                es.energySmear(dataChain, dataList, mcChain, FV, art)
                return

	if test:
                driftEvaluation(30, 1, True)
                return

                # Plot residual (data - mc)/mc as a function of angle theta
                # plotZtheta('dataStandoffSSTheta_z.root', 'mcStandoffSSTheta_z.root')

                # Draw slices of the detector volume and look at the residual
                # of data and mc
                # drawHistSlices(H_BINS, plotDir + 'SS/', SSoutDir, SSoutFile, True)

                # Randomly generate particles and save the resulting field lines
		# randomTest(H_BINS_RAND)
		randomTest(None, True)
		
		# pf.projectFit(H_BINS, SSoutDir, SSoutFile)

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
        print '---- Execution time: %s s ----' % (time.time() - start_time)

# ==== FUNCTIONS ====
def loadSettings(path, settings):
        sys.path.append(path)
	settingsLib = importlib.import_module(settings.split('.')[0])
	settingsMod = settingsLib.__dict__
	try:
		to_import = settingsLib.__all__
	except AttributeError:
		to_import = [name for name in settingsMod if not name.startswith('_')]
	globals().update({name: settingsMod[name] for name in to_import})

def writeSfepySettings(sigma_0, theta, z, ra, rb, sfepyDir):
	epsilon_0 = 8.854187817620e-12 # F/m
	charge = sigma_0/epsilon_0	# this is sigma

	# Circular charge density
	R = 0.1832356
	#if theta < 0:
	#	qx = 2*np.pi*R - theta*R
	#else:
	qx = theta*R
	qy = z
	qr = ra
	m = float(ra)/rb

	f = open(sfepyDir + '/settings.py', 'w')
	f.write('charge = %f\n' % charge)
	f.write('R = %f\n' % R)
	f.write('qx = %f\n' % qx)
	f.write('qy = %f\n' % qy)
	f.write('qr = %f\n' % qr)
	f.close()

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
	ap.add_argument('-so', '--standoff', help='Generate standoff distance plots', action='store_true')
	ap.add_argument('-t', '--test', help='Access to test-area', action='store_true')
        ap.add_argument('-e', '--energy', help='Energy plots', action='store_true')

	args = ap.parse_args()

	if not args.histogram and not args.plot and not args.energy and not args.drawcyl and not args.random and not args.standoff and not args.test:
		ap.error('Either --histogram, --energy, --drawcyl or --plot has to be executed.')
	if not args.multisite and not args.singlesite:
		ap.error('Either --multisite or --singlesite or both have to be selected.')

	return args.settings, args.histogram, args.multisite, args.singlesite, args.plot, args.show, args.cart, args.nodifference, args.drawcyl, args.random, args.standoff, args.test, args.energy

if __name__ == '__main__':
	main()

