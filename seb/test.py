#!/usr/bin/env python

# ==========================================
# File is used to create histograms used in
# the z-Standoff Plot for drifted MC vs 
# real Data using the command
# ./test.py -s <settings> -g
#
# Plot is shown using the command
# ./test.py -sh
# ==========================================

import argparse
import os
import importlib

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import sortData as sd
import plot_support as ps

def main():
	settings, generate, show = get_args()
	loadSettings(settings)

	if generate:
		t = ROOT.TFile.Open('driftData.root', 'READ')

		mcTree = t.Get('mcTree')

		# dataList is defined in the settings file
		dataChain = ps.getChain('dataTree', dataList)

		sd.sortDataSave(dataChain, mcTree, FV, 100, 'StandoffData')

	if show:
		hMCDict, hDataDict = sd.sortDataFile('mcStandoffData.root', 'dataStandoffData.root', nBins)
		hMCDictZ, hDataDictZ = sd.sortDataFile('mcStandoffDataZ.root', 'dataStandoffDataZ.root', nBins)
		hMCDictTheta, hDataDictTheta = sd.sortDataFile('mcStandoffDataTheta_z.root', 'dataStandoffDataTheta_z.root', nBins)

		xGrid, yGrid, valMCNorm, valDataNorm, valDiffNorm, valDiffNormError = sd.plotStandoff(hMCDictZ, hDataDictZ, nBins, showBin)
		print 'Creating valDiffShow...'
		valDiffShow = np.array(valDiffNorm)[:,showBin]

		print 'Creating valDiffErrorShow...'
		print len(valDiffNorm), len(valDiffNormError), len(valDiffNorm[0]), len(valDiffNormError[0])
		valDiffErrorShow = np.array(valDiffNormError)[:,showBin]

		print 'Done...'

		for i in range( len(valDiffShow) ):
			print valDiffShow[i], valDiffErrorShow[i]

		peak.findPeaks(xGrid, valDiffShow, valDiffErrorShow, True, '(data - mc)/data for Bin #0')

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
	
	ap.add_argument('-s', '--settings', type=str, help='Settings to load', required=True)
	ap.add_argument('-g', '--generate', help='Generate histograms', action='store_true', default=False)
	ap.add_argument('-sh', '--show', help='Show z-Standoff plot', action='store_true', default=False)

	args = ap.parse_args()

	if not args.generate and not args.show:
		ap.error('Either --generate or --show or both have to be selected.')

	return args.settings, args.generate, args.show

if __name__ == '__main__':
	main()

