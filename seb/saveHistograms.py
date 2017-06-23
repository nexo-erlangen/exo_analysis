import sys
import numpy as np

from plot_support import *
from plot_functions import *

import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')
ROOT.PyConfig.IgnoreCommandLineOptions = True

# ==== saveHistogramSS ====
def saveHistograms(H_BINS, FV, EN, preProcessDir, dataList, mcList, outDir, outFile, MS=True, nodiff=False, cyl=True):
	ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
	# ROOT.EXOFiducialVolume.SetUserRadialCut(*FV)

	dataChain = getChain('dataTree', dataList)
	mcChain = getChain('mcTree', mcList)

	# ==== GET HIST DIFFERENCE ====
	hist3d = ROOT.TH3D('hist3d', 'hist3d', *H_BINS)
	dict = {}

	if MS:
		name =  'MS'
		eventType = 'ms'
	else:
		name = 'SS'
		eventType = 'ss'

	dict['mcHist%s' % name] = fillHisto3d(mcChain, hist3d, getCut(True, True, eventType, True, *EN), 'mcHist%s' % name, cyl, True)
	dict['dataHist%s' % name] = fillHisto3d(dataChain, hist3d, getCut(False, True, eventType, False, *EN), 'dataHist%s' % name, cyl, True)
	if not nodiff:
		dict['diff%sRel' % name] = getHistoRelativeDiff3d(dict['dataHist%s' % name], dict['mcHist%s' % name], 'diff%sRel' % name)
		dict['diff%s' % name] = getHistoDiffInv(dict['dataHist%s' % name], dict['mcHist%s' % name], 'diff%s' % name)

	# Save histograms to file
	saveHistos(outDir + outFile, dict)
	