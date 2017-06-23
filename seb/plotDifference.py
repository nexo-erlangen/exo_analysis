import sys
import numpy as np

from plot_support import *
from plot_functions import *

import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')

# == FIDUCIAL VOLUME ==
APOTHEM = 200 #162
Z_MIN = 0 #5
Z_MAX = 192 #182

ROOT.EXOFiducialVolume.SetUserHexCut(APOTHEM, Z_MIN, Z_MAX)

# == ENERGY CUT ==
eMin = 700 #750
eMax = 3500

# == 1D HISTOGRAM PARAMETERS ==
H_BINS = 20
H_MIN = 0
H_MAX = 200

# == 3D HISTOGRAM PARAMETERS ==
H_BINS_ALL = 20
H_BIN_CART = [H_BINS_ALL, -200, 200, H_BINS_ALL, -200, 200, H_BINS_ALL, -200, 200]
# H_BIN_CYL = [H_BINS_ALL, 0, 180**2, H_BINS_ALL, -200, 200, H_BINS_ALL, -1, 1] 
H_BIN_CYL = [30, 0, 180**2, 40, -200, 200, 20, -1, 1]

# ==== MAIN ====
def main():
	'''
	if len(sys.argv):
		MC_FNAME = sys.argv[1]
		DATA_FNAME = sys.argv[2]
	else:
		print "Usage: MCFILE DATAFILE"

	# ==== TREES ====
	mcFile = ROOT.TFile.Open(MC_FNAME)
	mcTree = mcFile.Get('mcTree')
	dataFile = ROOT.TFile.Open(DATA_FNAME)
	dataTree = dataFile.Get('dataTree')
	'''

	preProcessDir = '/home/vault/capm/mppi025h/analysis/preprocess/'
	dataList = [preProcessDir + 'data_wThS5_denoised_camp%d.root' % i for i in range(3,6)] + [preProcessDir + 'prepro_ThS5_new_data_trial_camp%d.root' % i for i in range(1,6)]
	mcList = [preProcessDir + 'SourceS5_Th228.root']
	dataChain = getChain('dataTree', dataList)
	mcChain = getChain('mcTree', mcList)

	# ==== GET HIST DIFFERENCE ====
	# CYLINDRICAL
	hist3dCyl = ROOT.TH3D('hist3d', 'hist3d', *H_BIN_CYL)
	mcHistSS3dCyl = fillHisto3d(mcChain, hist3dCyl, getCut(True, True, 'ss', True, eMin, eMax), 'MC SS', True)
	dataHistSS3dCyl = fillHisto3d(dataChain, hist3dCyl, getCut(False, True, 'ss', False, eMin, eMax), 'Data SS', True)

	diffSS = getHistoRelativeDiff3d(mcHistSS3dCyl, dataHistSS3dCyl)
	#plotHistoSliceMulti(mcHistSS3dCyl, H_BIN_CYL, True, True, True)
	plotHisto2d(dataHistSS3dCyl.Project3D('yx'), '', '', '', True)
	plotHisto2d(mcHistSS3dCyl.Project3D('yx'), '', '', '', True)
	plotHisto2d(diffSS.Project3D('yx'), '', '', '', True)

	# plotHistoSliceMulti(dataHistSS3dCyl, H_BIN_CYL, True, True, True)
	#plotHistoSliceMulti(diffSS, H_BIN_CYL, True, True, True)

	# CARTESIAN
	'''
	hist3dCart = ROOT.TH3D('hist3d', 'hist3d', *H_BIN_CART)
	mcHistSS3dCart = fillHisto3d(mcChain, hist3dCart, getCut(True, True, 'ss', True, eMin, eMax), 'MC SS', False)
	dataHistSS3dCart = fillHisto3d(dataChain, hist3dCart, getCut(True, True, 'ss', False, eMin, eMax), 'Data SS', False)

	diffSS = getHistoRelativeDiff3d(mcHistSS3dCart, dataHistSS3dCart)
	diffSS2 = getHistoDiff(dataHistSS3dCart, mcHistSS3dCart)
	plotHistoSliceMulti(mcHistSS3dCart, H_BIN_CART, True, False, True)
	plotHistoSliceMulti(dataHistSS3dCart, H_BIN_CART, True, False, True)
	plotHistoSliceMulti(diffSS, H_BIN_CART, True, False, True)
	plotHistoSliceMulti(diffSS2, H_BIN_CART, True, False, True)
	'''

	'''
	mcProj = mcHistSS3dCyl.Project3D('yx')
	dataProj = dataHistSS3dCyl.Project3D('yx')
	mcProj.Scale(1./mcProj.Integral())
	dataProj.Scale(1./dataProj.Integral())
	#h = getHistoRelativeDiff2d(mcProj, dataProj)
	h = getHistoDiff(dataProj, mcProj)
	plotHisto2d(h, '', '', '', True)
	'''

def glOpacity(x, par):
	if(x[0] < 0.1*par[0]):
		return 0
	elif(x[0] < 0.3*par[0]):
		return 0.1
	elif(x[0] < 0.5*par[0]):
		return 0.3
	elif(x[0] < 0.8*par[0]):
		return 0.8
	else:
		return 1

if __name__ == '__main__':
	main()
