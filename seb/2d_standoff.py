import numpy as np
import sys

import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')

# important stuff
APOTHEM = 200 #162
Z_MIN = 0 #5
Z_MAX = 192 #182

# MC_FNAME = sys.argv[1]
# DATA_FNAME = sys.argv[2]

# set fiducial volume
ROOT.EXOFiducialVolume.SetUserHexCut(APOTHEM, Z_MIN, Z_MAX)

# r parameters
rBins = 40 
rMin = 0
rMax = 180

# z parameters
zBins = 40
zMin = -200
zMax = 200

# energy cut
eMin = 750
eMax = 3500

def main():
	mcHist = ROOT.TH2D('mc_cluster_hist', 'mc_cluster_hist', rBins, rMin, rMax, zBins, zMin, zMax)
	dataHist = ROOT.TH2D('data_cluster_hist', 'data_cluster_hist', rBins, rMin, rMax, zBins, zMin, zMax)

	preProcessDir = '/home/vault/capm/mppi025h/analysis/preprocess/'
	dataList = [preProcessDir + 'data_wThS5_denoised_camp%d.root' % i for i in range(3,6)] + [preProcessDir + 'prepro_ThS5_new_data_trial_camp%d.root' % i for i in range(1,6)]
	mcList = [preProcessDir + 'pre_1m_Th_real_diff_pre.root']
	dataTree = getChain('dataTree', dataList)
	mcTree = getChain('mcTree', mcList)

	'''
	mcFile = ROOT.TFile(MC_FNAME)
	dataFile = ROOT.TFile(DATA_FNAME)

	mcTree = mcFile.Get('mcTree')
	dataTree = dataFile.Get('dataTree')
	'''

	fillHisto(mcTree, mcHist, True, 'ss')
	fillHisto(dataTree, dataHist, False, 'ss')

	mcHist.Scale(1./mcHist.Integral())
	dataHist.Scale(1./dataHist.Integral())

	plotHisto(mcHist)
	plotHisto(dataHist)

	diffHist = dataHist.Clone()
	diffHist.Add(mcHist, -1)

	errorHist = diffHist.Divide(diffHist, mcHist)
	plotHisto(diffHist)

def fillHisto(Tree, Hist, MC=True, type='ms'):
	nEvents = Tree.GetEntries()
	NMOD = int(nEvents/10)

	for i in range(nEvents):
		Tree.GetEntry(i)
		es = Tree.EventSummary
		if(not (i % NMOD)):
			print 'Processed %d events' % i

		defaultCuts = es.isMissingPosition or es.multiplicity <= 0 or es.nsc != 1
		
		if(type == 'ms'):
			typeCut = abs(es.multiplicity - 1) < 0.1
			enCut = es.energy_ms > eMax or es.energy_ms < eMin
		else:
			typeCut = es.multiplicity > 1.1
			enCut = es.energy_ss > eMax or es.energy_ss < eMin

		if(MC):
			enCut = es.energy_mc > eMax or es.energy_mc < eMin
			dataCuts = False
		else:
			dataCuts = False #es.isNoise or (not es.isWithinDriftTime) or es.isDiagonallyCut() or es.isSolicitedTrigger

		sdCut = False # (es.standoff_distance < 0 or es.standoff_distance > 30) 
		fidCut = not es.isFiducial()

		# apply cuts
		if(defaultCuts or enCut or typeCut or dataCuts or sdCut or fidCut):
			continue

		eventPosX = np.array(es.cluster_x)
		eventPosY = np.array(es.cluster_y)
		eventPosZ = np.array(es.cluster_z)
		eventEnergy = np.array(es.cluster_energy)

		numCluster = len(eventPosX)
		if(numCluster <= 0 or numCluster != len(eventPosY) or numCluster != len(eventPosZ) or numCluster != len(eventEnergy)):
			continue

		# rPos = (eventPosX**2 + eventPosY**2)**(1./2)
		# rMaxEvent = np.amax(rPos * (rPos < 900.))
		rMaxEvent = es.GetMaxR()
		zMaxEvent = np.amax(eventPosZ * (abs(eventPosZ) < 900.))
		zMinEvent = np.amin(eventPosZ * (abs(eventPosZ) < 900.))
		if abs(zMaxEvent) < abs(zMinEvent):
			zMaxEvent = zMinEvent

		Hist.Fill(rMaxEvent**2, zMaxEvent)

def plotHisto(Hist):
	c = ROOT.TCanvas()
	Hist.SetStats(ROOT.kFALSE)
	Hist.SetTitle("Energy %i - %i [keV]" % (int(eMin), int(eMax)))
	Hist.GetXaxis().SetTitle("r [mm]")
	Hist.GetYaxis().SetTitle("z [mm]")
	Hist.Draw("")

	raw_input('end')

def getChain(treeName, fileList):
	ch = ROOT.TChain(treeName)
	if(not ch): return False

	for item in fileList:
		ch.Add(item)

	return ch

if __name__ == '__main__':
	main()
