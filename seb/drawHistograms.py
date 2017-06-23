import numpy as np
from plot_support import *
from plot_functions import *
import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')

def drawHistograms(H_BINS, plotDir, inDir, inFile, MS=True, cyl=True, nodiff=False, show=False):
	f = ROOT.TFile(inDir + inFile)
	hList = openAllHistos(f)

	projMultiDir = plotDir + 'ProjectionMulti/'
	projSingleDir = plotDir + 'ProjectionSingle/'
	sliceMultiDir = plotDir + 'SliceMulti/'
	sliceSingleDir = plotDir + 'SliceSingle/'
	sliceDiffDir = plotDir + 'SliceDiff/'
	projDiffDir = plotDir + 'ProjectionDiff/'

	for i in [projMultiDir, projSingleDir, sliceMultiDir, sliceSingleDir, sliceDiffDir, projDiffDir]:
		makeDir(i)
		print i

	for key, h in hList.iteritems():
		print key, h
		# Projection
		plotHistoSliceMulti(h, H_BINS, True, cyl, show, projMultiDir + key + '.pdf')
		plotMaximum3d(h, H_BINS, True, True, cyl, show, projSingleDir + key + '.pdf')

		# Slice
		# try, because may fail for diff-data
		getMaximum(h, H_BINS)
		try:
			plotHistoSliceMulti(h, H_BINS, False, cyl, show, sliceMultiDir + key + '.pdf')
		except:
			pass

		try:
			plotMaximum3d(h, H_BINS, False, True, cyl, show, sliceSingleDir + key + '.pdf')
		except:
			pass

	# ==== DIFF PLOTS ====
	# First the projections onto a TH2 are made and in
	# succession their difference is plotted.
	# For the other difference plots, the difference of 
	# the TH3 histograms is plotted.
	if(nodiff):
		return
	else:
		if MS:
			name =  'MS'
			eventType = 'ms'
		else:
			name = 'SS'
			eventType = 'ss'

		dataHist = hList['dataHist%s' % name]
		mcHist = hList['mcHist%s' % name]

		# Slice
		plotDiffSliceMulti(dataHist, mcHist, H_BINS, False, cyl, getHistoDiffInv, '', False, sliceDiffDir + 'sliceGetHistoDiff%s.pdf' % name)
		plotDiffSliceMulti(dataHist, mcHist, H_BINS, False, cyl, getHistoRelativeDiff3d, '', False, sliceDiffDir + 'sliceGetHistoRelativeDiff%s.pdf' % name)

		# Projection
		plotDiffSliceMulti(dataHist, mcHist, H_BINS, True, cyl, getHistoDiffInv, '', False, projDiffDir + 'projectGetHistoDiff%s.pdf' % name)
		plotDiffSliceMulti(dataHist, mcHist, H_BINS, True, cyl, getHistoRelativeDiff3d, '', False, projDiffDir + 'projectGetHistoRelativeDiff%s.pdf' % name)

 	f.Close()

def drawHistSlices(H_BINS, plotDir, inDir, inFile, show=False):
    f = ROOT.TFile(inDir + inFile)
    hList = openAllHistos(f)
    print hList
    h = hList['diffSS']

    plotDir += 'zSlices/'
    makeDir( plotDir ) 

    zBins = H_BINS[6]
    zBins = 60
    zArray = np.linspace(0, -200, zBins/2. + 2)[1:-1]

    for z in zArray:
        try:
            print 'Showing detector at %f.' % z
            plotHistoSliceSingle(h, H_BINS, 'z', z, False, True, False, plotDir + 'SS_at_%s.pdf' % str( z ).replace('.', '_'), (-1., 1.))
        except:
            pass

