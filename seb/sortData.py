import numpy as np
import re
import cPickle

import plot_support as ps
import plot_functions as pf
import peak_finder as peak
import generate_random as gr

import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')

# ==== CONSTANTS ====
CATHODE_ANODE_y_DISTANCE = 192.23657 # mm
REFLECTORINNERRAD = 183.2356 # mm

def sortDataSave(dataTree, mcTree, FV, bins=120, name='Standoff', art='ss', outDir='./'):
	ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
	ROOT.gROOT.cd()

	mcTreeCut = mcTree.CopyTree( ps.getCut(calibCut=True, energyCut=True, type=art, MC=True, sideCut=False), 'fast')
	dataTreeCut = dataTree.CopyTree( ps.getCut(calibCut=True, energyCut=True, type=art, MC=False, sideCut=False), 'fast')

	mcStandoff, mcStandoffZ, mcStandoffTheta = getSortData(mcTreeCut, 'mc', bins)
	dataStandoff, dataStandoffZ, dataStandoffTheta = getSortData(dataTreeCut, 'data', bins)

	saveToFile(mcStandoff, outDir + 'mc%s.root' % name)
	saveToFile(dataStandoff, outDir + 'data%s.root' % name)
	saveToFile(mcStandoffZ, outDir + 'mc%sZ.root' % name)
	print outDir + 'data%sZ.root' % name
	saveToFile(dataStandoffZ, outDir + 'data%sZ.root' % name)
	saveToFile(mcStandoffTheta, outDir + 'mc%sTheta_z.root' % name)
	saveToFile(dataStandoffTheta, outDir + 'data%sTheta_z.root' % name)

	# plotStandoff(mcStandoff, dataStandoff, bins)

def sortDataSaveRand(fName, fNameDrift, bins=120, name='Standoff'):
	standoff, standoffZ, standoffTheta = getSortData(fName, 'random', bins, True)
	driftStandoff, driftStandoffZ, driftStandoffTheta = getSortData(fNameDrift, 'randomDrift', bins, True)

	saveToFile(standoff, 'rand%s.root' % name)
	saveToFile(standoffZ, 'rand%sZ.root' % name)
	saveToFile(standoffTheta, 'rand%sTheta.root' % name)
	saveToFile(driftStandoff, 'randDrift%s.root' % name)
	saveToFile(driftStandoffZ, 'randDrift%sZ.root' % name)
	saveToFile(driftStandoffTheta, 'randDrift%sTheta.root' % name)

def sortDataFile(fNameMC, fNameData, bins):
	fMC = ROOT.TFile.Open(fNameMC)
	fData = ROOT.TFile.Open(fNameData)

	print fNameMC
	print fNameData
	hMC = ps.openAllHistos(fMC)
	hData = ps.openAllHistos(fData)
	
	hMCDict = correctDictNotation(hMC)
	hDataDict = correctDictNotation(hData)

	return hMCDict, hDataDict

def correctDictNotation(hDict):
	hCorrDict = {}
	for key in hDict:
		name = tuple( [float(item) for item in re.findall(r"[-+]?\d*\.\d+|\d+", key)] )
		hCorrDict[name] = hDict[key]

	return hCorrDict

def plotStandoff(mcStandoff, dataStandoff, bins, binN=0, binning=True, show=True, errorOut=False):
	xGrid = []
	yGrid = []
	# chi2 = []

	# Sum of all histograms
	hSumMC = ROOT.TH1D('sumMC', 'sumMC', 20, 0, 200)
	hSumData = ROOT.TH1D('sumData', 'sumData', 20, 0, 200)
	
	# hSumMC.Sumw2(True)
	# hSumData.Sumw2(True)

	# Bins of data and MC histograms
	valMC, valMCError = [], []
	valData, valDataError = [], []
	# Bins of the difference of data and MC histograms
	valDiff = []
	valDiffError = []

	Chi2 = []

	# ROOT.TH1.Sumw2(True)

	# Loop over all standoff distance histograms.
	# mcStandoff and dataStandoff have the same keys
	for key in mcStandoff:
		# Get the histograms from dictionaries
		mcH = mcStandoff[key]
		dataH = dataStandoff[key]
		
		# mcH.Sumw2(True)
		# dataH.Sumw2(True)

		'''
		if binning:
			# Scale the histograms
			try:
				mcH.Scale(1./mcH.Integral())
				dataH.Scale(1./dataH.Integral())
			except:
				pass
		'''

		# Add the histograms to get their total
		# sum after the loop is done
		hSumMC.Add(mcH)
		hSumData.Add(dataH)

		# Store the coordinates which are part 
		# of the key. Later used to do the plot
		if len(key) == 2:
			x, y = key
			xGrid.append(x)
			yGrid.append(y)
		else:
			x = key
			xGrid.append(x)
	
		# Get the distance of the histograms
		h = dataH.Clone()
		h.Add(mcH, -1)
		h.Divide(dataH)

		hEntries = dataH.GetEntries() - mcH.GetEntries()

		'''
		# Show the standoff histograms for
		# the position specified in key
		c = ROOT.TCanvas()
		mcH.Draw('hist')
		dataH.SetMarkerStyle(3)
		dataH.Draw('same P')
		raw_input('end')
		del c
		'''

		#pf.plotCompQQ(mcH, dataH, [20,0,200], 'mc', 'data', show=True, out='')
		
		# == CHI2 ==
		# res = np.array( [0]*20 )
		# mcH.Chi2Test(dataH, 'UU', res)
		# res[5] = sum( res[0:5] )
		# Chi2.append( res )
		
		# Loop over the first bins of the 
		# standoff histograms and store their
		# values in the val* lists
		vData, vDataError = [], []
		vMC, vMCError = [], []
		vDiff = []
		vDiffError = []
		for i in range(1, 10 + 1):
			print 'hEntries', h.GetEntries()
			dataBin = dataH.GetBinContent(i)
			mcBin = mcH.GetBinContent(i)

			try:
				dataBinErr = np.sqrt( dataBin ) / dataH.Integral()
			except:
				dataBinErr = 0

			try:
				mcBinErr = np.sqrt( mcBin ) / mcH.Integral()
			except:
				mcBinErr = 0
		
			# TEST
			try:
				dataBin /= dataH.Integral()
			except:
				pass
			try:
				mcBin /= mcH.Integral()
			except:
				pass
			
			# dataBinErr = np.sqrt( dataBin * (1 - dataBin/dataH.GetEntries()) )
			# mcBinErr = np.sqrt( mcBin * (1 - mcBin/mcH.GetEntries()) )

			print 'Bin#:', i
			print 'dataBin =', dataBin
			print 'mcBin =', mcBin
			print 'dataBinErr =', dataBinErr
			print 'mcBinErr =', mcBinErr

			# dataBinErr = dataH.GetBinError(i)
			# mcBinErr = mcH.GetBinError(i)
			try:
				diffBin = float(dataBin - mcBin)/(dataBin)
				diffBinError = abs(mcBin/dataBin) * np.sqrt( (mcBinErr/mcBin)**2 + (dataBinErr/dataBin)**2 )

			except:
				diffBin = 0
				diffBinError = 0

			print 'diffBin =', diffBin
			print 'diffBinError = ', diffBinError
			print

			vDiff.append( diffBin )
			vDiffError.append( diffBinError )
			vData.append( dataBin )
                        vDataError.append( dataBinErr )
			vMC.append( mcBin )
                        vMCError.append( mcBinErr )

		# Also append the sum of all bins as
		# last entry
		vData.append( sum(vData) )
		vMC.append( sum(vMC) )
		vDiff.append( sum(vDiff) )

		# Store the lists containing the single
		# bins in the lists containing the data
		# for every position
		valData.append( vData )
                valDataError.append( vDataError )
		valMC.append( vMC )
                valMCError.append( vMCError )
		valDiff.append( vDiff )
		valDiffError.append( vDiffError ) 

		# Delete the histograms to free memory
		del h
		del mcH
		del dataH

	# Sum over all histograms
	'''
	hSumMC.Scale(1./hSumMC.Integral())
	hSumData.Scale(1./hSumData.Integral())
	c1 = ROOT.TCanvas()
	hSumMC.Draw('HIST')
	hSumData.SetMarkerStyle(3)
	hSumData.Draw('same P')
	raw_input('end')
	'''
	
	if binning: 
                if show:
                    plotSplit(xGrid, yGrid, valMC, valData, valDiff, bins, binN)

		xGrid = [ item[0] for item in xGrid ]
		xMC, valDiff = sortLists(xGrid, valDiff)
		xMC, valDiffError = sortLists(xGrid, valDiffError)
                xMC, valMC = sortLists(xGrid, valMC)
                xMC, valMCError = sortLists(xGrid, valMCError)
                xMC, valData = sortLists(xGrid, valData)
                xMC, valDataError = sortLists(xGrid, valDataError)
                # print np.array(valDataError)[:,0]

		print 'Returning...'
                if errorOut:
                    return xMC, yGrid, valMC, valMCError, valData, valDataError, valDiff, valDiffError
                else:
                    return xMC, yGrid, valMC, valData, valDiff, valDiffError

	else:
		valMCNorm = normalizeValList(valMC, hSumMC)
		valDataNorm = normalizeValList(valData, hSumData)
		
		valDiffNorm = []
		for i in range(len(valMCNorm)):
			entryMC = np.array( valMCNorm[i] )
			entryData = np.array( valDataNorm[i] )

			valDiffNorm.append( (entryData - entryMC)/entryData )
			# valDiffNorm.append( abs(entryMC - entryData)/(entryMC + entryData) )

		valMCNorm = np.nan_to_num(valMCNorm)
		valDataNorm = np.nan_to_num(valDataNorm)
		valDiffNorm = np.nan_to_num(valDiffNorm)

                if show:
                    plotSplit(xGrid, yGrid, valMCNorm, valDataNorm, valDiffNorm, bins, binN)

		xGrid = [ item[0] for item in xGrid ]
		xMC, valDiffNorm = sortLists(xGrid, valDiffNorm)
		xMC, valDiffError = sortLists(xGrid, valDiffError)
		return xMC, yGrid, valMCNorm, valDataNorm, valDiffNorm, valDiffError

'''
def plotProjection(x, y, vMC, vData, vDiff, bins, binN):
	vMC = np.array(vDiff)[:,binN]
	lMC = ROOT.TGraph2D('MC', 'MC in bin %d' % binN, (bins-2)**2, np.array(x), np.array(y), np.array(vMC))

	c = ROOT.TCanvas()
	hMC = lMC.Project('x')
	hMC.SetMarkerStyle(3)
	hMC.Draw('P')
	raw_input('end')
'''

def plotSplit(x, y, vMC, vData, vDiff, bins, binN):
	# Select the standoff distance bin number
	dataData = np.array(vData)[:,binN]
	dataMC = np.array(vMC)[:,binN]
	dataDiff = np.array(vDiff)[:,binN]

	# Create Graphs
	if y:
		lMC = ROOT.TGraph2D('MC', 'MC in bin %d' % binN, (bins-2)**2, np.array(x), np.array(y), np.array(dataMC))
		lData = ROOT.TGraph2D('Data', 'Data in bin %d' % binN, (bins-2)**2, np.array(x), np.array(y), np.array(dataData))
		lDiff = ROOT.TGraph2D('Diff', 'Difference in bin %d' % binN, (bins-2)**2, np.array(x), np.array(y), np.array(dataDiff))
	else:
		x = [item[0] for item in x]
		xMC, dataMC = sortLists(x, dataMC)

		lMC = ROOT.TGraph(bins-2, np.array(xMC), np.array(dataMC))
		lMC.SetName('MC')
		lMC.SetTitle('MC in bin %d' % binN)

		xData, dataData = sortLists(x, dataData)
		lData = ROOT.TGraph(bins-2, np.array(xData), np.array(dataData))
		lData.SetName('Data')
		lData.SetTitle('Data in bin %d' % binN)

		xDiff, dataDiff = sortLists(x, dataDiff)
		lDiff = ROOT.TGraph(bins-2, np.array(xDiff), np.array(dataDiff))
		lDiff.SetName('Diff')
		lDiff.SetTitle('Difference in bin %d' % binN,)

	'''
	c = ROOT.TCanvas()
	lMC.Draw('AC*')
	raw_input('end')
	del c
	'''

	# Divide canvas in two pads in the upper
	# and one pad in the lower half
	# c1_1.1 | c1_1.2
	# ---------------
	#      c1.2
	c1 = ROOT.TCanvas('c1', 'Standoff Distance', 1)
	c1.Divide(1, 2)
	c1_1 = c1.cd(1)
	ROOT.gPad.Divide(2, 1)

	# Upper left
	c1_1.cd(1)
	if y:
		lMC.Draw('colz')
		lMC.GetYaxis().SetLimits(min(y)-1, max(y)+1)
	else:
		lMC.Draw('AC*')
	lMC.GetXaxis().SetLimits(min(x)-1, max(x)+1)
	c1.Update()

	# Upper right
	c1_1.cd(2)
	if y:
		lData.Draw('colz')
		lData.GetYaxis().SetLimits(min(y)-1, max(y)+1)
	else:
		lData.Draw('AC*')
	lData.GetXaxis().SetLimits(min(x)-1, max(x)+1)
	c1.Update()

	# Bottom
	c1.cd(2)
	if y:
		lDiff.Draw('colz')
		lDiff.GetYaxis().SetLimits(min(y)-1, max(y)+1)
	else:
		lDiff.Draw('AC*')
	lDiff.GetXaxis().SetLimits(min(x)-1, max(x)+1)
	c1.Update()
	
	raw_input('end')
	del c1
	del c1_1

def getSortData(tree, name, bins, rand=False):
	z = np.linspace(-200, 200, bins - 1)
	theta = np.linspace(-1.05, 1.05, 60)

	h = {}
	hZ = {}
	hTheta = {}

	meanBinsZ = []

	# === CREATE HISTOGRAMS AND MEAN BINS ===
	# Create z and (z, theta) histograms
	for i in range(len(z) - 1):
		zMin, zMax = z[i:i+2]
		zMean = zMin + float(zMax - zMin)/2
		meanBinsZ.append( zMean )
		
		hZ[zMean] = ROOT.TH1D('%s: z = %f' % (name, zMean), '%s: z = %f' % (name, zMean), 20, 0, 200)

		meanBinsTheta = []
		for j in range(len(theta) - 1):
			thetaMin, thetaMax = theta[j:j+2]
			thetaMean = thetaMin + float(thetaMax - thetaMin)/2
			meanBinsTheta.append( thetaMean )
	
			h[(zMean, thetaMean)] = ROOT.TH1D('%s: %f, %f' % (name, zMean, thetaMean), '%s: %f, %f' % (name, zMean, thetaMean), 20, 0, 200)

	# Create theta histograms
	for thetaMean in meanBinsTheta:
		hTheta[thetaMean] = ROOT.TH1D('%s: theta = %f' % (name, thetaMean), '%s: theta = %f' % (name, thetaMean), 20, 0, 200)

	# === FILL HISTOGRAMS ===
	if not rand:
		for i in range(tree.GetEntries()):
			tree.GetEntry(i)
			es = tree.EventSummary

			x, y, zVal = ps.getClusterPos(es)
			'''
			x = es.cluster_x[0]
			y = es.cluster_y[0]
			zVal = es.cluster_z[0]
			'''

			if np.isnan(x) or np.isnan(y) or np.isnan(zVal):
				continue

			if abs(zVal) >= 200.:
				continue

			thetaVal = np.arctan2(y, x) / np.pi
		
			#print (zVal, thetaVal), (np.digitize(zVal, z), np.digitize(thetaVal, theta))
			binTheta = meanBinsTheta[np.digitize(thetaVal, theta) - 1]
                        try:
			    binZ = meanBinsZ[np.digitize(zVal, z) - 1]
                        except:
                            print 'zVal:', zVal
                            continue
                            
                        #print (binZ, binTheta)
			#print z, theta

			d = getStandoff(x, y, zVal)
			hZ[binZ].Fill(d)

			# ATTENTION: Temporary cut!
                        # Used for theta only
			if zVal > (-35 - 15) and zVal < (-35 + 15): 
				hTheta[binTheta].Fill(d)

			h[(binZ, binTheta)].Fill(d)

		'''
		for key in h:
			try:
				h[key].Scale(1./h[key].Integral())
			except:
				pass
		'''

	# RANDOM DATA
	else:
		# Read list from file
		randFile = tree
		posList = cPickle.load(open(randFile, 'rb') )
                # FV = [162., 5., 182]
                # posList = gr.isFiducialList(posList, *list( np.array(FV)/1000. ))
		for item in posList:
			x, y, zVal = np.array( item ) * 1000
			thetaVal = np.arctan2(y, x)/np.pi

			try:
				binTheta = meanBinsTheta[np.digitize(thetaVal, theta) - 1]
			except:
				pass

			try:
				binZ = meanBinsZ[np.digitize(zVal, z) - 1]
			except:
				print 'Failed on z = %f' % zVal
				continue

			d = getStandoff(x, y, zVal)
			hZ[binZ].Fill(d)

			if zVal > (-35 - 15) and zVal < (-35 + 15):
				hTheta[binTheta].Fill(d)

			h[(binZ, binTheta)].Fill(d)

	return h, hZ, hTheta

# === PLOT ===
def getDictsReal(nBins, name = 'Standoff', outdir='./'):
        nBinsTheta = 60

	hMCDict, hDataDict = sortDataFile(outdir + 'mc%s.root' % name, outdir + 'data%s.root' % name, nBins)
	hMCDictZ, hDataDictZ = sortDataFile(outdir + 'mc%sZ.root' % name, outdir + 'data%sZ.root' % name, nBins)

	# ATTENTION: Test
	hMCDictTheta, hDataDictTheta = sortDataFile(outdir + 'mc%sTheta_z.root' % name, outdir + 'data%sTheta_z.root' % name, nBinsTheta)

	return hMCDict, hDataDict, hMCDictZ, hDataDictZ, hMCDictTheta, hDataDictTheta
	
def getDictsRand(nBins, name = 'Standoff'):
	# Load dictionaries containing the standoff
	# distance histograms from root files
	hMCDict, hDataDict = sortDataFile('rand%s.root' % name, 'randDrift%s.root' % name, nBins)
	hMCDictZ, hDataDictZ = sortDataFile('rand%sZ.root' % name, 'randDrift%sZ.root' % name, nBins)

	# ATTENTION: Test
        nBinsTheta = 60
	hMCDictTheta, hDataDictTheta = sortDataFile('rand%sTheta.root' % name, 'randDrift%sTheta.root' % name, nBinsTheta)

	return hMCDict, hDataDict, hMCDictZ, hDataDictZ, hMCDictTheta, hDataDictTheta

def getPlot(random=True, nBins=120, showBin=0, name='Standoff', outdir='./'):
	if random:
		hMCDict, hDataDict, hMCDictZ, hDataDictZ, hMCDictTheta, hDataDictTheta = getDictsRand(nBins, name)
	else:
		print name
		hMCDict, hDataDict, hMCDictZ, hDataDictZ, hMCDictTheta, hDataDictTheta = getDictsReal(nBins, name, outdir=outdir)

	# Plot the information of the standoff histograms in 
	# bin showBin for the specified dictionaries

	# == 2D ==
	# plotStandoff(hMCDict, hDataDict, 15, showBin)

	# == Z ==
	xGrid, yGrid, valMCNorm, valDataNorm, valDiffNorm, valDiffNormError = plotStandoff(hMCDictZ, hDataDictZ, nBins, showBin)

	if not random:
		# Store data in file
		out = {}
		out['x'] = xGrid
		out['mc'] = valMCNorm
		out['data'] = valDataNorm
		out['diff'] = valDiffNorm
		out['diff_err'] = valDiffNormError

		cPickle.dump( out, open('%s.p' % name, 'wb') )

	print 'Creating valDiffShow...'
	valDiffShow = np.array(valDiffNorm)[:,showBin]

	print 'Creating valDiffErrorShow...'
	print len(valDiffNorm), len(valDiffNormError), len(valDiffNorm[0]), len(valDiffNormError[0])
	valDiffErrorShow = np.array(valDiffNormError)[:,showBin]

	print 'Done...'

	'''
	for i in range( len(valDiffShow) ):
		print valDiffShow[i], valDiffErrorShow[i]
	'''

	# = Z-Standoff Plot =
	peak.findPeaks(xGrid, valDiffShow, valDiffErrorShow, True, '(data - mc)/data for Bin #%s' % showBin)

	# == THETA ==
	# xGrid, yGrid, valMCNorm, valDataNorm, valDiffNorm, valDiffError = plotStandoff(hMCDictTheta, hDataDictTheta, nBins, showBin)
        # peak.findPeaks(xGrid, np.array(valDiffNorm)[:,showBin], np.array(valDiffError)[:,showBin], False, '(data - mc)/data for Bin #%s' % showBin)

        nBinsTheta = 60
	hMCDictTheta, hDataDictTheta = sortDataFile('mcStandoffTheta.root', 'dataStandoffTheta.root', nBinsTheta)
	xGrid, yGrid, valMCNorm, valDataNorm, valDiffNorm, valDiffError = plotStandoff(hMCDictTheta, hDataDictTheta, nBinsTheta, showBin)
	return

# === SUPPORT ===
def sortLists(x, y):
	from operator import itemgetter
	return [list(x) for x in zip(*sorted(zip(x, y), key=itemgetter(0)))]

def normalizeValList(valList, sumHist):
	valListNorm = []
	for i in range(len(valList)):
		temp = []
		for j in range(len(valList[i])):
			# temp.append( valList[i][j] / sumHist.GetBinContent(j+1) )
			temp.append( valList[i][j] / sumHist.Integral(5, 20) )
		
		valListNorm.append( temp )
		
	return valListNorm

def saveToFile(hDict, fname):
	f = ROOT.TFile(fname, 'RECREATE')

	for key in hDict:
		hDict[key].Write()

	f.Close()

# Units in [mm], for [m] use same function in
# plot_support!
def getStandoff(x, y, z):
	r = np.sqrt( x**2 + y**2 )

	sd_r = REFLECTORINNERRAD - r
	sd_z = CATHODE_ANODE_y_DISTANCE - abs(z)

	if sd_r > sd_z:
		return sd_z
	else:
		return sd_r

# Save data of the format
# x, y, z0, z1, z2, ...
# to file
def saveToDat(fName, x, y, z):
	f = open(fName, 'w')
	for i in range(len(x)):
		f.write( '%f\t%f\t' % (x[i], y[i]) ) 
		for j in range(len(z[i])):
			f.write( '%f\t' % z[i][j] )
		f.write('\n')
	f.close()

