import numpy as np
import sys
import cPickle
import os

import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')
ROOT.PyConfig.IgnoreCommandLineOptions = True

# ==== CONSTANTS ====
CATHODE_ANODE_y_DISTANCE = 192.23657 # mm
REFLECTORINNERRAD = 183.2356 # mm

def getHistoDiff(h1, h2, name):
	h = h1.Clone()
	h.Add(h2, -1.0)
	h.Divide(h, h2)

	h.SetName(name)
	return h

def getHistoDiffInv(h1, h2, name):
	h = h1.Clone()
	h.Add(h2, -1.0)
	h.Divide(h, h1)

	h.SetName(name)
	return h

def getHistoRelativeDiff3d(h1, h2, name):
	h = h1.Clone()
	xBins = h.GetXaxis().GetNbins()
	yBins = h.GetYaxis().GetNbins()
	zBins = h.GetZaxis().GetNbins()

	for x in range(xBins + 1):
		for y in range(yBins + 1):
			for z in range(zBins + 1):
				b1 = h1.GetBinContent(x, y, z)
				b2 = h2.GetBinContent(x, y, z)
				try:
					# newVal = float(b1 - b2)/abs(b2)
					newVal = float(abs(b1 - b2))/((b1+b2)/2.) #max([b1, b2])
				except:
					newVal 	= 0

				h.SetBinContent(x, y, z, newVal)
	
	h.SetName(name)
	return h

def getHistoRelativeDiff2d(h1, h2):
	h = h1.Clone()
	xBins = h.GetXaxis().GetNbins()
	yBins = h.GetYaxis().GetNbins()

	for x in range(xBins):
		for y in range(yBins):
			b1 = h1.GetBinContent(x, y)
			b2 = h2.GetBinContent(x, y)
			try:
				# newVal = float(b1 - b2)/abs(b2)
				newVal = float(abs(b1 - b2))/((b1+b2)/2.)
			except:
				newVal 	= 0
			if(newVal > 1):
				print newVal

			h.SetBinContent(x, y, newVal)
	
	return h

def getHistoRelativeDiff(h1, h2):
	h = h1.Clone()
	xBins = h.GetXaxis().GetNbins()

	for x in range(xBins):
		b1 = h1.GetBinContent(x)
		b2 = h2.GetBinContent(x)
		try:
			# newVal = float(b1 - b2)/abs(b2)
			newVal = float(abs(b1 - b2))/((b1+b2)/2.)
		except:
			newVal 	= 0
		h.SetBinContent(x, newVal)
	
	return h

def fillHisto(tree, hist, cut, select, name):
	h = hist.Clone('h')
	cval = str(select) + '>>h'
	tree.Draw(cval, cut, 'goff')
	h.SetName(name)
	h.Sumw2()
	# h.Scale(1./h.Integral())
	return h

def fillHisto3d(tree, hist, cut, name, cyl=True, rsquare=False):
	h = hist.Clone('h')
	ROOT.gROOT.cd()
	#t = tree 
        t = tree.CopyTree(cut)
	nEntries = t.GetEntries()

	print 'Filling histogram %s with %d entries...' % (name, nEntries)
	for i in range(nEntries):
		t.GetEntry(i)
		es = t.EventSummary

		posX = np.array(es.cluster_x)
		posY = np.array(es.cluster_y)
		posZ = np.array(es.cluster_z)

		if(cyl):
			maxR, maxZ, maxTheta = getMaxCyl(posX, posY, posZ, es)
			if(not maxR or not maxZ):
				continue
			if(rsquare):
				h.Fill(maxR**2, maxZ, maxTheta)
			else:
				h.Fill(maxR, maxZ, maxTheta)
		else:
			try:
				maxX, maxY, maxZ = getMaxCart(posX, posY, posZ)
			except:
				continue

			'''
			# ATTENTION: Only the cluster of getMaxCart()
			# contributes to the standoff distance!
			for j in range(len(posX)):
				maxX = posX[j]
				maxY = posY[j]
				maxZ = posZ[j]
			'''

			if(not maxX and not maxY and not maxZ):
				continue

			h.Fill(maxX, maxY, maxZ)

		if not (i % int(nEntries/100.)):
			perc = int(100.*(i+1)/nEntries + 1)
			statusBar(perc, width='screen')
	print '\n'

	# h.Scale(1./h.GetMaximum())
	h.Scale(1./h.Integral())
	# h.Scale(1./nEntries)
	h.SetName(name)
	
	return h

def getMaxCyl(x, y, z, es):
	maxTheta = GetMaxTheta(x, y)
	maxZ = GetMaxZ(z)
	maxR = (es.GetMaxR() * (es.GetMaxR() < 900))
	return maxR, maxZ, maxTheta

def getMaxCart(x, y, z):
        x, y, z = np.array(x), np.array(y), np.array(z)

	# Return value with maximum radius
	r = np.array([ri for ri in np.sqrt((x**2 + y**2)) if abs(ri) < 900])

	if ( REFLECTORINNERRAD - max(r) ) > ( CATHODE_ANODE_y_DISTANCE - abs(GetMaxZ(z)) ):
		i = np.argmax(z)
	else:
		i = np.argmax(r)

	xMax = x[i]
	yMax = y[i]
	zMax = z[i]
	# xMax = np.amax(x * (np.absolute(x) < 900))
	# yMax = np.amax(y * (np.absolute(y) < 900)) 
	return xMax, yMax, zMax

def GetMaxZ(posZ):
        posZ = np.array([pZ for pZ in posZ if abs(pZ) < 900])
	zMax = np.amax(np.absolute(posZ))
	# zMin = np.amin(np.absolute(posZ))
	#if abs(zMax) < abs(zMin):
	#	return zMin
	#else:
	#	return zMax
        return zMax

# return theta in units of pi
def GetMaxTheta(posX, posY):
	xMax = np.amax(posX * (np.absolute(posX) < 900))
	yMax = np.amax(posY * (np.absolute(posY) < 900))
	return np.arctan2(xMax, yMax)/np.pi

def getCut(calibCut=True, energyCut=True, type='ms', MC=True, eMin=750, eMax=3500, sideCut=False, phase2=False):
        phase2 = True

	eventSum = ROOT.EXOEventSummary()
	prepTree = ROOT.EXOPreprocessedTreeModule() 
	prepTree.SetIsData(not MC)
	prepTree.SetEventSummary(eventSum)

	fFOVs = prepTree.GetFOVs()

	fFOVs.SetBooleanFlag('fIsDataOrMC', not MC)
        # Phase 1
        if not phase2:
            fFOVs.SetStringOption('fDiagonalCutDBFlavor', '2013-0nu-denoised')
        # Phase 2
        # else:
        #    fFOVs.SetStringOption('fDiagonalCutDBFlavor', 'phase2_v1')
	fFOVs.SetBooleanFlag('fCutMissingPosition', True)
	fFOVs.SetBooleanFlag('fVetoLikeCut', calibCut)
	fFOVs.SetBooleanFlag('fEqualizeThreshold', False)
	fFOVs.SetBooleanFlag('fApplyVetoes', True)
	fFOVs.SetBooleanFlag('flRandomTrigger', False)
	fFOVs.SetBooleanFlag('fApplyBadTimesVeto', True)
	fFOVs.SetBooleanFlag('fUse137XeVeto', False)

	# ATTENTION: It doesn't work to set the 
	# fiducial volume with flags
	# Use ROOT.EXOFiducialVolume.SetUserHexCut() 
	# in saveHistograms.py instead

        if sideCut:
            sCut = ' && (abs(atan2(cluster_y, cluster_x)) >= 0.52359877559829882) && (abs(atan2(cluster_y, cluster_x)) <= 2.6179938779914944) '
        else:
            sCut = ''

	if energyCut:
		if MC:
			enCut = " && (energy_mc > %d. && energy_mc < %d.) " % (eMin, eMax)
		else:
			enCut = " && (energy > %d. && energy < %d.) " % (eMin, eMax)
	else:
		enCut = ''

	if type == 'ms':
		typeCut = "&& multiplicity > 1.1 "
        elif type == 'ss':
		typeCut = "&& abs(multiplicity - 1) < 0.1"
        else:
                typeCut = ''

        cut = prepTree.GetDefaultCut()
        cut = cut.replace('&& !isSolicitedTrigger', '')
        if phase2:
            cut = cut.replace('!EventSummary.isDiagonallyCut() && ', '')
	return cut + enCut + typeCut + sCut

# Do no multiplicity cut!
def getCutCombined(calibCut=True, energyCut=True, MC=True, eMin=750, eMax=3500):
	eventSum = ROOT.EXOEventSummary()
	prepTree = ROOT.EXOPreprocessedTreeModule() 
	prepTree.SetIsData(not MC)
	prepTree.SetEventSummary(eventSum)

	fFOVs = prepTree.GetFOVs()

	fFOVs.SetBooleanFlag('fIsDataOrMC', not MC)
        '''
	fFOVs.SetStringOption('fDiagonalCutDBFlavor', '2013-0nu-denoised')
	fFOVs.SetBooleanFlag('fCutMissingPosition', True)
	fFOVs.SetBooleanFlag('fVetoLikeCut', calibCut)
	fFOVs.SetBooleanFlag('fEqualizeThreshold', False)
	fFOVs.SetBooleanFlag('fApplyVetoes', True)
	fFOVs.SetBooleanFlag('flRandomTrigger', False)
	fFOVs.SetBooleanFlag('fApplyBadTimesVeto', True)
	fFOVs.SetBooleanFlag('fUse137XeVeto', False)
        '''

	if energyCut:
		if MC:
			enCut = " && (energy_mc > %d. && energy_mc < %d.) " % (eMin, eMax)
		else:
			# enCut = " && (energy_ss > %d. && energy_ss < %d.) && (energy_ms > %d. && energy_ms < %d.) " % (eMin, eMax, eMin, eMax)
                        enCut = " && (energy > %d. && energy < %d.) " % (eMin, eMax)
	else:
		enCut = ''

	return prepTree.GetDefaultCut() + enCut

def getSlice(hist, BIN_SET, option='x', value=0):
	if(option is 'x'):
		BINM = BIN_SET[0:3]
		BINA = BIN_SET[3:6]
		BINB = BIN_SET[6:9]
	elif(option is 'y'):
		BINM = BIN_SET[3:6]
		BINA = BIN_SET[0:3]
		BINB = BIN_SET[6:9]
	elif(option is 'z'):
		BINM = BIN_SET[6:9]
		BINA = BIN_SET[0:3]
		BINB = BIN_SET[3:6]
	else:
		return False

	h = ROOT.TH2D('h', 'h', BINA[0], BINA[1], BINA[2], BINB[0], BINB[1], BINB[2])

	BIN_WIDTH_M = getBinWidth(BINM) 
	# BIN_WIDTH_M = (BINM[2] - BINM[1])/BINM[0]
	valBin = int((value - BINM[1])/BIN_WIDTH_M)

	for a in range(BINA[0]):
		for b in range(BINB[0]):
			# Function has to be changed according to option set
			if(option is 'x'):
				val = hist.GetBinContent(valBin, b, a)
			elif(option is 'y'):
				val = hist.GetBinContent(a, valBin, b)
			else:
				val = hist.GetBinContent(a, b, valBin)

			h.SetBinContent(a, b, val)

	h.SetEntries(BINA[0] * BINB[0])
	# h.Scale(1./h.Integral())
	return h

def getSlice2d(hist, option, value):
	h = hist.Clone()

	if(option == 'x'):
		binY = h.GetYaxis().FindBin(value)
		hOut = h.ProjectionY('', binY, binY).Clone()
	elif(option == 'y'):
		binX = h.GetXaxis().FindBin(value)
		hOut = h.ProjectionX('', binX, binX).Clone()
	else:
		return False

	return hOut

def getMaximum(hist, BIN_SET):
	BINS_X, BIN_X_MIN = BIN_SET[0:2]
	BINS_Y, BIN_Y_MIN = BIN_SET[3:5]
	BIN_W_X = getBinWidth(BIN_SET[0:3])
	BIN_W_Y = getBinWidth(BIN_SET[3:6])

	# bin = z * BINS_X * BINS_Y + y * BINS_X + x
	binM = hist.GetMaximumBin()
	z = int(binM / ((BINS_X+2) * (BINS_Y+2)))
	y = int((binM - z * ((BINS_X+2) * (BINS_Y+2))) / (BINS_X+2))
	x = binM - (z * (BINS_X+2) * (BINS_Y+2)) - (y * (BINS_X+2))

	xPos = hist.GetXaxis().GetBinCenter(x)
	yPos = hist.GetYaxis().GetBinCenter(y)
	zPos = hist.GetZaxis().GetBinCenter(z)

	return xPos, yPos, zPos

def getMaximum2d(hist, BIN_SET):
	BINS_X, BIN_X_MIN = BIN_SET[0:2]
	BIN_W_X = getBinWidth(BIN_SET[0:3])

	# bin = y * BINS_X + x
	binM = hist.GetMaximumBin()
	y = int(binM / (BINS_X + 2))
	x = binM - y * (BINS_X + 2)

	return hist.GetXaxis().GetBinCenter(x), hist.GetYaxis().GetBinCenter(y)

def getBinWidth(BIN_SET):
	BINS, BIN_MIN, BIN_MAX = BIN_SET
	return float(BIN_MAX - BIN_MIN)/BINS

def setShow(show):
	if(show):
		ROOT.gROOT.SetBatch(ROOT.kFALSE)
	else:
		ROOT.gROOT.SetBatch(ROOT.kTRUE)

def getChain(treeName, fileList):
	ch = ROOT.TChain(treeName)
	if(not ch): return False

	for item in fileList:
		ch.Add(item)

	return ch

def getFileList(dirName, fName):
	fList = os.listdir(dirName)
	fListFinal = []

	for item in fList:
		if fName in item:
			fListFinal.append( '%s/%s' % (dirName, item) )

	return fListFinal

def mergeLists(fList, outName, verbose=False):
	listTotal = []

	for f in fList:
		if verbose:
			print 'Merging %s...' % f
		listTotal += cPickle.load(open(f, 'rb') )

	if verbose:
		print 'Saving file %s' % outName
		print
	cPickle.dump(listTotal, open(outName, 'wb'))

def saveHistos(fName, hDict):
	f = ROOT.TFile(fName, 'recreate')
	for h in hDict:
		print 'Write %s to file...' % (h)
		hDict[h].Write()
	f.Close()
	print 'Saved histograms to file %s\n' % fName

def openAllHistos(f):
	if(not f):
		print 'File not found!'
		return False

	hNameList = f.GetListOfKeys()
	hDict = {}
	for h in hNameList:
		hName = h.GetName()
		print 'Read histogram %s' % (hName)
		hDict[hName] = f.Get(hName)

	return hDict

def openSingleHisto(f, hName):
	if(not f):
		print 'File %s not found!' % (fName)

	return f.Get(hName)

def makeDir(path):
	try:
		os.mkdir(path)
	except OSError:
		pass

def getHist(f):
	hList = openAllHistos(f)

	for key, h in hList.iteritems():
		print key, h
		if('data' in key):
			dataHist = hList[key]
		if('mc' in key):
			mcHist = hList[key]
		if('diff' in key and 'Rel' in key):
			diffRelHist = hList[key]
		elif('diff' in key):
			diffHist = hList[key]

	return mcHist, dataHist, diffHist, diffRelHist

def statusBar(perc, width='screen'):
	if width == 'screen':
		width = int(os.popen('tput cols', 'r').read()) - 8

	# width = int(width) - 8
	p = int(perc * float(width)/100)
	sys.stdout.write('\r')
	sys.stdout.write('[%-*s] %d%%' % (int(width), '='*p, perc))
	sys.stdout.flush()	

# isFiducial
# ==========
# For SS events only!
# Source: http://math.stackexchange.com/questions/41940/is-there-an-equation-to-describe-regular-polygons
def isFiducial(x, y, z, Apo, zMin, zMax):
		r = np.sqrt( x**2 + y**2 )
		theta = np.arctan2(y, x)
		#if theta < 0:
		#	theta = 2*np.pi - abs(theta)

		n = 6 # Hexagon
		R0 = (2./3) * np.sqrt(3) * Apo
		R = R0 * np.cos(np.pi/n) / ( np.cos((theta - np.pi/2) % (2*np.pi/n) - np.pi/n) )

		if (r <= R) and (abs(z) >= zMin) and (abs(z) <= zMax):
				return True
		else:
				return False

def getApothem(x, y, z):
    r = np.sqrt( x**2 + y**2 )
    theta = np.arctan2(y, x)
    n = 6

    A = np.sqrt(3)/2 * r * (1./np.cos(np.pi/n)) * np.cos((theta - np.pi/2) % (2*np.pi/n) - np.pi/n)
    return A

# ==== STANDOFF STUFF ====
# Units are [mm]!
def getClusterPos(es, cut=True):
	posX = np.array(es.cluster_x)
	posY = np.array(es.cluster_y)
	posZ = np.array(es.cluster_z)

        # if posX.any():
        #    print posX, posY, posZ

        if not cut:
            return posX, posY, posZ

        if not posX.all() or not posY.all() or not posZ.all() or abs(posX).any() > 900 or abs(posY).any() > 900:
            return np.nan, np.nan, np.nan 

        try:
	    maxX, maxY, maxZ = getMaxCart(posX, posY, posZ)
        except:
            return np.nan, np.nan, np.nan
	return maxX, maxY, maxZ

# Units have to be [m]!
def getStandoff(x, y, z):
	r = np.sqrt( x**2 + y**2 )

	sd_r = REFLECTORINNERRAD/1000. - r
	sd_z = CATHODE_ANODE_y_DISTANCE/1000. - abs(z)

	if sd_r > sd_z:
		return sd_z
	else:
		return sd_r

def getStandoffHisto(fName, tree, FV, type='ss', MC=True, name='default', out=None, bins=20):
	print 'Processing %s...' % name
	ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
	ROOT.gROOT.cd()

	select = 'standoff_distance'
	cut = getCut(calibCut=True, energyCut=True, type=type, MC=MC, eMin=980, eMax=9800)
        if MC:
            type = 'mc'
        # cut += ' && ( (energy_%s > 1000 && energy_%s < 1361) || (energy_%s > 1561 && energy_%s < 2000 ) )' % (type, type, type, type)
        # cut += ' && cluster_z > 0'

        print cut
	ROOT.EXOFiducialVolume.SetUserHexCut(*FV)

	f = ROOT.TFile.Open(fName, 'UPDATE')
	hist = ROOT.TH1D(name, name, bins, 0, 200)
        hStandoff = fillHisto(tree, hist, cut, select, name)
        if out:
            fOut = open(out, 'w')
            for i in range(1, int(hStandoff.GetSize())):
                data = hStandoff.GetBinContent(i)
                fOut.write('%d\n' % data)
            fOut.write('\n')
            fOut.close()

        hStandoff.Write('', ROOT.TObject.kOverwrite)
	f.Close()
        
def getStandoffDistMan(tree, FV, type='ss', MC=True, bins=20):
        RIR = 183.2356
        FV = [171, 0, 182]

	ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
	ROOT.gROOT.cd()
    
	cut = getCut(calibCut=True, energyCut=True, type=type, MC=MC, eMin=980, eMax=9800)
        print cut
	treeCut = tree.CopyTree(cut, 'fast') 	

        N = treeCut.GetEntries()
        
	posList = []
	for i in range( N ):
		treeCut.GetEntry(i)
		es = treeCut.EventSummary
                mul = es.multiplicity

		x, y, z = getClusterPos( es, True )
                if np.isnan(x) or np.isnan(y) or np.isnan(z):
                    continue

                x_, y_, z_ = x, y, z
                posList.append( list( np.array([x_, y_, z_]) )) # / 1000. ) )

        posList = np.array( posList )

        plotHex(posList[:,0], posList[:,1], FV[0], RIR, bins, r'$x$ [mm]', r'$y$ [mm]')
        plotHex(posList[:,0], posList[:,2], FV[0], FV[2], bins, r'$x$ [mm]', r'$z$ [mm]')
        plotHex(posList[:,1], posList[:,2], RIR, FV[2], bins, r'$y$ [mm]', r'$z$ [mm]')
        raw_input('')

def plotHex(x, y, rangeX, rangeY, bins, xlabel, ylabel):
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap
        from matplotlib import cm
        from copy import copy
        from scipy.stats import kendalltau
        import seaborn as sns
        # sns.set(style = 'ticks')
        # sns.set_palette()

        # cmapV = copy(cm.get_cmap('viridis'))
        # cmapV.set_under('w', alpha=0)
        # cmapV.set_bad('w', alpha=0)

        colorsData = sns.color_palette('viridis', 20) # sns.cubehelix_palette(20, start=.5, rot=-.75, dark=0.2, light=.7)
        colors = [(1., 1., 1.)] + colorsData # ['#FFFFFFFF'] + colorsData
        print colorsData
        cmap = ListedColormap(sns.color_palette(colors)) #.as_hex())
        print cmap

        colorsLine = [[c*1.2 if c*1.2 <= 1 else 1. for c in color] for color in colorsData]
        cmapLines = ListedColormap(sns.color_palette(colorsLine)) # sns.cubehelix_palette(20, start=.5, rot=-.75, dark=0.15, light=.65, as_cmap=True)

        #figXY, axXY = plt.subplots()
        # hb = axXY.hexbin(x, y, gridsize=bins, vmin=0.1, cmap=cmap)
        lm = sns.jointplot(x, y, kind='hex', stat_func=None, space=0, cmap=cmap) # , color='b')
        lm.plot_joint(sns.kdeplot, n_levels=6, cmap=cmapLines)
        lm.set_axis_labels(xlabel, ylabel)
        plt.xlim(-rangeX, rangeX)
        plt.ylim(-rangeY, rangeY)
        plt.subplots_adjust(left=0.15, bottom=0.1)

        # plt.subplots_adjust(left=0.2, right=0.8, top=0.8, bottom=0.2)
        # cbarAx = lm.fig.add_axes([.85, .25, .05, .4])
        # cb = plt.colorbar(cax=cbarAx)
        # cb.set_label('Counts')

        # lm.fig.xlabel(xlabel)
        # lm.fig.ylabel(ylabel)

        # axes = lm.axes
        # print axes
        # cb = figXY.colorbar(hb, ax=axXY)

        #axXY.set_xlim(-rangeX, rangeX)
        #axXY.set_ylim(-rangeY, rangeY)
        #axXY.set_xlabel(xlabel)
        #axXY.set_ylabel(ylabel)

        # figXY.show()
        plt.show()

def getStandoffHistoMan(fName, tree, FV, type='ss', MC=True, name='default', out=None, bins=20):
	import generate_random as gr
	ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
	ROOT.gROOT.cd()
    
	cut = getCut(calibCut=True, energyCut=True, type=type, MC=MC, eMin=980, eMax=9800)
        print cut
	treeCut = tree.CopyTree(cut, 'fast') 	

        N = treeCut.GetEntries()
        # N = 140000
    
	posList = []
        j = 0
        # while N > 0:
	for i in range( N ):
                # if not j % 1000:
                #    print N, j

		treeCut.GetEntry(i)
		# treeCut.GetEntry(j)
                # j += 1
		es = treeCut.EventSummary
                mul = es.multiplicity

                '''
                if type=='ss':
                    if mul > 1.1:
                        continue
                else:
                    if mul < 1.1:
                        continue
                '''

		x, y, z = getClusterPos( es, True )
                if np.isnan(x) or np.isnan(y) or np.isnan(z):
                    continue

                '''
                # Cut on angular segment in volume
                Apothem = 171
                theta = np.arctan2(x, y)
                thetaCutVal = 30.*np.pi/180 - np.arccos(float(Apothem)/REFLECTORINNERRAD)
                thetaOffset = 0 # -5.*np.pi/180
                if not (theta <= (thetaCutVal+thetaOffset) and theta >= (-thetaCutVal+thetaOffset)):
                    continue
                '''

                # if z > 0 or abs(z) > 160 or abs(z) < 40 or np.sqrt( x**2 + y**2 ) > 183.:
                #    continue

                # if x and y and z:
                # if x.any() and abs(x).any > 900:
                    # print x

                # for i in range( len(x) ):
                # x_, y_, z_ = x[i], y[i], z[i]
                x_, y_, z_ = x, y, z
                # if not isFiducial(x_, y_, z_, 162, 5, 182): # or z_ < 0 or abs(x_) > 900 or not (abs(z_) > 20 and abs(z_) < 130):
                #    continue
                posList.append( list( np.array([x_, y_, z_]) )) # / 1000. ) )
                N -= 1
	
	# fidVerts = gr.isFiducialList(posList, *list( np.array(FV)/1000. ))
        # fidVerts = posList

	f = ROOT.TFile.Open(fName, 'UPDATE')

        # Rewrote getStandoffList so it has a maximum of 200 instead of 0.2
	standoffCut = gr.getStandoffList(posList, name, bins, True)
	# standoffCut.Scale(1./standoffCut.Integral())
	standoffCut.Write('', ROOT.TObject.kOverwrite)
	f.Close()

	del treeCut

def getCompHisto(fName, treeData, treeMC, FV, type='ss', selectList=['standoff_distance'], name='default'):
	ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
	ROOT.gROOT.cd()

	cutData = getCut(calibCut=True, energyCut=True, type=type, MC=False, eMin=750, eMax=3500)
        cutMC = getCut(calibCut=True, energyCut=True, type=type, MC=True, eMin=750, eMax=3500)

        f = ROOT.TFile.Open(fName, 'UPDATE')
        for select in selectList:
            n = name + select
            nameData = n + 'Data' + type.upper()
            nameMC = n + 'MC' + type.upper()
            print 'Processing %s...' % nameData

            if select=='standoff_distance':
                binSettings = [20, 0, 200]
            elif select=='energy':
                binSettings = [200, 1000, 3000]
            elif select=='multiplicity':
                if type=='ms':
                    binSettings = [15, 1, 16]
                else:
                    continue
            elif select=='u_mst_metric' or select=='v_mst_metric':
                binSettings = [5, 0, 300] 
            elif select=='num_coll_wires' or select=='num_ind_wires':
                binSettings = [10, 0, 10]

            histData = ROOT.TH1D(nameData, nameData, *binSettings)
            histMC = ROOT.TH1D(nameMC, nameMC, *binSettings)

            hData = fillHisto(treeData, histData, cutData, select, nameData)
            if select=='energy':
                hMC = fillHisto(treeMC, histMC, cutMC, 'energy_mc', nameMC)
            else:
                hMC = fillHisto(treeMC, histMC, cutMC, select, nameMC)

            hData.Scale(1./hData.Integral())
            hMC.Scale(1./hMC.Integral())

            hRes = hData.Clone()
            hRes.SetNameTitle(n+'Res'+type.upper(), n+'Res'+type.upper())
            hRes.Add(hMC, -1)
            hRes.Divide(hMC)

            hData.Write('', ROOT.TObject.kOverwrite)
            hMC.Write('', ROOT.TObject.kOverwrite)
            hRes.Write('', ROOT.TObject.kOverwrite)

	f.Close()

def getZHisto(fName, tree, FV, type='ss', MC=True, name='default'):
    print 'Processing %s...' % name
    ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
    ROOT.gROOT.cd()

    select = 'cluster_z'
    cut = getCut(calibCut=True, energyCut=True, type=type, MC=MC, eMin=980, eMax=9800)
    print cut
    
    f = ROOT.TFile.Open(fName, 'UPDATE')
    hist = ROOT.TH1D(name, name, 30, -200, 200)
    hStandoff = fillHisto(tree, hist, cut, select, name)
    hStandoff.Write('', ROOT.TObject.kOverwrite)
    f.Close()

def getApothemHisto(fName, tree, FV, type='ss', MC=True, name='default'):
	ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
	ROOT.gROOT.cd()
    
	cut = getCut(calibCut=True, energyCut=False, type=type, MC=MC, eMin=980, eMax=9800)
        if MC:
            type = 'mc'
        cut += ' && ( (energy_%s > 1000 && energy_%s < 1361) || (energy_%s > 1561 && energy_%s < 2000 ) )' % (type, type, type, type)
        
	ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
        print cut
	treeCut = tree.CopyTree( cut )

        f = ROOT.TFile.Open(fName, 'UPDATE')
        histZ = ROOT.TH1D('z'+name, 'z'+name, 30, -200, 200)
        histA = ROOT.TH1D('a'+name, 'a'+name, 20, 0, 180)
        histEn = ROOT.TH1D('en'+name, 'en'+name, 300, 0, 2000)
        
        enList = []
        for j in range(treeCut.GetEntries()):
		treeCut.GetEntry(j)
		es = treeCut.EventSummary
		x, y, z = getClusterPos( es, False )
                if type == 'ss':
                    # print es.cluster_z[0], es.energy_ss
                    enList.append( es.energy_ss )
                    histEn.Fill( es.energy_ss )
                if type == 'mc':
                    enList.append( es.energy_mc )
                    histEn.Fill( es.energy_mc )

                histZ.Fill(z[0], 1.)

                A = getApothem(x, y, z)
                histA.Fill(A[0], 1.)

        if type == 'ss' or type == 'mc':
            enList = np.array( enList )
            print 'Mean: %f; Std: %f' % (np.mean(enList), np.std(enList))
            c = ROOT.TCanvas()
            histEn.Draw()
            raw_input('')

        histZ.Write('', ROOT.TObject.kOverwrite)
        histA.Write('', ROOT.TObject.kOverwrite)
        f.Close()

# NEvents only used for MC
def getApothemThetaHisto(fName, tree, FV, type='ss', MC=True, NEvents=0, name='default', N=3):
	ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
	ROOT.gROOT.cd()
    
	cut = getCut(calibCut=True, energyCut=True, type=type, MC=MC, eMin=980, eMax=9800)
        if MC:
            type = 'mc'
        # cut += ' && ( (energy_%s > 800 && energy_%s < 1200) )' % (type, type)
        # cut += ' && ( (energy_%s > 800 && energy_%s < 1400) || (energy_%s > 1530 && energy_%s < 1620 ) )' % (type, type, type, type)
        # cut += ' && (energy_%s > 800 && energy_%s < 2000) ' % (type, type)
        print cut

	ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
	treeCut = tree.CopyTree( cut )
        f = ROOT.TFile.Open(fName, 'UPDATE')

        Apothem = FV[0]
        if Apothem > 171:
            Apothem = 171
            FV[0] = Apothem

        # Maximum radius due to cut
        maxCoord = Apothem * 2./np.sqrt(3)

        # Empty result lists
        histZlist, histAlist, histSOlist, histRlist, histRZlist = [], [], [], [], []

        # 2d histogram bins
        # bin width in z should be a fraction of the field shaping ring distance
        dFSR = 16.87
        binWidthZ = dFSR / 3
        binZ = int( 400 / binWidthZ )
        zRange = 0.5 * binZ * binWidthZ

        binWidthA = 9. / 3
        binA = int( Apothem / binWidthA )

        binWidthSO = 9. / 3
        binSO = int( 200 / binWidthSO )

        binWidthR = 1. # 9. / 3
        binR = int( maxCoord / binWidthR ) # 20
        # binWidthZ = 400./binZ

        for i in range(N+1):
            histZlist.append( ROOT.TH1D('z'+name+str(i), 'z'+name+str(i), binZ, -200, 200) )
            histAlist.append( ROOT.TH1D('a'+name+str(i), 'a'+name+str(i), binA, 0, Apothem) )
            # histSOlist.append( ROOT.TH1D('so'+name+str(i), 'so'+name+str(i), int(200 // 9.4), 0, (200 // 9.4)*9.4) )
            histSOlist.append( ROOT.TH1D('so'+name+str(i), 'so'+name+str(i), binSO, 0, binWidthSO * binSO) )
            histRlist.append( ROOT.TH1D('r'+name+str(i), 'so'+name+str(i), binR, 0, binWidthR * binR) )
            histRZlist.append( ROOT.TH2D('rz'+name+str(i), 'rz'+name+str(i), binR, 0, maxCoord, binZ, -zRange, zRange) )

        l = [(2*n+1)*30*np.pi/180 for n in range(-3,3)] 

        # If data: get all events
        # If MC: use number of events specified, usually returned
        #        from this function applied to data
        # ATTENTION: MC file needs more events than data file for proper result!
        if not MC:
            NEvents = treeCut.GetEntries()
        else:
            nEvents = treeCut.GetEntries()
            if nEvents < NEvents:
                NEvents = nEvents

        for j in range( NEvents ):
		treeCut.GetEntry(j)
		es = treeCut.EventSummary
		x, y, z = getClusterPos( es, True )
                if np.isnan(x) or np.isnan(y) or np.isnan(z):
                    continue

                theta = np.arctan2(y, x)
                A = getApothem(x, y, z)
                r = np.sqrt(x**2 + y**2)

                idx = ( (np.digitize(theta, l) + N) % N )
                so = es.standoff_distance

                histZlist[idx].Fill(z, 1.)
                histAlist[idx].Fill(A, 1.)
                histSOlist[idx].Fill(so, 1.)
                histRlist[idx].Fill(r, 1.)
                histRZlist[idx].Fill(r, z, float( N )/histRZnorm(r, binWidthR, binWidthZ, Apothem))
                
                # Last entries in lists contain all points
                histZlist[-1].Fill(z, 1.)
                histAlist[-1].Fill(A, 1.)
                histSOlist[-1].Fill(so, 1.)
                histRlist[-1].Fill(r, 1.)
                histRZlist[-1].Fill(r, z, float( 1. )/histRZnorm(r, binWidthR, binWidthZ))

        for i in range(N+1):
            histZlist[i].Write('', ROOT.TObject.kOverwrite)
            histAlist[i].Write('', ROOT.TObject.kOverwrite)
            histSOlist[i].Write('', ROOT.TObject.kOverwrite)
            histRlist[i].Write('', ROOT.TObject.kOverwrite)
            histRZlist[i].Write('', ROOT.TObject.kOverwrite)

            # histRZlist[i].Draw('colz')
            # raw_input('')

        f.Close()

        return NEvents

def scatterEnergy(tree, FV, art='ss'):
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib import cm
    from copy import copy
    import seaborn as sns
    ROOT.EXOFiducialVolume.SetUserHexCut(*FV)
    ROOT.gROOT.cd()

    # cut = getCut(calibCut=True, energyCut=False, type=art)
    treeCut = tree # tree.CopyTree( cut )

    # h = ROOT.TH2D('h', 'h', 2000, 500, 2800, 2000, 0, 2800)

    en = []
    for i in range(treeCut.GetEntries()):
        treeCut.GetEntry(i)
        es = treeCut.EventSummary
        scint = es.e_scint * 0.984487
        charge = es.e_charge * 1.034517
        if (scint < 500 or scint > 3500) or (charge < 500 or charge > 3500):
            continue

        # h.Fill(scint, charge, 1)
        en.append( (scint, charge) )

    # c = ROOT.TCanvas()
    # h.Draw('colz')
    # raw_input('')

    colorsData = sns.color_palette('viridis', 256)
    colors = [(1., 1., 1.)] + colorsData
    cmap = ListedColormap(sns.color_palette(colors))

    colorsLine = [[c*1.8 if c*1.8 <= 1 else 1. for c in color] for color in colorsData]
    cmapLines = ListedColormap(sns.color_palette(colorsLine))

    scint, charge = np.array(en)[:,0], np.array(en)[:,1]
    print scint

    charge = np.array( [c for c in charge if c >= 500 and c <= 3500] )
    scint = np.array( [s for s in scint if s >= 500 and s <= 3500] )

    lm = sns.jointplot(charge, scint, kind='hex', stat_func=None, space=0, cmap=cmap, joint_kws={'gridsize': (200, 200)}, marginal_kws=dict(bins=500), xlim=(700, 3500), ylim=(700, 3500))
    # lm.plot_joint(sns.kdeplot, gridsize=200, n_levels=6, cmap=cmapLines, bw=.1)
    lm.set_axis_labels(r'$E_\mathrm{charge}$ [keV]', '$E_\mathrm{scint}$ [keV]')
    # plt.xlim(500, 3500)
    # plt.ylim(500, 3500)
    plt.subplots_adjust(left=0.15, bottom=0.1)
    plt.show()

    '''
    import matplotlib.pyplot as plt
    from scipy.stats import gaussian_kde
    f, ax = plt.subplots()

    x = np.array(en)[:,0]
    y = np.array(en)[:,1]
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    ax.scatter(x, y, c=z, s=100, edgecolor='')
    plt.show()
    raw_input('')
    '''

    return

# === GENERATE2NBB ===
# Generate random distribution of events within the TPC and apply
# fiducial cut. Process the events like in getApothemHisto.
# Fiducial cut volume of (162, 5, 182) is hardcoded.
def generate2nbb(fName, NEvents, FV, name='default', N=3):
    import random 
    f = ROOT.TFile.Open(fName, 'UPDATE')

    Apothem = FV[0]
    if Apothem > 171:
        Apothem = 171
        FV[0] = Apothem

    # Maximum radius due to cut
    maxCoord = Apothem * 2./np.sqrt(3)

    # Empty result lists
    histZlist, histAlist, histSOlist, histRlist, histRZlist = [], [], [], [], []

    # 2d histogram bins
    dFSR = 16.87
    binWidthZ = dFSR / 3
    binZ = int( 400 / binWidthZ )
    zRange = 0.5 * binZ * binWidthZ

    binR = 20
    binWidthR = maxCoord/binR

    for i in range(N+1):
        histZlist.append( ROOT.TH1D('z'+name+str(i), 'z'+name+str(i), 30, -200, 200) )
        histAlist.append( ROOT.TH1D('a'+name+str(i), 'a'+name+str(i), Apothem // 9, 0, Apothem) )
        histSOlist.append( ROOT.TH1D('so'+name+str(i), 'so'+name+str(i), int(200 // 9.4), 0, (200 // 9.4)*9.4) )
        histRlist.append( ROOT.TH1D('r'+name+str(i), 'r'+name+str(i), 20, 0, maxCoord) )
        histRZlist.append( ROOT.TH2D('rz'+name+str(i), 'rz'+name+str(i), binR, 0, maxCoord, binZ, -zRange, zRange) )

    l = [(2*n+1)*30*np.pi/180 for n in range(-3,3)] 

    while NEvents > 0:
        x, y, z = random.uniform(-maxCoord, maxCoord), random.uniform(-maxCoord, maxCoord), random.uniform(-maxCoord, maxCoord)

        if not isFiducial(x, y, z, *FV):
            continue
        else:
            NEvents -= 1

        theta = np.arctan2(y, x)
        A = getApothem(x, y, z)
        r = np.sqrt(x**2 + y**2)

        idx = ( (np.digitize(theta, l) + N) % N )
        so = getStandoff(x/1000., y/1000., z/1000.)*1000 

        histZlist[idx].Fill(z, 1.)
        histAlist[idx].Fill(A, 1.)
        histSOlist[idx].Fill(so, 1.)
        histRlist[idx].Fill(r, 1.)
        histRZlist[idx].Fill(r, z, 1./histRZnorm(r, binWidthR, binWidthZ))

        # Last entries in lists contain all points
        histZlist[-1].Fill(z, 1.)
        histAlist[-1].Fill(A, 1.)
        histSOlist[-1].Fill(so, 1.)
        histRlist[-1].Fill(r, 1.)
        histRZlist[-1].Fill(r, z, 1./histRZnorm(r, binWidthR, binWidthZ))

    for i in range(N+1):
        histZlist[i].Write('', ROOT.TObject.kOverwrite)
        histAlist[i].Write('', ROOT.TObject.kOverwrite)
        histSOlist[i].Write('', ROOT.TObject.kOverwrite)
        histRlist[i].Write('', ROOT.TObject.kOverwrite)
        histRZlist[i].Write('', ROOT.TObject.kOverwrite)

    f.Close()

# === NORMALIZATION ===
def circleHexagonIntersectArea(r, a):
    return 12 * (0.5*a*np.sqrt(r**2 - a**2) + 0.5*r**2*(np.pi/6. - np.arctan(np.sqrt((r/a)**2 - 1))))

# Normalize r-z-histogram to density. Therefore, volume
# of a tubular bin element has to be determined. Here,
# the intersection of a circle and a hexagon have to
# be considered
def histRZnorm(r, binWidthR, binWidthZ, a=162):
    # Maximum radius for apothem a
    maxR = a * 2./np.sqrt(3)

    # get lower edge of bin
    rl = (r // binWidthR) * binWidthR

    if (rl <= a) and (rl + binWidthR <= a):
        return np.pi*binWidthZ*((rl + binWidthR)**2 - rl**2)
    elif (rl <= a) and (rl + binWidthR >= a):
        return (circleHexagonIntersectArea(rl+binWidthR, a) - np.pi*rl**2) * binWidthZ
    elif (rl > a) and (rl + binWidthR <= maxR):
        return (circleHexagonIntersectArea(rl+binWidthR, a) - circleHexagonIntersectArea(rl, a)) * binWidthZ
    elif (rl > a) and (rl + binWidthR > maxR):
        return (circleHexagonIntersectArea(rl+binWidthR, a) - circleHexagonIntersectArea(maxR, a)) * binWidthZ

# === STATISTICS ===
def chiSquareTest(data1, data2):
    import scipy.stats
    N1, N2 = sum( data1 ), sum( data2 )
    print N1, N2
    print np.array( data1 ) / N1, np.array( data2 ) / N2

    dof = len( data1 )
    chiSum = 0
    for i in range( len(data1) ):
        if not data1[i] and not data2[i]:
            dof -= 1
            continue
        chiSum += (data1[i]/N1 - data2[i]/N2)**2 / (data1[i]/N1**2 + data2[i]/N2**2)

    return chiSum, scipy.stats.chisqprob(chiSum, dof-1)

def geometricTest(data1, data2):
    N1, N2 = sum( data1 ), sum( data2 )
    return sum( [data1[i]*data2[i]/(N1*N2) for i in range( len(data1) )] )**.5

def ksTest(data1, data2, alpha=0.05):
    N1, N2 = sum( data1 ), sum( data2 )
    data1, data2 = np.array(data1)/N1, np.array(data2)/N2

    def cumSum(l, i):
        return sum(l[:i]) 

    cumSumData1 = np.array( [cumSum(data1, i) for i in range( len(data1) )] )
    cumSumData2 = np.array( [cumSum(data2, i) for i in range( len(data2) )] )

    # Kolmogorov-Smirnov statistic
    D = max( abs( cumSumData1 - cumSumData2 ) )

    c = np.sqrt(-0.5 * np.log(alpha/2.))
    r = c * np.sqrt((N1 + N2) / (N1*N2))

    # Null hypothesis is rejected if D > r
    return D, r

