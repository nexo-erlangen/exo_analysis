import sys
import os
import numpy as np
import math as m
import pprint
import time
import warnings
from collections import deque
import multiprocessing
import logging
# logging.basicConfig(level=logging.DEBUG, format='%(levelname)s - %(message)s')
# logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')

from plot_support import *
import plot_functions as pf
import vtk_to_potential as vtp
import peak_finder as peak
import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')

# randomProcess
# =============
def randomProcess(randPos, tree, H_BINS_RAND, H_BINS, Apo, name, N, pot=False, potName=None, verbose=False):
	start = time.time()

	# Create random events using statistics of histogram
	# hist and save the vertices in a tree.
	# Create ROOT-file for the case RAM runs out
	# fRandTree = ROOT.TFile('/home/vault/capm/mppi025h/analysis/rand_tree.root', 'RECREATE')
	
	# ROOT.gROOT.cd()
	# tree = ROOT.TTree('random', 'random')
	# randPos = generateRandom(hist, tree, N)
	# fRandTree.Close()

	# Read the electrical field from file. The data is
	# stored in three TH3D histograms, one for each
	# vector component of the Efield. 

	#fEfield = ROOT.TFile('~/sfepy/efield_hist.root', 'READ')
	fEfield = ROOT.TFile('efield_hist.root')

	EFieldData = [ fEfield.Get('Ex'), fEfield.Get('Ey'), fEfield.Get('Ez') ]

	# Drift N electrons stored in tree along the Efield 
	# specified in EfieldData. An electron is drifted
	# within a bin, therefore the size of the bin stored
	# in H_BINS_RAND is needed.
	posHistory, posColl = driftElectronMulti(tree, EFieldData, H_BINS_RAND, N, pot, potName)

	fEfield.Close()

	end = time.time()
	print 'Elapsed time: %s s' % str(end - start)

	# Get last x, y and first z component of every
	# drifted particle in list posHistory
	randDriftPos = []
	for item in posHistory:
		x = item[-1][0]
		y = item[-1][1]
		z = item[0][2]

		randDriftPos.append( [x, y, z] )

	# random data after cut
	fiducialVerts = isFiducialList(randPos, *Apo)
	randHist = fillRandHist(fiducialVerts, H_BINS, name + ' (randHist)')
	randStandoff = getStandoffList(fiducialVerts, name + ' (randStandoff)')
	randHist.Scale(1./randHist.Integral())
	randStandoff.Scale(1./randStandoff.Integral())

	# random data after drift and cut
	fiducialDriftVerts = isFiducialList(randDriftPos, *Apo)	
	randHistDrift = fillRandHist(fiducialDriftVerts, H_BINS, name + ' (randHistDrift)')
	randDriftStandoff = getStandoffList(fiducialDriftVerts, name + ' (randDriftStandoff)')
	randHistDrift.Scale(1./randHistDrift.Integral())
	randDriftStandoff.Scale(1./randDriftStandoff.Integral())

	f = ROOT.TFile.Open('standoffDrift.root', 'UPDATE')
	randStandoff.Write()
	randDriftStandoff.Write()
	f.Close()

	# particles hitting the teflon 
	if pot and posColl:
		posCollCut = []
		for item in posColl:
			x = item[-1][0]
			y = item[-1][1]
			z = item[0][2]

			posCollCut.append( [x, y, z] )

		fiducialCollVerts = isFiducialList(posCollCut, *Apo)
		collHist = fillRandHist(fiducialCollVerts, H_BINS, name + ' (collHist)')
		collStandoff = getStandoffList(fiducialCollVerts, name + ' (collStandoff)')
		try:
			collHist.Scale(1./collHist.Integral())
			collStandoff.Scale(1./collStandoff.Integral())
		except:
			pass

		f = ROOT.TFile.Open('collision.root', 'UPDATE')
		collHist.Write()
		collStandoff.Write()
		f.Close()

	if verbose:
		plotParticleTrace(posHistory, 100)

		pf.plotCompTwo(randStandoff, randDriftStandoff, 'Standoff Distance [m]', 'Events / (0.01 m)', 'Random', 'Random drifted', True)

		c = ROOT.TCanvas()
		randHist.Project3D('yx').Draw('colz')
		c1 = ROOT.TCanvas()
		randHistDrift.Project3D('yx').Draw('colz')
		#c.SaveAs('randHist.pdf') 
		raw_input('end')

# driftElectronMulti
# ==================
# Drift N electrons stored in tree in the field
# specified by EFieldData. Returns the track history
# of every particle in a list posHistory.
def driftElectronMulti(tree, EFieldData, H_BINS_RAND, N, pot=False, potName=None):
	# Running the following loop multiple times is
	# faster than limiting the number of proccesses
	# run at the same time. Reason for this is not
	# known.
	# LEN specifices the number of threads to be run
	# at the same time. 

	# Set pointer from tree branches to coordinates
	# x, y, z. Access using tree.GetEntry(i) and x[0].
	x, y, z = readTreeToCoord(tree)
	
	if pot:
		potential = vtp.localRBF(*vtp.vtkToVec(potName))
		# potential = vtp.vtkToPotential(potName)

	print 'Drifting particles...'
	if N > 10:
		LEN = 10 #int(float(N)/STEP)
	else:
		LEN = N

	STEP = int( float(N)/LEN ) #10
	if not STEP:
		STEP = 1

	for j in range(STEP):
		manager = multiprocessing.Manager()
		out_q = manager.Queue() #multiprocessing.Queue()
		# pool = multiprocessing.Pool() #[]

		procs = []
		for i in range(j*LEN, (j+1)*LEN):
			tree.GetEntry(i)
			xVec = [x[0], y[0], z[0]]

			#pool.apply_async(driftElectron, args=(xVec, EFieldData, H_BINS_RAND, i, out_q,))
			if not pot:
				p = multiprocessing.Process(target=driftElectron, args=(xVec, EFieldData, H_BINS_RAND, i, out_q,))
			else:
				#print 'Starting Multiprocess'
				#print potName
				p = multiprocessing.Process(target=driftElectronPotential, args=(xVec, potential, H_BINS_RAND, i, out_q,))
			procs.append(p)
			p.start()

		# Dictionary containing the particle tracks.
		# Index specifies the particle number from
		# 0 to N.
		resultdict = {}
		# print 'Obtaining tracks for task %d...' % j
		for i in range(LEN): #range(j*LEN, (j+1)*LEN):
			try:
				resultdict.update(out_q.get())
			except:
				break

		p = int(100.*(j+1)/STEP)
		statusBar(p, 'screen')

		for p in procs:
			# print len(multiprocessing.active_children())
			p.join()

		#pool.close()
		#pool.join()

		# Read out the events from the dictionary.
		# Here, the try - except routine is needed because
		# if a collision took place, the entry in posHistory
		# with the same number doesn't exist.
		posHistory = []
		posColl = []

		for i in range(j*LEN, (j+1)*LEN):
			try:
				posHistory.append(resultdict[i])
			except:
				pass

		if pot:
			for i in range(j*LEN, (j+1)*LEN):
				try:
					posColl.append(resultdict[-i])
				except:
					pass
	print

	return posHistory, posColl

# generateRandom
# ==============
# Generate N random events following the distribution
# of histogram hist. Save events in specfied TTree tree
# and return a list of the coordinates of the generated
# events . 
def generateRandom(hist, tree, N):
	r_RID = 0.1832356		# m

	# important to set them to 0. 
	# otherwise the entries will have the 
	# wrong data type (int)
	x = np.array([0.])
	y = np.array([0.])
	z = np.array([0.])

	# randHist = ROOT.TH3F('random', 'random', *H_BINS)
	tree.Branch('x', x, 'x/D')
	tree.Branch('y', y, 'y/D')
	tree.Branch('z', z, 'z/D')

	randPos = []
	print 'Randomly generating %d particles...' % N
	for i in range(N):
		p = int(100.*(i+1)/N)
		statusBar(p, 'screen')

		hist.GetRandom3(x, y, z)
		# convert from mm to m
		x[0] *= 1.e-3
		y[0] *= 1.e-3
		z[0] *= 1.e-3

		# It is possible that events are generated
		# outside of the detector volume due to the
		# binning
		if (x[0]**2 + y[0]**2) > r_RID**2:
			continue

		#randHist.Fill(x[0], y[0], z[0])
		randPos.append( [x[0], y[0], z[0]] )
		tree.Fill()
	print 

	return randPos

def readTreeToCoord(tree):
	x = np.array([0.])
	y = np.array([0.])
	z = np.array([0.])

	tree.SetBranchAddress('x', x)
	tree.SetBranchAddress('y', y)
	tree.SetBranchAddress('z', z)

	return x, y, z

def driftElectronPotential(startPos, potential, H_BINS, eventNum=0, out_q=None):
	#print 'Starting Drift'

	v_DRIFT = 0.00171e6		# m/s
	t_STHB = 0.05e-6		# s
	x_CAD = 19.223657e-2	# m
	x_APD = 20.44065e-2		# m
	r_RID = 0.1832356		# m

	t = 0.005
	posHistory = [startPos]

	xVec = startPos
	dq = deque([0]*5)
	zLast = None
	collEvent = False

	while abs(xVec[2]) < x_CAD:
		if (xVec[0]**2 + xVec[1]**2) >= r_RID**2:
			posHistory.append( xVec ) 
			collEvent = True
			#print 'Coll'
			break

		if xVec[2] > 0.:
			xVec[2] = x_CAD
			posHistory.append( xVec )
			#print 'Positive'
			break 

		# grad = np.array( vtp.getGradient(xVec[0], xVec[1], xVec[2], potential) )
		grad = np.array( potential.getGradient(xVec, 0.05) )
		gradMag = np.sqrt( sum( [item**2 for item in grad] ) )
		grad /= gradMag
		#print 'Position:', xVec
		#print 'Gradient:', grad

		dq[0] = grad[2]
		# print dq
		meandq = sum(dq)/len(dq)
		# print grad[2], xVec[2]

		if abs( round(meandq, 3) ) == 1 or (zLast > 0 and grad[2] < 0) or (zLast < 0 and grad[2] > 0):
			# print 'zLast triggered'
			if xVec[2] < 0:
				xVec[2] = -x_CAD
			else:
				xVec[2] = x_CAD

			posHistory.append( list(xVec) )
			#print 'break'
			break

		else:
			dq.rotate(1)
			zLast = grad[2]
			xVec += grad * t # * v_DRIFT
			posHistory.append( list(xVec) )

	# print 'Have a Kit Kat'
	print posHistory
	if out_q:
		outdict = {}
		if collEvent:
			outdict[-eventNum] = posHistory
		else:
			outdict[eventNum] = posHistory
		out_q.put(outdict)
	else:
		return posHistory

def driftElectron(startPos, fEfieldData, H_BINS, eventNum=0, out_q=None):
	v_DRIFT = 0.00171e6		# m/s
	t_STHB = 0.05e-6		# s
	x_CAD = 19.223657e-2	# m
	x_APD = 20.44065e-2		# m

	bins, binMin, binMax, binWidth = getBinData(H_BINS)

	xVec = np.array(startPos)
	posHistory = [startPos]

	# Get the origin of the histogram bin for the first event
	cubeOrigin = np.array( [ roundBase(item, 10, binWidth[i]) for i, item in enumerate(xVec) ] )

	# Later used to fix collision problem
	topColl = False
	lastColl = [(False, False)] * 3

	breakCnt = 0
	while(True):
		if breakCnt > 100:
			if xVec[2] > 0.:
				xVec[2] = x_CAD
			else:
				xVec[2] = -x_CAD
			break

		if xVec[2] > 0.:
			xVec[2] = x_CAD
			posHistory.append( xVec )
			break
		
		Evec = getEvec(fEfieldData, xVec, H_BINS)
		# Subtract the bin origin from the particle position
		# to transform it to coordinates within the bin
		xVecCub = xVec - cubeOrigin

		#for i in range(3):
		#	if xVecCub[i] < 0.:
		#		xVecCub[i] += binWidth[i]

		logging.debug('CubeOrigin: %s' % pprint.pformat(cubeOrigin))
		logging.debug('\tEvec: %s' % pprint.pformat(Evec))
		logging.debug('\txVec: %s' % pprint.pformat(xVec))
		logging.debug('\txVecCub: %s' % pprint.pformat(xVecCub))

		newPosCub, coll = vectorBoxIntersection(xVecCub, Evec, binWidth)
		logging.debug('\tnewPosCub: %s' % pprint.pformat(newPosCub))

		# Check if bottom collision is followd by top
		# collision or vice versa. If this is the case,
		# the routine is traped. To solve the situation,
		# Add the vectors of the conflicting cells and
		# in the resulting direction
		# print 'Collision error!'

		'''
		if (coll[4] and topColl) or (coll[5] and not topColl):
			#print 'Evec:', Evec
			#print 'topColl:', topColl
			#print 'coll:', coll
			if coll[4]:
				zJump = -binWidth[2]
			else:
				zJump = binWidth[2]

			xJump = np.array(xVec) + np.array( [0, 0, zJump] )

			Evec = 0.5 * ( np.array(Evec) + np.array( getEvec(fEfieldData, xJump , H_BINS) ) )
			newPosCub, coll = vectorBoxIntersection(xVecCub, Evec, binWidth)
			cubeOrigin = np.array(cubeOrigin) + xJump
		'''

		xVec = newPosCub + cubeOrigin

		logging.debug('\tnewPos: %s' % pprint.pformat(xVec))
		logging.debug('')

		# Append new position to history
		posHistory.append(xVec)

		# Check if drift completed
		if(abs(xVec[2]) > x_CAD):
			xVec[2] = -x_CAD
			posHistory.append(xVec)
			break

		# Check collision to shift bin origin
		coll = [ (coll[2*i], coll[2*i+1]) for i in range(3) ]	# split list in list of tuples

		# print lastColl
		# print coll

		#if (np.array(lastColl)==np.array(coll)).all:
		#	break

		checkColl = np.array(lastColl) & np.array(coll)
		# print coll
		# i - index for coordinates
		for i, item in enumerate(checkColl):
			if item[0] and (Evec[i] < 0):
				cubeOrigin[i] -= binWidth[i]
			if item[1] and (Evec[i] > 0):
				cubeOrigin[i] += binWidth[i]

		lastColl = coll
		# Store last collision result for top/bottom
		# top/bottom = True/False
		tbColl = coll[-1]
		if tbColl[1]:
			topColl = True
		else:
			topColl = False

		breakCnt += 1
		
	if out_q:
		outdict = {}
		outdict[eventNum] = posHistory
		out_q.put(outdict)
	else:
		return posHistory

def getEField(fEfieldData, pos, H_BINS):
	bins, binMin, binMax, binWidth = getBinData(H_BINS)

	#EFieldX = fEfield.Get('Ex')
	#EFieldY = fEfield.Get('Ey')
	#EFieldZ = fEfield.Get('Ez')
	#EFieldList = [EFieldX, EFieldY, EFieldZ]

	pos = np.array( [ roundBase(item, prec=5, base=binWidth[i]) for i, item in enumerate(pos) ] )
	posBin = [ int( item ) for item in (pos - binMin) / binWidth ]
	x, y, z = posBin

	Evec = []
	for EEntry in fEfieldData:
		# If there is no bin at position x, y, z
		# an error is raised
		try:
			Evec.append( EEntry.GetBinContent(x, y, z) )
		except:
			Evec.append(0.)

	Ex, Ey, Ez = Evec
	# Ex = 0
	# Ey = 0
	# Ez = 1
	return Ex, Ey, Ez

def getEvec(fEfieldData, xVec, H_BINS):
	# Get the Efield vector at the current particle position
	Evec = list(getEField(fEfieldData, xVec, H_BINS))
	if not Evec[0] and not Evec[1] and not Evec[2]:
		if xVec[2] >= 0.:
			Evec[2] = 1.
		else:
			Evec[2] = -1.

	# Get Magnitude of Efield vector and normalize the vector
	Emag = np.sqrt( np.sum([item**2 for item in Evec]) )
	Evec = [ item/Emag for item in Evec ]
	return Evec

# ==== vectorBoxIntersection ====
# Get the coordinates of the intersection of vector vec
# located at position pos, defining a line in space, 
# and the bin boundaries binW. 
# Returns also which collision took place in the format
# coll = [Left, Right, Front, Back, Bottom, Top]
def vectorBoxIntersection(pos, vec, binW):
	pos = np.array(pos)
	vec = np.array(vec)
	binW = np.array(binW)

	tVec0 = []
	tVecW = []
	for i in range(3):
		if vec[i]:
			tVec0.append( (0 - pos[i]) / vec[i] )
			tVecW.append( (binW[i] - pos[i]) / vec[i] )
		else:
			tVec0.append( float('nan') )
			tVecW.append( float('nan') )

	with warnings.catch_warnings():
   		warnings.filterwarnings('ignore', r'All-NaN (slice|axis) encountered')
		test = np.nanmax([tVec0, tVecW], axis=0)
		testMin = np.nanmin([tVec0, tVecW], axis=0)
		fPos = pos + vec * np.nanmin(test)
		#if (fPos==pos).all:
		#	fPos = pos + vec * np.nanmin(tVec0)
		#	print 'Correction:', fPos
		#	print tVec0
		#	print tVecW

	coll = []
	for i, bw in enumerate(binW):
		if fPos[i] == 0:
			coll += [True, False]
		elif fPos[i] == bw:
			coll += [False, True]
		else:
			coll += [False, False]

	# print coll
	return fPos, coll

def getBinData(H_BINS):
	bins = np.array( [ H_BINS[0], H_BINS[3], H_BINS[6] ] )
	binMin = np.array( [ H_BINS[1], H_BINS[4], H_BINS[7] ] )
	binMax = np.array( [ H_BINS[2], H_BINS[5], H_BINS[8] ] )
	binWidth = ( binMax - binMin ) / bins 

	return bins, binMin, binMax, binWidth

def plotParticleTrace(posHistory, N=False):
	c = ROOT.TCanvas()
	view = ROOT.TView3D.CreateView(1)
	view.Centered3DImages()
	view.ShowAxis()
	view.SetRange(-0.2,-0.2,-0.2,0.2,0.2,0.2)
	l = np.array([ROOT.TPolyLine3D()]*len(posHistory))

	if len(posHistory) < N or not N:
		N = len(posHistory)

	for i in range(N):
		item = posHistory[i]
		item = np.array(item)
		# print item

		x = np.array(item[:,0])
		y = np.array(item[:,1])
		z = np.array(item[:,2])
		l[i] = ROOT.TPolyLine3D(len(item), x, y, z)
		l[i].Draw()
	c.Update()
	raw_input('end')

# fillRandHist
# ============
# history includes the track positions of all
# drifted particles. Because they arrived at the
# anode only the last entry of each particle
# excluding the z-component is needed

def fillRandHist(history, H_BINS, name='default'):
	H_BINS_COR = []

	# Change bin dimension: mm -> m
	for i in range(3):
		H_BINS_COR += [ H_BINS[3*i], H_BINS[3*i+1]/1000., H_BINS[3*i+2]/1000. ]

	randHist = ROOT.TH3D(name, name, *H_BINS_COR)

	for item in history:
		randHist.Fill(*item)

	return randHist

# rootToHist
# ==========
# Open root file fName and read stored histograms
# to dictionaries. Compare drifted with regular 
# data and store the results in a histogram. 
# The bin settings of the histogram are taken from
# spacing, which contains three np.linspaces used
# to set the parameters for drifting.
# ScaleBin specifies from which bin the integral
# starts, later used for normalisation.
def rootToHist(fNameRand, fNameTrue, spacing, scaleBin=5):
	# ==== RANDOM DATA ====
	# Get bin settings
	H_BINS = []
	for item in spacing:
		bins = len(item)
		binMin = item[0]
		binMax = item[-1]

		H_BINS += [bins, binMin, binMax]

	standoffHist = ROOT.TH3D('standoffStat', 'standoffStat', *H_BINS)

	# Separate drifted from regular data
	fRand = ROOT.TFile.Open(fNameRand)
	hDict = openAllHistos(fRand)

	hDictMC = {}
	hDictMCDrift = {}
	for key in hDict:
		nameList = []
		entry = key.replace('(', ' ')
		for elem in entry.split(' '):
			try:
				nameList.append( float(elem.split(',')[0]) )
			except:
				pass

		print nameList
		sigma0, theta, z = nameList
		if 'Drift' in key:
			hDictMCDrift[tuple(nameList)] = hDict[key]
		else:
			hDictMC[tuple(nameList)] = hDict[key]

	# ==== TRUE DATA ====
	fTrue = ROOT.TFile.Open(fNameTrue)
	hDict = openAllHistos(fTrue)
	hDataTrue = hDict['DataSS']
	hMcTrue = hDict['McSS']
		
	hDataTrue.Scale(1./hDataTrue.Integral(scaleBin, 20))
	hMcTrue.Scale(1./hMcTrue.Integral(scaleBin, 20))

	hDiffTrue = getHistoDiffInv(hDataTrue, hMcTrue, 'difference (true data)')
	hDiffTrue.SetBins(20, 0, 0.2)

	c2 = ROOT.TCanvas()
	#hDiffTrue.Draw('HIST')
	pf.plotCompTwo(hMcTrue, hDataTrue, 'Standoff Distance [m]', 'Events / ( 0.1 m )', 'Regular', 'Drifted', True)

	# ==== COMPARISON ====
	for key in hDictMC:
		hMC = hDictMC[key]
		hMCDrift = hDictMCDrift[key]

		hMC.Scale(1./hMC.Integral(scaleBin, 20))
		hMCDrift.Scale(1./hMCDrift.Integral(scaleBin, 20))

		hDiff = getHistoDiffInv(hMCDrift, hMC, 'difference')

		hDiffDiff = getHistoDiffInv(hDiff, hDiffTrue, 'hDiffDiff')

		rmsDiff = getRMSFromHist(hDiffDiff, 5)
		
		print np.square( np.array([hDiffDiff.GetBinContent(i) for i in range(1, 5 + 1)]) )
		print rmsDiff
		print

		if rmsDiff < 10:
			pf.plotCompTwo(hMC, hMCDrift, 'Standoff Distance [m]', 'Events / ( 0.1 m )', 'Regular', 'Drifted', True)

		#c = ROOT.TCanvas()
		#hDiffDiff.Draw('HIST')
		#hDiffTrue.SetMarkerStyle(3)
		#hDiffTrue.Draw('same')
		#raw_input('end')
		#del c

		#hMC.Scale(10000)
		#hMCDrift.Scale(10000)
		# pf.plotCompQQ(hMC, hMCDrift, [20, 0, 0.2], 'MC', 'MC drifted', True)

	fRand.Close()
	fTrue.Close()

def getRMSFromHist(hist, binN):
	rms = np.sqrt( sum( np.square( np.array([hist.GetBinContent(i) for i in range(1, binN + 1)]) ) ) )
	return rms

# ==== SUPPORT FUNCTIONS ====
def roundBase(x, prec=5, base=.5):
	# with int: always round down
	# with round: correct rounding
	if x >= 0.:
		return round(base * int(float(x)/base), prec)
	else:
		return np.floor(float(x) / base) * base
		# return round(base * int(float(x - base/2.)/base), prec)

def isFiducialList(ls, Apo, zMin, zMax):
	fidList = []
	for item in ls:
		x, y, z = item
		if isFiducial(x, y, z, Apo, zMin, zMax):
			fidList.append(item)
	return fidList

def getStandoffList(ls, name):
	standoffHist = ROOT.TH1F(name, name, 20, 0, 0.2)

	for item in ls:
		standoffHist.Fill( getStandoff(*item) )

	return standoffHist

# ==== TEST RANGE ====
def randomTest(H_BINS):
	fEfield = ROOT.TFile('efield_hist.root', 'READ')
	EFieldData = [ fEfield.Get('Ex'), fEfield.Get('Ey'), fEfield.Get('Ez') ]

	startPos = REFLECTORINNERRAD/1000. - 0.005
	endPos = REFLECTORINNERRAD/1000. - 0.001

	pot = vtp.vtkToPotential('exo200_wo_edges.vtk')
	posHistoryPot = vtp.driftParticle(pot, startPos, endPos, 20)

	f = open('raycastTest.dat', 'w')
	posHistoryTotal = []
	for x in np.linspace(startPos, endPos, 20):
		startPos = [x, 0., -0.01]
		print startPos
		posHistory = driftElectron(startPos, EFieldData, H_BINS, eventNum=0, out_q=None)

		for item in posHistory:
			f.write('%f\t%f\t%f\n' % (item[0], item[1], item[2]) )
		f.write('\n')
		posHistoryTotal.append(posHistory)

	f.close()

	c = ROOT.TCanvas()
	mgApprox = dataToGraph(posHistoryTotal, 'Approximation', 38)
	mg = dataToGraph(posHistoryPot, 'Potential')
	mg.Add(mgApprox)

	mg.Draw('AL')
	c.Update()
	raw_input('end')

	print list(reversed( posHistoryPot[-1] ))

	x = np.array( list(reversed( posHistoryPot[-1] )) )[:,2]*1000
	y = np.array( list(reversed( posHistoryPot[-1] )) )[:,0]
	peak.findPeaks(x, y, title='EField Drift')

	gammaSineFunc = ROOT.TF1('gammaSine', gammaSine, -160., -20, 6)
	gammaSineFunc.SetParName(0, 'alpha')
	gammaSineFunc.SetParName(1, 'theta')
	gammaSineFunc.SetParName(2, 'x0')
	gammaSineFunc.SetParName(3, 'a')
	gammaSineFunc.SetParName(4, 'b')
	gammaSineFunc.SetParName(5, 'A')

	gammaSineFunc.SetParameter(0, 2.)
	gammaSineFunc.SetParameter(1, 2.)
	gammaSineFunc.SetParameter(2, -80)
	gammaSineFunc.SetParameter(3, 0.001)
	gammaSineFunc.SetParameter(4, 0.001)
	gammaSineFunc.SetParameter(5, 1.)

	c = ROOT.TCanvas()
	gammaSineFunc.Draw('C')
	raw_input('end')

def gammaSine(x, par):
	return par[5] * ROOT.Math.fdistribution_pdf(-x[0], par[0], par[1], par[2]) + 0 * par[3]*np.sin(x[0] - par[4])

def dataToGraph(posHistory, name='default', color=1):
	mg = ROOT.TMultiGraph()
	mg.SetName(name)
	l = [0] * len(posHistory)
	for i, item in enumerate(posHistory):
		l[i] = ROOT.TGraph(len(item), np.array( np.array(item)[:,0] ), np.array( np.array(item)[:,2] ))
		l[i].SetLineColor(color)
		mg.Add(l[i])

	return mg
