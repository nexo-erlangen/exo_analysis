import sys
import os
import re
import cPickle
import numpy as np
import math as m
import pprint
import time
import warnings
import random
from collections import deque
import multiprocessing
import logging
from pympler import summary, muppy
# logging.basicConfig(level=logging.DEBUG, format='%(levelname)s - %(message)s')
# logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')

from plot_support import *
import plot_functions as pf
import vtk_to_potential as vtp
import vtk_to_vec as vtv
import projectFit as fit
import peak_finder as peak
import csv_to_vec as ctv
import artificialField as af
import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')

# randomProcess
# =============
# randPos used to get standard distribution;
# tree is used to generate drifted distribution
def randomProcess(randPos, tree, H_BINS_RAND, H_BINS, Apo, name, N, outName='standoffDrift.root', efield=False, efieldName=None, pot=False, potName=None, TwoDFlag=False, verbose=False):
	start = time.clock()

	# Create random events using statistics of histogram
	# hist and save the vertices in a tree.
	# Create ROOT-file for the case RAM runs out
	#fRandTree = ROOT.TFile('/home/vault/capm/mppi025h/analysis/rand_tree.root', 'RECREATE')
	
	#ROOT.gROOT.cd()
	#tree = ROOT.TTree('random', 'random')
	# randPos = generateRandom(hist, tree, N)
	# fRandTree.Close()

	# Read the electrical field from file. The data is
	# stored in three TH3D histograms, one for each
	# vector component of the Efield. 

	if efield:
		#fEfield = ROOT.TFile('~/sfepy/efield_hist.root', 'READ')
		# fEfield = ROOT.TFile('efield_hist.root')
		fEfield = ROOT.TFile(efieldName)

		EFieldData = [ fEfield.Get('Efield_histX'), fEfield.Get('Efield_histY'), fEfield.Get('Efield_histZ') ]
	else:
		EFieldData = None

	# Drift N electrons stored in tree along the Efield 
	# specified in EfieldData. An electron is drifted
	# within a bin, therefore the size of the bin stored
	# in H_BINS_RAND is needed.
	posHistory = driftElectronMulti(tree, efield, EFieldData, H_BINS_RAND, N, pot, potName, TwoDFlag)

	if efield:
		fEfield.Close()

	end = time.clock()
	print 'Elapsed time: %s s' % str(end - start)

	# Get last x, y and first z component of every
	# drifted particle in list posHistory
	randDriftPos = []
	for item in posHistory:
		x = item[-1][0]
		y = item[-1][1]
		z = item[0][2]

		randDriftPos.append( [x, y, z] )

	# random data wo cut
	randHist = fillRandHist(randPos, H_BINS, name + ' (randHist)')
	randStandoff = getStandoffList(randPos, name + ' (randStandoff)')
	randHist.Scale(1./randHist.Integral())
	randStandoff.Scale(1./randStandoff.Integral())

	# random data after cut
	fiducialVerts = isFiducialList(randPos, *Apo)
	randHistCut = fillRandHist(fiducialVerts, H_BINS, name + ' (randHistCut)')
	randStandoffCut = getStandoffList(fiducialVerts, name + ' (randStandoffCut)')
	randHistCut.Scale(1./randHistCut.Integral())
	randStandoffCut.Scale(1./randStandoffCut.Integral())

	# random data after drift wo cut
	randHistDrift = fillRandHist(randDriftPos, H_BINS, name + ' (randDriftHist)')
	randDriftStandoff = getStandoffList(randDriftPos, name + ' (randDriftStandoff)')
	try:
		randHistDrift.Scale(1./randHistDrift.Integral())
		randDriftStandoff.Scale(1./randDriftStandoff.Integral())
	except:
		pass

	# random data after drift and cut
	fiducialDriftVerts = isFiducialList(randDriftPos, *Apo)	
	randHistDriftCut = fillRandHist(fiducialDriftVerts, H_BINS, name + ' (randDriftHistCut)')
	randDriftStandoffCut = getStandoffList(fiducialDriftVerts, name + ' (randDriftStandoffCut)')
	try:
		randHistDriftCut.Scale(1./randHistDriftCut.Integral())
		randDriftStandoffCut.Scale(1./randDriftStandoffCut.Integral())
	except:
		pass

	# Write list to file
	outNameSplit = outName.split('/')
	if outNameSplit[0] == outName:
		outNameDir = './'
	else:
		outNameDir = outNameSplit[0] + '/'
		
	#try:
	print outNameDir
	print outName
	# print re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", outName)
	try:
		outNameNum = int( float( re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", outName)[0] ) )
	except:
		outNameNum = 0

	cPickle.dump(randPos, open(outNameDir + 'randPos_%d.p' % outNameNum, 'wb'))
	cPickle.dump(randDriftPos, open(outNameDir + 'randDriftPos_%d.p' % outNameNum, 'wb'))
	cPickle.dump(fiducialVerts, open(outNameDir + 'randPosCut_%d.p' % outNameNum, 'wb'))
	cPickle.dump(fiducialDriftVerts, open(outNameDir + 'randDriftPosCut_%d.p' % outNameNum, 'wb'))

	# Write histograms to root-file
	f = ROOT.TFile.Open(outName, 'RECREATE')
	
	randHist.Write()
	randHistCut.Write()
	randStandoff.Write()
	randStandoffCut.Write()

	randHistDrift.Write()
	randHistDriftCut.Write()
	randDriftStandoff.Write()
	randDriftStandoffCut.Write()
	
	f.Close()

	saveParticleTrace('particleTrace.dat', posHistory)
	
	if verbose:
		# plotParticleTrace(posHistory, 100)

		pf.plotCompTwo(randStandoffCut, randDriftStandoffCut, 'Standoff Distance [m]', 'Events / (0.01 m)', 'Random', 'Random drifted', True)

		c = ROOT.TCanvas()
		randHistCut.Project3D('yx').Draw('colz')
		c1 = ROOT.TCanvas()
		randHistDriftCut.Project3D('yx').Draw('colz')
		#c.SaveAs('randHist.pdf') 
		raw_input('end')

# driftElectronMulti
# ==================
# Drift N electrons stored in tree in the field
# specified by EFieldData. Returns the track history
# of every particle in a list posHistory.
def driftElectronMulti(tree, efield, EFieldData, H_BINS_RAND, N, pot=False, potName=None, TwoDFlag=False):
	if not H_BINS_RAND and efield:
		H_BINS_RAND = hist3dToH_BINS(EFieldData[0])
	# Running the following loop multiple times is
	# faster than limiting the number of proccesses
	# run at the same time. Reason for this is not
	# known.
	# LEN specifices the number of threads to be run
	# at the same time. 
	posHistory = []

	# Set pointer from tree branches to coordinates
	# x, y, z. Access using tree.GetEntry(i) and x[0].
	x, y, z = readTreeToCoord(tree)
	
	if pot:
		if TwoDFlag:
			posList = cPickle.load(open('efield_data/posList.p', 'rb'))
			vecList = cPickle.load(open('efield_data/vecList.p', 'rb'))

		else:
			posList = cPickle.load(open('posList.p', 'rb'))
			vecList = cPickle.load(open('vecList.p', 'rb'))

		potential = vtv.localRBF(posList, vecList)
		os.system('free -m')
		
		# potential = vtp.localRBF(*vtp.vtkToVec(potName))
		# potential = vtp.vtkToPotential(potName)

	print 'Drifting particles...'

	# View Memory contents
	#all_objects = muppy.get_objects()
	#sum1 = summary.summarize(all_objects)
	#summary.print_(sum1)
	
	if N > 40:
		LEN = 40 #int(float(N)/STEP)
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
			if efield and not pot:
				p = multiprocessing.Process(target=driftElectron, args=(xVec, EFieldData, H_BINS_RAND, i, out_q,))
			elif pot and not efield:
				if TwoDFlag:
					p = multiprocessing.Process(target=driftElectronPotential2d, args=(xVec, potential, None, i, out_q, False, False,))
				else:
					# print 'Starting Multiprocess'
					# print potName
					p = multiprocessing.Process(target=driftElectronPotential, args=(xVec, potential, H_BINS_RAND, i, out_q,))
			elif pot and efield:
				p = multiprocessing.Process(target=completeDrift, args=(xVec, potential, EFieldData, H_BINS_RAND, i, out_q, False))  

			procs.append(p)
			p.start()

		# Dictionary containing the particle tracks.
		# Index specifies the particle number from
		# 0 to N.
		resultdict = {}
		# print 'Obtaining tracks for task %d...' % j
		for i in range(LEN): #range(j*LEN, (j+1)*LEN):
			resultdict.update(out_q.get())

		p = int(100.*(j+1)/STEP)
		statusBar(p, 'screen')

		for p in procs:
			# print len(multiprocessing.active_children())
			p.join()

		#pool.close()
		#pool.join()

		for i in range(j*LEN, (j+1)*LEN):
			if resultdict[i]:
				posHistory.append(resultdict[i])
	print

	return posHistory

# generateRandomArt
# =================
# Generate N random events following a fitted distribution
def generateRandomArt(hist, tree, N):
	parList = (255.0, 64.96189303640774, 3.9, 18554.095519082686, -30.0, 104.85530141525176)

	x, y, z = np.array([0.]), np.array([0.]), np.array([0.])
	tree.Branch('x', x, 'x/D')
	tree.Branch('y', y, 'y/D')
	tree.Branch('z', z, 'z/D')

	print 'Randomly generating %d particles...' % N

	posList, h = fit.generateRandom(parList, N)	
	for pos in posList:
		x[0], y[0], z[0] = pos
		tree.Fill()
	print

	return posList

# generateRandomMC
# ================
# Get N cluster positions from real MC data
def generateRandomMC(hist, tree, N):
	x, y, z = np.array([0.]), np.array([0.]), np.array([0.])
	tree.Branch('x', x, 'x/D')
	tree.Branch('y', y, 'y/D')
	tree.Branch('z', z, 'z/D')

	randPos = []
	return

# generateRandom
# ==============
# Generate N random events following the distribution
# of histogram hist. Save events in specfied TTree tree
# and return a list of the coordinates of the generated
# events . 
def generateRandom(hist, tree, N, section='all'):
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

	# gRandom = ROOT.TRandom3(0)
	# ROOT.gRandom.SetSeed(random.randint(0, sys.maxint))
	ROOT.gRandom.SetSeed(0)

	# for i in range(N):
        i = 0
        j = N
        while j > 0:
		p = int(100.*(i+1)/N)
		statusBar(p, 'screen')

		hist.GetRandom3(x, y, z)
		if np.sqrt( x**2 + y**2 ) >= REFLECTORINNERRAD:
			continue

                if section == 'main':
                    theta = np.arctan2(y, x)
                    if not ( (abs(theta) <= np.radians(30)) or (abs(theta) >= np.radians(150)) ):
                        continue

                elif section == 'side':
                    theta = np.arctan2(y, x)
                    if ( (abs(theta) <= np.radians(30)) or (abs(theta) >= np.radians(150)) ):
                        continue

                # else: use whole volume

		# convert from mm to m
		x[0] *= 1.e-3
		y[0] *= 1.e-3
		z[0] *= 1.e-3
		#randHist.Fill(x[0], y[0], z[0])
		randPos.append( [x[0], y[0], z[0]] )
		tree.Fill()

                j -= 1
                i += 1
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

# ==== EXOAnalysis Functions ====
def getPotential_C(x, y, z, t):
	print 'Pot'
	posList = cPickle.load(open('/home/vault/capm/mppi025h/plot_scripts/efield_data/posList.p', 'rb'))
	vecList = cPickle.load(open('/home/vault/capm/mppi025h/plot_scripts/efield_data/vecList.p', 'rb'))
	potential = vtv.localRBF(posList, vecList)

	print 'Potential loaded'
	return potential

def driftElectronPotential2d_C(x, y, z, t, potential):
	startPos = list( np.array( [x, y, z] )/1000. )
	posHist, lengthCorrection = driftElectronPotential2d(startPos, potential, None, 0, None, False, False)

	if posHist:
		res = posHist[-1]
		x_, y_, z_ = list( np.array( res )*1000 )
                if z > 0:
			z__ = z - lengthCorrection
			if z__ < 0:
				x_, y_, z_ = 999., 999., 999.
			else:
				z_ = z__
		if z < 0:
			z__ = z + lengthCorrection
			if z__ > 0:
				x_, y_, z_ = 999., 999., 999.
			else:
				z_ = z__
	else:
		x_, y_, z_ = 999., 999., 999.

	return (x_, y_, z_, t)

# driftElectronPotential2d
# ========================
def driftElectronPotential2d(startPos, potential, H_BINS, eventNum=0, out_q=None, diffuse=False, saveHistory=False, useRK4=False, interpolate=False):
        from artificialField import reflectorRad as refRad

	# ATTENTION: Don't use saveHistory for large amount
	# of events, because memory will fill up quickly!

	x_CAD = 0.19223657
	r_RID = 0.1832356
	limit = 0.187

        if RK4:
            t = 0.0001
        else:
	    t = 0.000325

	xVec = list( startPos )
	posHistory = [ xVec ]
	# dq = deque([0]*5)
	zLast = None

	# ATTENTION: ONLY LOWER HALF OF DETECTOR IS USED
	upFlag = False
	if xVec[2] > 0:
		upFlag = True
		xVec[2] *= -1

	x, y, z = xVec
	theta = np.arctan2(y, x)
	r = np.sqrt( x**2 + y**2 )
        
        # Get reflector radius for angle theta
        reflectorRad2 = 0.04
        reflectorX2 = 0.215
        # r_RID = refRad(theta, r_RID, reflectorRad2, reflectorX2)

        # Height of the event in relation to limit
        lengthHeight = limit - abs(z)

	if diffuse:
		Ri = r
		Rvec = xVec[:2]

	xVecPol = np.array( [r, z] )

        length = 0
        while abs( xVecPol[1] ) < limit:
		# Use polar coordinates
                if xVecPol[0] > r_RID:
			posHistory = False
			print 'PTFE Hit!'
			break

		# print 'xVec:', xVec
		# print 'xVecPol:', xVecPol

                if useRK4:
                    grad = RK4(potential, [0., xVecPol[0], xVecPol[1]], 4*t, interpolate) 
                else:
                    if interpolate:
			grad = np.array( potential.getGradient([0., xVecPol[0], xVecPol[1]], t) )
                    else:
			grad = np.array( potential.getNearestGradient([0., xVecPol[0], xVecPol[1]]) )
    
		if np.isnan( grad ).any():
			posHistory = False
			break

		# print 'Grad:', grad
		gradMag = np.linalg.norm( grad ) # np.sqrt( sum( [item**2 for item in grad] ) )
		grad /= gradMag
		# print 'GradMag:', gradMag

		dr = grad[1]
		dz = grad[2]
		# polMag = np.sqrt( dr**2 + dz**2 )
		gradPol = np.array( [dr, dz] ) # /polMag
		# print 'GradPol:', gradPol

                '''
		dq[0] = dz
		meandq = sum(dq)/len(dq)

		if abs( round(meandq, 7) ) == 1:
			# print dq
			print
			print 'Field homogeneous'
			xVec = list( polToCart( xVecPol[0], xVecPol[1], theta ) )
			xVec[2] = -limit
			posHistory.append( xVec )
			break
                '''

                if False: # (zLast > 0 and gradPol[1] < 0) or (zLast < 0 and gradPol[1] > 0):
			print
			print 'zLast triggered'
			# NOTE: Currently the particle gets killed because
			# it probably is extremely close to the PTFE
			posHistory = False
			break

			# The z-component of the E-field flipped:
			# try to put the particle a bit beneath 
			# its current position and proceed the drift
			xVecPol[1] -= 0.1 * t
			if saveHistory:
				xVec = list( polToCart( xVecPol[0], xVecPol[1], theta ) )
				posHistory.append( xVec )

		else:
			# dq.rotate(1)
			zLast = dz
		        xVecPolOld = list( xVecPol )

			if not diffuse:
				xVecPol += np.array( gradPol ) * t 
			else:
				dx, dy, dz = getDiffusion( t )

				dR, dR_p = projectOnVector(dx, dy, Rvec[0], Rvec[1], Ri)
				dTheta = np.arctan2(dR_p, r) 

				'''
				print dx, dy, dz
				print dR, dR_p
				print dTheta
				print
				'''

				xVecPol += np.array( gradPol ) * t + np.array( [dR, dz] )
				theta += dTheta

			# print 'NewXVecPol', xVecPol
                        length += np.linalg.norm(xVecPol - xVecPolOld)

			if saveHistory:
				xVec = list( polToCart( xVecPol[0], xVecPol[1], theta ) )
				posHistory.append( xVec )
	else:
		if abs( xVecPol[1] ) > limit:
			# Fails if gradPol is not defined which is
			# the case for an event happening below limit
			try:
				deltaFrac = (abs(xVecPol[1]) - limit) / ( gradPol[1] * t )
				xVecPol += deltaFrac * t * np.array( gradPol )
                                length -= np.linalg.norm( deltaFrac * t * np.array( gradPol ) )

			except:
				pass

		posHistory.append( list( polToCart( xVecPol[0], xVecPol[1], theta ) ) )

	if posHistory and upFlag:
		posHistory = np.array( posHistory )
		posHistory[:,2] *= -1
		posHistory = posHistory.tolist()

	if out_q:
		# print posHistory
		outdict = {}
		outdict[eventNum] = posHistory
		out_q.put(outdict)
	else:
		return posHistory, abs(length - lengthHeight)

def polToCart(r, z, theta):
	x = r*np.cos( theta )
	y = r*np.sin( theta )

	return x, y, z

def getDiffusion(t):
	drift_velocity = 0.00171e6		# m / s 

	transDiffCoeff = 2.e-3			# m^2 / ns
	longDiffCoeff = 1.1e-2

	sigmaXY = np.sqrt( transDiffCoeff / drift_velocity * t )
	sigmaZ = np.sqrt( longDiffCoeff / drift_velocity * t )

	dx = np.random.normal(0, sigmaXY, 1)[0]
	dy = np.random.normal(0, sigmaXY, 1)[0]
	dz = np.random.normal(0, sigmaZ, 1)[0]

	return dx, dy, dz

# == projectOnVector ==
# Use the initial position of the particle
# to get xR, yR and R
def projectOnVector(x, y, xR, yR, R):
	# vec is vector to be projected
	vec = np.array( [x,  y] )

	# vecR points in direction of the radius
	vecR = 1/R * np.array( [xR, yR] )

	# vec_p is perpendicular to vec
	vecR_p = 1/R * np.array( [-yR, xR] )

	# projections
	projR = np.dot(vec, vecR) # * vecR 
	projR_p = np.dot(vec, vecR_p) # * vecR_p

	return projR, projR_p

# -- END - driftElectronPotential2d --

def driftElectronPotential(startPos, potential, H_BINS, eventNum=0, out_q=None, bottom=False, efield=False, saveHistory=False):
	# print 'Starting Drift'

	v_DRIFT = 0.00171e6		# m/s
	t_STHB = 0.05e-6		# s
	x_CAD = 19.53215e-2 # 19.223657e-2		# m
	x_APD = 20.44065e-2		# m
	r_RID = 0.180			# m
	CATHODE_ANODE_x_DISTANCE = 19.84065e-2
	CATHODE_ANODE_y_DISTANCE = 19.223657e-2

	WIRE_DIAMETER = 0.0127e-2

	t = v_DRIFT * t_STHB # 0.002
	
	xVec = startPos
	posHistory = [ xVec ]

	dq = deque([0]*10)
	zLast = None

	# ATTENTION: ONLY ONE HALF OF THE LOWER PART OF THE 
	# DETECTOR IS CURRENTLY USED
	leftFlag = False
	if xVec[0] < 0:
		leftFlag = True
		xVec[0] = -xVec[0]
	
	upFlag = False
	if xVec[2] > 0:
		upFlag = True
		xVec[2] = -xVec[2]

	# === MAIN LOOP ===
	if bottom:
		limit = CATHODE_ANODE_x_DISTANCE - WIRE_DIAMETER
	else:
		limit = x_CAD

	while abs(xVec[2]) < 0.199: # x_APD:
		# print xVec
		if (xVec[0]**2 + xVec[1]**2) >= r_RID**2:
			# posHistory.append( xVec ) 
			# collEvent = True
			# print 'Coll'

			posHistory = False
			break

		# grad = np.array( vtp.getGradient(xVec[0], xVec[1], xVec[2], potential) )

		if bottom:
			uvVec = [ -xVec[2]*1000 ] + list( uvPositionCorrection( *xy_to_uv( *(np.array(xVec[:2])*1000) ) ) )
			# print xVec, uvVec

			if efield:
				grad = getEvecMike(potential, uvVec, H_BINS)
				grad[0] *= -1
				grad[1] *= -1

				if not grad:
					break

			else:
				grad = np.array( potential.getGradient(uvVec, 0.05) )
		else:
			if efield:
				grad = getEvec(potential, xVec, H_BINS)
			else:
				grad = np.array( potential.getGradient(xVec, 0.05) )
		
		gradMag = np.sqrt( sum( [item**2 for item in grad] ) )
		grad /= gradMag
		# print 'Position:', xVec
		# print 'Gradient:', grad

		'''
		if (zLast > 0 and grad[2] < 0) or (zLast < 0 and grad[2] > 0):
			print 'zLast triggered'
			# grad = np.array( [0., 0., -1.] )
			if xVec[2] < 0:
				xVec[2] = -x_CAD
			else:
				xVec[2] = x_CAD
			posHistory.append( list(xVec) )
		'''

		dq[0] = xVec[2] #grad[2]
		meandq = sum(dq)/len(dq)
		# print dq, meandq, xVec[2]

		'''
		if bottom:
			# checkWireCollision(xVec)
			xVec += grad * t
			posHistory.append( list(xVec) )
		'''
		
		# Hit V- or U-wire
		if (abs( round(meandq, 5) ) == abs( round(xVec[2], 5) )) or (abs(xVec[2]) > limit): # or (zLast > 0 and grad[2] < 0) or (zLast < 0 and grad[2] > 0):
			# print 'Hit wire'

			# print 'zLast triggered'
			#if xVec[2] < 0:
			#	xVec[2] = -x_CAD
			#else:
			#	xVec[2] = x_CAD
			# xVec[2] = -x_CAD

			# Check if stuck in V-Wire
			if abs( xVec[2] ) < CATHODE_ANODE_y_DISTANCE + WIRE_DIAMETER:
				xVec[2] -= WIRE_DIAMETER
				posHistory.append( list(xVec) )
				print '\n=============\n'
				dq = deque([0]*10)

			# Hit U-wire!
			else:
				posHistory.append( list(xVec) )
				# print 'break'
				break

		else:
			dq.rotate(1)
			zLast = grad[2]
			xVec += grad * t # * v_DRIFT
			posHistory.append( list(xVec) )
	else:
		# Particle crossed APD plane -> kill it!
		# or reached end of Efield (z = 0.2), so 
		# xVec turns to [NaN] * 3

		if not np.isnan( np.array( xVec ) ).any():
			posHistory = False
			print
			print xVec
			print 'APD hit!'
		# else:
			# print 'End of Efield hit!'
			# Do nothing, use event

	# Quick and dirty fix for hitting a wire or 
	# reaching the end of the Efield
	if posHistory:
		if np.isnan( np.array(posHistory[-1]) ).any():
			posHistory.pop()

	if posHistory and (upFlag or leftFlag):
		posHistory = np.array( posHistory )
		if upFlag:
			posHistory[:,2] *= -1
		if leftFlag:
			posHistory[:,0] *= -1
		posHistory = posHistory.tolist()

	if not saveHistory and posHistory:
		posHistory = [posHistory[0], posHistory[-1]]

	# print 'Have a Kit Kat'
	if out_q:
		outdict = {}
		outdict[eventNum] = posHistory
		out_q.put(outdict)
	else:
		return posHistory

# === DRIFTELECTRON2D ===
# Using vector-box intersection
def driftElectron2d(startPos, potential=None, fEfieldData=None, PLANE_BINS=None, eventNum=0, out_q=None, saveHistory=False, verbose=False):
	x_CAD = 0.19223657
	r_RID = 0.1832356
	limit = 0.187

	# Bin settings for r and z
	bins = [ PLANE_BINS[0], PLANE_BINS[3] ]
	binMin = [ PLANE_BINS[1], PLANE_BINS[4] ]
	binMax = [ PLANE_BINS[2], PLANE_BINS[5] ]
	binWidth = [ (binMax[i] - binMin[i])/bins[i] for i in range(len(bins)) ]

        if verbose:
            print 'Starting drift...'
            print bins, binMin, binMax, binWidth

	posHistory = [startPos]
	xVec = list( startPos )
	posHistory = [ xVec ]

	upFlag = False
	if xVec[2] > 0:
		upFlag = True
		xVec[2] *= -1

	x, y, z = xVec
	theta = np.arctan2(y, x)
	r = np.sqrt( x**2 + y**2 )
	xVecPol = np.array( [r, z] )

        # Get reflector radius for angle theta
        reflectorRad2 = 0.04
        reflectorX2 = 0.215
        # r_RID = refRad(theta, r_RID, reflectorRad2, reflectorX2)

        # Height of the event in relation to limit
        lengthHeight = limit - abs(z)

        # Origin square for xVecPol
	sqOrigin = np.array( [roundBase(item, 10, binWidth[i]) for i, item in enumerate(xVecPol) ] )
        # Init last collision flags
	lastColl = [(False, False)] * 2

        length = 0
	while abs( xVecPol[1] ) < limit:
                if xVecPol[0] > r_RID:
                    posHistory = False
                    print 'PTFE Hit!'
                    break

		# getEvec needs cartesian bins, so a dummy xRange is used
                if potential:
                    Evec = np.array( potential.getNearestGradient([0., xVecPol[0], xVecPol[1]])[-2:] )
                    
                elif fEfieldData:
                    Evec = np.array( getEvec(fEfieldData, [0., sqOrigin[0], sqOrigin[1]], [1,0,1] + PLANE_BINS)[-2:] )
                    # Evec = getEvec(fEfieldData, [0., xVecPol[0], xVecPol[1]], [1,0,1] + PLANE_BINS)[-2:]

                if np.isnan( Evec ).any():
                    posHistory = False
                    break
                
                if verbose:
                    print 'Evec:', Evec
                    print 'xVecPol:', xVecPol
		
                xVecPolOld = list( xVecPol )
		xVecSq = xVecPol - sqOrigin

		newPosSq, coll = vectorBoxIntersection(xVecSq, Evec, binWidth)
		xVecPol = newPosSq + sqOrigin

                if saveHistory:
                    xVec = list( polToCart(xVecPol[0], xVecPol[1], theta) )
                    posHistory.append( xVec )

		coll = [ (coll[2*i], coll[2*i+1]) for i in range(2) ]
		checkColl = np.array(lastColl) % np.array(coll)

                if verbose:
                    print 'sqOrigin:', sqOrigin
                    print 'xVecSq:', xVecSq
                    print 'newPosSq:', newPosSq
                    print 'coll:', coll
                    print 'checkColl:', checkColl
                    print

		for i, item in enumerate(coll):
			if item[0] and (Evec[i] < 0):
				sqOrigin[i] -= binWidth[i]
			if item[1] and (Evec[i] > 0):
				sqOrigin[i] += binWidth[i]

		lastColl = coll
                length += np.linalg.norm(xVecPol - xVecPolOld)
        else:
                if abs( xVecPol[1] ) > limit:
                    t = 1
                    deltaFrac = (abs(xVecPol[1]) - limit) / ( Evec[1] * t )
                    xVecPol += deltaFrac * t * np.array( Evec )
                    length -= np.linalg.norm( deltaFrac * t * np.array( Evec ) )

                posHistory.append( list( polToCart( xVecPol[0], xVecPol[1], theta ) ) )

	if posHistory and upFlag:
		posHistory = np.array( posHistory )
		posHistory[:,2] *= -1
		posHistory = posHistory.tolist()

	if out_q:
		outdict = {}
		outdict[eventNum] = posHistory
		out_q.put(outdict)
	else:
		return posHistory, abs(length - lengthHeight)
	
def driftElectron(startPos, fEfieldData, H_BINS, eventNum=0, out_q=None):
	v_DRIFT = 0.00171e6		# m/s
	t_STHB = 0.05e-6		# s
	x_CAD = 19.223657e-2		# m
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

# == completeDrift ==
# Top part of the lower detector half is assumed to be radial 
# symmetric, so a 2d drift is done. Once a certain point is 
# reached, the 3d simulation near the V-lines is used to get the
# correct final event position
def completeDrift(startPos, potentialTop, efieldBottom, H_BINS, eventNum=0, out_q=None, saveHistory=False):
	# At top, do 2D-drift
	# topPosHistory, lengthCorrection = driftElectronPotential2d(startPos, potentialTop, None, eventNum, None, False, True, saveHistory)
	topPosHistory = af.artificialDrift(startPos, eventNum, None, False, saveHistory)

        if not topPosHistory:
            posHistory = False
        else:
            posHistory = topPosHistory

        '''
	if topPosHistory:
		endTop = topPosHistory[-1]
		if endTop:
			# At bottom, do 3D-drift
			bottomPosHistory = driftElectronPotential(endTop, efieldBottom, H_BINS, eventNum, None, True, True, saveHistory)
			
			if not bottomPosHistory:
				posHistory = topPosHistory # False
			else:
				if saveHistory:
					posHistory = topPosHistory[:-1] + bottomPosHistory[1:]
				else:
					posHistory = [topPosHistory[0], bottomPosHistory[-1]]
				#lastEntry = posHistory[-1]
				#lastEntry[2] = posHistory[0][2]
				#posHistory.append( lastEntry )

		else:
			posHistory = False
	else:
		posHistory = False
        '''

	# print 'eventNum:', eventNum
	# print 'posHistory:', posHistory
	# print

	outdict = {}
	outdict[eventNum] = posHistory
	out_q.put(outdict)

def uvPositionCorrection(u, v):
	UPPER_BOUND = 4.5
	LOWER_BOUND = 1.5

	return (LOWER_BOUND + u) % (UPPER_BOUND - LOWER_BOUND) + LOWER_BOUND, (LOWER_BOUND + v) % (UPPER_BOUND - LOWER_BOUND) + LOWER_BOUND

def checkWireCollision(xVec):
	x, y, z = xVec
	return

def getEField(fEfieldData, pos, H_BINS):
	bins, binMin, binMax, binWidth = getBinData(H_BINS)

	for i in range(3):
		if (pos[i] < binMin[i]) or (pos[i] > binMax[i]):
			return None

	#EFieldX = fEfield.Get('Ex')
	#EFieldY = fEfield.Get('Ey')
	#EFieldZ = fEfield.Get('Ez')
	#EFieldList = [EFieldX, EFieldY, EFieldZ]

	pos = np.array( [ roundBase(item, prec=10, base=binWidth[i]) for i, item in enumerate(pos) ] )
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
	print posBin, Evec
	return Ex, Ey, Ez

def getEvec(fEfieldData, xVec, H_BINS):
	# Get the Efield vector at the current particle position
	Evec = getEField(fEfieldData, xVec, H_BINS)
	
	if Evec:
		Evec = list( Evec )
	else:
		print xVec
		return None

	'''
	if not Evec[0] and not Evec[1] and not Evec[2]:
		if xVec[2] >= 0.:
			Evec[2] = 1.
		else:
			Evec[2] = -1.

		return Evec
	'''

	# Get Magnitude of Efield vector and normalize the vector
	Emag = np.sqrt( np.sum([item**2 for item in Evec]) )
	Evec = [ item/Emag for item in Evec ]
	return Evec

def getEvecMike(fEfieldData, xVec, H_BINS):
	bins, binMin, binMax, binWidth = getBinData(H_BINS)

	indexList = []
	for i, item in enumerate( xVec ):
		indexList.append( getIndex(item, binMin[i], binWidth[i] ) )
	field_index = indexList[0]*bins[1]*bins[2] + indexList[1]*bins[2] + indexList[2]
	# print 'field_index', field_index
	
	Evec = []
	for item in fEfieldData:
		Evec.append( item.GetBinContent( field_index ) )
	# print 'Evec', Evec

	Emag = np.sqrt( np.sum([item**2 for item in Evec]) )
	Evec = [ item/Emag for item in Evec ]
	return Evec

def getIndex(k, kMin, Dk):
	k -= kMin
	k /= Dk
	k += 0.5
	return int(k)

# ==== vectorBoxIntersection ====
# Get the coordinates of the intersection of vector vec
# located at position pos, defining a line in space, 
# and the bin boundaries binW. 
# Returns also which collision took place in the format
# 3D: coll = [Left, Right, Front, Back, Bottom, Top]
# 2D: coll = [Left, Right, Bottom, Top]
def vectorBoxIntersection(pos, vec, binW):
	pos = np.array(pos)
	vec = np.array(vec)
	binW = np.array(binW)

	tVec0 = []
	tVecW = []
	for i in range( len(pos) ):
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

def saveParticleTrace(fName, posHistory, N=None):
	f = open(fName, 'wb')

	if not N:
		N = len(posHistory)

	for i in range(N):
		item = np.array( posHistory[i] )
		x = np.array( item[:,0] )
		y = np.array( item[:,1] )
		z = np.array( item[:,2] )

		for j in range( len(x) ):
			f.write('%f\t%f\t%f\n' % (x[j], y[j], z[j]))
		f.write('\n\n')
	
	f.close()

def plotParticleTrace(posHistory, N=None):
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
	print hDict

	hDictMC = {}
	hDictMCDrift = {}
	for key in hDict:
		nameList = []
		if ' ' in key:
			entry = key.replace('(', ' ')
			for elem in entry.split(' '):
				try:
					nameList.append( float(elem.split(',')[0]) )
				except:
					pass

			print nameList
			sigma0, theta, z = nameList
		else:
			if 'Cut' in key and not 'Hist' in key:
				if 'Drift' in key:
					hDictMCDrift[tuple(nameList)] = hDict[key]
				else:
					hDictMC[tuple(nameList)] = hDict[key]

	# ==== TRUE DATA ====
	fTrue = ROOT.TFile.Open(fNameTrue)
	hDict = openAllHistos(fTrue)
	hDataTrue = hDict['DataSS']
	hMcTrue = hDict['McSS']
		
	# hDataTrue.Scale(1./hDataTrue.Integral(scaleBin, 20))
	# hMcTrue.Scale(1./hMcTrue.Integral(scaleBin, 20))

	hDataTrue.Scale(1./hDataTrue.Integral())
	hMcTrue.Scale(1./hMcTrue.Integral())

	hDiffTrue = getHistoDiffInv(hDataTrue, hMcTrue, 'difference (true data)')
	hDiffTrue.SetBins(20, 0, 0.2)

	# c2 = ROOT.TCanvas()
	# hDiffTrue.Draw('HIST')
	hMcTrue.SetTitle('True Data')
	pf.plotCompTwo(hMcTrue, hDataTrue, 'Standoff Distance [m]', 'Events / ( 0.1 m )', 'Regular', 'Drifted', True)

	# ==== COMPARISON ====
	for key in hDictMC:
		hMC = hDictMC[key]
		hMC.SetTitle('Random Data')
		hMCDrift = hDictMCDrift[key]

		# hMC.Scale(1./hMC.Integral(scaleBin, 20))
		# hMCDrift.Scale(1./hMCDrift.Integral(scaleBin, 20))

		hMC.Scale(1./hMC.Integral())
		hMCDrift.Scale(1./hMCDrift.Integral())

		hDiff = getHistoDiffInv(hMCDrift, hMC, 'difference')

		hDiffDiff = getHistoDiffInv(hDiff, hDiffTrue, 'hDiffDiff')

		rmsDiff = getRMSFromHist(hDiffDiff, 5)
		
		print np.square( np.array([hDiffDiff.GetBinContent(i) for i in range(1, 5 + 1)]) )
		print rmsDiff
		print

		# if rmsDiff < 10:
		pf.plotCompTwo(hMC, hMCDrift, 'Standoff Distance [m]', 'Events / ( 0.1 m )', 'Regular', 'Drifted', True)

		#c = ROOT.TCanvas()
		#hDiffDiff.Draw('HIST')
		#hDiffTrue.SetMarkerStyle(3)
		#hDiffTrue.Draw('same')
		#raw_input('end')
		#del c

		# hMC.Scale(10000)
		# hMCDrift.Scale(10000)
		pf.plotCompQQ(hMC, hMCDrift, [20, 0, 0.2], 'MC', 'MC drifted', True)

	fRand.Close()
	fTrue.Close()

def getRMSFromHist(hist, binN):
	rms = np.sqrt( sum( np.square( np.array([hist.GetBinContent(i) for i in range(1, binN + 1)]) ) ) )
	return rms

# ==== SUPPORT FUNCTIONS ====
def xy_to_uv(x, y):
	u = 0.5 * (x + np.sqrt(3.0) * y)
	v = 0.5 * (-x + np.sqrt(3.0) * y)
	return u, v

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

def hist3dToH_BINS(hist):
	axisList = [hist.GetXaxis(), hist.GetYaxis(), hist.GetZaxis()]
	
	H_BINS = []
	for axis in axisList:
		H_BINS += [int( axis.GetNbins() ), float( axis.GetXmin() ), float( axis.GetXmax() )]

	return H_BINS

# ==== TEST RANGE ====
def randomTest(H_BINS, TwoDFlag=False):
	if TwoDFlag:
		# get rid of x-range
		H_BINS = [1, 0, 1, 2000, -0.22724, 0.22724, 1000, -0.21788400000000002, 0.014380000000000004] # ctv.storeEField('./efield_data/', 'fieldExport_converted.csv', None)[3:]
		eFieldFn = './efield_data/efield_hist2d.root'
	else:
		eFieldFn = 'efield_hist.root'

	startPos = 0.16 # REFLECTORINNERRAD/1000. - 0.1
	endPos = REFLECTORINNERRAD/1000. - 0.0005
	height = -0.003
	N = 15

	driftTest(H_BINS, startPos, endPos, height, N, eFieldFn, TwoDFlag)

# def randomTest(H_BINS, startPos, endPos, height, N=20, TwoDFlag=False):
def driftTest(H_BINS, startPos, endPos, height, N, eFieldFn, TwoDFlag):
	# fEfield currently only used for Raycasting!
	fEfield = ROOT.TFile(eFieldFn, 'READ')
	EFieldData = [ fEfield.Get('Ex'), fEfield.Get('Ey'), fEfield.Get('Ez') ]

	# --- RBF INTERPOLATION ---

	# 2D interpolation
	if TwoDFlag:
		# Created by csv_to_vec
		posList = cPickle.load(open('efield_data/posList.p', 'rb'))
		vecList = cPickle.load(open('efield_data/vecList.p', 'rb'))
	else:
		# 3D interpolation data from vector field and potential
		# pot = vtp.vtkToPotential('exo200_wo_edges.vtk')

		posList = cPickle.load(open('posList.p', 'rb'))
		vecList = cPickle.load(open('vecList.p', 'rb'))
		# rbf = vtv.localRBF(*vtv.vtkToVec('exo200_wo_edges.vtk'))
	
	rbf = vtv.localRBF(posList, vecList)

        start = time.clock()
	for height in np.linspace(-0.01, -0.128, 5):
		f = open('gp_data/gradient2d_%d.dat' % (height*1000), 'w')

		for pos in np.linspace(startPos, endPos, N):
                        # RK4 with interpolation
			# posHistoryPot, lengthCorrection = driftElectronPotential2d([pos, 0, height], rbf, H_BINS, 0, None, False, True, True, True) 
                        posHistoryPot, lengthCorrection = driftElectron2d([pos, 0, height], rbf, None, H_BINS, 0, None, True)

                        # Old: uses vtk potential
                        # vtp.driftParticle(rbf, startPos, endPos, 'gp_data/gradient2d_%d.dat' % (height*1000), height, N)
			if posHistoryPot:	
				# print 'pos: %f, height: %f, length: %f' % (pos, height, lengthCorrection)
				for pos in posHistoryPot:
					x, y, z = pos
					f.write('%f\t%f\t%f\n' % (x, y, z))
				
				f.write('\n\n')
		f.close()
        print '-- Test Execution time: %s s --' % str(start - time.clock())

	# --- RAYCASTING ---
	'''
	print H_BINS

	f = open('raycastTest.dat', 'w')
	posHistoryTotal = []
	for x in np.linspace(startPos, endPos, N):
		startPos = [0, x, height]
		print startPos

		if TwoDFlag:
			posHistory = driftElectron2d(startPos, EFieldData, H_BINS, eventNum=0, out_q=None)
		else:
			posHistory = driftElectron(startPos, EFieldData, H_BINS, eventNum=0, out_q=None)

		for item in posHistory:
			f.write('%f\t%f\t%f\n' % (item[0], item[1], item[2]) )
		f.write('\n')
		posHistoryTotal.append(posHistory)

	f.close()
	'''

# === DRIFT EVALUATION ===
# Drift is performed in the x=0 area.
# N: Number of starting points in one dimension.
#    In total there are NxN points evaluated.
# reps: Number of repititions to perform. Each
#    repitition is evaluated and their mean 
#    time is put out in the end.
def driftEvaluation(N=30, reps=1, art=False):
    # Reflector radius
    r_RID = 0.1832356
    # Limit to switch to 3d drift
    limit = 0.187
    # Bin settings of the used 2d electric field
    H_BINS = [1, 0, 1, 2000, -0.22724, 0.22724, 1000, -0.21788400000000002, 0.014380000000000004]

    # Currently, only 2d vector field is supported
    if not art:
        posList = cPickle.load( open('efield_data/posList.p', 'rb') )
        vecList = cPickle.load( open('efield_data/vecList.p', 'rb') )
    else:
        # Use artificial field
        import artificialField as af

        rRange = H_BINS[4:6] + [H_BINS[3]]
        zRange = H_BINS[7:9] + [H_BINS[6]]
	R, Z = np.linspace(*rRange), np.linspace(*zRange)

        RR = np.array( list(R) * len(Z) )
        ZZ = np.array( [inner for outer in [len(R)*[z] for z in list(Z)] for inner in outer] )

        X = np.zeros( len(RR) )

        posList = zip(X, RR, ZZ)

	result = [np.array(af.getEField(r, z)) for (x, r, z) in posList]
        vecList = []
        for res in result:
            r, z = res
            norm = np.linalg.norm( res )
            vecList.append( (0., r/norm, z/norm) )

    # Create radial basis function
    rbf = vtv.localRBF(posList, vecList)

    # Get the parameter spaces
    rSpace = np.linspace(0.17, r_RID-0.0005, N)
    zSpace = np.linspace(-0.05, -limit+0.005, N) 

    # Return new radius after drift using posHistory
    def getRadius(posHistory):
        if posHistory:
            x, rRes, zRes = posHistory[-1]
        else:
            # Hit the PTFE
            rRes = np.nan

        return rRes

    # Compare to radii
    def compRadius(r1, r2):
        # print (r1, r2)
        if np.isnan(r1):
            if np.isnan(r2):
                return 0.
            else:
                return 1
        elif np.isnan(r2):
            if np.isnan(r1):
                return 0.
            else:
                return -1
        else:
            return r1 - r2

    # Arrays containing deviations
    RK4IntpTotal, RK4Total, EulerTotal, RBTotal = np.zeros(N*N), np.zeros(N*N), np.zeros(N*N), np.zeros(N*N)
    devList = [RK4IntpTotal, RK4Total, EulerTotal, RBTotal]

    # Arrays containing time information
    RK4IntpTime, RK4Time, EulerTime, RBTime = np.zeros(N*N), np.zeros(N*N), np.zeros(N*N), np.zeros(N*N)
    timeList = [RK4IntpTime, RK4Time, EulerTime, RBTime]

    # Loop over repititions
    for rep in range(reps):
        RK4Intp, RK4, Euler, RB = [], [], [], []
        RK4IntpT, RK4T, EulerT, RBT = [], [], [], []

        # Loop over z
        for z in zSpace:
            # Loop over r
            for r in rSpace:
                print (r, z)
                pos = [0., r, z]

                # Evaluate the different drifts.
                # Start with the best, and calculate the deviation
                # from its result for the other ones

                if art:
                    rComp = af.artificialDrift(pos, eventNum=0, out_q=None, diffuse=False, saveHistory=False)
                    if not rComp:
                        rComp = np.nan
                    else:
                        rComp = rComp[0][1]
                else:
                    rComp = r
                # print rComp
                    
                # def driftElectronPotential2d(startPos, potential, H_BINS, eventNum=0, out_q=None, diffuse=False, saveHistory=False, useRK4=False, interpolate=False)
                # RK4 with RBF interpolation
                startTime = time.clock()
                rRK4Intp = getRadius( driftElectronPotential2d(pos, rbf, H_BINS, 0, None, False, True, True, False)[0] )
                RK4Intp.append(compRadius(rComp, rRK4Intp))
                RK4IntpT.append(time.clock() - startTime)

                # RK4
                startTime = time.clock()
                rRK4 = getRadius( driftElectronPotential2d(pos, rbf, H_BINS, 0, None, False, True, True, False)[0] )
                RK4.append(compRadius(rComp, rRK4))
                RK4T.append(time.clock() - startTime)
                # Euler
                startTime = time.clock()
                rEuler = getRadius( driftElectronPotential2d(pos, rbf, H_BINS, 0, None, False, True, False, False)[0] )
                Euler.append(compRadius(rComp, rEuler))
                EulerT.append(time.clock() - startTime)
                # Ray-Box intersection
                startTime = time.clock()
                rRB = getRadius( driftElectron2d(pos, rbf, None, H_BINS, 0, None, True)[0] )
                RB.append(compRadius(rComp, rRB))
                RBT.append(time.clock() - startTime)

        # Add deviations to total
        RK4IntpTotal += np.array( RK4Intp )
        RK4Total += np.array( RK4 )
        EulerTotal += np.array( Euler )
        RBTotal += np.array( RB )

        # Add times to total
        RK4IntpTime += np.array( RK4IntpT )
        RK4Time += np.array( RK4T )
        EulerTime += np.array( EulerT )
        RBTime += np.array( RBT )

    # Normalize deviations to number of repititions
    devList = [item * 1./reps for item in devList]
    # Same for times
    timeList = [item * 1./reps for item in timeList]

    titleList = ['RK4 interpolate', 'RK4', 'Euler', 'Ray-Box']
    for i in range( len(devList) ):
        plotDeviation(rSpace, zSpace, devList[i])
        print '%s: %.3f s' % (titleList[i], sum( timeList[i] )/(N**2))

def plotDeviation(r, z, dev):
    import matplotlib.pyplot as plt
    N = len(r)
    
    r = np.array( N * list(r) ).reshape((N, N))
    z = np.array( [N * [i] for i in list(z)] ).reshape((N, N))
    dev = dev.reshape((N, N))

    cmap = plt.get_cmap('RdBu')

    fig, ax = plt.subplots()
    im = ax.pcolormesh(r, z, dev, vmin=-0.005, vmax=0.005, cmap=cmap)
    # ax.scatter(r, z, c=dev)
    fig.colorbar(im, ax=ax)

    fig.show()
    raw_input('')

# === RUNGE-KUTTA ===
# For vec = (0, r, z) 
def RK4(potential, vec, h=0.01, interpolate=False):
    # print vec
    vec = np.array( vec )
    # h: step-size

    def getGrad(vec, h, interpolate):
        # print vec,
        if interpolate:
            grad = np.array( potential.getGradient(vec, h) )
        else:
            grad = np.array( potential.getNearestGradient(vec) )
            # print 'Grad', grad

        if np.isnan(grad).any():
            return 0, 0

        dr = grad[1]
        dz = grad[2]
        # print (dr, dz)
        return dr, dz

    def getSlope(dr, dz):
        try:
            # slope = dz/dr
            slope = dr/dz
            if False: # abs(slope) > 5:
                return 0
                # return 0.01 * np.sign( slope )
            else:
                return slope

            '''
            slope = dz / abs( dr )
            if abs(slope) > 1:
                return np.sign(dz)
            else:
                return 0.01 * slope
            return 0.01 * dz / abs( dr )
            '''
        except:
            return 0.

    # k1
    dr, dz = getGrad(vec, h, interpolate)
    k1 = getSlope(dr, dz)

    # k2
    vR = np.array([0, 0.5*h*k1, 0])
    vZ = np.array([0, 0, -0.5*h])

    dr2, dz2 =  getGrad(vec + vR + vZ, h, interpolate) 
    k2 = getSlope(dr2, dz2)

    # k3
    # vZ remains the same
    vR = np.array([0, 0.5*h*k2, 0])

    dr3, dz3 = getGrad(vec + vR + vZ, h, interpolate)
    k3 = getSlope(dr3, dz3)

    # k4
    if False: # (dr > 0 and (dr2 < 0 or dr3 <0)) or (dr < 0 and (dr2 > 0 or dr3 > 0)):
        vR = np.array([0, 0, 0])
        dr4, dz4 = 0, 0
    else:
        vR = np.array([0, h*k3, 0.])
        vZ = np.array([0, 0, -h])
        dr4, dz4 = getGrad(vec + vR + vZ, h, interpolate)

    # k4 = getSlope(dr4, dz4)

    # k = 1./6 * (k1 + 2*k2 + 2*k3 + k4)
    
    # Calculate components from slope
    # x = 1 / np.sqrt( k**2 + 1 )
    # y = k / np.sqrt( k**2 + 1 ) 
    x = 1./6 * (dr + 2*dr2 + 2*dr3 + dr4)
    y = 1./6 * (dz + 2*dz2 + 2*dz3 + dz4)

    # Normalization
    n = np.sqrt(x**2 + y**2)
    x /= n
    y /= n
    # print (x, y)
    # print

    # print k1, k2, k3, k4, k
    # print (x, y), (x*h, y*h)
    return np.array( [0., x, y] )

# === SUPPORT ===
def dataToGraph(posHistory, name='default', color=1):
	mg = ROOT.TMultiGraph()
	mg.SetName(name)
	l = [0] * len(posHistory)
	for i, item in enumerate(posHistory):
		l[i] = ROOT.TGraph(len(item), np.array( np.array(item)[:,0] ), np.array( np.array(item)[:,2] ))
		l[i].SetLineColor(color)
		mg.Add(l[i])

	return mg

def rootListToRoot(treeName, dirName, outName):
	fList = getFileList(dirName, '.root')
	ch = getChain(treeName, fList)

	f = ROOT.TFile.Open(outName, 'RECREATE')
	f.Write(ch)
	f.Close()

