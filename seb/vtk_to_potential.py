try:
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
except:
    pass

import numpy as np
import random
from scipy import interpolate, spatial
from collections import deque
import time
# import matplotlib.pyplot as plt
# import matplotlib.cm as cm

import plot_support as ps
import ROOT

REFLECTORINNERRAD = 0.1832356
CATHODE_ANODE_y_DISTANCE = 0.19223657 # mm

def main():
	inFile = 'exo200_wo_edges.vtk'
	#inFile = '/home/hpc/capm/mppi025h/sfepy/exo200_wo_edges.vtk'
	# pot = vtkToPotential(inFile)
	# plotGradient(pot)
	# plotPotential(pot)
	# exportGradToDat(pot)
	# getEfieldHistogram(pot)
	# driftParticle(pot)

	# posArray, scalArray = vtkToVec(inFile)
	rbf = localRBF(*vtkToVec(inFile))
	# print rbf.getRBF([0, 0, -0.1], 0.01)
	# print rbf.getGradient([0.1, 0.1, 0.1], 0.1)

	# plotPotential(rbf)
	driftParticle(rbf, 0.162, 0.172, 10)

def driftParticle(potential, startPos, endPos, fName, height=-0.01, N=20):
	R = 0.1825

	lossCnt = 0

	posHistoryTotal = []
	t = 0.00001

	start = time.time()
	zLast = 0

	# for i in range(N):
	for i, x in enumerate( np.linspace(startPos, endPos, N) ):
		print 'Item %d' % i
		dq = deque([0]*5)
	
		xVec = [0, x, height]
		# xVec = [0.17, 0, -0.05] #getRandom()

		posHistory = [ list( xVec ) ]

		while abs(xVec[2]) < CATHODE_ANODE_y_DISTANCE:
			vecRsq = xVec[0]**2 + xVec[1]**2
			if vecRsq >= R**2:
				lossCnt += 1
				# break

			# grad = np.array( getGradient(xVec[0], xVec[1], xVec[2], potential) )
			grad = np.array( potential.getGradient(xVec, 0.05) )

			gradMag = np.sqrt( sum( [item**2 for item in grad] ) )
			grad /= gradMag

			if (zLast > 0 and grad[2] < 0) or (zLast < 0 and grad[2] > 0):
				grad = np.array( [0., 0., -1.] )

			dq[0] = grad[2]
			# print dq
			meandq = sum(dq)/len(dq)
			print grad[2], xVec[2]

			if abs( round(meandq, 7) ) == 1:
				if xVec[2] < 0:
					xVec[2] = -CATHODE_ANODE_y_DISTANCE
				else:
					xVec[2] = CATHODE_ANODE_y_DISTANCE

				posHistory.append( list(xVec) )
				# f.write('%f\t%f\t%f\n' % (xVec[0], xVec[1], xVec[2]))

				break

				#xVec = list(xVec)
				#xVec[2] = -0.18	
				#posHistory.append( xVec )
				#break
			else:
				dq.rotate(1)

				xVec += grad * t
				zLast = grad[2]
				posHistory.append( list(xVec) )

				# f.write('%f\t%f\t%f\n' % (xVec[0], xVec[1], xVec[2]))

		posHistoryTotal.append(posHistory)

	f = open(fName, 'w')
	for posHistory in posHistoryTotal:
		for item in posHistory:
			f.write('%f\t%f\t%f\n' % (item[0], item[1], item[2]))
		f.write('\n\n')
	f.close()
	
	end = time.time()
	print 'Elapsed time: %f s' % (end - start)
	print 'Mean time per event: %f s' % (float(end - start)/N)
	print 'Percent lost: %f' % round(float(lossCnt)/N, 2)

	return np.nan_to_num( posHistoryTotal )

def getEfieldHistogram(potential, H_BINS):
	print 'Generating electrical field...'

	N = H_BINS[0]
	# H_BINS = [N, -0.2, 0.2] * 3

	histPos = ROOT.TH3D('Epos', 'Epos', *H_BINS)

	histX = ROOT.TH3D('Ex', 'Ex', *H_BINS)
	histY = ROOT.TH3D('Ey', 'Ey', *H_BINS)
	histZ = ROOT.TH3D('Ez', 'Ez', *H_BINS)

	histList = [histX, histY, histZ]

	num = 0
	numTot = N**3
	for x in np.linspace(-0.2, 0.2, N):
		for y in np.linspace(-0.2, 0.2, N):
			for z in np.linspace(-0.2, 0.2, N):
				num += 1
				p = int(100.*(num+1)/numTot)
				ps.statusBar(p, 'screen')

				r = np.sqrt( x**2 + y**2 )
				if (r > REFLECTORINNERRAD) or (abs(z) > 0.2):
					continue

				vecDer = getGradient(x, y, z, potential)
				if np.isnan(vecDer[0]) or np.isnan(vecDer[1]) or np.isnan(vecDer[2]):
					continue

				histPos.Fill(*vecDer)

				for j, hist in enumerate(histList):
					hist.Fill(x, y, z, vecDer[j])

	print

	f = ROOT.TFile('/home/hpc/capm/mppi025h/sfepy/efield_hist.root', 'RECREATE')
	for hist in histList:
		hist.Write()
	f.Close()

def getRandom():
	r, z = 1, 1

	while (r > 0.18) and (abs(z) > 0.18):
		x = random.uniform(-0.2, 0.2)
		y = random.uniform(-0.2, 0.2)
		z = random.uniform(-0.2, 0.2)
		xVec = (x, y, z)
		r = np.sqrt( x**2 + y**2 )

	return x, y, z

def vtkToVec(inFile):
	print 'Reading data from %s...' % inFile,
	reader = vtk.vtkDataSetReader()
	reader.SetFileName(inFile)

	reader.ReadAllScalarsOn()
	reader.Update()

	data = reader.GetOutput()
	d = data.GetPointData()
	dScal = d.GetArray('t')

	posArray = []
	scalArray = []
	for i in range(data.GetNumberOfPoints()):
		posArray.append( data.GetPoint(i) )
		scalArray.append( dScal.GetValue(i) )

	scal = np.array( scalArray )
	print 'done'

	return posArray, scalArray

def vtkToPotential(inFile):
	posArray, scalArray = vtkToVec(inFile)

	print 'Interpolating potential...',
	potential = interpolate.LinearNDInterpolator(posArray, scalArray) # , rescale=True)
	# potential = interpolate.NearestNDInterpolator(posArray, scalArray)

	# Throws memory error
	# potential = interpolate.Rbf(x[:10000], y[:10000], z[:10000], scal[:10000], function='cubic')
	print 'done'

	return potential

def plotGradient(potential):
	x = np.linspace(0., 0.183, 100)
	y = 0
	z = np.linspace(-0.18, 0, 100)

	l = ROOT.TGraph2D( len(x)*len(z) )
	
	n = 0
	for i in x:
		for j in z:
			if (abs(j) > 0.18) or i > 0.18:
				continue
				
			grad = getGradient(i, y, j, potential)
			mag = np.sqrt( sum( [item**2 for item in grad] ) )

			l.SetPoint(n, i, j, grad[1]/mag)
			n += 1

	c = ROOT.TCanvas()
	l.Draw('colz')
	raw_input('end')

def exportGradToDat(potential):
	f = open('3dVector.dat', 'w')
	for i in np.linspace(0.13, 0.2, 10):
		for j in np.linspace(-0.1, 0.1, 10):
			for k in np.linspace(-0.1, 0, 10):
				grad = getGradient(i, j, k, potential)
				r = np.sqrt( i**2 + j**2 )
				if (abs(k) > 0.18) or r > 0.18:
					continue

				if np.isnan(grad[0]) and np.isnan(grad[1]) and np.isnan(grad[2]):
					continue

				grad = np.nan_to_num(grad)
				mag = np.sqrt( sum( [item**2 for item in grad] ) )

				f.write('%f\t%f\t%f\t%f\t%f\t%f\n' % (i, j, k, grad[0]/mag, grad[1]/mag, grad[2]/mag))

	f.close()

def plotPotential(potential):
	x = np.linspace(0., 0.01, 100)
	y = 0.
	z = np.linspace(-0.2, 0., 100)
	l = ROOT.TGraph2D(len(x) * len(z))

	n = 0
	val = [] # np.zeros((len(x), len(z)))

	nEntries = len(x) * len(z)
	print 'Interpolating entries...'
	for i, xi in enumerate(x):
		for j, zj in enumerate(z):
			rbfi = potential.getRBF([xi, y, zj], 1.)
			val.append( rbfi(xi, y, zj) )
			#val.append(potential(xi, y, zj))
			l.SetPoint(n, xi, zj, val[n])
			n += 1

			if not (n % int(nEntries/100.)):
				perc = int(100.*(n+1)/nEntries)
				ps.statusBar(perc, width='screen')
	print 
	print 'Done!'

	val = np.nan_to_num( np.array(val) )

	# x = np.unique(x)
	# z = np.unique(z)
	X, Z = np.meshgrid(x, z)
	VAL = val.reshape(len(x), len(z))

	plt.pcolormesh(X, Z, VAL)
	CS = plt.contour(X, Z, VAL, 30, linewidths=1, colors='k')
	plt.clabel(CS, inline=1, inline_spacing=3, fontsize=6, colors='k', use_clabeltext=1)
	# plt.colorbar()
	plt.show()

	'''
	startpos = REFLECTORINNERRAD - 0.01
	endPos = REFLECTORINNERRAD - 0.001

	posHistory = driftParticle(potential, startPos, endPos, 5)

	for item in posHistory:
		x = np.array(item)[:,0]
		z = np.array(item)[:,2]
		print x, z
		plt.plot(x, z)
	'''
	
	'''
	ROOT.gStyle.SetCanvasPreferGL(True)
	c = ROOT.TCanvas()
	l.Draw('colz')
	# c.SaveTo('plotPotential.pdf')
	raw_input('end')
	'''

def getGradient(x, y, z, potential, steps=100):
	xMin, xMax = (0., 0.19)
	yMin, yMax = (0., 0.19)
	zMin, zMax = (0., -0.2)

	DELTA = 0.01

	xRangeMin, xRangeMax = (x - DELTA if x - DELTA < xMin else xMin, x + DELTA if x + DELTA < xMax else xMax)
	yRangeMin, yRangeMax = (y - DELTA if y - DELTA < yMin else yMin, y + DELTA if y + DELTA < yMax else yMax)
	zRangeMin, zRangeMax = (z - DELTA if z - DELTA < zMin else zMin, z + DELTA if z + DELTA < zMax else zMax)

	spacingX = np.linspace(xRangeMin, xRangeMax, steps)
	spacingY = np.linspace(yRangeMin, yRangeMax, steps)
	spacingZ = np.linspace(zRangeMin, zRangeMax, steps)

	xVal = []
	yVal = []
	zVal = []

	for i in spacingX:
		xVal.append( potential( i, y, z ) )
	for i in spacingY:
		yVal.append( potential( x, i, z ) )
	for i in spacingZ:
		zVal.append( potential( x, y, i ) )

	splXDer = getDerivative(xVal, spacingX)
	splYDer = getDerivative(yVal, spacingY)
	splZDer = getDerivative(zVal, spacingZ)

	return splXDer(x), splYDer(y), splZDer(z)

	'''
	xDer = np.gradient(xVal)[len(spacingX)//2]
	yDer = np.gradient(yVal)[len(spacingY)//2]
	zDer = np.gradient(zVal)[len(spacingZ)//2]

	return xDer, yDer, zDer
	'''

def getDerivative(val, spacing):
	# set weights of nan values to zero
	# and also set the nan values themselves to zero
	w = np.isnan(val)
	val = np.nan_to_num(val)

	spl = interpolate.UnivariateSpline(spacing, val, w=~w, s=50000000)	# for 100: 100000

	splDer = spl.derivative()
	return splDer

# === RADIAL BASIS FUNCTION INTERPOLATION ===
class localRBF:
	def __init__(self, posArray, valArray, NN=500, intN=200):
		self.NN = NN		# number of nearest neighbors
		self.intN = intN	# number of interpolated points per axis
		self.R = 0.1832356
		self.zMax = 0.19223657

		self.x = np.array( np.array(posArray)[:, 0] )
		self.y = np.array( np.array(posArray)[:, 1] )
		self.z = np.array( np.array(posArray)[:, 2] )

		self.vals = np.array( valArray )

		xyz = np.c_[self.x, self.y, self.z]
		self.tree = spatial.cKDTree(xyz) # , leafsize=20)

	def getRBF(self, point, dist):
		point = np.array(point)
		distances, pIdx = self.tree.query(np.array( [point] ), k=self.NN, n_jobs=-1, distance_upper_bound=dist)
		distances = distances[:, 1:]
		self.distances = distances

		pIdxFilt = []
		for i, val in enumerate( distances[0] ):
			if val < dist and val != np.inf:
				pIdxFilt.append( pIdx[0][i] )

		xSelect = [ self.x[i] for i in pIdxFilt ]
		ySelect = [ self.y[i] for i in pIdxFilt ]
		zSelect = [ self.z[i] for i in pIdxFilt ]
		valSelect = [ self.vals[i] for i in pIdxFilt ]

		rbfi = interpolate.Rbf(xSelect, ySelect, zSelect, valSelect, function='multiquadric')
		return rbfi

	def getGradient(self, point, dist):
		rbfi = self.getRBF(point, dist)
		if dist > max( self.distances[0] ):
			dist = max( self.distances[0] )

		xInt, yInt, zInt = self.circleIntersection(point)
		xDist = list( self.correctBorders( xInt, (point[0] - dist, point[0] + dist) ) )
		yDist = self.correctBorders( yInt, (point[1] - dist, point[1] + dist) )
		zDist = list( self.correctBorders( zInt, (point[2] - dist, point[2] + dist) ) )

		# FIX! Because only half model is used
		if xDist[0] < 0:
			xDist[0] = 0
		#if zDist[1] > 0:
		#	zDist[1] = 0

		'''
		print xDist, xInt
		print yDist, yInt
		print zDist, zInt
		'''

		spacingX = list( np.linspace(xDist[0], point[0], self.intN//2) ) + list( np.linspace(point[0], xDist[1], self.intN//2)[1:] )
		spacingY = list( np.linspace(yDist[0], point[1], self.intN//2) ) + list( np.linspace(point[1], yDist[1], self.intN//2)[1:] )
		spacingZ = list( np.linspace(zDist[0], point[2], self.intN//2) ) + list( np.linspace(point[2], zDist[1], self.intN//2)[1:] )

		xVal = [ rbfi(x, point[1], point[2]) for x in spacingX ]
		yVal = [ rbfi(point[0], y, point[2]) for y in spacingY ]
		zVal = [ rbfi(point[0], point[1], z) for z in spacingZ ]

		# print xVal

		xDer = np.gradient(xVal)[self.intN//2 - 1]
		yDer = np.gradient(yVal)[self.intN//2 - 1]
		zDer = np.gradient(zVal)[self.intN//2 - 1]

		# print xDer
		# print

		return xDer, yDer, zDer

	def correctBorders(self, intPts, distLs):
		corrDist = list( distLs )
		if distLs[0] < intPts[0]:
			corrDist[0] = intPts[0]
		if distLs[1] > intPts[1]:
			corrDist[1] = intPts[1]
		return corrDist

	def solveCirc(self, c):
		x = np.sqrt(self.R**2 - c**2)
		return -x, x

	def circleIntersection(self, pos):
		# For x
		a, b = self.solveCirc(pos[1])	
		xInt = (a, b)
		# For y
		c, d = self.solveCirc(pos[0])
		yInt = (c, d)
		# For z
		zInt = (-self.zMax, self.zMax)

		return xInt, yInt, zInt

def createTree(posArray):
	# posArray, scalArray = vtkToVec(inFile)

	x = np.array( np.array(posArray)[:,0] )
	y = np.array( np.array(posArray)[:,1] )
	z = np.array( np.array(posArray)[:,2] )

	xyz = np.c_[x, y, z]
	tree = spatial.cKDTree(xyz)

	return tree, xyz

def getNeighbors(tree, xyz, point):
	point = np.array(point)

	distances, points = tree.query(point, k=2, n_jobs=-1)		# use multithreading
	# distances = distances[:, 1:]					# Remove k=1 distances
	# step = np.mean(distances)

	return distances, points

def getRBF(distances, points, dist):
	distancesFilt = []
	for val in distances:
		if val < dist:
			distancesFilt.append( val )

	#	interpolate.RBF(	

# === END - RBFI ===

def rk4(f, x0, y0, x1, n):
    vx = [0] * (n + 1)
    vy = [0] * (n + 1)
    h = (x1 - x0) / float(n)
    vx[0] = x = x0
    vy[0] = y = y0
    for i in range(1, n + 1):
        k1 = h * f(x, y)
        k2 = h * f(x + 0.5 * h, y + 0.5 * k1)
        k3 = h * f(x + 0.5 * h, y + 0.5 * k2)
        k4 = h * f(x + h, y + k3)
        vx[i] = x = x0 + i * h
        vy[i] = y = y + (k1 + k2 + k2 + k3 + k3 + k4) / 6
    return vx, vy

if __name__ == '__main__':
	main()

