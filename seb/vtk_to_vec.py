try:
    import vtk
except:
    pass

import numpy as np
import os, sys
import cPickle
import ROOT
from scipy import interpolate, spatial
# from settingsQxu_cart_rand import *

import plot_support as ps
import vtk_to_potential as vtp
ROOT.PyConfig.IgnoreCommandLineOptions = True

np.seterr(divide='ignore', invalid='ignore')

def main():
	inFile = '/exo200_wo_edges.vtk'
	# storeEField('.', inFile, H_BINS_RAND, True)

	# storeEFieldToFile(inFile)

	# Read list from file
	posList = cPickle.load(open('posList.p', 'rb') )
	vecList = cPickle.load(open('vecList.p', 'rb') )

	rbf = localRBF(posList, vecList)
	vtp.driftParticle(rbf, 0.17, 0.182, -0.02, 100)

def storeEFieldToFile(inFile):
	# Write list to file
	posList, vecList = storeEFieldToList('.', inFile)
	cPickle.dump(posList, open('posList.p', 'wb'))
	cPickle.dump(vecList, open('vecList.p', 'wb'))

def storeEField(inDir, inFile, H_BINS, verbose=False):
	print 'Reading vtk file...',
	# Read vtk-file and get stored vertices and vectors
	reader = vtk.vtkDataSetReader()
	reader.SetFileName(inDir + '/' + inFile)

	reader.ReadAllScalarsOn()
	reader.ReadAllVectorsOn()
	reader.Update()

	data = reader.GetOutput()

	CellArray = data.GetCells()
	Cells = CellArray.GetData()
	Vectors = data.GetCellData().GetVectors()

	# Bin settings of cartesian coordinates
	# H_BINS = [50, -0.2, 0.2, 50, -0.2, 0.2, 50, -0.2, 0.2]
	H_BIN_W = [ (H_BINS[2+3*i] - H_BINS[1+3*i]) / H_BINS[0+3*i] for i in range(0, 3) ]

	# Number of entries at position
	histPos = ROOT.TH3D('Epos', 'Epos', *H_BINS)

	# EField magnitude
	histMag = ROOT.TH3D('EMag', 'EMag', *H_BINS)

	# EField vector
	# The vector itself is in cartesian coordinates!
	histX = ROOT.TH3D('Ex', 'Ex', *H_BINS)
	histY = ROOT.TH3D('Ey', 'Ey', *H_BINS)
	histZ = ROOT.TH3D('Ez', 'Ez', *H_BINS)

	histList = [histX, histY, histZ]

	# Get maximum vector component to normalize values
	# vecMax = max( map( max, zip( *[ Vectors.GetTuple(i) for i in range(CellArray.GetNumberOfCells()) ] ) ) )

	for i in range(CellArray.GetNumberOfCells()):
		# Get the vertex-ids for a tetrahedron
		cellVals = [ Cells.GetValue(k) for k in range(i*5, (i+1)*5) ][1:] 
		# Calculate the center of mass position of the tetrahedron
		cellPos = np.sum( [np.array( data.GetPoint(entry) ) for entry in cellVals], axis=0 ) / 4.
		# Get the efield vector at the tetrahedron position
		cellVec = np.array( Vectors.GetTuple(i) )
		# Convert cell position to cylindrical coordinates
		r, t, z = cellPos
		#print '(%f, %f, %f):' % (r, t, z),
		#print cellVec

		# Vector magnitude
		vecMag = np.sqrt( sum([ item**2 for item in cellVec ]) )
		histMag.Fill(r, t, z, vecMag)

		# Increment count at cellPos 
		histPos.Fill(*cellPos)

		#x, y, z = cellPosBinNum
		for j, hist in enumerate(histList):
			# Increase histogram value at cellPos by
			# j-component of cellVec
			hist.Fill(r, t, z, cellVec[j])

	for hist in histList:
		hist.Divide(histPos)
	histMag.Divide(histPos)

	f = ROOT.TFile(inDir + '/efield_hist.root', 'recreate')
	for hist in histList:
		hist.Write()

	print 'done'

	#c = ROOT.TCanvas()
	#histX.Project3D('zy').Draw('SURF1 CYL')
	#c1 = ROOT.TCanvas()
	#histX.Project3D('xy').Draw('SURF1 POL')
	
	if verbose:
		c2 = ROOT.TCanvas()
		histPos.Project3D('zx').Draw('colz')
		c3 = ROOT.TCanvas()
		histPos.Project3D('zy').Draw('colz')
		c4 = ROOT.TCanvas()
		histPos.Project3D('yx').Draw('colz')

		raw_input('end')

	f.Close()
	
def storeEFieldToList(inDir, inFile):
	print 'Reading vtk file...',
	# Read vtk-file and get stored vertices and vectors
	reader = vtk.vtkDataSetReader()
	reader.SetFileName(inDir + '/' + inFile)

	reader.ReadAllScalarsOn()
	reader.ReadAllVectorsOn()
	reader.Update()

	data = reader.GetOutput()

	CellArray = data.GetCells()
	Cells = CellArray.GetData()
	Vectors = data.GetCellData().GetVectors()

	print 'Writing %d entries to list...' % CellArray.GetNumberOfCells()

	posList = []
	vecList = []
	for i in range(CellArray.GetNumberOfCells()):
		if not (CellArray.GetNumberOfCells() % (i+1)):
			p = int(100.*(i+1)/CellArray.GetNumberOfCells())
			ps.statusBar(p, 'screen')
		
		cellVals = [ Cells.GetValue(k) for k in range(i*5, (i+1)*5) ][1:]
		cellPos = np.sum( [np.array( data.GetPoint(entry) ) for entry in cellVals], axis=0) / 4.
		cellVec = np.array( Vectors.GetTuple(i) )

		posList.append( cellPos )
		vecList.append( cellVec )

	print 'Done!'
	return posList, vecList

def roundBase(x, prec=5, base=.05):
	return round(base * round(float(x)/base), prec)

def cartToCyl(x, y, z):
	R = np.sqrt(x**2 + y**2)
	theta = np.arctan2(float(y), x)
	return R, theta, z

# === RADIAL BASIS FUNCTION INTERPOLATION ===
class localRBF:
	def __init__(self, posArray, vecArray, NN=9):
		self.NN = NN		# number of nearest neighbors
		self.R = 0.1832356
		self.zMax = 0.19223657

		self.x = np.array( np.array(posArray)[:, 0] )
		self.y = np.array( np.array(posArray)[:, 1] )
		self.z = np.array( np.array(posArray)[:, 2] )

		self.vecs = np.array( vecArray )

		xyz = np.c_[self.x, self.y, self.z]
		self.tree = spatial.cKDTree(xyz, balanced_tree=True) # , leafsize=20)

	def getRBF(self, point, dist):
                # import time
                # start = time.clock()
                
		point = np.array(point)
                # print point
		distances, pIdx = self.tree.query(np.array( [point] ), k=self.NN, n_jobs=1, distance_upper_bound=dist)
		# distances = distances[:, 1:]
		self.distances = distances
		
                # stop1 = time.clock()
                # print 'Init:', stop1-start

                pIdxFilt = [ idx for idx in pIdx[0] if idx != len(self.x) ]
                # print pIdxFilt
                '''
		pIdxFilt = []
		for i, val in enumerate( distances[0] ):
			if val < dist and val != np.inf:
				pIdxFilt.append( pIdx[0][i] )
                '''

		xSelect = np.nan_to_num( np.array( [ self.x[i] for i in pIdxFilt ] ) )
		ySelect = np.nan_to_num( np.array( [ self.y[i] for i in pIdxFilt ] ) )
		zSelect = np.nan_to_num( np.array( [ self.z[i] for i in pIdxFilt ] ) )

		xVecSelect = np.nan_to_num( np.array( [ self.vecs[i][0] for i in pIdxFilt ] ) )
		yVecSelect = np.nan_to_num( np.array( [ self.vecs[i][1] for i in pIdxFilt ] ) )
		zVecSelect = np.nan_to_num( np.array( [ self.vecs[i][2] for i in pIdxFilt ] ) )

		'''
		print
		print np.mean( xSelect ), np.mean( xVecSelect )
		print np.mean( ySelect ), np.mean( yVecSelect )
		print np.mean( zSelect ), np.mean( zVecSelect )
		print
		'''

		rbfiX = interpolate.Rbf(xSelect, ySelect, zSelect, xVecSelect, function='multiquadric')
		rbfiY = interpolate.Rbf(xSelect, ySelect, zSelect, yVecSelect, function='multiquadric')
		rbfiZ = interpolate.Rbf(xSelect, ySelect, zSelect, zVecSelect, function='multiquadric')

                # print 'End:', time.clock()-stop1
		return rbfiX, rbfiY, rbfiZ

	def getGradient(self, point, dist):
		#if point[2] > 0:
		#	return [0., 0., 1]

		#if point[0] < 0:
		#	if point[2] > 0:
		#		return [0., 0., 1]
		#	else:
		#		return [0., 0., -1]

		rbfiX, rbfiY, rbfiZ = self.getRBF(point, dist)
                '''
		if dist > max( self.distances[0] ):
			dist = max( self.distances[0] )
                '''

		vecX, vecY, vecZ = [ rbfiX(*point), rbfiY(*point), rbfiZ(*point) ]
		# print vecX, vecY, vecZ
		mag = np.sqrt( vecX**2 + vecY**2 + vecZ**2 )
		return [ vecX/mag, vecY/mag, vecZ/mag ]

	def getNearestGradient(self, point):
		distances, pIdx = self.tree.query(np.array( [point] ))
		vecX, vecY, vecZ = self.vecs[pIdx][0]
		mag = np.sqrt( vecX**2 + vecY**2 + vecZ**2 )

		return [ vecX/mag, vecY/mag, vecZ/mag ]

if __name__ == '__main__':
	main()

