import numpy as np

def main():
	pos = (0.1, 0.1, 0.9)
	vec = (0.1, 0.1, -0.8)
	binW = 10
	vectorBoxIntersection(pos, vec, binW)

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

	tVec0 = (0 - pos) / vec
	tVecW = (binW - pos) / vec
	
	tAll = np.concatenate( (tVec0, tVecW)) 
	tPlus = tAll[tAll > 0.]
	tMinus = tAll[tAll < 0.]

	fPosPlus = pos + vec * min(tPlus)
	fPosMinus = pos + vec * max(tMinus)
	# print fPosPlus
	# print fPosMinus

	# ++++ OLD ++++
	'''
	test = []
	for i in range(3):
		if not (np.isnan(tVec0[i]) and np.isnan(tVecW[i])):
			print 'test', tVec0[i], tVecW[i]
			print max(tVec0[i], tVecW[i])
			test.append(max(tVec0[i], tVecW[i]))
	'''

	# ++++ OLD AND NOT WORKING ++++
	'''
	if (not vec[0] and not vec[1]) or (not vec[2] and not vec[0]) or (not vec[1] and not vec[2]):
		try:
			fPos = pos + vec * np.nanmin(test) #np.nanmax(tAll) #max(tPlus)
		except:
			print np.nanmin(tAll)
			fPos = pos + vec * np.nanmin(tVec0)
	else:
		fPos = pos + vec * np.nanmin(test) #np.nanmin(tAll)
	'''

	'''
	try:
		fPos = pos + vec * min(tPlus)
		print '\t\tfPosPlus:', pos + vec * min(tPlus)
	except:
		try:
			fPos = pos + vec * max(tMinus)
			print '\t\tfPosMinus:', pos + vec * max(tMinus)
		except:
			pass 
	'''

	coll = []
	for i, bw in enumerate(binWidth):
		if fPosPlus[i] is 0:
			coll += [True, False]
		elif fPosPlus[i] is bw:
			coll += [False, True]
		else:
			coll += [False, False]

	return fPosPlus, coll

def driftElectron(startPos, fEfield, H_BINS):
	v_DRIFT = 0.00171e6		# m/s
	t_STHB = 0.05e-6		# s
	x_CAD = 19.223657e-2	# m
	x_APD = 20.44065e-2		# m

	bins, binMin, binMax, binWidth = getBinData(H_BINS)

	xVec = np.array(startPos)
	posHistory = [startPos]

	# Get the origin of the histogram bin for the first event
	cubeOrigin = np.array( [ roundBase(item, 5, binWidth[i]) for i, item in enumerate(xVec) ] )

	while(True):
		# Get the Efield vector at the current particle position
		Evec = list(getEField(fEfield, *xVec))
		if not Evec[0] and not Evec[1] and not Evec[2]:
			Evec[2] = -1

		# Get Magnitude of Efield vector and normalize the vector
		Emag = np.sqrt( np.sum([item**2 for item in Evec]) )
		Evec = [ item/Emag for item in Evec ]

		# Subtract the bin origin from the particle position
		# to transform it to coordinates within the bin
		xVec -= cubeOrigin
		newPosCub, coll = vectorBoxIntersection(xVecCub, Evec, binWidth)
		newPos = newPosCub + cubeOrigin

		# Check if drift completed
		if(abs(newPos[2]) > x_CAD):
			break

		# Append new position to history
		posHistory.append(newPos)

		# Check collision to shift bin origin
		coll = [ (coll[2*i], coll[2*i+1]) for i in range(3) ]	# split list in list of tuples
		
		# i - index for coordinates
		for i, item in enumerate(coll):
			if item[0]:
				cubeOrigin[i] -= binWidth[i]
			if item[1]:
				cubeOrigin[i] += binWidth[i]

	return posHistory

if __name__ == '__main__':
	main()

def appendWorker()
	out_q = Queue()
	procs = []

	p = multiprocessing.Process(target=func, args=args)
	procs.append(p)
	p.start()

def worker(func, args, out_q):
	outdict = {}

def randomProcess(hist, H_BINS, N):
	start = time.time()

	tree = ROOT.TTree('random', 'random')
	randHist = generateRandom(hist, tree, H_BINS, N)
	x, y, z = readTreeToCoord(tree)
	fEfield = ROOT.TFile('~/sfepy/efield_hist.root', 'READ')

	N = tree.GetEntries()
	print 'Drifting %d particles...' % N

	posHistory = []
	for i in range(N+1):
		# TODO
		# Better code needed, because statusBar is 
		# printed everytime and if the percentage
		# is incremented
		p = int(100.*i/N)
		statusBar(p, os.popen('tput cols', 'r').read())

		tree.GetEntry(i)
		xVec = [x[0], y[0], z[0]]
		history = driftElectron(xVec, fEfield, H_BINS)
		for item in history:
			logging.info(item)
			#print item
		if history:
			posHistory.append(history)

	print 
	end = time.time()
	print 'Elapsed time: %s s' % str((end - start)/60.)

	plotParticleTrace(posHistory, 10)

	fEfield.Close()

def driftElectronOld(startPos, fEfield):
	v_DRIFT = 0.00171e6		# m/s
	t_STHB = 0.05e-6		# s
	x_CAD = 19.223657e-2	# m
	x_APD = 20.44065e-2		# m
	d_W = 0.0127e-2			# m

	STEP = 100				# append only every STEPth step

	# TODO: Need function to get electric field vector at position x, y, z
	# Vectors only known for specific points in space -> interpolate for positions between those points?

	xVec = np.array(startPos)

	posHistory = [startPos]
	s = 0
	while(True):
		Evec = list(getEField(fEfield, *xVec))
		if not Evec[0] and not Evec[1] and not Evec[2]:
			Evec[2] = -1

		Emag = np.sqrt( np.sum([item**2 for item in Evec]) )
		Evec = [item/Emag for item in Evec]
		print Evec

		dVec = [ 1 * v_DRIFT * t_STHB * item for item in Evec ]

		if(abs(xVec[2] + dVec[2]) > x_CAD):
			d_frac = min( [ 1, (x_CAD - abs(xVec[2])/dVec[2]) ] )
			dVec = [item * d_frac for item in dVec]
			break

		xVec += dVec
		if not s % STEP:
			posHistory.append(list(xVec))

		s += 1

	return xVec, posHistory

def randomProcess(hist, H_BINS, N):
	start = time.time()

	tree = ROOT.TTree('random', 'random')
	randHist = generateRandom(hist, tree, H_BINS, N)
	x, y, z = readTreeToCoord(tree)
	fEfield = ROOT.TFile('~/sfepy/efield_hist.root', 'READ')


	N = tree.GetEntries()
	print 'Drifting %d particles...' % N

	posHistory = []
	for i in range(N+1):
		# TODO
		# Better code needed, because statusBar is 
		# printed everytime and if the percentage
		# is incremented
		p = int(100.*i/N)
		statusBar(p, os.popen('tput cols', 'r').read())

		tree.GetEntry(i)
		xVec = [x[0], y[0], z[0]]
		history = driftElectron(xVec, fEfield, H_BINS)
		for item in history:
			logging.info(item)
			#print item
		if history:
			posHistory.append(history)

	print 
	end = time.time()
	print 'Elapsed time: %s s' % str((end - start)/60.)

	plotParticleTrace(posHistory, 10)

	fEfield.Close()
