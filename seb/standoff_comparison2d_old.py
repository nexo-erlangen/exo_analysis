def compStandoff2d(dataTree, mcTree, bins=100):
	print 'compStandoff2d'
	ROOT.gROOT.cd()
	
	mcTreeCut = mcTree.CopyTree( ps.getCut(calibCut=True, energyCut=True, type='ss', MC=True) )
	dataTreeCut = dataTree.CopyTree( ps.getCut(calibCut=True, energyCut=True, type='ss', MC=False) )

	posList = []
	valList = []
	print 'Beginning loop'

	z = np.linspace(-150, 150, bins - 1)
	theta = np.linspace(-1, 1, bins - 1)

	for i in range(len(z) - 1):
		for j in range(len(theta) - 1):
			zMin = z[i]
			zMax = z[i+1]

			thetaMin = theta[j]
			thetaMax = theta[j+1]

			print 'z in range %f to %f' % (zMin, zMax)
			print 'theta in range %f to %f' % (thetaMin, thetaMax)
			print

			mcHist = getCompStandoff2d(mcTreeCut, zMin, zMax, thetaMin, thetaMax, 'mcHist')
			dataHist = getCompStandoff2d(dataTreeCut, zMin, zMax, thetaMin, thetaMax, 'dataHist')

			h = dataHist.Clone()
			h.Add(mcHist, -1.)
			h.Divide(h, dataHist)

			x = zMin + float(zMax - zMin)/2 
			y = thetaMin + float(thetaMax - thetaMin)/2
			val = []	
			valErr = []		
			for k in range(1, 5 + 1):
				val.append( h.GetBinContent(k) )
				valErr.append( h.GetBinError(k) )
			val.append( sum(val) )

			posList.append( (x, y) )
			valList.append( val )

			del mcHist
			del dataHist
			del h

	for i in range(5):
		plotStandoffComp2d(posList, valList, i)

def plotStandoffComp2d(posLs, valLs, idx=2):
	posLs = np.array( posLs )
	valLs = np.array( valLs )

	#leg = ROOT.TLegend(.75, .8, .89, .89)
	c = ROOT.TCanvas()

	x = np.array(posLs[:,0])
	y = np.array(posLs[:,1])

	l = ROOT.TGraph2D(len(x), x, y, np.array(valLs[:,idx]))
	l.SetLineColor(1 + idx*2)
	l.Draw('surf1')

	#leg.AddEntry(l[i], 'Bin %d' % (i + 1), 'l')
	#leg.Draw()

	c.Update()

	c.SaveAs('plots/standoffComp2d_%d.pdf' % idx)
	raw_input('end')

def getCompStandoff2d(tree, zMin, zMax, thetaMin, thetaMax, name='default'):
	standoffHist = ROOT.TH1D(name, name, 20, 0, 200)
	
	# stupid solution to get rid of memory leak:
	# Create a temporary root file in which
	# the data is written. Deleting treeCut
	# doesn't work for some reason
	f = ROOT.TFile('temp.root', 'RECREATE')
	# ROOT.gROOT.cd()

	treeCut = tree.CopyTree( 'cluster_z > %f && cluster_z < %f' % (zMin, zMax) )
	for i in range(treeCut.GetEntries()):
		treeCut.GetEntry(i)
		es = treeCut.EventSummary

		posX = np.array(es.cluster_x)
		posY = np.array(es.cluster_y)
		posZ = np.array(es.cluster_z)

		tTheta = ps.GetMaxTheta(posX, posY)
		# tZ = ps.GetMaxZ(posZ)

		if (tTheta > thetaMin and tTheta < thetaMax):
			standoffHist.Fill(getStandoff(posX[0], posY[0], posZ[0]))

	try:
		standoffHist.Scale(1./standoffHist.Integral())
	except:
		pass

	del treeCut
	return standoffHist
