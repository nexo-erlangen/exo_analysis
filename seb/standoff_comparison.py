import numpy as np
import plot_support as ps

import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')

# ==== CONSTANTS ====
CATHODE_ANODE_y_DISTANCE = 192.23657 # mm
REFLECTORINNERRAD = 183.2356 # mm

def compStandoff(dataTree, mcTree, comp='Z', bins=10):
	print 'compStandoff'
	ROOT.gROOT.cd()
	mcTreeCut = mcTree.CopyTree( ps.getCut(calibCut=True, energyCut=True, type='ss', MC=True) )
	dataTreeCut = dataTree.CopyTree( ps.getCut(calibCut=True, energyCut=True, type='ss', MC=False) )

	x_ls = []
	y_ls = []
	print 'Beginning loop'

	if comp == 'Z':
		a = np.linspace(-150, 150, bins - 1)
	else:
		# Theta case
		a = np.linspace(-1, 1, bins - 1)

	for i in range(len(a) - 1):
		Min = a[i]
		Max = a[i+1]

		if comp == 'Z':
			print 'zRange: %f to %f' % (Min, Max)
			mcHist = getCompStandoffZ(mcTreeCut, Min, Max, 'mcHist')
			dataHist = getCompStandoffZ(dataTreeCut, Min, Max, 'dataHist')
		else:
			print 'thetaRange: %f to %f' % (Min, Max)
			mcHist = getCompStandoffTheta(mcTreeCut, Min, Max, 'mcHist')
			dataHist = getCompStandoffTheta(dataTreeCut, Min, Max, 'dataHist')

		h = dataHist.Clone()
		h.Add(mcHist, -1.)
		h.Divide(h, dataHist)
		# res = np.array([0] * int(mcHist.GetEntries()) )
		# dataHist.Chi2Test(mcHist, 'UU NORM P', res)

		x = []
		y = []
		yErr = []
		# h.GetEntries()

		# x.append( Min + float(Max - Min)/2 )
		x = Min + float(Max - Min)/2
		for j in range(1, 5 + 1):
			# x.append( Min + float(Max - Min)/2 )
			# y.append( res[j-1] )
			y.append( h.GetBinContent(j) )
			yErr.append( h.GetBinError(j) )
		y.append( sum(y) )
		print x
		print y

		x_ls.append(x)
		y_ls.append(y)

		'''
		c = ROOT.TCanvas()
		# h.Draw('HIST')
		# mcHist.SetMarkerStyle(4)
		mcHist.Draw('HIST')
		dataHist.SetMarkerStyle(3)
		dataHist.Draw('same P')
		raw_input('end')
		del c	
		'''

		del dataHist
		del mcHist
		del h
		
	plotStandoffComp(x_ls, y_ls)

def plotStandoffComp(x_ls, y_ls, num=6):
	x = np.array( x_ls ) # x_ls = np.array( x_ls )
	y_ls = np.array( y_ls )

	l = np.array([ROOT.TGraph()] * num)
	leg = ROOT.TLegend(.75, .8, .89, .89)
	c = ROOT.TCanvas()

	for i in range(num):
		# x = np.array(x_ls[:,i])
		y = np.array(y_ls[:,i])
		print x, len(x)
		print y, len(y)

		l[i] = ROOT.TGraph(len(x), x, y)
		l[i].SetLineColor(1 + i*2)
		l[i].GetYaxis().SetRangeUser(-1.3, 1.3)
		l[i].Draw()

		leg.AddEntry(l[i], 'Bin %d' % (i + 1), 'l')
		leg.Draw()

		c.Update()
		
	c.SaveAs('plots/standoffComp.pdf')
	raw_input('end')

def getCompStandoffZ(tree, Min, Max, name='default'):
	standoffHist = ROOT.TH1D(name, name, 20, 0, 200)

	temp = ROOT.TFile('temp.root', 'RECREATE')
	treeCut = tree.CopyTree( 'cluster_z > %f && cluster_z < %f' % (Min, Max) )
	for i in range(treeCut.GetEntries()):
		treeCut.GetEntry(i)
		es = treeCut.EventSummary

		posX = np.array(es.cluster_x)
		posY = np.array(es.cluster_y)
		posZ = np.array(es.cluster_z)		
		standoffHist.Fill(getStandoff(posX[0], posY[0], posZ[0]))

		# print posX[0], posY[0], posZ[0], getStandoff(posX[0], posY[0], posZ[0])
	
	#try:
	#	standoffHist.Scale(1./standoffHist.Integral())
	#except:
	#	pass

	del treeCut
	return standoffHist

def getCompStandoffTheta(tree, Min, Max, name='default'):
	temp = ROOT.TFile('temp.root', 'RECREATE')
	standoffHist = ROOT.TH1D(name, name, 20, 0, 200)

	for i in range(tree.GetEntries()):
		tree.GetEntry(i)
		es = tree.EventSummary

		posX = np.array(es.cluster_x)
		posY = np.array(es.cluster_y)
		posZ = np.array(es.cluster_z)

		t = ps.GetMaxTheta(posX, posY)

		if (t > Min and t < Max):
			standoffHist.Fill(getStandoff(posX[0], posY[0], posZ[0]))

	try:
		standoffHist.Scale(1./standoffHist.Integral())
	except:
		pass

	return standoffHist

def getStandoff(x, y, z):
	r = np.sqrt( x**2 + y**2 )

	sd_r = REFLECTORINNERRAD - r
	sd_z = CATHODE_ANODE_y_DISTANCE - abs(z)

	if sd_r > sd_z:
		return sd_z
	else:
		return sd_r

	return sd_r

