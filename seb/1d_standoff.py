import sys
import numpy as np

import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')

# important stuff
APOTHEM = 162
Z_MIN = 5
Z_MAX = 182

# set fiducial volume
ROOT.EXOFiducialVolume.SetUserHexCut(APOTHEM, Z_MIN, Z_MAX)

# histogram parameters
H_BINS = 20
H_MIN = 0
H_MAX = 200

# energy cut
eMin = 750
eMax = 3500

def main():
	if len(sys.argv):
		MC_FNAME = sys.argv[1]
		DATA_FNAME = sys.argv[2]
	else:
		print "Usage: MCFILE DATAFILE"

	mcFile = ROOT.TFile(MC_FNAME)
	mcTree = mcFile.Get('mcTree')
	dataFile = ROOT.TFile(DATA_FNAME)
	dataTree = dataFile.Get('dataTree')

	mcHistSS = fillHisto(mcTree, getCut(True, True, 'ss', True), 'MC SS')
	# mcHistSS = fillHisto(mcTree, getCut(True, True, 'ss', True), 'MC SS')
	# mcHistSS.SetMarkerStyle(3)
	dataHistSS = fillHisto(dataTree, getCut(True, True, 'ss', False), 'Data SS')

	# Set label
	ROOT.gStyle.SetOptTitle(0)	# Remove title
	ROOT.gStyle.SetOptStat(0)	# Remove stats

	mcHistSS.GetXaxis().SetTitle('Standoff Distance [mm]')
	dataHistSS.GetXaxis().SetTitle('Standoff Distance [mm')

	# Set legend
	leg = ROOT.TLegend(.75, .8, .89, .89)
	leg.AddEntry(mcHistSS, 'MC', 'l')
	leg.AddEntry(dataHistSS, 'Data', 'p')

	# Set style
	mcHistSS.SetLineStyle(1)
	mcHistSS.SetMarkerStyle(1)
	mcHistSS.SetFillColorAlpha(30, 0.5)

	dataHistSS.SetMarkerStyle(3)

	# Draw
	#c = ROOT.TCanvas()

	#mcHistSS.Draw('H')
	#dataHistSS.Draw('Psame')
	#leg.Draw('same')

	# mcHistSS.Draw('Psame')
	#c.SaveAs('./plots/1d_standoff_ss.pdf')

	compHisto(mcHistSS, 'MC - SS', dataHistSS, 'Data - SS')
	#raw_input('end')

def fillHisto(Tree, cut, name):
	hist = ROOT.TH1D('hist', 'hist', H_BINS, H_MIN, H_MAX)
	cval = 'standoff_distance>>hist'

	Tree.Draw(cval, cut, 'goff')	# goff -> graphics off
	hist.SetName(name)
	hist.Sumw2()
	hist.Scale(1./hist.Integral())
	return hist

def getCut(calibCut=True, energyCut=True, type='ms', MC=True):
	eventSum = ROOT.EXOEventSummary()
	prepTree = ROOT.EXOPreprocessedTreeModule() 
	prepTree.SetIsData(not MC)
	prepTree.SetEventSummary(eventSum)

	fFOVs = prepTree.GetFOVs()
	# fFOVs.SetBooleanFlag('fIsDataOrMC', True)
	fFOVs.SetBooleanFlag('fVetoLikeCut', True)
	fFOVs.SetBooleanFlag('fCutMissingPosition', False)
	fFOVs.SetBooleanFlag('fEqualizeThreshold', False)
	fFOVs.SetBooleanFlag('fApplyVetoes', False)
	fFOVs.SetStringOption('fDiagonalCutDBFlavor', 'no-cut')

	if energyCut:
		if MC:
			enCut = " && (energy_mc > %d. && energy_mc < %d.) " % (eMin, eMax)
		else:
			enCut = " && (energy_%s > %d. && energy_%s < %d.) " % (type, eMin, type, eMax)
	else:
		enCut = ''

	if type == 'ms':
		typeCut = "&& multiplicity > 1.1 "
	else:
		typeCut = "&& abs(multiplicity - 1) < 0.1"

	return prepTree.GetDefaultCut() + enCut + typeCut

def compHisto(hist1, title1, hist2, title2):
	res = np.array([0]*H_BINS)
	p =	hist1.Chi2Test(hist2, 'UU NORM P', res)
	print p, res

	x = np.array([0 + float(i)/H_BINS * H_MAX for i in range(H_BINS)])
	# x = np.array( [1 + i for i in range(H_BINS)] )
	resgr = ROOT.TGraph(H_BINS, x, res)
	f = ROOT.TF1('f', 'TMath::Gaus(x,0,1)',-10,10)
	qqplot = ROOT.TGraphQQ(H_BINS, res, f)

	c = ROOT.TCanvas()
	c.Divide(2, 2)
	ROOT.gStyle.SetOptTitle(1)

	c.cd(1)
	hist1.SetLineStyle(1)
	hist1.SetMarkerStyle(1)
	hist1.SetFillColor(0)
	
	hist1.SetTitle(title1)
	hist1.GetXaxis().SetTitle('Standoff Distance [mm]')
	hist1.Draw('HIST')

	c.cd(2)
	hist2.SetTitle(title2)
	hist2.GetXaxis().SetTitle('Standoff Distance [mm')
	hist2.Draw('HIST')

	c.cd(3)
	resgr.GetXaxis().SetRangeUser(H_MIN, H_MAX)
	resgr.SetTitle('Normalized Residuals')
	resgr.Draw('APL')
	line = ROOT.TLine(0, 0, H_MAX, 0)
	line.SetLineStyle(ROOT.kDashed)
	line.SetLineWidth(2)
	line.Draw()

	c.cd(4)
	qqplot.SetMarkerStyle(3)
	qqplot.SetTitle('Q-Q Plot of Normalized Residuals')
	qqplot.Draw('AP')

	c.cd(0)
	c.Update()
	c.SaveAs('./plots/1d_standoff_ss_comp.pdf')
	raw_input('end')

if __name__ == '__main__':
	main()
