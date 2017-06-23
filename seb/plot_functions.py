import numpy as np
from plot_support import *

import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')
ROOT.PyConfig.IgnoreCommandLineOptions = True

# plotHisto1d
# ===========
# Plot a TH1 histogram.

def plotHisto1d(h, xTitle, yTitle, title, show=True, out=''):
	setShow(show)
	c = ROOT.TCanvas()

	h.SetStats(ROOT.kFALSE)
	h.SetTitle(title)
	h.GetXaxis().SetTitle(xTitle)
	h.GetYaxis().SetTitle(yTitle)

	h.Draw('HIST')
	if(out):
		c.SaveAs(out)
	
	if(show):
		raw_input('end')

	c.Close()

# plotCompTwo
# ===========
# Plot two TH1 histograms h1 and h2 on top of each other.
# The titles 'title1' and 'title2' are used to represent
# the data in the legend. Titles 'xTitle' and 'yTitle' 
# are used to label the axes in the plot.

def plotCompTwo(h1, h2, xTitle, yTitle, title1, title2, show=True, out=''):
	setShow(show)
	c = ROOT.TCanvas()
	#ROOT.gStyle.SetTitle(ROOT.kFALSE)

	# Set legend
	leg = ROOT.TLegend(.75, .8, .89, .89)
	leg.AddEntry(h1, title1, 'l')
	leg.AddEntry(h2, title2, 'p')
	
	h1.SetStats(ROOT.kFALSE)
	#h1.SetTitle(title1)
	h1.GetXaxis().SetTitle(xTitle)
	h1.GetYaxis().SetTitle(yTitle)

	h2.SetStats(ROOT.kFALSE)
	#h2.SetTitle(title2)
	h2.GetXaxis().SetTitle(xTitle)
	h2.GetYaxis().SetTitle(yTitle)

	h1.Draw('HIST')
	h2.SetMarkerStyle(3)
	h2.Draw('Psame')
	leg.Draw('same')

	if(out):
		c.SaveAs(out)

	if(show):
		raw_input('end')

	c.Close()

# plotCompQQ
# ==========
# Plot the specified TH1 histograms hist1 and hist2
# on a 2x2 canvas together with the normalized residues
# and their QQ-plot

def plotCompQQ(hist1, hist2, hBins, title1, title2, show=True, out=''):
	setShow(show)
	H_BINS = hBins[0]
	H_MIN = hBins[1]
	H_MAX = hBins[2]

	res = np.array([0]*H_BINS)
	p = hist1.Chi2Test(hist2, 'UU NORM P', res)
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
	if(out):
		c.SaveAs(out)

	if(show):
		raw_input('end')

	c.Close()

# plotHisto2d
# ===========
# Plot a TH2 histogram.

def plotHisto2d(h, xTitle, yTitle, title, show=True, out=''):
	setShow(show)
	c = ROOT.TCanvas()

	h.SetStats(ROOT.kFALSE)
	h.SetTitle(title)
	h.GetXaxis().SetTitle(xTitle)
	h.GetYaxis().SetTitle(yTitle)

	h.Draw('colz')
	if(out):
		c.SaveAs(out)

	if(show):
		raw_input('end')

	c.Close()

# plotHisto3d
# ===========
# Plot a TH3 histogram h in a 3d-plot. The specified 
# function 'transFunc' is used to set the opacity of
# each bin according to their value.
# Uses GL to display canvas!

def plotHisto3d(h, xTitle, yTitle, zTitle, title, transFunc, show=True, out=''):
	setShow(show)
	ROOT.gStyle.SetCanvasPreferGL(ROOT.kTRUE)

	c = ROOT.TCanvas()

	lf = h.GetListOfFunctions()
	if(lf):
		tf = ROOT.TF1('TransferFunction', transFunc, 0, h.GetMaximum(), 1)
		tf.SetParameters(np.array([h.GetMaximum()]))
		lf.Add(tf)

	h.SetStats(ROOT.kFALSE)
	h.SetTitle(title)
	h.GetXaxis().SetTitle(xTitle)
	h.GetYaxis().SetTitle(yTitle)
	h.GetZaxis().SetTitle(zTitle)

	h.Draw('glcol')
	if(out):
		c.SaveAs(out)

	if(show):
		raw_input('end')

	c.Close()
	ROOT.gStyle.SetCanvasPreferGL(ROOT.kFALSE)

# plotHistoProject3d
# ==================
# Creates the projection of a TH3 histogram h onto the
# with option specified axis or surface.
def plotHistoProject3d(h, xTitle, yTitle, title, option='xy', show=True, out=''):
	setShow(show)
	c = ROOT.TCanvas()

	k = h.Project3D(option)

	k.SetStats(ROOT.kFALSE)
	k.SetTitle(title)
	k.GetXaxis().SetTitle(xTitle)
	k.GetYaxis().SetTitle(yTitle)

	k.Draw('colz')
	if(out):
		c.SaveAs(out)

	if(show):
		raw_input('end')

	c.Close()

# plotHistoSliceMulti
# ===================
# If proj is set, the projections of a TH3 histogram hist
# on each axis are performed and plotted on the same canvas.
# If unset, the slices at the found maxima of the 
# coordinates are done.
# cyl is set true, if the axis are depicted in cylindrical 
# coordinates.

def plotHistoSliceMulti(hist, BIN_SET, proj=False, cyl=True, show=True, out=''):
	setShow(show)
	if(proj):
		k = hist.Project3D('yx').Clone()
		l = hist.Project3D('zx').Clone()
		m = hist.Project3D('yz').Clone()

	else:
		x, y, z = getMaximum(hist, BIN_SET)
		k = getSlice(hist, BIN_SET, 'z', z).Clone()
		l = getSlice(hist, BIN_SET, 'y', y).Clone()
		m = getSlice(hist, BIN_SET, 'x', x).Clone()

	# k.Scale(1./k.Integral())
	# l.Scale(1./l.Integral())
	# m.Scale(1./m.Integral())
	plotHistoSliceMultiExecute(k, l, m, cyl, show, out)

def plotHistoSliceMultiExecute(k, l, m, cyl=True, show=False, out=''):
	if(cyl):
		k.SetTitle('z vs. r')
		k.GetXaxis().SetTitle('r [mm]')
		k.GetYaxis().SetTitle('z [mm]')

		l.SetTitle('theta vs. r')
		l.GetXaxis().SetTitle('r [mm]')
		l.GetYaxis().SetTitle('theta [pi]')

		m.SetTitle('z vs. theta')
		m.GetXaxis().SetTitle('theta [pi]')
		m.GetYaxis().SetTitle('z [mm]')
	else:
		k.SetTitle('y vs. x')
		k.GetXaxis().SetTitle('x [mm]')
		k.GetYaxis().SetTitle('y [mm]')

		l.SetTitle('z vs. x')
		l.GetXaxis().SetTitle('x [mm]')
		l.GetYaxis().SetTitle('z [mm]')

		m.SetTitle('y vs. z')
		m.GetXaxis().SetTitle('z [mm]')
		m.GetYaxis().SetTitle('y [mm]')

	c = ROOT.TCanvas()
	c.Divide(2, 2)

	c.cd(1)
	k.SetStats(ROOT.kFALSE)
	k.Draw('colz')
	#k.Draw('')

	c.cd(2)
	l.SetStats(ROOT.kFALSE)
	l.Draw('colz')
	#l.Draw('')

	c.cd(3)
	m.SetStats(ROOT.kFALSE)
	m.Draw('colz')
	#m.Draw('')

	c.cd(0)
	c.Update()
	if(out):
		c.SaveAs(out)

	if(show):
		raw_input('end')

	c.Close()

# plotHistoSliceSingle
# ====================
# Gets a slice or projection of a TH3 histogram as 
# a TH2 where option specifies the axis which is held 
# constant. 
# If proj2d is set, a projection onto TH2 is done, else
# the slice is performed. Next to the histogram, the TH1
# slices or projections, depending if proj1d is set or 
# not, are shown.

def plotHistoSliceSingle(hist, HBINS, option='x', value=0, proj2d=False, proj1d=False, show=True, out='', zRange=None):
	setShow(show)
	c = ROOT.TCanvas()
	c.Divide(2, 2)

	if(proj2d):
		if(option=='x'):
			hist2d = hist.Project3D('yz')
		elif(option=='y'):
			hist2d = hist.Project3D('zx')
		elif(option=='z'):
			hist2d = hist.Project3D('yx')
		else:
			return False
	else:
		hist2d = getSlice(hist, HBINS, option, value)

        hist2d.Scale(1./hist2d.Integral())
        if zRange:
            hist2d.SetMinimum(zRange[0])
            hist2d.SetMaximum(zRange[1])

	if(proj1d):
		hist1dX = hist2d.ProjectionX()
		hist1dY = hist2d.ProjectionY()
	else:
		x, y = getMaximum2d(hist2d, HBINS[3:9])
		hist1dX = getSlice2d(hist2d, 'y', x)
		hist1dY = getSlice2d(hist2d, 'x', y)

	hist1dX.Scale(1./hist1dX.Integral())
	hist1dY.Scale(1./hist1dY.Integral())
	hist1dX.SetFillColor(38)
	hist1dY.SetFillColor(38)

	# 2d histogram
	c.cd(1)
	hist2d.SetStats(ROOT.kFALSE)
	hist2d.Draw('colz')

	# first projection
	c.cd(2)
	hist1dY.SetStats(ROOT.kFALSE)
	hist1dY.Draw('HBAR4 l')

	# second projection
	c.cd(3)
	hist1dX.SetStats(ROOT.kFALSE)
	hist1dX.Draw('BAR4 l')

	c.cd(0)
	c.Update()
	if(out):
		c.SaveAs(out)

	if(show):
		raw_input('end')

	c.Close()

# plot3d
# ======
# Gets the maximum value of TH3 histogram hist and performs
# 'plotHistoSlice' for each of the coordinates.

def plotMaximum3d(hist, BIN_SET, proj2d=False, proj1d=False, cyl=True, show=True, out=''):
	x, y, z = getMaximum(hist, BIN_SET)
	print 'Maximum for ' + ('cylindrical ' if cyl else 'cartesian ') + 'coordinates at: (%f, %f, %f)' % (x, y, z)

	out = out.split('.pdf')[0]
	plotHistoSliceSingle(hist, BIN_SET, 'x', x, proj2d, proj1d, show, out + ('_r.pdf' if cyl else '_x.pdf'))
	plotHistoSliceSingle(hist, BIN_SET, 'y', y, proj2d, proj1d, show, out + ('_z.pdf' if cyl else '_y.pdf'))
	plotHistoSliceSingle(hist, BIN_SET, 'z', z, proj2d, proj1d, show, out + ('_theta.pdf' if cyl else '_z.pdf'))

def plotDiffSliceMulti(h1, h2, BIN_SET, proj=False, cyl=True, diffFunc=getHistoDiffInv, name='', show=False, out=''):
	setShow(show)
	if(proj):
		k1 = h1.Project3D('yz')
		k2 = h2.Project3D('yz')

		l1 = h1.Project3D('zx')
		l2 = h2.Project3D('zx')

		m1 = h1.Project3D('yx')
		m2 = h2.Project3D('yx')
	else:
		x, y, z = getMaximum(h1, BIN_SET)
		k1 = getSlice(h1, BIN_SET, 'z', z)
		k2 = getSlice(h2, BIN_SET, 'z', z)

		l1 = getSlice(h1, BIN_SET, 'y', y)
		l2 = getSlice(h2, BIN_SET, 'y', y)

		m1 = getSlice(h1, BIN_SET, 'x', x)
		m2 = getSlice(h2, BIN_SET, 'x', x)

	k = diffFunc(k1, k2, name)
	l = diffFunc(l1, l2, name)
	m = diffFunc(m1, m2, name)

	plotHistoSliceMultiExecute(k, l, m, cyl, show, out)

def plotStandoffHisto(fName, type='ss', show=True, out=''):
        f = ROOT.TFile.Open(fName)
        hDict = openAllHistos(f)
        if type == 'ss':
            hData = hDict['DataSS']
            hMc = hDict['McSS']
        else:
            hData = hDict['DataMS']
            hMc = hDict['McMS']

        # hDataSS.Scale(1./hDataSS.Integral(5, 20))
        # hMcSS.Scale(1./hMcSS.Integral(5, 20))
        hData.Scale(1./hData.Integral())
        hMc.Scale(1./hMc.Integral())

        plotCompTwo(hMc, hData, 'Standoff Distance [mm]', 'Normalized Events / (2 mm)', 'MC', 'Data', show, out)

        f.Close()

def hist2dToList(h):
    x_bins = h.GetNbinsX()
    y_bins = h.GetNbinsY()

    bins = np.zeros((x_bins, y_bins))
    x, y = [], []

    for y_bin in range( y_bins ):
        for x_bin in range( x_bins ):
            bins[x_bin, y_bin] = h.GetBinContent(x_bin + 1, y_bin + 1)

            x.append( h.GetXaxis().GetBinLowEdge(x_bin + 1) )
            y.append( h.GetYaxis().GetBinLowEdge(y_bin + 1) )

    x = list( sorted(set( x )) )
    y = list( sorted(set( y )) )
    x.append( x[-1] + x[1] - x[0] )
    y.append( y[-1] + y[1] - y[0] )

    return x, y, bins.T

# === LOW BACKGROUND ===
def histToList(h):
    N = h.GetNbinsX()
    x = []

    binValList = []
    binErrList = []
    for i in range(1, N+1):
        binValList.append( h.GetBinContent(i) )
        x.append( h.GetBinLowEdge(i) )
        
    return x, binValList

def getListNorm(binValList, normRange=[0,-1]):
    binValList = np.array( binValList )
    binErrList = np.sqrt( binValList )

    sumValList = sum(binValList) # sum(binValList[normRange[0]:normRange[1]])
    binValList = binValList / sumValList
    binErrList = binErrList / sumValList

    return binValList, binErrList

def histToListNorm(h, normRange=[0, 200]):
    x, binValList = histToList(h)
    normRange = [np.digitize(normRange[0], x), np.digitize(normRange[1], x)]
    print normRange
    binValList, binValListErr = getListNorm(binValList, normRange)
    return x, binValList, binValListErr

# use not normalized data
def shapeComparison(dataMean, data):
    import plot_support as ps
    # Not normalized
    x, dataMean = histToList(dataMean)
    x, data = histToList(data)
    
    # Chi-square
    print 'Chi-square'
    print ps.chiSquareTest(dataMean, data)
    print 'Bhattacharyya'
    print ps.geometricTest(dataMean, data)
    print 'KS'
    print ps.ksTest(dataMean, data)

def plotLbkg(fName, fit=False, art='ss', output='default.pdf', N=3):
    from matplotlib import pyplot as plt
    import matplotlib.ticker as ticker
    '''
    from matplotlib import rc
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    '''

    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(output)

    if fit:
        import scipy.optimize
        import scipy.stats
        def linear(x, m):
            return m * x

    f = ROOT.TFile.Open(fName, 'READ')
    hDict = openAllHistos(f)

    typeList = ['z', 'a', 'so', 'r']
    xLabelList = ['z [mm]', 'Apothem [mm]', 'Standoff distance [mm]', 'Radius [mm]']
    yLabel = 'Normalized counts'
    xRangeList = [(-200, 200), (0, 180), (0, 180), (0, 200)]
    yRangeList = [(0.02, 0.045), (0.0, 0.14), (0.0, 0.172), (0.0, 0.13)]
    # normList = [(-200, 200), (130, 180), (0, 50), (140, 200)]
    normList = xRangeList

    # colorList = [('#50ab6d', '#ab62c0'), ('#648ace', '#979e3d'), ('#c86f3e', '#ca5576')]
    if N == 3:
        colorList = [('#ca5576', '#ca5576'), ('#648ace', '#648ace'), ('#50ab6d', '#50ab6d'), ('black', 'blue')]
    elif N == 6:
        colorList = [('#ca5576', '#ca5576'), ('#648ace', '#648ace'), ('#50ab6d', '#50ab6d'), ('#ab62c0', '#ab62c0'), ('#979e3d', '#979e3d'), ('#c86f3e', '#c86f3e'), ('black', 'blue')]
    else:
        return

    for j, typ in enumerate(typeList):
        # == PLOT ==
        # fig1 = plt.figure()
        # ax1 = fig1.add_subplot(111)
        fig1, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False, gridspec_kw = {'height_ratios':[2.5, 1]})

        # needed for KS test
        dataMean = hDict['%sData%s%d' % (typ, art.upper(), N)]
        normRange = normList[j]
        x, dataMeanCont, dataMeanContErr = histToListNorm(dataMean, normRange)
        for i in range(N + 1):
            data = hDict['%sData%s%d' % (typ, art.upper(), i)]
            mc = hDict['%sMc%s%d' % (typ, art.upper(), i)]

            x, dataCont, dataErrCont = histToListNorm(data, normRange)
            x, mcCont, mcErrCont = histToListNorm(mc, normRange)
            x = np.array( x )
            binWidth = data.GetBinWidth(1)

            # = Residuals =
            res = (dataCont - mcCont) / dataCont
            resErr = abs(mcCont/dataCont) * np.sqrt( (mcErrCont/mcCont)**2 + (dataErrCont/dataCont)**2 )

            # = Plot =
            ax1.step(list(x)+[(x[-1]+binWidth)], list(mcCont)+[0.], where='post', c=colorList[i][0])
            ax1.errorbar(x+binWidth/2., mcCont, yerr=mcErrCont, fmt='.', c=colorList[i][0], label='MC#%d' % i)
            ax1.errorbar(x+binWidth/2., dataCont, xerr=binWidth/2., yerr=dataErrCont, c=colorList[i][1], fmt='x', label='Data#%d' % i)
            ax2.errorbar(x+binWidth/2., res, xerr=binWidth/2., yerr=resErr, fmt='x', color=colorList[i][1])

            ax2.grid()
            ax2.axhline(y=0, linewidth=0.5, color='k')
            print x+binWidth/2.
            print dataCont

            # = Fit =
            if fit and typ == 'a':
                popt, pcov = scipy.optimize.curve_fit(linear, x[:-2]+binWidth/2., dataCont[:-2]) 
                perr = np.sqrt(np.diag(pcov))
                print 'popt: ', popt
                print 'perr: ', perr
                print
                ax1.plot(np.array([0,162]), linear(np.array([0, 162]), *popt),c=colorList[i][1])

            if fit and typ == 'so':
                shapeComparison(dataMean, data)

                D, p = scipy.stats.ks_2samp(dataMeanCont, dataCont)
                print D, p
                print

        ax2.set_xlabel( xLabelList[j] )
        # ax2.set_ylabel(r"$\frac{n_\text{Data} - n_\text{MC}}{n_\text{Data}}$")
        if isinstance(binWidth, (int, long)):
            ax1.set_ylabel( yLabel + ' / (%d mm)' % binWidth )
        else:
            ax1.set_ylabel( yLabel + ' / (%.1f mm)' % binWidth )

        plotLeft, plotRight = xRangeList[j]
        # plotBottom, plotTop = yRangeList[j]
        plotBottom, plotTop = min(dataCont)*0.9, max(dataCont)*1.1
        ax1.set_xlim(left=plotLeft, right=plotRight)
        ax1.set_ylim(bottom=plotBottom, top=plotTop)
        # ax1.set_yscale('log')
        ax1.legend(loc='best')

        ax2.set_ylim(bottom=-0.5, top=0.5)
        ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.25))

        fig1.show()
        raw_input('')
        plt.close()

        fig1.savefig(pp, format='pdf')            
        plt.cla()

    # = 2D histograms =
    # r vs. z
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    fig, ax = plt.subplots(2, N+1, sharex=True, sharey=True)

    typList = ['Data', 'Mc']
    
    for j, typ in enumerate( typList ):
        for i in range(N + 1):
            hData = hDict['%s%s%s%d' % ('rz', typ, art.upper(), i)]
            # Normalize
            # hData.Scale(1./hData.GetEntries())

            # Calculate residuial
            # h = hData.Clone('h%d' % i)
            # h.Add(hMC, -1)
            # h.Divide(hData)

            r, z, bins = hist2dToList( hData )
            print 'Bin normalization %s #%d' % (typ, i)
            # maxDataHist = hDict['%s%s%s%d' % ('rz', typ, art.upper(), N)]
            # rMax, zMax, binsMax = hist2dToList( maxDataHist )
            # maxData = maxDataHist.GetBinContent(maxDataHist.GetMaximumBin())
            bins = np.array(bins) / np.sum( bins ) # np.mean( binsMax )
            print bins[0:5], np.sum(bins)
            print

            from matplotlib.colors import LogNorm
            im = ax[j][i].imshow(bins, interpolation='nearest', extent=[r[0], r[-1], z[0], z[-1]], aspect='auto', cmap=plt.get_cmap('inferno'), norm=LogNorm(vmin=0.0001))
            if i == N:
                ax[j][i].set_title('%s Total' % typ)
            else:
                ax[j][i].set_title('%s #%d' % (typ, i))

    fig.text(0.5, 0.02, 'r [mm]', ha='center')
    fig.text(0.02, 0.5, 'z [mm]', va='center', rotation='vertical')
    # plt.xlabel( 'r [mm]' )
    # plt.ylabel( 'z [mm]' )

    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.subplots_adjust(right=0.8)
    fig.colorbar(im, cax=cbar_ax)
    fig.show()

    raw_input('')
    fig.savefig(pp, format='pdf')

    f.Close()
    pp.close()

# === THETA PLOTS ===
def plotZtheta(fNdata, fNmc):
    import matplotlib.pyplot as plt
    dData = fileToDict(fNdata)
    dMC = fileToDict(fNmc)

    keyList, diffList = [], []
    for key, hData in dData.iteritems():
        keyList.append( key )

        try:
            hData.Scale(1./hData.Integral())
        except:
            pass

        hMC = dMC[key]
        try:
            hMC.Scale(1./hMC.Integral())
        except:
            pass

        diffList.append( getHistoDiffInv(hData, hMC, str(key)).GetBinContent(1) )

    diff = np.array( diffList )
    z = np.array( keyList )[:,0]
    try:
        theta = np.array( keyList )[:,1]
        tz = True
    except:
        theta = z
        tz = False

    if tz:
        print sorted( z )
        print len( set( z ) ), len( set( theta ) )

        print np.sqrt(len(theta))
        # diff = diff.reshape(38, 11)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(theta, z, diff)
        # plt.imshow(diff)
        # plt.colorbar()
        plt.show()

    else:
        plt.plot(theta, diff, 'o')
        plt.show()

# Store (z, theta) root-file to dict
def fileToDict(fN):
    import re
    f = ROOT.TFile.Open(fN, 'READ')
    hDict = openAllHistos(f)

    d = {}
    for key, h in hDict.iteritems():
        k = [float(i) for i in re.findall(r"[-+]?\d*\.\d+|\d+", key)]
        d[tuple(k)] = h

    return d

