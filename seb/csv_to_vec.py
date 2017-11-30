#!/usr/bin/env python
try:
	import matplotlib.pyplot as plt
	from matplotlib.colors import LogNorm
        from matplotlib import cm
        from matplotlib import rc
        from matplotlib.ticker import NullFormatter
        rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        rc('text', usetex=True)
        import seaborn as sns
        sns.setstyle()
except:
	pass

import numpy as np
import csv
import cPickle

import ROOT

import generate_random_new_ as gr

def main():
	# storeEField('./efield_data/', 'fieldExport_converted.csv', None)
	
	# posList, vecList = storeEFieldToList('./efield_data/', 'fieldExport_converted.csv')
	# posList, vecList = storeEFieldToList('./efield_data/', 'export2D_converted.csv')
	posList, vecList = storeEFieldToList('./efield_data/', 'export2D_AvgPotConverted.csv')

	#storeListToFile(posList, './efield_data/posList.p')
	#storeListToFile(vecList, './efield_data/vecList.p')

	plotEField(posList, vecList, getH_BINS(posList, vecList), component=None) 

def storeListToFile(ls, fName):
	cPickle.dump(ls, open(fName, 'wb'))

def storeEFieldToList(inDir, inFile):
	posList = []
	vecList = []
	with open(inDir + inFile, 'rb') as f:
		reader = csv.reader(f)
		for row in reader:
			if row[0].startswith('%'):
				continue
			
			row = [ float(item) for item in row ]

			pos = [ item for item in row[:3] ]
			# pos[2] -= 0.203796
			pos[2] -= 0.198
			posList.append( pos )

			vec = row[3:]
			vec = list( np.array(vec) * -1 )
			vecList.append( vec )

	return posList, vecList

def getH_BINS(posList, vecList):
        yLs = list( np.array( posList )[:,1] )
        zLs = list( np.array( posList )[:,2] )

        # print zLs
        print 'Elements to count:', (yLs[0], zLs[0])
        print 'Length of lists:', (len(yLs), len(zLs))

        q = yLs.count( yLs[0] )
        p = zLs.count( zLs[0] )
        H_BINS = [1, 0, 1, p, min( yLs ), max( yLs ), q, min( zLs ), max( zLs )]

        return H_BINS

def storeEField(inDir, inFile, H_BINS=None, verbose=False):
	posList, vecList = storeEFieldToList(inDir, inFile)

	noHbinFlag = False
	if not H_BINS:
		noHbinFlag = True	
		H_BINS = getH_BINS(posList, vecList)

	binW = gr.getBinData(H_BINS)[-1]

	print H_BINS
	print binW

	histPos = ROOT.TH3D('Epos', 'Epos', *H_BINS)
	histX = ROOT.TH3D('Ex', 'Ex', *H_BINS)
	histY = ROOT.TH3D('Ey', 'Ey', *H_BINS)
	histZ = ROOT.TH3D('Ez', 'Ez', *H_BINS)

	histList = [histX, histY, histZ]

	for i, hist in enumerate(histList):
		for j, pos in enumerate(posList):
			x, y, z = pos
			hist.Fill(x + binW[0]/2, y - binW[1]/2, z - binW[2]/2, vecList[j][i])

	if not noHbinFlag:
		for hist in histList:
			hist.Divide(histPos)
	
	f = ROOT.TFile(inDir + '/efield_hist2d.root', 'recreate')
	for hist in histList:
		hist.Write()

	f.Close()

	return H_BINS

def plotEField(posList, vecList, H_BINS, vmin=10.e3, component=None):
        # posList, vecList = zip(*[item for item in zip(posList, vecList) if (abs(item[0][2]) < 0.2 and item[0][2] < 0.)]) 
        print posList[0]
	print H_BINS
	ext = np.array( [H_BINS[4], H_BINS[5], H_BINS[7], H_BINS[8]] ) * 1000    # convert m to mm
        # ext = np.array( [-0.2, 0.2, -0.2, 0] ) * 1000
	q = int( H_BINS[3] )
	vecListChunk = list(chunks(vecList, q))

	vecList = []
	for item in vecListChunk:
                if not component:
                    vecList.append( [np.linalg.norm(x)/100 for x in item] )
                else:
                    idx = 0
                    if component == 'y':
                        idx = 1
                    elif component == 'z':
                        idx = 2

		    vecList.append( [ x[idx]/100 for x in item ] )

	fig, ax = plt.subplots(1)
	# im = ax.imshow(vecList, extent=ext, interpolation='bicubic', cmap='viridis', origin='lower') # , norm=LogNorm(vmin=vmin))
        import matplotlib as mpl
        logthresh = 1
        from matplotlib import cm
        cmap = cm.get_cmap('Blues', 11)
	im = ax.imshow(vecList, extent=ext, interpolation='bicubic', cmap=cmap, origin='lower', vmin=0, vmax=1000) # , norm=mpl.colors.SymLogNorm( logthresh ))
        # im = ax.imshow(vecList, extent=ext, interpolation='bicubic', cmap='viridis', origin='lower', vmin=-10.e3, vmax=0.)
        # im, cb = imshow_symlog(fig, ax, vecList, ext, -800, 800, 1)
        ax.set_xlabel('x [mm]')
        ax.set_ylabel('z [mm]')
        ax.set_xlim(0, 195)
        ax.set_ylim(-200, 0)

        vmin, vmax = 800, 800
        # maxlog=int(np.ceil( np.log10(vmax) ))
        # minlog=int(np.ceil( np.log10(-vmin) ))
        # tick_locations=([-(10**x) for x in xrange(minlog,-logthresh-1,-1)] +[0.0] + [(10**x) for x in xrange(-logthresh,maxlog+1)] )
        cb = fig.colorbar(im)
        cb.set_label(r'$E$ [V/cm]')
	#cbar = fig.colorbar(im)
	plt.show()

	raw_input('')
	return 

def imshow_symlog(fig, ax, my_matrix, extent, vmin, vmax, logthresh=5):
    import matplotlib as mpl
    img = ax.imshow( my_matrix , extent=extent,
        vmin = float(vmin), vmax=float(vmax),
        norm = mpl.colors.SymLogNorm(10**-logthresh),
        interpolation = 'bicubic')

    maxlog=int(np.ceil( np.log10(vmax) ))
    minlog=int(np.ceil( np.log10(-vmin) ))

    #generate logarithmic ticks 
    tick_locations=([-(10**x) for x in xrange(minlog,-logthresh-1,-1)] +[0.0] + [(10**x) for x in xrange(-logthresh,maxlog+1)] )

    cb = fig.colorbar(img, ticks=tick_locations)
    return img, cb

def chunks(l, n):
	for i in xrange(0, len(l), n):
		yield l[i:i + n]

if __name__ == '__main__':
	main()

