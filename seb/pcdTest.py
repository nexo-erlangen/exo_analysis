#!/usr/bin/env python
import numpy as np
import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')

import plot_functions as pf

H_BINS_ALL = 100
H_BINS = [H_BINS_ALL, -200, 200, H_BINS_ALL, -200, 200, H_BINS_ALL, -200, 200] 

def main():
    f = ROOT.TFile('pcdDistro.root', 'READ')
    drawHist(f, H_BINS)
    f.Close()
    return

    L = 0.16
    N = None # 100000

    fName = '/home/vault/capm/mppi025h/analysis/preprocess/SourceS5_Th228.root'
    posList = getSSPositions(fName, L, N)

    f = ROOT.TFile('pcdDistro0_16.root', 'RECREATE')
    h = ROOT.TH3F('h', 'h', *H_BINS)

    for pos in posList:
        x, y, z = pos
        h.Fill(x, y, z, 1)

    # h.Project3D('yz').Draw('colz')
    h.Write()

    # drawHist( f, H_BINS )

    f.Close()

def getSSPositions(fName, L, N):
    f = ROOT.TFile.Open(fName, 'READ')
    t = f.Get('mcTree')
    if not N:
        N = int( t.GetEntries() * 0.1 )

    print 'Filtering SS data from %d events' % N
    print '======'

    ssEvents = []

    for i in range(N):
        t.GetEntry( i )
        es = t.EventSummary
        
        eventLs = []
        for j in range( len( np.array( es.pcd_x ) ) ):
            eventLs.append( [np.array( es.pcd_x )[j], np.array( es.pcd_y )[j], np.array( es.pcd_z )[j]] )

        # print eventLs
        for j in range( len(eventLs) - 1):
            x = np.array( eventLs[j] )
            y = np.array( eventLs[j+1] )

            d = np.linalg.norm( x - y )
            if d > L:
                break
        else:
            try:
                ssEvents.append( list( np.sum( [ np.array( item ) for item in eventLs ], axis=0 ) / len( eventLs ) ) )
            except:
                pass
            # print 'SS detected!'
        # print

    #for item in ssEvents:
    #    print item

    print
    print '%d Events for length %.2f' % (len(ssEvents), L)
    print 'SS-Fraction: %.2f' % (float( len(ssEvents) ) / N)

    f.Close()

    return ssEvents

def drawHist(f, H_BINS):
    h = f.Get('h')

    c = ROOT.TCanvas()
    h.Project3D('zx').Draw('COLZ')

    '''
    pf.plotHistoSliceMulti(h, H_BINS, proj=True, cyl=False, show=True, out='')
    pf.plotHistoSliceSingle(h, H_BINS, option='x', value=0, proj2d=False, proj1d=True, show=True, out='')
    pf.plotHistoSliceSingle(h, H_BINS, option='y', value=0, proj2d=False, proj1d=True, show=True, out='')
    pf.plotHistoSliceSingle(h, H_BINS, option='z', value=0, proj2d=False, proj1d=True, show=True, out='')
    '''

    raw_input('')

if __name__ == '__main__':
    main()

