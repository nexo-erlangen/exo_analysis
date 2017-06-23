#!/usr/bin/env python
import numpy as np
import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')

import artificialField as af
import plot_support as ps

REFLECTORINNERRAD = 183.2356
CATHODE_ANODE_y_DISTANCE = 192.23657

def main():
    fName = '/home/vault/capm/mppi025h/analysis/preprocess/pre_testSingle_pre.root'
    fNameSSonly = '/home/vault/capm/mppi025h/analysis/preprocess/pre_testSingleDrift_pre.root'
    # fNameFid = '/home/vault/capm/mppi025h/analysis/preprocess/pre_testSingleFid_pre.root'
    # fNameSSonlyFid = '/home/vault/capm/mppi025h/analysis/preprocess/pre_testSingleDriftFid_pre.root'
    
    # fName = '/home/vault/capm/mppi025h/analysis/preprocess/S5Th/Xe134_run_4009_tree.root'
    fName = '/home/vault/capm/mppi025h/analysis/preprocess/pre_artDrift_pre.root'
    fNameSSonly = '/home/vault/capm/mppi025h/analysis/preprocess/SourceS5_Th228.root'

    # msCluster( fName )
    # driftCluster( fName )

    '''
    hPcd = pcd2d( fName )
    hPcdDrift = pcd2d( fNameSSonly )
    hPcd.Scale(1./hPcd.Integral())
    hPcdDrift.Scale(1./hPcdDrift.Integral())

    c = ROOT.TCanvas()
    hDiff = ps.getHistoDiffInv(hPcd, hPcdDrift, '')
    hDiff.Draw('colz')
    # hPcd.Draw('colz')
    # c1 = ROOT.TCanvas()
    # hPcdDrift.Draw('colz')

    # hPcdFid = pcd2d( fNameFid )
    # hPcdDriftFid = pcd2d( fNameSSonlyFid )
    # hPcdFid.Scale(1./hPcdFid.Integral())
    # hPcdDriftFid.Scale(1./hPcdDriftFid.Integral())

    # hDiffFid = ps.getHistoDiff(hPcdFid, hPcdDriftFid, '')

    # c2 = ROOT.TCanvas()
    # hDiffFid.Draw('colz')

    # hPcdFid.Draw('colz')
    # c3 = ROOT.TCanvas()
    # hPcdDriftFid.Draw('colz')
    raw_input('')
    return
    '''

    hArt = pcdStandoff( fName )
    hArtSSonly = pcdStandoff( fNameSSonly )
    # hArtFid = pcdStandoff( fNameFid )
    # hArtSSonlyFid = pcdStandoff( fNameSSonlyFid )

    hArtSSonly.SetMarkerStyle(2)
    # hArtSSonlyFid.SetMarkerStyle(2)
    
    hArt.Scale(1./hArt.Integral())
    hArtSSonly.Scale(1./hArtSSonly.Integral())
    # hArtFid.Scale(1./hArtFid.Integral())
    # hArtSSonlyFid.Scale(1./hArtSSonlyFid.Integral())

    c = ROOT.TCanvas()
    ROOT.gStyle.SetOptStat(0)
    leg = ROOT.TLegend(.7, .8, .9, .9)
    hArt.Draw()
    hArtSSonly.Draw('Psame')

    leg.AddEntry(hArt, 'normal')
    leg.AddEntry(hArtSSonly, 'drifted')
    leg.Draw()

    '''
    c2 = ROOT.TCanvas()
    leg2 = ROOT.TLegend(.7, .8, .9, .9)
    hArtFid.Draw('')
    hArtSSonlyFid.Draw('Psame')

    leg2.AddEntry(hArtFid, 'normal fid')
    leg2.AddEntry(hArtSSonlyFid, 'drifted fid')
    leg2.Draw()
    '''

    raw_input('')

def msCluster(fName):
    f = ROOT.TFile.Open(fName, 'READ')
    to = f.Get('mcTree')
    ROOT.gROOT.cd()

    cut = ps.getCut(calibCut=True, energyCut=True, type='ms', MC=True, eMin=750, eMax=3500)
    t = to.CopyTree( cut )

    N = t.GetEntries()
    print 'Processing %d events...' % N

    clusterPos = []
    for i in range(N):
        t.GetEntry( i )
        es = t.EventSummary

        mul = es.multiplicity
        mulO = mul
        print i, mul,
        if mul > 1.1:
            posX, posY, posZ, energy = np.array( es.cluster_x ), np.array( es.cluster_y ), np.array( es.cluster_z ), np.array( es.cluster_energy )
            
            '''
            if len( posX ) == 1:
                print 'What?'
                clusterPos.append( [posX[0], posY[0], posZ[0]] )
            '''

            clusterSurvive = []
            for i in range( len(posX) ):
                x, y, z, en = posX[i], posY[i], posZ[i], energy[i]
                x_, y_, z_, t_ = af.artificialDrift_C(x, y, z, 0 )

                if x_ >= 900:
                    mul -= 1
                else:
                    clusterSurvive.append( [x_, y_, z_] )

            '''
            if mul != mulO:
                print 'Hit!'
            else:
                print
            '''

            if mul == 1:
                try:
                    print clusterSurvive[0]
                    clusterPos.append( clusterSurvive[0] )
                except:
                    pass
            else:
                print

    h = ROOT.TH1F('h', 'h', 50, 0, 200)
    for item in clusterPos:
        r = np.sqrt( item[0]**2 + item[1]**2 )
        z = item[2]

        print (REFLECTORINNERRAD-r, z)
        d_r = REFLECTORINNERRAD - r
        d_z = CATHODE_ANODE_y_DISTANCE - abs(z)

        if d_r < d_z:
            h.Fill( d_r )
        else:
            h.Fill( d_z )

    c = ROOT.TCanvas()
    h.Draw()

    raw_input('')

    f.Close()

def driftCluster(fName):
    f = ROOT.TFile.Open(fName, 'READ')
    to = f.Get('mcTree')
    ROOT.gROOT.cd()

    cut = ps.getCut(calibCut=True, energyCut=False, type=None, MC=True, eMin=750, eMax=3500)
    t = to.CopyTree( cut )

    N = t.GetEntries()
    h = ROOT.TH1F('h', 'h', 40, 0, 200)

    for i in range( N ):
        t.GetEntry(i)
        es = t.EventSummary

        x, y, z, en = np.array( es.cluster_x )[0], np.array( es.cluster_y )[0], np.array( es.cluster_z )[0], np.array( es.cluster_energy )[0]
        x_, y_, z_, t_ = af.artificialDrift_C(x, y, z, 0 )

        if abs( x_ ) < 900 and z_ < 0:
            r = np.sqrt( x_**2 + y_**2 )
            d_r = REFLECTORINNERRAD - r
            d_z = CATHODE_ANODE_y_DISTANCE - abs(z_)
    
            if d_r < d_z:
                h.Fill( d_r )
            else:
                h.Fill( d_z )

    c = ROOT.TCanvas()
    h.Scale(1./h.Integral())
    h.Draw()
    raw_input('')

def pcdStandoff(fName):
    f = ROOT.TFile.Open(fName, 'READ')
    to = f.Get('mcTree')
    if not to:
        to = f.Get('dataTree')

    h = ROOT.TH1F('h', 'h', 40, 0, 200)

    ROOT.gROOT.cd()
    cut = '' # ps.getCut(calibCut=True, energyCut=True, type=None, MC=True, eMin=750, eMax=3500)
    t = to # to.CopyTree( cut )

    N = 100000 # t.GetEntries()
    if N > t.GetEntries():
        N = t.GetEntries()
    j = 0
    
    # for i in range(N):
    while N > 0:
        t.GetEntry( j )
        if not j % 10000:
            print N, j
        j += 1

        es = t.EventSummary

        if es.multiplicity > 1.1: # or es.isMissingPosition:
           continue

        pcdX, pcdY, pcdZ = np.array( es.pcd_x ), np.array( es.pcd_y ), np.array( es.pcd_z )
        # pcdX, pcdY, pcdZ = np.array( es.cluster_x ), np.array( es.cluster_y ), np.array( es.cluster_z )
        for i in range( len(pcdX) ):
            x, y, z = pcdX[i], pcdY[i], pcdZ[i]
            N -= 1
            r = np.sqrt( x**2 + y**2 )

            if (r > REFLECTORINNERRAD) or not ps.isFiducial(x, y, z, 162, 5, 182) or abs(x) > 900: # or z > 0:
                continue

            d_r = REFLECTORINNERRAD - r
            d_z = CATHODE_ANODE_y_DISTANCE - abs(z)
            if d_r < d_z:
                h.Fill( d_r )
            else:
                h.Fill( d_z )
        
    return h

def pcd2d(fName):
    f = ROOT.TFile.Open(fName, 'READ')
    to = f.Get('mcTree')
    if not to:
        to = f.Get('dataTree')
    h = ROOT.TH2F('h', 'h', 50, -190, 190, 50, -190, 190)

    ROOT.gROOT.cd()
    cut = '' # ps.getCut(calibCut=True, energyCut=False, type='ss', MC=True, eMin=750, eMax=3500)
    t = to # to.CopyTree( cut )

    N = 1000000 # t.GetEntries()
    if N > t.GetEntries():
        N = t.GetEntries()
    j = 0
    
    # for i in range(N):
    while N > 0:
        t.GetEntry( j )
        if not j % 10000:
            print N, j
        j += 1
        es = t.EventSummary

        # pcdX, pcdY, pcdZ = np.array( es.pcd_x ), np.array( es.pcd_y ), np.array( es.pcd_z )
        pcdX, pcdY, pcdZ = np.array( es.cluster_x ), np.array( es.cluster_y ), np.array( es.cluster_z )
        for i in range( len(pcdX) ):
            x, y, z = pcdX[i], pcdY[i], pcdZ[i]
            if not ps.isFiducial(x, y, z, 162, 5, 182): # or z > 0:
                continue

            N -= 1
            h.Fill(x, y)

            '''
            r = np.sqrt( x**2 + y**2 )

            d_r = REFLECTORINNERRAD - r
            d_z = CATHODE_ANODE_y_DISTANCE - abs(z)
            if d_r < d_z:
                if d_r > 10 and d_r < 30:
                    h.Fill(x, y)
            else:
                if d_z > 10 and d_z < 30:
                    h.Fill(x, y)
            '''
    c = ROOT.TCanvas()
    return h

if __name__ == '__main__':
    main()

