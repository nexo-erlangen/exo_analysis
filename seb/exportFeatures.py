import numpy as np
import ROOT

def main():
    getCountsData('./standoffArtDrift_neg3Uncut/dataStandoffSSZ.root', 'dataValuesSS.dat')
    getCountsData('./standoffArtDrift_neg3Uncut/mcStandoffSSZ.root', 'mcValuesSS.dat')
    return

    f = ROOT.TFile.Open('realNodiffFeatures.root', 'READ')

    hData = f.Get('standoff_distanceDataSS')
    hMC = f.Get('standoff_distanceMCSS')

    fOut = open('exportFeatures.dat', 'w')
    fOut.write('# data\t mc\t (data-mc)/data\n')
    for i in range(1, int(hData.GetSize()) - 2):
        data = hData.GetBinContent(i)
        mc = hMC.GetBinContent(i)
        fOut.write('%f\t%f\t%f\n' % (data, mc, (data-mc)/data))
    fOut.write('\n')
    fOut.close()

def getCountsData(fName, fOutName):
    f = ROOT.TFile.Open(fName, 'READ')
    
    hNameList = f.GetListOfKeys()
    entryList = []

    for h in hNameList:
        hName = h.GetName()
        print 'Read histogram %s' % (hName)
        z = float( hName.split(' = ')[-1] )
        hist = f.Get(hName)

        N = hist.GetBinContent(1)
        try:
            hist.Scale(1./hist.Integral())
        except:
            pass
        n = hist.GetBinContent(1)

        entryList.append( (z, int(N), n) )

    entryList = sorted(entryList, key=lambda x:x[0])
    fOut = open(fOutName, 'w')
    fOut.write('z\tN\tn\n')
    for entry in entryList:
        fOut.write('%f\t%d\t%f\n' % entry)
    fOut.write('\n')
    fOut.close()

if __name__ == '__main__':
    main()

