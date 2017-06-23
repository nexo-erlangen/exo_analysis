#!/usr/bin/env python
import ROOT

def main():
    f = ROOT.TFile.Open('standoffDrift.root', 'OPEN')

    h = f.Get('(randHist)')
    hCut = f.Get('(randHistCut)')
    hDrift = f.Get('(randDriftHist)')
    hDriftCut = f.Get('(randDriftHistCut)')

    c = ROOT.TCanvas()
    h.Project3D('yx').Draw('colz')
    cCut = ROOT.TCanvas()
    hCut.Project3D('yx').Draw('colz')
    cDrift = ROOT.TCanvas()
    hDrift.Project3D('yx').Draw('colz')
    cDriftCut = ROOT.TCanvas()
    hDriftCut.Project3D('yx').Draw('colz')

    raw_input('')
    f.Close()

if __name__ == '__main__':
    main()

