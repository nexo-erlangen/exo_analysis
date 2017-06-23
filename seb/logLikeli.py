#!/usr/bin/env python
import numpy as np
import sortData as sd
import plot_support as ps
import cPickle
import os

# Fiducial volume
FV = [162, 5, 182]

# zLimits = [(-160, -20), (20, 160)]
zLimitsBottom = (-160, -20)
zLimitsTop = (20, 160)

llN = 2 # Number of bins to get logLikelihood value from

# Paths
preProcessDir = '/home/vault/capm/mppi025h/analysis/preprocess/'
dataList = []
for fn in os.listdir(preProcessDir + 'S5Th/'):
    #runNum = [int(s) for s in fn.split('_') if s.isdigit()][0] 
    #if runNum >= 4500 and runNum <= 4600:
    dataList.append( preProcessDir + 'S5Th/%s' % fn )

# TODO: set
# mcList = [preProcessDir + 'SourceS5_Th228.root']
# mcList = [preProcessDir + 'pre_artDrift_neg_pre.root']
mcList = [preProcessDir + 'pre_nelderMead_pre.root']

def main():
    cnt = 0
    exoRootOut = '/home/vault/capm/mppi025h/plot_scripts/NMExoAnalysis/'
    print getLogLikelihood(mcList, outFn='Top%d' % cnt, outDir=exoRootOut, art='ss')
    # print getLogLikelihood(mcList, outFn='Top%d' % cnt, outDir=exoRootOut, art='ms')

def getLogLikelihood(mcFile, outFn='Standoff', outDir='./', art='ss'):
    nBins = 120
    
    dataChain = ps.getChain('dataTree', dataList)
    mcChain = ps.getChain('mcTree', mcList)

    # getPlot creates data dictionaries
    sd.sortDataSave(dataChain, mcChain, FV, nBins, outFn + '%s' % art.upper(), art, outDir)
    hMCDict, hDataDict, hMCDictZ, hDataDictZ, hMCDictTheta, hDataDictTheta = sd.getDictsReal(nBins, outFn + '%s' % art.upper(), outDir)
    xGrid, yGrid, valMCNorm, valDataNorm, valDiffNorm, valDiffNormError = sd.plotStandoff(hMCDictZ, hDataDictZ, nBins, 0, True, False)
    # sd.getPlot(False, nBins, 0, outFn + '%s' % art.upper(), outDir)
    # d = cPickle.load(open('Standoff%s.p' % art.upper(), 'rb')) 

    llSumTop, llSumBottom = 0, 0
    for i in range(llN):
        z = np.array( xGrid )

        mcTop, mcBottom = cutData(z, np.array( valMCNorm )[:,i], zLimitsTop), cutData(z, np.array( valMCNorm )[:,i], zLimitsBottom)
        dataTop, dataBottom = cutData(z, np.array( valDataNorm )[:,i], zLimitsTop), cutData(z, np.array( valDataNorm )[:,i], zLimitsBottom)

        llResTop, llResBottom = poissonLogLikelihood(dataTop, mcTop), poissonLogLikelihood(dataBottom, mcBottom)
        llSumTop += llResTop 
        llSumBottom += llResBottom

    return (float(llSumTop) / llN), (float(llSumBottom) / llN)

def cutData(z, data, limits):
    dataCut = []
    a, b = limits
    for i, d in enumerate(data):
        if z[i] >= a and z[i] <= b:
            dataCut.append( d )

    return dataCut

def poissonLogLikelihood(data, mc):
    return sum( np.nan_to_num( [data[i] - mc[i]*np.log( data[i] ) if (data[i] and mc[i]) else 0 for i in range( len(data) )] ) )

if __name__ == '__main__':
    main()

