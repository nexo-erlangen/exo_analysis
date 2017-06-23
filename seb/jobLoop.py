#!/usr/bin/env python
import subprocess
import os
import sys
import time
import shutil
import cPickle
import fileinput
import itertools
import argparse
import numpy as np

import gen_ran_qsub as grq
import sortData as sd
import peak_finder as peak
import plot_support as ps
import curve_evaluation as ce

# === Constants ===
# N = int( 1.e6 )        # Number of events
N = int( 10.e6 )

settingsF = 'settings_mcSource.py'
outF = 'standoffDrift.root'
nBins = 120

# === EXO Script ===
exoScript = '/home/vault/capm/mppi025h/EXOoffline/exoout_/macros/S5.exo'
exoAnalysisOut = '/home/vault/capm/mppi025h/run_data/MC/nelderMead/'
exoRootOut = '/home/vault/capm/mppi025h/plot_scripts/NMExoAnalysis/'

USER = 'mppi025h'

# === Folder ===
resultDir = 'result/'
pdfDir = resultDir + 'pdf/'
dataDir = resultDir + 'data/'
gpDir = resultDir + 'gp/'

# Used by Job Loop
resumeFile = resultDir + 'jobLoopResume.p'
resultFile = resultDir + 'results.dat'

# Used by Nelder Mead
resultDB = resultDir + 'resultDB.p'
resultDBTop = resultDir + 'resultDBTop.p'
resultDBBottom = resultDir + 'resultDBBottom.p'
resultFileNM = resultDir + 'resultsNM.dat'
resultFileNMTop = resultDir + 'resultsNMTop.dat'
resultFileNMBottom = resultDir + 'resultsNMBottom.dat'

def main():
    global lima
    resume, nmFlag, lima, exo = get_args()

    # Check if folders exist, if not create
    makeDir( resultDir )
    makeDir( pdfDir )
    makeDir( dataDir )
    makeDir( gpDir )

    # == JOB LOOP ==
    # Loop over all parameters given in parDictList
    if not nmFlag:
        # Parameter ranges
        if resume and os.path.isfile( resumeFile ):
            parDictList = cPickle.load( open(resumeFile, 'rb') )
        else:
            parDictList = []

        if not parDictList:
            parList = [ ['c_m', 1.1, 1.1, 1], ['c_A', 5, 13, 3], ['c_m2', 0.04, 0.06, 2] ]
            parDictList = getParamDictList( parList )

        parDictListFinished = []        # Holds already processed parameters

        # Make sure all jobs are finished before starting new ones
        qstatWait()

        if not resume:
            f = open(resultFile, 'w')
            for item in sorted( list( parDictList[0].keys() ) ):
                f.write( '%s\t' % item )
            f.write( '\n' )

        for parDict in parDictList:
            outName, outNameCut = evalParams( parDict )
            result = ce.evalPara( outName, 'real/real_60.p' ) 
            resultCut = ce.evalPara( outNameCut, 'real/realCut_60.p', True )
            result = 0.25 * (result + 3*resultCut)

            # Store results to file
            storeResult( f, parDict )

            # Store to be processed paramaters in file
            parDictListFinished.append( parDict )
            dDump = []
            for d in parDictList:
                if not d in parDictListFinished:
                    dDump.append( d )
            cPickle.dump( dDump, open(resumeFile, 'wb') )

        f.close()

    # == NELDER-MEAD ==
    else:
        if not exo:
            # initialGuessOld = [ ('c_m', 1.1), ('c_m2', 0.06), ('s_m', 0.175), ('c_A', 5) ]
            # initialGuessOld2 = [ ('c_m', 1.084), ('c_m2', 0.079), ('s_m', 0.222), ('c_A', 5.047) ]
            # initialGuess = [ ('c_m', 1.1), ('c_m2', 0.06), ('s_m', 0.175), ('c_A', 13.) ]
            # initialGuess = [ ('c_m', 1.5), ('c_m2', 0.01), ('s_m', 0.175), ('c_A', 20.), ('t_h', 0.05) ]
            # initialGuess = [ ('c_m', 1.5), ('c_m2', 0.01), ('s_m', 0.175), ('c_A', 20.), ('t_h', 0.2) ]
            
            # Initial Guess for artDrift_neg4
            # initialGuess = [ ('c_m', -2.), ('c_A', 6.25), ('t_h', 1.2) ]

            # Initial Guess for using only data on the side
            initialGuess = [ ('c_m', -1.59), ('c_A', 7.67), ('t_h', 1.177) ]

            bounds = {}
            # bounds['s_m'] = (0.172, 0.183)
            bounds['t_h'] = (0., 2.)
            bounds['c_m'] = (-5., 0.)
            # bounds['c_m2'] = (-0.02, 0.02)
            bounds['c_A'] = (0., 30.)

            # initialGuess = [ ('c_m', 1.1), ('c_A', 5) ]
            initialDict = getParamDict(initialGuess)
            # qstatWait()

            # Main loop
            nelderMead(initialDict, bounds)

        # NELDER-MEAD & EXOANALYSIS
        else:
            # Boundary conditions
            bounds = {}
            bounds['t_h'] = (0., 2.)
            bounds['c_m'] = (-5., 0.)
            bounds['c_A'] = (0., 30.)

            # Initial values
            initialGuessTop = [ ('c_A', 7.671374), ('c_m', -1.593037), ('t_h', 1.177416) ]
            initialDictTop = getParamDict(initialGuessTop)

            initialGuessBottom = [ ('c_A', 7.671374), ('c_m', -1.593037), ('t_h', 1.177416) ]
            initialDictBottom = getParamDict(initialGuessBottom)

            # Main
            nelderMeadEXOAnalysis(initialDictTop, initialDictBottom, bounds)

def evalParamsEXOAnalysis(parDictTop, parDictBottom):
    import create_jobs as cj
    print 'Creating job scripts...'

    # Clear old files
    if os.path.isdir(exoAnalysisOut):
        shutil.rmtree(exoAnalysisOut)

    # Write parameter to artifical drift
    writeArtificialParam(parDictTop, 'artificialField_paramEXOTop.py')
    writeArtificialParam(parDictBottom, 'artificialField_paramEXOBottom.py')
    freeNodes = getFreeNodes()

    # Set number of jobs
    if lima:
        job = 50
        jobFac = 0.6
    else:
        job = 20
        jobFac = 0.6

    if freeNodes*jobFac > job:
        nJobs = int( jobFac * freeNodes )
    else:
        nJobs = int( job )

    # Number of events N is defined in header
    cj.createJobs(exoScript, N, nJobs, exoAnalysisOut, 'nelderMead', lima)

    # Submit scripts and wait
    print '\nSubmitting Jobs...'
    print exoAnalysisOut
    jobNr, err = sendJobs('%s/run_submit.sh' % exoAnalysisOut)
    if err:
        print err
        return
    qstatWait(jobNr, 60, True)

    # Once EXOAnalysis has finished, the MC files 
    # need to be preprocessed. Since lima gets an error,
    # switch to woodycap for execution
    cmd = ("ssh -Y %s@woodycap.rrze.uni-erlangen.de -t '%s'" % (USER, exoAnalysisOut + '/preprocess.sh'))
    out, err = sendCmd( cmd )
    if err:
        print err
        return

    # The resulting file lies in $VAULT/analysis/preprocess
    # and is named pre_nelderMead_pre.root
    return '$VAULT/analysis/preprocess/pre_nelderMead_pre.root'

def evalParams( parDict ):
    # Set filenames according to parameter
    outName = 'so'
    for key, value in parDict.iteritems():
        val = ('%.3f' % value).replace('.', '_')
        outName += '-%s%s' % (key, val)
    print 'Processing %s...' % outName
    outNameCut = outName + 'Cut'

    # Write parameter to artificial drift 
    writeArtificialParam(parDict, 'gp_data/artificialField_param.py')

    freeNodes = getFreeNodes()

    # Set number of jobs to be used
    if lima:
        job = 60
        jobFac = 0.3
    else:
        job = 20
        jobFac = 0.6

    if freeNodes > job*jobFac:
        nJobs = int( jobFac * freeNodes )
    else:
        nJobs = int( job*jobFac )

    # Get number of events per job
    nEventsPerJob = int( float(N) / nJobs )

    # Generate qsub-scripts and output folder
    # Clear old data first
    if lima:
        if os.path.isdir('drift_data_lima'):
            shutil.rmtree('drift_data_lima')
    else:
        if os.path.isdir('drift_data'):
            shutil.rmtree('drift_data')
    
    print 'Generating qsub scripts for %d jobs with %d events per job' % (nJobs, nEventsPerJob)
    grq.execGen(settingsF, outF, nEventsPerJob, nJobs, lima)

    # Submit the scripts and wait till they're finished
    print '\nSubmitting scripts...'

    # grq.genRand(1000, 'standoffDrift.root')
    jobNr, err = sendJobs('ran_qsub.sh')
    if err:
        print err
        return
    qstatWait(jobNr, 60, True)

    # Once jobs are finished, combine results
    grq.getStat(True, lima)
    # os.rename('randPos_0.p', 'randPos.p')
    # os.rename('randDriftPos_0.p', 'randDriftPos.p')

    # = WITHOUT CUT =
    # Create files for z-standoff plot
    sd.sortDataSaveRand('randPos.p', 'randDriftPos.p', nBins)

    # Create the plot and get data of resulting curves
    soData = plotStandoffZ(pdfDir + '%s.pdf' % outName, dataDir + '%s.p' % outName)

    # Write data to .dat files
    exportDictToDat(soData, gpDir + '%s.dat' % outName)

    # = WITH CUT =
    # Create files for z-standoff plot
    sd.sortDataSaveRand('randPosCut.p', 'randDriftPosCut.p', nBins)

    # Create the plot and get data of resulting curves
    soDataCut = plotStandoffZ(pdfDir + '%s.pdf' % outNameCut, dataDir + '%s.p' % outNameCut)

    # Write data to .dat files
    exportDictToDat(soDataCut, gpDir + '%s.dat' % outNameCut)

    # After one run is finished, store remaining parameter
    # sets to file which can be used to resume an interrupted loop
    print 'Done!'
    print

    return outName + '.p', outNameCut + '.p'

def sendCmd(cmd):
    p = subprocess.Popen(cmd , stdout=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    return out, err

# Returns number of free cluster nodes
def getFreeNodes():
    cmd = ('pbsnodes | grep -o "= free" | wc -l')
    return int( sendCmd( cmd )[0] )
    
# Waits until qstat gives no output
def qstatWait(jobNr=None, delay=60, verbose=False):
    cmd = 'qstat'
    out, err = sendCmd( cmd )

    if verbose:
        startTime = time.time()
        print 'Started at:', time.strftime('%H:%M:%S', time.localtime(time.time()))

    while out or err:
        time.sleep( delay )
        if verbose:
            sys.stdout.write('Running for %s min\r' % int((time.time() - startTime)/60) )
            sys.stdout.flush()

        out, err = sendCmd( cmd )
        if jobNr:
            runJobs = [ item.split('.')[0] for item in out.split('\n')[2:] ][:-1]
            for job in jobNr:
                if job in runJobs:
                   break 
            else:
                break

    if verbose:
        print 'Finished at:', time.strftime('%H:%M:%S', time.localtime(time.time()))

def sendJobs(script):
    cmd = ( 'bash %s' % script )
    print cmd
    out, err = sendCmd( cmd  )

    print out

    jobNr = [ item.split('.')[0] for item in out.split('\n') ][:-1]

    # if err:
    #    return err

    return jobNr, err

def plotStandoffZ(pdfOut, dataOut):
    hMCDictZ, hDataDictZ = sd.sortDataFile('randStandoffZ.root', 'randDriftStandoffZ.root', nBins)
    
    # Get values of the curves
    # yGrid is empty! Only used vor z - theta plot
    xGrid, yGrid, valMCNorm, valDataNorm, valDiffNorm, valDiffNormError = sd.plotStandoff(hMCDictZ, hDataDictZ, nBins, 0, True, False)
   
    # Prepare data for plot
    valDiffShow = np.array( valDiffNorm )
    valDiffErrorShow = np.array( valDiffNormError )

    # Create plots
    pdfOutName = pdfOut.split('.')[0]
    peak.findPeaks(xGrid, valDiffShow[:,0], valDiffErrorShow[:,0], True, '(data - mc)/data for Bin #%s' % 0, pdfOutName + '-0.pdf', False)
    peak.findPeaks(xGrid, valDiffShow[:,1], valDiffErrorShow[:,1], True, '(data - mc)/data for Bin #%s' % 1, pdfOutName + '-1.pdf', False)

    # Arrange data in dictionary
    out = {}
    out['x'] = xGrid
    out['mc'] = valMCNorm
    out['data'] = valDataNorm
    out['diff'] = valDiffNorm
    out['diff_err'] = valDiffNormError

    # Store data in file
    cPickle.dump( out, open(dataOut, 'wb') )

    return out

def writeArtificialParam(paramDict, outFn): # e.g.: c_m=1.68588137007714, d_A=1.68588137007714
    # Get main file containing all important parameters
    inFileName = 'gp_data/artificialField_paramRoot.py'

    # Copy file and change parameters there
    outFileName = outFn # './artificialField_param.py'
   
    outFile = open(outFileName, 'w')

    keyList, valList = [], []
    for key, val in paramDict.iteritems():
        keyList.append(key)
        valList.append(val)

    for line in open(inFileName, 'r'):
        for key, value in paramDict.iteritems():
            if line.strip().startswith('%s =' % key):
                outFile.write( '%s = %.8f\n' % (key, value) )
                break
        else:
            outFile.write( line )
    
    outFile.close()

# Used in Nelder Mead
def getParamDict(parList):
    entryDict = {}
    for entry in parList:
        name, val = entry
        entryDict[name] = val

    return entryDict

# parList has the form
# [ [entry1], [entry2], ... ] 
# with each entry of the form
# entryi = ['name', start, end, num]
def getParamDictList(parList):
    parDictList = []

    nameList = []
    entryList = []
    for entry in parList:
        name, start, end, num = entry
        nameList.append( name )
        entryList.append( list( np.linspace(start, end, num) ) )

    for lst in list( itertools.product(*entryList) ):
        entryDict = {}
        for i, item in enumerate( lst ):
            entryDict[nameList[i]] = item

        parDictList.append( entryDict )

    return parDictList

# Write values in dictionary to ascii file 
def exportDictToDat(d, fName):
    f = open(fName, 'w')

    # Get x-values, which are the same for each bin
    dataDict = dict(d)
    x = dataDict['x']
    del dataDict['x']

    keyList = []
    valList = []

    for key, val in dataDict.iteritems():
        keyList.append( key )
        valList.append( val )

    f.write( '# x\t%s\n' % ('\t'.join(keyList)) )
    
    for binN in range( 5 ):
        binList = []
        for i, key in enumerate( keyList ):
            binList.append( np.array( valList[i] )[:,binN] )

        for i, xElm in enumerate( x ):
            f.write( '%f\t%s\n' % (xElm, '\t'.join(str(item) for item in np.array( binList )[:,i]) ) )

        # keyEntry = valList[i]
            # row = [] + list( np.array( keyEntry )[:,binN] )
            #f.write( '%s\n' % ('\t'.join(str(item) for item in row)) )
        f.write('\n\n')
    
    f.close()

def storeResult( f, parDict ):
    keyList = []
    valList = []

    for key, val in parDict.iteritems():
        keyList.append( key )
        valList.append( val )

    valList = [ x for (y, x) in sorted( zip(keyList, valList) ) ]
    keyList = sorted( keyList )

    for item in valList:
        f.write( '%s\t' % item )
    f.write( '\n' )

# === NELDER MEAD ===
def parDict2vec(parDict):
    keyList = []
    valList = []

    for key, val in parDict.iteritems():
        keyList.append( key )
        valList.append( val )

    valList = [ x for (y, x) in sorted( zip(keyList, valList) ) ]
    keyList = sorted( keyList )

    return keyList, valList

def vec2parDict(keyList, valList):
    d = {}
    # print keyList, valList

    for i in range(len(keyList)):
        d[keyList[i]] = list(valList)[i]

    return d
            
def oneHotVector(n, i):
    return np.array( list( np.zeros(i) ) + [1] + list( np.zeros(n-i-1) ) )

def h(x, i):
    if x[i] != 0.:
        # return 0.05
        return 0.25
    else:
        return 0.00025

# xList is a list of n vectors
def centroid(xList):
    n = len(xList)
    return 1./n * np.sum( xList, axis=0 )

def parCut(keyList, x, bounds):
    # print keyList
    # print x
    # Check if a value in x is out of bounds
    # and return 0 if so
    for key, tup in bounds.iteritems():
        val = x[keyList.index( key )]
        # print key, val, 

        low, high = tup
        # print low, high
        if val < low or val > high:
            print 'Result out of bounds, returning 0'
            return 0.

    # keyList must not contain res!
    parDict = vec2parDict(keyList, x)
    return parDict

def searchParResult(parDict, parDictList, keyList):
    for d in parDictList:
        # len+1 because result was added
        if len(d) != len(parDict)+1:
            continue

        try:
            for key in keyList:
                if not d[key] == parDict[key]:
                    break
            else:
                res = d['res'] # numpy.random.uniform() 
                print 'Found existing result: %f' % res
                print d
                return res

        except:
            continue
    else:
        return None

def getRes(keyList, x, bounds):
    import numpy.random
    parDict = parCut(keyList, x, bounds)

    # Store results in dictionaries
    try:
        parDictList = cPickle.load( open(resultDB, 'rb') )
    except:
        parDictList = []

    # Check if file was already progressed

    parRes = searchParResult(parDict, parDictList, keyList)
    if parRes:
        return -1 * parRes
    
    else:
        outName, outNameCut = evalParams(parDict)

        # res = ce.evalPara( outName, 'real/real_60.p' )[-1] 
        # resCut = ce.evalPara( outNameCut, 'real/realCut_60.p', True )[-1]
        resCut = ce.evalPara( outNameCut, 'real/realSide_60.p', True )[-1]
        res = resCut
        # res = 0.25 * (res + 3*resCut)

        parDict['res'] = res
        parDictList.append( parDict )
        cPickle.dump(parDictList, open(resultDB, 'wb'))

    return -1 * res

def getResEXOAnalysis(keyListTop, xTop, keyListBottom, xBottom, bounds, cnt=0):
    import logLikeli as ll

    parDictTop, parDictBottom = parCut(keyListTop, xTop, bounds), parCut(keyListBottom, xBottom, bounds)

    # Top
    try:
        parDictListTop = cPickle.load( open(resultDBTop, 'rb') )
    except:
        parDictListTop = []

    # Bottom
    try:
        parDictListBottom = cPickle.load( open(resultDBBottom, 'rb') )
    except:
        parDictListBottom = []

    parResTop, parResBottom = searchParResult(parDictTop, parDictListTop, keyListTop), searchParResult(parDictBottom, parDictListBottom, keyListBottom)

    if parResTop and parResBottom:
        return parResTop, parResBottom
    else:
        # Creates MC file
        print 'Create MC files...',
        mcFile = evalParamsEXOAnalysis(parDictTop, parDictBottom)
        # mcFile = '$VAULT/analysis/preprocess/pre_nelderMead_pre.root'
        print 'Done!'

        print 'Get log likelihood...',
        resTopSS, resBottomSS = ll.getLogLikelihood([mcFile], outFn='%d' % cnt, outDir=exoRootOut, art='ss')
        resTopMS, resBottomMS = ll.getLogLikelihood([mcFile], outFn='%d' % cnt, outDir=exoRootOut, art='ms')
        print 'Done!'

        resTop = 0.5 * (resTopSS + resTopMS)
        resBottom = 0.5 * (resBottomSS + resBottomMS)

        parDictTop['res'], parDictBottom['res'] = resTop, resBottom
        parDictListTop.append( parDictTop )
        parDictListBottom.append( parDictBottom )

        cPickle.dump(parDictListTop, open(resultDBTop, 'wb'))
        cPickle.dump(parDictListBottom, open(resultDBBottom, 'wb'))
        return resTop, resBottom

def logResults(vecList, resList, logFile=resultFileNM):
    with open(logFile, 'a') as f:
        for i in range(len(vecList)):
            for vec in vecList[i]:
                f.write('%f\t' % vec)
            f.write('%f\n' % (-1*resList[i]))
        f.write('\n\n')

def nelderMeadEXOAnalysis(initialGuessTop, initialGuessBottom, bounds):
    import logLikeli as ll

    with open(resultFileNMTop, 'a') as f:
        f.write('\n% ============\n')
    with open(resultFileNMBottom, 'a') as f:
        f.write('\n% ============\n')

    if len(initialGuessTop) != len(initialGuessBottom):
        print 'Need same initial parameters for top and bottom!'
        return False

    n = len( initialGuessTop )

    # Create initial points
    keyListTop, initVecsTop = getInitialVectors(initialGuessTop) 
    keyListBottom, initVecsBottom = getInitialVectors(initialGuessBottom)

    # Evaluate the n+1 vertices
    print 'Evaluate initial vectors...',
    resListTop, resListBottom = [], []
    loopCnt = 0
    for i in range(len(initVecsTop)):
        resTop, resBottom = getResEXOAnalysis(keyListTop, initVecsTop[i], keyListBottom, initVecsBottom[i], bounds, loopCnt)
        resListTop.append( resTop )
        resListBottom.append( resBottom )

        loopCnt += 1
        
    print 'Done!'

    # Create Nelder Mead generators
    nmTop, nmBottom = nelderMeadGen(initVecsTop, keyListTop, resListTop, bounds, resultFileNMTop), nelderMeadGen(initVecsBottom, keyListBottom, resListBottom, bounds, resultFileNMBottom)

    iterations = 40
    xTop, xBottom = nmTop.next(), nmBottom.next()
    while loopCnt < iterations:
        # Get current vector
        print 'Top & Bottom:'
        print xTop, xBottom

        # Evaluate
        resTop, resBottom = getResEXOAnalysis(keyListTop, xTop, keyListBottom, xBottom, bounds, loopCnt)
        print resTop, resBottom
        print

        # Send results to generators
        nmTop.next()
        xTop = nmTop.send( resTop )
        nmBottom.next()
        xBottom = nmBottom.send( resBottom )

        loopCnt += 1

def nelderMeadGen(vecList, keyList, resList, bounds, logFile=resultFileNM):
    # Dimension of the problem
    n = len( keyList )

    # Used parameters
    alpha = 1.
    beta = 1 + 2./n
    gamma = 0.75 - 1./(2*n)
    delta = 1 - 1./n
    
    # = START =
    while True:
        # 1. Sort
        # Sort results and vectors,
        # ordering is from best to worst
        vecList = [list(vec) for vec in vecList]
        vecList = [ x for (y, x) in sorted( zip(resList, vecList) ) ]
        resList = sorted( resList )
        
        # Log the results in resultFileNM
        # print vecList
        print resList
        logResults(vecList, resList, logFile)

        # 2. Computing Centroid
        # Get centroid of all vectors except worst
        x_ = centroid( vecList[:-1] )

        # 3. Transformation
        # a. Reflection
        x_r = x_ + alpha * (x_ - vecList[n])
        print
        print x_
        print x_r
        yield x_r
        res_r = yield # getRes(keyList, x_r, bounds)
        print res_r
        
        if resList[0] <= res_r and res_r < resList[n-1]:
            resList[n] = res_r
            vecList[n] = x_r
            print 'Reflection done'
            continue

        # b. Expansion
        if res_r < resList[0]:
            x_e = x_ + beta*(x_r - x_)
            yield x_e
            res_e = yield # getRes(keyList, x_e, bounds)

            if res_e < res_r:
                vecList[n] = x_e
                resList[n] = res_e
            else:
                vecList[n] = x_r
                resList[n] = res_r

            print 'Expansion done'
            continue

        # c. Outside Contraction
        if resList[n-1] <= res_r and res_r < resList[n]:
            x_oc = x_ + gamma*(x_r - x_)
            yield x_oc
            res_oc = yield # getRes(keyList, x_oc, bounds)
            
            if res_oc <= res_r:
                vecList[n] = x_oc
                resList[n] = res_oc

                print 'Outside Contraction done'
                continue
            else:
                # e. Shrink
                for i in range(1, n+1):
                    vecList[i] = list( np.array(vecList[0]) + delta*(np.array(vecList[i]) - np.array(vecList[0])) )
                    yield vecList[i]
                    resList[i] = yield # getRes(keyList, vecList[i], bounds)
                print 'Shrink done'
                continue


        # d. Inside Contraction
        if res_r >= resList[n]:
            x_ic = x_ - gamma*(x_r - x_)
            yield x_ic
            res_ic = yield # getRes(keyList, x_ic, bounds)

            if res_ic < resList[n]:
                vecList[n] = x_ic
                resList[n] = res_ic
        
                print 'Inside Contraction done'
                continue

            else:
                # e. Shrink
                for i in range(1, n+1):
                    vecList[i] = list( np.array(vecList[0]) + delta*(np.array(vecList[i]) - np.array(vecList[0])) )
                    yield vecList[i]
                    resList[i] = yield # getRes(keyList, vecList[i], bounds)
                print 'Shrink done'
                continue

def getInitialVectors(initialGuess):
    n = len(initialGuess)

    # Initial vectors
    keyList, xI, = parDict2vec( initialGuess )
    xI = np.array( xI )

    # Generate n+1 initial points,
    # where x[0] is the guessed value
    vecList = [xI]
    for i in range(n):
        # x := x_(i+1)
        x = xI + xI * h(xI, i) * oneHotVector(n, i)
        vecList.append( list(x) )

    return keyList, vecList
    
def nelderMead(initialGuess, bounds, limit=0.9, iterations=30):
    with open(resultFileNM, 'a') as f:
        f.write('\n% ============\n')

    # n is the dimension of the problem
    n = len( initialGuess )

    # Used parameters
    alpha = 1.
    beta = 1 + 2./n
    gamma = 0.75 - 1./(2*n)
    delta = 1 - 1./n

    # Initial vectors
    keyList, xI = parDict2vec( initialGuess )
    xI = np.array( xI )

    # Generate n+1 initial points,
    # where x[0] is the guessed value
    vecList = [xI]
    for i in range(n):
        # x := x_(i+1)
        x = xI + xI * h(xI, i) * oneHotVector(n, i)
        vecList.append( list(x) )
    
    resList = []

    # = START =
    # Evaluate the n+1 vertices
    print 'Evaluate initial vectors...',
    for vec in vecList:
        res = getRes(keyList, vec, bounds)
        resList.append( res )
    print 'Done!'
        
    loopCnt = 0
    while True:
        loopCnt += 1

        # 1. Sort
        # Sort results and vectors,
        # ordering is from best to worst
        vecList = [ x for (y, x) in sorted( zip(resList, vecList) ) ]
        resList = sorted( resList )
        
        # Log the results in resultFileNM
        logResults(vecList, resList)

        # Break the loop if one of the conditions is fulfilled
        if -resList[0] > limit or loopCnt > iterations:
            print 'BREAK'
            break

        # 2. Computing Centroid
        # Get centroid of all vectors except worst
        x_ = centroid( vecList[:-1] )

        # 3. Transformation
        # a. Reflection
        x_r = x_ + alpha * (x_ - vecList[n])
        print
        print x_
        print x_r
        res_r = getRes(keyList, x_r, bounds)
        print resList
        
        if resList[0] <= res_r and res_r < resList[n-1]:
            resList[n] = res_r
            vecList[n] = x_r
            print 'Reflection done'
            continue

        # b. Expansion
        if res_r < resList[0]:
            x_e = x_ + beta*(x_r - x_)
            res_e = getRes(keyList, x_e, bounds)

            if res_e < res_r:
                vecList[n] = x_e
                resList[n] = res_e
            else:
                vecList[n] = x_r
                resList[n] = res_r

            print 'Expansion done'
            continue

        # c. Outside Contraction
        if resList[n-1] <= res_r and res_r < resList[n]:
            x_oc = x_ + gamma*(x_r - x_)
            res_oc = getRes(keyList, x_oc, bounds)
            
            if res_oc <= res_r:
                vecList[n] = x_oc
                resList[n] = res_oc

                print 'Outside Contraction done'
                continue
            else:
                # e. Shrink
                for i in range(1, n+1):
                    vecList[i] = list( np.array(vecList[0]) + delta*(np.array(vecList[i]) - np.array(vecList[0])) )
                    resList[i] = getRes(keyList, vecList[i], bounds)
                continue
                print 'Shrink done'


        # d. Inside Contraction
        if res_r >= resList[n]:
            x_ic = x_ - gamma*(x_r - x_)
            res_ic = getRes(keyList, x_ic, bounds)

            if res_ic < resList[n]:
                vecList[n] = x_ic
                resList[n] = res_ic
        
                print 'Inside Contraction done'
                continue

            else:
                # e. Shrink
                for i in range(1, n+1):
                    vecList[i] = list( np.array(vecList[0]) + delta*(np.array(vecList[i]) - np.array(vecList[0])) )
                    resList[i] = getRes(keyList, vecList[i], bounds)
                print 'Shrink done'
                continue

    # Good enough solution found!
    print 'Success - Solution found'
    print '========================'
    print 'Used vectors: ', vecList
    print 'Results: ', resList

    return

# === SUPPORT FUNCTIONS ===
def makeDir(name):
    if not os.path.exists(name):
        os.makedirs(name)

def killAllJobs():
    print 'All running jobs are killed'
    cmd = ( 'qselect -u mppi025h | xargs qdel' )
    sendCmd( cmd )

# Argument parser
def get_args():
    ap = argparse.ArgumentParser(description='')

    ap.add_argument('-r', '--resume', help='Resume from last finished parameter set', required=False, action='store_true', default=False)
    ap.add_argument('-nm', '--nelder_mead', help='Use Nelder-Mead simplex algorithm to find solution', required=False, action='store_true', default=False)
    ap.add_argument('-l', '--lima', help='Run on lima cluster', required=False, action='store_true', default=False)
    ap.add_argument('-e', '--exo', help='Use EXOAnalysis', required=False, action='store_true', default=False)

    args = ap.parse_args()
    if args.exo and not args.nelder_mead:
        ap.error('Currently, only Nelder Mead is implemented for EXOAnalysis')

    return args.resume, args.nelder_mead, args.lima, args.exo

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass
        # killAllJobs()

