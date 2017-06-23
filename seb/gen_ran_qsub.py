#!/usr/bin/env python

import argparse
import os
import importlib
import re
import cPickle

import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')
ROOT.PyConfig.IgnoreCommandLineOptions = True

import generate_random as gr
import plot_support as ps

def main():
	settings, outfile, generate, N, Nj, getstat, gethist = get_args()
        loadSettings(settings)
        
        if generate:
            execGen(settings, outfile, N, Nj)
        elif getstat:
            getStat()
        elif gethist:
            f = ROOT.TFile.Open('standoffDrift.root', 'RECREATE')

            ApoRand = [162./1000, 5./1000, 182./1000]
            getHistFiducial(f, 'randPos.p', H_BINS, ApoRand, 'rand')
            getHistFiducial(f, 'randDriftPos.p', H_BINS, ApoRand, 'randDrift')

            f.Close()
	else:
	    genRand(N, outfile)

def execGen(settings, outfile, N, Nj, lima=False):
	loadSettings(settings)

        if lima:
            outdir = 'drift_data_lima/'
        else:
            outdir = 'drift_data/'
	scriptDir = outdir + 'drift_scripts/'

        f = open('ran_qsub.sh', 'w')
        f.write('#!/bin/bash\n\n')

        fileNumbers = []
        if not os.path.exists(scriptDir):
                os.makedirs(outdir)
                os.makedirs(scriptDir)
                startIdx = 0

        else:
                for file in os.listdir(outdir):
                        if file.endswith('.root'):
                                fileNumbers.append( int( re.findall(r'\d+', file)[0] ) )

                if fileNumbers:
                        startIdx = max(fileNumbers)
                else:
                        startIdx = 0

        for i in range(startIdx, startIdx + Nj):
                outf = outfile.split('.')
                scriptfile = 'gen_rand_%d.sh' % i
                genRandScript(scriptDir + scriptfile, settings, N, Nj, outdir + outf[0] + '_%d.root' % i, lima)
                f.write('qsub %s\n' % (scriptDir + scriptfile))
        f.write('\n')

        f.close()

def genRandScript(scriptfile, settings, N, Nj, outfile, lima=False):
	f = open(scriptfile, 'w')

	# Write header
	f.write('#!/bin/bash -l\n')
	f.write('#\n')
	f.write('#PBS -V -j oe\n')
        if lima:
            f.write('#PBS -l nodes=1:ppn=24,walltime=03:00:00\n')
        else:
	    f.write('#PBS -l nodes=1:ppn=4,walltime=03:00:00\n')
	f.write('#PBS -N qxu_comp\n')
	f.write('#PBS -o /home/vault/capm/mppi025h/plot_scripts/log/gen_ran_qsub.log\n')
	f.write('#\n')
	f.write('cd $PBS_O_WORKDIR\n')

        # if run on lima, source specific settings
        if lima:
            f.write('source $VAULT/.exorc_lima\n')

	# Write command
	f.write('python /home/vault/capm/mppi025h/plot_scripts/gen_ran_qsub.py -s %s -N %d -J %d -o %s' % (settings, N, Nj, outfile))

	f.close()

def genRand(N, outfile):
	ApoRand = [162.e-3, 5.e-3, 182.e-3]

	# if ms:
	# 	f = ROOT.TFile(MSoutDir + MSoutFile)
	# else:

	f = ROOT.TFile(SSoutDir + SSoutFile)
	mcHist, dataHist, diffHist, diffRelHist = ps.getHist(f) 

        '''
        f = ROOT.TFile('pcdDistro.root', 'READ')
        mcHist = f.Get('h')
        '''

	ROOT.gROOT.cd()
	randTree = ROOT.TTree('random', 'random')
	dumpTree = ROOT.TTree('dump', 'dump')

	print 'Random positions:'
	randPos = gr.generateRandom(mcHist, dumpTree, N, 'side')
	# randPos = gr.generateRandom(mcHist, randTree, N)
		
	print 'Random drifted positions:'
	gr.generateRandom(mcHist, randTree, N, 'side')
	del dumpTree

	name = 'StandoffDrift'
	inFile = None		# only when used with potential

	efieldName = '/home/vault/capm/mppi025h/g/exo_data6/groups/3DFieldMaps/3Dmaxwell/3D_Efield.root'
	# Set Bins to None, so H_BINS are read from efield-file 
	H_BINS_RAND = None

	# Set to use completeDrift
	gr.randomProcess(randPos, randTree, H_BINS_RAND, H_BINS, ApoRand, name, N, outfile, True, efieldName, True, inFile, True, False)

	f.Close()

def loadSettings(settings):
        import sys
        sys.path.append('settings/')

	settingsLib = importlib.import_module(settings.split('.')[0])
	settingsMod = settingsLib.__dict__
	try:
		to_import = settingsLib.__all__
	except AttributeError:
		to_import = [name for name in settingsMod if not name.startswith('_')]
	globals().update({name: settingsMod[name] for name in to_import})

def getStat(merge=True, lima=False):
        if lima:
            outdir = 'drift_data_lima'
        else:
            outdir = 'drift_data'

	if merge:
		ps.mergeLists( ps.getFileList(outdir, 'randPos_'), 'randPos.p', True )
		ps.mergeLists( ps.getFileList(outdir, 'randDriftPos_'), 'randDriftPos.p', True )
		ps.mergeLists( ps.getFileList(outdir, 'randPosCut_'), 'randPosCut.p', True )
		ps.mergeLists( ps.getFileList(outdir, 'randDriftPosCut_'), 'randDriftPosCut.p', True )

	# Load positions from file
	randPosList = cPickle.load(open('randPos.p', 'rb'))
	randDriftPosList = cPickle.load(open('randDriftPos.p', 'rb'))
	randPosCutList = cPickle.load(open('randPosCut.p', 'rb'))
	randDriftPosCutList = cPickle.load(open('randDriftPosCut.p', 'rb'))
	print 'Opened position files'
	
	# out = ROOT.TFile.Open('standoff.root', 'RECREATE')
	outDrift = ROOT.TFile.Open('standoffDriftData.root', 'RECREATE')

	# Create tree structure
	t = ROOT.TTree('mcTree', 'mcTree')
	tDrift = ROOT.TTree('dataTree', 'dataTree')

	# Fill trees with positions
	listToRoot(randPosList, t)
	listToRoot(randDriftPosList, tDrift)

	# Store to TFile
	outDrift.Write()
	outDrift.Close()

        # Same for cut data
        outDriftCut = ROOT.TFile.Open('standoffDriftDataCut.root', 'RECREATE')
	tCut = ROOT.TTree('mcTree', 'mcTree')
	tDriftCut = ROOT.TTree('dataTree', 'dataTree')
	listToRoot(randPosCutList, tCut)
	listToRoot(randDriftPosCutList, tDriftCut)
        outDriftCut.Write()
        outDriftCut.Close()

def getHistFiducial(f, fPosName, H_BINS, Apo, name):
	posList = cPickle.load(open(fPosName, 'rb'))

	fiducialVerts = gr.isFiducialList(posList, *Apo)
	h = gr.fillRandHist(posList, H_BINS, '(%sHist)' % name)
	h1 = gr.fillRandHist(fiducialVerts, H_BINS, '(%sHistCut)' % name)
	so = gr.getStandoffList(fiducialVerts, '(%sStandoffCut)' % name)

	h.Scale(1./h.Integral())
        h1.Scale(1./h1.Integral())
	so.Scale(1./so.Integral())
 
	h.Write()
        h1.Write()
	so.Write()

def listToRoot(l, t):
	es = ROOT.EXOEventSummary()
	t.Branch('EventSummary', es)

	for pos in l:
		tempVecList = [ ROOT.std.vector('double')(1), ROOT.std.vector('double')(1), ROOT.std.vector('double')(1) ]
		
		for i, coord in enumerate( pos ):
			tempVecList[i][0] = coord * 1000

		# print tempVecList[0][0], tempVecList[1][0], tempVecList[2][0]
		es.cluster_x = tempVecList[0]
		es.cluster_y = tempVecList[1]
		es.cluster_z = tempVecList[2]

		t.Fill()

def get_args():
	ap = argparse.ArgumentParser(description=' ')

	ap.add_argument('-s', '--settings', type=str, help='Settings to load', required=True)
	ap.add_argument('-o', '--outfile', type=str, help='name of output-file', required=True)
	ap.add_argument('-N', '--number', type=int, help='Number of events per job', required=False)
	ap.add_argument('-J', '--jobs', type=int, help='Number of jobs', required=False)
	ap.add_argument('-g', '--generate', action='store_true', default=False)
	ap.add_argument('-gs', '--getstat', action='store_true', default=False)
	ap.add_argument('-gh', '--gethist', action='store_true', default=False)

	args = ap.parse_args()
	if not args.getstat and not args.gethist and not args.number and not args.jobs:
		ap.error('Number of jobs and events per job needed')

	return args.settings, args.outfile, args.generate, args.number, args.jobs, args.getstat, args.gethist

if __name__ == '__main__':
	main()

