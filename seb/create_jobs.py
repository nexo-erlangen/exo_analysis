#/usr/bin/env python

import os
import sys
import stat
import random
from datetime import timedelta
from shutil import copyfile

TIME_PER_RUN = timedelta(hours=23, minutes=59, seconds=00)

# ~2:34 for 1000 events
# ~10:24 for 5000 events

def main():
	print "Creating multiple .exo scripts for job execution"
	print "================================================"
	fname = raw_input("Script to use: ")
	if not fname:	
		print "Need a script!"
		sys.exit()

	n = input("Number of events to generate: ")
	m = input("Number of jobs to create: ")
	outdir = raw_input("Specify output dir (blank for '.') ")
	if not outdir:
		outdir = './'
	outfile = raw_input("Name of output files (wo extension): ")

        createJobs(fname, n, m, outdir, outfile)

def createJobs(fname, n, m, outdir, outfile, lima=False):
	f = open(fname, 'r')
	events_per_job = float(n)/m

	# create out folder
	makedir(outdir)
	makedir(outdir + '/root_files')

	script = []
	for line in f:
		script.append(line)
	f.close()

	out_line = 0
	include = 0
	initial_seed = 0
	for i, line in enumerate(script):
		if 'maxevents' in line:
			script[i] = 'maxevents ' + str(int(events_per_job)) + '\n'
		if 'exosim/initial_seed' in line:
			initial_seed = i
		if 'toutput/file' in line:
			out_line = i
		if 'exosim/macro' in line:
			include = i

	if out_line == 0:
		print "Error: no output file set!"
		sys.exit()

	if include != 0 and outdir != './':
		inc_file = ((script[include].split(' ')[1]).split('/')[-1]).rstrip()
		copyfile('/home/vault/capm/mppi025h/EXOoffline/exoout_/macros/' + inc_file, outdir + '/' + inc_file)

	f = open(outdir + '/run_submit.sh', 'w')
	f.write('#!/bin/bash\n')
	for i in range(m):
		outname = outfile + '_' + str(i)
		script[out_line] = '/toutput/file ' + outdir + '/root_files/' + outname + '.root' + '\n'
		script[initial_seed] = '/exosim/initial_seed ' + str(random.randint(0, sys.maxint)) + '\n'
                script[include] = '/exosim/macro %s' % ('/home/vault/capm/mppi025h/EXOoffline/exoout_/macros/' + inc_file + '\n')
		list_to_file(script, outdir, outname)
		writeJobScript(outdir, outname, lima)
		f.write('qsub ' + outdir + '/job_scripts/submit_job_' + outname + '.sh\n')
	f.write('\n')
	os.chmod(outdir + '/run_submit.sh', stat.S_IRWXO | stat.S_IRWXG | stat.S_IRWXU)
	f.close()

	writeDatFile(outdir, outfile, m)
	writePreprocessScript(outdir, outfile)

def list_to_file(script, outdir, outname):
	f = open(outdir + '/' + outname + '.exo', 'w')
	for line in script:
		f.write('%s' % line)
	f.close()

def writeJobScript(outdir, outname, lima=False):
	global TIME_PER_RUN
	
	# Create folder for log files
	try:
		os.mkdir(outdir + '/log')
		os.mkdir(outdir + '/job_scripts')
	except OSError:
		pass

	# Write script
	fname = outdir + '/job_scripts/submit_job_' + outname + '.sh'
	f = open(fname, 'w')
	f.write('#!/bin/bash -l\n')
	f.write('#\n')
	f.write('#PBS -V -j oe\n')
        if lima:
            f.write('#PBS -l nodes=1:ppn=24,walltime=' + str(TIME_PER_RUN) + '\n')
        else:
            f.write('#PBS -l nodes=1:ppn=4,walltime=' + str(TIME_PER_RUN) + '\n')
	f.write('#PBS -N ' + outname + '\n')
	f.write('#PBS -o ' + outdir + '/log/' + outname + '.log' + '\n')
	f.write('#\n')
	f.write('cd $PBS_O_WORKDIR\n')
        if lima:
            f.write('source $VAULT/.exorc_lima\n')
        else:
            f.write('source $VAULT/.exorc\n')
	f.write('EXOAnalysis ' + outdir + '/' + outname + '.exo\n\n')
	f.close()

def writeDatFile(outdir, outname, n_jobs):
	fulldir = outdir + '/root_files'
	f = open(fulldir + '/' + outname + '.dat', 'w')
	for i in range(n_jobs):
		f.write(fulldir + '/' + outname + '_' + str(i) + '.root\n')
	f.close()

def writePreprocessScript(outdir, outname):
	fname = outdir + '/preprocess.sh'
	if not os.path.isfile(fname):
		f = open('/home/vault/capm/mppi025h/EXOoffline/exoout_/macros/preprocess.sh', 'r')
		contents = f.readlines()
		f.close()

		contents.insert(1, 'NAME=' + outname + '\n')
		contents.insert(1, 'OUTFILE=' + outname + '_pre.root' + '\n')
		contents.insert(1, 'OUTDIR=' + outdir + '\n')
		contents.insert(1, 'INFILE=' + outname + '.dat\n')
		contents.insert(1, 'INDIR=' + outdir + '/root_files\n')

		f = open(fname, 'w')
		f.writelines(contents)
		f.close()
	os.chmod(fname, stat.S_IRWXO | stat.S_IRWXG | stat.S_IRWXU)

def makedir(path):
	try:
		os.mkdir(path)
	except OSError:
		pass

if __name__ == '__main__':
	main()

