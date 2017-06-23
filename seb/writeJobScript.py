from datetime import timedelta
import os

# Writes a job script to run a specified python script
def main():
	outdir = '/home/vault/capm/mppi025h/plot_scripts'
	outname = 'qxu_comp'
	writeJobScript(outdir, outname, '/home/vault/capm/mppi025h/plot_scripts/saveHistograms.py')

def writeJobScript(outdir, outname, scriptName, h=12, m=0, s=0):
	TIME_PER_RUN = timedelta(hours=h, minutes=m, seconds=s)
	
	# Create folder for log files
	try:
		os.mkdir(outdir + '/log')
		# os.mkdir(outdir + '/job_scripts')
	except OSError:
		pass

	# Write script
	fname = outdir + '/submit_job_' + outname + '.sh'
	f = open(fname, 'w')
	f.write('#!/bin/bash -l\n')
	f.write('#\n')
	f.write('#PBS -V -j oe\n')
	f.write('#PBS -l nodes=1:ppn=4,walltime=' + str(TIME_PER_RUN) + '\n')
	f.write('#PBS -N ' + outname + '\n')
	f.write('#PBS -o ' + outdir + '/log/' + outname + '.log' + '\n')
	f.write('#\n')
	f.write('cd $PBS_O_WORKDIR\n')
	f.write('python ' + scriptName + '\n\n')
	f.close()

if __name__ == '__main__':
	main()

