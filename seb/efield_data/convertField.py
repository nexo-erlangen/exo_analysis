#!/usr/bin/env python
import numpy as np
import csv

def main():
	fInName = 'export2D_AvgPot.txt'
	fOutName = 'export2D_AvgPotConverted.csv'
	convert(fInName, fOutName)

def convert(fInName, fOutName):
	fOut = open(fOutName, 'w')

	with open(fInName, 'rb') as f:
		for line in f:
			if line.startswith('%'):
				continue
			row = list( np.nan_to_num( [float(item) for item in (line.strip()).split()] ) )
			fOut.write('%f,%f,%f,%f,%f,%f\n' % (0, row[0], row[1], 0, row[2], row[3]))
			
	fOut.close()

def convertCSV(fInName, fOutName):
	fOut = open(fOutName, 'w')

	with open(fInName, 'rb') as f:
		reader = csv.reader(f)
		for row in reader:
			if row[0].startswith('#'):
				continue

			row = list( np.nan_to_num( [ float(item) for item in row ] ) )
			fOut.write('%f,%f,%f,%f,%f,%f\n' % (0, row[0], -row[2], row[4], row[3], -row[5]))

	fOut.close()

if __name__ == '__main__':
	main()

