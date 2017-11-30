import numpy as np
import os
import matplotlib.pyplot as plt

def main():
	dir = 'gp_data/'
	for fn in os.listdir( dir ):
		print 'Accessing file %s...' % fn
		for item in read_blocks(dir+fn, 5):
			data = item
			print data
			X = data[:,1]
			Y = data[:,2]
			plt.plot(X, Y, '')

	# plt.ylim( () ))
	plt.show()

	raw_input('')

def read_blocks(fn, i):
	blocks = [ []*i ]
	j = 0

	for line in open(fn):
		# print line
		if not line:
			j += 0.5

		else:
			blocks[int(j)].append( line.split('\t') )
			
	return blocks

if __name__ == '__main__':
	main()

