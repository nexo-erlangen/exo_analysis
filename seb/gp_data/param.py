#!/usr/bin/env python

import numpy as np

def main():
	fInName = 'param.dat'
	fOutName = 'param_out.dat'
	
	fOut = open(fOutName, 'w')
        fOut.write('# r\tm\tt\tA\tb\n')

        REFLECTORINNERRAD = 183.2356
        startPos = 0.16 # REFLECTORINNERRAD/1000. - 0.1
        endPos = REFLECTORINNERRAD/1000. - 0.0005

        r = np.linspace(startPos, endPos, 20)

	with open(fInName) as f:
		lines = [ x.strip() for x in f.readlines() ]
		for i, line in enumerate( lines ):
                        if line:
                            if line[0] == '#':
                                    continue
                                    
                            content = line.split(' ')[:-1]
                            content.insert(0, r[i-1])
                            print content
                            for entry in content:
                                fOut.write('%f\t' % float( entry ))
                            fOut.write('\n')

if __name__ == '__main__':
	main()

