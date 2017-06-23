#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def main():
    fN = 'alpha_pos.dat'
    with open(fN) as f:
        content = [x.strip().split('\t') for x in f.readlines()][:-1]
    content = [[float(i) for i in j] for j in content]

    theta = np.array(content)[:,0]
    z = np.array(content)[:,1]

    # 2d histogram
    heatmap, thetaEdges, zEdges = np.histogram2d(theta, z, bins=40)
    extent = [thetaEdges[0], thetaEdges[-1], zEdges[0], zEdges[-1]]

    thetaPos = []
    thetaNeg = []
    for zVal, thetaVal in zip(z, theta):
        if zVal < 0:
            thetaNeg.append( thetaVal )
        else:
            thetaPos.append( thetaVal )

    plt.clf()
    plt.hist(thetaPos, 40, normed=0)
    plt.hist(thetaNeg, 40, normed=0, alpha=0.5)
    plt.show()
    raw_input('')

    plt.clf()
    plt.imshow(heatmap.T, origin='lower')
    plt.show()

if __name__ == '__main__':
    main()

