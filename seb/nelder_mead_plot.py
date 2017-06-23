#!/usr/bin/env python
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

import cPickle
from operator import itemgetter

def main():
    # plotElement('c_m')
    # plt.show()
    plotNM()
    return

def plotNM():
    d = cPickle.load( open('resultDB.p', 'rb') )
    dCut = []
    for item in d:
        if len(item) == 4:
            dCut.append( item )

    keyList = dCut[0].keys()
    keyList.remove('res')
    keyList.append('res')
    print 'Using keys: ', keyList

    # Get tetraheder
    tetList = []
    for item in dCut:
        tet = []
        for key in keyList:
            tet.append(item[key])
        tetList.append( tet )

    print tetList

    # Create plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('$' + keyList[0] + '$')
    ax.set_ylabel('$' + keyList[1] + '$')
    ax.set_zlabel('$' + keyList[2] + '$')

    vertList = []
    tet = tetList[0:4]
    for i in range(4, len(tetList)-4):
        x = np.array( tet )[:,0]
        y = np.array( tet )[:,1]
        z = np.array( tet )[:,2]
        print x
        print y
        print z

        X = [x[j] for j in (0, 1, 0, 2, 0, 3, 1, 2, 3)]
        Y = [y[j] for j in (0, 1, 0, 2, 0, 3, 1, 2, 3)]
        Z = [z[j] for j in (0, 1, 0, 2, 0, 3, 1, 2, 3)]

        vertList.append( (X, Y, Z) )
        tet = sorted(tet, key=itemgetter(3), reverse=True)[0:3] + [tetList[i]]
        print tet, tetList[i]

    for i, item in enumerate(vertList):
    # for item in vertList[0:3]:
        X, Y, Z = item
        ax.scatter(X, Y, Z)
        ax.plot3D(X, Y, Z)
        plt.savefig('nm_plot/nm%d.pdf' % i)
        plt.pause(0.1)
    plt.show()

def plotElement(elm='c_A'):
    dList = cPickle.load( open('resultDB.p', 'rb') )

    valList = []
    for d in dList:
        res = d['res']
        val = d[elm]
        valList.append( (val, res) )

    valList = sorted(valList, key=lambda x:x[1])

    plt.plot(np.array(valList)[:,1], np.array(valList)[:,0])

if __name__ == '__main__':
    main()

