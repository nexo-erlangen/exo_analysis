import numpy as np 
import cPickle
import sys
sys.path.insert(0, '..')

import csv_to_vec as ctv

def main():
	posList, vecList = ctv.storeEFieldToList('./', 'Efield_yz-plane_spreadsheet.csv')

	# Correct coordinates
	for pos in posList:
		pos[2] -= 0.203796
	
	newPosList = []
	newVecList = []
	for i, pos in enumerate( posList ):
		x, y, z = pos

		if z < getFieldRingPos(7) and z > getFieldRingPos(8):
			newPosList.append( pos )
			newVecList.append( vecList[i] )

	zList = np.array( newPosList )[:,2]
	zDiff = max(zList) - min(zList)
	zMin = -min(zList)

	print zDiff, zMin
	remain = (-min(zList) % zDiff)
	times = int(-min(zList)/zDiff)
	print times, remain

	finalPosList = []
	finalVecList = []
	
	for i, pos in enumerate( posList ):
		x, y, z = pos
		if (z < 0 and z > -remain):
			finalPosList.append( [x, y, z] )
			finalVecList.append( vecList[i] )

		elif (z < -zMin):
			finalPosList.append( [-x, -y, z] )
			finalVecList.append( vecList[i] )

	for i in range(times):
		for j, pos in enumerate( newPosList ):
			x, y, z = pos
			z += i * zDiff

			finalPosList.append( [x, y, z] )
			finalVecList.append( newVecList[j] )

	#for i in range( len(finalPosList) ):
	#	print finalPosList[i], finalVecList[i]

	print len(newPosList), len(newPosList)
	print len(finalPosList), len(finalVecList)

	storeListToFile(finalPosList, 'posListTest.p')
	storeListToFile(finalVecList, 'vecListTest.p')

def storeListToFile(ls, fName):
	cPickle.dump(ls, open(fName, 'wb'))
	
def getFieldRingPos(n):
	return -0.0182753 - n*0.0168656

if __name__ == '__main__':
	main()

