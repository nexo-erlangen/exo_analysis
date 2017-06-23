#!/usr/bin/env python
import numpy as np
import generate_random as gr
import importlib
try:
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
except:
	pass

from scipy.optimize import fsolve
import artificialField_paramTest as par
import plot_support as ps

# Bulge in PTFE at resistor chain
reflectorRad1 = par.REFLECTORINNERRAD
reflectorRad2 = 0.04
reflectorX2 = 0.215

def main():
        # getCharge()
        # return

	# drawVectorField(*generateEField([0.14, 0.1832, 500], [-0.003, -0.19, 500]))
	# return

        '''
        theta = np.linspace(-1., 1., 5000)*np.pi
        r = [reflectorRad(t, reflectorRad1, reflectorRad2, reflectorX2) for t in theta]

        ax = plt.subplot(111, projection='polar')
        ax.plot(theta, r)
        ax.grid(True)
        plt.show()
        return
        '''

	fName = 'gp_data/artDriftLin2.dat'
	fN = open(fName, 'w')
	Nr = 40
	Nz = 10

	lossCnt = 0
	# for startX in np.linspace(par.REFLECTORINNERRAD - 0.1, par.REFLECTORINNERRAD - 0.0005, Nr):
	for startX in np.linspace(0., 0., 1):
		for startY in np.linspace(0.17, par.REFLECTORINNERRAD - 0.0005, Nr):
			for startZ in np.linspace(-0.003, -0.02, Nz):
				if np.sqrt(startX**2 + startY**2) >= par.REFLECTORINNERRAD:
					continue
				startPos = [startX, startY, startZ]
				# out = artificialDrift_C(0., startR, startZ, 250.) 
				out = artificialDrift(startPos, 0, None, False, True)
				l, hl = driftLength(startPos[1], startPos[2])
                                if round(l, 3) != round(hl, 3):
                                    print 'l = %f, hl = %f' % (l, hl)

				'''
				x, y, z, t = out
				print x, y, z
				'''

				if out:
					for item in out:
						x, y, z = item
						fN.write( '%.8f\t%.8f\t%.8f\n' % (x, y, z) )
					fN.write('\n\n')
				else:
					lossCnt += 1
	fN.close()
	print 'Percent lost:', float( lossCnt )/(Nr * Nz)

	for n in reversed(range(-10, 10)):    
		r_ = getR(0.18, -0.038)
		print n, rad(n, r_), zed(n, r_)

	drawDrift( fName, False )

def artificialTest():
    zStart = -0.01
    startPos = 0.175 # REFLECTORINNERRAD - 0.1
    endPos = par.REFLECTORINNERRAD - 0.0005

    fN = open('gp_data/artTest.dat', 'w')
    for r in np.linspace(startPos, endPos, 20):
        for zStart in np.linspace(-0.003, -0.17, 5):
            rNew = getR(r, zStart)
            # print r, rNew, f(rNew, zStart)

            '''
            isInt = fIntersected(rNew)
            print isInt
            
            if isInt and (isInt < zStart):
                zEnd = isInt
            else:
                zEnd = -0.18
            '''

            zEnd = -0.18
	    res = rNew
            for z in np.linspace(zStart, zEnd, 1000):
                res = f(rNew, z)
                #if res < REFLECTORINNERRAD:
                fN.write( '%.8f\t%.8f\t%.8f\n' % (r, z, res) )
                #else:
                #    break
		# r = getR(res, z)

            fN.write('\n')
        fN.write('\n')

    fN.close()

def artificialDrift_C(x, y, z, t):
    # print 'Initial Position:', (x, y, z)
    startPos = [coord/1000. for coord in [x, y, z]]
    r = np.sqrt( startPos[0]**2 + startPos[1]**2 )
    if r >= par.REFLECTORINNERRAD:
        x_, y_, z_ = 999, 999, 999
        return (x_, y_, z_, t)

    posHist = artificialDrift(startPos, 0, None, False, False)

    if posHist:
        res = posHist[-1]

        l, lh = driftLength(r, startPos[2])
        # print 'l = %d, lh = %d' % (l, lh)
        x_, y_, z_ = list( [coord*1000. for coord in res] )
        if z > 0:
            z_ = ( z/1.e3 - abs( lh - l ) ) * 1000
        else:
            z_ = ( z/1.e3 + abs( lh - l ) ) * 1000

    else:
        x_, y_, z_ = 999, 999, 999

    # print 'Final position:', (x_, y_, z_)

    # return x_, y_, z_, t
    return (x_, y_, z_, t)

def artificialDrift(startPos, eventNum=0, out_q=None, diffuse=False, saveHistory=False):
    # Important constants
    x_CAD = 0.19223657
    r_RID = 0.1832356
    limit = 0.187

    # Cartesian and polar vector
    xVec = list( startPos )
    x, y, z = xVec
    theta = np.arctan2(y, x)
    r = np.sqrt( x**2 + y**2 )

    # Flip if in wrong section of volume
    upFlag = False
    if z > 0:
        upFlag = True
        z *= -1

    # Get correct radius
    r_ = getR(r, z)
    xVecPol = np.array( [r, z] )

    posHistory = totalIntersect(r, r_, z, theta)
    if not posHistory:
        return False

    if saveHistory:
	    # Number of steps
	    N = 1000

	    for i, z in enumerate( np.linspace(z, -limit, N) ):
                res = total(r, r_, z)
		
		if not N % 20:
		    xVec = list( gr.polToCart( res, z, theta ) )
		    posHistory.append( xVec )

		'''
		if res > r_RID:
		    posHistory = False
		    print 'Intersect'
		    break
		'''

	    else:
                fPos = list( gr.polToCart(res, -limit, theta) )
                #if ps.isFiducial(fPos[0], fPos[1], fPos[2], 0.174, 0., 200.):
                posHistory.append( fPos )
                #else:
                #    posHistory = False

    else:
	    res = h(r_, -limit)
            fPos = list( gr.polToCart(res, -limit, theta) )
            #if ps.isFiducial(fPos[0], fPos[1], fPos[2], 0.174, 0., x_CAD):
            posHistory.append( fPos )
            #else:
            #    posHistory = False

    if posHistory and upFlag:
        posHistory = np.array( posHistory )
        posHistory[:, 2] *= -1
        posHistory = posHistory.tolist()

    if out_q:
        outdict = {}
        outdict[eventNum] = posHistory
        out_q.put(outdict)
    else:
        return posHistory
    
# === REFLECTOR RADIUS ===
def reflectorRad(theta, r1, r2, x2):
    # Intersection points of circles 1 and 2
    x_ = (r1**2 - r2**2 + x2**2) / (2*x2)
    y_ = np.sqrt( r1**2 - x_**2 )
    theta_ = np.arctan2(y_, x_)

    if ((theta > 0 and theta <=np.pi) and theta >= theta_) or ((theta <= 0 and theta >= -np.pi) and theta <= -theta_):
        return par.REFLECTORINNERRAD
    else:
        return (x2 - np.sqrt(-(x2 * np.tan(theta))**2 + (r2 * np.tan(theta))**2 + r2**2 )) / (np.tan(theta)**2 + 1)

# === FUNCTION DEFINITION ===
# = COMBINED =
def total(r, r_, z):
    if par.z2 <= z and z <= par.z1:
        res = f(r_, z)
    elif z > par.z1:
        res = f(r_, z) #g(r, r_, z)
    elif z < par.z2:
        res = h(r_, z)

    return res

def totalIntersect(r, r_, z, theta):
    R = reflectorRad(theta, reflectorRad1, reflectorRad2, reflectorX2)

    # Check for intersection
    if par.z2 < z and z < par.z1:
	interRes = fIntersect(r_, z, R)
    elif z >= par.z1:
	interRes = gIntersect(r, r_, z, R)
    elif z <= par.z2:
	interRes = hIntersect(r, r_, z, R)

    if interRes:
        posHistory = False
	return posHistory
    else:
	posHistory = [ list( gr.polToCart( r, z, theta ) ) ]
        if par.z2 <= z and z <= par.z1:
            posHistory += [ list( gr.polToCart( f(r_, z), z, theta ) ) ]
	elif z > par.z1:
            # posHistory = [ list( gr.polToCart( g(r, r_, z), z, theta ) ) ]
	    posHistory += [ list( gr.polToCart( f(r_, z), z, theta ) ) ]
        elif z < par.z2:
            posHistory += [ list( gr.polToCart( h(r_, z), z, theta ) ) ]
        return posHistory

def total_der(r, r_, z):
	if par.z2 < z and z < par.z1:
		res = f_der(r_, z)
	elif z >= par.z1:
		res = f_der(r_, z) # False
	elif z <= par.z2:
		res = h_der(r_, z)
	return res

# = F-SECTION = for z2 < z < z1
def f(r, z): 
    # return m(r) * z + A(r) * S(z) + t(r)
    return m(r) * z + A(r) * S(z) + t(r)

def m(r):
    # return np.exp( c_m * (r - d_m) )
    # return 1 + c_m * (r - d_m) + (c_m * (r - d_m))**2

    # limit = (par.d_m2 - par.d_m) / (par.c_m - par.c_m2)

    if r < par.s_m:
        return par.c_m2 * r + par.d_m2
    else:
        return par.c_m * r + par.d_m

def t(r):
    # return m_t * r + t_t
    return par.m_t * r

def A(r):
    # return -np.exp( c_A * (r - d_A) )
    # return -( 1 + c_A * (r - d_A) + (c_A * (r - d_A))**2 )

    if r < par.d_A:
        return par.c_A2 * (r - par.d_A) + par.d_A2
        # return par.d_A2 / par.d_A * r
    else:
        return -par.c_A * (r - par.d_A)**2 + par.d_A2

def S(z):
	#return np.sin( -(z - z_0)/b - np.pi )
	return hull(z) * np.sin( 2*np.pi * (z - par.z_0)/par.b )

def S_der(z):
	return 2*np.pi/par.b * np.cos(2*np.pi*(z - par.z_0)/par.b) * hull(z) + np.sin(2*np.pi*(z-par.z_0)/par.b)*(par.t_h - 1)/par.h_z0

def hull(z):
	return (par.t_h - 1)/par.h_z0 * (z - par.h_z0) + par.t_h

def f_der(r, z):
	return m(r) + A(r) * S_der(z)

def fIntersect(r, z, R):
    if par.c_m >= 0:
        t_ = (par.t_h - 1) / par.h_z0
        z_ = (R + A(r) - t(r)) / (m(r) - A(r)*t_)
        n_ = (z - par.z_0)/par.b + 0.25

        res = f(r, par.b*(np.floor(n_) - 0.25) + par.z_0)
        print res, np.floor(n_)
        if res >= R:
            return True
        else:
            return False
    else:
        n = int( np.floor( (-0.187 - par.z_0)/par.b + 0.25 ) )
        for ni in range(n, 0):
            res = f(r, par.b*(np.floor(ni) - 0.25) + par.z_0)
            if res >= R:
                return True
        else:
            return False

    '''
    n = (R + A(r) - par.m_t * r) / (2*np.pi*par.b * m(r)) - 0.75
    #z_ = (R + A(r) - m_t * r) / m(r)
    # z_ = (R + A(r)*np.sin(np.arccos(b*m(r)/A(r))) - m_t*r) / m(r)
    z_ = (R - par.m_t*r + A(r)*np.sqrt( 1 - ( (par.b*m(r)) / (2*np.pi*A(r)) )**2 )) / m(r)
    
    if z_ <= 0: # and z_ >= -0.18: 
	    if z >= z_:
		res = True
	    else:
		res = False
    else:
	    res = False

    # print z, n, z_, res
    return res
    '''

# = G-SECTION = for z >= z1
def g(r, r_, z):
    c = f(r_, par.z1)
    return C(r) * np.exp(-(z + gShift(r))**2/(2*sigma(r)**2)) + c

def sigma(r):
    return r * par.s_g 
    # return 1e-5 * 1 / (REFLECTORINNERRAD - r)**(1)

def C(r):
    lim = par.REFLECTORINNERRAD - (-par.c_C / par.a_C)**( 1./par.b_C )
    if r >= lim:
        return par.a_C * (-r+par.REFLECTORINNERRAD)**(par.b_C) + par.c_C
    else:
        return 0

def gShift(r):
    # return 0.019/REFLECTORINNERRAD * r
    return par.A_gs * np.exp(5 * (r - par.REFLECTORINNERRAD))

def gIntersect(r, r_, z, R):
    c = f(r_, par.z1)
    # if C(r) + c >= par.REFLECTORINNERRAD and z >= -gShift(r): # and z >= z1:
    if g(r, r_, z) >= R:
	    return True
    else:
	    return False
    
    '''
    if ( C(r) < 0 ) or ( r_ >= REFLECTORINNERRAD ):
        return False
    else:
	tmp = np.log( C(r) / (REFLECTORINNERRAD - c) )
	if tmp < 0:
		return False

    z_ = -sigma(r) * np.sqrt( 2*np.log(C(r)/(REFLECTORINNERRAD - c)) ) + z1
    print np.log( C(r) / (REFLECTORINNERRAD - c) ), z_,

    if z >= z_:
	    print True
	    return True
    else:
	    print False
	    return False
    '''

# = H-SECTION = for z >= -0.18
def a(r):
	return ( m(r) + A(r)*S_der(par.z2) ) / ( 2*(par.z2 - par.z0) )

def off(r):
	return -a(r) * (par.z2 - par.z0)**2 + f(r, par.z2)

def B(r):
    m = 0.966359
    t = 0.00446672
    return m*r + t

def h(r, z):
    # off = -a(r)*(par.z2 - par.z0)**2 + f(r, par.z2)
    # return a(r) * (z - par.z0)**2 + off

    '''
    limit = (par.d_m2 - par.d_m) / (par.c_m - par.c_m2)
    if r < limit:
        c_m = par.c_m2
        d_m = par.d_m2
    else:
        c_m = par.c_m
        d_m = par.d_m

    if r <= par.d_A:
        a = ( c_m * r + d_m + ( par.c_A2 * (r - par.d_A + par.d_A2) ) * S_der(par.z2)) / ( 2*(par.z2 - par.z0) )
    else:
        a = (c_m * r + d_m + ( par.d_A2 - par.c_A*(r - par.d_A)**2 ) * S_der(par.z2) ) / ( 2*(par.z2 - par.z0) )
    '''
    return a(r) * (z - par.z0)**2 + off(r)

def h_der(r, z):
	return 2 * a(r) * (z - par.z0)	

def hIntersect(r, r_, z, R):
    if r_ > R:
        return True
    return False

# === getR ===
def getR(r_x, z_x):
    # (r + z*(c_m*d_m - 1) + S(z) * (1 - c_A*d_A) - t_t) / (c_m * z - c_A * S(z) + m_t)
    # return r / m_t

    '''
    func = lambda r : -np.exp(c_A * (r - d_A)) * S(z_x) + np.exp(c_m * (r - d_A))*z_x + m_t*r - r_x
    r_init = 0.16
    r_sol = fsolve(func, r_init)

    return r_sol[0]
    '''

    b = -( (par.c_m * z_x + par.m_t) / (par.c_A * S(z_x)) + 2*par.d_A )
    c = (r_x - par.d_m * z_x) / (par.c_A * S(z_x)) - par.d_A2 / par.c_A + par.d_A**2

    r1 = ( -b + np.sqrt( b**2 - 4*c ) ) / 2
    r2 = ( -b - np.sqrt( b**2 - 4*c ) ) / 2

    if r1 < 0. or r1 > 0.2:
            res_r = r2
    else:
            res_r = r1

    if res_r > par.d_A:
        return res_r
    else:
	return (r_x + par.c_A2 * par.d_A * S(z_x) - par.d_m2 * z_x - par.d_A2 * S(z_x)) / (par.c_m2 * z_x + par.c_A2 * S(z_x) + par.m_t)

def rad(n, r):
    z = zed(n, r)
    return m(r)*z - A(r) + t(r), f(r, zed(n, r))

def zed(n, r):
    #return 2*np.pi*b*(n - 0.25)
    #return b*(2*np.pi*n - np.arccos(-b*m(r)/A(r)))
    return par.z_0 - (par.b*np.arccos(-par.b*m(r)/(2*np.pi*A(r))))/(2*np.pi) + par.b*n

# =======================
def driftLength(r, z):
	rStart = float( r )
	zStart = -float( abs( z ) )

	r_ = getR(r, zStart)
	N = int( 10.e2 )
	z = np.linspace(zStart, -0.187, N)

	r = [ total(rStart, r_, zCoord) for zCoord in z ]
	l = 0.
	for i in range( len(r) - 1 ):
        	l += np.sqrt( (r[i] - r[i+1])**2 + (z[i] - z[i+1])**2 )

	dl = np.sqrt( (r[0] - r[-1])**2 + (z[0] - z[-1])**2 )
	hl = abs( z[0] - z[-1] )

	# print 'Total length =', l
	# print 'Direct length =', dl
	# print 'Height length =', hl

	return l, hl

# === Charge Distribution ===
def getCharge():
        import csv_to_vec as ctv

        # H_BINS = [1, 0, 1, 2000, -0.22724, 0.22724, 1000, -0.21788400000000002, 0.014380000000000004]

        # Choose H_BINS that posListSim = posListArt
        posListSim, vecListSim = readEField()
        vecListSim = np.array( [np.array(vec)/np.linalg.norm(vec) for vec in vecListSim] )

        H_BINS = ctv.getH_BINS(posListSim, vecListSim)
        print H_BINS
        posListArt, vecListArt = getEFieldValues(H_BINS)

        posList = posListSim
        vecList = vecListSim - vecListArt

        for i in range(1000, 1100):
            print posListSim[i], vecListSim[i], posListArt[i], vecListArt[i], (posList[i], vecList[i])

        ctv.plotEField(posListSim, abs(vecListSim), H_BINS, 1.e-5) 
        ctv.plotEField(posListArt, abs(vecListArt), H_BINS, 1.e-5) 
        ctv.plotEField(posList, abs(vecList), H_BINS, 1.e-5) 

# === EField Plot ===
def generateEField(rRange, zRange):
	r = np.linspace( *rRange )
	z = np.linspace( *zRange )

	R, Z = np.meshgrid(r, z)
	res = np.array( [getEField(r, z) for (r, z) in zip(np.ravel(R), np.ravel(Z))] )
	resR = res[:,0]
	resZ = res[:,1]

	RESR, RESZ = resR.reshape(R.shape), resZ.reshape(Z.shape)

	return R, Z, RESR, RESZ

def drawVectorField(r, z, vecR, vecZ):
	# plt.figure()
	# Q = plt.quiver(R, Z, RESR, RESZ, units='inches')
	f, ax = plt.subplots()
	strm = ax.streamplot(r, z, vecR, vecZ, color=vecR, linewidth=1, density=2, arrowstyle='->', cmap=plt.cm.autumn)
	f.colorbar(strm.lines)

	plt.show()
	return

def readEField():
        import cPickle
	posList = cPickle.load(open('efield_data/posList.p', 'rb'))
	vecList = cPickle.load(open('efield_data/vecList.p', 'rb'))

	return np.array(posList), np.array(vecList)

def getEField(r, z):
	r_ = getR(r, z)
	rDer = np.nan_to_num( -total_der(r, r_, z) )
	zDer = -1

	if rDer:
		mag = 1./np.sqrt(rDer**2 + zDer**2)
		return rDer / mag, zDer / mag
		# Evec = 1./np.sqrt(rDer**2 + zDer**2) * np.array([rDer, zDer])
	else:
		return 0., 0.

def getEFieldValues(H_BINS):
        rRange = H_BINS[4:6] + [H_BINS[3]]
        zRange = H_BINS[7:9] + [H_BINS[6]]
	R, Z = np.linspace(*rRange), np.linspace(*zRange)

        RR = np.array( list(R) * len(Z) )
        ZZ = np.array( [inner for outer in [len(R)*[z] for z in list(Z)] for inner in outer] )

        X = np.zeros( len(RR) )

        posList = zip(X, RR, ZZ)

	result = [np.array(getEField(abs(r), z)) for (x, r, z) in posList]
        vecList = []
        for res in result:
            r, z = res
            norm = np.linalg.norm( res )
            vecList.append( (0., r/norm, z/norm) )

        return np.array(posList), np.array(vecList)

# === Plot Drift Lines ===
def getLines(fn):
    l = []
    with open(fn) as f:
        content = [x.strip() for x in f.readlines()]
        ls = []
        for item in content:
            if item.strip() == '':
                if ls:
                    l.append( ls )
                ls = []
            else:
                ls.append( [float(x) for x in item.split( '\t' )] )

    return l

def drawDrift(fn, threeD=False):
    l = getLines(fn)

    fig = plt.figure()
    if threeD:
        ax = fig.add_subplot(111, projection='3d')
    else:
	ax = fig.add_subplot(111)
	plt.grid()

    for line in l:
        vert = np.array( line )
        x, y, z = vert[:,0], vert[:,1], vert[:,2]

        if threeD:
            ax.plot(x, y, z, color='b')
	else:
	    ax.plot(y, z, color='b')

    fig.show()

    raw_input('--- Press any key ---')

if __name__ == '__main__':
    main()

