#!/usr/bin/env python
import numpy as np
import generate_random as gr
from scipy.optimize import fsolve

# Fit parameter
b = 0.00261
z_0 = np.pi # 1. - b/2.

# m(r)
c_m = 1.68588137007714
d_m = -0.2924038847091
c_m2 = 0.0408885891201855 # -0.0177060260082422
d_m2 = 0. # 0.0099479324509953

# t(r)
m_t = 1.00109952487923

# A(r)
c_A = 6.01494860061481
d_A = 0.173172199655505
# c_A2 = 0.000143566044902527 # -0.00292957115068489
# d_A2 = 0. # -3.19103731564881e-05

c_A2 = -0.00292957115068489
d_A2 = -3.19103731564881e-05

REFLECTORINNERRAD = 183.2356e-3

def main():
    zStart = -0.01
    startPos = 0.16 # REFLECTORINNERRAD - 0.1
    endPos = REFLECTORINNERRAD - 0.0005

    fN = open('gp_data/artTest.dat', 'w')
    for r in np.linspace(startPos, endPos, 20):
        for zStart in np.linspace(-0.003, -0.17, 5):
            rNew = getR(r, zStart)
            # print r, rNew, f(rNew, zStart)

            '''
            isInt = isIntersected(rNew)
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

def artificialDrift(startPos, eventNum=0, out_q=None, diffuse=False, saveHistory=False):
    # Important constants
    x_CAD = 0.19223657
    r_RID = 0.1832356
    limit = 0.190

    # Cartesian and polar vector
    xVec = list( startPos )
    posHistory = [ xVec ]
    x, y, z = xVec
    theta = np.arctan2(y, x)
    r = np.sqrt( x**2 + y**2 )

    # Flip if in wrong section of volume
    upFlag = False
    leftFlag = False
    if xVec[2] > 0:
        upFlage = True
        xVec[2] *= -1

    # Get correct radius
    r_ = getR(r, z)
    xVecPol = np.array( [r, z] )

    for z in np.linspace(z, limit, 1000):
        res = f(r_, z)

        if saveHistory:
            xVec = list( gr.polToCart( res, z, 0. ) )
            posHistory.append( xVec )

        if res > r_RID:
            posHistory = False
            break

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
    
def f(r, z): 
    # return m(r) * z + A(r) * S(z) + t(r)
    
    return m(r) * z + A(r) * S(z) + t(r)

def m(r):
    # return np.exp( c_m * (r - d_m) )
    # return 1 + c_m * (r - d_m) + (c_m * (r - d_m))**2

    limit = (d_m2 - d_m) / (c_m - c_m2)
    if r < limit:
        return c_m2 * r + d_m2
    else:
        return c_m * r + d_m

def t(r):
    # return m_t * r + t_t
    return m_t * r

def A(r):
    # return -np.exp( c_A * (r - d_A) )
    # return -( 1 + c_A * (r - d_A) + (c_A * (r - d_A))**2 )

    if r < d_A:
        # return c_A2 * (r - d_A) + d_A2
        return d_A2 / d_A * r
    else:
        return -c_A * (r - d_A)**2 + d_A2

def S(z):
    return np.sin(-(z - z_0)/b)
    # return 0

def getR(r_x, z_x):
    # (r + z*(c_m*d_m - 1) + S(z) * (1 - c_A*d_A) - t_t) / (c_m * z - c_A * S(z) + m_t)
    # return r / m_t

    '''
    func = lambda r : -np.exp(c_A * (r - d_A)) * S(z_x) + np.exp(c_m * (r - d_A))*z_x + m_t*r - r_x
    r_init = 0.16
    r_sol = fsolve(func, r_init)

    return r_sol[0]
    '''

    if r_x > d_A:
        b = -( (c_m * z_x + m_t) / (c_A * S(z_x)) + 2*d_A )
        c = (r_x - d_m * z_x) / (c_A * S(z_x)) - d_A2 / c_A + d_A**2

        r1 = ( -b + np.sqrt( b**2 - 4*c ) ) / 2
        r2 = ( -b - np.sqrt( b**2 - 4*c ) ) / 2

        if r1 < 0. or r1 > 0.2:
                return r2
        else:
                return r1
    else:
        return (r_x + c_A2 * d_A * S(z_x) - d_m2 * z_x - d_A2 * S(z_x)) / (c_m2 * z_x + c_A2 * S(z_x) + m_t)

def isIntersected(r):
    if r < d_A:
        func = lambda z : (c_m2 * r + d_m2) * z + (c_A2 * (r - d_A) + d_A2) * S(z) + m_t*r - REFLECTORINNERRAD
    else:
        func = lambda z : (c_m * r + d_m) * z + (-c_A * (r - d_A)**2 + d_A2) * S(z) + m_t*r - REFLECTORINNERRAD
    
    #func = lambda z : -np.exp(c_A * (r - d_A)) * S(z) + np.exp(c_m * (r - d_A))*z + m_t*r - REFLECTORINNERRAD
    z_init = -0.
    z_sol = fsolve(func, z_init)

    if z_sol > -0.18 and z_sol < 0.:
        return z_sol
    else:
        return False

if __name__ == '__main__':
    main()

