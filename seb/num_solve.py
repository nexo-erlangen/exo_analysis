#!/usr/bin/env python
import numpy as np
import scipy.integrate
import scipy.special

def main():
    # Calculation for an energy of 1MeV
    vD = 170.e3
    X = 380.        # Electric field
    u = vD / X      # Mobility
    u /= 2
    alpha = 1.240   # Recombination coefficient
    b = 0.002       # Initial ionization distribution shape
    Dt = 55.        # Diffusion constant
    Dl = 0.1 * Dt
    d = 0.234       # Length of the column in cm
    N0 = 6.4e4 / d  # Number of initially present ions in a column of 1 cm length

    # Stopping power
    density = 3.057 # g/cm^3
    x, y = np.loadtxt('stopping_power.dat', delimiter='\t', usecols=range(2), unpack=True)
    x = np.array( x ) * 1.e3        # Convert MeV to keV
    y = np.array( y ) * density * 1.e3    # Multiply with density

    # Mean diffusion constant
    meanDiff = meanDiffusion(Dt, Dl)
    print 'Mean diffusion constant:', meanDiff
    
    # Filter by energy
    dataFilt = np.array( [a for a in zip(x, y) if a[0] <= 1000] )
    x, y = dataFilt[:,0], 1./np.array( dataFilt[:,1] )

    d = scipy.integrate.simps(y, x, even='avg')
    print 'Integration of stopping power:', d

    plot(x, [y], axisLabel=[r'Energy $E$ [keV]', r'Stopping power $\frac{1}{S(E)}$ [cm]'])

    # Data
    bRange = np.linspace(0., 0.0003, 1000) 
    # Y_perp_values_low = [Y_perp(alpha, 1000000, Dl, z(b, u, X, Dl)) for b in bRange]
    # Y_perp_values_high = [Y_perp(alpha, 10000000, Dl, z(b, u, X, Dl)) for b in bRange]
    
    Y_perp_values = [Y_perp(alpha, N0, meanDiff, z(b, u, X, meanDiff)) for b in bRange]
    Y_par_values = [Y_par(alpha, b, N0, u, X, Dt, d) for b in bRange]
    X = 567
    vD = 185.e3
    u = vD / X
    Y_perp_p2_values = [Y_perp(alpha, N0, meanDiff, z(b, u, X, meanDiff)) for b in bRange]
    Y_par_p2_values = [Y_par(alpha, b, N0, u, X, Dt, d) for b in bRange]
    # plt.fill_between(bRange, Y_perp_values_low, Y_perp_values_high, alpha=.5)
    # plt.plot(bRange, [li(y2(alpha, b, N0, u, X, D, d)) for b in bRange])

    # plot(bRange, [Y_perp_values, Y_par_values], labelList=['perpendicular', 'parallel'], axisLabel=[r'$b$ [cm]', r'$\mathcal R$'])
    plot(bRange, [Y_perp_values, Y_par_values, Y_perp_p2_values, Y_par_p2_values], labelList=['perpendicular', 'parallel', 'perpendicular (Phase 2)', 'parallel (Phase 2)'], axisLabel=[r'$b$ [cm]', r'$\mathcal R$'])

    raw_input('')

def z(b, u, X, D):
    return b**2 * u**2 * X**2 / (2 * D**2)

def S(z):
    return 1./np.sqrt(np.pi) * np.array( scipy.integrate.quad(lambda s: S_int(s, z), 0, np.inf) )[0]

def S_int(x, z):
    return np.exp(-x) / np.sqrt( x*(1 + x/z) )

def y1(alpha, N0, D):
    return 8*np.pi*D / (alpha * N0)

def y2(alpha, b, N0, u, X, D, d):
    return y1(alpha, N0, D) + np.log((4*D*T(u, X, d) + b**2)/(b**2))

def T(u, X, d):
    return d/(2*u*X)

def li(x):
    try:
        # return scipy.integrate.quad(lambda u: np.exp(u)/u, -np.inf, x)[0]
        return scipy.special.expi(x) 
    except:
        return np.nan

# The field is parallel to the colomuns
def Y_par(alpha, b, N0, u, X, D, d):
    return u/(2*D) * X * b**2/d * y1(alpha, N0, D) * np.exp( -y1(alpha, N0, D) ) * (li(y2(alpha, b, N0, u, X, D, d)) - li(y1(alpha, N0, D)))

# The field is perpendicular to the columns
def Y_perp(alpha, N0, D, z):
    return 1. / (1 + alpha*N0 / (8*np.pi*D) * np.sqrt(np.pi / z) * S(z))

def ellipseRadius(x, a, b):
    return a*b / np.sqrt(a**2*np.cos(x)**2 + b**2*np.sin(x)**2)

def meanDiffusion(Dx, Dy):
    return 1./(2*np.pi) * scipy.integrate.quad(lambda theta: ellipseRadius(theta, Dx, Dy), 0, 2*np.pi)[0]

def plot(x, yList, labelList=None, axisLabel=None):
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick
    from matplotlib import rc
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)

    if not labelList:
        labelList = [None] * len(yList)
    if not axisLabel:
        axisLabel = [None, None]

    f = plt.figure()
    ax = f.add_subplot(111)

    for i, y in enumerate(yList):
        ax.plot(x, y, label=labelList[i])
        # ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
        ax.get_xaxis().get_major_formatter().set_powerlimits((0, 0))
        ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))
    
    ax.set_xlabel(axisLabel[0])
    ax.set_ylabel(axisLabel[1])
    ax.legend(loc='best')
    f.show()

if __name__ == '__main__':
    main()

