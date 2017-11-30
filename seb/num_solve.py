#!/usr/bin/env python
import numpy as np
import scipy.integrate
import scipy.special

def main():
    # Calculation for an energy of 1MeV
    vD = 170.5e3 #  / 2
    X = 380.        # Electric field
    u = vD / X      # Mobility
    # u /= 2
    alpha = 0.001240   # Recombination coefficient
    b = 2.34e-6 # 0.002       # Initial ionization distribution shape in cm^2
    Dt = 55.        # Diffusion constant
    Dl = 0.1 * Dt
    d = 0.05 # 0.233790 # 1.e-2 # .1*0.234       # Length of the column in cm

    # Stopping power
    density = 3.057 # g/cm^3
    x, y = np.loadtxt('stopping_power.dat', delimiter='\t', usecols=range(2), unpack=True)
    x = np.array( x ) * 1.e3        # Convert MeV to keV

    y = np.array( y ) * density * 1.e3    # Multiply with density

    # Mean diffusion constant
    meanDiff = meanDiffusion(Dt, Dl)
    print 'Mean diffusion constant:', meanDiff
    # meanDiff = 40
    
    # Filter by energy
    # dataFilt = np.array( [a for a in zip(x, y) if a[0] <= 1000] )
    x, y = zip(*[item for item in zip(x, y) if item[0] <= 1000])
    y = 1. / np.array( y )

    # Length of track in cm
    dTrack = scipy.integrate.simps(y, x, even='avg')
    print 'Integration of stopping power: %f cm' % dTrack
    # d = dTrack

    # Number of initially present ions in column of 1 cm length
    N0 = 5000 / d # 6.4e4 / d

    # Plot of the stopping power
    plot(x, [y], axisLabel=[r'Energy $E$ [keV]', r'Inverse stopping power $\frac{1}{S(E)}$ [cm]'], despine=True)

    # Data
    bRange = np.linspace(0., 3.e-4, 10000) 
    Y_perp_values = [Y_perp(alpha, N0, meanDiff, z(b, u, X, meanDiff)) for b in bRange]
    Y_par_values = [Y_par(alpha, b, N0, u, X, Dt, d) for b in bRange]
    X = 567
    vD = 185.e3
    u = vD / X
    Y_perp_p2_values = [Y_perp(alpha, N0, Dl, z(b, u, X, Dl)) for b in bRange]
    Y_par_p2_values = [Y_par(alpha, b, N0, u, X, Dt, d) for b in bRange]
    # plt.fill_between(bRange, Y_perp_values_low, Y_perp_values_high, alpha=.5)
    # plt.plot(bRange, [li(y2(alpha, b, N0, u, X, D, d)) for b in bRange])

    # plot(bRange, [Y_perp_values, Y_par_values], labelList=['perpendicular', 'parallel'], axisLabel=[r'$b$ [cm]', r'$\mathcal R$'])

    sigmaRange = 1./np.sqrt(2) * np.array(bRange)
    Ylist = [Y_perp_values, Y_par_values] # Y_perp_p2_values, Y_par_p2_values]
    plot(sigmaRange, Ylist, labelList=['perpendicular', 'parallel', 'perpendicular (Phase 2)', 'parallel (Phase 2)'], axisLabel=[r'$\sigma$ [cm]', r'Diffusion loss'], vline=None, despine=True)

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
    return 1. / (1 + alpha*float(N0) / (8*np.pi*D) * np.sqrt(np.pi / z) * S(z))

def ellipseRadius(x, a, b):
    return a*b / np.sqrt(a**2*np.cos(x)**2 + b**2*np.sin(x)**2)

def meanDiffusion(Dx, Dy):
    return 1./(2*np.pi) * scipy.integrate.quad(lambda theta: ellipseRadius(theta, Dx, Dy), 0, 2*np.pi)[0]

def plot(x, yList, labelList=None, axisLabel=None, vline=None, despine=False):
    from plot_functions import CurvedText
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick
    from matplotlib import rc
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)

    import seaborn as sns
    sns.set_style('whitegrid', {'axes.grid' : False})
    sns.set(style='ticks')

    if not labelList:
        labelList = [None] * len(yList)
    if not axisLabel:
        axisLabel = [None, None]

    f = plt.figure(figsize=(7, 3))
    ax = f.add_subplot(111)

    for i, y in enumerate(yList):
        ax.plot(x, y, label=labelList[i])
        # if labelList[i]:
        #    text = CurvedText(x=x, y=y, text=labelList[i], va='bottom', axes=ax)
        # ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
        ax.get_xaxis().get_major_formatter().set_powerlimits((0, 0))
        ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))

    for i in range(len(yList)/2):
        y1, y2 = np.array(yList[2*i]), np.array(yList[2*i+1])
        idx = np.argwhere(np.diff(np.sign(y1 - y2)) != 0).reshape(-1) + 0
        ax.plot(x[idx], y1[idx], 'ro')
        print (x[idx], y1[idx])

    if vline:
        ax.axvline(x=vline, ls='--')
    
    ax.set_xlabel(axisLabel[0])
    ax.set_ylabel(axisLabel[1])
    ax.legend(loc='best', frameon=False)

    ax.set_xlim(xmin=0)
    ax.set_ylim(ymin=0)

    if despine:
        # sns.despine(fig=f, ax=ax, left=True)
        sns.despine(fig=f, ax=ax)
        # ax.set_yticks([])
    f.show()

if __name__ == '__main__':
    main()

