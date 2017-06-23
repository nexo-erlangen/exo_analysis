from plot_support import *
from plot_functions import *
import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')

def drawHistoCyl(hist):
	setShow(True)
	c = ROOT.TCanvas()
	hist.Project3D('yz').Draw('SURF1 CYL')
	#diffRelHist.Project3D('xz').Draw('SURF1 POL')
	raw_input('end')
	c.Close()
