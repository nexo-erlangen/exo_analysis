import os

# Default parameters
#H_BINS = [30, 0, 180**2, 40, -200, 200, 40, -1, 1]
H_BINS_ALL = 60
H_BINS = [H_BINS_ALL, -200, 200, H_BINS_ALL, -200, 200, H_BINS_ALL, -200, 200]
# === RANDOM BINS ===
H_BINS_RAND = [300, -0.2, 0.2, 300, -0.2, 0.2, 300, -0.2, 0.2]

# Set apothem value high so no radial cut is done
FV = [162, 20, 162] 
# FV = [300, 5, 182] 
EN = [980, 9800]

# Paths
preProcessDir = '/home/vault/capm/mppi025h/analysis/preprocess/'
dataList = [preProcessDir + 'data_Run2abc.root']
# dataList = [preProcessDir + 'pre_lbkg_pre.root']
# mcList = [preProcessDir + 'bb2n.root']
mcList = [preProcessDir + 'pre_artDrift_bb2n_pre.root']

# SS paths
SSoutDir = '/home/vault/capm/mppi025h/analysis/histo/'
SSoutFile = 'SS_lbkg.root'

# MS paths
MSoutDir = '/home/vault/capm/mppi025h/analysis/histo/'
MSoutFile = 'MS_lbgk.root'

# Plot paths
plotDir = '/home/vault/capm/mppi025h/plot_scripts/plots/cart_lbgk/'

