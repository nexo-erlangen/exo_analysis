import os

# Default parameters
#H_BINS = [30, 0, 180**2, 40, -200, 200, 40, -1, 1]
H_BINS_ALL = 60
H_BINS = [H_BINS_ALL, -200, 200, H_BINS_ALL, -200, 200, H_BINS_ALL, -200, 200]
# === RANDOM BINS ===
H_BINS_RAND = [300, -0.2, 0.2, 300, -0.2, 0.2, 300, -0.2, 0.2]

# Set apothem value high so no radial cut is done
FV = [162, 5, 182]
EN = [700, 3500]

# Paths
preProcessDir = '/home/vault/capm/mppi025h/analysis/preprocess/'
dataList = []
for fn in os.listdir(preProcessDir + 'S5Th/'):
    dataList.append( preProcessDir + 'S5Th/%s' % fn )

mcList = [preProcessDir + 'pre_drift0_172_pre.root']

# SS paths
SSoutDir = '/home/vault/capm/mppi025h/analysis/histo/'
SSoutFile = 'SS_S5Th.root'

# MS paths
MSoutDir = '/home/vault/capm/mppi025h/analysis/histo/'
MSoutFile = 'MS_S5Th.root'

# Plot paths
plotDir = '/home/vault/capm/mppi025h/plot_scripts/plots/cart_S5Th/'
