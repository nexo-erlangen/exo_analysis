# Default parameters
#H_BINS = [30, 0, 180**2, 40, -200, 200, 40, -1, 1]
H_BINS_ALL = 200
H_BINS = [H_BINS_ALL, -200, 200, H_BINS_ALL, -200, 200, H_BINS_ALL, -200, 200]
# === RANDOM BINS ===
H_BINS_RAND = [200, -0.2, 0.2, 200, -0.2, 0.2, 200, -0.2, 0.2]

# Set apothem value high so no radial cut is done
FV = [300, 0, 300] #[162, 5, 182]
EN = [700, 3500]

# Paths
preProcessDir = '/home/vault/capm/mppi025h/plot_scripts/'
dataList = [preProcessDir + 'standoffDriftData.root'] 
mcList = [preProcessDir + 'standoffDriftData.root']

# SS paths
SSoutDir = '/home/vault/capm/mppi025h/analysis/histo/'
SSoutFile = 'SSstandoffDrift.root'

# MS paths
MSoutDir = '/home/vault/capm/mppi025h/analysis/histo/'
MSoutFile = 'MSstandoffDrift.root'

# Plot paths
plotDir = '/home/vault/capm/mppi025h/plot_scripts/plots/standoffDrift/'
