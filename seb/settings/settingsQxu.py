# Default parameters
#H_BINS = [30, 0, 180**2, 40, -200, 200, 20, -1, 1]
H_BINS = [30, 0, 180**2, 40, -200, 200, 40, -1, 1]
FV = [162, 5, 182]
EN = [700, 3500]

# Paths
preProcessDir = '/home/vault/capm/mppi025h/analysis/preprocess/'
dataList = [preProcessDir + 'data_wThS5_denoised_camp%d.root' % i for i in range(3,6)] + [preProcessDir + 'prepro_ThS5_new_data_trial_camp%d.root' % i for i in range(1,6)]
mcList = [preProcessDir + 'SourceS5_Th228.root']

# SS paths
SSoutDir = '/home/vault/capm/mppi025h/analysis/histo/'
SSoutFile = 'SSqxu_comp.root'

# MS paths
MSoutDir = '/home/vault/capm/mppi025h/analysis/histo/'
MSoutFile = 'MSqxu_comp.root'

# Plot paths
plotDir = '/home/vault/capm/mppi025h/plot_scripts/plots/'
