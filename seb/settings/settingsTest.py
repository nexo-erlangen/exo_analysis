import os

# Default parameters
#H_BINS = [30, 0, 180**2, 40, -200, 200, 40, -1, 1]
H_BINS_ALL = 60
H_BINS = [H_BINS_ALL, -200, 200, H_BINS_ALL, -200, 200, H_BINS_ALL, -200, 200]
# === RANDOM BINS ===
H_BINS_RAND = [300, -0.2, 0.2, 300, -0.2, 0.2, 300, -0.2, 0.2]

# Set apothem value high so no radial cut is done
FV = [162, 5, 182] # [300, 0, 300] 
EN = [700, 3500]

# Paths
preProcessDir = '/home/vault/capm/mppi025h/analysis/preprocess/'
dataList = []
for fn in os.listdir(preProcessDir + 'S5Th/'):
    dataList.append( preProcessDir + 'S5Th/%s' % fn )

# dataList = [preProcessDir + 'prepro_ThS5_new_data_trial_camp1.root']
dataList = [preProcessDir + 'ThS5/Xe134_run_3877_tree.root']

# dataList = [preProcessDir + 'ThS5/' 'ljcalib_2014v1_20130nudenoised_2014v1_run_%s_tree.root' % i for i in [5056, 4449, 4450, 4451, 5057, 5396]]
mcList = [preProcessDir + 'pre_testOld_pre.root']

# SS paths
SSoutDir = '/home/vault/capm/mppi025h/analysis/histo/'
SSoutFile = 'SS_test.root'

# MS paths
MSoutDir = '/home/vault/capm/mppi025h/analysis/histo/'
MSoutFile = 'MS_test.root'

# Plot paths
plotDir = '/home/vault/capm/mppi025h/plot_scripts/plots/test/'
