#!/usr/bin/env python
import os
import ROOT
ROOT.gSystem.Load('/home/hpc/capm/mppi025h/EXOoffline/EXOFitting/EXO_Fitting/lib/libEXOFitting.so')

preProcessDir = '/home/vault/capm/mppi025h/analysis/preprocess/'

t = ROOT.TChain('dataTree')
# t.Add( preProcessDir + 'prepro_ThS5_new_data_trial_camp3.root' )
t.Add( preProcessDir + 'prepro_ThS5_new_data_trial_camp3.root' )

'''
for fn in os.listdir(preProcessDir + 'S5Th/'):
    runNum = [int(s) for s in fn.split('_') if s.isdigit()][0] 
    if runNum >= 4700 and runNum <= 4800:
        t.Add( preProcessDir + 'S5Th/%s' % fn )
'''

# cut = ROOT.TCut('energy_ms >= 700 && energy_ms < 3000 && multiplicity > 1 && !isMissingPosition && (Min$(abs(cluster_z)) > 20) && Max$(abs(cluster_z)) < 160')
geoCut = ROOT.TCut('abs(cluster_x)<900 && abs(cluster_y)<900 && (Min$(abs(cluster_z)) > 20 && Max$(abs(cluster_z)) < 160 && Max$(abs(cluster_x)*(abs(cluster_x)<900)) < 183 && Max$(abs( 0.5*( ((cluster_z>0) - (cluster_z<=0))*cluster_x*(abs(cluster_x)<900) + sqrt(3.0)*cluster_y*(abs(cluster_y)<900))  )) < 183 && Max$(abs( 0.5*( ((cluster_z<=0) - (cluster_z>0))*cluster_x*(abs(cluster_x)<900) + sqrt(3.0)*cluster_y*(abs(cluster_y)<900))  )) < 183 && Max$(cluster_x*cluster_x*(abs(cluster_x)<900) + cluster_y*cluster_y*(abs(cluster_y)<900)) < 33575.3)')
energyCut = ROOT.TCut('energy_ms >= 1400 && energy_ms <= 1800')
energyCut2 = ROOT.TCut('energy_ms >= 2000 && energy_ms <= 2300')

elseCut = ROOT.TCut('multiplicity > 1 && !isMissingPosition && nsc == 1 && !EventSummary.isDiagonallyCut()')

c = ROOT.TCanvas()
# c.SetLogy()
t.Draw('sqrt(cluster_x*cluster_x + cluster_y*cluster_y)', geoCut + energyCut + elseCut)

c2 = ROOT.TCanvas()
# c2.SetLogy()
t.Draw('sqrt(cluster_x*cluster_x + cluster_y*cluster_y)', geoCut + energyCut2 + elseCut)

raw_input('')

