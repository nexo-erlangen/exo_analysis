import ROOT

preProcessDir = '/home/vault/capm/mppi025h/analysis/preprocess/'
dataList = [preProcessDir + 'data_wThS5_denoised_camp%d.root' % i for i in range(3,6)] + [preProcessDir + 'prepro_ThS5_new_data_trial_camp%d.root' % i for i in range(1,6)]
# mcList = [preProcessDir + 'SourceS5_Th228.root']

def main():
	print 'Getting chain...'
	Tree = getChain('dataTree', dataList)
	f = ROOT.TFile(preProcessDir + 'data_ThS5.root', 'RECREATE')
	Tree.CloneTree(-1, 'fast')

	f.Write()
	f.Close()

	print 'Done'

def getChain(treeName, fileList):
	ch = ROOT.TChain(treeName)
	if(not ch): 
		return False

	for item in fileList:
		ch.Add(item)

	return ch

if __name__ == '__main__':
	main()

