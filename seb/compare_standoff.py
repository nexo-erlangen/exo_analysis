#!/usr/bin/env python
import numpy as np
import plot_support as ps
import plot_functions as pf

# Compare data1 with data2
def main():
    nBins = 120
    sortDataOutDir = '/home/vault/capm/mppi025h/plot_scripts/standoffFiles'

    file1 = 'S5fid'
    file1Dir = '%s/%s/' % (sortDataOutDir, file1)
    file2 = 'ArtDriftOld'
    file2Dir = '%s/%s/' % (sortDataOutDir, file2)
    outFile = '%s_vs_%s' % (file1, file2)

    for art in ['ss', 'ms']:
        # Not working, because MC vs. data is read from same file
        # plotStandoffHistoFancy('standoffRoot/' + fName, art, False, 'standoffPlots/%s/' % sortDataOutFn + fName.split('.')[0] + '%s%d.pdf' % (art.upper(), resBins))
        for region in ['all', 'top', 'bottom']:
            for b in range(3):
                pf.plotStandoffZ(file1Dir + 'mc%s%sZ.root' % (file1, art.upper()), file2Dir + 'mc%s%sZ.root' % (file2, art.upper()), nBins, b, region, False, 'standoffPlots/%s%s%dzStandoff%sBin%d.pdf' % (outFile, art.upper(), nBins, region.title(), b))

if __name__ == '__main__':
    main()

