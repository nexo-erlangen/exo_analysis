#!/usr/bin/env python
import numpy as np
from StringIO import StringIO
import gzip
import urllib2
import xmltodict
import ast

# Phase I
# =======
# Run range 2464 -> 6370
# Run2a  2464-->3564
# Run2b  3567-->4481
# Run2c  4482-->5599
# Run2d  5600--> 6328
# source: https://confluence.slac.stanford.edu/display/exo/EXOProcessData+--+Phase+I+and+Phase+II+Data+Locations

# runRange = (2464, 6370)
runRange = (6385, 8208)
sourceType = 'Th-228:weak'

def main():
    url = 'http://exo-data.slac.stanford.edu/ExoDatacat/rest/runs'

    request = urllib2.Request( url )
    request.add_header('Accept', 'application/xml')
    request.add_header('Accept-encoding', 'gzip')

    response = urllib2.urlopen( request )

    if response.info().get('Content-Encoding') == 'gzip':
        uzdata = response.read()
        buf = StringIO( uzdata )
        f = gzip.GzipFile(fileobj=buf)
        data = f.read()
    else:
        return

    data = xmltodict.parse( data )
    for run in data['runs']['run']:
        runNum = int(run['run'])
        if run['runType'] != 'Data-Source calibration' or not (runNum >= runRange[0] and runNum <= runRange[1]): # or run['sourceType'] != sourceType:
            continue

        runFilterList = []
        # Sometimes no comments are given
        try:
            if 'S5+1"' in run['comment']:
                print runNum, run['comment']
                runFilterList.append( runNum )
        except:
            continue

    print runFilterList

if __name__ == '__main__':
    main()

