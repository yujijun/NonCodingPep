import sys, os, traceback
import re

statFile = '/home/binsrv5/encode/stats.xls'
outName = os.path.split(statFile)[1].replace('.xls', '.txt')
outFile = '/u/laumontc/tsaPaper/encode/' + outName

reSRR = re.compile('SRR[0-9]+')
reSpotCount = re.compile('spot_count="[0-9]+"')


with open(statFile, 'r') as sF, open(outFile, 'w') as oF :
    
    oF.write('\t'.join(['SRR', 'spotCount', 'nbReads']) + '\n')
    
    for line in sF:
        line = line.strip('\n')
        srr = re.search(reSRR, line)
        spotCount = re.search(reSpotCount, line)

        # print line
        if srr is not None :
            srr = srr.group()
            spotCount = spotCount.group().split('"')[1]
            nbReads = int(spotCount) * 2
            oF.write('\t'.join([srr, spotCount, str(nbReads)]) + '\n')
            # print srr, spotCount, nbReads
    
