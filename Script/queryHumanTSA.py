import sys, os, traceback
import pandas as pd
import glob

sys.path.append('/u/laumontc/tsaPaper/scripts/lib/km/')
from Jellyfish import Jellyfish


# Useful functions
def kmerSet(seq, k) :
    ks = [seq[i:i+k] for i in range(0, len(seq)-k+1)]
    return ks

def getCount(kmSet, jf) :
    kc = [jf.query(kS) for kS in kmSet]
    return kc

# def meanCount(countList) :
#     if sum(True if b == 0 else False for b in countList) > 0 :
#         return 0 # only consider a PCR detected when all k-mers in this PCR are detected in the sample
#     else :
#         return int(round((float(sum(countList)) / len(countList)), 0))

def cumCount(countList) :
    if sum(True if b == 0 else False for b in countList) > 0 :
        return 0 # only consider a PCR detected when all k-mers in this PCR are detected in the sample
    else :
        return sum(countList)

def minCount(countList) :
    return min(countList)

def createPythonDic3(path_file, sep = '\t', index = None) :
    dic = {}
    dup = False
    with open(path_file, 'r') as i :
        for lines in i :
            sl = [s.strip('\n') for s in lines.split(sep)]
            if sl[index] not in dic :
                dic[sl[index]] = sl[0:index] + sl[index +1 :]
            else :
                dup = True
                break
        if dup :
            return -1
        else :
            return dic


# Code starts here !!!
agFile = sys.argv[1]
statFile = sys.argv[2]
pathJf = sys.argv[3]
path_out = sys.argv[4]

print >> sys.stdout, pathJf

#Extract tissue and sample id
sampleID = os.path.split(pathJf)[0].split('/')[-1]
tissue = os.path.split(pathJf)[0].split('/')[-2]

# Load jf db
print >> sys.stdout, '%s, %s, %s' %('tissue', 'sampleId', 'totalCount')
jfDb = Jellyfish(pathJf)
dbTot = jfDb.tot()
print >> sys.stdout, '%s, %s, %i' %(tissue, sampleID, dbTot)


# Load agFile
agDf = pd.read_csv(agFile, sep = '\t', header = None)
agDf.columns = ['cancerSampleID', 'mcsID', 'pep', 'mcs', 'cancerPCRCount', 'cancerPCRFreq']

# Create dic from stat file
statDic = createPythonDic3(statFile, sep = '\t', index = 0)

# Queries and compute normalized freq
mCountList = []
cCountList = []
ncCountList = []
mcCountList = []

for index, row in agDf.iterrows() :
    mcs = row['mcs']
    kSet = kmerSet(mcs, jfDb.k)
    
    kCount = getCount(kSet, jfDb)
    mCount = minCount(kCount)
    mCountList.append(mCount)
    mcCount = (float(mCount) * 10**8) / float(statDic[sampleID][1])
    mcCountList.append(mcCount)
    
    cCount = cumCount(kCount)
    cCountList.append(cCount)
    ncCount = (float(cCount) / dbTot)*10**9
    ncCountList.append(ncCount)

    

agDf['agCumCount'] = pd.Series(cCountList)
agDf['agCumNorm'] = pd.Series(ncCountList)
agDf['agMinCount'] = pd.Series(mCountList)
agDf['agMinNorm_100M'] = pd.Series(mcCountList)
agDf['tissue'] = pd.Series([tissue] * len(cCountList))
agDf['SRR'] = pd.Series([sampleID] * len(cCountList))
agDf['nbReads'] = pd.Series([statDic[sampleID][1]] * len(cCountList))

agDf.to_csv('%s/%s_%s.txt' %(path_out, tissue, sampleID), sep = '\t', index = False)
