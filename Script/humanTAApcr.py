import sys, os, traceback
from pyGeno.Genome import *


ref = Genome(name = 'GRCh38.88')

inFile = '/u/laumontc/tsaPaper/gtex/natMed/humanTAAsShort.txt'
outFile = os.path.split(inFile)[0] + '/'+ 'taaCandidates_pcrs.txt'

with open(inFile, 'r') as iF, open(outFile, 'w') as oF :
    iF.next() # skip header
    for line in iF :
        sl = line.rstrip().split('\t')

        pep = sl[0]
        ensg = sl[1]

        # print ensg, pep

        g = ref.get(Gene, id = ensg)[0]

        for p in g.get(Protein) :
            pID = p.id
            pSEQ = p.sequence

            pStart = pSEQ.find(pep)
            coding = 'NA'
            if pStart > -1 :
                pEnd = pStart + len(pep)
                # print p.sequence[pStart:pEnd]
                
                tStart = pStart * 3
                tEnd = pEnd * 3

                dnac = p.transcript.cDNA

                coding = dnac[tStart : tEnd]
                
                break
                
        if coding != 'NA' :
            # print pep, coding
            oF.write('\t'.join(['human_TAA', pep + '_1', pep, coding, 'NA', 'NA']) + '\n')
        

