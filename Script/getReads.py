import sys, os, traceback
from argparse import ArgumentParser
import pandas as pd
import gzip

from itertools import islice
from string import maketrans



# ###########################################################################
# Init function
def is_valid_file(parser, p_file):
    p_file = os.path.abspath(p_file)
    if not os.path.exists(p_file):
        parser.error("The file %s does not exist!" % p_file)
    else:
        return p_file

def outDirCheck(parser, outDir) :
    outDir = os.path.abspath(outDir)
    if not os.path.isdir(outDir) :
        parser.error('The specified output directory %s does not exist!' %(outDir))
    else :
        return outDir

def get_parser():
    parser = ArgumentParser()

    parser.add_argument("-tsa",
                        dest="TSA",
                        help="List of TSA candidates (one per line, no header, .txt extension required)",
                        default=None,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-maps",
                        dest="MAPS",
                        help="Annotated file of all MAPs (capaMHC)",
                        default=None,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("--notpaired",
                        dest="PAIRED",
                        help="Use when working with unpaired RNA-seq data",
                        action='store_false')
    parser.add_argument("-fastq1",
                        dest="FASTQ1",
                        help="Path to R1 fastq file (fastq.gz)",
                        default=None,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-fastq2",
                        dest="FASTQ2",
                        help="Path to R2 fastq file (fastq.gz)",
                        default=None,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-o",
                        dest="DIR_OUT",
                        help="Path to the output directory",
                        type=lambda x: outDirCheck(parser, x))
    parser.set_defaults(PAIRED=True)

    return parser

# ################################################################################################
# Useful functions
def complement (seq) :
    '''
    Complements a DNA sequence, returning the reverse complement.
    '''
    tb = maketrans("ACGTRYMKWSBDHVNacgtrymkwsbdhvn",
                   "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn")        
    return seq[-1::-1].translate (tb)


def codingSeq(tsaFile, mapsAnnotated, outPath) :
    outFile = outPath + '/' + os.path.split(tsaFile)[1].replace('.txt', '_pcrs.txt')
    df = pd.read_csv(mapsAnnotated, header = 0, sep = '\t')

    with open(tsaFile, 'r') as tF, open(outFile, 'w') as oF :
        for line in tF :
            pep = line.rstrip()

            pepID = pep

            fdf = df[df['Peptide sequence'] == pep]
            pcrs = set(fdf['PCR'].values)

            i = 1
            for p in pcrs :
                pCount = pd.unique(df[df['PCR'] == p]['cancerPCRCount'])
                pFreq = pd.unique(df[df['PCR'] == p]['cancerPCRFreq'])
                nLine = '\t'.join([pepID + '_' + str(i), pep, p, str(pCount[0]), str(pFreq[0])]) + '\n'
                oF.write(nLine)
                i+=1

        return(outFile)


def tsaReads(tsaPcrs, paired, fastq1, fastq2, outPath) :
    # Generate dictionnaries with open file for both fwd (R2) and rc (R1)
    kDic = {}
    krevDic = {}
    with open(tsaPcrs, 'r') as tP :
        for line in tP :
            kID, pep, PCR = line.rstrip().split('\t')[0:3]
            rcPCR = complement(PCR)

            if paired :
                f1 = outPath + '/' + kID + '_R1_paired.fasta'
                f2 = outPath + '/' + kID + '_R2_paired.fasta'

                if PCR not in kDic :
                    kDic[PCR] = [open(f2, 'w')]
                    krevDic[rcPCR] = [open(f1, 'w')]
                else :
                    print '%s already exists in dic' %(PCR)

            else :
                f1 = outPath + '/' + kID + '_R1_unpaired.fasta'
                f2 = outPath + '/' + kID + '_R2_unpaired.fasta'

                if PCR not in kDic :
                    kDic[PCR] = [open(f2, 'w')]
                    krevDic[rcPCR] = [open(f1, 'w')]
                else :
                    print '%s already exists in dic' %(PCR)

    # Strating looking for TSA-coding reads
    pcrList = kDic.keys()
    revpcrList = krevDic.keys()
    i = 0

    try :
        if paired :
            with gzip.open(fastq1, 'r') as r1, gzip.open(fastq2, 'r') as r2 : # r1

                while True :
                    i += 1
                    if i % 1000 == 0 :
                        print 'Analyzed %s lines...' %(str(i*2048))
                    next_n_lines1 = list(islice(r1, 2048))
                    next_n_lines2 = list(islice(r2, 2048))

                    if not next_n_lines1:
                        break

                    cptL = 1
                    while cptL < len(next_n_lines1):
                        read1 = next_n_lines1[cptL]
                        qual1 = '>' + next_n_lines1[cptL-1]
                        read2 = next_n_lines2[cptL]
                        qual2 = '>' + next_n_lines2[cptL-1]

                        posPcr = [kDic[pcr][0].write(qual2 + read2) for pcr in pcrList if read2.find(pcr) > -1] # r2
                        posrevPCR = [krevDic[revpcr][0].write(qual1 + read1) for revpcr in revpcrList if read1.find(revpcr) > -1] #r1

                        cptL += 4

        else :
            fileList = [fastq1, fastq2]

            for fl in fileList :
                i = 0
                print fl

                readType = fl.find('R1')

                with gzip.open(fl, 'r') as r :

                    while True :
                        i += 1
                        if i % 1000 == 0 :
                            print 'Analyzed %s lines...' %(str(i*2048))

                        next_n_lines = list(islice(r, 2048))

                        if not next_n_lines :
                            break

                        cptL = 1
                        while cptL < len(next_n_lines):
                            read = next_n_lines[cptL]
                            qual =  '>' + next_n_lines[cptL-1]

                            if readType > -1 :
                                posrevPCR = [krevDic[revpcr][0].write(qual + read) for revpcr in revpcrList if read.find(revpcr) > -1]
                            else :
                                posPcr = [kDic[pcr][0].write(qual + read) for pcr in pcrList if read.find(pcr) > -1]

                            cptL += 4
        return(0)

    except:
        print("Error in openning file")
        traceback.print_exc()
        return(-1)


    close_r1 = [kDic[pcr][0].close() for pcr in pcrList]
    close_r2 = [krevDic[revpcr][0].close() for revpcr in revpcrList]


# ###########################################################################
# Main function
def main():
    print("\n---------------------------------------------------------------------")
    print("getReads.py: Return fastas containing RNA-seq reads")
    print("                        encoding TSA candidates (stranded data only)")
    print("---------------------------------------------------------------------\n")

    global args
    args = get_parser().parse_args()

    # Retrieve coding sequences
    print 'Retrieving all PCRs for each TSA candidate...'
    tsaPcr = codingSeq(args.TSA, args.MAPS, args.DIR_OUT)
    print 'Done!'
    print ''

    # Retrieve reads for each pcrs
    print 'Retrieving reads associated to each TSA/pcr...'
    tsaReads(tsaPcr, args.PAIRED, args.FASTQ1, args.FASTQ2, args.DIR_OUT)
    print 'Done!'

main()
