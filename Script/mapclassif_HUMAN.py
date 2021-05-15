import sys, os, traceback
# from string import maketrans
from argparse import ArgumentParser

sys.path.append('/u/laumontc/tsaPaper/scripts/lib/km/')
from Jellyfish import Jellyfish

sys.path.append('/u/laumontc/tsaPaper/scripts/lib/nektar/nektar_py/')
import utils.common as uc
import utils.kmer as uk

import imp
fm = imp.load_source('fmindex', '/u/laumontc/tsaPaper/scripts/lib/nektar/fmindex/FmIndex.py')

import pandas as pd
import re



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
        parser.error('The specificied output directory %s does not exist!' %(outDir))
    else :
        return outDir

def get_parser():
    parser = ArgumentParser()

    parser.add_argument("-pep",
                        dest="PEP",
                        help="Path to identification file from capaMHC (.csv)",
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-dbFm",
                        dest="DB_FM",
                        help="Path to fmindex of the database used for MS identification (.fm)",
                        default=None,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-dbI",
                        dest="DB_I",
                        help="Path to index of the database used for MS identification (.index)",
                        default=None,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-ct",
                        dest="CONTIGT",
                        help="Sample's assembly.tab from NEKTAR's assembly (Default: None)",
                        default=None,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-pt",
                        dest="PERSO_TRANSCR",
                        help="Personalized transcriptome generated with cDNAperso.py (.txt)",
                        default=None,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-nFM",
                        dest="NORMAL_FM",
                        help="Path to fm index of the normal personalized proteome (.fm)",
                        default=None,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-cFM",
                        dest="CANCER_FM",
                        help="Path to fm index of the cancer personalized proteome (.fm)",
                        default=None,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-nJF",
                        dest="NORMAL_JF",
                        help="Path to kmer database of the normal sample (.jf)",
                        default=None,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-cJF",
                        dest="CANCER_JF",
                        help="Path to kmer database of the cancer sample (.jf)",
                        default=None,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-sgJF",
                        dest="DBSNPG_JF",
                        help="Path to kmer database of the last dbSNP database - genome (.jf)",
                        default=None,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-stJF",
                        dest="DBSNPT_JF",
                        help="Path to kmer database of the last dbSNP database - transcriptome (.jf)",
                        default=None,
                        type=lambda x: is_valid_file(parser, x))
    return parser

# ##############################################################################
# File related functions
def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = ''.join(s.strip() for s in faiter.next())
        yield header, seq


def createPythonDic3(path_file, sep = '\t', index = None) :
    '''
    given a tabulated file, returns a python dictionnary using the user-specified index (0-based) as key
    Therefore, no duplicated values are allowed in the 'index' column.
    '''
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


def load_fmindex(path_to_fm, b_format = 'Nuc4'):
    ''' Takes path to an fnindex file and loads it
        b_format : info on bit storage -- 'Nuc4' for fm done on DNA/RNA | 'AA8' for fm done on amino acids
    '''
    fmindex = fm.FmIndex(b_format)
    fmindex.Load(path_to_fm)
    return fmindex

    
# ##########################################################################
# k-mer related functions

def kmerSet(seq, k) :
    ks = [seq[i:i+k] for i in range(0, len(seq)-k+1)]
    return ks

def kmer_find(fmindex, seq, max_match = 10000000):
    ''' Takes an fmindex and a sequence
        Returns its position in the given fmindex
    '''
    pos_seq = fmindex.SearchSeq(seq, max_match)
    return pos_seq

def kmer_check(fmindex, seq):
    ''' Takes an fmindex and a sequence
        Returns if it was found or not in the given fmindex
    '''
    occ_seq = int(fmindex.CheckSeq(seq))
    return occ_seq

def getCount(kmSet, jf) :
    kc = [jf.query(kS) for kS in kmSet]
    return kc

def cumCount(countList) :
    if sum(True if b == 0 else False for b in countList) > 0 :
        return 0 # only consider a PCR detected when all k-mers in this PCR are detected in the sample
    else :
        return sum(countList)

def countNorm(count, total) :
    res = (float(count) / total) * 10**9
    return res



# ##########################################################################
# ID related functions

def convertIndex(index, indexDic) :
    res = indexDic[index]
    if len(res) == 1 :
        return(res[0])
    else :
        return(res)

def getCodingSeq(peptide, seq, stranded, aa_format) :
    frame_dic = {0 : 'f1', 1 : 'f2', 2 : 'f3', 3 : 'r1', 4 : 'r2', 5 : 'r3' }

    if stranded :
        tr_seq = threeFrameTrans(seq, aa_format)
    else :
        tr_seq = sixFrameTrans(seq, aa_format)
        
    info = [[i, tr_seq[i].find(peptide)] for i in range(0, len(tr_seq)) if tr_seq[i].find(peptide) > -1]
    frame = frame_dic.get(info[0][0], '-1')
    pep_pos = info[0][1] * 3
    
    if frame != -1 :
        seq_in_frame = seqTrans(seq, frame)
        coding_seq = seq_in_frame[pep_pos : pep_pos + (len(peptide) * 3)]
    else :
        coding_seq = None # peptide not found in six-frame translation
        
    return (frame, coding_seq)
    
polym_nt_var = {
    'R' : 'A/G', 'Y' : 'C/T', 'M': 'A/C',
    'K' : 'T/G', 'W' : 'A/T', 'S' : 'C/G',
    'B': 'C/G/T', 'D' : 'A/G/T', 'H' : 'A/C/T',
    'V' : 'A/C/G', 'N': 'A/C/G/T'
}

def return_all_variants(pep) :
    ''' Takes a sequence (nt or aa) with variants and return a list of all possible variants
        ART/KM --> [ARTM, ARKM]
    '''
    # Look for the first variant encoded [A-Z]/[A-Z]/...
    pep = [pep]
    m = re.search('([A-Z| \*]/)+', pep[0]) 

    while m != None :
        # Determine position and identity of the variant in pep
        pattern_pos = (m.start(), m.end() + 1)
        pattern = pep[0][pattern_pos[0] : pattern_pos[1]].split('/')

        # Generate all possible pep as a function of all possible variants
        tmp = []
        for p in pep :
            for pat in pattern :
                tmp.append(p[: pattern_pos[0]] + pat + p[pattern_pos[1] :])
        pep = tmp

        # Update m to determine if there is another variant to consider
        m = re.search('([A-Z| \*]/)+', pep[0])

    # Return list of all possible variants for a given pep
    return pep

def is_valid_seq_nt(seq):
    ''' Takes a nt seq as input and checks if it contains authorized characters only (ATCGN)
    Possible to include one letter nomenclature for polymorphic nts if code the portion for the other functions to handle it
    '''
    s = seq.upper().strip('\n')
    return set(s).issubset({'A', 'T', 'C', 'G', 'U', 'N'})
    # return set(s).issubset({'A', 'T', 'C', 'G', 'U', 'R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N'})


codon_to_aa = {
    'TTT' : ['Phe', 'F'], 'TTC' : ['Phe', 'F'], 'TTA' : ['Leu', 'L'], 'TTG' : ['Leu', 'L'],
    'TCT' : ['Ser', 'S'], 'TCC' : ['Ser', 'S'], 'TCA' : ['Ser', 'S'], 'TCG' : ['Ser', 'S'],
    'TAT' : ['Tyr', 'Y'], 'TAC' : ['Tyr', 'Y'], 'TAA' : ['Stop', '*'], 'TAG' : ['Stop', '*'],
    'TGT' : ['Cys', 'C'], 'TGC' : ['Cys', 'C'], 'TGA' : ['Stop', '*'], 'TGG' : ['Trp', 'W'],
    'CTT' : ['Leu', 'L'], 'CTC' : ['Leu', 'L'], 'CTA' : ['Leu', 'L'], 'CTG' : ['Leu', 'L'],
    'CCT' : ['Pro', 'P'], 'CCC' : ['Pro', 'P'], 'CCA' : ['Pro', 'P'], 'CCG' : ['Pro', 'P'],
    'CAT' : ['His', 'H'], 'CAC' : ['His', 'H'], 'CAA' : ['Gln', 'Q'], 'CAG' : ['Gln', 'Q'],
    'CGT' : ['Arg', 'R'], 'CGC' : ['Arg', 'R'], 'CGA' : ['Arg', 'R'], 'CGG' : ['Arg', 'R'],
    'ATT' : ['Ile', 'I'], 'ATC' : ['Ile', 'I'], 'ATA' : ['Ile', 'I'], 'ATG' : ['Met', 'M'],
    'ACT' : ['Thr', 'T'], 'ACC' : ['Thr', 'T'], 'ACA' : ['Thr', 'T'], 'ACG' : ['Thr', 'T'],
    'AAT' : ['Asn', 'N'], 'AAC' : ['Asn', 'N'], 'AAA' : ['Lys', 'K'], 'AAG' : ['Lys', 'K'],
    'AGT' : ['Ser', 'S'], 'AGC' : ['Ser', 'S'], 'AGA' : ['Arg', 'R'], 'AGG' : ['Arg', 'R'],
    'GTT' : ['Val', 'V'], 'GTC' : ['Val', 'V'], 'GTA' : ['Val', 'V'], 'GTG' : ['Val', 'V'],
    'GCT' : ['Ala', 'A'], 'GCC' : ['Ala', 'A'], 'GCA' : ['Ala', 'A'], 'GCG' : ['Ala', 'A'],
    'GAT' : ['Asp', 'D'], 'GAC' : ['Asp', 'D'], 'GAA' : ['Glu', 'E'], 'GAG' : ['Glu', 'E'],
    'GGT' : ['Gly', 'G'], 'GGC' : ['Gly', 'G'], 'GGA' : ['Gly', 'G'], 'GGG' : ['Gly', 'G'],
}


def seqTrans(seq, frame = 'f1') :
    ''' Takes a nt sequence as input (seq) and format it to be translated in the indicated frame (frame) by the translate_seq
        Returns the formated sequence.
        frame can be: f1, f2, f3 or r1, r2, r3 for the respective forward and reverse translation
    '''

    assert is_valid_seq_nt(seq), 'Invalid sequence'
    seq = seq.upper().replace('U', 'T').strip('\n')
    
    if frame == 'f1' :
        seqTrans = seq
        return seqTrans
    elif frame == 'f2' :
        seqTrans = seq[1:]
        return seqTrans
    elif frame == 'f3' :
        seqTrans = seq[2:]
        return seqTrans
    elif frame == 'r1' :
        seqTrans = get_reverse_complement_seq(seq)
        return seqTrans
    elif frame == 'r2' :
        seqTrans = get_reverse_complement_seq(seq)[1:]
        return seqTrans
    elif frame == 'r3' :
        seqTrans = get_reverse_complement_seq(seq)[2:]
        return seqTrans
    else :
        raise ValueError('Unknown reading frame: %s, should be one of the following: f1, f2, f3, r1, r2, r3' %frame)


def translate_seq(seqTrans, aa_format = 1) :
    ''' Takes a nt sequence and translate it
        aa_format: 0 --> 3 letter code | 1 --> 1 letter code
        In translation output : '*' == stop codon, '-' == unknown Aa because of at least one 'N' in the queried codon
    '''

    assert is_valid_seq_nt(seqTrans), 'Invalid sequence'
    assert aa_format == 0 or aa_format == 1, 'aa_format should be 0 (3 letter code) or 1 (1 letter code)'
    # assert frame == 'f1' or frame == 'f2' or frame == 'f3' or frame == 'r1' or frame == 'r2' or frame == 'r3', 'Unknown reading frame: %s, it should be f1, f2, f3, r1, r2, r3)' %frame

    seqTrans = seqTrans.upper().replace('U', 'T').strip('\n')
    i = 0
    prot = ''
    while len(seqTrans) - (len(seqTrans) - (i + 3)) <= len(seqTrans) :
            # print i
            try :
                ## Avoid throwing unecessarry KeyErrors
                ## Unknown codons are replaced by '-' (can only contain one or more Ns)
                prot += codon_to_aa.get(seqTrans[i:i+3], '-')[aa_format]
                i += 3
            except IndexError :
                i += 3

    return prot


def threeFrameTrans(seq, aa_format = 1) :
    ''' Takes a nt seq and return a tuple containing its 3-frame translation (f1, f2, f3)
        seq : str
    '''
    assert is_valid_seq_nt(seq), 'Invalid sequence'
    assert aa_format == 0 or aa_format == 1, 'aa_format should be 0 (3-letter code) or 1 (1-letter code)'
    
    prot = (translate_seq(seqTrans(seq, 'f1'), aa_format), translate_seq(seqTrans(seq, 'f2'), aa_format),
            translate_seq(seqTrans(seq, 'f3'), aa_format))

    return prot


# cPP, nPP, cPCR, nPCR
ppDic = {'[1, 1, 1, 1]' : 'No',
         '[1, 0, 1, 1]' : 'Yes',
         '[1, 1, 1, 0]' : 'No',
         '[1, 0, 1, 0]' : 'Yes'}
speDic = {'[1, 1, 1, 1]' : 'No',
          '[0, 1, 1, 1]' : 'No',
          '[1, 0, 1, 1]' : 'Yes',
          '[1, 1, 1, 0]' : 'No',
          '[0, 0, 1, 1]' : 'Maybe',
          '[0, 1, 1, 0]' : 'No',
          '[1, 0, 1, 0]' : 'Yes',
          '[0, 0, 1, 0]' : 'Yes'}



def annotation(capaDf, headerFull, cDicSPE, cPersoTranscr, dbFm, dbFmIndex, cPP, nPP, cT, totcT, nT, totnT, dbsnpG, dbsnpT, out) :
    pattern = re.compile('ENST[0-9]+') #ENSMUST for murine samples

    with open(out, 'w') as o :
        o.write('\t'.join(headerFull) + '\n')
        for index, row in capaDf.iterrows() :
            # t1 = time.time()
            pep = row['Peptide sequence']
            # print pep

            # Get peptide count in cancer and normal personalized proteome
            cPPcount = kmer_check(cPP, pep)
            nPPcount = kmer_check(nPP, pep)
            # print cPPcount, nPPcount


            accList = [entry.split(':') for entry in kmer_find(dbFm, pep).split(',')[1:]]
            # print accList

            for acc in accList :
                accIndex, protPos = acc
                accId = convertIndex(accIndex, dbFmIndex)
                cdnaPos = (int(protPos.split('-')[0])-1) * 3

                if accId.find('PP') > -1 :
                    dbOrigin = 'PP'
                    for a in accId.split('|') :
                        pcr = None
                        idT = re.search(pattern, a).group()
                        cDNA = cPersoTranscr[idT][0]

                        pPCR = cDNA[cdnaPos : cdnaPos + (len(pep) * 3)]

                        mutPCR = ''.join([polym_nt_var[c] if c in polym_nt_var else c for c in pPCR])
                        allMutPCR = return_all_variants(mutPCR)
                        pcr = [allM for allM in allMutPCR if pep == translate_seq(allM, aa_format = 1)]
                        # print a, pcr

                        if pcr is not None :
                            for p in pcr :
                                pcrSet = kmerSet(p, cT.k)
                                cCount = cumCount(getCount(pcrSet, cT))
                                cFreq = countNorm(cCount, totcT)
                                nCount = cumCount(getCount(pcrSet, nT))
                                nFreq = countNorm(nCount, totnT)
                                agCode = [cPPcount, nPPcount, cCount, nCount]
                                agBin = [1 if aC > 0 else 0 for aC in agCode]

                                dbsnpGcount = cumCount(getCount(pcrSet, dbsnpG))
                                dbsnpTcount = cumCount(getCount(pcrSet, dbsnpT))

                                # print pep, a, p, agCode, agBin, ppDic.get(str(agBin), 'NA')

                                line = list(row) + [dbOrigin, a, p] + agCode + [cFreq, nFreq] + [str(agBin)] + [ppDic.get(str(agBin), 'NA')] + [dbsnpGcount, dbsnpTcount]
                                strLine = [str(el) for el in line]
                                o.write('\t'.join(strLine) + '\n')
                        else :
                            print 'No pcr associated to %s...' %(pep)
                            # return(-1)


                elif accId.find('SPE') > -1 :
                    dbOrigin = 'SPE'
                    for a in accId.split('|') :
                        pcr = None
                        # seq = cDicSPE[a.split('_')[0]][-1] # 7 for EL4 cells, 6 for CT26 but always last element of the list
                        seq = cDicSPE[a.split('-')[-1].split('_')[0]][-1]
                        pcr = [getCodingSeq(pep, seq, stranded = True, aa_format = 1)[1]] # stranded --> True: 3 frame trans, False: 6 frame trans
                        # print a, pcr

                        if pcr is not None :
                            for p in pcr :
                                pcrSet = kmerSet(p, cT.k)
                                cCount = cumCount(getCount(pcrSet, cT))
                                cFreq = countNorm(cCount, totcT)
                                nCount = cumCount(getCount(pcrSet, nT))
                                nFreq = countNorm(nCount, totnT)
                                agCode = [cPPcount, nPPcount, cCount, nCount]
                                agBin = [1 if aC > 0 else 0 for aC in agCode]
                                
                                dbsnpGcount = cumCount(getCount(pcrSet, dbsnpG))
                                dbsnpTcount = cumCount(getCount(pcrSet, dbsnpT))
                                
                                # print pep, a, p, agCode, agBin, speDic.get(str(agBin), 'NA')

                                line = list(row) + [dbOrigin, a, p] + agCode + [cFreq, nFreq] + [str(agBin)] + [speDic.get(str(agBin), 'NA')] + [dbsnpGcount, dbsnpTcount]
                                strLine = [str(el) for el in line]
                                o.write('\t'.join(strLine) + '\n')
                        else :
                            print 'No pcr associated to %s...' %(pep)
                            # return(-1)

                else :
                    print 'Pb with peptide %s associated to %s id' %(pep, accId)
                    # return(-1)



# ###########################################################################
# Main function
def main():
    print("\n------------------------------------------------------------------")
    print("mapclassif.py: Flag TSA candidates within a list of MAPs (capaMHC)")
    print("------------------------------------------------------------------\n")

    global args
    args = get_parser().parse_args()


    print 'Loading MAPs data...'
    capaDf = pd.read_csv(args.PEP, header = 1, skiprows = 1)
    header = capaDf.columns.values.tolist()
    headerFull = header + ['dbOrigin', 'protAccession', 'PCR', 'cancerProteomeCount', 'normalProteomeCount', 'cancerPCRCount', 'normalPCRcount', 'cancerPCRFreq', 'normalPCRFreq', 'agBinaryCode', 'immStatus', 'dbSNPgenomeCount', 'dbSNPtranscriptomeCount']
    cDicSPE = createPythonDic3(args.CONTIGT, sep = '\t', index = 0)
    cPersoTranscr = createPythonDic3(args.PERSO_TRANSCR, sep = '\t', index = 0)
    print ''

    print 'Loading fm indexes...'
    dbFm = load_fmindex(args.DB_FM, b_format = 'AA8')
    dbFmIndex = createPythonDic3(args.DB_I, sep = '\t', index = 0)
    cPP = load_fmindex(args.CANCER_FM, b_format = 'AA8')
    nPP = load_fmindex(args.NORMAL_FM, b_format = 'AA8')
    print ''

    print 'Loading jellyfish databases...'
    cT = Jellyfish(args.CANCER_JF) # to create Jellyfish class and be able to query in python
    totcT = cT.tot()
    nT = Jellyfish(args.NORMAL_JF)
    totnT = nT.tot()
    dbsnpG = Jellyfish(args.DBSNPG_JF)
    dbsnpT = Jellyfish(args.DBSNPT_JF)
    print ''

    out = args.PEP.replace('.csv', '_annotated.txt')

    print out
    print 'Annotating file...'
    annotation(capaDf, headerFull, cDicSPE, cPersoTranscr, dbFm, dbFmIndex, cPP, nPP, cT, totcT, nT, totnT, dbsnpG, dbsnpT, out)
    print ''

main()
