import sys, os, traceback
from itertools import groupby
from argparse import ArgumentParser
from string import maketrans
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
        parser.error('The specified output directory %s does not exist!' %(outDir))
    else :
        return outDir

def get_parser():
    parser = ArgumentParser()

    parser.add_argument("-n",
                        dest="NAME",
                        help="Sample name (Default: Cancer)",
                        default='Cancer',
                        type=str)
    parser.add_argument("-pp",
                        dest="PP",
                        help="Sample's personalized proteome from capaMHC (Default: None)",
                        default=None,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-cf",
                        dest="CONTIGF",
                        help="Sample's assembly.fasta from NEKTAR's assembly (Default: None)",
                        default=None,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-ct",
                        dest="CONTIGT",
                        help="Sample's assembly.tab from NEKTAR's assembly (Default: None)",
                        default=None,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("--nostrand",
                        dest="STRAND",
                        help="Use when working with unstranded contigs",
                        action='store_false')
    parser.add_argument("-cl",
                        dest="CLEN",
                        help="Minimal contig length for translation (>=, Default: 34)",
                        default=34,
                        type=int)
    parser.add_argument("-pl",
                        dest="PLEN",
                        help="Minimal protein length to be included in MS database (>=, Default: 8)",
                        default=8,
                        type=int)
    parser.add_argument("-o",
                        dest="DIR_OUT",
                        help="Path to the output directory",
                        type=lambda x: outDirCheck(parser, x))
    parser.set_defaults(STRAND=True)

    return parser

# ##############################################################################
# File parser functions
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


# ###############################################################################
# Nucleotide seq functions
def is_valid_seq_nt(seq):
    ''' Takes a nt seq as input and checks if it contains authorized characters only (ATCG)
    Possible to include one letter nomenclature for polymorphic nts if code the portion for the other functions to handle it
    '''
    s = seq.upper().strip('\n')
    return set(s).issubset({'A', 'T', 'C', 'G'})
    

def complement (seq):
    '''
    Complements a DNA sequence, returning the reverse complement.
    '''
    tb = maketrans("ACGTRYMKWSBDHVNacgtrymkwsbdhvn",
                          "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn")        
    return seq[-1::-1].translate (tb)


# ####################################################################################
# Translation functions and variables
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
    seq = seq.upper().strip('\n')
    
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
        seqTrans = complement(seq)
        return seqTrans
    elif frame == 'r2' :
        seqTrans = complement(seq)[1:]
        return seqTrans
    elif frame == 'r3' :
        seqTrans = complement(seq)[2:]
        return seqTrans
    else :
        raise ValueError('Unknown reading frame: %s, should be one of the following: f1, f2, f3, r1, r2, r3' %frame)

        
def translate_seq(seqTrans, aa_format = 1) :
    ''' Takes a nt sequence and translate it
        aa_format: 0 --> 3 letter code | 1 --> 1 letter code
        In translation output : '*' == stop codon, '-' == unknown Aa because of at least one 'N' in the queried codon
    '''

    assert is_valid_seq_nt(seqTrans), 'Invalid sequence'
    assert aa_format == 0 or aa_format == 1, 'aa_format should be 0 (3-letter code) or 1 (1-letter code)'
    # assert frame == 'f1' or frame == 'f2' or frame == 'f3' or frame == 'r1' or frame == 'r2' or frame == 'r3', 'Unknown reading frame: %s, it should be f1, f2, f3, r1, r2, r3)' %frame

    seqTrans = seqTrans.upper().strip('\n')
    i = 0
    prot = ''
    while len(seqTrans) - (len(seqTrans) - (i + 3)) <= len(seqTrans) :
            # print i
            try :
                ## Avoid throwing unecessarry KeyErrors
                ## Unknown codons are replaced by '-'
                prot += codon_to_aa.get(seqTrans[i:i+3], '-')[aa_format]
                i += 3
            except IndexError :
                i += 3

    return prot


def sixFrameTrans(seq, aa_format = 1) :
    ''' Takes a nt seq and return a tuple containing its 6-frame translation (f1, f2, f3, r1, r2, r3)
        seq : str
    '''
    assert is_valid_seq_nt(seq), 'Invalid sequence'
    assert aa_format == 0 or aa_format == 1, 'aa_format should be 0 (3-letter code) or 1 (1-letter code)'
    
    prot = (translate_seq(seqTrans(seq, 'f1'), aa_format), translate_seq(seqTrans(seq, 'f2'), aa_format),
            translate_seq(seqTrans(seq, 'f3'), aa_format), translate_seq(seqTrans(seq, 'r1'), aa_format),
            translate_seq(seqTrans(seq, 'r2'), aa_format), translate_seq(seqTrans(seq, 'r3'), aa_format))

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
    

def filter_translation_output(prot, length_threshold = 8) :
    ''' Takes a tuple of protein sequences as input (prot), split each of them by stop ('*')
        Returns a list of fragments (filtered_prot) having a length >= to length_threshold
    '''
    filtered_prot = []
    l = 0
    
    for p in prot :
        exp_p = p.split('*')
        l += len(exp_p)
        
        for ep in exp_p  :
            if len(ep) >= length_threshold :
                filtered_prot.append(ep)
                
    return [l, filtered_prot]


# ###########################################################################
# Database creation functions

def ppReformat(name, path_to_pp, path_out) :
    ppGen = fasta_iter(path_to_pp)
    
    pattern = re.compile('ENS[A-Z]+[0-9]+')
    protDic = {}
    index = 1


    for entry in ppGen :
        header = entry[0].rstrip()
        protSeq = entry[1].rstrip()

        subHeader = '%s-PP_%s_' %(name, str(index)) + '_'.join(re.findall(pattern, header))
        index += 1

        if protSeq not in protDic :
            protDic[protSeq] = [subHeader]
        else :
            protDic[protSeq].append(subHeader)

    ppName = os.path.basename(path_to_pp).split('.')[0]
    outFile = path_out + '/' + name + '_' + ppName + '_reformat.fasta'

    with open(outFile, 'w') as o :
        for entry in protDic :
            o.write('>' + '|'.join(protDic[entry]) + '\n' + entry + '\n')

    return(outFile)


def contigTrans(name, contigTab, contigFasta, contigLen, strandness, proteinLen, path_out) :
    contigDic = createPythonDic3(contigTab, sep = '\t', index = 0)
    contigSeq = fasta_iter(contigFasta)

    nbF = '3-frame' if strandness else '6-frame'
    protDic = {}
    counter = 0
    outContig = contigTab.replace('.tab','_translatedContigs.tab')

    with open(outContig, 'w') as oC :
        oC.write('ID_SEQ' + '\t' + '\t'.join(contigDic['ID_SEQ']) + '\t' + 'SEQ' + '\n')

        print 'Start contig translation: '
        for entry in contigSeq :
            header = entry[0]
            seq = entry[1]

            counter += 1
            if counter % 50000 == 0 :
                print '    Translated %i contigs...' %(counter)

            if len(seq) >= contigLen :
                if strandness :
                    trans = threeFrameTrans(seq, aa_format = 1)
                else :
                    trans = sixFrameTrans(seq, aa_format = 1)

                transFiltered = filter_translation_output(trans, length_threshold = proteinLen)

                i = 0
                for t in transFiltered[1] :
                    if t not in protDic :
                        protDic[t] = [name + '-SPE_' + 'contig-' + header + '_' + str(i)]
                    else :
                        protDic[t].append(name + '-SPE_' + 'contig-' + header + '_' + str(i))
                    i += 1

                line = header + '\t' + '\t'.join(contigDic[header]) + '\t' + seq + '\n'
                oC.write(line)
        print 'Done with contig translation:'
        print ''

    print 'Writing output file...'
    fastaName = path_out + '/' + name + '_' + '_'.join([str(contigLen), nbF, str(proteinLen)]) + '.fasta'
    with open(fastaName, 'w') as fN :
        for entry in protDic :
            fN.write('>' + '|'.join(protDic[entry]) + '\n' + entry + '\n')
    print 'Done!'
    print ''

    return (outContig, fastaName)


# ###########################################################################
# Main function
def main():
    print("\n----------------------------------------------------------------")
    print("msDbGeneration.py: Generate cancer database to identify TSAs by MS")
    print("------------------------------------------------------------------\n")

    global args
    args = get_parser().parse_args()

    
    if args.PP is not None :
        ppStuff = ppReformat(args.NAME, args.PP, args.DIR_OUT)
        print 'Personalized proteome for sample %s is here: %s' %(args.NAME, ppStuff)
        print ''
        
    elif args.CONTIGF is not None and args.CONTIGT is not None :
        cStuff = contigTrans(args.NAME, args.CONTIGT, args.CONTIGF, args.CLEN, args.STRAND, args.PLEN, args.DIR_OUT)
        print 'File containing translated contigs for sample % s is here: %s' %(args.NAME, cStuff[0])
        print 'Translation of contigs for sample %s is here: %s' %(args.NAME, cStuff[1])
        print '    Parameters used:'
        print '        Minimal contig length: %s' %(str(args.CLEN))
        print '        Contig translation: %s' %('3-frame' if args.STRAND else '6-frame')
        print '        Minimal peptide length: %s' %(str(args.PLEN))
        print ''

    else :
        print 'Specify a file for the -pp or -cf and -ct keys...'
        print ''

main()





        
