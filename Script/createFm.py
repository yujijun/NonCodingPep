import sys, os, traceback

import subprocess
from itertools import groupby
from argparse import ArgumentParser


# #########################################################################
# Init functions
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

    parser.add_argument("-f",
                        dest="FASTA",
                        help="Fasta file to use (.fasta)",
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-nbChar",
                        dest="NB_CHAR",
                        help="Number of character per line (Default: 70)",
                        default=70,
                        type=int)
    parser.add_argument("--noprot",
                        dest="PROT",
                        help="Use when working with DNA/RNA",
                        action='store_false')
    parser.add_argument("--nostrand",
                        dest="STRAND",
                        help="Use when working with unstranded data",
                        action='store_false')
    parser.add_argument("-o",
                        dest="DIR_OUT",
                        help="Path to the output directory",
                        type=lambda x: outDirCheck(parser, x))
    parser.set_defaults(PROT=True)
    parser.set_defaults(STRAND=True)

    return parser


# #########################################################################
# File parser
def fasta_iter(fasta_name):
    fh = open(fasta_name)

    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == '>'))

    for header in faiter:
        header = header.next()[1:].strip()
        seq = ''.join(s.strip() for s in faiter.next())
        yield header, seq


# ##########################################################################
# Other functions
def insert_newlines(string = '', every = 0):
    return '\n'.join(string[i : i + every] for i in xrange(0, len(string), every))

def processFasta(fastaGen, char, outFile) :
    try :
        out = outFile.replace('.fasta', '_%schar.fasta' %(str(char)))
        out2 = outFile.replace('.fasta', '_%schar.index' %(str(char)))

        i = 1
        with open(out, 'w') as o, open(out2, 'w') as o2 :
            o2.write('\t'.join(['index', 'id']) + '\n')
            for entry in fastaGen :
                header = entry[0].rstrip()
                protSeq = entry[1].rstrip()

                protSeqN = insert_newlines(protSeq, every = char)

                newLine =  '\n'.join(['>' + str(i), protSeqN]) + '\n'

                o.write(newLine)
                o2.write('\t'.join([str(i), header]) + '\n')

                i += 1

        return out

    except :
        print 'Unexpected error: ', sys.exc_info()[0]
        raise

def fastaToFm(fastaProc, prot, strand) :
    outDir = fastaProc.strip('.fasta')

    if prot :
        cmd = 'python /u/laumontc/tsaPaper/scripts/lib/nektar/nektar_py/nektar.py create_fmindex -f %s -o %s -prot -rev' %(fastaProc, outDir)
        print cmd
        process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        print output

    else :
        if strand :
            cmd = 'python /u/laumontc/tsaPaper/scripts/lib/nektar/nektar_py/nektar.py create_fmindex -f %s -o %s -rev' %(fastaProc, outDir)
            process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
            print output
        else :
            cmd = 'python /u/laumontc/tsaPaper/scripts/lib/nektar/nektar_py/nektar.py create_fmindex -f %s -o %s ' %(fastaProc, outDir)
            process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
            print output


# #########################################################################
# Main function

def main():
    print("\n------------------------------------------------------------------------")
    print("createFm.py: Pre-process and create an fm index for nuc/protein fasta file")
    print("------------------------------------------------------------------------\n")

    global args
    args = get_parser().parse_args()

    outName = args.DIR_OUT + '/' + os.path.basename(args.FASTA)
    fG = fasta_iter(args.FASTA)

    print 'Pre-processing fasta file...'
    processedFasta = processFasta(fG, args.NB_CHAR, outName)
    # print processedFasta
    print 'Done!'
    print ''

    print 'Generating fm-index...'
    fastaToFm(processedFasta, args.PROT, args.STRAND)
    print 'Done!'
    print ''

main()
