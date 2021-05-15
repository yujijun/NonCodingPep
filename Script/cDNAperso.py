import sys, os, traceback

from argparse import ArgumentParser

from pyGeno.Genome import *
from pyGeno.SNPFiltering import SNPFilter, SequenceSNP
Transcript.ensureGlobalIndex('id')



# ################################################################
# Init functions
def outDirCheck(parser, outDir) :
    outDir = os.path.abspath(outDir)
    if not os.path.isdir(outDir) :
        parser.error('The specified output directory %s does not exist!' %(outDir))
    else :
        return outDir

def get_parser():
    parser = ArgumentParser()

    parser.add_argument("-v",
                        dest="REF",
                        help="Version of the reference to be used (Default: GRCm38.87)",
                        default='GRCm38.87',
                        type=str)
    parser.add_argument("-s",
                        dest="SAMPLE",
                        help="Name of the sample (Default: Cancer)",
                        default='Cancer',
                        type=str)
    parser.add_argument("-snp",
                        dest="SNP_SET",
                        help="Name of the snp set of be used (Default: None)",
                        default=None)
    parser.add_argument("-qual",
                        dest="SNP_QUALITY",
                        help="Minimun snp quality (Default: 20)",
                        default=20,
                        type=float)
    parser.add_argument("-o",
                        dest="DIR_OUT",
                        help="Path to the output directory",
                        type=lambda x: outDirCheck(parser, x))

    return parser


# ######################################################################################
# Functions
class Qual_filter(SNPFilter):

    def __init__(self, threshold):
        self.threshold = threshold

    def filter(self, chromosome, **kwargs):
        for snp_set, snp in kwargs.iteritems():
            if float(snp.quality) > self.threshold:
                return SequenceSNP(snp.alt.replace(',', ''))
		# return SequenceSNP(snp.alt)

        return None

def createPersoGenome(ref, snp, qual) :
    try :
        pG = Genome(name = ref, SNPs = snp, SNPFilter = Qual_filter(qual))
        # pG.get(Transcript, id = 'ENSMUST00000070533')[0]
        return pG
    except KeyError:
        return -1
    except ValueError:
        return -2
        
    
def exportcDNA(pG, name, outDir, ref, snp, qual) :
    outFileName = outDir + '/' + '_'.join([name, ref, snp, str(qual)]) + '.txt'
    exclFileName = outFileName.replace('.txt', '_nocDNA.txt')

    try :
        with open(outFileName, 'w') as o, open(exclFileName, 'w') as e :

            header = '\t'.join(['id', 'cDNAseq']) + '\n'
            o.write(header)
            e.write(header)

            # print pG.iterGet(Transcript)
            for t in pG.iterGet(Transcript) :
                # print t
                idT = t.id
                cDNA = t.cDNA

                if len(cDNA) > 0 :
                    line = '\t'.join([idT, cDNA]) + '\n'
                    o.write(line)
                else :
                    l = '\t'.join([idT, 'None']) + '\n'
                    e.write(l)

    except :
        print 'Unexpected error: ', sys.exc_info()[0]
        raise
        


# ###########################################################################
# Main function
def main():
    print("\n---------------------------------------------------------------------")
    print("cDNAperso.py: Create a .txt file containing cDNA personalized versions")
    print("                             Requires everything installed in pyGeno\n")
    print("-----------------------------------------------------------------------")

    global args
    args = get_parser().parse_args()

    # print args.SAMPLE, args.DIR_OUT, args.REF, args.SNP_SET, args.SNP_QUALITY

    print 'Creating personalized genome...'
    persoGenome = createPersoGenome(args.REF, args.SNP_SET, args.SNP_QUALITY)

    if persoGenome == -1 :
        print '    Genome %s is not installed in pyGeno' %(args.REF)
    elif persoGenome == -2 :
        print '    Snp set %s has not been imported in pyGeno' %(args.SNP_SET)
    else :
        print '    Starting to export personalized transcripts...'
        exportcDNA(persoGenome, args.SAMPLE, args.DIR_OUT, args.REF, args.SNP_SET, args.SNP_QUALITY)
        print '    Done!'

main()






