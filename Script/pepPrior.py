import sys, os, traceback
import pandas as pd
import re
import numpy as np

capa = sys.argv[1]


outNo = capa.replace('.txt', '.excluded')
outYes = capa.replace('.txt', '.yes')
outMaybeYes = capa.replace('.txt', '.maybeYes')
outMaybeNo = capa.replace('.txt', '.maybeNo')
outAmbYes = capa.replace('.txt', '.ambYes')
outAmbNo = capa.replace('.txt', '.ambNo')

df = pd.read_csv(capa, sep = '\t', na_values = 'NA')

pepList = list(set(df['Peptide sequence']))

# print df


ms_aa = {
    'I' : 'I/L', 'L' : 'L/I'
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



with open(outNo, 'a') as oN, open(outYes, 'a') as oY, open(outMaybeYes, 'a') as oMY, open(outMaybeNo, 'a') as oMN, open(outAmbYes, 'a') as oAY, open(outAmbNo, 'a') as oAN :
    header = list(df.columns.values) + ['existingVar']
    oN.write('\t'.join(header) + '\n')
    oY.write('\t'.join(header) + '\n')
    oMY.write('\t'.join(header) + '\n')
    oMN.write('\t'.join(header) + '\n')
    oAY.write('\t'.join(header) + '\n')
    oAN.write('\t'.join(header) + '\n')
    
    for p in pepList :
        sub = df[(df['Peptide sequence'] == p)]
        iS = list(set(sub['immStatus']))

        pVar = ''.join([ms_aa[aa] if aa in ms_aa else aa for aa in p])
        pVarList = [r for r in return_all_variants(pVar) if r != p]

        pVarListFilt = '|'.join([pvl for pvl in pVarList if len(df[(df['Peptide sequence'] == pvl)]) > 0])
        finalpVar = 'NA' if len(pVarListFilt) == 0 else pVarListFilt
        sub['existingVar'] = np.repeat(finalpVar, len(sub))


        if len(iS) == 1 :
            if iS[0] == 'No' :
                sub.to_csv(oN, header = False, sep = '\t', index = False, doublequote = False, mode = 'a')

            elif iS[0] == 'Maybe' :
                ratios = sub['cancerPCRFreq'] / sub['normalPCRFreq']
                if sum([r >= 10.0 for r in ratios]) == len(ratios) :
                    sub.to_csv(oMY, header = False, sep = '\t', index = False, doublequote = False, mode = 'a', na_rep = 'NA')
                    # maybeYes
                else :
                    sub.to_csv(oMN, header = False, sep = '\t', index = False, doublequote = False, mode = 'a', na_rep = 'NA')
                    # maybeNo
            elif iS[0] == 'Yes' :
                sub.to_csv(oY, header = False, sep = '\t', index = False, doublequote = False, mode = 'a', na_rep = 'NA')
                # yes
            else :
                print set(sub['Peptide sequence'])

        else :
            if len([i for i in iS if i == 'NA']) > 0 :
                print set(sub['Peptide sequence'])
            elif len([i for i in iS if i == 'No']) > 0 :
                sub.to_csv(oAN, header = False, sep = '\t', index = False, doublequote = False, mode = 'a', na_rep = 'NA')
                # ambNo
            elif len([i for i in iS if i == 'Maybe']) > 0 :
                subMaybe = sub[(sub['immStatus'] == 'Maybe')]
                ratios2 = subMaybe['cancerPCRFreq'] / subMaybe['normalPCRFreq']
                if sum([r2 >= 10.0 for r2 in ratios2]) == len(ratios2) :
                    sub.to_csv(oAY, header = False, sep = '\t', index = False, doublequote = False, mode = 'a', na_rep = 'NA')
                    # ambYEs
                else :
                    sub.to_csv(oAN, header = False, sep = '\t', index = False, doublequote = False, mode = 'a', na_rep = 'NA')
                    # ambNo
            else :
                print sub
        

    
    


