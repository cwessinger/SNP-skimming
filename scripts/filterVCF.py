"""
author: Carrie Wessinger

Usage: filterVCF.py infile.vcf outfile.filtered.vcf nsamples minsamples

This script will filter a vcf file based on mapping quality, missing data,
and proportion of reads that match the ref allele (across all individuals)

It also removes sites with multiple ALT alleles.
"""

import sys
import numpy as np

##############################################################
def skipHeader(line):
    """ For a line in VCF file, if line is not a header line, extract columns into list """
    cols = line.replace('\n','').split('\t')
    if len(cols) < 2:
        return 'header'
    elif cols[0] == '#CHROM':
        return 'header'
    else:
        return cols

def extractVCFfields(sampleData):
    """ Extract data from sample-level fields in VCF file """
    if sampleData != './.':
        fields = sampleData.split(':')
        if (len(fields) > 4):
            alleleDepths = fields[1].split(',')
            totDepth = fields[2]
            phreds = fields[4].split(',')
            return [totDepth, alleleDepths, phreds]
        else:
            return 'missing'
    else:
        return 'missing'

###################################################################

vcf = open(sys.argv[1], 'rU')
outvcf = open(sys.argv[2], 'w')

nsamples = int(sys.argv[3])
minsamples = int(sys.argv[4])
minMQ = 20 # minimum mapping quality
min_proportionref = 0.20
max_proportionref = 0.80

indivRange = range(9, 9 + nsamples)
for line in vcf:
    cols = skipHeader(line)
    if cols != 'header':
        INFO = cols[7].split(';')
        listpos = len(INFO) - 1
        MQscore = 0

        while listpos > 0: # search for mapping quality score in INFO field
            if (INFO[listpos].split('='))[0] == 'MQ': # found it!
                MQscore = float(INFO[listpos].split('=')[1])
                listpos = 0
            else:
                listpos -= 1 # keep searching

        # filter based on mapping quality
        if MQscore >= minMQ:
            alt_base = cols[4]

            # skip lines with more than one alt base
            if len(alt_base) > 1:
                pass

            else:
                totREF = 0
                totALT = 0
                calls  = 0

                for j in indivRange:
                    annot = extractVCFfields(cols[j])
                    if annot != 'missing':
                        calls += 1
                        alleleDepths = annot[1]
                        refDepth = int(alleleDepths[0])
                        altDepth = int(alleleDepths[1])
                        totREF += refDepth
                        totALT += altDepth

                proportionref = float(totREF)/float(totREF + totALT)

                if totREF == 0:
                    pass

                elif (max_proportionref > proportionref > min_proportionref) and calls >= minsamples:
                    outvcf.write(line)
outvcf.close()