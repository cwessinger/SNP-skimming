"""
author: Carrie Wessinger

Usage: python MLEq.reest.py input.vcf input.qs.txt input.tks.txt output.reestimated.qs.txt nsamples

Generates an output file that contains Tk values for each individual for k = [1-15]
"""

import sys
import os
from math import log
from scipy import optimize

#############################################################

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

def sumLLq(qGuess, qGenoData): # sums LL for q, given genotype data across all indivs at a site
    """ sums log likelihood for a given value of q """
    q = qGuess
    totLL = 0.0 # we will increment this

    for k in range(len(qGenoData)):
        if qGenoData[k] != 'missing':
            tk = qGenoData[k][0]
            phreds = qGenoData[k][1]
            if phreds[0] == 0:
                ldata = q*q + q*(1-q)*(1-tk)
            elif phreds[1] == 0:
                ldata = 2*q*(1-q)*tk
            elif phreds[2] == 0:
                ldata = (1-q)*(1-q) + q*(1-q)*(1-tk)

            if ldata > 0.0:
                totLL += log(ldata)
            else:
                return WORST_LL
    return totLL

def scipyLLq(qvalue):
    """ function to minimize using optimize.brent """
    global qGenoData
    negLnL = -sumLLq(qvalue, qGenoData)
    return negLnL

#####################################################

vcf = open(sys.argv[1], 'rU')
AFs = open(sys.argv[2], 'rU')
tks = open(sys.argv[3], 'rU')
outfile = open(sys.argv[4], 'w')
nsamples = int(sys.argv[5])

indivRange = range(9, 9 + nsamples)
min_kvalue = 1
max_kvalue = 15
sm_val = 0.00001 # small deviation
WORST_LL = float('-inf')

outfile.write('contig\tposition\tindivs\tinitial_q\treestimate_q\n')

qlist = []
for site in AFs:
    cols = site.replace('\n', '').split('\t')
    if cols[0] != 'contig':
        qlist.append(float(cols[3]))

tklist = [0]
for dat in tks:
    cols = dat.replace('\n', '').split('\t')
    if cols[0] != 'depth':
        tklist.append(float(cols[3]))

sitecounter = 0
newlyfixed = 0
for line in vcf:
    cols = skipHeader(line)
    if cols != 'header':
        qEst = qlist[sitecounter]
        contig = cols[0]
        position = cols[1]
        qGenoData = []
        count = 0
        for j in indivRange:
            annot = extractVCFfields(cols[j])
            if annot != 'missing':
                count += 1
                sitedepth = int(annot[0])
                if sitedepth <= max_kvalue:
                    tk = tklist[sitedepth]
                else:
                    tk = 1.0
                phreds = annot[2]
                qGenoData.append([tk, [float(phreds[0]), float(phreds[1]), float(phreds[2])]])
            else:
                qGenoData.append('missing')
        q0sumLL = sumLLq(0.0 + sm_val, qGenoData)
        q1sumLL = sumLLq(1.0 - sm_val, qGenoData)
        qMsumLL = sumLLq(qEst, qGenoData)

        if (qMsumLL > q0sumLL) and (qMsumLL > q1sumLL):
            bracket = (0.0, qEst, 1.0)
            MLE = optimize.brent(scipyLLq, brack=bracket)

        elif q0sumLL >= qMsumLL:
            MLE = 0.0
            newlyfixed += 1

        elif q1sumLL >= qMsumLL:
            MLE = 1.0
            newlyfixed += 1
        negLL = -sumLLq(MLE, qGenoData)

        sitecounter += 1
        outfile.write(contig+'\t'+position+'\t'+str(count)+'\t'+str(qEst)+'\t'+str(MLE)+'\n')
print('sites that now appear fixed: ', str(newlyfixed))
outfile.close()