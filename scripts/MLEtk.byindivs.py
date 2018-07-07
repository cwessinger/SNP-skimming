"""
author: Carrie Wessinger

Usage: python MLEtk.allindivs.py input.vcf input.qs.txt output.tks.txt nsamples

Generates an output file that contains Tk values for each individual for k = [1-15]
"""

# This is a script to estimate T_k based, using previously estimated values of q_i
import sys
import os
from math import log
from scipy import optimize
from scipy.stats import chi2
import numpy as np

##########################################################################

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

def genoCounts(data):
    """ genotype counts from compiled data for a given read depth """
    RR = 0
    RA = 0
    AA = 0
    for d in range(len(data)):
        genos = data[d][0]
        if genos[0] == '0':
            RR += 1
        elif genos[1] == '0':
            RA += 1
        elif genos[2] == '0':
            AA += 1
    tot = float(RR + RA + AA)
    if tot == 0:
        pRR = 0
        pRA = 0
        pAA = 0
    else:
        pRR = RR/tot
        pRA = RA/tot
        pAA = AA/tot
    return pRR, pRA, pAA, tot


def sumLLtk(tkGuess, data):
    """ sums log likelihood for a given value of tk, given genotype data across all indivs and all sites with read depth k"""
    tk = tkGuess
    totLL = 0.0  # we will increment this
    for d in range(len(data)):
        q = data[d][1]
        phreds = data[d][0]
        if phreds[0] == '0':
            ldata = q * q + q * (1 - q) * (1 - tk)
        elif phreds[1] == '0':
            ldata = 2 * q * (1 - q) * tk
        elif phreds[2] == '0':
            ldata = (1 - q) ** 2 + q * (1 - q) * (1 - tk)

        if ldata > 0.0:
            totLL += log(ldata)
        else:
            return WORST_LL
    return totLL


def testBoundaries(data):
    """ check whether highest likelihood is at either of the boundaries of [0, 1] """
    testvals = [0.0+sm_val, 0.0+2*sm_val, 1.0-2*sm_val, 1.0-sm_val]
    testSumLLs = []
    for v in testvals:
        testSumLLs.append(sumLLtk(v, data))
    if testSumLLs[0] < testSumLLs[1] and testSumLLs[2] > testSumLLs[3]:
        if testSumLLs[1] > testSumLLs[2]:
            testBest = 0.0 + 2*sm_val
        elif testSumLLs[1] < testSumLLs[2]:
            testBest = 1.0 - 2*sm_val
    elif testSumLLs[0] > testSumLLs[1] and testSumLLs[2] > testSumLLs[3]:
        testBest = 0.0 + sm_val
    elif testSumLLs[0] < testSumLLs[1] and testSumLLs[2] < testSumLLs[3]:
        testBest = 1.0 - sm_val
    else:
        testBest = 'error'
    return testBest
        
def scipyLLtk(tkvalue):
    """ function to minimize using optimize.brent """
    global currData
    negLL = -sumLLtk(tkvalue, currData)        
    return negLL       

def likeRatioTest(data, sumLLbest):
    """ is the optimized tk value a signif. better fit than tk = 1? """
    LLA = sumLLbest
    LL0 = sumLLtk(1.0, data)
    LRT = 2*(LLA-LL0)
    pval = 1.0 - chi2.cdf(LRT, 1)
    return(LRT, pval)

def interpretResults(best, currData):
    if best == 'error':
        MLE = 'NA'
        negLL = 'NA'
        LRT = 'NA'
        pval = 'NA'
    elif best == 1 - sm_val:
        MLE = 1.0
        negLL = -sumLLtk(1.0-sm_val, currData)
        LRT = 0.0
        pval = 1.0
    elif best == 0+sm_val:
        MLE = 0.0
        negLL = -sumLLtk(0.0+sm_val, currData)
        LRT, pval = likeRatioTest(currData, -negLL)
    else:
        MLE = optimize.brent(scipyLLtk, brack = (0.0+sm_val, best, 1.0-sm_val))
        negLL = -sumLLtk(MLE, currData)
        LRT, pval = likeRatioTest(currData, -negLL)
    return MLE, negLL, LRT, pval

#######################################################################

vcf = open(sys.argv[1], 'rU')
AFs = open(sys.argv[2], 'rU')
outfile = open(sys.argv[3], 'w')
nsamples = int(sys.argv[4])

indivRange = range(9, 9 + nsamples)
min_kvalue = 1
max_kvalue = 15
sm_val = 0.0001 # small deviation
WORST_LL = float('-inf')

outfile.write('indiv\tdepth\ttotsites\tpropRA\tMLE_tk\n')

qlist = []
for site in AFs:
    cols = site.replace('\n', '').split('\t')
    if cols[0] != 'contig':
        qlist.append(float(cols[3]))

gatkData = [[[] for k in range(max_kvalue)] for j in indivRange]
sitecounter = 0
for line in vcf:
    cols = skipHeader(line)
    if cols != 'header':
        for j in indivRange:
            annot = extractVCFfields(cols[j])
            q = qlist[sitecounter]
            if annot != 'missing':
                totDepth = int(annot[0])
                phreds = annot[2]
                if min_kvalue <= totDepth <= max_kvalue:
                    gatkData[j-9][totDepth - min_kvalue].append([phreds, q])
                else:
                    pass
            else:
                pass # since lists are based on read depth, can't do anything with missing data
        sitecounter += 1

for j in indivRange:
    for k in range(min_kvalue, max_kvalue+1):
        currData = gatkData[j-9][k-min_kvalue]
        testBest = testBoundaries(currData)
        MLE, NegLL, LRT, Pval = interpretResults(testBest, currData)
        RR, RA, AA, tot = genoCounts(currData)

        outfile.write(str(j-9)+'\t'+str(k)+'\t'+str(tot)+'\t'+str(RA)+'\t'+str(MLE)+'\n')

outfile.close()
        
