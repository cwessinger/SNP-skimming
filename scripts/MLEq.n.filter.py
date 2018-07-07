"""
author: Carrie Wessinger

Usage: python MLEq.n.filter.py input.vcf output_prefix nsamples

Finds maximum likelihood estimate for frequency of the ref allele (q).

Generates four output files that begin with output_prefix.

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

def sumLLq(qGuess, qGenoData):
    """ sums log likelihood for a given value of q, given genotype data across all indivs at a site """
    q = qGuess
    totLL = 0.0 # we will increment this
    for k in range(len(qGenoData)):
        if qGenoData[k] != 'missing':
            # "unphred" GATK's genotype likelihoods
            LRR = 10**(-qGenoData[k][0]/10.0)
            LRA = 10**(-qGenoData[k][1]/10.0)
            LAA = 10**(-qGenoData[k][2]/10.0)
        else:
            LRR = 1.0
            LRA = 1.0
            LAA = 1.0
        ldata = q*q*LRR + 2*q*(1-q)*LRA + (1-q)*(1-q)*LAA
        if ldata > 0.0:
            totLL += log(ldata)
        else:
            return WORST_LL
    return totLL

def testBoundaries(data):
    """ check whether highest likelihood is at either of the boundaries of [0, 1] """
    testvals = [0.0, 0.0+sm_val, 1.0-sm_val, 1.0]
    testSumLLs = []
    for v in testvals:
        testSumLLs.append(sumLLq(v, data))            
    if testSumLLs[0] < testSumLLs[1] and testSumLLs[2] > testSumLLs[3]:
        if testSumLLs[1] > testSumLLs[2]:
            testBest = 0.0 + sm_val
        elif testSumLLs[1] < testSumLLs[2]:
            testBest = 1.0 - sm_val        
    elif testSumLLs[0] > testSumLLs[1] and testSumLLs[2] > testSumLLs[3]:
        testBest = 0.0
    elif testSumLLs[0] < testSumLLs[1] and testSumLLs[2] < testSumLLs[3]:
        testBest = 1.0
    else:
        testBest = 'error'
    return testBest

def scipyLLq(qvalue):
    """ function to minimize using optimize.brent """
    global qGenoData
    negLnL = -sumLLq(qvalue, qGenoData)  
    return negLnL

#####################################################

vcf = open(sys.argv[1], 'rU')
prefix = sys.argv[2]
nsamples = int(sys.argv[3])

#output files
polyvcf = open("{0}.poly.vcf".format(prefix), 'w')
polyAFs = open("{0}.poly.Qs.txt".format(prefix), 'w')
fixedvcf = open("{0}.fixed.vcf".format(prefix), 'w')
fixedAFs = open("{0}.fixed.Qs.txt".format(prefix), 'w')

# modify this if needed:
indivRange = range(9, 9 + nsamples) # The VCF columns where the samples appear

sm_val = 0.00001 # small deviation
WORST_LL = float('-inf')

polyAFs.write('contig\tposition\tindivs\tMLE_q\n')
fixedAFs.write('contig\tposition\tindivs\tMLE_q\n')

for line in vcf:
    cols = skipHeader(line)
    if cols != 'header':
        contig = cols[0]
        position = cols[1]
        qGenoData = []
        count = 0
        for j in indivRange:
            annot = extractVCFfields(cols[j])
            if annot != 'missing':
                count += 1
                phreds = annot[2]                
                qGenoData.append([float(phreds[0]), float(phreds[1]), float(phreds[2])])        
        best = testBoundaries(qGenoData)
        if 0<best<1:
            bracket = (0.0, best, 1.0)
            MLE = optimize.brent(scipyLLq, brack = bracket)
            negLL = -sumLLq(MLE, qGenoData)
            polyvcf.write(line)
            polyAFs.write(contig+'\t'+position+'\t'+str(count)+'\t'+str(MLE)+'\n')
        else:
            MLE = best
            negLL = -sumLLq(MLE, qGenoData)
            fixedvcf.write(line)
            fixedAFs.write(contig+'\t'+position+'\t'+str(count)+'\t'+str(MLE)+'\n')
polyAFs.close()
fixedAFs.close()
polyvcf.close()
fixedvcf.close()