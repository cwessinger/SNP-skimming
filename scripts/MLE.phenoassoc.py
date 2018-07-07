"""
author: Carrie Wessinger

This will take vcf and pheno data and calculate associations

Usage: MLE.phenoassoc.py input.vcf allelefreqs.txt input_allelefreqs.txt input_tks.txt input_phenos.txt nsamples vcfID
"""

import sys
import os
import numpy as np
import pandas as pd
from math import log
from scipy import optimize
from scipy.stats import chi2

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

def scipyLL(params):
    """ function to minimize using optimize.brent """
    global model
    global genoData
    global focalPheno
    global q
    mu = params[0]
    sigma = params[1]
    if model == 'H0':
        negLnL = -H0sumLL(q, mu, sigma, genoData, focalPheno)
    elif model == 'H1':
        a = params[2]
        negLnL = -H1sumLL(q, mu, sigma, a, genoData, focalPheno)
    return negLnL
    
def genoLs(indivData):
    """ assign genotype likelihoods """
    if indivData != 'missing':
        tk = indivData[0]
        phreds = indivData[1]
        if phreds[0] == '0':
            LRR = 1.0
            LRA = (1-tk)/2.0
            LAA = 10 ** (-float(phreds[2])/10.0)
        elif phreds[1] == '0':
            LRR = 10 ** (-float(phreds[0])/10.0)
            LRA = 1.0
            LAA = 10 ** (-float(phreds[2])/10.0)
        elif phreds[2] == '0':
            LRR = 10 ** (-float(phreds[0])/10.0)
            LRA = (1-tk)/2.0
            LAA = 1.0
    else:
        LRR = 1.0
        LRA = 1.0
        LAA = 1.0
    return LRR, LRA, LAA

def H0sumLL(q, mu, sigma, genoData, phenoData):
    """ sum log likelihood under null model that locus does not affect phenotype """
    totLL = 0
    for j in range(len(genoData)):
        LRR, LRA, LAA = genoLs(genoData[j])
        if np.isnan(phenoData[j]) == False:
            z = phenoData[j]
            Ppheno = (1/(np.sqrt(2*sigma*np.pi)))*np.exp((-(z-mu)**2)/(2*sigma))
        else:
            Ppheno = 1.0
        Pdata = Ppheno*LRR*q*q + Ppheno*LRA*2*q*(1-q) + Ppheno*LAA*(1-q)*(1-q)
        if Pdata > 0.0:
            totLL += log(Pdata)
        else:
            return WORST_LL
    return totLL

def H1sumLL(q, mu, sigma, a, genoData, phenoData):
    """ sum log likelihood under alt model that locus affects phenotype """
    muRR = mu
    muRA = mu + a
    muAA = mu + (2*a)
    totLL = 0
    for j in range(len(genoData)):
        LRR, LRA, LAA = genoLs(genoData[j])
        if np.isnan(phenoData[j]) == False:
            z = phenoData[j]
            PdataRR = (1/(np.sqrt(2*sigma*np.pi)))*np.exp((-(z-muRR)**2)/(2*sigma)) * LRR*q*q
            PdataRA = (1/(np.sqrt(2*sigma*np.pi)))*np.exp((-(z-muRA)**2)/(2*sigma)) * LRA*2*q*(1-q)
            PdataAA = (1/(np.sqrt(2*sigma*np.pi)))*np.exp((-(z-muAA)**2)/(2*sigma)) * LAA*(1-q)*(1-q)
            Pdata = PdataRR + PdataRA + PdataAA
        else:
            Pdata = LRR*q*q + LRA*2*q*(1-q) + LAA*(1-q)*(1-q)
        if Pdata > 0.0:
            totLL += np.log(Pdata)
        else:
            return WORST_LL
    return totLL


##########################################################################
vcf = open(sys.argv[1], 'rU')
AFs = open(sys.argv[2], 'rU')
tks = open(sys.argv[3], 'rU')
phenos = pd.read_table(sys.argv[4])
nsamples = int(sys.argv[5])
phenoName = sys.argv[6]
vcfID = sys.argv[7]

indivRange = range(9, 9 + nsamples)
min_kvalue = 1
max_kvalue = 15
sm_val = 0.0001
WORST_LL = float('-inf')

outfile = open('assocs.{0}.{1}.txt'.format(phenoName, vcfID), 'w')
if vcfID == '0':
    outfile.write('contig\tposition\tcount\tq\th0_mu\th0_sigma\th1_muRR\th1_a\th1_sigma\tLRT\tpval\n')

qlist = []
for site in AFs:
    cols = site.replace('\n', '').split('\t')
    if cols[0] != 'contig':
        qlist.append(float(cols[4]))

tklist = [0]
for dat in tks:
    cols = dat.replace('\n', '').split('\t')
    if cols[0] != 'depth':
        tklist.append(float(cols[3]))

focalPheno = phenos[phenoName]

sitecounter = 0
for line in vcf:
    cols = skipHeader(line)
    if cols != 'header':
        contig = cols[0]
        position = cols[1]
        genoData = []
        count = 0
        q = qlist[sitecounter]

        for j in indivRange:

            annot = extractVCFfields(cols[j])
            if annot != 'missing':
                count += 1
                readdepth = int(annot[0])
                phreds = annot[2]
                if readdepth <= max_kvalue:
                    tk = tklist[readdepth]
                else:
                    tk = 1.0
                genoData.append([tk, phreds])
            else:
                genoData.append('missing')

        muMin = focalPheno.min()
        muEst = focalPheno.mean()
        muMax = focalPheno.max()

        sigmaMin = 0.05*focalPheno.var()
        sigmaEst = focalPheno.var()
        sigmaMax = 2*focalPheno.var()

        aMin = -(focalPheno.max() - focalPheno.min())
        aEst = 0.0
        aMax = focalPheno.max() - focalPheno.min()

        # first find likelihood of null model:
        model = 'H0'
        H0InitialValues = np.array([muEst, sigmaEst])
        H0Bounds = [(muMin,muMax), (sigmaMin,sigmaMax)]
        H0MLEs, H0NegLL, d0 = optimize.fmin_l_bfgs_b(scipyLL, x0=H0InitialValues, bounds=H0Bounds, approx_grad=True)

        # next find likelihood of alt model:
        model = 'H1'
        H1InitialValues = np.array([muEst, sigmaEst, aEst])
        H1Bounds = [(muMin,muMax), (sigmaMin,sigmaMax), (aMin,aMax)]
        H1MLEs, H1NegLL, d1 = optimize.fmin_l_bfgs_b(scipyLL, x0=H1InitialValues, bounds=H1Bounds, approx_grad=True)

        H0mu = H0MLEs[0]
        H0sigma = H0MLEs[1]

        H1muRR = H1MLEs[0]
        H1sigma = H1MLEs[1]
        H1a = H1MLEs[2]

        H0LL = -H0NegLL
        H1LL = -H1NegLL
        LRT = 2*(H1LL - H0LL)
        pval = 1.0 - chi2.cdf(LRT, 1)

        sitecounter += 1
        outfile.write(contig+'\t'+position+'\t'+str(count)+'\t'+str(q)+'\t'+str(H0mu)+'\t'+str(H0sigma)+'\t'+str(H1muRR)+'\t'+str(H1a)+'\t'+str(H1sigma)+'\t'+str(LRT)+'\t'+str(pval)+'\n')

outfile.close()
            