'''
Simulates a genome of 10010 loci, 10 of which have sparse effects on phenotype,
runs SNP-skimming likelihood method to find associations.
Also outputs a bimbam file for use in gemma and outputs the simulated phenotypes.
'''

import sys
import os
import random
import numpy as np
from math import log
from scipy import optimize
from scipy.stats import chi2
import pandas as pd


#############################################################

def scipyLL(params):
    global model
    global genoData
    global phenos
    q = params[0]
    mu = params[1]
    sigma = params[2]
    if model == 'H0':
        negLnL = -H0sumLL(q, mu, sigma, genoData, phenos)
    elif model == 'H1':
        a = params[3]
        negLnL = -H1sumLL(q, mu, sigma, a, genoData, phenos)
    return negLnL


def genoLs(indivData):
    if indivData[1] != 'missing':
        tk = indivData[0]
        obs = indivData[1]
        if obs == 'R':
            LRR = 1.0
            LRA = (1 - tk) / 2.0
            LAA = 0.0
        elif obs == 'H':
            LRR = 0.0
            LRA = 1.0
            LAA = 0.0
        elif obs == 'A':
            LRR = 0.0
            LRA = (1 - tk) / 2.0
            LAA = 1.0
    else:
        LRR = 1.0
        LRA = 1.0
        LAA = 1.0
    return LRR, LRA, LAA


def H0sumLL(q, mu, sigma, genoData, phenos):
    totLL = 0
    for j in range(len(genoData)):
        LRR, LRA, LAA = genoLs(genoData[j])
        if np.isnan(phenos[j]) == False:
            z = phenos[j]
            Ppheno = (1 / (np.sqrt(2 * sigma * np.pi))) * np.exp((-(z - mu) ** 2) / (2 * sigma))
        else:
            Ppheno = 1.0
        Pdata = Ppheno * LRR * q * q + Ppheno * LRA * 2 * q * (1 - q) + Ppheno * LAA * (1 - q) * (1 - q)
        if Pdata > 0.0:
            totLL += log(Pdata)
        else:
            return WORST_LL
    return totLL


def H1sumLL(q, mu, sigma, a, genoData, phenos):
    muRR = mu
    muRA = mu + a
    muAA = mu + (2 * a)
    totLL = 0
    for j in range(len(genoData)):
        LRR, LRA, LAA = genoLs(genoData[j])
        if np.isnan(phenos[j]) == False:
            z = phenos[j]
            PdataRR = (1 / (np.sqrt(2 * sigma * np.pi))) * np.exp((-(z - muRR) ** 2) / (2 * sigma)) * LRR * q * q
            PdataRA = (1 / (np.sqrt(2 * sigma * np.pi))) * np.exp((-(z - muRA) ** 2) / (2 * sigma)) * LRA * 2 * q * (
                    1 - q)
            PdataAA = (1 / (np.sqrt(2 * sigma * np.pi))) * np.exp((-(z - muAA) ** 2) / (2 * sigma)) * LAA * (1 - q) * (
                    1 - q)
            Pdata = PdataRR + PdataRA + PdataAA
        else:
            Pdata = LRR * q * q + LRA * 2 * q * (1 - q) + LAA * (1 - q) * (1 - q)
        if Pdata > 0.0:
            totLL += np.log(Pdata)
        else:
            return WORST_LL
    return totLL


###############################################################

outfile = open(sys.argv[1], 'w') # outputs the LRT values (our method)
outfile.write('rep\tq\tmuRR\tsigma\ta\tLRT\n')

bimbam = open(sys.argv[2], 'w') #  the continuous value genotypes for gemma, in bimbam format
phenofile = open(sys.argv[3], 'w') # the phenotypic values for gemma


nsamples = 291
meandepth = 0.8384035
indivRange = range(9, nsamples + 9)
true_q = 0.3376 # observed allele frequency in amplicon step
nsparse = 10 # number  of sites with sparse effects on phenotype
n_other = 10000 # number of noncausal (neutral) sites. 
sm_val = 0.0001
WORST_LL = float('-inf')

#----------------------------------------
# Set up genotypes for first 10 genotypes

sparseGenoData = [[] for i in range(nsparse)]
qEst_vals = [] # proportion ref alleles observed for each indiv, based on read depth. used as a guess of q for optimizer

true_nref = [0 for i in range(ntot)] # number of ref alleles in each indiv. Used to assign phenos.
Vg = 0.0 # this gets incremented according to the inclusion of each sparse site
for i in range(nsparse):
    refCounter = 0
    dataCounter = 0
    for j in range(ntot):
        k = np.random.poisson(meandepth)
        if k > 0:
            tk = 1 - 2 * (0.5) ** k
            dataCounter += 1
        else:
            tk = 0

        x = random.random()
        if x <= true_q*true_q:
            true_nref[j] += 1
            if k > 0:
                obs = 'R'
                refCounter += 1.0
            else:
                obs = 'missing'

        elif x <= (true_q*true_q) + (2*(1-true_q)*true_q):
            if k > 0:
                x2 = random.random()
                if x2 <= (0.5)**k:
                    obs = 'R'
                    refCounter += 1.0
                elif x2 <= 2*(0.5)**k:
                    obs = 'A'
                else:
                    obs = 'H'
                    refCounter += 0.5

            else:
                obs = 'missing'

        else:
            true_nref[j] -= 1
            if k > 0:
                obs = 'A'
            else:
                obs = 'missing'

        sparseGenoData[i].append([tk, obs])
    qEst_vals.append(refCounter / float(dataCounter))
    Vg += 2*true_q*(1-true_q) # variance in pheno explained by all sites

# in this case, we want the proportion of total variance explained
# by our sparse sites to be 50% (each explains 5%)
Vp = Vg*2
Vr = Vp - Vg # residual variance not explained by site
#-----------------------------------------------------

# based on nref, assign phenotypes
phenos = []

for j in range(ntot):
    ph = random.normalvariate(true_nref[j], np.sqrt(Vr))
    x = phenos.append(ph)
    phenofile.write(str(ph) + '\n')
phenofile.close()


# extract some useful info for phenos that I'll need for finding assocs.
muMin = np.min(phenos)
muEst = np.mean(phenos)
muMax = np.max(phenos)
sigmaMin = 0.05 * np.var(phenos)
sigmaEst = np.var(phenos)
sigmaMax = 2 * np.var(phenos)
aMin = -(np.max(phenos) - np.min(phenos))
aEst = 0.0
aMax = np.max(phenos) - np.min(phenos)

qMin = 0.0 + sm_val
qMax = 1.0 - sm_val
#--------------------------------------------

# First need to run calculations for the sparse sites.

sitecounter = 0

for i in range(nsparse):
    genoData = sparseGenoData[i]
    bimbam.write(str(sitecounter) + ', T, A')
    qEst = qEst_vals[i]

    # find likelihood of null model
    model = 'H0'
    H0InitialValues = np.array([qEst, muEst, sigmaEst])
    H0Bounds = [(qMin, qMax), (muMin, muMax), (sigmaMin, sigmaMax)]
    H0MLEs, H0NegLL, d0 = optimize.fmin_l_bfgs_b(scipyLL, x0=H0InitialValues, bounds=H0Bounds, approx_grad=True)

    # find likelihood of alt model
    model = 'H1'
    H1InitialValues = np.array([qEst, muEst, sigmaEst, aEst])
    H1Bounds = [(qMin, qMax), (muMin, muMax), (sigmaMin, sigmaMax), (aMin, aMax)]
    H1MLEs, H1NegLL, d1 = optimize.fmin_l_bfgs_b(scipyLL, x0=H1InitialValues, bounds=H1Bounds, approx_grad=True)

    H0q = H0MLEs[0]
    H0mu = H0MLEs[1]
    H0sigma = H0MLEs[2]

    H1q = H1MLEs[0]
    H1muRR = H1MLEs[1]
    H1sigma = H1MLEs[2]
    H1a = H1MLEs[3]

    H0LL = -H0NegLL
    H1LL = -H1NegLL
    LRT = 2 * (H1LL - H0LL)
    outfile.write(str(sitecounter) + '\t' +
                  str(H1q) + '\t' +
                  str(H1muRR) + '\t' +
                  str(H1sigma) + '\t' +
                  str(H1a) + '\t' +
                  str(LRT) + '\n')

    for g in range(len(genoData)):
        LRR, LRA, LAA = genoLs(genoData[g])
        val = (1.0 * LRA * 2 * H0q * (1 - H0q) + 2.0 * LRR * H0q * H0q) / (
                (1 - H0q) * (1 - H0q) * LAA + LRA * 2 * H0q * (1 - H0q) + LRR * H0q * H0q)
        bimbam.write(', ' + str(val))
    bimbam.write('\n')
    sitecounter += 1


# Now simulate genotypes for the non-causal sites

for i in range(n_other):
    bimbam.write(str(sitecounter) + ', T, A')
    genoData = []

    refCounter = 0
    dataCounter = 0

    for j in range(ntot):
        k = np.random.poisson(meandepth)
        if k > 0:
            tk = 1 - 2 * (0.5) ** k
            dataCounter += 1
        else:
            tk = 0

        x = random.random()
        if x <= true_q*true_q:
            if k > 0:
                obs = 'R'
                refCounter += 1.0
            else:
                obs = 'missing'

        elif x <= (true_q*true_q) + (2*(1-true_q)*true_q):
            if k > 0:
                x2 = random.random()
                if x2 <= (0.5)**k:
                    obs = 'R'
                    refCounter += 1.0
                elif x2 <= 2*(0.5)**k:
                    obs = 'A'
                else:
                    obs = 'H'
                    refCounter += 0.5

            else:
                obs = 'missing'

        else:
            if k > 0:
                obs = 'A'
            else:
                obs = 'missing'

        genoData.append([tk, obs])
    qEst = refCounter/dataCounter

    # find likelihood of null model
    model = 'H0'
    H0InitialValues = np.array([qEst, muEst, sigmaEst])
    H0Bounds = [(qMin, qMax), (muMin, muMax), (sigmaMin, sigmaMax)]
    H0MLEs, H0NegLL, d0 = optimize.fmin_l_bfgs_b(scipyLL, x0=H0InitialValues, bounds=H0Bounds, approx_grad=True)

    # find likelihood of alt model
    model = 'H1'
    H1InitialValues = np.array([qEst, muEst, sigmaEst, aEst])
    H1Bounds = [(qMin, qMax), (muMin, muMax), (sigmaMin, sigmaMax), (aMin, aMax)]
    H1MLEs, H1NegLL, d1 = optimize.fmin_l_bfgs_b(scipyLL, x0=H1InitialValues, bounds=H1Bounds, approx_grad=True)

    H0q = H0MLEs[0]
    H0mu = H0MLEs[1]
    H0sigma = H0MLEs[2]

    H1q = H1MLEs[0]
    H1muRR = H1MLEs[1]
    H1sigma = H1MLEs[2]
    H1a = H1MLEs[3]

    H0LL = -H0NegLL
    H1LL = -H1NegLL
    LRT = 2 * (H1LL - H0LL)

    outfile.write(str(sitecounter) + '\t' +
                  str(H1q) + '\t' +
                  str(H1muRR) + '\t' +
                  str(H1sigma) + '\t' +
                  str(H1a) + '\t' +
                  str(LRT) + '\n')

    for g in range(len(genoData)):
        LRR, LRA, LAA = genoLs(genoData[g])
        val = (1.0 * LRA * 2 * H0q * (1 - H0q) + 2.0 * LRR * H0q * H0q) / (
                (1 - H0q) * (1 - H0q) * LAA + LRA * 2 * H0q * (1 - H0q) + LRR * H0q * H0q)
        bimbam.write(', ' + str(val))
    bimbam.write('\n')
    sitecounter += 1

outfile.close()
bimbam.close()












