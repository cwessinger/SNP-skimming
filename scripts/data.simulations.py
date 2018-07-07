"""
author: Carrie Wessinger

Usage: python data.simulations.py outfile.txt

Simulates shallow sequence data for a locus with a true effect on
phenotype and estimates phenotypic effects observed from dataset.

"""
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
    """ function to minimize using optimize.brent """
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
    """ assign genotype likelihoods """
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
    """ sum log likelihood under null model that locus does not affect phenotype """
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
    """ sum log likelihood under alt model that locus affects phenotype """
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

outfile = open(sys.argv[1], 'w')
outfile.write('rep\tq\tmuRR\tsigma\ta\tLRT\n')

nreps = 100 # number of replicate simulations
nsamples = 291 # number of individuals sampled from the population
meandepth = 0.8384035 # mean sequence read depth modeled

true_q = 0.3376 # population allele frequency of the simulated locus
true_a = 0.2916 # additive effect of the simulated locus (effect of alt allele)
true_muRR = 20.611 #  mean phenotype of homozygous (ref allele) individuals
Vp = 2.2082 # phenotypic variance

sm_val = 0.0001
WORST_LL = float('-inf')

Vg = 2*true_q*(1-true_q)*true_a*true_a # variance in phenotype explained by site
Vr = Vp - Vg # residual variance not explained by site

for r in range(nreps):

    # set up lists of read depths drawn from Poisson,
    # 'true' genotypes (drawn from HWE)
    #  observed data (based on true genotypes and drawn read depths
    #  and simulated phenotypes for each individual

    depths = []
    trueGenos = []
    obsGenos = []
    phenos = []

    for i in range(nsamples):
        k = np.random.poisson(meandepth)
        depths.append(k)
        x = random.random()
        if x <= true_q*true_q:
            trueGenos.append('R')
            if k > 0:
                obsGenos.append('R')
            else:
                obsGenos.append('missing')
            phenos.append(random.normalvariate(true_muRR, np.sqrt(Vr)))

        elif x <= (true_q*true_q) + (2*(1-true_q)*true_q):
            trueGenos.append('H')
            if k > 0:
                x2 = random.random()
                if x2 <= (0.5)**k:
                    obsGenos.append('R')
                elif x2 <= 2*(0.5)**k:
                    obsGenos.append('A')
                else:
                    obsGenos.append('H')
            else:
                obsGenos.append('missing')
            phenos.append(random.normalvariate(true_muRR + true_a, np.sqrt(Vr)))

        else:
            trueGenos.append('A')
            if k > 0:
                obsGenos.append('A')
            else:
                obsGenos.append('missing')
            phenos.append(random.normalvariate(true_muRR + 2 * true_a, np.sqrt(Vr)))

    # extract some useful info for phenos that I'll need for finding assocs.
    muMin = np.min(phenos)
    muEst = np.mean(phenos)
    muMax = np.max(phenos)
    sigmaMin = 0.05*np.var(phenos)
    sigmaEst = np.var(phenos)
    sigmaMax = 2*np.var(phenos)
    aMin = -(np.max(phenos) - np.min(phenos))
    aEst = 0.0
    aMax = np.max(phenos) - np.min(phenos)

    # make lists of genotype and phenotype for random subsamples
    genoData = []
    refCounter = 0.0
    dataCounter = 0
    for i in range(nsamples):
        k = depths[i]
        if k > 0:
            dataCounter += 1
            tk = 1-2*(0.5)**k
            if obsGenos[i] == 'R':
                refCounter += 1.0
            elif obsGenos[i] == 'H':
                refCounter += 0.5
        else:
            tk = 0
        genoData.append([tk, obsGenos[i]])

    qMin = 0.0 + sm_val
    qEst = refCounter / dataCounter
    qMax = 1.0 - sm_val

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
    outfile.write(str(r) + '\t' +
                  str(H1q) + '\t' +
                  str(H1muRR) + '\t' +
                  str(H1sigma) + '\t' +
                  str(H1a) + '\t' +
                  str(LRT) + '\n')

outfile.close()











