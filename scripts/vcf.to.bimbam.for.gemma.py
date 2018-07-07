import sys
import os
# import pandas as pd


vcf = open(sys.argv[1], 'rU')
AFs = open(sys.argv[2], 'rU')
tyx = open(sys.argv[3], 'rU')
bimbam = open(sys.argv[4], 'w')


def skipHeader(line):
    cols = line.replace('\n','').split('\t')
    if len(cols) < 2:
        return 'header'
    elif cols[0] == '#CHROM':
        return 'header'
    else:
        return cols

def extractVCFfields(sampleData):
    if sampleData != './.':
        fields = sampleData.split(':')
        alleleDepths = fields[1].split(',')
        totDepth = fields[2]
        phreds = fields[4].split(',')
        return [totDepth, alleleDepths, phreds]
    else:
        return 'missing'

def genoLs(annot):
    if annot == 'missing':
        LRR = 1.0
        LRA = 1.0
        LAA = 1.0
    else:
        readdepth = int(annot[0])
        phreds = annot[2]
        if readdepth <= max_kvalue:
            tk = tks[readdepth]
        else:
            tk = 1.0
        if phreds[0] == '0':
                LRR = 1.0
                LRA = (1.0-tk)/2.0
                LAA = 10**(-float(phreds[2])/10.0)
        elif phreds[1] == '0':
                LRR = 10**(-float(phreds[0])/10.0)
                LRA = 1.0
                LAA = 10**(-float(phreds[2])/10.0)
        elif phreds[2] == '0':
                LRR = 10**(-float(phreds[0])/10.0)
                LRA = (1.0-tk)/2.0
                LAA = 1.0
    return LRR, LRA, LAA

def calcGenoValue(annot, q):
    LRR, LRA, LAA = genoLs(annot)
    val = ( 1.0*LRA*2*q*(1-q) + 2.0*LRR*q*q)/((1-q)*(1-q)*LAA + LRA*2*q*(1-q) + LRR*q*q )
    return round(val, 4)


############################################################################

min_kvalue = 1
max_kvalue = 15

indivRange = range(9, 300)


tks = [0]
for idx1, line in enumerate(tyx):
	colx = line.replace('\n', '').split('\t') 
	if idx1>0:
		tks.append(float(colx[3]))

q_list = []
for site in AFs:
    cols = site.replace('\n', '').split('\t')
    if cols[0] != 'contig':
        q_list.append(round(float(cols[4]),4))

sitecounter = 0
for line in vcf:
    cols = skipHeader(line)

    if cols != 'header':
	if sitecounter % 1000 == 0:
		print sitecounter
        tig = cols[0]
        pos = cols[1]
        bimbam.write('{0}_{1}'.format(tig, pos) + ', T, A') # gemma does not use these dummy values

        q = q_list[sitecounter]

        for j in indivRange:

            annot = extractVCFfields(cols[j])
            val = calcGenoValue(annot, q)
            bimbam.write(', '+str(val))

        bimbam.write('\n')
        sitecounter += 1

bimbam.close()


