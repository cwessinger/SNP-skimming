"""
author: Carrie Wessinger

Usage: python pyLDx.phase1.py input.vcf input.sam output.haplos.txt

Generates a file listing all instances of VCF sites that occur on same read
"""

import re
import sys
import os

vcf = open(sys.argv[1], 'rU')
sam = open(sys.argv[2], 'rU')
outfile = open(sys.argv[3], 'w')

###########################################################

### First fill in a dictionary of snps of interest. 
### that has structure snpDict[contig][position] = [RefBase, AltBase]

snpDict = {}
for vline in vcf:
    vcols = vline.replace('\n', '').split('\t')
    snpTig = vcols[0]
    snpPos = vcols[1]
    snpRef = vcols[3]
    snpAlt = vcols[4]
    if snpDict.has_key(snpTig):
        snpDict[snpTig][snpPos] = [snpRef, snpAlt]
    else:
        snpDict[snpTig] = {}
        snpDict[snpTig][snpPos] = [snpRef, snpAlt]

#------------------------------------------------------------------------------------------

cigar_regex = re.compile(r'(\d+[MIDSHP])') # regex search for operations within cigar
cigarpart_regex = re.compile(r'(\d+)([MIDSHP])') # regex search to split up individual cigar ops

for line in sam:
    
    cols = line.replace('\n', '').split('\t')
    contig = cols[2]
    start = int(cols[3])
    cigar = cols[5]
    SEQ = cols[9]
    QUAL = cols[10]

    # use cigar string to create a final read sequence that can be related to the reference coordinates.
    
    cigarOps = cigar_regex.findall(cigar) # divides cigar string into components
    readsite = 0 # keeps track of where we are on the read (SEQ string)
    readseq = '' # string of base calls for read
    readqual = '' # string of quality scores for read

    for opp in cigarOps:
        opLength = int(cigarpart_regex.search(opp).group(1)) # the number part
        opCode = cigarpart_regex.search(opp).group(2) # the letter part

        if opCode == 'H': # if the operation is a hard clip, don't record
            pass

        elif opCode == 'I' or opCode == 'S': # if the operation is an insertion or soft clip
            for i in range(opLength):
                readsite += 1 # don't record information to readseq or readqual.

        if opCode == 'M': # if the operation is a matching alignment
            for i in range(opLength):
                readseq = readseq + SEQ[readsite] 
                readqual = readqual + QUAL[readsite]
                readsite += 1

        elif opCode == 'D': # if the operation is a deletion
            for i in range(opLength):
                readseq = readseq + '.' # record blanks to readseq and readqual at this location
                readqual = readqual + '.'

    ##### Now we look for the significant SNPs
    LocusList = []
    AlleleList = []
    
    for j in range(len(readseq)):
        currSite = start + j
        
        try:
            refBase = snpDict[contig][str(currSite)][0]
            altBase = snpDict[contig][str(currSite)][1]
            
            if readseq[j] == refBase and (ord(readqual[j])-33) >= 20:
                LocusList.append(currSite)
                AlleleList.append('R')
                
            elif readseq[j] == altBase and (ord(readqual[j])-33) >= 20:
                LocusList.append(currSite)
                AlleleList.append('A')

        except KeyError:
            pass
        
    if len(LocusList) >= 1:
        outfile.write(contig)
        for k in range(len(LocusList)):
            outfile.write('\t'+str(LocusList[k])+'\t'+str(AlleleList[k]))
        outfile.write('\n')
outfile.close()