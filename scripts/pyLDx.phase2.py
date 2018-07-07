"""
author: Carrie Wessinger

Usage: python pyLDx.phase2.py input.haplos.txt output.LDresults.txt

Compile haplotype information from phase1 to find LD between pairs of
sites sequenced on the same reads.
"""
import sys
import os

infile = open(sys.argv[1], 'rU')
outfile = open(sys.argv[2], 'w')
outfile.write('contig\t' + 'site1\t' +
              'site2\t' +
              'distance\t' +
              'freq_site1\t' +
              'freq_site2\t' +
              'xRR\t' +
              'xRA\t' +
              'xAR\t' +
              'xAA\tr2\n')

siteDict = {}
pairDict = {}

for idx, line in enumerate(infile):

    cols = line.replace('\n', '').split('\t')
    currcontig = cols[0]

    if idx == 0:
        pass

    elif currcontig != prevcontig:

        for i in pairDict:
            sites = i.split('_')
            site1 = sites[0]
            site2 = sites[1]
            distance = int(site2) - int(site1)

            xRR = int(pairDict[i][0])
            xRA = int(pairDict[i][1])
            xAR = int(pairDict[i][2])
            xAA = int(pairDict[i][3])

            psite1 = float(xRR + xRA) / float(xRR + xRA + xAR + xAA)
            psite2 = float(xRR + xAR) / float(xRR + xRA + xAR + xAA)
            if psite1 == 0.0 or psite1 == 1.0 or psite2 == 0.0 or psite2 == 1.0:
                pass

            else:

                r2 = (((xRR / float(xRR + xRA + xAR + xAA)) - (psite1 * psite2)) ** 2) / (
                        psite1 * psite2 * (1 - psite1) * (1 - psite2))

                outfile.write(str(prevcontig) + '\t' + str(site1) + '\t' +
                              str(site2) + '\t' +
                              str(distance) + '\t' +
                              str(psite1) + '\t' +
                              str(psite2) + '\t' +
                              str(xRR) + '\t' +
                              str(xRA) + '\t' +
                              str(xAR) + '\t' +
                              str(xAA) + '\t' +
                              str(r2) + '\n')

        siteDict = {}
        pairDict = {}

    else:
        pass

    numSites = (len(cols) - 1) / 2
    positionList = []
    callList = []

    # first get the lists for positions and calls.
    for i in range(numSites):
        positionList.append(cols[i * 2 + 1])
        callList.append(cols[i * 2 + 2])

    for j in range(len(positionList)):

        position = positionList[j]
        call = callList[j]

        if siteDict.has_key(position) == False:
            siteDict[position] = [0, 0]
        if call == 'R':
            siteDict[position][0] += 1
        elif call == 'A':
            siteDict[position][1] += 1

    for k in range(len(positionList)):

        step = 1

        while (k + step) <= range(len(positionList))[-1]:

            if pairDict.has_key('{0}_{1}'.format(positionList[k], positionList[k + step])) == False:
                pairDict['{0}_{1}'.format(positionList[k], positionList[k + step])] = [0, 0, 0, 0]
            if callList[k] == 'R' and callList[k + step] == 'R':
                pairDict['{0}_{1}'.format(positionList[k], positionList[k + step])][0] += 1
            elif callList[k] == 'R' and callList[k + step] == 'A':
                pairDict['{0}_{1}'.format(positionList[k], positionList[k + step])][1] += 1
            elif callList[k] == 'A' and callList[k + step] == 'R':
                pairDict['{0}_{1}'.format(positionList[k], positionList[k + step])][2] += 1
            if callList[k] == 'A' and callList[k + step] == 'A':
                pairDict['{0}_{1}'.format(positionList[k], positionList[k + step])][3] += 1

            step += 1

    prevcontig = currcontig

outfile.close()
