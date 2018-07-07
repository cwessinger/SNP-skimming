"""
author: Carrie Wessinger

Usage: python divideVCF.py infile.vcf nlines_per_file

Splits a given VCF file into multiple "subVCF" files that each contain
'nlines_per_file' lines.
"""

import sys
import os

infile = open(sys.argv[1], 'rU')

divideby = int(sys.argv[2]) # to divide into files containing this number of lines
count = 0
for idx, line in enumerate(infile):

    if idx == 0:
        outfile = open('./subVCF_0.txt', 'w')
        outfile.write(line)

    else:

        if idx % divideby != 0:
            outfile.write(line)

        else:
            count += 1
            outfile.close()
            outfile = open('./subVCF_{0}.txt'.format(count), 'w')
            outfile.write(line)

outfile.close()