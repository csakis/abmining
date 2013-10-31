#!/usr/bin/python
__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'
import sys
from collections import defaultdict
from operator import itemgetter
"""This file searches for barcodes in a trimmed fasta file"""

if len(sys.argv) < 2:
  fasta_in = raw_input('Please input the fasta file: ')
else:
  fasta_in = sys.argv[1]
file_in = open(fasta_in, 'rU')  # get the file output by cdr3_cluster
sample_name = fasta_in.split('.')[0]
file_out = sample_name + '.barcodes' #set the output file name
f_out = open(file_out, 'w')
barcodes = defaultdict(int)
for line in file_in:
    seq = line[:10]
    barcodes[seq] +=1
for barcode, occurrence in sorted(barcodes.items(), key=itemgetter(1), reverse=True):
  f_out.write('%s, %i\n' % (barcode, occurrence))
f_out.close()