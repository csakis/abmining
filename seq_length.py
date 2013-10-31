#!/usr/bin/python
__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'
"""This script uses a trimmed fasta file and counts lengths of the sequences. """
import sys
from collections import defaultdict
from operator import itemgetter

if len(sys.argv) < 2:
    file_in = raw_input('Please input the fasta file name: ')
else:
    file_in = sys.argv[1]
f_in = open(file_in, 'rU')
file_lengths = defaultdict(int)
for line in f_in:
    file_lengths[len(line)] += 1
for length, number in sorted(file_lengths.items(), key=itemgetter(1), reverse=True):  # order the CDR3s by abundance
    print 'Sequence length: %d \tOccurrences: %d' % (length, number) # the csv file contains unique seq, count