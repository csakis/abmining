#!/usr/bin/python
__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'
"""This script uses the cdr3.dna file and search for the DNA sequences for the CDR3s provided by a text file containing the CDR3 peptide seqeunces. """
from collections import defaultdict
from operator import itemgetter
import sys

if len(sys.argv) < 2:
  file_in = raw_input('Please input the dna.cdr3 file name: ')
else:
  file_in = sys.argv[1]
f_in = open(file_in, 'rU')
sample_name = ''.join(file_in.split('.')[:-2])
pattern_file_name = raw_input('Please input the CDR3 pattern file name: ')
pattern_file = open(pattern_file_name, 'rU')
for cdr3_pattern in pattern_file:  # the for loop going through each CDR3 in the pattern file
  cdr3_pattern = cdr3_pattern[:-1]  # remove new line from the end
  print cdr3_pattern
  file_out_name =  sample_name + cdr3_pattern + '.dna.txt'
  file_out = open(file_out_name, 'w')
  cdr3_dna_list = []
  unique_cdr3_dict = defaultdict(int)
  for dna_cdr3_line in f_in:
    cdr_protein, cdr3_dna = dna_cdr3_line.split(',')
    if cdr_protein == cdr3_pattern:
      unique_cdr3_dict[cdr3_dna[:-1]] += 1
  for cdr, index in sorted(unique_cdr3_dict.items(), key=itemgetter(1), reverse=True):  # order the CDR3s by abundance
    file_out.write('%s, %s\n' % (cdr, index)) # the csv file contains unique seq, count
  file_out.close()
  f_in.seek(0)







