#!/usr/bin/python
__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'

import sys
from string import maketrans

def reverse_comp(line, revcomp_key):
  line = line[:-1]  # remove new line character from end
  line = line[::-1]  # reverse sequence
  revcomp = line.translate(revcomp_key)  # complement sequence
  return(revcomp)
"""This script takes two paired end MiSeq fasta files and combines them into one paired end read file.
The script takes two command arguments: the first file should be the 5'-reads and the second argument is the fasta
file containing the '3-end reads. The second file's reads will be reverse-complemented and added to the first file's reads."""

if len(sys.argv) < 3:
  file_in_1 = raw_input('Please enter the name of the fasta file containing the forward reads: ')
  file_in_2 = raw_input('Please enter the name of the fasta file containing the reverse reads: ')
else:
  file_in_1 = sys.argv[1]
  file_in_2 = sys.argv[2]
f_in_1 = open(file_in_1,'rU')
f_in_2 = open(file_in_2,'rU')
f_out_name = '%s_%s_paired_end.fasta' % (file_in_1, file_in_2)
f_out = open(f_out_name, 'w')
revcomp_key = maketrans('ACGTRYMKBDHVacgtrymkbdhv', 'TGCAYRKMVHDBtgcayrkmvhdn')
file_1_lines = f_in_1.readlines()
file_2_lines = f_in_2.readlines()
for i in range(len(file_1_lines)):
  if i%10000 == 0:
    print 'So far %d sequences have been processed.' % (i / 2)
  if i%2 == 0:
    f_out.write(file_1_lines[i])
  if i%2 == 1:
    sequence_1 = file_1_lines[i][:-1]
    sequence_2 = reverse_comp(file_2_lines[i], revcomp_key)
    full_seq = '%s%s\n' % (sequence_1, sequence_2)
    f_out.write(full_seq)
print 'The %s file has been created successfully.' % f_out_name



