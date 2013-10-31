#!/usr/bin/python
__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'
"""This scripts takes the fasta output of mothur and removes the identifiers, 
reverse complements the DNA seqeunce and outputs the file with .rev extension. """
import re
import sys
from string import maketrans

if len(sys.argv) < 2:
  fasta_in = raw_input('Please input the sample name (ie. sd11.trim.fasta): ')
else:
  fasta_in = sys.argv[1]
rev_out = fasta_in + '.rev'
revcomp_key = maketrans('ACGTRYMKBDHVacgtrymkbdhv', 'TGCAYRKMVHDBtgcayrkmvhdn')
f_in = open(fasta_in, 'rU')
f_out = open(rev_out, 'w')
for line in f_in.readlines():
  if not re.match('^>', line):
    line = line[:-1] # remove new line character from end
    line = line[::-1]  #reverse sequence
    revcomp = line.translate(revcomp_key) + '\n' #complement sequence
    f_out.write(revcomp)
f_in.close()
print 'The sequences have been reverse complemented. The %s file has been created' % rev_out
f_out.close()


