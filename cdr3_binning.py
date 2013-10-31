#!/usr/bin/python
__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'
""" This script's input is the file created by the cdr3search.py python script.
It searches for the unique cdr3s and count to occurrences of them.
It writes the data out into a csv file containing, the unique sequences and their counts."""
import sys
from operator import itemgetter

if len(sys.argv) < 2:
  cdr_in = raw_input('Please input the sample name sd6.cdr3: ')
else:
  cdr_in = sys.argv[1]
bins_out = 'bins_' + cdr_in
f_in = open(cdr_in, 'rU')
cdr3_dict = {}
for line in f_in.readlines(): #read CDR3s line by line
  line = line[:-1] #remove \n character from lines
  if line in cdr3_dict: #check if the cdr3 is unique
    cdr3_dict[line] += 1
  else:
    cdr3_dict[line] = 1
  if not len(cdr3_dict) % 1000:
    print '%d of CDR3s have been binned.' % len(cdr3_dict)
f_in.close()
f_out = open(bins_out, 'w')
for cdr, index in sorted(cdr3_dict.items(), key=itemgetter(1), reverse=True): #order the CDR3s by abundance
  f_out.write('%s, %s\n' % (cdr, index)) #the csv file contains unique seq, count
f_out.close()
print '\nThere are %d unique CDR3s found.' % len(cdr3_dict)
