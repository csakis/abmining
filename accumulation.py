#!/usr/bin/python
__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'
import sys
from blist import blist

"""This file estimates the library size using accumulation using the .cdr3 output of the pipeline"""

if len(sys.argv) < 2:
  cdr3_in = raw_input('Please input the file containing all the HCDR3s: ')
else:
  cdr3_in = sys.argv[1]
file_in = open(cdr3_in, 'rU')  # get the file output by cdr3_cluster
sample_name = cdr3_in.split('.')[0]
file_out = sample_name + '.accumulation' #set the output file name
f_out = open(file_out, 'w')
# the count for the actual CDR3 sample size
# shuffle three times just in case
cdr3_list = blist(file_in.readlines())
ordered_cdr3_list = sorted(cdr3_list, key = lambda str: len(str), reverse= True)
no_of_cdr3 = len(cdr3_list) #the number of CDR3s
cdr3_count = 0
unique_cdr3_no = 0
cdr3_dict = {} # the dictionary of all the CDR3s
unique_cdr3s = {} # the dictionary of unique CDR3 counts
for cdr3 in ordered_cdr3_list:
  if cdr3[:-1] in cdr3_dict: #check if the cdr3 is unique
    cdr3_dict[cdr3[:-1]] += 1
  else:
    cdr3_dict[cdr3[:-1]] = 1
    unique_cdr3_no += 1
  cdr3_count += 1
  if cdr3_count % 10000 == 0:
    print 'So far %d CDR3s have been checked.' % cdr3_count
  unique_cdr3s[cdr3_count] = unique_cdr3_no
for count, unique_no in sorted(unique_cdr3s.items()): #order the CDR3s by abundance
  f_out.write('%s, %s\n' % (count, unique_no)) #the csv file contains unique seq, count
f_out.close()
print 'The %s file has been created successfully.' % file_out