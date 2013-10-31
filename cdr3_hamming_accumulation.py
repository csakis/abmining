#!/usr/bin/python
import random
from blist import blist
from collections import defaultdict
import sys
import operator
from itertools import imap

__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'
"""This script's input (.cdr3) is the file created by the cdr3_pipeline.py python script, each line
containing a single CDR3 sequence.
It searches for the unique cdr3s and count to occurrences of them.
After that, the script bins the CDR3s at a given Hamming distance."""


def dist(str1, str2):
  ne = operator.ne
  hamming = sum(imap(ne, str1, str2))
  return hamming


if len(sys.argv) < 2:
  cdr_in = raw_input('Please input the sample name sd6.cdr3: ')
else:
  cdr_in = sys.argv[1]
bins_out = 'bins_' + cdr_in
f_in = open(cdr_in, 'rU')
cdr3_dict = {}
full_cdr_list_dict = defaultdict(list) # a dict containing all the CDR3s by length
cdr3_length_list = []
no_of_all_cdr3 = 0
for line in f_in.readlines(): #read CDR3s line by line
  line = line[:-1] #remove \n character from lines
  no_of_all_cdr3 +=1
  cdr3_length = len(line)
  full_cdr_list_dict[cdr3_length].append(line)
  if cdr3_length not in cdr3_length_list:
    cdr3_length_list.append(cdr3_length)
  if line in cdr3_dict: #check if the cdr3 is unique
    cdr3_dict[line] += 1
  else:
    cdr3_dict[line] = 1

print 'There are %d CDR3s to cluster.' % no_of_all_cdr3
print 'There are %d unique CDR3s in the file.' % len(cdr3_dict)

unique_cdr3_name_count = {}
unique_cdr3list_by_count = defaultdict(blist)
for cdr, index in sorted(cdr3_dict.items(), key=operator.itemgetter(1), reverse=True):
  unique_cdr3_name_count[cdr] = 0
  unique_cdr3list_by_count[len(cdr)].append(cdr) # unique cdr3 list descending count by length

Hamming_distance = int(raw_input('Please enter the desired Hamming distance: '))
sample_name = cdr_in.split('.')[0]
file_out_cluster_name = 'bins_Hamming_%d_%s.csv' % (Hamming_distance, sample_name) #file that contains the CDR3 clusters
file_out_accumulation_name = 'accumulation_Hamming_%d_%s.csv' % (Hamming_distance, sample_name)

hamming_accumulation_dict = {}  # contains the accumulation numbers at a given Hamming distance
hamming_accumulation_counter = 0
unique_no = 0
for cdr3_length in sorted(cdr3_length_list, reverse=True):
  all_cdr_list = blist(full_cdr_list_dict[cdr3_length])
  random.shuffle(all_cdr_list)
  unique_cdrlist = blist(unique_cdr3list_by_count[cdr3_length])
  for cdr3 in all_cdr_list:
    hamming_accumulation_counter += 1
    if hamming_accumulation_counter % 10000 == 0:
      print "%d CDR3s have been checked so far, %d to go." % (hamming_accumulation_counter, (no_of_all_cdr3 - hamming_accumulation_counter))
    if unique_cdr3_name_count[cdr3] == 0: # if we have not found this CDR3 yet
      for unique_cdr in unique_cdrlist:
        if dist(cdr3, unique_cdr) <= Hamming_distance:
          if unique_cdr3_name_count[unique_cdr] == 0:
            unique_no += 1
            unique_cdr3_name_count[unique_cdr] = 1
          else:
            unique_cdr3_name_count[unique_cdr] += 1
          break
    else:
      unique_cdr3_name_count[cdr3] += 1
    hamming_accumulation_dict[hamming_accumulation_counter] = unique_no
cdr3_count = 0
f_out = open(file_out_cluster_name, 'w')
for cdr3, no_of_cdr3 in sorted(unique_cdr3_name_count.items(), key=operator.itemgetter(1), reverse=True):
  if int(no_of_cdr3) > 0:
    cdr3_count += 1
    f_out.write('%s, %s\n' % (cdr3, no_of_cdr3)) #the csv file contains unique seq, count
f_out.close()

f_out = open(file_out_accumulation_name, 'w')
for cdr3, unique_no in sorted(hamming_accumulation_dict.items(), key=operator.itemgetter(1)):
  f_out.write('%s, %s\n' % (cdr3, unique_no)) #the csv file contains unique seq, count
f_out.close()
print 'There are %d unique clusters at Hamming distance %d' % (cdr3_count, Hamming_distance)
print 'The %s file has been created successfully.' % file_out_accumulation_name
print 'The %s file has been created successfully.' % file_out_cluster_name
