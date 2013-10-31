#!/usr/bin/python


__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'

""" The script needs the output of theDNA bin output of the cdr3_pipeline.py as the input.
This script calculates the minimal Hamming distances of a dataset of unique CDR3s. First it creates subsets of CDR3s of equal length, then
calculates the Hamming distances each possible pairs. The output of the script is a csv files containing the number of occurrences of Hamming distances
for each CDR3 length.
In order to make this script run as fast as possible, nothing is written on the screen when it's running.
These calculations can be computationally demanding and depending on the size of the dataset, the script could run for days."""

from itertools import imap
import operator
import sys

def dist(str1, str2):
  ne = operator.ne
  hamming = sum(imap(ne, str1, str2))
  return hamming


if len(sys.argv) < 2:
  file_in = raw_input('Please enter the name of the file containing the unique CDR3s: ')
else:
  file_in = sys.argv[1]

cdr3_file = open(file_in, 'rU')
file_out = 'minimal_hamming_dna.' + file_in  #the file containing the hamming distances
dist_file = open(file_out, 'w')
cdrdict = {} #the dictionary that contains the list of CDR3s grouped by their lengths as the key.
cdrlendict = {} # the library of CDR3 lengths

# building the CDR3 libraries
for line in cdr3_file:
  cdr = line.split(',')[0]
  if len(cdr) in cdrdict:
    cdrdict[len(cdr)].append(cdr)
  else:
    cdrdict[len(cdr)] = [cdr]
  if len(cdr) in cdrlendict:
    cdrlendict[len(cdr)] += 1
  else:
    cdrlendict[len(cdr)] = 1
cdr3_file.close()
#CDR3 libraries done.

#Calculating the Hamming distances
hammingdict = {} #The library contains the composite Hamming distances of all CDR3s
hammingdict[1] = 0
count_check = 0
for length in cdrlendict.keys():  #the claculation is done for each lengths of CDR3s.
  cdrlist = cdrdict[length] #the cdrlist contains all the CDR3s of equal lengths
  length_of_cdrlist = len(cdrlist)
  for i in range(0, length_of_cdrlist - 1):
    min_hamming_dist = 1000 # just an arbitrary big number that we start with
    for j in range(i + 1, length_of_cdrlist):
      count_check += 1
      if count_check % 50000000 == 0:
        print '%d cdr3s have been checked' % count_check
      hamming_dist = dist(cdrlist[i], cdrlist[j])
      if hamming_dist < min_hamming_dist:
        min_hamming_dist = hamming_dist
      if min_hamming_dist == 1:
        hammingdict[1] += 1
        break
    if min_hamming_dist > 1:
      if min_hamming_dist in hammingdict:
        hammingdict[min_hamming_dist] += 1
      else:
        hammingdict[min_hamming_dist] = 1

#writing the results into a file
text_out = 'Hamming Distance, Frequency\n'
dist_file.write(text_out)
for ind, dis in sorted(hammingdict.items()):
  text_out = '%d, %d\n' % ( ind, dis)
  dist_file.write(text_out)
