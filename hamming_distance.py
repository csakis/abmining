#!/usr/bin/python

__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'

""" This script calculates the Hammming distances of a dataset of CDR3s. First it creates subsets of CDR3s of equal length, then
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
  file_in = raw_input('Please enter the name of the file containing the CDR3s: ')
else:
  file_in = sys.argv[1]

cdr3_file = open(file_in, 'rU')
file_out = 'hamming_' + file_in  #the file containing the hamming distances
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
for length in cdrlendict.keys():  #the calculation is done for each lengths of CDR3s.
  cdrlist = cdrdict[length] #the cdrlist contains all the CDR3s of equal lengths
  length_of_cdrlist = len(cdrlist)
  hammingdict = {} #The library contains the number of occurences of Hamming distances
  for i in range(0, length_of_cdrlist):
    for j in range(i + 1, length_of_cdrlist):
      hamming_dist = dist(cdrlist[i], cdrlist[j])
      if hamming_dist in hammingdict:
        hammingdict[hamming_dist] += 1
      else:
        hammingdict[hamming_dist] = 1

  #writing the results into a file
  text_out = '\nHamming distance of %d length\n' % length
  dist_file.write(text_out)
  for  ind, dis in sorted(hammingdict.items()):
    text_out = '%d, %d, %d\n' % (length, ind, dis)
    dist_file.write(text_out)
