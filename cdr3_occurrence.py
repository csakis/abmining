#!/usr/bin/python
__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'
import sys

"""This file counts the occurrences of CDR3s in a dataset. It needs the binned CDR3 file and puts out the occurrence file
that contains the number of occurrences of CDR3s."""

if len(sys.argv) < 2:
  bin_in = raw_input('Please input the binned CDR3 file name: ')
else:
  bin_in = sys.argv[1]
file_in = open(bin_in, 'rU')  # get the file output by cdr3_binning
file_out = bin_in + '.occurrence' #set the output file name
print 'We will go through each CDR3 and its frequencies and count how many times it is found.'
counter = 0
count_dict = {} # the dictionary contains the count pairs
for line in file_in.readlines():
  counter += 1
  if counter % 1000 == 0:
    print 'So far %d CDR3s have been checked.' % counter
  value_list = line.split(',') #split the count values and the unique sequence
  cdr3_count = int(value_list[1])
  if cdr3_count in count_dict:
    count_dict[cdr3_count] += 1
  else:
    count_dict[cdr3_count] = 1
file_in.close()
f_out = open(file_out, 'w')
for cdr3_count, occurrence in sorted(count_dict.items()):
  f_out.write('%i, %i\n' % (cdr3_count, occurrence))
f_out.close()
print 'The %s file has been created successfully.' % file_out


