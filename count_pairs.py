#!/usr/bin/python
__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'
import sys

"""This file counts the occurrences of CDR3s in two populations.
It need a file that was clustered using cdr3_cluster.py.
The file should contain a cluster of two sequencing outputs."""

if len(sys.argv) < 2:
  cluster_in = raw_input('Please input the cluster name : ')
else:
  cluster_in = sys.argv[1]
sample_name = cluster_in.split('.')[0]
file_in = open(cluster_in, 'rU')  # get the file output by cdr3_cluster
file_out = sample_name + '.count' #set the output file name
count_dict = {} # the dictionary contains the count pairs
line_count = 0
for line in file_in.readlines():
  if line_count != 0:
    value_list = line.split(',') #split the count values and the unique sequence
    count1 = int(value_list[1])
    count2 = int(value_list[2][:-1])
    if (count1, count2) in count_dict:
      count_dict[count1, count2] += 1
    else:
      count_dict[count1, count2] = 1
  else:
    value_list = line.split(',')
    sample1 = value_list[1]
    sample2 = value_list[2][:-1]
  line_count += 1
file_in.close()
f_out = open(file_out, 'w')
f_out.write('%s, %s, Occurrence\n' % (sample1, sample2))
for counts, occurrence in sorted(count_dict.items()):
  count1, count2 = counts
  f_out.write('%i, %i, %i\n' % (count1, count2, occurrence))
f_out.close()
print 'The %s file has been created successfully.' % file_out