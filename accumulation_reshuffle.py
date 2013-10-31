#!/usr/bin/python
from blist import blist
from collections import defaultdict
import random
"""This script recreates the "random" distribution of accumulation of a Hamming accumulation set"""

__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'
import sys
if len(sys.argv) < 2:
  accum_in = raw_input('Please input the accumulation file name: ')
else:
  accum_in = sys.argv[1]
cdr_in = raw_input('Please enter the cdr3 file name: ')
f_cdr_in = open(cdr_in, 'rU')

cdr_length_dict = defaultdict(int)  # a dict the number of CDR3 for a defined length
cdr_length_list = []  # contains all the lengths of CDR3s
for line in f_cdr_in.readlines():  #read CDR3s line by line
  line = line[:-1] #remove \n character from lines
  cdr_length = len(line)
  if cdr_length not in cdr_length_list:
    cdr_length_list.append(cdr_length)
  cdr_length_dict[cdr_length] +=1
f_cdr_in.close()
cdr_length_list.sort(reverse=True)
print cdr_length_list
f_in = open(accum_in, 'rU')
accum_out = '%s_.reshuffle' % (accum_in)
f_out = open(accum_out,'w')
cdr3_list = f_in.readlines()
list_length = len(cdr3_list)
increment_list =blist() # the increment list contains the derivative of each neighboring CDR3
increment_list = [1]
for i in range(0,list_length-1):
  read1, unique1 = cdr3_list[i].split(',')
  read2, unique2 = cdr3_list[i+1].split(',')
  if (int(unique2[:-1]) - int(unique1[:-1])) == 0:
    increment_list.append(0)
  else:
    increment_list.append(1)

cdr_dict_by_length = defaultdict(blist)
"""this dictionary contains the ordered derivative list for each CDR3 lengths"""

for cdr_length in cdr_length_list:
  cdr_dict_by_length[cdr_length] = increment_list[0: cdr_length_dict[cdr_length]]
  del increment_list[0:cdr_length_dict[cdr_length]]

"""The accumulation list has been chopped up into pieces for each CDR3 lengths"""


"""Now we have to create a new randomized new accumulation list.
I'll create a random CDR3 length and pop the first number from the cdr_dict_by_length
and add it to the new accumulation list."""
new_accumulation_list = blist()
while len(cdr_length_list) > 0:
  random_cdr_length = cdr_length_list[random.randrange(0, len(cdr_length_list))]
  new_accumulation_list.append(cdr_dict_by_length[random_cdr_length].pop())
  if len(cdr_dict_by_length[random_cdr_length]) == 0:
    cdr_length_list.remove(random_cdr_length)
new_accumulation_list.reverse()
# new_accumulation list is done

count = 1
accumulation = 1
for number in new_accumulation_list:
  f_out.write(('%d, %d\n') % (count, accumulation))
  count +=1
  accumulation += number

