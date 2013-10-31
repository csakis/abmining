#!/usr/bin/python
__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'
""""""
import re
import sys
#from collections import defaultdict
if len(sys.argv) < 2:
  file_in = raw_input('Please input the sample name sd6.cdr3: ')
else:
  file_in = sys.argv[1]
file_out = file_in +'.out'
f_in = open(file_in, 'rU')
f_out = open(file_out, 'w')
aa_dict = {}
unique_list = []
for line in f_in:
  match = re.search('^>', line)
  if match:
    LCDR3 = line.split(';')[5]
  elif LCDR3 != '':
    extended_cdr = '.....' + LCDR3 + '.....'
    match2 = re.search(extended_cdr, line)
    if match2:
      if match2.group() not in unique_list:
        unique_list.append(match2.group())

for cdr in unique_list:
  for i in range(0,7):
    if (i, cdr[i]) in aa_dict:
      aa_dict[(i, cdr[i])] += 1
    else:
      aa_dict[(i, cdr[i])] = 1
    if (i+7, cdr[-(i+1)]) in aa_dict:
      aa_dict[(i+7, cdr[-(i+1)])] += 1
    else:
      aa_dict[(i+7, cdr[-(i+1)])] = 1

for key in sorted(aa_dict.keys()):
  position, aa = key
  count = aa_dict[key]
  print position, aa, count







