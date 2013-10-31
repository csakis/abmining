#!/usr/bin/python
__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'
import re
import sys

if len(sys.argv) < 2:
  file_in = raw_input('Please input the file containing the unique CDR3s')
else:
  file_in = sys.argv[1]
f_in = open(file_in, 'rU')
#creating a cdr3 dictionary.
#key is the length of the CDR3
#items the equal length CDR3s as a list
cdr3_dict = {}
for line in f_in.readlines():
  cdr3 = line.split(',')[0]
  cdr3_length = len(cdr3)
  if cdr3_length in cdr3_dict:
    cdr3_dict[cdr3_length].append(cdr3)
  else:
    cdr3_dict[cdr3_length] = [cdr3]
#cdr3 dictionary is built
f_in.close()
aa_dict = {}
for length, cdrlist in sorted(cdr3_dict.items()):
  for cdr3 in cdrlist:
    codon_list = re.findall('.', cdr3)
    position = 1
    for codon in codon_list:
      if (length, position, codon) in aa_dict:
        aa_dict[length, position, codon] += 1
      else:
        aa_dict[length, position, codon] = 1
      position += 1
#The aa dictionary is built
#Now we will write it into  file
file_out = file_in + '.aa'
f_out = open(file_out, 'w')
f_out.write('Length, Position, aa_codon, aa_count\n')
for key, aa_count in sorted(aa_dict.items()):
  length, position, codon = key
  text_out = ('%i, %i, %s, %i\n') % (length, position, codon, aa_count)
  f_out.write(text_out)
f_out.close()
print 'The %s file has been created successfully.' % file_out

