#!/usr/bin/python
import sys

__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'
if len(sys.argv) < 2:  # check if the user entered any command arguments
  """If there's not command line argument, ask for the file names."""
  bin_file = raw_input('Please enter the name of the file containing unique CDR3s (ie. bins_sd6.cdr3.csv): ')
else:
  bin_file = sys.argv[1]
cdr3_file = open(bin_file, 'rU')
file_out = 'length_' + bin_file
f_out = open(file_out, 'w')
cdrlendict = {}
for line in cdr3_file:
  cdr = line.split(',')[0]
  if len(cdr) in cdrlendict:
    cdrlendict[len(cdr)] += 1
  else:
    cdrlendict[len(cdr)] = 1
cdr3_file.close()
for length, number in sorted(cdrlendict.items()):
  text_out = '%d, %d\n' % (length, number)
  f_out.write(text_out)
f_out.close()
print 'The %s file has been created successfully.' % file_out