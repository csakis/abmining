#!/usr/bin/python
__author__ = 'csaba'
__email__ = 'csakis[at]lanl[dot]gov'
import sys
if len(sys.argv) < 2:
  accum_in = raw_input('Please input the accumulation file name: ')
else:
  accum_in = sys.argv[1]

f_in = open(accum_in, 'rU')
resample = int(raw_input('Please enter the resampling depth (i.e. 100): '))
accum_out = '%s_%s.resample' % (accum_in, resample)
f_out = open(accum_out,'w')
for line in f_in.readlines()[::resample]:
  f_out.write(line)

