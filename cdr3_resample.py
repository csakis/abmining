#!/usr/bin/python
__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'
import sys
import random
"""This file resamples the input file to the desired depth"""
if len(sys.argv) < 2:
  file_in = raw_input('Please input the bins file name : ')
else:
  file_in = sys.argv[1]
f_in = open(file_in, 'rU')
lines = f_in.readlines()
go_back = True
while   go_back:
  sampling = int(raw_input('Please enter the number of sequences you need: '))
  if sampling > len(lines):
    go_back = True
    print 'The number you entered is too big'
    print "There are %d sequences in the file" % len(lines)
  else:
    go_back = False

out_filename = file_in + '.sample'

f_out = open(out_filename, 'w')

lines[-1] = lines[-1] + '\n'
print 'shuffle 1'
random.shuffle(lines)
print 'shuffle 2'
random.shuffle(lines)
print 'shuffle 3'
random.shuffle(lines)
for item in lines[:sampling]:
  f_out.write(item)
print 'The %s file has been created successfully.' % out_filename