#!/usr/bin/python
__author__ = 'csaba'
__email__ = 'csakis[at]lanl[dot]gov'
import sys
import os

"""This script calculates the mean hydophobicity of CDR3sw. It uses the unique CDR3 file (bins_*.cdr3) output.
The calculation is based on the CCS scale (consensus hydrophobicity scale)
[A. Tossi, L. Sandri, A. Giangaspero (2002) New consensus hydrophobicity scale extended to non-proteinogenic amino acids.
In Peptides 2002: Proceedings of the twenty-seventh European peptide symposium. Edizioni Ziino, Napoli, Italy. pp. 416-417]"""

# Building the hydrophobicity library
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
if not os.path.isfile(os.path.join(__location__, 'hydro.csv')):
  print 'You need to download the hydro.csv file from the ionpython website!'
  print 'http://sourceforge.net/projects/ionpython/'
  sys.exit(
    'Please run the program after you copied the hydro.csv file into the directory where all the other python scripts are located!')
hydro_file = open(os.path.join(__location__, 'hydro.csv')) # the csv file containing the translation table
hydro_dict = {}
for line in hydro_file.readlines():
  aa = line.split(',')
  hydro_dict[aa[0]] = float(aa[1][:-1])
# hydrophobicity is done.

if len(sys.argv) < 2:
  file_in = raw_input('Please enter the name of the file containing the CDR3s: ')
else:
  file_in = sys.argv[1]
sample_name = file_in.split('.')[0][5:]
file_out_name = sample_name + '.hydro'
f_in = open(file_in, 'rU')
f_out = open(file_out_name,'w')
cdr3_counter = 0
for line in f_in.readlines():  # read CDR3s line by line
  cdr3_counter +=1
  if cdr3_counter % 10000 == 0:
    print 'So far %d CDR3s have been checked.' % cdr3_counter
  line = line[:-1]  # remove \n character from lines
  line = line.split(',')
  hydro = 0
  for char in line[0]:
    hydro += hydro_dict[char]
  hydro = hydro / len(line[0])
  f_out.write('%s, %s, %.2f\n' % (line[0], line[1], hydro))
f_in.close()
f_out.close()
print 'The %s file has been created successfully.' % file_out_name



