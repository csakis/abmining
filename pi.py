#!/usr/bin/python
__author__ = 'csaba'
__email__ = 'csakis[at]lanl[dot]gov'
import sys
import os
from collections import defaultdict


"""This script calculates the mean pI of CDR3s. It uses the unique CDR3 file (bins_*.cdr3) output as the input file.
The calculation is based on

    Bjellqvist, B.,Hughes, G.J., Pasquali, Ch., Paquet, N., Ravier, F.,
    Sanchez, J.-Ch., Frutiger, S. & Hochstrasser, D.F. The focusing
    positions of polypeptides in immobilized pH gradients can be predicted
    from their amino acid sequences. Electrophoresis 1993, 14, 1023-1031.

    MEDLINE: 8125050 """

# Building the pI library
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
if not os.path.isfile(os.path.join(__location__, 'pi.csv')):
  print 'You need to download the pi.csv file from the ionpython website!'
  print 'http://sourceforge.net/projects/ionpython/'
  sys.exit(
    'Please run the program after you copied the hydro.csv file into the directory where all the other python scripts are located!')
pi_file = open(os.path.join(__location__, 'pi.csv'))  # the csv file containing the translation table
pi_dict = {}
for line in pi_file.readlines():
  line = line[:-1]  # removing \n
  pi_list = line.split(',')
  pi_dict[pi_list[0]] = [float(pi_list[1]), float(pi_list[2]), float(pi_list[3])]
# the  pI library is done.

if len(sys.argv) < 2:
  file_in= 'test.txt'
  # file_in = raw_input('Please enter the name of the file containing the CDR3s: ')
else:
  file_in = sys.argv[1]
sample_name = file_in.split('.')[0][5:]
file_out_name = sample_name + '.pi'
f_in = open(file_in, 'rU')
f_out = open(file_out_name,'w')
cdr3_counter = 0
for line in f_in.readlines():  # read CDR3s line by line
  cdr3_counter +=1
  if cdr3_counter % 10000 == 0:
    print 'So far %d CDR3s have been checked.' % cdr3_counter
  line = line[:-1]  # remove \n character from lines
  line = line.split(',')
  cdr3_seq = line[0]  # The particular CDR3
  cdr3_count = line[1]  # The count of the particular CDR3
  nterm = cdr3_seq[0]
  cterm = cdr3_seq[-1]

  aa_count_dict = defaultdict(int)  # this dictionary holds the number occurrences of each aa in the CDR3
  for aa in cdr3_seq[1:-1]:
    aa_count_dict[aa] += 1
  # pI calculation begins here
  pHmin = 0.0 # the lowest pH
  pHmax = 14.0  # the highest pH
  maxiteration = 2000 # max iteration value
  epsi = 0.0001  # desired precision
  i = 0  # the iteration counter
  charge = 1
  while (i < maxiteration and (pHmax - pHmin > epsi)):
    pHmid = pHmin + (pHmax - pHmin) / float(2)
    pHmid_exp = 10**(-pHmid)
    cter = 10**(-pi_dict[cterm][0]) / (10**(-pi_dict[cterm][0]) + pHmid_exp)
    nter = pHmid_exp / (10**(-pi_dict[nterm][1]) + pHmid_exp)
    carg = chis = clys = casp = cglu = ccys = ctyr = 0
    if aa_count_dict['R']>0:
      carg = aa_count_dict['R'] * pHmid_exp / (10**(-pi_dict['R'][2]) + pHmid_exp)
    if aa_count_dict['H']>0:
      chis = aa_count_dict['H'] * pHmid_exp / (10**(-pi_dict['H'][2]) + pHmid_exp)
    if aa_count_dict['K']>0:
      clys = aa_count_dict['K'] * pHmid_exp / (10**(-pi_dict['K'][2]) + pHmid_exp)
    if aa_count_dict['D']>0:
      casp = aa_count_dict['D'] * 10**(-pi_dict['D'][2]) / (10**(-pi_dict['D'][2]) + pHmid_exp)
    if aa_count_dict['E']>0:
      cglu = aa_count_dict['E'] * 10**(-pi_dict['E'][2]) / (10**(-pi_dict['E'][2]) + pHmid_exp)
    if aa_count_dict['C']>0:
      ccys = aa_count_dict['C'] * 10**(-pi_dict['C'][2]) / (10**(-pi_dict['C'][2]) + pHmid_exp)
    if aa_count_dict['Y']>0:
      ctyr = aa_count_dict['Y'] * 10**(-pi_dict['Y'][2]) / (10**(-pi_dict['Y'][2]) + pHmid_exp)


    charge = carg + clys + chis + nter - (casp + cglu + ctyr + ccys + cter)  # charge at pHmid
    if charge > 0:
      pHmin = pHmid
    else:
      pHmax = pHmid
    i += 1

  f_out.write('%s, %s, %.3f\n' % (line[0], line[1], pHmid))
f_in.close()
f_out.close()
print 'The %s file has been created successfully.' % file_out_name



