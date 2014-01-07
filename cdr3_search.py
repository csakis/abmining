#!/usr/bin/python
import os

__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'
"""This scripts takes the fasta output of mothur and removes the identifiers, 
reverse complements the DNA sequence and outputs the file with .rev extension.
It also finds CDR3s based on the cdr3_pattern regular expression pattern and translates
the found DNA seqeunces into peptide seqeunces and saves the list as .cdr file."""
import sys
import re
from string import maketrans

def protein_translate(aa_dict, seq):
  cdr3 = '' #empty protein sequence
  codon_list = re.findall('...', seq) #separates the sequence by triplets
  for codon in codon_list:
    cdr3 = cdr3 + aa_dict[codon] #use the triplet list to translat to aa
  return cdr3


def stars():
  print '\n', '*' * 80
  return


def reverse_complement(line):
  """
  This function reverse-complements a DNA sequence that contains a \n character at the end
  """
  revcomp_key = maketrans('ACGTRYMKBDHVacgtrymkbdhv', 'TGCAYRKMVHDBtgcayrkmvhdn') # the DNA complementation key
  line = line[:-1] # remove new line character from end
  line = line[::-1]  #reverse sequence
  rev_comp_line = line.translate(revcomp_key) + '\n' #complement sequence and add new line
  return rev_comp_line

#build the AA translation dictionary
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
if not os.path.isfile(os.path.join(__location__, 'aa_table.csv')):
  stars()
  print 'You need to download the aa_table.csv file from the ionpython website!'
  print 'http://sourceforge.net/ionpython/'
  stars()
  sys.exit(
    'Please run the program after you copied the aa_table.csv file into the directory where all the other python scripts are located!')
aa_file = open(os.path.join(__location__, 'aa_table.csv')) # the csv file containing the translation table
aa_dict = {}
for line in aa_file.readlines():
  aa = line.split(',')
  aa_dict[aa[1][:-1]] = aa[0]
#AA dictionary ready

"""Searching for CDR3s."""
if len(sys.argv) < 2:
  fasta_in = raw_input('Please input the sample file name (i.e.: sd6.trim.rev) ')
else:
  fasta_in = sys.argv[1]
cdr_out = fasta_in + '.cdr'

cdr_file_name = fasta_in + '.cdr3' #The file name that will contain the HCDR3 sequences.

# *** Here is the regular expression pattern. ***
cdr3_pattern = r'(TT[TC]|TA[CT])(TT[CT]|TA[TC]|CA[TC]|GT[AGCT]|TGG)(TG[TC])(([GA][AGCT])|TC)[AGCT]([ACGT]{3}){5,32}TGGG[GCT][GCT]'
#improved CDR3 regex pattern using the Quake dataset
stars()
print 'Now we are looking for CDR3 sequences'
seq_count = 0
cdr3_count = 0
fasta_file = open(fasta_in, 'rU')
fasta_list = fasta_file.readlines()
cdr_file_out = open(cdr_out, 'w')
#CDR3 search loop
for DNA_seq in fasta_list: #Go through each line of the trimmed fasta file
  seq_count += 1
  if re.match('^[ACGTN]', DNA_seq): #check if the line contains seq ID or DNA seq
    match = re.search(cdr3_pattern, reverse_complement(DNA_seq)) #check the reverse complement first
    if match:
      cdr3_count += 1
      cdr_file_out.write(protein_translate(aa_dict, match.group())[2:-1] + '\n') #write the CDR3 output file
    else:
      rev_match = re.search(cdr3_pattern, DNA_seq) #check if the seq contains a CDR3
      if rev_match:
        cdr3_count += 1
        cdr_file_out.write(protein_translate(aa_dict, rev_match.group())[2:-1] + '\n') #write the CDR3 output file
    if (seq_count % 250000) == 0:
      print 'So far %d sequences have been checked.' % seq_count
#End CDR3 search loop
stars()
print '%d CRD3s were found.' % cdr3_count
fasta_file.close() #close trimmed fasta file
cdr_file_out.close() #close cdr file
print 'The %s file has been created and contains all the CDR3s.' % cdr_file_name
stars()
# End of CDR3 search


