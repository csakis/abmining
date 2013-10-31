#!/usr/bin/python
__author__ = 'Csaba Kiss'
__email__ = 'csakis[at]lanl[dot]gov'
import sys
import os
import re
from string import maketrans
from operator import itemgetter
from Bio import SeqIO
import datetime


def stars():
  print '\n', '*' * 80
  return


def window(sequence, winSize, step):
  numOfChunks = ((len(sequence) - winSize) / step) + 1
  for i in range(0, numOfChunks * step, step):
    yield sequence[i:i + winSize]


def homopol(seq, homop):
  homo_regex = r'(A{%i,}|C{%i,}|G{%i,}|T{%i,})' % (homop, homop, homop, homop)
  return bool(re.search(homo_regex, seq))


def protein_translate(aa_dict, seq):
  """
  This function takes a DNA sequence and translates into a protein
  """
  cdr3 = '' #empty protein sequence
  codon_list = re.findall('...', seq)  # separates the sequence by triplets
  for codon in codon_list:
    cdr3 = cdr3 + aa_dict[codon]  # use the triplet list to translat to aa
  return cdr3


def reverse_complement(line):
  """
  This function reverse-complements a DNA sequence that contains a \n character at the end
  """
  revcomp_key = maketrans('ACGTRYMKBDHVacgtrymkbdhv', 'TGCAYRKMVHDBtgcayrkmvhdn') # the DNA complementation key
  line = line[:-1]  # remove new line character from end
  line = line[::-1]  # reverse sequence
  rev_comp_line = line.translate(revcomp_key) + '\n' #complement sequence and add new line
  return rev_comp_line


def barcode_search(cdr3, barcodes):
  barcode_found = False
  for barcode_seq, barcode_sample_name in barcodes:
    barcode_pattern = r'^' + barcode_seq
    match = re.search(barcode_pattern, cdr3) #check the reverse complement first
    if match:
      barcode_found = True
      break
  return barcode_found, cdr3, barcode_sample_name

# Logging all screen outputs to a log file

class Logger(object):
  def __init__(self):
    log_file_name = 'barcode_' + sample_name + '_' + datetime.datetime.now().strftime('%m-%d-%Y-%H-%M') + '.log'
    self.terminal = sys.stdout
    self.log = open(log_file_name, "w")

  def write(self, message):
    self.terminal.write(message)
    self.log.write(message)


#-----------------------------------------------------------------------------------------------------------------------
#Main program starts here.
stars()
print '\t\t*** Script Description ***\n'
print 'This script uses raw data (sff file) from 454 or Ion torrent sequncers or\n' \
      'fastq files from 454, Ion Torrent or MiSeq (fastq extracted using Casava 1.8).\n' \
      'After entering the quality trimming requirements, the script finds the HCDR3 sequences in the\n' \
      'DNA sequences and then translates them. Finally it bins them according to their abundance.\n' \
      'The script produces three or (four) output files for each barcode.\n' \
      '(If the input was the raw sff file, it creates the fastq file).\n' \
      'The trim.fasta file contains all the quality filtered DNA files.\n' \
      'The .cdr3 file contains all the HCDR3 sequences.\n' \
      'The .dna.cdr3 file contains all the DNA HCDR3 sequences.\n' \
      'The  bins...cdr3.csv file contains the HCDR3 sequences with their abundances.\n' \
      'The  bins...dna.cdr3.csv file contains the DNA sequences\n encoding the HCDR3 peptides with their abundances.\n'
stars()
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

while_condition = True
while while_condition: #sequence file check
  if len(sys.argv) < 2:
    raw_file_name = raw_input('Please enter the sff or fastq file name (ie. sd6.sff or sd6.fastq): ')
  else:
    raw_file_name = sys.argv[1]
  sample_name = raw_file_name.split('.')[0]
  file_extension = raw_file_name.split('.')[1]
  if len(raw_file_name.split('.')) != 2 or (file_extension != 'sff'
                                            and file_extension != 'fastq') or not os.path.isfile(raw_file_name):
    print 'You entered an incorrect file name or the file does not exist.'
    print 'Please enter the correct file name!'
    sys.exit()
  else:
    while_condition = False

sys.stdout = Logger()

#*** Barcode checking ***
while_condition = True
while while_condition:
  while_condition = False
  barcode_file_name = raw_input('Please enter the barcode file name: ')
  if not os.path.isfile(barcode_file_name):
    print 'The file does not exit. Try again'
    while_condition = True
  else:
    barcode_file = open(barcode_file_name, 'rU')
    barcodes = [] # the list containing the barcode sequences and names in a tuple format (seq, sample_name)
    for line in barcode_file.readlines():
      barcode_list = line.split('\t')
      if len(barcode_list) == 3 and barcode_list[0] == 'barcode':
        barcodes.append((barcode_list[1], barcode_list[2][:-1]))
      else:
        print 'The barcode file is not the correct format.'
        print 'One line should look like the oligo format mothur uses:'
        print 'barcode[TAB]CTAAAGGTAAC[TAB]sample_name'
        print 'Please load the correct file!'
        while_condition = True
        break
# *** Barcodes loaded successfully ***
# *** they are stored in the barcodes list as (seq, samplename)


"""Sequence trimming settings"""

#minimum length
while_condition = True
while while_condition:
  stars()
  min_length = raw_input("Please enter the minimum required DNA length (Default=50): ")
  if not min_length:
    min_length = 50
    while_condition = False
  else:
    min_length = int(min_length)
  if min_length < 1 or min_length > 300:
    print 'You entered a wrong value.'
    print 'Try again!'
  else:
    while_condition = False

#Qwindow size
while_condition = True
while while_condition:
  stars()
  qwindowsize = raw_input('Please enter the quality window size (Default=50): ')
  if not qwindowsize:
    qwindowsize = 50
    while_condition = False
  else:
    qwindowsize = int(qwindowsize)
  if qwindowsize > min_length or qwindowsize < 2:
    print 'You entered a wrong value.'
    print 'Try again!'
  else:
    while_condition = False

# Max homopolymer setting
while_condition = True
while while_condition:
  stars()
  homop = raw_input('Please enter the max homopolymer value (Default=8): ')
  if not homop:
    homop = 8
    while_condition = False
  else:
    homop = int(homop)
  if homop > 20 or homop < 2:
    print 'You entered a wrong value.'
    print 'Try again!'
  else:
    while_condition = False

#window step setting
while_condition = True
stars()
print 'The step stetting influences the speed and sensitivity quality trimming.'
print 'The value of 1 is the slowest and most accurate. If speed is an issue,'
print 'it is recommended to experiment with settings up to 10.'
while while_condition:
  step = raw_input('Please enter the step value (Default=1): ')
  if not step:
    step = 1
    while_condition = False
  else:
    step = int(step)
  if step > 10 or step < 1:
    print 'You entered a wrong value.'
    print 'Try again!'
  else:
    while_condition = False

#Quality average setting
while_condition = True
stars()
print 'This is the most important setting for quality trimming.'
print 'For Ion Torrent reads the optimal value is 20-22.'
print 'For 454 sequencing 18 is a good compromise.'
print 'Too low of a value will increase sequencing read errors.'
print 'Too high of a value will eliminate too many sequences.'
while while_condition:
  qwindow_average = raw_input('Please enter the quality average value (Default=20): ')
  if not qwindow_average:
    qwindow_average = 20
    while_condition = False
  else:
    qwindow_average = int(qwindow_average)
  if qwindow_average > 40 or qwindow_average < 0:
    print 'You entered a wrong value.'
    print 'Try again!'
  else:
    while_condition = False

#sff sequence extraction
if file_extension == 'sff':  # extract the sequences from the sff file
  fastq_file = sample_name + '.fastq'
  stars()
  print 'The sff file is now being converted into a fastq file.'
  print 'This could take a while...'
  sff_seq_count = SeqIO.convert(raw_file_name, "sff-trim", fastq_file, "fastq")
  print '%i sequences have been converted to fastq format.' % sff_seq_count
  print '%s file has been created' % fastq_file
  stars()

if file_extension == 'fastq':
  fastq_file = raw_file_name

#*** DNA sequence quality trimming ***
stars()
print 'Now we trim the sequences using the desired quality settings.'
print 'Depending on the settings, this could take a while... (10-20 minutes)'
print 'Please be patient!'
stars()
good_reads = []
good_seq_count = 0
for rec in SeqIO.parse(fastq_file, 'fastq'):
  seq_length = len(rec.seq) # length of the raw sequence
  DNA_seq = str(rec.seq) # the raw DNA sequence
  quality_list = rec.letter_annotations['phred_quality']  #raw quality list
  keep_seq = True

  #trimming the left and right ends
  left_trim = 0
  right_trim = seq_length - 1
  while quality_list[right_trim] < qwindow_average and right_trim > 0:
    right_trim -= 1
  while quality_list[left_trim] < qwindow_average and left_trim < seq_length - 1:
    left_trim += 1

  # check if sequence is still longer than minimum legth required.
  if (seq_length - (left_trim + seq_length - right_trim)) <= min_length:
    keep_seq = False

  # check for homopolymers
  if keep_seq:
    trimmed_seq = DNA_seq[left_trim:right_trim]
    if homopol(trimmed_seq, homop):
      keep_seq = False

  # check for quality
  if keep_seq:
    current_qual_list = rec.letter_annotations["phred_quality"][left_trim:right_trim]
    for win in window(current_qual_list, qwindowsize, step):
      current_qual = sum(win) / qwindowsize
      if current_qual < qwindow_average:
        keep_seq = False
        break
  if keep_seq:
    good_seq_count += 1
    if (good_seq_count % 100000) == 0:
      print 'So far %d good sequences have been found.' % good_seq_count
    good_reads.append(trimmed_seq)  # good_reads list contains all the quality sequences

#Searching for barcodes
barcode_dict = {}  # the barcode dictionary will contain the barcode separated sequence_lists usign the barcode as key
for barcode_seq, barcode_name in barcodes:
  barcode_dict[barcode_name] = []

for cdr3_sequence in good_reads:
  barcode_found, cdr3_found, barcode_sample_name = barcode_search(cdr3_sequence, barcodes)
  if barcode_found:
    barcode_dict[barcode_sample_name].append(
      cdr3_found)  # add the found cdr3 sequence to the appropriate list in the dictionary
# barcode library is done

for barcode, seq_list in barcode_dict.items():
  trim_fasta_file_name = barcode + '_' + sample_name + '.trim.fasta'
  trim_fasta_file = open(trim_fasta_file_name, 'w')
  for cdr3 in seq_list:
    trim_fasta_file.write('%s\n' % cdr3)
  stars()
  print '*** Sample: %s ***' % barcode
  print 'The %s file has been created and contains %d trimmed DNA sequences' % (trim_fasta_file_name, len(seq_list))
  """Searching for CDR3s."""
  cdr_file_name = barcode + '_' + sample_name + '.cdr3' #The file name that will contain the HCDR3 sequences.
  cdr_dna_file_name = barcode + '_' + sample_name + '.dna.cdr3' #The file name that will contain the HCDR3 sequences.

  # *** Here is the regular expression pattern. ***
  cdr3_pattern = r'TA[CT](TT[CT]|TA[TC]|CA[TC]|GT[AGCT])TG[TC][GA][AGCT]([ACGT]{3}){5,32}[AGCT]TGGG[GCT][GCT]' #CDR3 search pattern
  # cdr3_pattern = r'(TT[TC]|TA[CT])(TT[CT]|TA[TC]|CA[TC]|GT[AGCT]|TGG)(TG[TC])(([GA][AGCT])|TC)[AGCT]([ACGT]{3}){5,32}TGGG[GCT][GCT]' #CDR3 search pattern
  if not os.path.isfile(trim_fasta_file_name):
    sys.exit('Please check that you have the trimmed fasta file in the directory')
  print 'Now we are looking for CDR3 sequences in the %s sample' % barcode
  seq_count = 0
  cdr3_count = 0
  fasta_file = open(trim_fasta_file_name, 'rU')
  fasta_list = fasta_file.readlines()
  cdr_file_out = open(cdr_file_name, 'w')
  cdr_dna_file_out = open(cdr_dna_file_name, 'w')
  #CDR3 search loop
  for DNA_seq in fasta_list: #Go through each line of the trimmed fasta file
    seq_count += 1
    if re.match('^[ACGTN]', DNA_seq): #check if the line contains seq ID or DNA seq
      match = re.search(cdr3_pattern, reverse_complement(DNA_seq)) #check the reverse complement first
      if match:
        cdr3_count += 1
        cdr3_peptide = protein_translate(aa_dict, match.group())[2:-1]
        cdr_file_out.write(cdr3_peptide + '\n')  # write the CDR3 output file
        cdr_dna_text = '%s, %s\n' % (cdr3_peptide, match.group())
        cdr_dna_file_out.write(cdr_dna_text)  # write the CDR3 DNA output file
      else:
        rev_match = re.search(cdr3_pattern, DNA_seq)  # check if the seq contains a CDR3
        if rev_match:
          cdr3_count += 1
          cdr3_peptide = protein_translate(aa_dict, rev_match.group())[2:-1]
          cdr_file_out.write(cdr3_peptide + '\n')  # write the CDR3 output file
          cdr_dna_text = '%s, %s\n' % (cdr3_peptide, rev_match.group())
          cdr_dna_file_out.write(cdr_dna_text)  # write the CDR3 DNA output file
      if (seq_count % 250000) == 0:
        print 'So far %d sequences have been checked.' % seq_count
        #End CDR3 search loop
  print '\n%d CDR3s were found in the %s sample.' % (cdr3_count, barcode)
  fasta_file.close() #close trimmed fasta file
  cdr_file_out.close() #close cdr file
  cdr_dna_file_out.close() #close cdr DNA file
  print 'The %s file has been created and contains all the CDR3s.' % cdr_file_name
  # End of CDR3 search

  # CDR3 clustering begins here
  print '\nNow we will find each unique CDR3s and their abundance.'
  bins_out = 'bins_' + barcode + '_' + sample_name + '.cdr3.csv'
  cdr_in = open(cdr_file_name, 'rU')
  cdr_list = cdr_in.readlines()
  cdr3_count = 0
  cdr3_dict = {}
  for line in cdr_list: #read CDR3s line by line
    cdr3_count += 1
    line = line[:-1] #remove \n character from lines
    if line in cdr3_dict: #check if the cdr3 is unique
      cdr3_dict[line] += 1
    else:
      cdr3_dict[line] = 1
    if not cdr3_count % 100000:
      print '%d of CDR3s have been searched so far.' % cdr3_count
  cdr_in.close()
  print 'There are %d unique CDR3s found in the %s sample.' % ((len(cdr3_dict)), barcode)
  f_out = open(bins_out, 'w')
  for cdr, index in sorted(cdr3_dict.items(), key=itemgetter(1), reverse=True):  # order the CDR3s by abundance
    f_out.write('%s, %s\n' % (cdr, index))  # the csv file contains unique seq, count
  f_out.close()
  #end clustering

  # CDR3 DNA clustering begins here
  print '\nNow we will find each unique CDR3 DNA sequences and their abundance.'
  bins_dna_out = 'bins_' + barcode + '_' + sample_name + '.dna.cdr3.csv'
  cdr_dna_in = open(cdr_dna_file_name, 'rU')
  cdr_list = cdr_dna_in.readlines()
  cdr3_count = 0
  cdr3_dict = {}
  for line in cdr_list: # read CDR3s line by line
    cdr3_count += 1
    line = line[:-1] # remove \n character from lines
    if line in cdr3_dict: # check if the cdr3 is unique
      cdr3_dict[line] += 1
    else:
      cdr3_dict[line] = 1
    if not cdr3_count % 100000:
      print '%d of CDR3s have been searched so far.' % cdr3_count
  cdr_dna_in.close()
  print 'There are %d unique CDR3 DNAs found in the %s sample.' % ((len(cdr3_dict)), barcode)
  f_out = open(bins_dna_out, 'w')
  for cdr, index in sorted(cdr3_dict.items(), key=itemgetter(1), reverse=True):  # order the CDR3s by abundance
    f_out.write('%s, %s\n' % (cdr, index))  # the csv file contains unique seq, count
  f_out.close()
  #end CDR3 DNA clustering

stars()
print 'Sequence processing has been finished.'
print '\t\t *** Summary of the results ***\n'
if file_extension == 'sff':
  print 'The %s file has been created containing %i sequences.' % (fastq_file, sff_seq_count)

for barcode_seq, barcode_name in barcodes:
  trim_fasta_file_name = barcode_name + '_' + sample_name + '.trim.fasta'
  cdr_file_name = barcode_name + '_' + sample_name + '.cdr3'
  cdr3_count = 0
  cdr3_in = open(cdr_file_name, 'rU')
  for line in cdr3_in.readlines():
    cdr3_count += 1
  cdr3_in.close()
  bins_out = 'bins_' + barcode_name + '_' + sample_name + '.cdr3.csv'
  unique_count = 0
  bins_in = open(bins_out, 'rU')
  for line in bins_in.readlines():
    unique_count += 1
  bins_in.close()
  bins_dna_out = 'bins_' + barcode_name + '_' + sample_name + '.dna.cdr3.csv'
  unique_dna_count = 0
  bins_dna_in = open(bins_dna_out, 'rU')
  for line in bins_dna_in.readlines():
    unique_dna_count += 1
  bins_dna_in.close()
  print '*** Sample: %s ***' % barcode_name
  print 'The %s file has been created containing %i trimmed sequences.' % (
    trim_fasta_file_name, len(barcode_dict[barcode_name]))
  print 'The %s file has been created containing %i CDR3 sequences.' % (cdr_file_name, cdr3_count)
  print 'The %s file has been created containing %i unique CDR3 peptide sequences.' % (bins_out, unique_count)
  print 'The %s file has been created containing %i unique CDR3 DNA sequences.' % (bins_dna_out, unique_dna_count)
  stars()
