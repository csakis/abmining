#!/usr/bin/python
"""This script's inputs is the file created by cdr3_binning.py.
It clusters any number of samples containing unique CDR3 sequences and their counts.
First the user enters the number of tables and their names.
The script builds a dictionary containing all the sequences with their counts.
Then creates another dictionary that has the unique sequences as keys and the counts of the sequences as values in a tuple."""
import csv, sys

file_name_list = []
if len(sys.argv) < 2:  # check if the user entered any command arguments
  """If there's not command line argument, ask for the file names."""
  no_of_tables = int(raw_input('Please enter the number of files you would like to cluster? '))
  for i in range(no_of_tables):
    file_name_list.append(raw_input('Enter the name of file No: %s: ' % str(i + 1)))
else:
  """Get the file names from the command line arguments."""
  print "You have entered %d file names" % (len(sys.argv) - 1 )
  for i in range(1, len(sys.argv)):
    file_name_list.append(sys.argv[i])
out_file = raw_input('Please enter the cluster file name: ') + '.csv'
seq_dictionary = {}
for file_name in file_name_list:
  f = open(file_name, 'rU')  #open each sequence file
  cdr3_count = 0
  for line in f.readlines():
    cdr3_count += 1
    values = line.split(',')
    seq_dictionary[(file_name, values[0])] = values[1][
                                             :-1] #build a dictionary containing (name of table, sequence) = count of sequences
  print '%s table contains %d unique CDR3s.' % (file_name, cdr3_count)
#Build the cluster dictionary containing the unique sequence as the key, and a tuple containing the number of count is each table
cluster_dict = {}
for keys, count in seq_dictionary.items():
  table_name, uni_seq = keys
  temp_count = []
  if uni_seq not in cluster_dict:
    for table in file_name_list:
      if (table, uni_seq) in seq_dictionary:
        temp_count.append(int(seq_dictionary[(table, uni_seq)]))
      else:
        temp_count.append(0)
    cluster_dict[uni_seq] = temp_count
print 'The cluster contains %d sequences.' % len(cluster_dict)
file_out = csv.writer(open(out_file, 'wb'), delimiter=',')
file_name_list_2 =[]
for filename in file_name_list:
  if filename.startswith('bins_') == True and filename.endswith('.cdr3.csv') == True:
    file_name_list_2.append(filename[5:-9])
  else:
    file_name_list_2.append(filename)
file_name_list_2.insert(0, 'uni_seq')
file_out.writerow(file_name_list_2)
for seq, counts in cluster_dict.items():
  counts.insert(0, seq)
  file_out.writerow(counts)
