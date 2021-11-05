##########
# 
# PROGRAM: make_GloVe_input.py
# DESCRIPTION: Program reads in the filtered final sequence table and creates an ASV co-occurence file.
#               Each line in the ASV co-occurence file contains the full length ASV sequence of all ASVs 
#               in one sample. 
# INPUT: Program reads in the filtered final sequence table that was produced in the DADA2 pipeline. The 
#               sequence table contains 50425 columns (ASVs) and 15706 row entries (samples).
# OUTPUT: Program outputs a GloVe input file. This is a an ASV co-occurrence file: each line in the 
#               file contains the full length ASV sequence of all ASVs in one sample. The final file 
#               contained 15,706 lines, one for every sample.  
#        
# NOTE: When running on the command line this script must be run using "python3" not "python"
#
##########

#!/usr/bin/env python

##---------------------------
## import necessary libraries
##---------------------------
import csv
import numpy as np
import random
import pandas as pd

##--------------------------------
## Read in file and open output file
##---------------------------------
data_dir = "/nfs3/PHARM/David_Lab/austin_e/filtered_final_seqtab/"
f = open(data_dir + "V2_transposed_filtered_final_seqtab.csv", 'r')
outfile = open("/nfs3/PHARM/David_Lab/austin_e/filtered_final_seqtab/glove_input.txt", mode = 'w')
print("V2_transposed_filtered_final_seqtab.csv")
writer = csv.writer(outfile, delimiter = "\t", quoting = csv.QUOTE_NONE, escapechar = '')

##-------------------------------------
## Get list of all taxa (ASV sequences)
##-------------------------------------
taxa_names = f.readline()
taxa_names = taxa_names.strip().split(",")
taxa_names = taxa_names[1:]


i = 0 # set counter
#random.seed(0)
##---------------------------------------------------------------------------------------
## Loop through each line of the sequence table (not including the 1st line/column names)
##---------------------------------------------------------------------------------------
for line in f:
        vec = line.split(",")
        sample_id = vec[0]

        present = [float(i) > 0 for i in vec[1:]] # get list of columns in the current row that have non-zero entry
        writer.writerow(np.array(taxa_names)[present]) # add column names (ASVs) that are present to output file 
        print(i, end = '\t')
        i = i + 1

##-------------
## close files
##-------------
f.close()
outfile.close()
