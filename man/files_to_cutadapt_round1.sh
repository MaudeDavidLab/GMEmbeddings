##################
#
# This script performs the first round of cutadapt to the fastq files
#       Here we just remove the specified adapters from each read 
#
##################

#!/bin/bash

cd all_numbered_fastq_files
for file in *.fastq
do
        echo cutadapt -b GTGYCAGCMGCCGCGGTAA -b GGACTACNVGGGTWTCTAAT -b GTGCCAGCMGCCGCGGTAA -b GGACTACHVGGGTWTCTAAT -b AATGATACGGCGACCACCGAGATCTACACGCT -b CAAGCAGAAGACGGCATACGAGAT -o CA_"$file" $file

        ##### CUTADAPT parameters: #####
        # -b --> removes the specified adapter from each read
        #
        # -o --> output file name to write to
        #       - here we just append "CA_" to the beginning of the original file name

done
