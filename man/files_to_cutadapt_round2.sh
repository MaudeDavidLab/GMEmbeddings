##################
#
# This script performs the second round of cutadapt to the fastq files
#       Here we just remove the last 5 bases from each read 
#
# NOTE: This could have been added to the first cutadapt file, but we didn't realize until later
#
##################

#!/bin/bash

# loops through all fastq files in the folder
for file in *.fastq
do
        echo cutadapt -u -5 -o trimmed_"$file" $file

        ##### CUTADAPT parameters: #####
        # -u  --> removes bases from each read 
        #       - here we remove the last 5 bases from each read by specifying it as negative (-5)
        #
        # -o --> output file name to write to
        #       - here we just append "trimmed_" to the beginning of the original file name

done
