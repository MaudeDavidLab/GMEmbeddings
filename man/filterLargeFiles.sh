##################
#
#       This script goes through all fastq files in current folder and checks 
#       if they contain over 5000 reads. If they do, fastq file is moved to 
#       the newly created folder "FilesWithOver5000Reads". The leftover files
#       in the current folder are bad (contain less than 5000 reads) and will
#       not be used in further analyses.
#
##################

#!/bin/bash

mkdir FilesWithOver5000Reads # folder to output "good" files to

#This for loop tells you how many reads are in each fastq file
#for file in *.fastq
#do
#       echo $file
#       echo $(cat $file|wc -l)/4|bc
#done


# loop over all fastq files in folder
for file in *.fastq
do
        # get number of reads in current file
        declare -i x="$(cat $file|wc -l)/4|bc"

        # checks if current files has over 5000 reads 
        # if true, prints file name and number of reads 
        # contained. It then moves the file to "good" folder
        if [ "$x" -gt "5000" ]
        then
                echo $file 
                echo $x
                mv $file FilesWithOver5000Reads
        fi
done
