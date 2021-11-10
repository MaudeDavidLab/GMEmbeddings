##################
#
# This script filters BLAST ouput files based on evalue(low), pident(high), and length(high)
#
# INPUT: BLAST output file
# OUTPUT: BLAST output file containing only the 'best' hits
#
##################

#!/bin/bash

## Find out how many different unique sequences are in the BLAST output (filters the 'qseqid' column)
max_lines=$(cat Schirmer_blast_hits.tsv | awk '{print $1}' | sort -u -k1,1 --merge | wc -l)

echo "number of different sequences:" $max_lines

## Loop through the unique sequences in the 'qseqid' column
for n in $( eval echo {1..$max_lines} )
do
        ## Sort the current unique 'qseqid':
        ## 1. make empty file to move all BLAST output matching current unique 'qseqid'
        ## 2. move all current 'qseqid' matches to the new file 'temp_file.tsv' 
        echo "filtering seq$n"
        cat Schirmer_blast_hits.tsv | grep -w "seq$n" >> temp_file.tsv
        LINES=$(cat temp_file.tsv)

        ## Get the filtering parameters of the current unique 'qseqid'
        ## 1. Get lowest evalue
        ## 2. Get largest pident (percentage of identical matches)
        ## 3. Get largest length (alignment length (sequence overlap))
        evalue=$(cat temp_file.tsv | sort -k5,5gr --merge | sed -n 1p | awk '{print $5}')
        pident=$(cat temp_file.tsv | sort -k8,8n --merge | sed -n 1p | awk '{print $8}')
        length=$(cat temp_file.tsv | sort -k7,7n --merge | sed -n 1p | awk '{print $7}')

        ## set counter
        line=1

        ## While loop to go though the newly made file
        ##      It loops through each line of in the file
        ## the parameters in this while loop allow us to keep leading whitespace, 
        ##      doesn't interpret backslash sequences, or skip the last line if it's 
        ##      missing a terminating linefeed
        while IFS="" read -r p || [ -n "$p" ] 
        do
                ## Get the parameters of the current line:
                ## 1. evalue
                ## 2. pident
                ## 3. length
                evalue2=$(echo $p | awk '{print $5}')
                pident2=$(echo $p | awk '{print $8}')
                length2=$(echo $p | awk '{print $7}')

                ## Check if the parameters on the current line match what was found at the beginning of 
                ## the for loop
                ##      If parameters match, pipe the BLAST hit to output file
                if [[ "$evalue2" == "$evalue" && "$pident2" == "$pident" && "$length2" -eq "$length" ]]
                then
                        echo $p >> Schirmer_filtered_best_hits.tsv
        #       else
        #               echo "bad sequence"
        #               echo $evalue "vs" $evalue2
        #               echo $pident "vs" $pident2
        #               echo $length "vs" $length2
                fi
                line=$(($line + 1))
        done <temp_file.tsv

        ## delete the temporary file
        rm -rf temp_file.tsv
done
