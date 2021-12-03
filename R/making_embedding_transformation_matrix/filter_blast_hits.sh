#!/bin/bash

export TMPDIR=/nfs3/PHARM/David_Lab/christine/gut_microbiome_embeddings/
dir=$1
echo $dir/
rm $dir/best_hits.tsv
rm $dir/temp_file.tsv

unique_ASVs=$(cat $dir/blast_hits.tsv | awk '{print $1}' | sort -u -k1,1 --merge)
for id in $unique_ASVs;
do
	START1="$(date +%s)"
	echo "filtering $id"
	cat $dir/blast_hits.tsv | awk -v var="$id" '$1==var {print $0}'| sort -k5,5g  -k8,8nr -k7,7nr > $dir/temp_file.tsv
	evalue=$(awk  'FNR <= 1' $dir/temp_file.tsv | awk '{print $5}')
	pident=$(awk  'FNR <= 1' $dir/temp_file.tsv | awk '{print $8}')
	length=$(awk  'FNR <= 1' $dir/temp_file.tsv | awk '{print $7}')
	echo $evalue
	echo $pident
	echo $length	
	line=1
	while IFS="" read -r p || [ -n "$p" ]
	do
		evalue2=$(echo $p | awk '{print $5}')
		pident2=$(echo $p | awk '{print $8}')
		length2=$(echo $p | awk '{print $7}')
		if [[ "$evalue2" == "$evalue" && "$pident2" == "$pident" && "$length2" -eq "$length" ]]
		then
			echo $p >> $dir/best_hits.tsv
		else
			break
	#		echo "bad sequence"
	#		echo $evalue "vs" $evalue2
	#		echo $pident "vs" $pident2
	#		echo $length "vs" $length2
		fi
		line=$(($line + 1))
	done < $dir/temp_file.tsv
	rm -rf $dir/temp_file.tsv
	END1="$(date +%s)"
	DURATION1=$[${END1} - ${START1}]
	#echo "Filtering one ASV  took: ${DURATION1}"
 
done
