#!bin/bash

# This script changes headers of trinity assembly files (or other fasta files) to include the expected number of reads (from RSEM output) as needed for ContamCheck. Each assembly is found in a different folder, and the script runs across all subfolders in the directory you are currently in.

# Requires Seqkit, which can be installed using Conda


curdir=$(pwd)
for folder in ./*; do 
  [ -d "$folder" ] && cd "$folder" && 

	# Get headers from RSEM output, remove header line, and sort:
	cat RSEM.isoforms.results | cut -f 1 | tr -s '\t' '_' | tail -n +2 | sort  > RSEM_headers


	# sort assembly file
	for file in *.fas; do
 	 seqkit sort -2 -n $file > $file.sorted;
	done

	#get headers from assembly fasta file without >
	for file in *.sorted; do
	  grep '>' $file | sed 's/>//' > assmbl_headers;
	done

	# Check if headers from RSEM output and assembly are identical and in the same order
	DIFF=$(diff RSEM_headers assmbl_headers)
	if [ "$DIFF" = "" ] 
	then 
	    # Combine headers and expected counts, remove header of resulting file and sort
	    cat RSEM.isoforms.results | cut -f 1,4 | tr -s '\t' '_' | tail -n +2 | sort > headers_with_expected_count

	    #Make key value file
	    paste assmbl_headers headers_with_expected_count > key_values

	    # Replace header with key value
	       for file in *.sorted; do
	         seqkit replace -p '(.+)' --replacement '{kv}' --kv-file key_values $file > $file.renamed.fasta;
	       done
	else 
	    echo 'Error: Order of header files in RSEM output and assembly do not match' > ERROR.txt;
	fi

done
cd $curdir

