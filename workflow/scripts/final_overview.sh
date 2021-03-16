#!/bin/bash
#build a matrix which species is in which cluster...	

#usage:
# bash final_overview.sh -i input/dir/ -o output-file 

#set language so that grep is working...
LANG="en_IN.utf8"
export LANG

#get options

while getopts i:o: flag
do
    case "${flag}" in
        i) IN_DIR=${OPTARG};;
        o) OUT_FILE=${OPTARG};;
    esac
done


#initialize first line of outputfile with species names
printf '' >${OUT_FILE}
for SPEC in $(grep '>' ${IN_DIR%/}/*.fas  | cut -f2 -d':'  | cut -f1,2,3 -d'_' | sort | uniq | sed 's/>//') ; do 
printf "\t${SPEC}"  >>${OUT_FILE}
done
printf "\n" >>${OUT_FILE}

#go through each file
for FILE in $(ls ${IN_DIR%/}/*.fas );do 

printf "${FILE}\t">>${OUT_FILE} #print file name for first column

echo "${FILE}"

for SPEC in $(grep '>' ${IN_DIR%/}/*.fas  | cut -f2 -d':'  | cut -f1,2,3 -d'_' | sort | uniq | sed 's/>//') ; do 	
HIT=''
	HIT=$(grep "${SPEC}" ${FILE});
	if [[ ! -z "${HIT}" ]]; then
	printf "1\t">>${OUT_FILE}
	else printf "0\t">>${OUT_FILE}
	fi
	done
printf "\n">>${OUT_FILE}

done
