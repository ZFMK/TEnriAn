#!/bin/bash

## Modify for setup


my_ogs_version='version_1'
MY_ORTHO_SET='test_Orthoset'
MY_ORTHOGRAPH_DB='test_Orthograph_DB.sqlite'

CONTAINER=/PATH/TO/SINGULARITY-CONTAINERS/TEnriAn_container.sif

COG_file='resources/orthograph/lepi-tabfile-exons.txt'
INIT_config='resources/orthograph/orthograph_init.conf'

#Prepare OGS sequence files
#make sure header is matching the COG file

#create an initial config file
echo "created orthograph_init.conf file with the initial configurations for Orthograph"

cat <<EOF > ${INIT_config}
species-name = Initial_spec
input-file = resources/orthograph/Initialize_sets.fas
output-directory = Init_spec

database-backend   = sqlite
sqlite-program     = /usr/bin/sqlite3
sqlite-database    = $(pwd)/resources/orthograph/${MY_ORTHOGRAPH_DB}
sets-dir           = $(pwd)/resources/orthograph/sets
ortholog-set       = ${MY_ORTHO_SET}
make-set           = 1

alignment-program    = mafft-linsi
hmmbuild-program     = hmmbuild
makeblastdb-program  = makeblastdb
translate-program    = fastatranslate
hmmsearch-program    = hmmsearch
blast-program        = blastp
exonerate-program    = exonerate

EOF


printf "\nrunning orthograph-manager to setup the sqlite-DB\nName of the sqlite-DB is ${MY_ORTHOGRAPH_DB}\n\n"
printf "\n\n#####\n\nFor a more detailed descrition on how to setup orthograph please visit:\n\thttps://github.com/mptrsen/Orthograph/ \n\n#####\n\n"

printf "y\n" | singularity exec $CONTAINER orthograph-manager --create -c ${INIT_config}

#All OGS sequence files must be located in the folder INPUT_OGS, must be in fasta format, with .fas-ending and the file name must be the species name

echo "uploading all OGS files"

if [ -z "$(ls -A resources/orthograph/INPUT_OGS)" ]; then
	echo "the folder INPUT_OGS is empty. It should contain the OGS files...";
else
	echo "folder not empty..."

	for file in $(ls resources/orthograph/INPUT_OGS/*.fas) ; do 
		printf "\nusing file ${file}\n" ; 
		SPEC=$(basename ${file%.*} ) ; 
		echo "species name: ${SPEC} " ; 
	
		singularity exec $CONTAINER orthograph-manager --load-ogs-peptide ${file} --ogs-version ${my_ogs_version} --ogs-taxon-name ${SPEC} -c ${INIT_config} ; 	
	done
fi


printf "\nadding the tab delimited file with the COG information\n"
printf "################################\n"
printf "\nPlease enter the the requested information and confirm the correct assignment of your official gene sets\n\n"
printf "Please make sure that you provied the same set name that you provide here in the config file of the workflow!\n"
printf "Name of your Ortholog set :\t${MY_ORTHO_SET}\n"
printf "################################\n\n"

singularity exec $CONTAINER orthograph-manager ${COG_file} -c ${INIT_config}



printf "\n\nif the sets are not present they we will now start an initial run to produce them...\n\n"
printf "\n#####\nnow running an initial analysis on a very short set ... \n\tresources/orthograph/Initialize_sets.fas\n\nfeedback will be written to initial_orthograph_run.log\n#####\n\n"

singularity exec $CONTAINER orthograph-analyzer -c ${INIT_config} &> initial_orthograph_run.log

printf "\n\nfinished initial run. Now cleaning the intermediate output ... \n"

#clean up 

rm -r Init_spec

printf "DONE\n"

