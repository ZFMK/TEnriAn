# The TEnriAn (Target ENrichment ANalysis) workflow


In this workflow, Target Enrichment (TE) data is assembled, assigned to orthologous clusters, aligned against reference Hidden Markov Models (HMMs) and filtered to produce Multiple Sequences Alignments (MSAs). The MSAs can be used for phylogenetic analysis of a taxonomic group.
    • Trimming with FASTP
    • Assembly with TRINITY
    • Cross contamination filtering
    • Orthology assignment with ORTHOGRAPH
    • Alignment with HMMALIGN
    • Filtered by multiple criteria (Aliscore, Outlier, Coverage)

## Prerequisites

In order to use this workflow, you will need singularity and anaconda installed on your system.
Detailed instructions can be found here:
singularity installation
conda installation
This workflow was built using Snakemake, which is a powerful workflow management system tool to create reproducible and scalable data analyses. For a detailed introduction please visit the snakemake documentation.

## Clone repository

In order to download the workflow to your current directory, please use the following git command:
git clone https://datacenter.zfmk.de/gitlab/smartin/lepi_te_workflow.git

## Setup

### Conda

To run the workflow we need to install snakemake and some other basic utilities. This can be done with the environment.yml file, which is part of this repository. It provides all the versions and packages that need to be installed. Create the conda environment with the following command:
conda env create -f environment.yml
After the environment was successfully created you need to activate it.
conda activate TEnriAn

### Singularity

The workflow uses a singularity container with a lot of software installed. This container brings its own operating system (Ubuntu 18.04), which makes you independent from your local operating system and ensures reproducibility of your analyses. The singularity container needs to be built from a recipe file. This can be done with the helper script.
You can build the container on your local system, but you will need sudo rights. A recipe file is located in workflow/container/EnGeTrAl_container.def. The singularity command to build the container is: sudo singularity build workflow/container/EnGeTrAl_container.{sif,def}
Building the container can take some time.

## Configuration and resources

### Config file

The singularity container is expected to be located in workflow/container/. If the container is located in a different directory, this can be specified in the workflow/config/config.yml file.
Parameters for the final alignment filter steps can be adjusted for each dataset.

### Raw data example

This workflow takes raw sequencing data from different samples as input. These fastq.gz files need to be located in resources/Raw_data.
As input we expect paired end data, where the file name consists of samplenameR1.fastq.gz and samplenameR2.fastq.gz.
The sample name needs to be specified in the workflow/config/species_samples.tsv file. Here a corresponding species name needs to be specified.
The sample name has to consist of three parts separated by an underscore, e.g. Family_Genus_spec-1. The species name is also used to name output files and folders. If you have multiple samples from the same species you need to make sure that you do not overwrite results, by choosing an appropriate identifier. During the contamination check this format is important. Here the species name is extracted from the first three parts of the sequence header, with the separators '_' AND '-'. Sequences from the same species are not evaluated for contamination. The third part of the species name of closely related subspecies or samples from the same species can then be further subdivided into parts that are separated by a '-'.

sample
species
9_S17
Epicopeiidae_Parabraxas_davidi
18_S33
Epicopeiidae_Psychostrophia_nymphidiaria-1
18_S34
Epicopeiidae_Psychostrophia_nymphidiaria-2
...
...







### External genomes or transcriptomes

External data from genomes or transcriptomes can be added to the workflow. They need to be located in workflow/config/species_external.tsv file.

species
sequence
Mecoptera_Panorpa_vulgaris
Mecoptera_Panorpa_vulgaris_transcriptome.fas
Noctuidae_Trichoplusia_ni
Noctuidae_Trichoplusia_ni_cds.fas
...
...

### Prepare Orthograph sets

In order to predict orthologous sequences from the assembled target enrichment sequences, Orthograph is used. For detailed information please visit the Orthograph GITHUB repository. To run this workflow you need to place the Orthograph sets and the sqlite database in the resources/orthograph/ directory and specify the name in the config file. If you start your analysis without a set and database, please have a look at the Orthograph_setup.sh script. You will need to make some adjustments in the beginning of the script, so that it works with your dataset. You will need the following input data:
    • Tab-delimited file that contains information about your orthologous groups, with the following three columns: 'name of ortho-group' ‘gene id' 'taxon_name' (Please have a look at the example file resources/orthograph/lepi-tabfile-exons.txt)
    • Official gene sets (OGS) that contain all predicted protein sequences for each taxon. Each OGS needs to be in a separate sequence file. The sequence headers need to be identical to the gene_id from the information file about the orthologous groups. Not all sequences in the OGS need to be included in the tab-delimited file, but the OGS needs to be complete in order to identify best reciprocal hits.
    • OGS sequences need to be in FASTA format and have a sequence header with only alphanumeric signs and '-' or '_' . Any special character can interfere with the software, which uses some characters for field delimitation.
    • OGS sequence files need to be placed in resources/orthograph/INPUT_OGS/ and the sequence file name needs to have the format 'TAXON_name.fas'

The Orthograph_setup.sh script will create an initial config file for Orthograph, then create the orthograph-db and upload all OGS files into this database. Then it will also upload the orthologous groups file. In this step some user interaction is needed. You have to input the name of the orthology set and confirm the OGS that you want to use. Next it will start an initial Orthograph run with resources/orthograph/Initialize_sets.fas. This may take some time because Orthograph will now create the sets. Feedback will be written to initial_orthograph_run.log file.

## Run the analysis

The workflow is divided into three parts:
    1. assembly (join raw reads, trim by quality, decompress reads, de-novo assemble with trinity, estimate abundance, contamination check)
    2. orthology (create config file for Orthograph, run orthograph-analyzer and orthograph-reporter, collect orthologous clusters)
    3. alignment_and_filter (hmmalign, convert Stockholm format, align corresponding nucleotide sequences, trim hmm-alignment, remove outliers, identify Target Enrichment taxa, filter alignment, Aliscore)

Each part must be executed one after the other:
snakemake --until run_assembly --use-singularity --use-conda --cores 12
snakemake --until run_orthology --use-singularity --use-conda --cores 12
snakemake --until run_alignment_filtering --use-singularity --use-conda --cores 12
Snakemake will create the environment automatically at the beginning of each run. It is possible to create them beforehand with the following command: snakemake --use-conda --use-singularity --conda-create-envs-only --cores 1. We need to use the option --use-singularity in order to activate the singularity container that provides all the major tools for this pipeline. Also for some steps, a specific conda environment needs to be activated so we have to specify the --use-conda option. Each rule will check that the previous rule was run and look for flag-files in the results folder. If you want to rerun a specific rule you can also use the --forcerun option. Before you start your actual run it is recommended to use the --dry-run option, to check if everything is set up correctly. The results of each step will be written to the corresponding directory in the results folder. For more detailed information on how to run snakemake, please see the snakemake documentation.

## Citation

OUR PAPER

## References

    • Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33. https://doi.org/10.12688/f1000research.29032.1
    

singularity Orthograph fastp hmmalign trinity
