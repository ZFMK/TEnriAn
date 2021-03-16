####################################
## rules for orthology prediction ##
####################################

import pandas as pd
configfile: "workflow/config/config.yaml"
container: config["container"]["singularity"]

samples = pd.read_table(config["sample_info"] , dtype=str).set_index(["sample"], drop=False)
external = pd.read_table(config["external_evidence"] , dtype=str).set_index(["species"], drop=False)

## functions

def get_ext_species(wildcards):
    #transform dataframe to dictionary
    area_dict = dict(zip(external['sequence'],external['species']))
    species = area_dict[wildcards.file]
    print('the result for the species extraction {}'.format(species))
    return species

def get_ext_file(wildcards):
    #transform dataframe to dictionary
    area_dict = dict(zip(external['species'],external['sequence']))
    file = area_dict[wildcards.species]
    print('the result for the file extraction {}'.format(file))
    return 'resources/External_data/' + file

def return_list_of_ext_orthoIn():
    mylist=[]
    for dataset in external['species'].tolist():
        mylist.append("results/02_orthology/input/{dataset}.fas".format(dataset=dataset))
    return set(mylist)

def return_list_of_orthoOut():
    mylist=[]
    for dataset in external['species'].tolist():
        mylist.append("results/02_orthology/results/{species}/summary.txt".format(species=dataset))
    for dataset in samples['species'].tolist():
        mylist.append("results/02_orthology/results/{species}/summary.txt".format(species=dataset))
        
    return set(mylist)

rule add_Tax2headers_external:
    input:
        IN=get_ext_file
    output:
        OUT="results/02_orthology/input/{species}.fas" 
    shell:
        "cp {input.IN} {output.OUT} ; "
        "sed -i 's/ .*//g' {output.OUT} ; "
        "if grep -q '>{wildcards.species}_' {output.OUT} ; then echo {output.OUT} already contains species name: {wildcards.species} ; "
        " else sed -i 's/>/>{wildcards.species}_/g' {output.OUT} ;fi"

rule check_files_complete_external:
    input: 
        return_list_of_ext_orthoIn()
    output: 
        "results/02_orthology/prep_ext_input.done"
    shell:
        "touch {output}"

rule orthograph_config:
    input:
        check="results/01_assembly-step.complete",
        check2="results/02_orthology/prep_ext_input.done"
    output:
        conf="results/02_orthology/config/{species}_orthograph.conf"
    params:
        assembly="results/02_orthology/input/{species}.fas" ,
        out_dir="results/02_orthology/results/{species}",
        sets_dir=config["orthograph"]["sets_dir"],
        ortho_set=config["orthograph"]["ortho_set"],
        database=config["orthograph"]["database"]
    threads:
        13
    shell:
        "printf \"species-name = {wildcards.species} \ninput-file = {params.assembly}\noutput-directory = {params.out_dir} \n\" > {output.conf} ;"
        "printf \"database-backend = sqlite \nsqlite-program = sqlite3\nsqlite-database = {params.database}\n\" >> {output.conf} ;"
        "printf \"sets-dir = {params.sets_dir}\n\" >> {output.conf} ;"
        "printf \"ortholog-set = {params.ortho_set}\nalignment-program = mafft-linsi \nhmmbuild-program = hmmbuild\nmakeblastdb-program  = makeblastdb\n\" >> {output.conf};"
        "printf \"translate-program = fastatranslate\nhmmsearch-program = hmmsearch\nblast-program = blastp\nexonerate-program    = exonerate\n\" >> {output.conf};"
        "printf \"num-threads = {threads}\nbrh-only = 1\n\" >> {output.conf}"


rule orthograph_analyzer:
    input:
        conf="results/02_orthology/config/{species}_orthograph.conf",
    output:
        #directory("results/02_orthology/results/{species}/"),
        "results/02_orthology/results/{species}/summary.txt"
    log:
        analyzer="results/logs/02_orthology/{species}_orthograph_run.log",
        reporter="results/logs/02_orthology/{species}_orthograph_run_report.log"
    threads:
        13
    shell:
        "orthograph-analyzer -c {input.conf} &> {log.analyzer} ; "
        "orthograph-reporter -c {input.conf} &> {log.reporter} ; "

rule check_orthograph_results:
    input: 
        return_list_of_orthoOut()
    output: 
        "results/02_orthology/orthograph.done"
    shell:
        "touch {output}"


rule Cluster_collection:
    input:
        "results/02_orthology/orthograph.done"
    output:
        dir=directory("results/02_orthology/Cluster_collection"),
        check="results/02_orthology/orthograph_collection.done"
    conda:
        "../envs/bio.yaml"
    shell:
        "mkdir {output.dir} ; "
        "python workflow/scripts/Cluster_collection.py ; touch {output.check} "


rule orthograph_cleanup:
    input:
        "results/02_orthology/orthograph_collection.done"
    output:
        check="results/02_orthology/orthograph_cleanup.done"
    params:
        ortho_dir="results/02_orthology/results"
    conda:
        "../envs/bio.yaml"
    shell:
        "rm -r {params.ortho_dir} ; touch {output.check} "



