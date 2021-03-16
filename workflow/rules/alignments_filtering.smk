########################################
## rules for alignment and filtering: ##
########################################

#configfile: "workflow/config/config.yaml"
#container: config["container"]["singularity"]

samples = pd.read_table(config["sample_info"] , dtype=str).set_index(["sample"], drop=False)
IDS, = glob_wildcards("results/02_orthology/Cluster_collection/{id}.aa.fas")

## functions
def get_list_of_final():
    mylist=[]
    for ID in IDS :
        mylist.append("results/03_alignments/final_aa/{EOG}.log".format(EOG=ID))
    return set(mylist)


rule Remove_stopCodons:
    input:
        aa="results/02_orthology/Cluster_collection/{EOG}.aa.fas",
        nt="results/02_orthology/Cluster_collection/{EOG}.nt.fas",
    output:
        aa="results/02_orthology/Cluster_collection/{EOG}.aa.clean.fas",
        nt="results/02_orthology/Cluster_collection/{EOG}.nt.clean.fas"
    conda:
        "../envs/bio.yaml"
    log:
        "results/logs/03_alignments/cleaned_collection/{EOG}.log"
    shell:
        "python ./workflow/scripts/Remove_stopCodons_seqPairs.py {input.aa} {input.nt}"


rule hmmalign:
    input:
        aa="results/02_orthology/Cluster_collection/{EOG}.aa.clean.fas",
        #hmm="resources/orthograph/sets/Lepidoptera_orthoDB9_extended/hmms/{EOG}.hmm"
    output:
        "results/03_alignments/hmmaligned/{EOG}.ali.sth"
    params:
        sets_dir=config["orthograph"]["sets_dir"],
        ortho_set=config["orthograph"]["ortho_set"]
    log:
        "results/logs/03_alignments/hmmalign/{EOG}.log"
    shell:
        "( hmmalign -o {output} --amino  --informat FASTA --outformat Stockholm"
        " {params.sets_dir}/{params.ortho_set}/hmms/{wildcards.EOG}.hmm {input.aa} ) &> {log}"


rule trans_stockholm:
    input: 
        "results/03_alignments/hmmaligned/{EOG}.ali.sth"
    output: 
        "results/03_alignments/hmmaligned/{EOG}.ali.fas"
    conda:
        "../envs/bio.yaml"
    log:
        "results/logs/03_alignments/trans_stockholm/{EOG}.log"
    shell:
        "./workflow/scripts/convert_sth2fasta.py {input} {output} &> {log}"


rule nt_ali:
    input: 
        aa_ali="results/03_alignments/hmmaligned/{EOG}.ali.fas",
        nt="results/02_orthology/Cluster_collection/{EOG}.nt.clean.fas"
    output: 
        "results/03_alignments/hmmaligned_nt/{EOG}.nt_ali.fas"
    log:
        "results/logs/03_alignments/nt_ali/{EOG}.log"
    shell:
        "( perl ./workflow/scripts/pal2nal.mod.pl -output fasta {input.aa_ali} {input.nt} > {output} ) &> {log}"


rule trim_hmmali:
    input:
        "results/03_alignments/hmmaligned/{EOG}.ali.sth",
        "results/03_alignments/hmmaligned_nt/{EOG}.nt_ali.fas"
    output:
        "results/03_alignments/trimmed_hmmali/{EOG}.aa.fas",
        "results/03_alignments/trimmed_hmmali/{EOG}.nt.fas"
    log:
        "results/logs/03_alignments/trim_hmmali/{EOG}.log"
    shell:
        "perl ./workflow/scripts/hmmalign_cut_mod.pl {input} {output} &> {log}"


rule remove_outlier:
    input:
        aa="results/03_alignments/trimmed_hmmali/{EOG}.aa.fas",
        nt="results/03_alignments/trimmed_hmmali/{EOG}.nt.fas"
    output:
        aa_out="results/03_alignments/trimmed_hmmali/{EOG}.aa_orm.fas",
        nt_out="results/03_alignments/trimmed_hmmali/{EOG}.nt_orm.fas"
        #this nt_out is produces by default from the program, but needs to be specified in order to be picked up later
    log:
        "results/logs/03_alignments/remove_outlier/{EOG}.log"
    shell:
        "./workflow/scripts/remove_outlier_sequences-v0.9.5-dist/Remove-Outlier-Sequences-v0.9.5 --remove-gap-ambig-lowerCase-sites --windowsSize 15 --IQR-factor 2.0 "
        "-i {input.aa} -o {output.aa_out} --corresponding-nuc-fasta-file-name {input.nt} &> {log}"
#sometimes it can happen, that no outlier sequences are found and no output is generated... 
#double check that output is present and otherwise
#copy input to output...


rule get_HE_taxa_list:
    output:
        HEtax="workflow/config/HE-taxa_list.txt"
    run:
        f = open( output.HEtax, 'w')
        for dataset in samples['species'].tolist():
            f.write("{dataset}\n".format(dataset=dataset))


rule filter_1:
    input:
        aa="results/03_alignments/trimmed_hmmali/{EOG}.aa_orm.fas",
        nt="results/03_alignments/trimmed_hmmali/{EOG}.nt_orm.fas",
        HEtax="workflow/config/HE-taxa_list.txt"
    output:
        aa="results/03_alignments/trimmed_hmmali/{EOG}.aa_orm_filtered.fas",
        nt="results/03_alignments/trimmed_hmmali/{EOG}.nt_orm_filtered.fas"
        #output files are named after input file with '_filtered' in the end
    params:
        HE=config["filtering"]["HE_per_position"],
        taxa=config["filtering"]["taxa_per_position"],
        min_aa=config["filtering"]["min_aa"]
    conda:
        "../envs/bio.yaml"
    log:
        "results/logs/03_alignments/filter_1/{EOG}.log"
    shell:
        "(python ./workflow/scripts/filter_high-zones_options.py --aa_input {input.aa} --nt_input {input.nt} "
        "--HE_list {input.HEtax} --simple_header --HE_per_position {params.HE} --taxa_per_position {params.taxa} --min_aa {params.min_aa} "
        "; if [ -s {output.aa} ]; then skipredundant -sequences {output.aa} -mode 1 -threshold 99 -minthreshold 10 -maxthreshold 99 -gapopen 10 -gapextend 0.5 -outseq {output.aa}.check -redundantoutseq {output.aa}.red ;"
        " if [ $(grep -c '>' {output.aa}.check ) -eq 1 ] ; then rm {output.aa} {output.nt} ; fi ; fi "
        "; if [ ! -f {output.aa} ]; then touch {output.aa}; fi"
        "; if [ ! -f {output.nt} ]; then touch {output.nt}; fi ) &> {log}"

rule pairwise_sequence_distance:
    input:
        aa="results/03_alignments/trimmed_hmmali/{EOG}.aa_orm_filtered.fas",
        nt="results/03_alignments/trimmed_hmmali/{EOG}.nt_orm_filtered.fas"
    output:
        aa="results/03_alignments/trimmed_hmmali/{EOG}.aa_orm_filtered_pairwise.fas",
        nt="results/03_alignments/trimmed_hmmali/{EOG}.nt_orm_filtered_pairwise.fas"
    params:
        aa="results/03_alignments/trimmed_hmmali/{EOG}.aa_orm_filtered_pairwise.fas.o"
    log:
         "results/logs/03_alignments/pairwise_seq_dist/{EOG}.log"
    shell:
        "( if [ -s {input.aa} ] ; then ./workflow/scripts/Pairwise-sequence-distance-using-existing-alignment-v1.3-dist_tenrian/Pairwise-alignment-distances-identity-v1.3 --fasta-input-file-name {input.aa} --data-type PROTEIN --minimum-pairs-in-window 10 &> {params.aa} ; "
        " Num_below=$(grep 'Number defined below 0.50:' {params.aa} | cut -f2 -d':' | tr -d ' ') ; "
        " if [ ! -z $Num_below ] ; then "
        " if [ $Num_below -lt 10 ] ; then cp {input.aa} {output.aa} ; cp {input.nt} {output.nt} ; echo 'Number of Sequences below 0.5 : $Num_below' ; "
        " else touch {output.aa} ; touch {output.nt} ; echo 'Number of Sequences below 0.5 : $Num_below' ; fi ; fi " 
        "   else echo {input.aa} is empty; touch {output.aa}; touch {output.nt} ; fi)&> {log}"

rule Aliscore:
    input:
        aa="results/03_alignments/trimmed_hmmali/{EOG}.aa_orm_filtered_pairwise.fas",
        nt="results/03_alignments/trimmed_hmmali/{EOG}.nt_orm_filtered_pairwise.fas",
        HEtax="workflow/config/HE-taxa_list.txt"
    output:
        aa="results/03_alignments/trimmed_hmmali/{EOG}.aa_orm_filtered_pairwise_filtered.fas",
        nt="results/03_alignments/trimmed_hmmali/{EOG}.nt_orm_filtered_pairwise_filtered.fas"
    params:
        HE=config["filtering"]["HE_per_position"],
        taxa=config["filtering"]["taxa_per_position"],
        min_aa=config["filtering"]["min_aa"],
        min_ali=config["filtering"]["min_ali_length"]
    conda:
        "../envs/bio.yaml"
    log:
        "results/logs/03_alignments/aliscore/{EOG}.log"
    shell:
        "( if [ -s {input.aa} ]; then perl ./workflow/scripts/ALISCORE_v2.0/Aliscore.02.2.pl -r 1000000 -i {input.aa}  ;"
        "rm {input.aa}.profiles.svg ;"
        "python ./workflow/scripts/filter_high-zones_options.py --aa_input {input.aa} --nt_input {input.nt} "
        "--filter_aliscore {input.aa}_List_random.txt --min_ali_length {params.min_ali} --HE_list {input.HEtax} --HE_per_position {params.HE} --taxa_per_position {params.taxa} --min_aa {params.min_aa} ; "
        " if [ ! -f {output.aa} ]; then touch {output.aa}; fi ;"
        " if [ ! -f {output.nt} ]; then touch {output.nt}; fi; "
        "else echo {input.aa} is empty; touch {output.aa}; touch {output.nt} ; fi)&> {log}"


rule final_collection:
    input:
        aa="results/03_alignments/trimmed_hmmali/{EOG}.aa_orm_filtered_pairwise_filtered.fas",
        nt="results/03_alignments/trimmed_hmmali/{EOG}.nt_orm_filtered_pairwise_filtered.fas"
    output:
        aa=temp("results/03_alignments/final_aa/{EOG}.log"), 
        nt=temp("results/03_alignments/final_nt/{EOG}.log"),
    params:
        aa="results/03_alignments/final_aa/{EOG}.fas",
        nt="results/03_alignments/final_nt/{EOG}.fas"
    log:
        "results/logs/03_alignments/final/{EOG}.log"
    shell:
        "( if [ -s {input.aa} ]; then cp {input.aa} {params.aa} ; cp {input.nt} {params.nt}; "
        " echo 'filtering was successful' ; "
        "else echo 'nothing left after filtering to use for final analysis' ; fi ; touch {output.aa} ; touch {output.nt} ) &> {log}"


rule check_filter_complete:
    input: 
        get_list_of_final()
    output: 
        "results/03_alignments/03_alignments_filtering.complete"
    shell:
        "touch {output}"

rule abs_pres_matrix:
    input:
        check="results/03_alignments/03_alignments_filtering.complete"
    output:
        check="results/03_alignments/03_alignments_filtering-step.complete",
        file="results/03_alignments/abs-pres-matrix.txt"
    params:
        dir="results/03_alignments/final_aa",
    log:
        "results/logs/03_alignments/abs-pres-matrix.log"
    shell:
        "(./workflow/scripts/final_overview.sh -i {params.dir} -o {output.file} ; touch {output.check}) &> {log} " 
