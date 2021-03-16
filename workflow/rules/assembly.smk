###################################
## rules for assembly of TE data ##
###################################

import pandas as pd

configfile: "workflow/config/config.yaml"
# container: config["container"]["singularity"]

samples = pd.read_table(config["sample_info"] , dtype=str).set_index(["sample"], drop=False)
#samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])  # enforce str in index

## functions:

def get_sample_species(sample_id):
    #transform dataframe to dictionary
    area_dict = dict(zip(samples['sample'],samples['species']))
    species = area_dict[sample_id]
    print('the result for the species extraction {}'.format(species))
    return species

def return_list_of_sample_output():
    mylist=[]
    for dataset in samples['sample'].tolist():
        mylist.append("results/01_assembly/abundance/{dataset}/{dataset}.renamed.fas".format(dataset=dataset))
    return set(mylist)


rule join_raw_reads:
    output:
        R1=temp("results/01_assembly/Raw_merged/{sample}_R1.fq.gz"),
        R2=temp("results/01_assembly/Raw_merged/{sample}_R2.fq.gz"), 
    shell:
        "cat resources/Raw_data/{wildcards.sample}*R1_001.fastq.gz >> {output.R1} ;"
        "cat resources/Raw_data/{wildcards.sample}*R2_001.fastq.gz >> {output.R2} ;"


rule Trim_raw2:
    input:
        R1="results/01_assembly/Raw_merged/{sample}_R1.fq.gz",
        R2="results/01_assembly/Raw_merged/{sample}_R2.fq.gz",
    output:
        R1=temp("results/01_assembly/Raw_trimmed/{sample}_R1.fq.gz"),
        R2=temp("results/01_assembly/Raw_trimmed/{sample}_R2.fq.gz")
    params:
        spec=lambda wildcards: get_sample_species('{sample_id}'.format(sample_id=wildcards.sample))
    conda:
        "../envs/trim.yaml"
    log:
        "results/logs/01_assembly/trim_raw/{sample}.log"
    shell:
        "(fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} "
        " --length_required 100 --low_complexity_filter --detect_adapter_for_pe "
        " --json results/logs/01_assembly/trim_raw/{wildcards.sample}.fastp.json"
        " --html results/logs/01_assembly/trim_raw/{wildcards.sample}.fastp.html  ; echo '{params.spec}') &> {log}"


rule decompress_reads:
    input:
        R1="results/01_assembly/Raw_trimmed/{sample}_R1.fq.gz",
        R2="results/01_assembly/Raw_trimmed/{sample}_R2.fq.gz"
    output:
        R1="results/01_assembly/Raw_trimmed/{sample}_R1.fastq",
        R2="results/01_assembly/Raw_trimmed/{sample}_R2.fastq"
    shell:
        "zcat {input.R1} > {output.R1} ; zcat {input.R2} > {output.R2}"


rule trinity_assembly:
    input:
         R1="results/01_assembly/Raw_trimmed/{sample}_R1.fastq",
         R2="results/01_assembly/Raw_trimmed/{sample}_R2.fastq"
    output:
        out="results/01_assembly/assembly/{sample}_trinity.Trinity.fas"
    params:
        outdir="results/01_assembly/assembly/{sample}_trinity",
        max_ram=config["assembly"]["Trinity_Max_RAM"]
    threads: 55
    conda:
        "../envs/trinity.yaml"
    log:
        "results/logs/01_assembly/assembly/{sample}_trinity.log"
    shell:
        "(printf 'starting Trinity assembly ...\n\n' ; Trinity --seqType fq --left {input.R1} --right {input.R2} "
        "--CPU {threads} --SS_lib_type FR --no_version_check --max_memory {params.max_ram} --output {params.outdir} --full_cleanup ; "
        "sleep 100 ; cp {params.outdir}.Trinity.fasta {output.out} ) &> {log} " 


rule add_Tax2headers:
    input:
        IN="results/01_assembly/assembly/{sample}_trinity.Trinity.fas"
    output:
        OUT="results/01_assembly/assembly_renamed/{sample}.fas"
    params:
        spec=lambda wildcards: get_sample_species('{sample_id}'.format(sample_id=wildcards.sample))
    shell:
        "cp {input.IN} {output.OUT} ; "
        "sed -i 's/ .*//g' {output.OUT} ; "
        "sed -i 's/>/>{params.spec}_/g' {output.OUT}"


rule estimate_abundance:
    input:
        R1="results/01_assembly/Raw_trimmed/{sample}_R1.fastq",
        R2="results/01_assembly/Raw_trimmed/{sample}_R2.fastq",
        TE="results/01_assembly/assembly_renamed/{sample}.fas"
    output:
        rsem="results/01_assembly/abundance/{sample}/RSEM.isoforms.results",
        fas="results/01_assembly/abundance/{sample}/{sample}.fas"
    params:
        dir="results/01_assembly/abundance/{sample}"
    threads: 
        11
    container: 
        config["container"]["singularity"]
    log:
        "results/logs/01_assembly/abundance/{sample}_rsem.log"
    shell:
        "cp {input.TE} {output.fas} ; "
        "align_and_estimate_abundance.pl --transcripts {output.fas}"
        " --seqType fq --left {input.R1} --right {input.R2} --est_method RSEM --aln_method bowtie2 --trinity_mode "
        " --prep_reference --thread_count {threads} --output_dir {params.dir} --coordsort_bam &> {log} "


rule add_abundance_into_to_header:
    input:
        Abundance="results/01_assembly/abundance/{sample}/RSEM.isoforms.results",
        FASTA="results/01_assembly/abundance/{sample}/{sample}.fas"
    output:
        Renamed="results/01_assembly/abundance/{sample}/{sample}.renamed.fas"
    params:
        RSEM_headers="results/01_assembly/abundance/{sample}/{sample}.RSEM_headers",
        assmbl_headers="results/01_assembly/abundance/{sample}/{sample}.assmbl_headers",
        fasta_sorted="results/01_assembly/abundance/{sample}/{sample}.sorted.fas",
        new_headers="results/01_assembly/abundance/{sample}/{sample}.new_headers",
        header_pairs="results/01_assembly/abundance/{sample}/{sample}.header_pairs"
    conda:
        "../envs/seqkit.yaml"
    shell:
        "cat {input.Abundance} | cut -f 1 | tr -s '\t' '_' | tail -n +2 | sort  > {params.RSEM_headers} ; "
        "seqkit sort -2 -n {input.FASTA} > {params.fasta_sorted} ; "
        "grep '>' {params.fasta_sorted} | sed 's/>//' > {params.assmbl_headers} ; "
        
	    "if [ \"$(diff {params.RSEM_headers} {params.assmbl_headers})\" = '' ] ; then "
        "cat {input.Abundance} | cut -f 1,4 | tr -s '\t' '_' | tail -n +2 | sort > {params.new_headers} ; "
        " paste {params.assmbl_headers} {params.new_headers} > {params.header_pairs} ; "
        " seqkit replace -p '(.+)' --replacement '{{kv}}' --kv-file {params.header_pairs} {params.fasta_sorted} > {output.Renamed}; fi"


rule check_files_complete:
    input: 
        return_list_of_sample_output()
    output: 
        "results/01_assembly/contam_check/check_files.done"
    #message: 
    #    "Confirming that required files exist. Contamcheck can now proceed..."
    shell:
        "touch {output}"


rule join_assemblies_for_contamcheck:
    input:
        "results/01_assembly/contam_check/check_files.done"
    output:
        "results/01_assembly/contam_check/all_assemlies.fas"
    params:
        assembly_list=' '.join(return_list_of_sample_output())
    shell:
        "cat {params.assembly_list} > {output}"


rule run_contamcheck:
    input:
        "results/01_assembly/contam_check/all_assemlies.fas"
    output:
        "results/01_assembly/contam_check/run_complete.check"
    params:
        out_dir="results/01_assembly/contam_check/"
    conda:
        "../envs/contam.yaml"
    log:
        "results/logs/01_assembly/contam_check/run_contamcheck.log"
    shell:
        "python workflow/scripts/After_ContamCheck/contamcheck5_float_Trinity_mod.py {input} {params.out_dir} &> {log} ; "
        "touch {output} &>> {log} "


rule prepare_sequences_for_orthograph:
    input:
        "results/01_assembly/contam_check/run_complete.check"
    output:
        check="results/02_orthology/data_assembly_complete.check",
    params:
        out_dir=directory("results/02_orthology/input/")
    shell:
        "[ ! -d {params.out_dir} ] && mkdir -p {params.out_dir} ; cp results/01_assembly/contam_check/free_by_species/*.fas {params.out_dir}  ; touch {output.check}"
