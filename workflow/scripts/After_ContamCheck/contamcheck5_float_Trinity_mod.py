#!/usr/bin/etc python3

# prerequisites:
# samtools, python 3.5+, result from ncbi blast
# add identifier for each sample as first part of fasta header, e.g. >species1_contig1...
# ATTENTION: read numbers should be in the 8th part of the fasta header (seperated by underscore)
# combine all samples in one multifasta-file
# before using the skript perform a blast all-vs-all with options:
# makeblastdb -in MYSEQUENCES.fasta -out MY_SEQ_DB -dbtype nucl
# blastn -query MYSEQUENCES.fasta -db MY_SEQ_DB -evalue 1e-50 -perc_identity 99 -outfmt 6 -out blast_all_vs_all

import subprocess
import sys
from os import path

if len(sys.argv) == 1:
    print("please provide filename with all sequences")
    quit()

Sequences = sys.argv[1] # contains all sequences of all taxa

Out_dir = sys.argv[2] #path to the output directory

if path.exists(Sequences):
    print("file exists")
else:
    exit("file not found")

Blastresult=Out_dir + "blast_all_vs_all" # name of the output file from blast all versus all (-evalue 1e-50, -outfmt 6)

if path.exists(Out_dir + "blast_all_vs_all"):
    print("blastresult provided, start filtering")
else:
    subprocess.run("makeblastdb -in "+Sequences+" -out "+ Out_dir +"mydb -dbtype nucl",shell=True)
    subprocess.run("blastn -query "+Sequences+" -db "+ Out_dir +"mydb -evalue 1e-50 -out "+Blastresult+" -outfmt 6 -perc_identity 99 -num_threads 8",shell=True)

Column = 9   # the readnumber is x-th element of the fasta header, when seperated by underscores (~ abundance)
Column = Column-1 # e.g. element number 8 is number 7 in python list

# Filters hits from the same species

BLAST = open(Blastresult,"r")
FILTERED = open(Out_dir +"filtered_blast_result.txt","w")

for Line in BLAST:
    Items = Line.split("\t")
    Query = Items[0]
    Subj = Items[1]
    Quid = Query.split("_")
    Sbid = Subj.split("_")
    Qspec = Quid[0]+"_"+Quid[1]+"_"+Quid[2]
    Sspec = Sbid[0]+"_"+Sbid[1]+"_"+Sbid[2]
    if Qspec != Sspec:
        FILTERED.write(Line)

BLAST.close()
FILTERED.close()


# Makes black list with potential contaminations to be removed
# As well makes a list of most frequent species pairs

FILTERED = open(Out_dir +"filtered_blast_result.txt","r")
BLACKLIST = open(Out_dir +"blacklist.txt","w")
PAIRS = open(Out_dir +"blast_pairs.txt","w")

for Line in FILTERED:
    Items = Line.split("\t")
    Query = Items[0]
    Subj = Items[1]
    Quid = Query.split("_")
    Sbid = Subj.split("_")
    if float(Sbid[Column]) == 0 or float(Quid[Column])/float(Sbid[Column]) > 2: 
        BLACKLIST.write(Subj)
        BLACKLIST.write("\n")
    elif float(Quid[Column])/float(Sbid[Column]) < 0.5:
        BLACKLIST.write(Query)
        BLACKLIST.write("\n")
    else:
        BLACKLIST.write(Subj)
        BLACKLIST.write("\n")
        BLACKLIST.write(Query)
        BLACKLIST.write("\n")

    PAIRS.write("_".join([Quid[0],Quid[1],Quid[2],Sbid[0],Sbid[1],Sbid[2],"\n"]))

FILTERED.close()
BLACKLIST.close()
PAIRS.close()

subprocess.run("sort "+ Out_dir +"blacklist.txt | uniq > "+ Out_dir +"blacklist_uniq.txt", shell=True)
subprocess.run("cut -f1 -d'_' "+ Out_dir +"blacklist_uniq.txt | uniq -c > "+ Out_dir +"how_many_contam_per_species.txt", shell=True)
subprocess.run("sort "+ Out_dir +"blast_pairs.txt | uniq -c | sort -nr > "+ Out_dir +"species_pairs_by_occurrence.txt", shell=True)
subprocess.run("rm "+ Out_dir +"blacklist.txt "+ Out_dir +"blast_pairs.txt", shell=True)

# makes the white list (contamination free sequences)

subprocess.run("grep '>' "+Sequences+" | sed 's/>//' > "+ Out_dir +"all_headers.txt", shell=True)
subprocess.run("cat "+ Out_dir +"blacklist_uniq.txt "+ Out_dir +"all_headers.txt | sort | uniq -u > "+ Out_dir +"whitelist.txt", shell=True)
subprocess.run("rm "+ Out_dir +"all_headers.txt",shell=True)

#subprocess.run("samtools faidx "+Sequences, shell=True)
#subprocess.run("xargs -a blacklist_uniq.txt -i samtools faidx "+Sequences+" {} > seq_contaminated.fasta", shell=True)
#subprocess.run("xargs -a whitelist.txt -i samtools faidx "+Sequences+" {} > seq_contam_free.fasta", shell=True)

subprocess.run("python workflow/scripts/After_ContamCheck/seqout.py "+ Out_dir +"blacklist_uniq.txt "+Sequences+" "+ Out_dir +"seq_contaminated.fasta", shell=True)
subprocess.run("python workflow/scripts/After_ContamCheck/seqout.py "+ Out_dir +"whitelist.txt "+Sequences+" "+ Out_dir +"seq_contam_free.fasta", shell=True)

#sorts sequences into single species files

subprocess.run("mkdir "+ Out_dir +"free_by_species", shell=True)
Oldname = "   "
flag = 1
FILE = open(Out_dir +"seq_contam_free.fasta","r")
for Line in FILE:
    Line = Line.strip()
    if ">" in Line:
        Items = Line.split("_")
        Name = Items[0][1:]+"_"+Items[1]+"_"+Items[2] # +"_"+Items[3]
        if Name == Oldname:
            NEWSPEC.write(Line+"\n")
        else:
            if flag == 1:
                NEWSPEC=open(Out_dir +"free_by_species/"+Name+".fas","w")
                NEWSPEC.write(Line+"\n")
                Oldname=Name
                flag = 2
            else:
                NEWSPEC.close()
                NEWSPEC=open(Out_dir +"free_by_species/"+Name+".fas","w")
                NEWSPEC.write(Line+"\n")
                Oldname=Name
    else:
        NEWSPEC.write(Line+"\n")
NEWSPEC.close()
FILE.close()


# sorts also the contaminated ones

subprocess.run("mkdir " + Out_dir +"contam_by_species", shell=True)

Oldname = "   "
flag = 1
FILE = open(Out_dir +"seq_contaminated.fasta","r")
for Line in FILE:
    Line=Line.rstrip()
    if ">" in Line:
        Items = Line.split("_")
        Name = Items[0][1:]+"_"+Items[1]+"_"+Items[2] #+"_"+Items[3]
        if Name == Oldname:
            NEWSPEC.write(Line+"\n")
        else:
            if flag == 1:
                NEWSPEC=open(Out_dir +"contam_by_species/"+Name+".fas","w")
                NEWSPEC.write(Line+"\n")
                Oldname=Name
                flag = 2
            else:
                NEWSPEC.close()
                NEWSPEC=open(Out_dir +"contam_by_species/"+Name+".fas","w")
                NEWSPEC.write(Line+"\n")
                Oldname=Name
    else:
        NEWSPEC.write(Line+"\n")
NEWSPEC.close()
FILE.close()


###subprocess.run("rm seq_contam_free.fasta seq_contaminated.fasta",shell=True)

