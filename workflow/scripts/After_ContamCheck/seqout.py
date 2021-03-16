#!/usr/bin/etc python3

from Bio import SeqIO
import sys

Reallist=sys.argv[1]
Alltaxa=sys.argv[2]
Outfile=sys.argv[3]

LIST = open(Reallist,"r") ### list of hmm files
ALLSEQ = open(Alltaxa,"r") ### fasta file with all sequences from all taxa
Allfasta_dict = SeqIO.to_dict(SeqIO.parse(ALLSEQ, 'fasta'))
ALLSEQ.close()

MYOUTFILE = open(Outfile,"w")
Wanted = [line.rstrip() for line in LIST]
for Header in Wanted:
    SeqIO.write(Allfasta_dict[Header],MYOUTFILE,'fasta')
MYOUTFILE.close()

