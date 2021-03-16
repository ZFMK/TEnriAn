#!/usr/bin/env python3
#python script to identify stop-codons in aminoacid sequences and replace them with X. Also replace corresponding nucleotides...

import sys
import re

from Bio import SeqIO
from Bio.Seq import Seq

input_AA   =   sys.argv[1]
input_NUC   =   sys.argv[2]

f_aa= open(re.sub('.fas', '.clean.fas', input_AA),"w")
f_nuc= open(re.sub('.fas', '.clean.fas', input_NUC),"w")

#load files in dictionary
SequenceFileDict_AA=SeqIO.to_dict(SeqIO.parse( input_AA, "fasta"))
SequenceFileDict_NUC=SeqIO.to_dict(SeqIO.parse( input_NUC, "fasta"))

for title in SequenceFileDict_AA.keys():

    if '*' in str(SequenceFileDict_AA[title].seq) :
        
        print('we have a stop codon in {}'.format(title))
        
        AAs     = list(SequenceFileDict_AA[title].seq)
        NUCs    = list(SequenceFileDict_NUC[title].seq)
        
        for pos, item in enumerate(AAs): #not good python style...
            if item == '*': #find position in list that has a stop codon

                AAs[pos]='X'
                
                for j in range((pos*3),(pos*3)+3): #also replace in nucleotide sequences...

                    NUCs[j]='N'
                
        f_aa.write('>{}\n{}\n'.format( title,''.join(AAs) ))
        f_nuc.write('>{}\n{}\n'.format( title,''.join(NUCs) ))
        
        
    else:
        f_aa.write('>{}\n{}\n'.format( title,str(SequenceFileDict_AA[title].seq) ))
        f_nuc.write('>{}\n{}\n'.format( title, str(SequenceFileDict_NUC[title].seq) ))
