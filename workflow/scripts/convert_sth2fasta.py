#!/usr/bin/env python3

from Bio import SeqIO
import os
import re
import sys

seqfile = sys.argv[1]
output = sys.argv[2]

path = os.path.dirname(seqfile)

input_handle = open(seqfile, "r")


#output = os.path.splitext(seqfile)[0] +'.fas'

#print(seqfile)
#if os.path.exists(output):
#    os.remove(output)

output_handle = open(output, "w")

alignments = SeqIO.parse(input_handle, "stockholm")

SeqIO.write(alignments, output_handle, "fasta")

output_handle.close()
input_handle.close()
 
