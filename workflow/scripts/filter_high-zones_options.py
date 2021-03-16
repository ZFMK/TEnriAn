#!/usr/bin/env python3
import pdb

import sys
import time
import re
import os
import subprocess
import argparse

from Bio import SeqIO
from Bio.Seq import Seq


def simplify_header(Seq_Dict):
    #SequenceFileDict=SeqIO.to_dict(SeqIO.parse(seqfile, "fasta")) #load sequence in dictionary....
    mod_seq_header={}

    for key in Seq_Dict.keys():
        parts=key.split("_")

        new_header = '_'.join( [parts[0],parts[1],parts[2]] )    
        mod_seq_header[new_header] = Seq_Dict[key].seq 
       
    return(mod_seq_header)

def reformat(Seq_Dict):
    #SequenceFileDict=SeqIO.to_dict(SeqIO.parse(seqfile, "fasta")) #load sequence in dictionary....
    mod_seq_header={}

    for key in Seq_Dict.keys():
            
        mod_seq_header[key] = Seq_Dict[key].seq 
       
    return(mod_seq_header)

def get_HE_seqs(headerlist,Seq_Dict):
    HE_seqs = {}
    for key in Seq_Dict.keys():
        if (key in headerlist):
            HE_seqs[key]= Seq_Dict[key]
    if(bool(HE_seqs)):
        return(HE_seqs)
    else:
        print('the file does not contain any of the Target Enrichment species')
        sys.exit(0)


def count_he_per_position(he_seqs):
#
    total_length = len(next(iter(he_seqs.values())))
    PositionCount = {}

    #count HEs for each position
    for i in range(0, total_length, 1):
        PositionCount[i]=0
        for title in he_seqs.keys():
            if (he_seqs[title][i] != '-' and he_seqs[title][i] != 'X') :
                PositionCount[i] += 1

    return(PositionCount)

def count_taxa_per_position(sequences):
#
    total_length = len(next(iter(sequences.values())))  #first entry of the dictionary....
    PositionCount = {}

    #count HEs for each position
    for i in range(0, total_length, 1):
        PositionCount[i]=0
        for title in sequences.keys():
            if (sequences[title][i] != '-' and sequences[title][i] != 'X') :
                PositionCount[i] += 1

    return(PositionCount)


def filter_problematic_position(limit, counts_per_position, aa_seq_dict, nt_seq_dict):

    problematic_pos=[]
    for position in counts_per_position.keys():
        if counts_per_position[position] < limit :
            
            problematic_pos.append(position) #add one so that we can start counting with position 1
            
    counter = 0
    for pos in problematic_pos:
        #print(len(next(iter(aa_seq_dict.values()))))
        
        for record in aa_seq_dict:
        
            aa_seq_dict[record] = aa_seq_dict[record][:pos-counter] + aa_seq_dict[record][pos +1 -counter:]
            nt_seq_dict[record] = nt_seq_dict[record][:pos*3-counter*3] + nt_seq_dict[record][pos*3 +3 -counter*3:]

        counter +=1
    if (next(iter(aa_seq_dict.values()))):
        return(aa_seq_dict, nt_seq_dict)
    else:
        print('\tno sequence left after all problematic positions were removed ...')
        sys.exit(0)


def filter_aliscore(problematic_pos, aa_seq_dict, nt_seq_dict):
      
    counter = 0
    for pos in problematic_pos:
        #print(len(next(iter(aa_seq_dict.values()))))
        
        for record in aa_seq_dict:
        
            aa_seq_dict[record] = aa_seq_dict[record][:pos-counter] + aa_seq_dict[record][pos +1 -counter:]
            nt_seq_dict[record] = nt_seq_dict[record][:pos*3-counter*3] + nt_seq_dict[record][pos*3 +3 -counter*3:]

        counter +=1
    if (next(iter(aa_seq_dict.values()))):
        return(aa_seq_dict, nt_seq_dict)
    else:
        print('\tno sequence left after all problematic positions were removed ...')
        sys.exit(0)



def filter_short(limit,aa_seq_dict, nt_seq_dict): ###filter out short sequences

    print('\tfiltering short sequences with limit {}'.format(limit))

    SequenceFileDict_long={}
    SequenceFileDict_nt_long={}

    for record in aa_seq_dict.keys() :
        seq = str(aa_seq_dict[record])
        if seq.replace("-", "") and len(seq.replace("-", "")) > limit:
            SequenceFileDict_long[record] = str(aa_seq_dict[record])
            SequenceFileDict_nt_long[record] = str(nt_seq_dict[record])
        else:
            print('\tthis sequences will be removed:\t{}'.format(record))

    return(SequenceFileDict_long,SequenceFileDict_nt_long)

def check_if_aligned(seqs):
    #check if all sequences have the same length
    total_length = len(next(iter(seqs.values())))
    Aligned = True

    for header in seqs.keys():
        if total_length != len(seqs[header]):
            Aligned = False

    if Aligned:
        return total_length
    else:
        return False

def get_arguments():
    parser = argparse.ArgumentParser()
    
    optional = parser._action_groups.pop() 
    required = parser.add_argument_group('required arguments')

    required.add_argument("--aa_input", "-a", help = "aminoacid input filename", required=True)
    required.add_argument("--nt_input", help = "nucleotide input filename", required=True)
    
    
    optional.add_argument("--HE_list", help = "list of taxa that belong to the Hybrid Enrichment set")
    optional.add_argument('--HE_per_position', dest='HE_pos_limit' ,type=int, help = "number of HE taxa per position in order ot keep position")

    optional.add_argument('--simple_header', action='store_true', help = "only use the first three parts of the header (family_genus_species)")
    optional.add_argument('--taxa_per_position', dest='taxa_pos_limit' ,type=int, help = "number of taxa per position in order ot keep position")
    optional.add_argument('--min_aa', type=int, help = "minimal number of aminoacids in a sequence")
    optional.add_argument('--min_ali_length' , type=int , help="minimal length of the alignment after all other filtering steps")
    optional.add_argument('--filter_aliscore', dest='aliscore_file' ,help = "use aliscore results for aminoacid sequence to filter aa and nt sequences")

    parser._action_groups.append(optional) # added this line
    return(parser.parse_args())




if __name__ == "__main__":

    args = get_arguments() 

    print('Aminoacid file: {}\nNucleotide file: {}\nHE-Taxa: {}\n'.format(args.aa_input,args.nt_input,args.HE_list))
    
    SequenceFileDict=SeqIO.to_dict(SeqIO.parse( args.aa_input , "fasta"))
    SequenceFileDict_nt=SeqIO.to_dict(SeqIO.parse( args.nt_input , "fasta"))
    
    output_aa = args.aa_input.replace('.fas', '_filtered.fas')
    output_nt = args.nt_input.replace('.fas', '_filtered.fas')
    

    #pdb.set_trace()
    alignemnt_length_before=''
    alignemnt_length_after=''
    seq_number_before=''
    seq_number_after=''

    if (check_if_aligned(SequenceFileDict) and check_if_aligned(SequenceFileDict_nt)):
        alignemnt_length_before=check_if_aligned(SequenceFileDict)
        seq_number_before=len(SequenceFileDict.keys())
        print('\tThe {} input files are aligned ...\n\ttotal length of alignment:\t{}'.format(seq_number_before,alignemnt_length_before))

    elif(check_if_aligned(SequenceFileDict)):

        print('\n\tThe nucleotide file is not aligned \n Please only use alignements')
        sys.exit(0)
    elif(check_if_aligned(SequenceFileDict_nt)):
        print('\n\tThe aminoacid file is not aligned')
        sys.exit(0)
    else:
        print('\n\tBoth input files are not aligned...\nPlease only use alignements')
        sys.exit(0)


    if(args.simple_header):

        #simplify headers so that they are only Family_Genus_Species - style
        #print('simplify_header')
        SequenceFileDict = simplify_header(SequenceFileDict)
        SequenceFileDict_nt = simplify_header(SequenceFileDict_nt)
    else:
        
        SequenceFileDict = reformat(SequenceFileDict)
        SequenceFileDict_nt = reformat(SequenceFileDict_nt)

    if ( args.aliscore_file ):
        positions = []
        for line in open(args.aliscore_file ,'r'):
            for pos in line.split():
                positions.append(int(pos)-1)
        #print(positions)
        [SequenceFileDict, SequenceFileDict_nt] = filter_aliscore(positions, SequenceFileDict, SequenceFileDict_nt)


    if ( args.HE_list ):
        hetitles = [line.rstrip('\n') for line in open(args.HE_list)]


    #if no he sequences are present we will exit...
        he_seqs = get_HE_seqs(hetitles,SequenceFileDict)
        print('\t{} HE sequences'.format(len(he_seqs.keys())))

        positional_count = count_he_per_position(he_seqs)

        if(args.HE_pos_limit):
            print('\tthe minimal number of HE taxa per position is set to {}'.format(args.HE_pos_limit))
            limit = args.HE_pos_limit
        
        else:
            Number_of_HE_seqs = len(he_seqs.keys())
            limit = int(Number_of_HE_seqs * 0.5 )
            print('\tno limit was set for minimal number of taxa per position. We will continue with minimum of half of the included HE-taxa: {} taxa'.format(limit))


        [SequenceFileDict, SequenceFileDict_nt] = filter_problematic_position(limit ,positional_count, SequenceFileDict, SequenceFileDict_nt )
        print('\tNumber of sequences after filtering low HE regions:\t{}'.format(len(SequenceFileDict.keys())))

    if(args.taxa_pos_limit):

        positional_count_all = count_taxa_per_position(SequenceFileDict)

        [SequenceFileDict, SequenceFileDict_nt] = filter_problematic_position(args.taxa_pos_limit ,positional_count_all, SequenceFileDict, SequenceFileDict_nt )

    #filter short sequences
    if(args.min_aa):
        [SequenceFileDict, SequenceFileDict_nt] = filter_short(args.min_aa, SequenceFileDict, SequenceFileDict_nt)

    if bool(SequenceFileDict):
        alignemnt_length_after=check_if_aligned(SequenceFileDict)
        seq_number_after=len(SequenceFileDict.keys())
        print('\t#Input\tsequences before/after\t\talignment length before/after\n\t##stat\t{}\t{}\t{}\t{}\t{}\n'.format(args.aa_input,seq_number_before,seq_number_after,alignemnt_length_before,alignemnt_length_after))
    else:
        print('no sequences left\n')
        sys.exit(0)

    if(args.min_ali_length and bool(SequenceFileDict) ):
        if (args.min_ali_length > check_if_aligned(SequenceFileDict) ):
            print('total alignment size is not large enough')
            sys.exit(0)

    ##write to output
    with open(output_aa, "w") as output_handle:
        for key in SequenceFileDict.keys():
            output_handle.write('>{}\n{}\n'.format(key,SequenceFileDict[key]))

    with open(output_nt, "w") as nt_output_handle:
        for key in SequenceFileDict_nt.keys():
            nt_output_handle.write('>{}\n{}\n'.format(key,SequenceFileDict_nt[key]))
    
    
