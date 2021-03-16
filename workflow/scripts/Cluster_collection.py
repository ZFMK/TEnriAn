#!/usr/bin/env python3

from Bio import SeqIO
import os
import re

#collect sequences from the different orthograph results and join them in files containing all members of a cluster.


##get all the subdirectories with os.walk(directory)
directory_list = next(os.walk('results/02_orthology/results/'))[1]
#print(directory_list)

#files_list = next(os.walk('.'))[2]
#print(files_list)

#loop through all subdirectories and collect data from the files in aa-subdir and nt-subdir
for directory in directory_list:
    if 'aa' in next(os.walk('results/02_orthology/results/'+directory))[1] and 'nt' in next(os.walk('results/02_orthology/results/'+directory))[1] : #check if aa and nt foleder are there
        print('the path is: {}/aa/'.format(directory))      #print some feedback
        
        #for the AAs
        files_list_sub = next(os.walk('results/02_orthology/results/' + directory +'/aa/'))[2]
        for file in files_list_sub:
            if re.match('.*.fa', file): #do sth for each .fa file...
                #print('\tthis is file: {}'.format(file))
                
                SequenceFileDict=SeqIO.to_dict(SeqIO.parse('results/02_orthology/results/' + directory + '/aa/' + file , "fasta")) #load sequence in dictionary....
                for key in SequenceFileDict.keys():
                    parts=key.split("|")
                    if re.match('.*' + directory + '.*' , parts[1]) and re.match('.*translate.*', key) : #if the sequence corresponds to the species, do sth and it is not from the original Aminoacids...
                        #print('cluster {} sequence {}'.format(parts[0], parts[2]))
                        
                        species_start=parts[1] + '_'
                        if ( parts[2].startswith(species_start)):

                            myfile = open("results/02_orthology/Cluster_collection/"+parts[0] + '.aa.fas', "a+")
                            myfile.write('>{}\n{}\n'.format(parts[2],SequenceFileDict[key].seq))
                            myfile.close()

                        else:
                            myfile = open("results/02_orthology/Cluster_collection/"+parts[0] + '.aa.fas', "a+")

                            myfile.write('>{}\n{}\n'.format(species_start+parts[2],SequenceFileDict[key].seq))
                            myfile.close()
#for the NTs
        print('the path is: {}/nt/'.format(directory))
        files_list_sub = next(os.walk('results/02_orthology/results/' + directory +'/nt/'))[2]
        for file in files_list_sub:
            if re.match('.*.fa', file): #do sth for each .fa file...
                #print('\tthis is file: {}'.format(file))
                
                SequenceFileDict=SeqIO.to_dict(SeqIO.parse('results/02_orthology/results/' +directory + '/nt/' + file , "fasta")) #load sequence in dictionary....
                for key in SequenceFileDict.keys():
                    parts=key.split("|")
                    if re.match('.*' + directory + '.*' , parts[1]): #if the sequence corresponds to the species, do sth
                        #print('cluster {} sequence {}'.format(parts[0], parts[2]))
                        
                        species_start=parts[1] + '_'
                        if ( parts[2].startswith(species_start)):

                            myfile = open("results/02_orthology/Cluster_collection/"+parts[0] + '.nt.fas', "a+")
                            myfile.write('>{}\n{}\n'.format(parts[2],SequenceFileDict[key].seq))
                            myfile.close()
                        else:
                            myfile = open("results/02_orthology/Cluster_collection/"+parts[0] + '.nt.fas', "a+")

                            myfile.write('>{}\n{}\n'.format(species_start+parts[2],SequenceFileDict[key].seq))
                            myfile.close()
                        
#                print(len(SequenceFileDict))
