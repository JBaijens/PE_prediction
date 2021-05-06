# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 13:35:26 2020

@author: Jan Baijens, janbaijens@live.nl

This script creates a table of human pathogenic mutations for the quantification of predicted prime editing efficiencies.
Current input is downloaded from the NCBI ClinVar database, by pathogenic mutations under 51 bp length.
"""
#Import libraries:
import os
from pathlib import Path
import numpy as np
import pandas as pd
import collections
import matplotlib.pyplot as plt
#from scipy import stats
from Bio import SeqIO, Seq
plt.style.use('ggplot')
#Setup working environment and paths:
current_dir = os.getcwd() #Make sure to run from Scripts directory.
PE_efficiency_dir = Path(os.path.split(current_dir)[0]) #Move up one directory
Figures_dir = PE_efficiency_dir / 'Figures'
Files_dir = PE_efficiency_dir / 'Files'
Data_dir = PE_efficiency_dir / 'Data'

#Define functions:
def get_edit_type(mutated, healthy):
    if len(mutated) == len(healthy) == 1:
        return 'SNP'
    elif len(mutated) < len(healthy):
        return 'insertion'
    elif len(mutated) > len(healthy):
        return 'deletion'
    else:
        return 'other'
    
def lookup_mutation(chromosome, row, left_flank_len, right_flank_len, chromosome_len):
    healthy = row['healthy_nucleotide']
    mutated = row['mutated_nucleotide']
    position = int(row['mutation_location'])
    
    #Check if sequence location corresponds to mutation location:
    healthy_seq = chromosome[position:position+len(healthy)]
    if healthy_seq != healthy:
        print('Mismatch for row {}, sequence {}, {}.'.format(row.index, healthy_seq, healthy))
        return
    
    #Get left and right flanks of mutation:
    left = max(position - left_flank_len, 0)
    right = min(position + len(healthy) + right_flank_len, chromosome_len)
    
    left_flank = chromosome[left:position]
    right_flank = chromosome[position+len(healthy):right]
    
    #Store healthy and mutated sequences with flanks, mutation location and description:
    row['seq_healthy'] = chromosome[left:right]
    row['seq_mutated_f'] = left_flank + mutated + right_flank
    row['mut_loc_f'] = len(left_flank)
    row['seq_mutated_r'] = Seq.reverse_complement(row['seq_mutated_f'])
    row['mut_loc_r'] = len(right_flank)
    row['edit_f'] = mutated + ' to ' + healthy
    row['edit_r'] = Seq.reverse_complement(mutated) + ' to ' + Seq.reverse_complement(healthy)
    row['edit_type'] = get_edit_type(mutated, healthy)
    return row

def func(pct, allvals):
    absolute = int(pct/100.*np.sum(allvals))
    return "{:d}\n({:.1f}%)".format(absolute, pct)

#Read input:
input_file = str(Data_dir / '12-23_clinvar_pathogenic.txt')
pathogenic = []
bad_lines = []
with open(input_file, 'r') as handle:
    header = handle.readline().strip().split('\t')
    for line in handle:
        data = line.strip().split('\t')
        if len(data) == len(header):
            pathogenic.append(data)
        else:
            bad_lines.append(data)

pathogenic_df = pd.DataFrame(pathogenic, columns = header)

#Only keep mutations that are mutated in a single location, filter out others:
incorrect_SPDI = set([SPDI for SPDI in list(pathogenic_df['Canonical SPDI']) if len(SPDI.split(':')) != 4])
pathogenic_df = pathogenic_df.loc[~pathogenic_df['Canonical SPDI'].isin(incorrect_SPDI)].copy()
pathogenic_df.reset_index(inplace = True, drop = True)

#List sequence identifiers of involved mutations using SPDI (https://pubmed.ncbi.nlm.nih.gov/31738401/):
seq_ids, locations, healthy, mutated = list(zip(*[i.split(':') for i in list(pathogenic_df['Canonical SPDI'])])) 
pathogenic_df['RefSeq_accession_number'] = seq_ids
pathogenic_df['mutation_location'] = locations
pathogenic_df['healthy_nucleotide'] = healthy
pathogenic_df['mutated_nucleotide'] = mutated
unique_seq_ids = list(set(seq_ids))

#Write sequence identifiers to file:
with open(str(Files_dir / '10-31seq_ids'), 'w') as handle:
    for i in unique_seq_ids:
        handle.write(i + '\n')
        
#Sequences may now be downloaded to Data directory using Batch Entrez (https://www.ncbi.nlm.nih.gov/sites/batchentrez).
#Sequence fasta file is stored as refseq_human_genome.fasta in PE_efficiency/Data/ folder.
        
#Loop over clinvar mutations, obtain mutation sequences for both strands, along with mutation info and location:
#Get indexes of mutations for each chromosome:
by_acc = pathogenic_df.groupby('RefSeq_accession_number')
chromosome_mutation_indexes = {i : list(by_acc.get_group(i).index) for i in unique_seq_ids}    

#Add new columns to table:
new_cols = ['seq_healthy', 'seq_mutated_f', 'seq_mutated_r', 'mut_loc_f', 'mut_loc_r', 'edit_f', 'edit_r', 'edit_type']
for col in new_cols:
    pathogenic_df[col] = [''] * len(pathogenic_df)

#Set parameters for flanking sequence lengths on each side of mutations:
left_flank_len = 50
right_flank_len = 50

#Loop over fasta file, lookup mutations for each chromosome:
for record in SeqIO.parse(str(Data_dir / 'refseq_human_genome.fasta'), 'fasta'):
    chromosome = record.id
    sequence = str(record.seq)
    sequence_len = len(sequence)
    print('Working on {}, length {}'.format(chromosome, sequence_len))
    if chromosome in chromosome_mutation_indexes:
        for index in chromosome_mutation_indexes[chromosome]:
            row = pathogenic_df.iloc[index].copy()
            pathogenic_df.iloc[index] = lookup_mutation(sequence, row, left_flank_len, right_flank_len, sequence_len)

#Save table to file:
pathogenic_df.to_csv(str(Data_dir / '01-03Clinvar_mutations.csv'), index = False)

#Plot edit types:
edit_types = dict(collections.Counter(list(pathogenic_df.edit_type)))
plt.pie(edit_types.values(), labels = edit_types.keys(), autopct = lambda pct: func(pct, list(edit_types.values())))
plt.axis('equal')
plt.title('Edit type for Clinvar pathogenic mutations <51 bp')
plt.savefig(str(Figures_dir / 'edit_types.png'), dpi = 600)
plt.show()


