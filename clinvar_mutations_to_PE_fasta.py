# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 23:22:33 2020

@author: Jan Baijens, janbaijens@live.nl

This script converts Pathogenic ClinVar mutations into input for DeepPE.
The sequence around the edit is scanned for PAM presence, targets are stored in a table and an input file for DeepPE is created.
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
import matplotlib.pyplot as plt
from time import gmtime, strftime
plt.style.use('ggplot')
date = strftime('%m-%d', gmtime())
#Setup working environment and paths:
current_dir = os.getcwd() #Make sure to run from Scripts directory.
PE_efficiency_dir = Path(os.path.split(current_dir)[0]) #Move up one directory
Figures_dir = PE_efficiency_dir / 'Figures'
Files_dir = PE_efficiency_dir / 'Files'
Data_dir = PE_efficiency_dir / 'Data'

#Define functions:
def insertion_edits(mut, hel, pos):
    #Used to convert SPDI format to edit format from DeepPE
    insertions = {'A', 'C', 'G', 'T', 'AG', 'AGGAA', 'AGGAATCATG'} #Not used right now, could check in place if edit will be predicted.
    edits = []
    if mut == '':
        edit = '+' + str(pos) + ', ' + hel + ' insertion'
        edits.append(edit)
        return edits
    ins_len = len(hel) - len(mut)
    for i in range(len(hel)):
        edited = hel[:i] + hel[i+ins_len:]
        if edited == mut: # The deleted part in hel should be added in mut to get hel.
            edit =  '+' + str(pos + i) + ', ' + hel[i:i+ins_len] + ' insertion'
            edits.append(edit)
    return edits

def deletion_edits(mut, hel, pos):
    #Used to convert SPDI format to edit format from DeepPE
    edits = []
    if hel == '':
        edit = '+' + str(pos) + ', ' + str(len(mut)) + ' bp deletion'
        edits.append(edit)
        return edits
    del_len = len(mut) - len(hel)
    for i in range(len(mut)):
        edited = mut[:i] + mut[i + del_len:]
        if edited == hel:
            edit = '+' + str(pos + i) + ', ' + str(del_len) + ' bp deletion'
            edits.append(edit)
    return edits

def get_desired_edit(edit_type, edit, position):
    #Used to convert SPDI format to edit format from DeepPE
    desired_edits = set()
    if edit_type == 'SNP':
        desired_edit = '+' + str(position) + ', ' + str(edit)
        desired_edits.add(desired_edit)
            
    elif edit_type == 'deletion':
        mut, hel = edit.split(' to ')
        desired_edits = desired_edits.union(deletion_edits(mut, hel, position))
    
    elif edit_type == 'insertion':
        mut, hel = edit.split(' to ')
        desired_edits = desired_edits.union(insertion_edits(mut, hel, position))
    return desired_edits

def find_pams(sequence, begin, end, pams, pam_len = 3):
    '''
    Search the sequence from begin to end, return list of all found pam locations.
    '''
    found_pam_locations = []
    for i in range(begin, end):
        scanner = sequence[i:i+pam_len]
        if scanner in pams:
           found_pam_locations.append(i)
    return found_pam_locations

def sequence_to_PE_input(sequence, pam_loc, edit_loc):
    #Given a sequence with a known PAM location and edit location:
    #Return a 46 bp sequence where the PAM has a fixed position and return the relative position of the edit.
    start = pam_loc - 24
    end = pam_loc + 23
    input_seq = sequence[start : end]
    edit_position = 4 + (edit_loc - pam_loc)
    return input_seq, edit_position

def transversion(nucl):
    transversions = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    if nucl in transversions:
        return transversions[nucl]
    
def three_edits(nucl, pos):
    edits = set()
    for other in ['A', 'C', 'G', 'T']:
        edit = '+{}, {} to {}'.format(pos, nucl, other)
        if nucl != other:
            edits.add(edit)
    return edits

def get_predictable_edits(seq):
    #Determine which types of edits PE_Position and PE_Type will be able to predict for a given input sequence.
    SNP_edit_pos = set([1,2,3,4,6,7,8,9,11,14])
    #Add insertions and deletions:
    predictable_edits = {'+1, 1 bp deletion',
'+1, 10 bp deletion',
'+1, 2 bp deletion',
'+1, 5 bp deletion',
'+1, A insertion',
'+1, AG insertion',
'+1, AGGAA insertion',
'+1, AGGAATCATG insertion',
'+1, C insertion',
'+1, G insertion',
'+1, T insertion'}
    #Add all edits on position +1:
    predictable_edits = predictable_edits.union(three_edits(seq[21], 1))
    #Add all edits from PE_position:
    pos = 0
    for nucl in seq[21:35]:
        pos += 1
        if pos in SNP_edit_pos:
            trans = transversion(nucl)
            edit = '+{}, {} to {}'.format(pos, nucl, trans)
            predictable_edits.add(edit)
    return predictable_edits

#Load data:
mutations_file = str(Data_dir / '01-03Clinvar_mutations.csv')
editable_df_file = str(Data_dir / '03-14editable.csv')
edits_file = str(Data_dir / '03-13Clinvar_edits.csv')
mutations = pd.read_csv(mutations_file, low_memory = False)

#Change Clinvar mutations to editable mutations, for both strands:
editable_muts = []
for i in range(len(mutations)):
    if i%10000 == 0:
        print('Working on row: {}'.format(i))
    row = mutations.iloc[i].copy()
    row['strand'] = 0
    row['edit'] = ''
    row['mut_loc'] = 0
    #Forward strand:
    editables = get_desired_edit(row.edit_type, row.edit_f, row.mut_loc_f)
    for editable in editables:
        pos, edit = editable.split(', ')
        new_row = row.copy()
        new_row.strand = 'forward'
        new_row.mut_loc = int(pos.replace('+', ''))
        new_row.edit = edit
        editable_muts.append(new_row)
    #Reverse strand:
    editables = get_desired_edit(row.edit_type, row.edit_r, row.mut_loc_r)
    for editable in editables:
        pos, edit = editable.split(', ')
        new_row = row.copy()
        new_row.strand = 'reverse'
        new_row.mut_loc = int(pos.replace('+', ''))
        new_row.edit = edit
        editable_muts.append(new_row)
    
editable_df = pd.concat(editable_muts, axis = 1).T

#Save editable_df to file:
editable_df.to_csv(editable_df_file, index = False)

#Add non editable mutations, we do want to count their PAMs.
editable_SPDI_set = set(editable_df['Canonical SPDI'])
mutations_SPDI_set = set(mutations['Canonical SPDI'])
non_editable_mutations = mutations_SPDI_set.difference(editable_SPDI_set)

new_edits = []
for mut in non_editable_mutations:
    row = mutations.loc[mutations['Canonical SPDI'] == mut].iloc[0].copy()
    row['strand'] = 'forward'
    row['edit'] = np.nan
    row['mut_loc'] = row.mut_loc_f
    new_edits.append(row.copy())
    
    row['strand'] = 'reverse'
    row['edit'] = np.nan
    row['mut_loc'] = row.mut_loc_r
    new_edits.append(row.copy())
    
new_edits_df = pd.concat(new_edits, axis = 1).T
editable_df = editable_df.append(new_edits_df)

#Set PAM search window parameters:
begin = -10 #Relative to starting location of mutation in sequence.
end = 4 #Relative to starting location of mutation in sequence. #TRY WITH 4.
pam_names = ['NGG', 'NAN', 'NGN']
pams = [{'AGG', 'CGG', 'GGG', 'TGG'}, #PAMs that will be accepted.
        {'AAA','AAC','AAG','AAT','CAA','CAC','CAG','CAT','GAA','GAC','GAG','GAT','TAA','TAC','TAG','TAT'}, 
        {'AGA','AGC','AGG','AGT','CGA','CGC','CGG','CGT','GGA','GGC','GGG','GGT','TGA','TGC','TGG','TGT'}] 
edit_positions_dict = dict()
#Loop over mutations, store all found PAM locations:
for i in range(len(pam_names)):
    pam = pams[i]
    pam_name = pam_names[i]
    
    edit_positions = []
    for i in range(len(editable_df)):
        if i%33333 == 0:
            print('Working on row: {}'.format(i))
            
        row = editable_df.iloc[i].copy()
        if row.strand == 'forward':
            #Search PAMs in forward strand:
            found_pam_locations = find_pams(row.seq_mutated_f, row.mut_loc + begin, row.mut_loc + end, pam)
            for pam_loc in found_pam_locations:
                row['strand'] = 'forward'
                row['pam_location'] = pam_loc
                row['input_seq'], row['edit_position'] = sequence_to_PE_input(row.seq_mutated_f, pam_loc, row.mut_loc)
                edit_positions.append(row.copy())
        elif row.strand == 'reverse':
            #Search PAMs in reverse strand:
            found_pam_locations = find_pams(row.seq_mutated_r, row.mut_loc + begin, row.mut_loc + end, pam)
            for pam_loc in found_pam_locations:
                row['strand'] = 'reverse'
                row['pam_location'] = pam_loc
                row['input_seq'], row['edit_position'] = sequence_to_PE_input(row.seq_mutated_r, pam_loc, row.mut_loc)
                edit_positions.append(row.copy())
    edit_positions_dict[pam_name] = pd.concat(edit_positions, axis = 1).T

#Save edit_positions to file for each PAM search:
for pam in pam_names:
    edit_positions_df = edit_positions_dict[pam]
    file_name = '{}edits_{}.csv'.format(date, pam)
    edit_positions_df.to_csv(str(Data_dir / file_name), index = False)
    
#edit_positions_df = pd.concat(edit_positions, axis = 1).T

##Reload data: TEMP
#editable_df = pd.read_csv(editable_df_file, low_memory = False)
#edit_positions_df = pd.read_csv(str(Data_dir / '03-14edits_{}.csv'.format('NGG')), low_memory = False)

##Sanity check on input sequences:
#input_seq_lens = dict(collections.Counter([len(s) for s in list(edit_positions_df.input_seq)])) #Should be length 47
#incorrect_input_seq_alphabet = set([s for s in list(edit_positions_df.input_seq) if len(set(s)) > 4])
#edit_positions_df = edit_positions_df.loc[~edit_positions_df.input_seq.isin(incorrect_input_seq_alphabet)].copy()

##Check for duplicates caused by corrected mut_loc, multiple PAM search:
#unique_edits = set()
#duplicate_i = []
#for i in range(len(edit_positions_df)):
#    row = edit_positions_df.iloc[i]
#    edit_data = (row.input_seq, row.edit_position, row.edit)
#    if edit_data in unique_edits:
#        duplicate_i.append(i)
#    else:
#        unique_edits.add(edit_data)
        
#Add desired_edit, predictable_edit, predictable cols:
desired_edits = []
predictable_edits = []
predictables = []
for i in range(len(edit_positions_df)):
    row = edit_positions_df.iloc[i]
    desired = '+{}, {}'.format(row.edit_position, row.edit)
    predictable_edit = get_predictable_edits(row.input_seq)
    if desired in predictable_edit:
        predictable = 1
    else:
        predictable = 0
    desired_edits.append(desired)
    predictable_edits.append(predictable_edit)
    predictables.append(predictable)

edit_positions_df['desired_edit'] = desired_edits
edit_positions_df['predictable_edits'] = predictable_edits   
edit_positions_df['predictable'] = predictables 
      
#Reindex dataframe:
edit_positions_df.reset_index(inplace = True, drop = True)
edit_positions_df['input_name'] = edit_positions_df.index
edit_positions_df.to_csv(str(Data_dir / '03-14edits_NGG.csv'), index = False)

#Get predictable df:
edits_predictable = edit_positions_df.loc[edit_positions_df.predictable == 1].copy()
edits_predictable.to_csv(str(Data_dir / '03-14edits_NGG_predictable.csv'), index = False)

##Reload data: TEMP
#editable_df = pd.read_csv(editable_df_file, low_memory = False)
#edit_positions_df = pd.read_csv(str(Data_dir / '03-14edits_{}.csv'.format('NGG')), low_memory = False)
#edits_predictable = pd.read_csv(str(Data_dir / '03-14edits_NGG_predictable.csv'), low_memory = False)

##Analyse number of found PAMs per mutation
#by_mutation = edit_positions_df.groupby('Canonical SPDI')
#mutation_counts = by_mutation.count()['Name']
#nr_PAMs_frequency = dict(collections.Counter(list(mutation_counts)))
#nr_PAMs_frequency[0] = len(set(mutations['Canonical SPDI']).difference(set(edit_positions_df['Canonical SPDI'])))
#
#xy = sorted(list(nr_PAMs_frequency.items()))
#plt.bar(x = [e[0] for e in xy], height = [e[1] for e in xy])
#plt.xticks([0,2,4,6,8,10,12,14,16])
#plt.xlabel('Nr found NGN PAMs')
#plt.ylabel('Frequency')
#plt.title('Number of NGN PAMs per Clinvar pathogenic mutation')
##plt.savefig(str(Figures_dir / '02-26PAM_count_NGN.png'), dpi = 600)
#plt.show()
#
#dict(collections.Counter(edit_positions_df.edit_position))
#plt.hist(edit_positions_df.edit_position)

##Reload data:
#edit_positions_df = pd.read_csv(edits_file, low_memory = False)

##Create random samples, write input to fasta files:
#random_subset = np.random.choice(range(len(edits_predictable)), 250, replace = False)
#random_subset_SNP = np.random.choice(list(edits_predictable.loc[edits_predictable.edit_type == 'SNP'].index), 250, replace = False)
#random_subset_del = np.random.choice(list(edits_predictable.loc[edits_predictable.edit_type == 'deletion'].index), 250, replace = False)
#random_subset_ins = np.random.choice(list(edits_predictable.loc[edits_predictable.edit_type == 'insertion'].index), 250, replace = False)
#fasta_file = str(Data_dir / '04-26predictable_edits_sample250_del.fasta')
#with open(fasta_file, 'a+') as handle:
##    for i in range(len(edit_positions_df)):
#    for i in random_subset_del:
#        row = edits_predictable.iloc[i]
##        if row['edit_type'] == 'SNP':
#        id_ = str(int(row['input_name']))
#        seq = row['input_seq']
#        handle.write('>' + id_ + '\n' + seq + '\n')

#http://deepcrispr.info/DeepPE/page/user/a55a5d370501291e045ae57dde164e4d9bacd20d #Result random sample
#http://deepcrispr.info/DeepPE/page/user/96f2f87aa47c3e82a04ee4d42267258d22952cd2 #Result SNP sample
#http://deepcrispr.info/DeepPE/page/user/fa2f891416bc7e94b7217dc3c8666c9319df0ebe #Result ins sample
#http://deepcrispr.info/DeepPE/page/user/bf4668cc31e42ed2d50478d5968212de685e7fb7 #Result del sample
        
#Create fasta files with chunks of input:
df = edits_predictable
n_chunks = 10
chunksize = int(len(df) / n_chunks)
for i in range(1, n_chunks + 1): #Create 1-index
    fasta_file = str(PE_efficiency_dir / 'input_chunks\\input_chunk_{}.fa'.format(i))
    print('Writing file: {}'.format(fasta_file))
    n_commands = 0
    with open(fasta_file, 'a+') as handle:
        for j in range(chunksize*int(i-1), chunksize*int(i)):
            row = df.iloc[j]
            id_ = str(int(row['input_name']))
            seq = row['input_seq']
            handle.write('>' + id_ + '\n' + seq + '\n')
            n_commands += 1
    print('Wrote {} entries to file.'.format(n_commands))

        

