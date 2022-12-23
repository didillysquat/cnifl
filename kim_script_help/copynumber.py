
import os
import pandas as pd
import numpy as np
import pysam

# #read in bamfile as alignment file, 'rb' stands for 'reading' + 'bam'
# bamfile = pysam.AlignmentFile("amplicons/barcode08/barcode08_sorted.bam", "rb" )

# #read in shannon.txt file as df with pandas
# SE = pd.read_csv("amplicons/barcode08/barcode08_CC7000007l_2618834-2620854_CC7000007l.shannon.txt", sep = "\t")
# print(SE)

# #only select rows where the Shannon entropy > 0.8 and positions are within blast hit for exon
# SNP = SE.loc[(SE["shannon_entropy"] > 0.8) & (SE["position"] >= 2618834) & (SE["position"] <= 2620854)]
# print(SNP)

# # create a list from dataframe SNP with all positions between 2618834-2620854 where SE > 0.8
# b = SNP["position"].values.tolist()
# print(b)

# #get coverage and bases in reads at given position
# #indicate contig/region, e.g. "CC7000007l" and range
# #set truncate to true to avoid getting all reads which overlap the region
# #python is 0 based, so base 1 would be base 0 in py (b = base in bam, so b-1 is correct base in py)

# filename = "ben_py/b8_pileup.txt"
# outfile = open(filename, "w")

# for i in b:
#     for pileupcolumn in bamfile.pileup("CC7000007l", (i-1), (i), truncate = True, min_base_quality = 0):
#         for pileupread in pileupcolumn.pileups:
#             if not pileupread.is_del and not pileupread.is_refskip:
#                 # query position is None if is_del or is_refskip is set.
#                 outfile.write('\n%s \t %s \t %s' %
#                     (pileupread.alignment.query_name,
#                     pileupread.alignment.query_sequence[pileupread.query_position], pileupcolumn.pos))

# bamfile.close()
# outfile.close()
filename = "/Users/benjaminhume/Downloads/gapdh_pileup.txt"

# create pandas dataframe from outfile with column headers
infile = pd.read_csv(filename,  sep = "\t", header=None, names=["read", "base", "position"])

# # transform base column into row with position column as corresponding column header
# # save dataframe as txt file
# # for barcode data that returns an error code due to 'duplication' use pivot_table
# barcode = infile.pivot_table(index = 'read', columns='position', values='base', aggfunc="first").reset_index()
# pivoted = infile.pivot(index='read', columns='position', values='base').reset_index()

# #pivoted.columns.name=None
# pivoted = barcode
# pivoted.to_csv("ben_py/b8_pivoted.txt", sep='\t')


# transform base column into row with position column as corresponding column header
# save dataframe as txt file

from collections import Counter

pivoted = infile.pivot(
    index='read', columns='position', values='base'
    ).reset_index()
pivoted.set_index("read", inplace=True, drop=True)
pivoted_no_nan = pivoted[~pivoted.isnull().any(axis=1)]
#base_n = len(pivoted_no_nan)
# The no_nan_list will only contain sequences that are 'full length'.
# I.e. it will not contain sequences that don't have the * at the beginning of end.
no_nan_list_long = pivoted_no_nan.apply(lambda x: ''.join(x.values.tolist()).replace(" ", ""), axis=1).values.tolist()
base_n = len(set(no_nan_list_long))
# counter_object = Counter(no_nan_list_long)
# keys = counter_object.keys()
# base_n = len(keys)
pivoted = pivoted.replace(np.nan, "*")

# remove the gaps from the ends
old_seq_list = pivoted.apply(lambda x: ''.join(x.values.tolist()).replace(" ", "").strip("*"), axis=1).values.tolist()

n = base_n

from itertools import combinations_with_replacement
all_combi_seqs = []
# sequences that only have * at the begining or end
no_nan_list_short = []
base_seqs = list(set(no_nan_list_long))
# This is a list of all of the sequences that we have 'counted'
# as a haplotype
counted_sequences = []
# Add the no_nan_list_long sequences to the counted_sequences list
counted_sequences.extend(base_seqs)
# A dictionary to relate each of the asterix containing seuqences
# to their combination sequences
all_combi_seqs_dict = {}
for seq in old_seq_list:
    print(f"Processing {seq}")
    
    # remove the current sequence from the list
    # and check to see if the current sequence is found twice in the list
    # current_sequence = "".join(ser.values.tolist()).replace(" ", "")
    if seq in counted_sequences:
        continue # Then the seuqence is found more than once and we can't be sure that it is different form the other sequences.
    # We are only interested in sequences that had NAs in them
    # 20221122 No, this is not true
    # The sequences that don't have "*" in them, and are not the original base sequences
    # will need to be searched for in the list of sequences that is generated from 
    # the combinations work we do below.

    # For the sequences that don't have *s in them
    if not "*" in seq:
        no_nan_list_short.append(seq)
        continue # Deal with these sequences after we have generated all the combination sequences
    
    # # If not then make a new list without the current seuqences in it so that we can compare to this list.
    # 20221122 This logic is flawed. If two sequences are the same, then we should still count 1 of them.
    # temp_old_seq_list = [_ for _ in old_seq_list if _ != seq]

    # get the number of nan (only those in the middle of the sequence) You are not concered with the ends of the sequences.
    num_nas = seq.count("*")
   
    # Generate all the combinations that those na positions could be
    combins = combinations_with_replacement(["*", "A", "G", "C", "T"], num_nas)
    new_seq_list = []
    # For each combination generate a sequence where the *s are replaced by the combination
    # Add this to a new list that we can then compare to the counted_sequences list
    # If none of the combins sequences are in the counted_sequences list, then we can add this sequence
    # (in its raw form i.e. with the *).

    # Get the indices of the asterixes
    ast_indices = [ind for ind, val in enumerate(seq) if val == "*"]
    for comb in combins:
        new_seq = list(seq)
        for i in range(len(comb)):
            # replace the correct asterix
            new_seq[ast_indices[i]] = comb[i]
        new_seq_list.append("".join(new_seq))
    
    # Add the generated combi seqs to the collection
    all_combi_seqs.extend(new_seq_list)
    all_combi_seqs_dict[seq] = new_seq_list
    all_combi_seqs = list(set(all_combi_seqs))

    # Now need to check if our sequences are subset or superset of the counted_sequences
    # If match was true, then we want to make sure that we retain the longest version of the two sequences that made the match
    # else, we could end up with a sequences like 'A' that will match to every sequences and ruin our ability to count
    matches = []
    for new_seq in new_seq_list:
        for counted_seq in counted_sequences:
            if new_seq in counted_seq or counted_seq in new_seq:
                matches.append(counted_seq)
    
    if not matches:
        n += 1
        counted_sequences.append(seq)
    else:
        # Make sure that the longest sequences in the matches is the one that is retained in the
        # counted sequences list
        if len(seq) > len(max(counted_sequences, key=len)):
            # Then the current seq is longer than the seq that was matched and it should replace it
            counted_sequences[counted_sequences.index(matches[0])] = seq



# Now check the no_nan_list_short seqs whether they are found in the counted_sequences
for no_nan_seq in no_nan_list_short:
    matches = []
    # Check the combi_seq list
    for combi_key, combi_list in all_combi_seqs_dict.items():
        for combi_seq in combi_list:
            if no_nan_seq in combi_seq or combi_seq in no_nan_seq:
                matches.append(combi_key)

    if matches:
        # If match was true, then we want to make sure that we retain the longest version of the two sequences that made the match
        # else, we could end up with a sequences like 'A' that will match to every sequences and ruin our ability to count
        # Make sure that the longest sequences in the matches is the one that is retained in the
        # counted sequences list
        if len(no_nan_seq) > len(max(all_combi_seqs, key=len)):
            # Then the current seq is longer than the seq that was matched and it should replace it
            counted_sequences[counted_sequences.index(max(matches, key=len))] = seq
        continue
    # Also check to see that the sequence is not found in the counted_sequences
    else:
        for counted_seq in counted_sequences:
            if no_nan_seq in counted_seq or counted_seq in no_nan_seq:
                matches.append(counted_seq)        

    if not matches:
        # Then we can count this seq as unique
        n += 1
        counted_sequences.append(seq)
    else:
        if len(seq) > len(max(counted_sequences, key=len)):
            # Then the current seq is longer than the seq that was matched and it should replace it
            counted_sequences[counted_sequences.index(max(matches, key=len))] = seq
        continue



print(f"n is {n}")



foo = "bar"



