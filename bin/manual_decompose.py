"""
The MED implementation of the MED algorithm is not working (I think because there is too much error)
so that we only have unique sequences.

I still think that the algorithm in principal should work.

I have dived into their source code and pulled out the entropy calculation formula.
We will manually calculate the entropies and visualize their distributions.

We should have a directory that contains the fastas for each of the pcr_origin and filtered to only contain the target amplicon.
We will aim to make a multi-part figure where each alignment has 2 figures (to start with).
1 - the entropy along the alinment that is essentially a bar chart
2 - a histogram of the column entropies

There will be 19 fastas to work with so we will aim to split the figure into 4 columns and 5 rows
"""

import numpy
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import math
import compress_pickle
class ManMed:
    def __init__(self, long_short="long"):
        if long_short == "long":
            self.path_to_alignments = "/home/humebc/projects/cnifl/results/aligned_multi_sample_single_product_fastas"
        elif long_short == "short":
            raise NotImplementedError
            self.path_to_alignments = "/home/humebc/projects/cnifl/results/aligned_multi_sample_single_product_fastas"
        else:
            raise NotImplementedError
        self.long_short = long_short
        fastas = os.listdir(self.path_to_alignments)
        self.pcrOrigin_to_fasta_dict = {"_".join(_.split("_")[:4]):_ for _ in fastas}
        self.valid_chars = set(['A', 'T', 'C', 'G', '-'])
        nrows = 5
        ncols = 8
        # self.fig, self.ax_arr = plt.subplots(nrows=nrows, ncols=ncols)
        # Currently 2 but will get bigger
        self.num_plots_per_alignment = 2

        # Get all of the fastas into a dict
        self.pcrOrigin_to_fasta_list = {}
        i = 0
        for name, fasta in self.pcrOrigin_to_fasta_dict.items():
            fig, ax = plt.subplots(nrows=1, ncols=1)
            # col_index_ax_1 = (i%ncols)*self.num_plots_per_alignment
            # col_index_ax_2 = col_index_ax_1 + 1
            # row_index_ax_1_2 = math.floor((i*2)/ncols)
            # ax_1 = self.ax_arr[row_index_ax_1_2, col_index_ax_1]
            # ax_2 = self.ax_arr[row_index_ax_1_2, col_index_ax_2]
            self._do_the_work(name, fasta, ax)
            i += 1


    def _do_the_work(self, name, fasta, ax_1):
        fasta_cache_path = f"cache/{name}.fasta_df.gz"
        if os.path.exists(fasta_cache_path):
            print(f"Reading in cache: {fasta_cache_path}")
            self.a_df = compress_pickle.load(fasta_cache_path)
            print("Done")
        else:
            print(f"Reading fasta: {fasta}")
            with open(os.path.join(self.path_to_alignments, fasta), "r") as f:
                self.current_fasta_list = self.convert_interleaved_to_sequencial_fasta([_.rstrip() for _ in f])
            print("Done")
            self.names = [_[1:] for _ in self.current_fasta_list if _.startswith(">")]
            print("Creating df")
            self.a_df = pd.DataFrame([list(_) for _ in self.current_fasta_list if not _.startswith(">")])
            self.a_df.index = self.names
            # remove columns with only "-"
            self.a_df = self.a_df.loc[:, (self.a_df != "-").any(axis=0)]
            print("Done")
            print(f"Cacheing out to: {fasta_cache_path}")
            compress_pickle.dump(self.a_df, fasta_cache_path)
            print("Done")
        
        
        
        # for each columns/position of the alignment calc the entropy
        entropy_cache_path = f"cache/{name}.columnsEntropies.gz"
        if os.path.exists(entropy_cache_path):
            print(f"reading entropy cache: {entropy_cache_path}")
            entropy_list = compress_pickle.load(entropy_cache_path)
            print("Done")
        else:
            entropy_list = []
            print("Calculating entropies")
            indi_entropies = []
            i = 0
            tot_cols = len(list(self.a_df))
            for col in list(self.a_df):
                i += 1
                indi_entropy, col_entropy = self.entropy(list(self.a_df[col]))
                entropy_list.append(col_entropy)
                indi_entropies.extend(indi_entropy)
                if i%100 == 0:
                    print(f"{i} out of {tot_cols}")
            print("Done")
            print(f"Cacheing out to: {entropy_cache_path}")
            compress_pickle.dump(entropy_list, entropy_cache_path)
            print("Done")

        foo = "bar"
        entropy_cutoff = 1.5
        # now plot a histogram

        print("plotting entropy histogram")
        ax_1.hist(entropy_list, bins=int(len(list(self.a_df))/50))
        ax_1.set_yscale('log')
        ax_1.set_title(f"{name} entropy histogram", fontsize=6)
        ax_1.set_xlim(0,2)
        # ax_1.xaxis.set_visible(False)
        print("Done")

        # print("plotting entropy across alignment")
        # ax_2.bar(list(self.a_df), height=entropy_list, width=1)
        # ax_2.set_title(f"{name} across alignment entropy", fontsize=6)
        # ax_2.set_xlabel("BP in alignment")
        # # ax_2.xaxis.set_visible(False)
        # # ax_2.hlines(entropy_cutoff, xmin=ax_2.get_xlim()[0], xmax=ax_2.get_xlim()[1], color="red")
        # print("Done")
        # # plt.tight_layout()
        print("saving")
        plt.savefig(f"testing_figures/{name}_entropies.{self.long_short}.png", figsize=(10, 5))
        print("Done")

        # # Now make a corrected alignment and write this out using the cutoff provided above
        # # we can do this either based on the column entropies, or we can do it on the character entropies
        # # TODO think about implemeneting correction based on character entropies.
        # # However, this is not simple as what do you correct the to? We don't
        # # want to artifically inflate the number of haplotypes by correcting error bases
        # # to the wrong base.
        # indices_to_revert_to_consensus = [index for index, entr in enumerate(entropy_list) if entr <= entropy_cutoff]
        
        # print(f"Converting {len(indices_to_revert_to_consensus)} out of {len(entropy_list)} positions back to consensus")
        # print(f"{len(entropy_list) - len(indices_to_revert_to_consensus)} positions remain")

        # # We also want to know how many unique sequences there are before and after the correction
        # pre_correction_seqs = []
        # for samp, ser in self.a_df.iterrows():
        #     pre_correction_seqs.append("".join(ser.values).replace("-", ""))

        # # # Run an exercise to see how the number of unique sequences decrease as a function of the entropy cutoff
        # # for ent_cutoff in np.arange(1.5,2.5,0.1):
        # #     temp_df = self.a_df
        # #     cor_inds = [index for index, entr in enumerate(entropy_list) if entr <= ent_cutoff]
        # #     for ind in cor_inds:
        # #         temp_df[ind] = self.a_df[ind].value_counts().idxmax()
        # #     post_cor = []
        # #     for samp, ser in temp_df.iterrows():
        # #         post_cor.append("".join(ser.values).replace("-", ""))
        # #     print(f"At entropy cutoff correction at {ent_cutoff} there are {len(set(post_cor))} unique seqs out of {len(pre_correction_seqs)} sequences")

        

        # print(f"Before correction there are {len(set(pre_correction_seqs))} unique seqs out of {len(pre_correction_seqs)} sequences")

        # for ind in indices_to_revert_to_consensus:
        #     self.a_df[ind] = self.a_df[ind].value_counts().idxmax()
        
        # # Remove columns that are only "-"
        # self.a_df = self.a_df.loc[:, (self.a_df != "-").any(axis=0)]
        # # Then write back out.
        # with open(f"{sys.argv[1].replace('.fasta', f'.errorCorrected.{entropy_cutoff}.fasta')}", "w") as f:
        #     for ind, ser in self.a_df.iterrows():
        #         f.write(f">{ind}\n")
        #         f.write(f"{''.join(ser.values)}\n")


        # post_correction_seqs = []
        # for samp, ser in self.a_df.iterrows():
        #     post_correction_seqs.append("".join(ser.values).replace("-", ""))
        # print(f"After correction there are {len(set(post_correction_seqs))} unique seqs out of {len(pre_correction_seqs)} sequences")
            


    def convert_interleaved_to_sequencial_fasta(self, fasta_as_list):
        new_fasta = []
        temp_seq_string_list = []
        for i, fasta_line in enumerate(fasta_as_list):
            if fasta_line.startswith('>'):
                if temp_seq_string_list:
                    new_fasta.append(''.join(temp_seq_string_list))
                    temp_seq_string_list = []
                    new_fasta.append(fasta_line)
                else:
                    new_fasta.append(fasta_line)
            else:
                temp_seq_string_list.append(fasta_line)
        new_fasta.append(''.join(temp_seq_string_list))
        return new_fasta

    def entropy(self, l):
        """
        calulates the entropy for a single position in analignment.
        l: list of characters found at the position
        returns: returns a float that is the entropy
        """
        l = [_.upper() for _ in l]

        character_entropies = []
        for char in self.valid_chars:
            character_proportion = (l.count(char) / len(l)) + 0.0000000000000000001
            character_entropies.append((character_proportion * numpy.lib.scimath.log2(character_proportion)) * -1)
    
        return character_entropies, (sum(character_entropies))

ManMed()