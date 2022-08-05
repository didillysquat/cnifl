"""
The MED implementation of the MED algorithm is not working (I think because there is too much error)
so that we only have unique sequences.

I still think that the algorithm in principal should work.

I have dived into their source code and pulled out the entropy calculation formular.
I propose that manually do the decomposisiont.
First I was to plot up with the entropy looks like across the alignment just as he does in his outputs
SO that is what I am going to do here.
"""

import numpy
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

class ManMed:
    def __init__(self, collapse_by="character"):
        self.collapse_by=collapse_by
        with open(sys.argv[1], "r") as f:
            self.alignment_list = self.convert_interleaved_to_sequencial_fasta([_.rstrip() for _ in f])
        self.names = [_ for _ in self.alignment_list if _.startswith(">")]
        self.a_df = pd.DataFrame([list(_) for _ in self.alignment_list if not _.startswith(">")])
        self.a_df.index = self.names
        # remove columns with only "-"
        self.a_df = self.a_df.loc[:, (self.a_df != "-").any(axis=0)]

        self.valid_chars = set(['A', 'T', 'C', 'G', '-'])

    def calculate_entropies_of_alignment(self):
        # for each columns/position of the alignment calc the entropy
        entropy_list = []
        indi_entropies = []
        for col in list(self.a_df):
            indi_entropy, col_entropy = self.entropy(list(self.a_df[col]))
            entropy_list.append(col_entropy)
            indi_entropies.extend(indi_entropy)

        foo = "bar"
        entropy_cutoff = 1.5
        # now plot a histogram
        ax = plt.subplot()
        ax.hist(entropy_list, bins=int(len(list(self.a_df))/50))
        ax.set_yscale('log')
        ax.set_title(f"{sys.argv[1].split('/')[-1]} entropy histogram columns")
        plt.savefig(f"{sys.argv[1].split('/')[-1].split('.')[1]}test_histogram.columns.png")
        plt.close()
        print("Done plotting columns hist")

        # Now plot a histogram of the inidividual entropies
        ax = plt.subplot()
        ax.hist(indi_entropies, bins=int(len(list(self.a_df))/50))
        ax.set_yscale('log')
        ax.set_title(f"{sys.argv[1].split('/')[-1]} entropy histogram characters")
        plt.savefig(f"{sys.argv[1].split('/')[-1].split('.')[1]}test_histogram.character.png")
        plt.close()
        print("Done plotting character hist")
        


        # We also want to plot up the entropies along the alignment
        ax = plt.subplot()
        ax.bar(list(self.a_df), height=entropy_list, width=1)
        ax.set_title(f"{sys.argv[1].split('/')[-1]} entropy along alignment", fontsize=6)
        ax.set_xlabel("BP in alignment")
        ax.hlines(entropy_cutoff, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], color="red")
        plt.tight_layout()
        plt.savefig(f"{sys.argv[1].split('/')[-1].split('.')[1]}_entropy.png")

        
        # Now make a corrected alignment and write this out using the cutoff provided above
        # we can do this either based on the column entropies, or we can do it on the character entropies
        # TODO think about implemeneting correction based on character entropies.
        # However, this is not simple as what do you correct the to? We don't
        # want to artifically inflate the number of haplotypes by correcting error bases
        # to the wrong base.
        indices_to_revert_to_consensus = [index for index, entr in enumerate(entropy_list) if entr <= entropy_cutoff]
        
        print(f"Converting {len(indices_to_revert_to_consensus)} out of {len(entropy_list)} positions back to consensus")
        print(f"{len(entropy_list) - len(indices_to_revert_to_consensus)} positions remain")

        # We also want to know how many unique sequences there are before and after the correction
        pre_correction_seqs = []
        for samp, ser in self.a_df.iterrows():
            pre_correction_seqs.append("".join(ser.values).replace("-", ""))

        # # Run an exercise to see how the number of unique sequences decrease as a function of the entropy cutoff
        # for ent_cutoff in np.arange(1.5,2.5,0.1):
        #     temp_df = self.a_df
        #     cor_inds = [index for index, entr in enumerate(entropy_list) if entr <= ent_cutoff]
        #     for ind in cor_inds:
        #         temp_df[ind] = self.a_df[ind].value_counts().idxmax()
        #     post_cor = []
        #     for samp, ser in temp_df.iterrows():
        #         post_cor.append("".join(ser.values).replace("-", ""))
        #     print(f"At entropy cutoff correction at {ent_cutoff} there are {len(set(post_cor))} unique seqs out of {len(pre_correction_seqs)} sequences")

        

        print(f"Before correction there are {len(set(pre_correction_seqs))} unique seqs out of {len(pre_correction_seqs)} sequences")

        for ind in indices_to_revert_to_consensus:
            self.a_df[ind] = self.a_df[ind].value_counts().idxmax()
        
        # Remove columns that are only "-"
        self.a_df = self.a_df.loc[:, (self.a_df != "-").any(axis=0)]
        # Then write back out.
        with open(f"{sys.argv[1].replace('.fasta', f'.errorCorrected.{entropy_cutoff}.fasta')}", "w") as f:
            for ind, ser in self.a_df.iterrows():
                f.write(f">{ind}\n")
                f.write(f"{''.join(ser.values)}\n")


        post_correction_seqs = []
        for samp, ser in self.a_df.iterrows():
            post_correction_seqs.append("".join(ser.values).replace("-", ""))
        print(f"After correction there are {len(set(post_correction_seqs))} unique seqs out of {len(pre_correction_seqs)} sequences")
            


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

ManMed().calculate_entropies_of_alignment()