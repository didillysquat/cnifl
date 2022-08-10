"""This script will read in the mmseqs table and us it to send each of the seuqences in the fq.gz to a new fasta"""

import sys
import gzip
from Bio import SeqIO
from collections import defaultdict
import matplotlib.pyplot as plt

class SplitFasta:
    def __init__(self):
        sample = sys.argv[3]
        origin = sys.argv[4]
        with open(sys.argv[1], "r") as f:
            self.alignment_file = f.read().splitlines() 

        # go through the alignment file and pull out the best score for each of the lines
        best_match_dict = {}
        best_score_dict = {}
        pid_score_dict = defaultdict(int)
        len_match_dict = defaultdict(int)
        for line in self.alignment_file:
            items = line.split("\t")
            q_seq = items[0]
            score = int(items[11])
            if score >= 500:
                if q_seq in best_score_dict:
                    if best_score_dict[q_seq] < score:
                        best_match_dict[q_seq] = items[1]
                        best_score_dict[q_seq] = score
                        pid_score_dict[q_seq] = float(items[2])
                        len_match_dict[q_seq] = abs(int(items[7]) - int(items[6]))
                    else:
                        continue
                else:
                    best_match_dict[q_seq] = items[1]
                    best_score_dict[q_seq] = score
                    pid_score_dict[q_seq] = float(items[2])
                    len_match_dict[q_seq] = abs(int(items[7]) - int(items[6]))

        # Get some stats for the scores associated with each of the references
        plt.hist(best_score_dict.values(), bins=100)
        plt.title(f"{sample}_{origin} BIT score distribution from best matches")
        plt.savefig("bit_score_hist.png")
        plt.close()

        plt.hist(pid_score_dict.values(), bins=100)
        plt.title(f"{sample}_{origin} PID score distribution from best matches")
        plt.savefig("pid_hist.png")
        plt.close()
        
        # We should have 5 mapping for D_cDNA, and we do see 5 fastas with the highest amount of sequences mapped
        # but there are also mappings to the other genes. So... I think we need to have a look to see if there are
        # discrete differences in the distribution of PIDs and BIT scores for each of the target genes
        for match in list(set(best_match_dict.values())):
            bit_scores = [v for k, v in best_score_dict.items() if best_match_dict[k] == match]
            plt.hist(bit_scores, bins=100)
            plt.title(f"{sample}_{origin}_{match} BIT score distribution from best matches")
            plt.savefig(f"{match}_bit_hist.png")
            plt.close()
            pid_scores = [v for k, v in pid_score_dict.items() if best_match_dict[k] == match]
            plt.hist(pid_scores, bins=100)
            plt.title(f"{sample}_{origin}_{match} PID score distribution from best matches")
            plt.savefig(f"{match}_pid_hist.png")
            plt.close()
            len_scores = [v for k, v in len_match_dict.items() if best_match_dict[k] == match]
            plt.hist(len_scores, bins=100)
            plt.title(f"{sample}_{origin}_{match} length score distribution from best matches")
            plt.savefig(f"{match}_len_hist.png")
            plt.close()

        # Now we have a dictionary that has a the sequences mapped to the best match
        fasta_dict = defaultdict(list)
        with gzip.open(sys.argv[2], "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                foo = "bar"
                try:
                    fasta_dict[best_match_dict[record.id]].extend([f">{record.id}", str(record.seq)])
                except KeyError:
                    fasta_dict["missing"].extend([f">{record.id}", record.seq])

        # Here we have the fastas ready for output
        for k, v in fasta_dict.items():
            with gzip.open(f"{sample}_{origin}_{k}.fa.gz", 'wb') as f:
                for line in v:  
                    f.write(f"{line}\n".encode())

        foo = "bar"

SplitFasta()