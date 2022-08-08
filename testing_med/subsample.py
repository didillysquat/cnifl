"""
Sample the input fasta keeping an equal number of sequences per sample in the file
"""

import sys
from collections import defaultdict
path_to_fasta = sys.argv[1]
with open(path_to_fasta, "r") as f:
    fasta_lines = [_.rstrip() for _ in f]

samples = [_[1:].split("_")[0] for _ in fasta_lines if _.startswith(">")]

target_seqs = 5000

sample_count_dict = defaultdict(int)

subbed_fasta = []
for i, line in enumerate(fasta_lines):
    if line.startswith(">"):
        sample = line[1:].split("_")[0]
        if sample_count_dict[sample] <= target_seqs:
            # Then we are still collecting sequences for this samples
            subbed_fasta.extend([line, fasta_lines[i+1]])
            sample_count_dict[sample] += 1
        else:
            pass

with open(path_to_fasta.replace(".fasta", ".subsampled.fasta"), "w") as f:
    for line in subbed_fasta:
        f.write(f"{line}\n")

