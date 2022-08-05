"""
So it turned out that there are multiple products in the PCR we are looking at and I'm going to try
to map them out to the aptasia genome that I have down loaded from here. https://www.ncbi.nlm.nih.gov/assembly/GCF_001417965.1/
I've mapped the 100 subsampled fasta to play with using mapPacBio.sh.

I want to make a histogram of the different contigs where the sequences map
"""
import sys
sam_path = sys.argv[1]


import re
from collections import Counter, defaultdict
from Bio.Seq import Seq

with open(sam_path, "r") as f:
    sam_list = [_.rstrip() for _ in f]

mapped_contigs = []
scaff_to_fasta_list_dict = defaultdict(list)
for line in sam_list:
    seq_name = line.split("\t")[0]
    seq_seq = line.split("\t")[9]
    bit_flag = int(line.split("\t")[1])
    if bit_flag not in [0,16]:
        continue

    scaff = re.findall("scaffold\d+", line)[0]
    if bit_flag == 0:
        scaff_to_fasta_list_dict[scaff].extend([f">{seq_name}", seq_seq])
    elif bit_flag == 16: # The seq needs reverse complementing
        scaff_to_fasta_list_dict[scaff].extend([f">{seq_name}", str(Seq(seq_seq).reverse_complement())])

    mapped_contigs.append(scaff)

for scaff, scaff_fasta_list in scaff_to_fasta_list_dict.items():
    with open(sam_path.replace(".no.header.sam", f".{scaff}.fasta"), "w") as f:
        for line in scaff_fasta_list:
            f.write(f"{line}\n")

# Now write out split
# We know that the reference sequence mapped to scaffold 273 This is the 3rd most abundant match in the sub 100
# so not great but not a complete loss.


print(Counter(mapped_contigs))