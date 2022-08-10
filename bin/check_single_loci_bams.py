"""
Here I want to check that the long read bams that have been intersected with the reference CDS mapped reads
only contain 1 loci. I.e. that the PCRs that were done didn't amplify other loci.

I will try to do this using pysam
"""

import pysam
import os
from collections import defaultdict
class SingleLoci:
    def __init__(self):
        
        self.path_to_mapped_longs = "/home/humebc/projects/cnifl/results/long_intersect_ref"
        for bam in [_ for _ in os.listdir(self.path_to_mapped_longs) if _.endswith("bam")]:
            ref_scaf_dict = defaultdict(int)
            ref_scaf_pos_dict = defaultdict(int)
            i = 0
            # print(f"\n\n\nChecking bam: {bam}")
            
            # with pysam.AlignmentFile(os.path.join(self.path_to_mapped_longs, bam), "rb") as bam_f:
            with pysam.AlignmentFile("/home/humebc/projects/cnifl/reference/aip.ficolin.cds.sorted.bam", "rb") as bam_f:
                
                tot = bam_f.mapped
                for bam_entry in bam_f:
                    print(bam_entry.qname)
                    ref_name = bam_f.header["SQ"][bam_entry.reference_id]["SN"]
                    ref_scaf_dict[ref_name] += 1
                    print(f"{ref_name}:{bam_entry.reference_start}-{bam_entry.reference_end}")
                    i += 1
                    # if i%100 == 0:
                    #     print(f"{i}/{tot}")
                    # # if i > 1000:
                    # #     break
                    # ref_name = bam_f.header["SQ"][bam_entry.reference_id]["SN"]
                    # ref_scaf_dict[ref_name] += 1
                    # ref_scaf_pos_dict[f"{ref_name}:{bam_entry.reference_start}-{bam_entry.reference_end}"] += 1
                    # foo = "bar"
                print("\n\n Ref scaff dict:")
                for k, v in ref_scaf_dict.items():
                    print(f"{k}: {v}")
                
                print("\n\n Ref scaff post dict:")
                for k, v in ref_scaf_pos_dict.items():
                    print(f"{k}: {v}")
                print("\n\n\n\n")
                foo = "bar"
            foo = "bar"
            
SingleLoci()
"""
Getting a bit of a strange result with the intersect with the reference.
"""