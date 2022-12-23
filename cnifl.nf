#!/usr/bin/env nextflow

// After the in silico PCR, and the screen for size, the next issue was that there were still multiple
// products in the PCR that mapped to different regions of the genome.
// So to further break down the PCR into the different amplicons we will map against the aptasia
// genome and create separate fasta file based on the scaffold that they map to.
// We can also work with a 'target scaffold' which is the scaffold that the reference sequence maps to.

// The sequences are really quite messy I will try to do an insilico PCR. This will remove the need for porechop

// This script will handle the inital processing of the cnifl data.
// The sequencing has been done on the minion by Matt and Kim.
// They have used an 8 adapter set to index the four samples, and the gDNA vs cDNA amplicons.
// These have been split into separate collections of fastq file. There are 4000 reads per fastq file
// For each of the indices the PCRs have been pooled (13 different genes.)
// As such we will need to map these out into separate files so that we can run MED on them.

// We want to end up with 1 fastq.gz file for each sample/origin (i.e. gDNA or cDNA)/gene combination
// This should be roughly 4 X 2 X 13 different fastq files. ~ 104 files
// We will then want to combine the samples into the same files for MED, so this will mean that we have
// 2 X 13 fastq.gz files, one for the gDNA and one for the cDNA for each of the 13 genes.


seq_ch = Channel.fromPath('/home/humebc/projects/cnifl/raw_seq_data/concatenated_files/*.fastq.gz').map{[it.getName().split("_")[0] + "_" + it.getName().split("_")[1].replaceAll(".fastq.gz", ""), it]}
ref = file("/home/humebc/projects/cnifl/reference/ref")
ref_bam = file("/home/humebc/projects/cnifl/reference/aip.ficolin.cds.sorted.bam")
ref_bai = file("/home/humebc/projects/cnifl/reference/aip.ficolin.cds.sorted.bam.bai")
split_by_sam = file("/home/humebc/projects/cnifl/testing_med/sam.py")
subsample = file("/home/humebc/projects/cnifl/testing_med/subsample.py")
// individual_pcrs = Channel.from([["A_gdna", "cnifl_4936_1"], ["A_gdna", "cnifl_6313_1"], ["D_cdna", "cnifl_4936_1"]])
// primers = Channel.from([[["TTCACTCCTGCTTGGTGTCC", "CGAGCCAACCCTCAATCTGT"], "cnifl_4936_1", "cnifl3"], [["CATATCTGTGCTCCCGGTCC", "AGGACTTGTGCTTACCGTGG"], "cnifl_6313_1", "cnifl9"]])
// I have included a cutoff that we will use for the minimum size that is half of the product supposed size.
primers = Channel.from([
    [["TTCACTCCTGCTTGGTGTCC",	"CGAGCCAACCCTCAATCTGT"],  "cnifl_4936_1" , 813],
    [["CATATCTGTGCTCCCGGTCC",	"AGGACTTGTGCTTACCGTGG"],  "cnifl_6313_1" , 682],	
    [["CATATCTGTGCTCCCGGTCC",	"CATTGTACCACCAGGCACCT"],  "cnifl_6313_2" , 655],	
    [["GTGCAAGAGGAAAGTCAGGC",	"TGGAGTTCGGAACAGTCTTTTG"],"cnifl_6313_3" , 457],	
    [["TTGGTGGATTTGAATGCCGC",	"TAGTCGTGAACGCCATGTCG"],  "cnifl_11648_1", 850],
    [["TTCGTCCGCTACTTACGCTC",	"GCGTTTGGGCGGATTTTCAT"],  "cnifl_13073_1", 1007],
    [["TGCCATCACTCTATCGGTCG",	"CCACCAGGCACCGGTATAAG"],  "cnifl_15897_1", 919],
    [["GGTCCTTCCTTAGCAAGTCCA",	"CTGATTTGTAAAGTTCAGCGCA"],"cnifl_15897_2", 404],
    [["AACGATGGCCATCAGTTCGA",	"ACCGGTCCATGATGTAGCAC"],  "cnifl_15899_1", 268],
    [["CCTACGCGCTCTTCATCAGT",	"CCAGCTGTACCCGAATACCC"],  "cnifl_18899_1", 915],
    [["TGTTCACCATCGTGGCTTCA",	"TGGTGGCAAAAGCCCTGTTA"],  "cnifl_18899_2", 926],
    [["TCGCGAGTTTGTTGTTGTGC",	"CTCTCGTCGATGCTCTTGCT"],  "cnifl_19581_1", 406],
    [["GCCATTCCTGCTTTGTGTCA",	"AATACGAGCCAACCCTCAGC"],  "cnifl_27171_1", 851]
    ]
    )

individual_pcrs = Channel.from([
    "A_gdna_cnifl_4936_1",  "A_gdna_cnifl_6313_1",  "A_gdna_cnifl_6313_2",  "A_gdna_cnifl_11648_1",
    "A_gdna_cnifl_13073_1", "A_gdna_cnifl_15897_1",	"A_gdna_cnifl_15899_1",	"A_gdna_cnifl_18899_1",
    "A_gdna_cnifl_18899_2", "A_gdna_cnifl_19581_1",	"A_gdna_cnifl_27171_1",	"A_gdna_cnifl_15897_2",
    "A_gdna_cnifl_6313_3",  "B_gdna_cnifl_4936_1",	"B_gdna_cnifl_6313_1",	"B_gdna_cnifl_6313_2",
    "B_gdna_cnifl_11648_1",	"B_gdna_cnifl_13073_1", "B_gdna_cnifl_15897_1",	"B_gdna_cnifl_15899_1",
    "B_gdna_cnifl_18899_1",	"B_gdna_cnifl_18899_2",	"B_gdna_cnifl_19581_1",	"B_gdna_cnifl_27171_1",
    "B_gdna_cnifl_15897_2",	"B_gdna_cnifl_6313_3",  "C_gdna_cnifl_4936_1",	"C_gdna_cnifl_6313_1",
    "C_gdna_cnifl_6313_2",  "C_gdna_cnifl_11648_1",	"C_gdna_cnifl_13073_1", "C_gdna_cnifl_15897_1",
    "C_gdna_cnifl_15899_1",	"C_gdna_cnifl_18899_1",	"C_gdna_cnifl_18899_2",	"C_gdna_cnifl_19581_1",
    "C_gdna_cnifl_27171_1", "C_gdna_cnifl_15897_2",	"C_gdna_cnifl_6313_3",  "D_gdna_cnifl_4936_1",
    "D_gdna_cnifl_6313_1",	"D_gdna_cnifl_6313_2",  "D_gdna_cnifl_11648_1",	"D_gdna_cnifl_13073_1",
    "D_gdna_cnifl_15897_1",	"D_gdna_cnifl_15899_1",	"D_gdna_cnifl_18899_1",	"D_gdna_cnifl_18899_2",
    "D_gdna_cnifl_19581_1",	"D_gdna_cnifl_27171_1", "D_gdna_cnifl_15897_2",	"D_gdna_cnifl_6313_3",
    "A_cdna_cnifl_6313_2",  "A_cdna_cnifl_15897_2", "A_cdna_cnifl_6313_3",  "B_cdna_cnifl_18899_1",
    "B_cdna_cnifl_15897_2", "B_cdna_cnifl_6313_3",  "C_cdna_cnifl_4936_1",  "C_cdna_cnifl_6313_2",
    "C_cdna_cnifl_15897_1", "C_cdna_cnifl_18899_1", "C_cdna_cnifl_15897_2", "C_cdna_cnifl_6313_3",
    "D_cdna_cnifl_4936_1",  "D_cdna_cnifl_6313_2",  "D_cdna_cnifl_18899_1", "D_cdna_cnifl_15897_2",
    "D_cdna_cnifl_6313_3"
    ]).map{[it.split("_")[0] + "_" + it.split("_")[1], it.split("_")[2] + "_" + it.split("_")[3] + "_" + it.split("_")[4]]}

// A map to let us identify which scaffold the amplicon we were aiming for is from
target_scaffold_map = Channel.from([
    ["cnifl_4936_1" , "scaffold273"],
    ["cnifl_6313_1" , "scaffold39"],
    ["cnifl_6313_2" , "scaffold39"],
    ["cnifl_6313_3" , "scaffold39"],
    ["cnifl_11648_1", "scaffold30"],
    ["cnifl_13073_1", "scaffold1"],
    ["cnifl_15897_1", "scaffold337"],
    ["cnifl_15897_2", "scaffold337"],
    ["cnifl_15899_1", "scaffold337"],
    ["cnifl_18899_1", "scaffold493"],
    ["cnifl_18899_2", "scaffold493"],
    ["cnifl_19581_1", "scaffold416"],
    ["cnifl_27171_1", "scaffold279"]
])

// A map to let us identify which scaffold the amplicon we were aiming for is from
// and I have also included the bp range of where the reference CDS maps. I did this by opening
// up the reference mapped CDS file in py sam and looking at the reference_start and reference_ed
target_scaffold_map_full_name = Channel.from([
    [ "cnifl_4936_1_cdna" , "LJWW01000272.1 Exaiptasia pallida isolate CC7 scaffold273, whole genome shotgun sequence:128264-130134"],
    [ "cnifl_6313_1_cdna" , "LJWW01000039.1 Exaiptasia pallida isolate CC7 scaffold39, whole genome shotgun sequence:327291-329484"],
    [ "cnifl_6313_2_cdna" , "LJWW01000039.1 Exaiptasia pallida isolate CC7 scaffold39, whole genome shotgun sequence:327291-329484"],
    [ "cnifl_6313_3_cdna" , "LJWW01000039.1 Exaiptasia pallida isolate CC7 scaffold39, whole genome shotgun sequence:327291-329484"],
    ["cnifl_11648_1_cdna", "LJWW01000030.1 Exaiptasia pallida isolate CC7 scaffold30, whole genome shotgun sequence:61078-62941"],
    ["cnifl_13073_1_cdna", "LJWW01000001.1 Exaiptasia pallida isolate CC7 scaffold1, whole genome shotgun sequence:182325-184457"],
    ["cnifl_15897_1_cdna", "LJWW01000336.1 Exaiptasia pallida isolate CC7 scaffold337, whole genome shotgun sequence:35544-37588"],
    ["cnifl_15897_2_cdna", "LJWW01000336.1 Exaiptasia pallida isolate CC7 scaffold337, whole genome shotgun sequence:35544-37588"],
    ["cnifl_15899_1_cdna", "LJWW01000336.1 Exaiptasia pallida isolate CC7 scaffold337, whole genome shotgun sequence:213282-214977"],
    ["cnifl_18899_1_cdna", "LJWW01000492.1 Exaiptasia pallida isolate CC7 scaffold493, whole genome shotgun sequence:54787-56966"],
    ["cnifl_18899_2_cdna", "LJWW01000492.1 Exaiptasia pallida isolate CC7 scaffold493, whole genome shotgun sequence:54787-56966"],
    ["cnifl_19581_1_cdna", "LJWW01000415.1 Exaiptasia pallida isolate CC7 scaffold416, whole genome shotgun sequence:151158-153256"],
    ["cnifl_27171_1_cdna", "LJWW01000278.1 Exaiptasia pallida isolate CC7 scaffold279, whole genome shotgun sequence:166092-168459"],

    [ "cnifl_4936_1_gdna" , "LJWW01000272.1 Exaiptasia pallida isolate CC7 scaffold273, whole genome shotgun sequence:128264-130134"],
    [ "cnifl_6313_1_gdna" , "LJWW01000039.1 Exaiptasia pallida isolate CC7 scaffold39, whole genome shotgun sequence:327291-329484"],
    [ "cnifl_6313_2_gdna" , "LJWW01000039.1 Exaiptasia pallida isolate CC7 scaffold39, whole genome shotgun sequence:327291-329484"],
    [ "cnifl_6313_3_gdna" , "LJWW01000039.1 Exaiptasia pallida isolate CC7 scaffold39, whole genome shotgun sequence:327291-329484"],
    ["cnifl_11648_1_gdna", "LJWW01000030.1 Exaiptasia pallida isolate CC7 scaffold30, whole genome shotgun sequence:61078-62941"],
    ["cnifl_13073_1_gdna", "LJWW01000001.1 Exaiptasia pallida isolate CC7 scaffold1, whole genome shotgun sequence:182325-184457"],
    ["cnifl_15897_1_gdna", "LJWW01000336.1 Exaiptasia pallida isolate CC7 scaffold337, whole genome shotgun sequence:35544-37588"],
    ["cnifl_15897_2_gdna", "LJWW01000336.1 Exaiptasia pallida isolate CC7 scaffold337, whole genome shotgun sequence:35544-37588"],
    ["cnifl_15899_1_gdna", "LJWW01000336.1 Exaiptasia pallida isolate CC7 scaffold337, whole genome shotgun sequence:213282-214977"],
    ["cnifl_18899_1_gdna", "LJWW01000492.1 Exaiptasia pallida isolate CC7 scaffold493, whole genome shotgun sequence:54787-56966"],
    ["cnifl_18899_2_gdna", "LJWW01000492.1 Exaiptasia pallida isolate CC7 scaffold493, whole genome shotgun sequence:54787-56966"],
    ["cnifl_19581_1_gdna", "LJWW01000415.1 Exaiptasia pallida isolate CC7 scaffold416, whole genome shotgun sequence:151158-153256"],
    ["cnifl_27171_1_gdna", "LJWW01000278.1 Exaiptasia pallida isolate CC7 scaffold279, whole genome shotgun sequence:166092-168459"]
])

// NB we did the initial QC of the short reads

// get the unique and name pairs to run the PCRs on.

// Turns out that uniqueing the long reads is almost pointless because of the higher error rate so I won't bother.
// Plus it is breaking the pipeline because sometimes names files aren't produced.
process mothur_fastq_to_fasta{
    tag "${sample_origin}"
    container "didillysquat/mothur:latest"
    cpus 1

    input:
    tuple val(sample_origin), path(fastq) from seq_ch

    output:
    tuple val(sample_origin), path("${sample_origin}.fasta") into unique_fasta_ch

    script:
    """
    gzip -d -c $fastq > ${sample_origin}.fastq
    echo 'fastq.info(fastq=${sample_origin}.fastq)' >> batch_file.txt
    mothur batch_file.txt
    """
}

// Here we run a each PCR on each 
process mothur_in_silico_pcr_for_real{
    tag "${sample_origin}_${pcr_name}"
    container "didillysquat/mothur:latest"
    publishDir "/home/humebc/projects/cnifl/results/in_silico_pcr/${sample_origin}_${pcr_name}", mode: 'copy'
    cpus 50

    input:
    tuple val(pcr_name), val(sample_origin), path(unique_fasta), val(primers), val(screen_len) from individual_pcrs.combine(unique_fasta_ch, by:0).combine(primers, by:1)

    output:
    tuple val(pcr_name), val(sample_origin), path("${sample_origin}_${pcr_name}.pcr.minsize.fasta") into med_in_ch

    script:
    """
    echo forward ${primers[0]} >> oligo.txt
    echo reverse ${primers[1]} >> oligo.txt
    echo 'pcr.seqs(fasta=${unique_fasta}, oligos=oligo.txt, pdiffs=3, rdiffs=3, processors=${task.cpus})' >> batch_file.txt
    echo 'screen.seqs(fasta=${sample_origin}.pcr.fasta, minlength=$screen_len)' >> batch_file.txt
    mothur batch_file.txt
    rm batch_file.txt
    mv ${sample_origin}.pcr.good.fasta ${sample_origin}_${pcr_name}.pcr.minsize.fasta
    rm *scrap*
    """
}

// create val of sample, origin, pcr, fasta and modify the header line to add the sample name and a count
// med_in_ch.map{[it[1].split("_")[0], it[1].split("_")[1], it[0], it[2]]}
process modify_header{
    tag "${sample_origin}_${pcr_name}"
    cpus 1

    input:
    tuple val(pcr_name), val(sample_origin), path(fasta) from med_in_ch

    output:
    tuple val(pcr_name), val(sample_origin), path("${sample_origin}_${pcr_name}.pcr.minsize.modifiedHeader.fasta") into modify_header_out_ch

    script:
    sample = sample_origin.split("_")[0]
    origin = sample_origin.split("_")[1]
    """
    awk 'BEGIN{count=0} {if(\$0 ~ />/){print ">${sample}_${origin}_${pcr_name}_" count; count = count + 1} else{print \$0}}' $fasta > ${sample_origin}_${pcr_name}.pcr.minsize.modifiedHeader.fasta
    """
}

modify_header_out_ch.map{[it[0] + "_" + it[1].split("_")[1], it[2]]}.groupTuple().into{map_out_multiple_products_ch_in; mpac_bio_bam}
// modify_header_out_ch.map{[it[0] + "_" + it[1].split("_")[1], it[2]]}.groupTuple().view()

//group by pcr_origin so that we have all the fastqs for a given pcr_origin
// At this point we need to map out the different products that are in the PCRs
process mapPacBio{
    tag "${pcr_origin}"
    container "staphb/bbtools:latest"
    publishDir "/home/humebc/projects/cnifl/results/multi_sample_multi_product_sams", mode: 'symlink'
    cpus 50

    input:
    tuple val(pcr_origin), path(fatas) from map_out_multiple_products_ch_in
    file ref

    output:
    tuple val(pcr_origin), path("${pcr_origin}_multiSample_fasta.mapped.no.header.sam") into med_out_ch

    script:
    """
    cat *.fasta > ${pcr_origin}_multiSample_fasta.fasta
    mapPacBio.sh in=${pcr_origin}_multiSample_fasta.fasta out=${pcr_origin}_multiSample_fasta.mapped.sam build=1 threads=${task.cpus} nodisk
    samtools view ${pcr_origin}_multiSample_fasta.mapped.sam > ${pcr_origin}_multiSample_fasta.mapped.no.header.sam
    rm ${pcr_origin}_multiSample_fasta.mapped.sam
    """
}

// Here we produce multiple fastas according to which scaffold the products mapped to
process split_by_mapping{
    tag "${pcr_origin}"
    conda "python=3 numpy pandas biopython"
    publishDir "/home/humebc/projects/cnifl/results/multi_sample_single_product_fastas", mode: 'symlink'

    input:
    tuple val(pcr_origin), path(sam) from med_out_ch
    file split_by_sam

    output:
    tuple val(pcr_origin), path("*.fasta") into split_alignments_out_ch

    script:
    """
    python3 $split_by_sam $sam
    """
}

// Some of the alignments have a large ammount of sequencs in them (i.e. 300 000) and this is causing alignment to be extremely slow, even with mafft
// We don't need that many seuqneces in order to be able to accurately cout the haplotype so we will subsample this down to 5000 sequeneces per sample per alignment
// and work with that. If there are less than 5000 samples then we will just take what we have
// We'll do this in python rather than awk because it will be quicker for me to write
// From this point on we are also only working with the target scaffold amplicon and we are not continuing with the rest
process subsample{
    tag "${pcr}_${origin}"
    conda "python=3"

    input:
    tuple val(pcr), val(origin), path(fastas), val(target_scaffold) from split_alignments_out_ch.map{[it[0].replaceAll("_" + it[0].split("_").last(), ""), it[0].split("_").last(), it[1]]}.combine(target_scaffold_map, by:0)
    file subsample

    output:
    tuple val(pcr), val(origin), val(target_scaffold), path("${pcr}_${origin}_multiSample_fasta.mapped.${target_scaffold}.subsampled.fasta") into subsample_out_ch

    script:
    """
    python3 $subsample ${pcr}_${origin}_multiSample_fasta.mapped.${target_scaffold}.fasta
    """
}

// Then need to align the mapped out fastas in preparation for analysis of the entropy
// We will filter to only work with those 'target' products according to the target_scaffold_map
// example of an input fasta cnifl_15899_1_gdna_multiSample_fasta.mapped.scaffold94.fasta
process align_mapped_out{
    tag "${pcr}_${origin}"
    conda "mafft"
    publishDir "/home/humebc/projects/cnifl/results/aligned_multi_sample_single_product_fastas", mode: 'symlink'

    input:
    tuple val(pcr), val(origin), val(target_scaffold), path(subbed_fasta) from subsample_out_ch

    output:
    tuple val(pcr), val(origin), path("${pcr}_${origin}_multiSample_fasta.mapped.${target_scaffold}.subsampled.aligned.fasta") into align_mapped_out_ch, long_aligned_out

    script:
    """
    mafft ${pcr}_${origin}_multiSample_fasta.mapped.${target_scaffold}.subsampled.fasta > ${pcr}_${origin}_multiSample_fasta.mapped.${target_scaffold}.subsampled.aligned.fasta
    """
}

raw_seqs = Channel.fromFilePairs("/home/humebc/projects/cnifl/raw_seq_data/short_read_data/*{R1,R2}*.fastq.gz")

// Clean adapters from the short read input seqs
process fastp{
    tag "${sample}"
    conda "fastp"

    input:
    tuple val(sample), path(reads) from raw_seqs

    output:
    file "*.html" into fastp_out_no_dedup_ch
    tuple val(sample), file("${sample}_1.clean.fq.gz"), file("${sample}_2.clean.fq.gz") into fastp_out_ch

    script:
    """
    fastp -q 20 -i ${reads[0]} -I ${reads[1]} -o ${sample}_1.clean.fq.gz -O ${sample}_2.clean.fq.gz
    """
}

// Map the short reads to the aiptasia genome
process bbmap{
    tag "${sample}"
    container "staphb/bbtools:latest"
    publishDir "/home/humebc/projects/cnifl/results/short_read_mapped", mode: 'symlink'
    cpus 50

    input:
    tuple val(sample), file(fwd), file(rev) from fastp_out_ch
    file ref

    output:
    tuple val(sample), path("${sample}.mapped.sorted.bam"), path("${sample}.mapped.sorted.bam.bai") into bbmap_out

    script:
    """
    bbmap.sh in1=$fwd in2=$rev outm=${sample}.mapped.bam build=1 threads=${task.cpus}
    samtools view -h ${sample}.mapped.bam | awk -F'\\t' 'substr(\$0,1,1)=="@" || (\$9>= 1 && \$9<=500) || (\$9<=-1 && \$9>=-500)' | samtools view -b -f2 > ${sample}.mapped.filtered.bam
    samtools sort ${sample}.mapped.filtered.bam > ${sample}.mapped.sorted.bam
    samtools index ${sample}.mapped.sorted.bam
    rm ${sample}.mapped.bam ${sample}.mapped.filtered.bam
    """
}

// We want to investigate whether the short reads support the sequence diversity we are seeing in the long reads.
// To do this, we will aim to do two sets of intersects using BED tools so that we only end up with the short reads
// that cover the regions of the long reads that we are interested in.
// The first intersect will be the long reads with the mapped reference CDSs.
// Once we have done this first set of interects I want to check that we only have one amplicon
// i.e. 1 loci per long read PCR, else something unexpected is happening. ( I checked and we have multiple amplicons so we have to filter as done below)
// From here, we can do an intersect
// of the short reads with the reference-intersected long reads. We have a large number
// of short read sequences, that we should split up into cDNA and gDNA, but other than that
// we don't know which cnifls or samples the short reads refer to so we will do a pairwise
// intersect where we map all short read reads to all long read reads.
// This will end up with a whole load of short reads bams per long read bam so we will then need to 
// merge the short read bams. In this way we should end up with one long read and possibly an associated short read.
// We can then look at these in a fasta to see if the diversity is supported or something similar.

// Remap the long reads this time into bam format with headers and sort and index
process mapPacBio_long_bam{
    tag "${pcr_origin}"
    container "staphb/bbtools:latest"
    publishDir "/home/humebc/projects/cnifl/results/long_reads_mapped_to_ref", mode: 'symlink'
    cpus 50

    input:
    tuple val(pcr_origin), path(fatas) from mpac_bio_bam
    file ref

    output:
    tuple val(pcr_origin), path("${pcr_origin}_multiSample_fasta.mapped.sorted.bam"), path("${pcr_origin}_multiSample_fasta.mapped.sorted.bam.bai") into intersect_long_ref_ch

    script:
    """
    cat *.fasta > ${pcr_origin}_multiSample_fasta.fasta
    mapPacBio.sh in=${pcr_origin}_multiSample_fasta.fasta out=${pcr_origin}_multiSample_fasta.mapped.bam build=1 threads=${task.cpus} nodisk
    samtools sort ${pcr_origin}_multiSample_fasta.mapped.bam > ${pcr_origin}_multiSample_fasta.mapped.sorted.bam
    samtools index ${pcr_origin}_multiSample_fasta.mapped.sorted.bam
    rm ${pcr_origin}_multiSample_fasta.mapped.bam
    rm ${pcr_origin}_multiSample_fasta.fasta
    """
}


// Intersect the long read mapped reads with the reference cds reads that have also been mapped
process intersect_long_ref{
    tag "${pcr_origin}"
    container "didillysquat/bedtools_samtools"
    publishDir "/home/humebc/projects/cnifl/results/long_intersect_ref", mode: 'symlink'
    cpus 50

    input:
    tuple val(pcr_origin), path(long_bam), path(long_bai) from intersect_long_ref_ch
    file ref_bam
    file ref_bai

    output:
    tuple val(pcr_origin), path("${pcr_origin}.longMappedToRef.sorted.bam"), path("${pcr_origin}.longMappedToRef.sorted.bam.bai") into intersect_long_ref_out_ch

    script:
    """
    bedtools intersect -a $long_bam -b $ref_bam  > ${pcr_origin}.longMappedToRef.bam
    samtools sort ${pcr_origin}.longMappedToRef.bam > ${pcr_origin}.longMappedToRef.sorted.bam
    samtools index ${pcr_origin}.longMappedToRef.sorted.bam
    """
}

// Here we should check that there is only one amplicon in each of the bam files. Perhaps we can check this in python with pysam.
// We have looked and this assumption holds true for some of the samples but not for others.
// So we will once again have to extract out the scaffold and therefore amplicons that we are interested/targeted in each of the bams
// before these bams will be a suitable template for running the intersect with the short read seqs.
process filter_by_target_scaff{
    container "didillysquat/bedtools_samtools"
    tag "${pcr_origin}"
    publishDir "/home/humebc/projects/cnifl/results/long_intersect_ref_filtered", mode: 'symlink'
    cpus 50

    input:
    tuple val(pcr_origin), path(bam), path(bai), val(target_coords) from intersect_long_ref_out_ch.combine(target_scaffold_map_full_name, by: 0)

    output:
    tuple val(pcr_origin), path("${pcr_origin}.longMappedToRef.filteredToTarget.sorted.bam"), path("${pcr_origin}.longMappedToRef.filteredToTarget.sorted.bam.bai") into filter_by_target_scaff_out_ch, merge_with_short_bam_in_ch, merge_with_short_bam_subsample_in_ch

    script:
    """
    samtools view -hb ${pcr_origin}.longMappedToRef.sorted.bam \"$target_coords\" > ${pcr_origin}.longMappedToRef.filteredToTarget.bam
    samtools sort ${pcr_origin}.longMappedToRef.filteredToTarget.bam > ${pcr_origin}.longMappedToRef.filteredToTarget.sorted.bam
    samtools index ${pcr_origin}.longMappedToRef.filteredToTarget.sorted.bam
    rm ${pcr_origin}.longMappedToRef.filteredToTarget.bam
    """
}

// Finally we want to intersect the short read sequences with these filtered long read bams.
// We will interesect each short read bam file with each long read bam file.
// Then we will need to merge all of the short read bams that are associated to a given long read bam.
// These files are then what we will want to analyses by simply looking at them in the alignment viewer.

// split the cdna short read samples into cdna and gdna channels
bbmap_out.map{
    if (it[0].contains("gDNA")){
        ["gdna", it[0], it[1], it[2]]
    }else if (it[0].contains("cDNA")){
        ["cdna", it[0], it[1], it[2]]
    }
}.branch {
        short_cdna: it[0] == "cdna"
        short_gdna: it[0] == "gdna"
    }.set { short_result }


// Do the same for the long reads
filter_by_target_scaff_out_ch.map{
    if (it[0].contains("gdna")){
        ["gdna", it[0], it[1], it[2]]
    }else if (it[0].contains("cdna")){
        ["cdna", it[0], it[1], it[2]]
    }
}.branch{
    long_cdna: it[0] == "cdna"
    long_gdna: it[0] == "gdna"
}.set {long_result}

// Intersect the short reads with the long reads.
process intersect_short_long{
    tag "${pcr_origin}_${sample}"
    container "didillysquat/bedtools_samtools"
    publishDir "/home/humebc/projects/cnifl/results/short_intersect_long_individual", mode: 'symlink'
    cpus 1

    input:
    tuple val(origin), val(pcr_origin), path(filtered_long_bam), path(filtered_long_bai), val(origin), val(sample), path(short_bam), path(short_bai) from long_result.long_cdna.combine(short_result.short_cdna).mix(long_result.long_gdna.combine(short_result.short_gdna))

    output:
    tuple val(pcr_origin), path("${pcr_origin}.${sample}.short_intersect_long.sorted.bam"), path("${pcr_origin}.${sample}.short_intersect_long.sorted.bam.bai") into intersect_short_long_out_ch

    script:
    """
    bedtools intersect -a $short_bam -b $filtered_long_bam  > ${pcr_origin}.${sample}.short_intersect_long.bam
    samtools sort ${pcr_origin}.${sample}.short_intersect_long.bam > ${pcr_origin}.${sample}.short_intersect_long.sorted.bam
    samtools index ${pcr_origin}.${sample}.short_intersect_long.sorted.bam
    rm ${pcr_origin}.${sample}.short_intersect_long.bam
    """
}

// Finally we want to merge all of the individual short intersect long bams by pcr_origin
// There will be some empty bams so we will screen for these at the input by requiring > 100kb in size for the bam.
process merge_short_intersect_long{
    tag "${pcr_origin}"
    container "didillysquat/bedtools_samtools"
    publishDir "/home/humebc/projects/cnifl/results/short_intersect_long_merged", mode: 'symlink'
    cpus 1

    input:
    tuple val(pcr_origin), path(bams), path(bais) from intersect_short_long_out_ch.filter{it[1].size() > 100000}.groupTuple()

    output:
    tuple val(pcr_origin), path("${pcr_origin}.short_intersect_long.merged.sorted.bam"), path("${pcr_origin}.short_intersect_long.merged.sorted.bam.bai") into merge_short_intersect_long_out_ch, merge_short_intersect_long_subsample_out_ch

    script:
    """
    samtools merge -o ${pcr_origin}.short_intersect_long.merged.bam *.bam
    samtools sort ${pcr_origin}.short_intersect_long.merged.bam > ${pcr_origin}.short_intersect_long.merged.sorted.bam
    samtools index ${pcr_origin}.short_intersect_long.merged.sorted.bam
    rm ${pcr_origin}.short_intersect_long.merged.bam
    """
}


// rather than work with fasta
// lets try to stay in the bam format so that we are relative to the reference.
// let's merge the short reads with the corresponding long reads and look at them in IGV.
process merge_long_and_short_bams{
    tag "${pcr_origin}"
    container "didillysquat/bedtools_samtools"
    publishDir "/home/humebc/projects/cnifl/results/short_and_long_merged_by_pcr_origin_bams", mode: 'symlink'
    cpus 1

    input:
    tuple val(pcr_origin), path(short_bam), path(short_bai), path(long_bam), path(long_bai) from merge_short_intersect_long_out_ch.combine(merge_with_short_bam_in_ch, by:0)

    output:
    tuple val(pcr_origin), path("${pcr_origin}.short_and_long_reads.sorted.bam"), path("${pcr_origin}.short_and_long_reads.sorted.bam.bai") into merge_long_and_short_bams_ch

    script:
    """
    samtools merge -o ${pcr_origin}.short_and_long_reads.bam $short_bam $long_bam
    samtools sort ${pcr_origin}.short_and_long_reads.bam > ${pcr_origin}.short_and_long_reads.sorted.bam
    samtools index ${pcr_origin}.short_and_long_reads.sorted.bam
    rm ${pcr_origin}.short_and_long_reads.bam
    """
}

// The problem is that we want to be able to visualize some long and some short in the same alignment.
// to acheive that I want to subsample both of the bams before merging them.
process merge_long_and_short_bams_subsample{
    tag "${pcr_origin}"
    container "didillysquat/bedtools_samtools"
    publishDir "/home/humebc/projects/cnifl/results/subsampled_short_and_long_merged_by_pcr_origin_bams", mode: 'symlink'
    cpus 1

    input:
    tuple val(pcr_origin), path(short_bam), path(short_bai), path(long_bam), path(long_bai) from merge_short_intersect_long_subsample_out_ch.combine(merge_with_short_bam_subsample_in_ch, by:0)

    output:
    tuple val(pcr_origin), path("${pcr_origin}.short_and_long_reads.subsampled.sorted.bam"), path("${pcr_origin}.short_and_long_reads.subsampled.sorted.bam.bai") into merge_long_and_short_bams_subsample_ch

    script:
    """
    # Make the 100 subsampled shorts
    samtools view -H $short_bam > short.headers.sam
    samtools view $short_bam | shuf | head -n 400 > short.reads.sam
    cat short.headers.sam short.reads.sam > short.subsampled.sam
    samtools view -bh short.subsampled.sam > short.subsampled.bam
    samtools sort short.subsampled.bam > short.subsampled.sorted.bam
    samtools index short.subsampled.sorted.bam

    # Make the 100 subsampled longs
    samtools view -H $long_bam > long.headers.sam
    samtools view $long_bam | shuf | head -n 100 > long.reads.sam
    cat long.headers.sam long.reads.sam > long.subsampled.sam
    samtools view -bh long.subsampled.sam > long.subsampled.bam
    samtools sort long.subsampled.bam > long.subsampled.sorted.bam
    samtools index long.subsampled.sorted.bam

    # Finally do the merge of the short and long subsampled reads
    samtools merge -o ${pcr_origin}.short_and_long_reads.subsampled.bam short.subsampled.sorted.bam long.subsampled.sorted.bam
    samtools sort ${pcr_origin}.short_and_long_reads.subsampled.bam > ${pcr_origin}.short_and_long_reads.subsampled.sorted.bam
    samtools index ${pcr_origin}.short_and_long_reads.subsampled.sorted.bam
    """
}
