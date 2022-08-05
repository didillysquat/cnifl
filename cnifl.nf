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
split_by_sam = file("/home/humebc/projects/cnifl/testing_med/sam.py")
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

// get the unique and name pairs to run the PCRs on.
// NB ignore the name :)
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

modify_header_out_ch.map{[it[0] + "_" + it[1].split("_")[1], it[2]]}.groupTuple().set{map_out_multiple_products_ch_in}
// modify_header_out_ch.map{[it[0] + "_" + it[1].split("_")[1], it[2]]}.groupTuple().view()

//group by pcr_origin so that we have all the fastqs for a given pcr_origin
// At this point we need to map out the different products that are in the PCRs
process mapPacBio{
    tag "${pcr_origin}"
    container "staphb/bbtools:latest"
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

// Then need to align the mapped out fastas in preparation for analysis of the entropy
// We will filter to only work with those 'target' products according to the target_scaffold_map
// example of an input fasta cnifl_15899_1_gdna_multiSample_fasta.mapped.scaffold94.fasta
process align_mapped_out{
    tag "${pcr}_${origin}"
    conda "mafft"
    publishDir "/home/humebc/projects/cnifl/results/aligned_multi_sample_single_product_fastas", mode: 'symlink'

    input:
    tuple val(pcr), val(origin), path(fastas), val(target_scaffold) from split_alignments_out_ch.map{[it[0].replaceAll("_" + it[0].split("_").last(), ""), it[0].split("_").last(), it[1]]}.combine(target_scaffold_map, by:0)

    output:
    tuple val(pcr), val(origin), path("${pcr}_${origin}_multiSample_fasta.mapped.${target_scaffold}.aligned.fasta") into align_mapped_out_ch

    script:
    """
    mafft ${pcr}_${origin}_multiSample_fasta.mapped.${target_scaffold}.fasta > ${pcr}_${origin}_multiSample_fasta.mapped.${target_scaffold}.aligned.fasta
    """
}

// At this point we are ready to look at the entropy
// It is perhaps good if we put this into a compound figure.

// seq_ch = Channel.fromPath('/home/humebc/projects/cnifl/raw_seq_data/concatenated_files/*.fastq.gz').map{[it.getName().split("_")[0], it.getName().split("_")[1].replaceAll(".fastq.gz", ""), it]}
// aipgene11648_ref = file("/home/humebc/projects/cnifl/reference/aipgene11648.fa")
// aipgene13073_ref = file("/home/humebc/projects/cnifl/reference/aipgene13073.fa")
// aipgene13580_ref = file("/home/humebc/projects/cnifl/reference/aipgene13580.fa")
// aipgene15890_ref = file("/home/humebc/projects/cnifl/reference/aipgene15890.fa")
// aipgene15897_ref = file("/home/humebc/projects/cnifl/reference/aipgene15897.fa")
// aipgene15899_ref = file("/home/humebc/projects/cnifl/reference/aipgene15899.fa")
// aipgene18899_ref = file("/home/humebc/projects/cnifl/reference/aipgene18899.fa")
// aipgene19581_ref = file("/home/humebc/projects/cnifl/reference/aipgene19581.fa")
// aipgene27171_ref = file("/home/humebc/projects/cnifl/reference/aipgene27171.fa")
// aipgene27771_ref = file("/home/humebc/projects/cnifl/reference/aipgene27771.fa")
// aipgene4936_ref = file("/home/humebc/projects/cnifl/reference/aipgene4936.fa")
// aipgene6313_ref = file("/home/humebc/projects/cnifl/reference/aipgene6313.fa")
// aipgene7322_ref = file("/home/humebc/projects/cnifl/reference/aipgene7322.fa")

// ref_fa = file("/home/humebc/projects/cnifl/reference/mmseqs/ref.fa")

// split_script = file("/home/humebc/projects/cnifl/bin/split_fasta_by_gene.py")



// // porechop
// process porechomp{
//     tag "${sample} ${origin}"
//     container "biocontainers/porechop:v0.2.4dfsg-1-deb_cv1"
//     cpus 20

//     input:
//     tuple val(sample), val(origin), path(fastq) from seq_ch

//     output:
//     tuple val(sample), val(origin), path("${sample}_${origin}.chopped.fastq.gz") into pore_chop_out_ch

//     script:
//     """
//     porechop -i $fastq -o ${sample}_${origin}.chopped.fastq.gz --threads ${task.cpus}
//     """
// }

// // Once we have the reads clean we will want to map them against the reference 13 sequences
// // We will try doing it with bbsplit 
// // process bbsplit{
// //     tag "${sample} ${origin}"
// //     container "staphb/bbtools:latest"
// //     cpus 10

// //     input:
// //     tuple val(sample), val(origin), path(fastq) from pore_chop_out_ch
// //     file aipgene11648_ref
// //     file aipgene13073_ref
// //     file aipgene13580_ref
// //     file aipgene15890_ref
// //     file aipgene15897_ref
// //     file aipgene15899_ref
// //     file aipgene18899_ref
// //     file aipgene19581_ref
// //     file aipgene27171_ref
// //     file aipgene27771_ref
// //     file aipgene4936_ref 
// //     file aipgene6313_ref 
// //     file aipgene7322_ref 

// //     output:
// //     path("${sample}_${origin}_*") into bbsplit_out_ch

// //     script:
// //     """
// //     bbsplit.sh -Xmx800g build=1 threads=${task.cpus} in1=$fastq \
// //     ref_aipgene11648=$aipgene11648_ref \
// //     ref_aipgene13073=$aipgene13073_ref \
// //     ref_aipgene13580=$aipgene13580_ref \
// //     ref_aipgene15890=$aipgene15890_ref \
// //     ref_aipgene15897=$aipgene15897_ref \
// //     ref_aipgene15899=$aipgene15899_ref \
// //     ref_aipgene18899=$aipgene18899_ref \
// //     ref_aipgene19581=$aipgene19581_ref \
// //     ref_aipgene27171=$aipgene27171_ref \
// //     ref_aipgene27771=$aipgene27771_ref \
// //     ref_aipgene4936=$aipgene4936_ref \
// //     ref_aipgene6313=$aipgene6313_ref \
// //     ref_aipgene7322=$aipgene7322_ref \
// //     out_aipgene11648=${sample}_${origin}_aipgene11648.fq.gz \
// //     out_aipgene13073=${sample}_${origin}_aipgene13073.fq.gz \
// //     out_aipgene13580=${sample}_${origin}_aipgene13580.fq.gz \
// //     out_aipgene15890=${sample}_${origin}_aipgene15890.fq.gz \
// //     out_aipgene15897=${sample}_${origin}_aipgene15897.fq.gz \
// //     out_aipgene15899=${sample}_${origin}_aipgene15899.fq.gz \
// //     out_aipgene18899=${sample}_${origin}_aipgene18899.fq.gz \
// //     out_aipgene19581=${sample}_${origin}_aipgene19581.fq.gz \
// //     out_aipgene27171=${sample}_${origin}_aipgene27171.fq.gz \
// //     out_aipgene27771=${sample}_${origin}_aipgene27771.fq.gz \
// //     out_aipgene4936=${sample}_${origin}_aipgene4936.fq.gz \
// //     out_aipgene6313=${sample}_${origin}_aipgene6313.fq.gz \
// //     out_aipgene7322=${sample}_${origin}_aipgene7322.fq.gz \
// //     outu=${sample}.unmapped.fq.gz
// //     """
// // }

// process mmseqs{
//     tag "${sample} ${origin}"
//     container "soedinglab/mmseqs2:latest"
//     cpus 10

//     input:
//     tuple val(sample), val(origin), path(fastq) from pore_chop_out_ch
//     file ref_fa

//     output:
//     tuple val(sample), val(origin), path("${sample}_${origin}_alnResult.m8"), path(fastq) into mmseqs_out_ch

//     script:
//     """
//     mmseqs easy-search $fastq $ref_fa  ${sample}_${origin}_alnResult.m8 tmp --search-type 3
//     """
// }

// process split_fastq_using_mmseqs_matches{
//     tag "${sample} ${origin}"
//     conda "python3 conda-forge::biopython conda-forge::gzip"

//     input:
//     tuple val(sample), val(origin), path(mmseqs_alignment), path(fastq) from mmseqs_out_ch
//     file split_script

//     output:
//     path("*fa.gz") into split_fastq_out

//     script:
//     """
//     python3 $split_script $mmseqs_alignment $fastq $sample $origin
//     """
// }

// // Convert the fastq to fasta and modify the header line so that it includes the sample name and a count
// // process fastq_to_fasta{
// //     tag "${unique_name}"
// //     container 'biocontainers/seqtk:v1.3-1-deb_cv1'

// //     input:
// //     tuple val(sample), val(origin), val(gene), file(fastq) from bbsplit_out_ch.collect().map{[it.getName().replaceAll("fastq.gz", "").split("_")[0], it.getName().replaceAll("fastq.gz", "").split("_")[1], it.getName().replaceAll("fastq.gz", "").split("_")[2], it]}

// //     output:
// //     tuple val("${origin}_${gene}"), file("${sample}_${origin}_${gene}.sample_headers.fa.gz") into fastq_to_fasta_out_ch

// //     script:
// //     """
// //     seqtk seq -A $fastq | awk 'BEGIN{count=0} {if(\$0 ~ />/){print ">${sample}_" count; count = count + 1} else{print \$0}}' - | gzip > ${sample}_${origin}_${gene}.sample_headers.fa.gz
// //     """
// // }

// // // Then group the fasta files by gene origin combinations i.e. so that we have collections of 4 fasta files that are the four samples
// // // And create a single fasta from these
// // // Then use this as input to med.
// // process med{
// //     tag "${origin_gene}"
// //     container "meren/oligotyping:latest"
// //     publishDir "results/med/${origin_gene}", mode: 'symlink'

// //     input:
// //     tuple val(origin_gene), file(fastas) from fastq_to_fasta_out_ch.groupTuple()

// //     output:
// //     file("*") into med_out_ch

// //     script:
// //     """

// //     """


// // }
// // fastq_to_fasta_out_ch.groupTuple()
