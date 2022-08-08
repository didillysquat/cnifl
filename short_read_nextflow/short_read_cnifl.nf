#!/usr/bin/env nextflow

// Chris has magically produced a load of short read data that may overlap with some of the amplicons that Matt has sequenced
// with the long read sequencing.

// There are a load of different samples worth of sequencing and we're not sure how they relate to different amplicons etc.

// So the easiest way forwards that is robust is to map all of the short reads to the aiptasia genome
// and then look at the regions that have been mapped to that overlap with the mapping positions of the 
// long read sequences.

raw_seqs = Channel.fromFilePairs("/home/humebc/projects/cnifl/raw_seq_data/short_read_data/*{R1,R2}*.fastq.gz")
ref = file("/home/humebc/projects/cnifl/reference/ref")

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

// Now we want to map to the aiptasia genome. We can do this using bbmap so that it is compatible with the long
// read mapping
process bbmap{
    tag "${sample}"
    container "staphb/bbtools:latest"
    publishDir "/home/humebc/projects/cnifl/results/short_read_mapped", mode: 'symlink'
    cpus 50

    input:
    tuple val(sample), file(fwd), file(rev) from fastp_out_ch
    file ref

    output:
    tuple path("${sample}.mapped.sorted.bam"), path("${sample}.mapped.sorted.bam.bai") into bbmap_out

    script:
    """
    bbmap.sh in1=$fwd in2=$rev outm=${sample}.mapped.bam build=1 threads=${task.cpus}
    samtools sort ${sample}.mapped.bam > ${sample}.mapped.sorted.bam
    samtools index ${sample}.mapped.sorted.bam
    rm ${sample}.mapped.bam
    """
}

// bbmap_out.collect().view()

// process merge_bams{
//     tag "${sample}"
//     container "staphb/bbtools:latest"


//     input:
//     tuple path(bams) from bbmap_out.collect()


//     output:
//     tuple path("shortreads.sorted.bam"), path("shortreads.sorted.bam.bai") into bbmap_out

//     script:
//     """
//     samtools merge -o shortreads.bam *.bam
//     """
// }

// Need to get the reads that share an interval with the long read mappings.