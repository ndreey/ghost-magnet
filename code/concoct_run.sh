#!/bin/bash

#SBATCH --job-name=bowtie2_test          # name that will show up in the queue
#SBATCH --output=slurm-bowtie2_build_test-%j.out       # filename of the output; the %j is equal to jobID
#SBATCH --error=slurm-bowtie2_build_test-%j.err        #
#SBATCH --partition=cpuqueue        #
#SBATCH --ntasks=1                  # number of tasks (analyses) to run
#SBATCH --cpus-per-task=2           # the number of threads allocated to each task
#SBATCH --mem-per-cpu=4G           # memory per cpu-core
#SBATCH --time=01:30:00             # time for analysis (day-hour:min:sec)
#SBATCH --mail-type=ALL             # send all type of email
#SBATCH --mail-user=andre.bourbonnais@sund.ku.dk


# I. Define directory names [DO NOT CHANGE]
# =========================================
TMPDIR=/tmp/00_concoct_tmp/
submitdir=${pwd}
workdir=${TMPDIR}
path_gfr=/mnt/c/Users/andbo/thesis_andbo

echo "$(date)    Accessed ${workdir}"

# 1. Lock and load module and data
# ============================================

repli="x700c"


for i in {2..11}
do
    jobid=$i

    # Different host-contamination level
    hc_level=("00" "01" "02" "03" "04" "05" "06" "07" "08" "090" "095")
    
    # Define prefix based on array id
    hc_prefix=${hc_level[$jobid-1]}
    
    assembly=${repli}/${hc_prefix}_700c.contigs.fa.gz
    reads=/home/andbo/02_final_approach/bin_refs/${repli}/reads

    echo "$(date)    Replica: ${repli} Sample: ${hc_prefix} ENGAGED"
    echo "$(date)    Bowtie2 Engaged!"
    # # Builds index of assembly
    bowtie2-build --threads 4 ${assembly} ${hc_prefix}_index
    
    # Mapping reads
    bowtie2 -p 4 -q \
        -x ${hc_prefix}_index \
        -1 ${reads}/${hc_prefix}_filt_R1.fastq.gz \
        -2 ${reads}/${hc_prefix}_filt_R2.fastq.gz  \
        -S ${hc_prefix}_map.sam
    echo "$(date)    Bowtie2 Disengaged!"
##
    echo "$(date)    SAMtools Engaged!"    
    samtools view -@ 4 -b -S ${hc_prefix}_map.sam > ${hc_prefix}_map.bam
    echo "$(date)    SAMtools view complete" 
    samtools sort -@ 4 ${hc_prefix}_map.bam -o ${hc_prefix}_map_sorted.bam
    echo "$(date)    SAMtools sort complete" 
    samtools index -@ 4 ${hc_prefix}_map_sorted.bam ${hc_prefix}_map_sorted.bam.bai
    echo "$(date)    SAMtools index complete" 
    echo "$(date)    SAMtools Disengaged!"
##
    #
    echo "$(date)    concoct Engaged!"
    zcat ${assembly} > ${hc_prefix}_tmp.fa
    #
    #
    # Contigs are split into lengths of 10k while the .bed file specifies
    # which contig the sub contigs belong to.
    #
    # -c, --chunk_size    CHUNK_SIZE
    # -o, --overlap_size  OVERLAP_SIZE
    # -m, --merge_last    Concatenate final part to last contig
    # -b, --bedfile       BEDfile to be created with exact regions of the
    # original contigs corresponding to the newly created contigs
    cut_up_fasta.py ${hc_prefix}_tmp.fa -c 10000 -o 0 --merge_last \
        -b ${hc_prefix}_contigs_10K.bed > ${hc_prefix}_contigs_10K.fa
    echo "$(date)    Contigs are cut into subcontigs"
    
    # Generates the coverage table
    concoct_coverage_table.py ${hc_prefix}_contigs_10K.bed \
        ${hc_prefix}_map_sorted.bam > ${hc_prefix}_coverage_table.tsv
    echo "$(date)    Coverage file is generated"
    
    rm *.sam *bt2
    # Unsupervised binning
    concoct --composition_file ${hc_prefix}_contigs_10K.fa \
        --coverage_file ${hc_prefix}_coverage_table.tsv \
        --threads 4 \
        --length_threshold 500 \
        --read_length 79 \
        -b ${hc_prefix}_concoct/
    echo "$(date)    Unsupervised binning is complete"
    
    merge_cutup_clustering.py ${hc_prefix}_concoct/clustering_gt500.csv > ${hc_prefix}_concoct/clustering_merged.csv
    echo "$(date)    Clusters are merged"
    
    mkdir ${hc_prefix}_concoct/fasta_bins
    
    extract_fasta_bins.py ${hc_prefix}_tmp.fa ${hc_prefix}_concoct/clustering_merged.csv \
        --output_path ${hc_prefix}_concoct/fasta_bins
    echo "$(date)    Bins are stored as individual FASTA"
    
    rm *tmp.fa
    echo "$(date)    Replica: ${repli}$ Sample: ${hc_prefix} COMPLETE"

done   