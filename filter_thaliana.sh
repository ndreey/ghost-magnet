#!/usr/bin/env bash

# ----------Issue_3-----------
# Filter away the A. thaliana (taxid: 3702) reads

# loop through all sample directories
for sample_dir in data/raw/simulation_short_read/*; do

  echo "Processing ${sample_dir}..."

  # get the sample number from the directory name
  sample_num=$(basename "$sample_dir" | cut -d '_' -f 4)

  # create the "clean" directory
  clean_dir="data/processed/clean_sample_${sample_num}"
  mkdir -p "$clean_dir"

  # set input file paths
  reads="${sample_dir}/reads/anonymous_reads.fq.gz"
  map="${sample_dir}/reads/reads_mapping.tsv.gz"

  # set output file paths
  filtered_fastq="${clean_dir}/filtered_reads.fq.gz"
  thaliana_reads="${clean_dir}/thaliana_reads.txt"

  # extract reads_mapping.tsv.gz to reads_mapping.tsv
  gunzip -c "$map" > reads_mapping.tsv

  # filter out A. thaliana reads using the tax_id column
  awk -F $'\t' '$3 == 3702 {print $1}' reads_mapping.tsv > "$thaliana_reads"

  # filter A. thaliana reads in each sample directory
  zcat "$reads" | awk -v RS="@S" '$2 !~ /taxid:3702/ {print "@"$0}' | gzip > "$filtered_fastq"

  # cleanup intermediate files
  rm reads_mapping.tsv 

  echo "Finished processing ${sample_dir}"
  echo " "

done