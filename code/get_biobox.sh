#!/bin/bash

# Variable to call the .py module
path=/home/andbo/mambaforge/envs/pip_amber/bin/src/utils
path_bin=/home/andbo/concoct_arena/01_bench

for i in 10 9 8 7 6 5 4 3 2 1
do
    jobid=$i

    # Different host-contamination level
    hc_level=("00" "01" "02" "03" "04" "05" "06" "07" "08" "090" "095")
    
    # Define prefix based on array id
    hc_prefix=${hc_level[$jobid-1]}

    echo "$(date)    SAMPLE ${hc_prefix} STARTED"    
    
    mkdir ${hc_prefix}

    # Generate the biobox format 
    python ${path}/convert_fasta_bins_to_biobox_format.py \
        ${path_bin}/${hc_prefix}_concoct/fasta_bins/* \
        -o ${hc_prefix}/${hc_prefix}_bin_map.tsv

    echo "$(date)    BIN MAP COMPLETE" 

    # Store the .tsv file
    bin_map=""${hc_prefix}"/"${hc_prefix}"_bin_map.tsv"
    
    # Store the gsb
    gsb="${hc_prefix}_gsa_mapping.tsv"

    # remove the @version and @sampleid lines of the file
    sed -i '2,3d' $bin_map
    
    # Adds @Version and @SampleID that will match with GSA. 
    # Added these lines manually to GSA beforehand
    sed -i "1i @Version:0.9.0\n@SampleID:"${hc_prefix}"" $bin_map

    sed -i "1i @Version:0.9.0\n@SampleID:"${hc_prefix}"" $gsb
    
    echo "$(date)    SAMPLE ID ADDED" 

    # Run AMBER
    amber.py -g ${hc_prefix}_gsa_mapping.tsv \
        ${hc_prefix}/*_bin_*.tsv -o ${hc_prefix}/results

    echo "$(date)    SAMPLE ${hc_prefix} COMPLETE"     
done    