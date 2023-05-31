#!/bin/bash

for i in {2..11}
do
    jobid=$i

    # Different host-contamination level
    hc_level=("00" "01" "02" "03" "04" "05" "06" "07" "08" "090" "095")
    
    # Define prefix based on array id
    hc_prefix=${hc_level[$jobid-1]}
    while read line
    do
        fa_file=$(echo "$line" | awk '{print $1}')
        bin_name=$(echo "$line" | awk '{print $3}')

        if [[ "$bin_name" == "Platanthera_zijinensis_chr" && "$fa_file" =~ \.fa$ ]]
        then
            cat "$fa_file" >> bin_refs/700c_bin/${hc_prefix}_700c.fa
        fi

    done < ${hc_prefix}/results/bin_metrics.tsv
done


# So i dont forget
# We loop through all the AMBER results and access the bin_metrics.tsv
# We extract 1 column adn 3 column values (bin id, abundant genome)
# If the bin id ends with .fa we know its not from golden standard.
# if abundant genome is == to host, we access that bin and concatenate it
# to the samples GHOST-MAG.