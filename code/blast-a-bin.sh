#!/bin/bash


                                           
echo " _   _         _       _____     _   _"
echo "| |_| |___ ___| |_ ___|  _  |___| |_|_|___ "
echo "| . | | .'|_ -|  _|___|     |___| . | |   |"
echo "|___|_|__,|___|_|     |__|__|   |___|_|_|_|"
echo "$(date)      ndreey"
echo ""
echo ""
                                           
# The different replicants
for repli in "x1000c" "x700c" "maxbin"; do

    
    for i in 1 # 1..11}
    do
        
        # Different host-contamination level
        hc_level=("00" "01" "02" "03" "04" "05" "06" "07" "08" "090" "095")
        
        # Define prefix based on array id
        hc_prefix=${hc_level[$i-1]}
        
        # Set the path to the bins
        binpath=/home/andbo/concoct_arena/${repli}/${hc_prefix}_concoct/fasta_bins
    
        if [ ! -d "${repli}" ]; then
            mkdir "${repli}_X00"
        fi
    
        touch ${repli}_X00/X${hc_prefix}_blast.tsv

        echo "######################################################"
        echo "################# blast-a-bin ########################"
        echo "------------------------------------------------------"
        echo "$(date)    Accessed: ${repli} samples"
        echo "$(date)    Sample: X${hc_prefix} ENGAGED"
        echo "######################################################"
        echo ""
    
        for bin in $(ls ${binpath}); do
            echo ""
            echo "$(date)        BLASTn of ${repli}_X00/X${hc_prefix}/${bin}"
    
            blastn -db database/mockdb \
                -query ${binpath}/${bin} \
                -outfmt "6 qseqid sseqid pident evalue bitscore" \
                -out tmp/tmp.tsv \
                -evalue 1e-100 \
                -num_threads 4
            
            # Adding the bin variable
            awk -v bin=b${bin} 'BEGIN {OFS="\t"} {$6=bin; print}' tmp/tmp.tsv > tmp/tmp2.tsv
            
            # Adding to the sample summary file
            echo "$(date)        Writing results to ${repli}_X00/X${hc_prefix}_blast.tsv"        
            cat tmp/tmp2.tsv >> ${repli}_X00/X${hc_prefix}_blast.tsv
            rm tmp/*
    
            echo "$(date)        BLASTn of ${repli}_X00/X${hc_prefix}/${bin} COMPLETE"
            echo ""               
        done
        echo "######################################################"
        echo "################# blast-a-bin ########################"
        echo "------------------------------------------------------"
        echo "$(date)    Sample: X${hc_prefix} COMPLETE"
        echo "######################################################"
        echo ""
    done
        echo "$(date)    ${repli} samples COMPLETE"
        echo ""
done

echo " _   _         _       _____     _   _"
echo "| |_| |___ ___| |_ ___|  _  |___| |_|_|___ "
echo "| . | | .'|_ -|  _|___|     |___| . | |   |"
echo "|___|_|__,|___|_|     |__|__|   |___|_|_|_|"
echo "$(date)      ndreey"

                                           