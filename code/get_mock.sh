#!/usr/bin/env bash

# ----------Issue_1-----------
#        Get Mock Data

# set the output directory for the downloads
output_dir="data/raw"

# Grabs the ten first mock samples using wget.
for i in {0..9}; do
    link="https://frl.publisso.de/data/frl:6425521/plant_associated/short_read/rhimgCAMI2_sample_"$i"_reads.tar.gz"
    wget -P $output_dir $link
done
