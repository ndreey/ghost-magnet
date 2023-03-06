#!/usr/bin/env bash

# ----------Issue_1-----------
#       Unpack Mock Data

# set the output directory for the downloads
output_dir="data/raw"

for i in {0..3}; do
    file_path="data/raw/rhimgCAMI2_sample_"$i"_reads.tar.gz"
    tar -zxvf $file_path -C $output_dir
done