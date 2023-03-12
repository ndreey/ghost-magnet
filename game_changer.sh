#!/bin/bash

# Initialize sum_logs variable
sum_logs=0

# Read in abundance1.tsv file and loop over rows
while IFS=$'\t' read -r col1 col2; do
    # Check if value in column 2 is non-zero
    if (( $(echo "$col2 != 0" | bc -l) )); then
        # Add the logarithm of the non-zero value to the sum_logs variable
        sum_logs=$(echo "$sum_logs + l($col2)" | bc -l)
    fi
done < "data/references/simulation_short_read/abundance2.tsv"

# Print out the final sum_logs value
echo $sum_logs
