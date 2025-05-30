#!/bin/bash

# List of filenames
filenames=("crater_300K_f24_600ps_100eV.dump" "crater_300K_f25_600ps_100eV.dump" "crater_300K_f26_600ps_100eV.dump" "explin_300K_f25_900ps_100eV_r70A_incr5A.dump" "explin_300K_f26_150ps_100eV_r70A_incr20A.dump")

# List of integers
integers=(10 10 10 5 5)

# Get the number of items
len=${#filenames[@]}

# Loop through the arrays using indices
for (( i=0; i<$len; i++ )); do
  filename="${filenames[$i]}"
  integer="${integers[$i]}"
  echo "Running with $filename and $integer"
  python3 pythonscripts/helpers/ovito_temp.py "$filename" "$integer"
done