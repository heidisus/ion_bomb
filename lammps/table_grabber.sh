#!/bin/bash

# Exit if no input file is given
if [ -z "$1" ]; then
  echo "Usage: ./table_grabber inputfile outputfile"
  exit 1
fi

# Take in the cmd line args
inputfile=$1
outputfile=$2

# Check the input file exists and is readable
if [ ! -r "$inputfile" ]; then
  echo "Error: Cannot read file '$inputfile'"
  exit 1
fi

# Get the header
grep '^\s*Step' "$inputfile" | head -1 > $outputfile

# Get the values
awk '/^[[:space:]]*[0-9]+[[:space:]]+[0-9.]+[[:space:]]/' "$inputfile" >> $outputfile

# Get the final printed values:
grep '^Area:' "$inputfile" >> $outputfile
grep '^Inserted' >> $outputfile
grep '^Deposited' >> $outputfile
