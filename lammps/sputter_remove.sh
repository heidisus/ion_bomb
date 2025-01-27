#!/bin/bash

for i in {1..49}    
do
    rm Cu_sputter_loops/data/data.sputter_loop_$i    
    rm Cu_sputter_loops/dump/sputter_loop_$i.dump    
    rm Cu_sputter_loops/output/output_sputter$i.tsv    
done
