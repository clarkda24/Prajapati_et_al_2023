#! /bin/bash

> compute_Cuts_and_Occ_jobs

# Process all BAM files
ls *.bam | while read x
do 
	exp_name="${x:0:${#x}-4}"
	echo "/data/$USER/scripts/sacCer3/run_Compute_Cuts_and_Occ_single_sacCer3_BAM_file_amin_amax.sh /usr/local/matlab-compiler/v94 $PWD/${exp_name}.bam $1 $2 " >> compute_Cuts_and_Occ_jobs
done

# Run swarm
swarm -g 64 -t 32 -f compute_Cuts_and_Occ_jobs --partition quick,norm --time=4:00:00
