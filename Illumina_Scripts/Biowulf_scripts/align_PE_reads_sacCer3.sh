#!/bin/bash
> bowtie2_jobs

for file in *_R1.fastq.gz;
do
	file=${file##*/}
	exp_name="${file:0:${#file}-12}"
	echo "bowtie2 -X 5000 --very-sensitive --no-discordant --no-mixed --no-unal \\"$'\n'$'\t'$'\t'"-p 31 \\"$'\n'$'\t'$'\t'"-x bowtie2Indexes/Saccharomyces_cerevisiae/sacCer3/genome \\"$'\n'$'\t'$'\t'"-1 ${exp_name}_R1.fastq.gz \\"$'\n'$'\t'$'\t'"-2 ${exp_name}_R2.fastq.gz \\"$'\n'$'\t'$'\t'" 2>${exp_name}.bowtie2_statistics.log | \\"$'\n'$'\t'"samtools view -bS -f 0x2 -F 0x300 -q 1 - | samtools sort -T tmp_${exp_name} -@ 30 -o ${exp_name}.bam - ; samtools index ${exp_name}.bam"$'\n' >> bowtie2_jobs

	echo >> bowtie2_PE_jobs_sacCer3
done
# Run swarm
swarm -g 64 -t 32 -f bowtie2_jobs --module bowtie,samtools --time=4:00:00
