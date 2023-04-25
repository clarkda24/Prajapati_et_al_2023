#!/bin/bash
for i in *_S[0-9]_R1*.fastq.gz; do mv "$i" "${i/_S[0-9]_R1_001.fastq.gz/_R1.fastq.gz}"; done
for i in *_S[0-9][0-9]_R1*.fastq.gz; do mv "$i" "${i/_S[0-9][0-9]_R1_001.fastq.gz/_R1.fastq.gz}"; done
for i in *_S[0-9]_R2*.fastq.gz; do mv "$i" "${i/_S[0-9]_R2_001.fastq.gz/_R2.fastq.gz}"; done
for i in *_S[0-9][0-9]_R2*.fastq.gz; do mv "$i" "${i/_S[0-9][0-9]_R2_001.fastq.gz/_R2.fastq.gz}"; done
