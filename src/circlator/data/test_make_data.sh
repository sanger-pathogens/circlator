#!/usr/bin/env bash
# This is the script that was used to generate the
# test data, which is used when running 'circlator test'
set -e
set -x

ref=test_ref.fa
for_reads=test_for_reads.fa
reads=test_reads.fq.gz
contigs=test_contigs.fa

fastaq make_random_contigs --seed 42 1 300000 $ref
samtools faidx $ref

echo ">1.twice" > $for_reads.$$
samtools faidx $ref 1 1 | grep -v ">" >> $for_reads.$$
fastaq to_fasta $for_reads.$$ $for_reads
rm $for_reads.$$
fastaq to_perfect_reads --seed 42 $for_reads $reads 16000 1 20 8000


samtools faidx test_ref.fa 1:500-148000 > $contigs
samtools faidx test_ref.fa 1:150000-299500 >> $contigs
