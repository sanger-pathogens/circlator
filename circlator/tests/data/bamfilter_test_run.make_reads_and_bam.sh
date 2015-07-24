#!/usr/bin/env bash -e
reads=bamfilter_test_run.reads.fa
ref=bamfilter_test_run.ref.fa
samtools faidx bamfilter_test_run.ref.fa contig1:1-100     > $reads
samtools faidx bamfilter_test_run.ref.fa contig1:301-400  >> $reads
samtools faidx bamfilter_test_run.ref.fa contig2:1-100    >> $reads
samtools faidx bamfilter_test_run.ref.fa contig4:1-100    >> $reads
samtools faidx bamfilter_test_run.ref.fa contig4:251-350  >> $reads
samtools faidx bamfilter_test_run.ref.fa contig4:201-350  >> $reads
samtools faidx bamfilter_test_run.ref.fa contig4:401-500  >> $reads
samtools faidx bamfilter_test_run.ref.fa contig4:651-750  >> $reads
samtools faidx bamfilter_test_run.ref.fa contig4:651-800  >> $reads
samtools faidx bamfilter_test_run.ref.fa contig4:651-801  >> $reads
samtools faidx bamfilter_test_run.ref.fa contig4:851-950  >> $reads
samtools faidx bamfilter_test_run.ref.fa contig4:901-1000 >> $reads
fastaq make_random_contigs 1 100 - | awk 'BEGIN{getline; print ">unmapped_read"} 1' >> $reads

smalt index -k 9 -s 1 $ref $ref
smalt map $ref $reads | samtools view -bS - > tmp.$$.bam
samtools sort tmp.$$.bam bamfilter_test_run
samtools index bamfilter_test_run.bam
rm tmp.$$.bam $ref.sm{a,i}


