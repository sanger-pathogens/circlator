import argparse
import sys
import circlator

def run():
    parser = argparse.ArgumentParser(
        description = 'Map reads using bwa mem',
        usage = 'circlator mapreads [options] <reference.fasta> <reads.fasta> <out.bam>')
    parser.add_argument('--bwa_opts', help='BWA options, in quotes [%(default)s]', default='-x pacbio', metavar='STRING')
    parser.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1, metavar='INT')
    parser.add_argument('--verbose', action='store_true', help='Be verbose')
    parser.add_argument('ref', help='Name of input reference FASTA file', metavar='reference.fasta')
    parser.add_argument('reads', help='Name of corrected reads FASTA file', metavar='reads.fasta')
    parser.add_argument('bam', help='Name of output BAM file', metavar='out.bam')
    options = parser.parse_args()

    circlator.mapping.bwa_mem(
      options.ref,
      options.reads,
      options.bam,
      threads=options.threads,
      bwa_options=options.bwa_opts,
      verbose=options.verbose,
    )
