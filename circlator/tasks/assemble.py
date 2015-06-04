import argparse
import sys
import pyfastaq
import circlator

def run():
    parser = argparse.ArgumentParser(
        description = 'Assemble reads using SPAdes',
        usage = 'circlator assemble [options] <in.reads.fasta> <out_dir>')
    parser.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1, metavar='INT')
    parser.add_argument('--verbose', action='store_true', help='Be verbose')
    parser.add_argument('--spades_k', help='Comma separated list of kmers to use when running SPAdes. Max kmer is 127 and each kmer should be an odd integer [%(default)s]', default='127,121,111,101,95,91,85,81,75,71', metavar='k1,k2,k3,...')
    parser.add_argument('reads', help='Name of input reads FASTA file', metavar='in.reads.fasta')
    parser.add_argument('out_dir', help='Output directory (must not already exist)')
    options = parser.parse_args()

    a = circlator.assemble.Assembler(
        options.reads,
        options.out_dir,
        threads=options.threads,
        spades_kmers=options.spades_k,
        verbose=options.verbose
    )
    a.run()
