import argparse
import sys
import pyfastaq
import circlator

def run(args=None):
    parser = argparse.ArgumentParser(
        description = 'Assemble reads using SPAdes',
        usage = 'circlator assemble [options] <in.reads.fasta> <out_dir>')
    parser.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1, metavar='INT')
    parser.add_argument('--verbose', action='store_true', help='Be verbose')
    parser.add_argument('reads', help='Name of input reads FASTA file', metavar='in.reads.fasta')
    parser.add_argument('out_dir', help='Output directory (must not already exist)')
    options = parser.parse_args(args)

    a = circlator.assemble.Assembler(
        options.reads,
        options.out_dir,
        threads=options.threads,
        verbose=options.verbose
    )
    a.run()

