import argparse
import sys
import pyfastaq
import circlator

def run():
    parser = argparse.ArgumentParser(
        description = 'Make reads from mapping to be reassembled',
        usage = 'circlator bam2reads [options] <in.bam> <outprefix>')
    parser.add_argument('--discard_unmapped', action='store_true', help='Use this to not keep unmapped reads')
    parser.add_argument('--only_contigs', help='File of contig names (one per line). Only reads that map to these contigs are kept (and unmapped reads, unless --discard_unmapped is used).', metavar='FILENAME')
    parser.add_argument('--length_cutoff', type=int, help='All reads mapped to contigs shorter than this will be kept [%(default)s]', default=100000, metavar='INT')
    parser.add_argument('--min_read_length', type=int, help='Minimum length of read to output [%(default)s]', default=250, metavar='INT')
    parser.add_argument('--verbose', action='store_true', help='Be verbose')
    parser.add_argument('bam', help='Name of input bam file', metavar='in.bam')
    parser.add_argument('outprefix', help='Prefix of output filenames')
    options = parser.parse_args()

    bam_filter = circlator.bamfilter.BamFilter(
        options.bam,
        options.outprefix,
        length_cutoff=options.length_cutoff,
        min_read_length=options.min_read_length,
        contigs_to_use=options.only_contigs,
        discard_unmapped=options.discard_unmapped,
        verbose=options.verbose,
    )
    bam_filter.run()

