import argparse
import sys
import pyfastaq
import circlator

def run(args=None):
    parser = argparse.ArgumentParser(
        description = 'Make reads from mapping to be reassembled',
        usage = 'circlator bam2reads [options] <in.bam> <outprefix>')
    parser.add_argument('--short', type=int, help='All reads mapped to contigs shorter than this will be kept [%(default)s]', default=5000, metavar='INT')
    parser.add_argument('--long', type=int, help='All reads mapped to contigs of length between --short and --long will be broken [%(default)s]', default=100000, metavar='INT')
    parser.add_argument('--min_read_len', type=int, help='Minimum length of read to keep [%(default)s]', default=250, metavar='INT')
    parser.add_argument('--end_dist', type=int, help='Fro long contigs, distance to end of contig for which to keep reads [%(default)s]', default=45000, metavar='INT')
    parser.add_argument('bam', help='Name of input bam file', metavar='in.bam')
    parser.add_argument('outprefix', help='Prefix of output filenames')
    options = parser.parse_args(args)

    bam_filter = circlator.bamfilter.BamFilter(
        options.bam,
        options.outprefix,
        min_length_short=options.short,
        min_length_long=options.long,
        min_read_length=options.min_read_len,
        min_remove_end_distance=options.end_dist
    )
    bam_filter.run()

