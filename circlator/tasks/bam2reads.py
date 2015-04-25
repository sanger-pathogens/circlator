import argparse
import sys
import pyfastaq
import circlator

def run():
    parser = argparse.ArgumentParser(
        description = 'Make reads from mapping to be reassembled',
        usage = 'circlator bam2reads [options] <in.bam> <outprefix>')
    parser.add_argument('--length_cutoff', type=int, help='All reads mapped to contigs shorter than this will be kept [%(default)s]', default=100000, metavar='INT')
    parser.add_argument('bam', help='Name of input bam file', metavar='in.bam')
    parser.add_argument('outprefix', help='Prefix of output filenames')
    options = parser.parse_args()

    bam_filter = circlator.bamfilter.BamFilter(
        options.bam,
        options.outprefix,
        length_cutoff=options.length_cutoff
    )
    bam_filter.run()

