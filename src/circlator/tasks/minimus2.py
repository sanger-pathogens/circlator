import argparse
import os
import sys
import pyfastaq
import circlator

class Error (Exception): pass


def run():
    parser = argparse.ArgumentParser(
        description = 'Runs minimus2 circularisation pipeline, see https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Circularizing-and-trimming ... this script is a modified version of that protocol. It first runs minimus2 on the input contigs (unless --no_pre_merge is used). Then it tries to circularise each contig one at a time, by breaking it in the middle and using the two pieces as input to minimus2. If minimus2 outputs one contig, then that new one is assumed to be circularised and is kept, otherwise the original contig is kept.',
        usage = 'circlator minimus2 [options] <assembly.fasta> <output prefix>')
    parser.add_argument('--no_pre_merge', action='store_true', help='Do not do initial minimus2 run before trying to circularise each contig')
    parser.add_argument('assembly', help='Name of original assembly', metavar='assembly.fasta')
    parser.add_argument('outprefix', help='Prefix of output files', metavar='output prefix')
    options = parser.parse_args()

    m = circlator.minimus2.Merger(
        options.assembly,
        options.outprefix,
        merge = not options.no_pre_merge
    )

    m.run()
