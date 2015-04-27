import argparse
import sys
import pyfastaq
import circlator

def run():
    parser = argparse.ArgumentParser(
        description = 'Fix SNPs and indels in FASTA using BAM file of mapped corrected reads',
        usage = 'circlator fixvars [options] <in.fasta> <in.bam> <outprefix>')
    parser.add_argument('--verbose', action='store_true', help='Be verbose')
    parser.add_argument('fasta_in', help='Name of input FASTA file', metavar='in.fasta')
    parser.add_argument('bam_in', help='Name of input BAM file', metavar='in.bam')
    parser.add_argument('outprefix', help='Prefix of output files')
    options = parser.parse_args()

    fixer = circlator.fixvars.VariantFixer(
        fasta_in=options.fasta_in,
        bam_in=options.bam_in,
        outprefix=options.outprefix,
        verbose=options.verbose
    )
    fixer.run()
