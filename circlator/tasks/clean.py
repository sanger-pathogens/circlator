import argparse
import circlator

def run():
    parser = argparse.ArgumentParser(
        description = 'Clean contigs',
        usage = 'circlator clean [options] <in.fasta> <outprefix>')
    parser.add_argument('--min_length', type=int, help='Minimum contig length to keep [%(default)s]', default=2000, metavar='INT')
    parser.add_argument('--keep', help='File of contig names to keep in output file', metavar='FILENAME')
    parser.add_argument('fasta_in', help='Name of input FASTA file', metavar='in.fasta')
    parser.add_argument('outprefix', help='Prefix of output files')
    options = parser.parse_args()

    cleaner = circlator.clean.Cleaner(
        options.fasta_in,
        options.outprefix,
        min_length=options.min_length,
        keepfile=options.keep
    )
    cleaner.run()
