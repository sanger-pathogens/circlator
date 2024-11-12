import argparse
import circlator

def run():
    parser = argparse.ArgumentParser(
        description = 'Clean contigs',
        usage = 'circlator clean [options] <in.fasta> <outprefix>')
    parser.add_argument('--min_contig_length', type=int, help='Contigs shorter than this are discarded (unless specified using --keep) [%(default)s]', default=2000, metavar='INT')
    parser.add_argument('--min_contig_percent', type=int, help='If length of nucmer hit is at least this percentage of length of contig, then contig is removed. (unless specified using --keep) [%(default)s]', default=95, metavar='FLOAT')
    parser.add_argument('--diagdiff', type=int, help='Nucmer diagdiff option [%(default)s]', metavar='INT', default=25)
    parser.add_argument('--min_nucmer_id', type=float, help='Nucmer minimum percent identity [%(default)s]', metavar='FLOAT', default=95)
    parser.add_argument('--min_nucmer_length', type=int, help='Minimum length of hit for nucmer to report [%(default)s]', metavar='INT', default=500)
    parser.add_argument('--breaklen', type=int, help='breaklen option used by nucmer [%(default)s]', metavar='INT', default=500)
    parser.add_argument('--keep', help='File of contig names to keep in output file', metavar='FILENAME')
    parser.add_argument('--verbose', action='store_true', help='Be verbose')
    parser.add_argument('fasta_in', help='Name of input FASTA file', metavar='in.fasta')
    parser.add_argument('outprefix', help='Prefix of output files')
    options = parser.parse_args()

    cleaner = circlator.clean.Cleaner(
        options.fasta_in,
        options.outprefix,
        min_contig_length=options.min_contig_length,
        min_contig_percent_match=options.min_contig_percent,
        nucmer_diagdiff=options.diagdiff,
        nucmer_min_id=options.min_nucmer_id,
        nucmer_min_length=options.min_nucmer_length,
        nucmer_breaklen=options.breaklen,
        keepfile=options.keep,
        verbose=options.verbose
    )

    cleaner.run()
