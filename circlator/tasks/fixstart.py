import argparse
import circlator

def run():
    parser = argparse.ArgumentParser(
        description = 'Change start point of each sequence in assembly',
        usage = 'circlator fixstart [options] <assembly.fasta> <outprefix>')
    parser.add_argument('--genes_fa', help='FASTA file of genes to search for to use as start point. If this option is not used, a built-in set of dnaA genes is used', metavar='FILENAME')
    parser.add_argument('--ignore', help='Absolute path to file of IDs of contigs to not change', metavar='FILENAME')
    parser.add_argument('--mincluster', type=int, help='The -c|mincluster option of promer. If this option is used, it overrides promer\'s default value', metavar='INT')
    parser.add_argument('--min_id', type=float, help='Minimum percent identity of promer match between contigs and gene(s) to use as start point [%(default)s]', default=70, metavar='FLOAT')
    parser.add_argument('--verbose', action='store_true', help='Be verbose')
    parser.add_argument('assembly_fa', help='Name of input FASTA file', metavar='assembly.fasta')
    parser.add_argument('outprefix', help='Prefix of output files')
    options = parser.parse_args()

    fixer = circlator.start_fixer.StartFixer(
        options.assembly_fa,
        options.outprefix,
        min_percent_identity=options.min_id,
        promer_mincluster=options.mincluster,
        genes_fa=options.genes_fa,
        ignore=options.ignore,
        verbose=options.verbose,
    )
    fixer.run()
