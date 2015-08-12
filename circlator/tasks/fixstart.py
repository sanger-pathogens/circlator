import argparse
import circlator

def run():
    parser = argparse.ArgumentParser(
        description = 'Change start point of each sequence in assembly',
        usage = 'circlator fixstart [options] <assembly.fasta> <outprefix>')
    parser.add_argument('--genes_fa', help='FASTA file of genes to search for to use as start point', metavar='FILENAME')
    parser.add_argument('--ignore', help='Absolute path to file of IDs of contigs to not change', metavar='FILENAME')
    parser.add_argument('--min_id', type=float, help='Minimum percent identity of promer match to dnaA gene [%(default)s]', default=70, metavar='FLOAT')
    parser.add_argument('assembly_fa', help='Name of input FASTA file', metavar='assembly.fasta')
    parser.add_argument('outprefix', help='Prefix of output files')
    options = parser.parse_args()

    fixer = circlator.fixstart.StartFixer(
        options.assembly_fa,
        options.outprefix,
        min_percent_identity=options.min_id,
        genes_fa=options.genes_fa,
        ignore=options.ignore
    )
    fixer.run()
