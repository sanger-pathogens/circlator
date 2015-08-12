import argparse
import circlator

class Error (Exception): pass

def run():
    parser = argparse.ArgumentParser(
        description = 'Downloads and filters a file of dnaA genes from uniprot',
        usage = 'circlator get_dnaa [options] <output prefix>')
    parser.add_argument('--min_length', type=int, help='Minimum length in amino acids [%(default)s]', default=333)
    parser.add_argument('--max_length', type=int, help='Maximum length in amino acids [%(default)s]', default=500)
    parser.add_argument('outprefix', help='Prefix of output files', metavar='output prefix')
    options = parser.parse_args()

    dna_getter = circlator.dnaa.UniprotDownloader(
        min_gene_length=options.min_length,
        max_gene_length=options.max_length
    )

    dna_getter.run(options.outprefix)
