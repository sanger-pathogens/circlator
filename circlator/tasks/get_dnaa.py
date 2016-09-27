import argparse
import circlator

class Error (Exception): pass

def run():
    parser = argparse.ArgumentParser(
        description = 'Downloads and filters a file of dnaA (or other) genes from uniprot',
        usage = 'circlator get_dnaa [options] <output prefix>')
    parser.add_argument('--min_length', type=int, help='Minimum length in amino acids [%(default)s]', default=333, metavar='INT')
    parser.add_argument('--max_length', type=int, help='Maximum length in amino acids [%(default)s]', default=500, metavar='INT')
    parser.add_argument('--uniprot_search', help='Uniprot search term [%(default)s]', default='dnaa', metavar='STRING')
    parser.add_argument('--name_re', help='Each sequence name must match this regular expression [%(default)s]', default='Chromosomal replication initiat(or|ion) protein.*dnaa', metavar='STRING')
    parser.add_argument('--name_re_case_sensitive', action='store_true', help='Do a case-sensitive match to regular expression given by --name_re. Default is to ignore case.')

    parser.add_argument('outprefix', help='Prefix of output files', metavar='output prefix')
    options = parser.parse_args()

    dna_getter = circlator.dnaa.UniprotDownloader(
        min_gene_length=options.min_length,
        max_gene_length=options.max_length,
        uniprot_search=options.uniprot_search,
        header_regex=options.name_re,
        header_regex_ignorecase=(not options.name_re_case_sensitive)
    )

    dna_getter.run(options.outprefix)
