import argparse
import sys
import pyfastaq
import circlator

def run():
    parser = argparse.ArgumentParser(
        description = 'Assemble reads using SPAdes/Canu/Racon',
        usage = 'circlator assemble [options] <in.reads.fasta> <out_dir>')
    parser.add_argument('--not_careful', action='store_true', help='Do not use the --careful option with SPAdes (used by default)')
    parser.add_argument('--not_only_assembler', action='store_true', help='Do not use the --assemble-only option with SPAdes (used by default)')
    parser.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1, metavar='INT')
    parser.add_argument('--verbose', action='store_true', help='Be verbose')
    parser.add_argument('--spades_k', help='Comma separated list of kmers to use when running SPAdes. Max kmer is 127 and each kmer should be an odd integer [%(default)s]', default='127,117,107,97,87,77', metavar='k1,k2,k3,...')
    parser.add_argument('--spades_use_first', action='store_true', help='Use the first successful SPAdes assembly. Default is to try all kmers and use the assembly with the largest N50')
    parser.add_argument('--assembler', choices=circlator.common.allowed_assemblers, help='Assembler to use for reassemblies [%(default)s]', default='spades')
    parser.add_argument('--data_type', choices=circlator.common.allowed_data_types, help='String representing one of the 4 type of data analysed (only used for Canu and Racon) [%(default)s]', default='pacbio-raw')
    parser.add_argument('reads', help='Name of input reads FASTA file', metavar='in.reads.fasta')
    parser.add_argument('out_dir', help='Output directory (must not already exist)')
    options = parser.parse_args()

    a = circlator.assemble.Assembler(
        options.reads,
        options.out_dir,
        threads=options.threads,
        careful=not options.not_careful,
        only_assembler=not options.not_only_assembler,
        spades_kmers=options.spades_k,
        spades_use_first_success=options.spades_use_first,
        assembler=options.assembler,
        data_type=options.data_type,
        verbose=options.verbose
    )
    a.run()
