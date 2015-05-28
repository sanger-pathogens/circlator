import argparse
import os
import sys
import pyfastaq
import circlator

class Error (Exception): pass


def print_message(m, opts):
    if opts.verbose:
        print(m)


def run():
    parser = argparse.ArgumentParser(
        description = 'Run mapreads, bam2reads, assemble, merge, clean, fixstart',
        usage = 'circlator all [options] <assembly.fasta> <reads.fasta> <output directory>')
    parser.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1, metavar='INT')
    parser.add_argument('--verbose', action='store_true', help='Be verbose')
    parser.add_argument('assembly', help='Name of original assembly', metavar='assembly.fasta')
    parser.add_argument('reads', help='Name of corrected reads FASTA file', metavar='reads.fasta')
    parser.add_argument('outdir', help='Name of output directory (must not already exist)', metavar='output directory')

    mapreads_group = parser.add_argument_group('mapreads options')
    mapreads_group.add_argument('--bwa_opts', help='BWA options, in quotes [%(default)s]', default='-x pacbio', metavar='STRING')

    bam2reads_group = parser.add_argument_group('bam2reads options')
    bam2reads_group.add_argument('--b2r_length_cutoff', type=int, help='All reads mapped to contigs shorter than this will be kept [%(default)s]', default=100000, metavar='INT')

    # no assemble options (yet?)
    #assemble_group = parser.add_argument_group('assemble_group')

    merge_group = parser.add_argument_group('merge options')
    merge_group.add_argument('--merge_diagdiff', type=int, help='Nucmer diagdiff option [%(default)s]', metavar='INT', default=25)
    merge_group.add_argument('--merge_min_id', type=float, help='Nucmer minimum percent identity [%(default)s]', metavar='FLOAT', default=95)
    merge_group.add_argument('--merge_min_length', type=int, help='Minimum length of hit for nucmer to report [%(default)s]', metavar='INT', default=500)
    merge_group.add_argument('--merge_min_length_merge', type=int, help='Minimum length of nucmer hit to use when merging [%(default)s]', metavar='INT', default=4000)
    merge_group.add_argument('--merge_min_spades_circ_pc', type=float, help='Min percent of contigs needed to be covered by nucmer hits to spades circular contigs [%(default)s]', metavar='FLOAT', default=95)
    merge_group.add_argument('--merge_breaklen', type=int, help='breaklen option used by nucmer [%(default)s]', metavar='INT', default=500)
    merge_group.add_argument('--merge_ref_end', type=int, help='max distance allowed between nucmer hit and end of input assembly contig [%(default)s]', metavar='INT', default=15000)
    merge_group.add_argument('--merge_reassemble_end', type=int, help='max distance allowed between nucmer hit and end of reassembly contig [%(default)s]', metavar='INT', default=1000)
    merge_group.add_argument('--no_pair_merge', action='store_true', help='Do not merge pairs of contigs when running merge task')

    clean_group = parser.add_argument_group('clean options')
    clean_group.add_argument('--clean_min_length', type=int, help='Minimum contig length to keep [%(default)s]', default=2000, metavar='INT')

    fixstart_group = parser.add_argument_group('fixstart options')
    fixstart_group.add_argument('--genes_fa', help='FASTA file of genes to search for to use as start point', metavar='FILENAME')

    options = parser.parse_args()


    print_message('{:_^79}'.format(' Checking external programs '), options)
    circlator.external_progs.check_all_progs(verbose=options.verbose)

    original_assembly = os.path.abspath(options.assembly)
    original_reads = os.path.abspath(options.reads)

    try:
        os.mkdir(options.outdir)
    except:
        print('Error making output directory', options.outdir, file=sys.stderr)
        sys.exit(1)

    os.chdir(options.outdir)

    bam = '01.mapreads.bam'
    filtered_reads_prefix = '02.bam2reads'
    filtered_reads =  filtered_reads_prefix + '.fasta'
    assembly_dir = '03.assemble'
    reassembly = os.path.join(assembly_dir, 'contigs.fasta')
    merge_prefix = '04.merge'
    merged_fasta = merge_prefix + '.fasta'
    clean_prefix = '05.clean'
    clean_fasta = clean_prefix + '.fasta'
    fixstart_prefix = '06.fixstart'
    fixstart_fasta = fixstart_prefix + '.fasta'


    #-------------------------------- mapreads -------------------------------
    print_message('{:_^79}'.format(' Running mapreads '), options)
    circlator.mapping.bwa_mem(
      original_assembly,
      original_reads,
      bam,
      threads=options.threads,
      bwa_options=options.bwa_opts,
      verbose=options.verbose,
    )


    #-------------------------------- bam2reads ------------------------------
    print_message('{:_^79}'.format(' Running bam2reads '), options)
    bam_filter = circlator.bamfilter.BamFilter(
        bam,
        filtered_reads_prefix,
        length_cutoff=options.b2r_length_cutoff
    )
    bam_filter.run()


    #-------------------------------- assemble -------------------------------
    print_message('{:_^79}'.format(' Running assemble '), options)
    a = circlator.assemble.Assembler(
        filtered_reads,
        assembly_dir,
        threads=options.threads,
        verbose=options.verbose
    )
    a.run()


    #-------------------------------- merge ----------------------------------
    print_message('{:_^79}'.format(' Running merge '), options)
    if not options.no_pair_merge:
        merge_reads = filtered_reads
    else:
        merge_reads = None
        options.merge_opts.extend(['--reads', filtered_reads])

    m = circlator.merge.Merger(
        original_assembly,
        reassembly,
        merge_prefix,
        nucmer_diagdiff=options.merge_diagdiff,
        nucmer_min_id=options.merge_min_id,
        nucmer_min_length=options.merge_min_length,
        nucmer_min_length_for_merges=options.merge_min_length_merge,
        min_spades_circular_percent=options.merge_min_spades_circ_pc,
        nucmer_breaklen=options.merge_breaklen,
        ref_end_tolerance=options.merge_ref_end,
        qry_end_tolerance=options.merge_reassemble_end,
        threads=options.threads,
        verbose=options.verbose,
        reads=merge_reads
    )
    m.run()


    #-------------------------------- clean ----------------------------------
    merge_log = merge_prefix + '.circularise.log'
    contigs_to_keep = []
    contigs_to_not_fix_start = []
    with open(merge_log) as f:
        for line in f:
            if not line.startswith('[merge circularised]\t'):
                continue
            if line.rstrip() == '\t'.join(['[merge circularised]', '#Contig', 'repetitive_deleted', 'circl_using_nucmer', 'circl_using_spades', 'circularised']):
                continue

            x, name, u, y, z, circularised = line.rstrip().split('\t')
            if circularised == '1':
                contigs_to_keep.append(name)
            else:
                contigs_to_not_fix_start.append(name)

    clean_keep_file = clean_prefix + '.contigs_to_keep'
    with open(clean_keep_file, 'w') as f:
        if len(contigs_to_keep) > 0:
            print('\n'.join(contigs_to_keep), file=f)

    not_fix_start_file = fixstart_prefix + '.contigs_to_not_change'
    with open(not_fix_start_file, 'w') as f:
        if len(contigs_to_not_fix_start) > 0:
            print('\n'.join(contigs_to_not_fix_start), file=f)

    print_message('{:_^79}'.format(' Running clean '), options)

    cleaner = circlator.clean.Cleaner(
        merged_fasta,
        clean_prefix,
        min_length=options.clean_min_length,
        keepfile=clean_keep_file
    )
    cleaner.run()


    #-------------------------------- fixstart -------------------------------
    print_message('{:_^79}'.format(' Running fixstart '), options)
    fixer = circlator.fixstart.StartFixer(
        clean_fasta,
        fixstart_prefix,
        genes_fa=options.genes_fa,
        ignore=not_fix_start_file
    )
    fixer.run()
