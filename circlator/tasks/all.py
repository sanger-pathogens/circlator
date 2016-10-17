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
    parser.add_argument('--unchanged_code', type=int, help='Code to return when the input assembly is not changed [%(default)s]', default=0, metavar='INT')
    parser.add_argument('assembly', help='Name of original assembly', metavar='assembly.fasta')
    parser.add_argument('reads', help='Name of corrected reads FASTA file', metavar='reads.fasta')
    parser.add_argument('outdir', help='Name of output directory (must not already exist)', metavar='output directory')

    mapreads_group = parser.add_argument_group('mapreads options')
    mapreads_group.add_argument('--bwa_opts', help='BWA options, in quotes [%(default)s]', default='-x pacbio', metavar='STRING')

    bam2reads_group = parser.add_argument_group('bam2reads options')
    bam2reads_group.add_argument('--b2r_discard_unmapped', action='store_true', help='Use this to not keep unmapped reads')
    bam2reads_group.add_argument('--b2r_only_contigs', help='File of contig names (one per line). Only reads that map to these contigs are kept (and unmapped reads, unless --b2r_discard_unmapped is used). Note: the whole assembly is still used as a reference when mapping', metavar='FILENAME')
    bam2reads_group.add_argument('--b2r_length_cutoff', type=int, help='All reads mapped to contigs shorter than this will be kept [%(default)s]', default=100000, metavar='INT')
    bam2reads_group.add_argument('--b2r_min_read_length', type=int, help='Minimum length of read to output [%(default)s]', default=250, metavar='INT')

    assemble_group = parser.add_argument_group('assemble options')
    parser.add_argument('--assemble_spades_k', help='Comma separated list of kmers to use when running SPAdes. Max kmer is 127 and each kmer should be an odd integer [%(default)s]', default='127,117,107,97,87,77', metavar='k1,k2,k3,...')
    parser.add_argument('--assemble_spades_use_first', action='store_true', help='Use the first successful SPAdes assembly. Default is to try all kmers and use the assembly with the largest N50')
    parser.add_argument('--assemble_not_careful', action='store_true', help='Do not use the --careful option with SPAdes (used by default)')
    parser.add_argument('--assemble_not_only_assembler', action='store_true', help='Do not use the --assemble-only option with SPAdes (used by default)')

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
    clean_group.add_argument('--clean_min_contig_length', type=int, help='Contigs shorter than this are discarded (unless specified using --keep) [%(default)s]', default=2000, metavar='INT')
    clean_group.add_argument('--clean_min_contig_percent', type=int, help='If length of nucmer hit is at least this percentage of length of contig, then contig is removed. (unless specified using --keep) [%(default)s]', default=95, metavar='FLOAT')
    clean_group.add_argument('--clean_diagdiff', type=int, help='Nucmer diagdiff option [%(default)s]', metavar='INT', default=25)
    clean_group.add_argument('--clean_min_nucmer_id', type=float, help='Nucmer minimum percent identity [%(default)s]', metavar='FLOAT', default=95)
    clean_group.add_argument('--clean_min_nucmer_length', type=int, help='Minimum length of hit for nucmer to report [%(default)s]', metavar='INT', default=500)
    clean_group.add_argument('--clean_breaklen', type=int, help='breaklen option used by nucmer [%(default)s]', metavar='INT', default=500)

    fixstart_group = parser.add_argument_group('fixstart options')
    fixstart_group.add_argument('--genes_fa', help='FASTA file of genes to search for to use as start point. If this option is not used, a built-in set of dnaA genes is used', metavar='FILENAME')
    fixstart_group.add_argument('--fixstart_mincluster', type=int, help='The -c|mincluster option of promer. If this option is used, it overrides promer\'s default value', metavar='INT')
    fixstart_group.add_argument('--fixstart_min_id', type=float, help='Minimum percent identity of promer match between contigs and gene(s) to use as start point [%(default)s]', default=70, metavar='FLOAT')

    options = parser.parse_args()

    print_message('{:_^79}'.format(' Checking external programs '), options)
    if options.verbose:
        circlator.versions.get_all_versions(sys.stdout, raise_error=True)
    else:
        circlator.versions.get_all_versions(None, raise_error=True)


    files_to_check = [options.assembly, options.reads]
    if options.b2r_only_contigs:
        files_to_check.append(options.b2r_only_contigs)
        options.b2r_only_contigs = os.path.abspath(options.b2r_only_contigs)

    if options.genes_fa:
        files_to_check.append(options.genes_fa)

    circlator.common.check_files_exist(files_to_check)

    if options.genes_fa:
        options.genes_fa = os.path.abspath(options.genes_fa)

    original_assembly = os.path.abspath(options.assembly)
    original_reads = os.path.abspath(options.reads)


    try:
        os.mkdir(options.outdir)
    except:
        print('Error making output directory', options.outdir, file=sys.stderr)
        sys.exit(1)

    os.chdir(options.outdir)

    with open('00.info.txt', 'w') as f:
        print(sys.argv[0], 'all', ' '.join(sys.argv[1:]), file=f)
        circlator.versions.get_all_versions(f)

    original_assembly_renamed = '00.input_assembly.fasta'
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

    pyfastaq.tasks.to_fasta(
        original_assembly,
        original_assembly_renamed,
        strip_after_first_whitespace=True,
        check_unique=True
    )

    #-------------------------------- mapreads -------------------------------
    print_message('{:_^79}'.format(' Running mapreads '), options)
    circlator.mapping.bwa_mem(
      original_assembly_renamed,
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
        length_cutoff=options.b2r_length_cutoff,
        min_read_length=options.b2r_min_read_length,
        contigs_to_use=options.b2r_only_contigs,
        discard_unmapped=options.b2r_discard_unmapped,
        verbose=options.verbose,
    )
    bam_filter.run()


    #-------------------------------- assemble -------------------------------
    print_message('{:_^79}'.format(' Running assemble '), options)
    a = circlator.assemble.Assembler(
        filtered_reads,
        assembly_dir,
        threads=options.threads,
        careful=not options.assemble_not_careful,
        only_assembler=not options.assemble_not_only_assembler,
        spades_kmers=options.assemble_spades_k,
        spades_use_first_success=options.assemble_spades_use_first,
        verbose=options.verbose
    )
    a.run()


    #------------------------------ filter original assembly -----------------
    if options.b2r_only_contigs:
        print_message('{:_^79}'.format(' --b2r_only_contigs used - filering contigs '), options)
        assembly_to_use = merge_prefix + '.00.filtered_assembly.fa'
        pyfastaq.tasks.filter(original_assembly_renamed, assembly_to_use, ids_file=options.b2r_only_contigs)
    else:
        assembly_to_use = original_assembly_renamed


    #-------------------------------- merge ----------------------------------
    print_message('{:_^79}'.format(' Running merge '), options)
    if not options.no_pair_merge:
        merge_reads = filtered_reads
    else:
        merge_reads = None

    m = circlator.merge.Merger(
        assembly_to_use,
        reassembly,
        merge_prefix,
        nucmer_diagdiff=options.merge_diagdiff,
        nucmer_min_id=options.merge_min_id,
        nucmer_min_length=options.merge_min_length,
        nucmer_min_length_for_merges=options.merge_min_length_merge,
        min_spades_circular_percent=options.merge_min_spades_circ_pc,
        spades_kmers=options.assemble_spades_k,
        spades_use_first_success=options.assemble_spades_use_first,
        spades_careful=not options.assemble_not_careful,
        spades_only_assembler=not options.assemble_not_only_assembler,
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
        min_contig_length=options.clean_min_contig_length,
        min_contig_percent_match=options.clean_min_contig_percent,
        nucmer_diagdiff=options.clean_diagdiff,
        nucmer_min_id=options.clean_min_nucmer_id,
        nucmer_min_length=options.clean_min_nucmer_length,
        nucmer_breaklen=options.clean_breaklen,
        keepfile=clean_keep_file,
        verbose=options.verbose
    )
    cleaner.run()


    #-------------------------------- fixstart -------------------------------
    print_message('{:_^79}'.format(' Running fixstart '), options)
    fixer = circlator.start_fixer.StartFixer(
        clean_fasta,
        fixstart_prefix,
        min_percent_identity=options.fixstart_min_id,
        promer_mincluster=options.fixstart_mincluster,
        genes_fa=options.genes_fa,
        ignore=not_fix_start_file,
        verbose=options.verbose
    )
    fixer.run()

    #-------------------------------- summary -------------------------------
    print_message('{:_^79}'.format(' Summary '), options)
    number_of_input_contigs = pyfastaq.tasks.count_sequences(original_assembly_renamed)
    final_number_of_contigs = pyfastaq.tasks.count_sequences(fixstart_fasta)
    number_circularized = len(contigs_to_keep)
    print_message('Number of input contigs: ' + str(number_of_input_contigs), options)
    print_message('Number of contigs after merging: ' + str(final_number_of_contigs), options)
    print_message(' '.join(['Circularized', str(number_circularized), 'of', str(final_number_of_contigs), 'contig(s)']), options)

    with open(fixstart_prefix + '.ALL_FINISHED', 'w') as f:
        pass

    if number_of_input_contigs == final_number_of_contigs and number_circularized == 0:
        sys.exit(options.unchanged_code)

