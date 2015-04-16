import argparse
import os
import sys
import pyfastaq
import circlator

class Error (Exception): pass

class OptSplitAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        print('OptSplitAction:', values)
        a = values.split(',')
        if len(a) % 2 != 0:
            raise Error('Error in options:', a)

        for i in range(0, len(a), 2):
            if len(a[i]) == 1:
                a[i] = '-' + a[i]
            else:
                a[i] = '--' + a[i]
        
        setattr(namespace, self.dest, a)


def add_verbose(l, opts):
    if opts.verbose and '--verbose' not in l:
        l = ['--verbose'] + l
    return l


def add_threads(l, opts):
    if opts.threads and '--threads' not in l:
        l = ['--threads', str(opts.threads)] + l
    return l


def print_message(m, opts):
    if opts.verbose:
        print(m)


def run():
    parser = argparse.ArgumentParser(
        description = 'Run mapreads, bam2reads, assemble, merge',
        usage = 'circlator map [options] <assembly.fasta> <reads.fasta> <output directory>')
    parser.add_argument('--no_pair_merge', action='store_true', help='Do not merge pairs of contigs when running merge task')
    parser.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1, metavar='INT')
    parser.add_argument('--verbose', action='store_true', help='Be verbose')
    parser.add_argument('--mapreads_opts', action=OptSplitAction, help='mapreads options', default=[])
    parser.add_argument('--bam2reads_opts', action=OptSplitAction, help='bam2reads options', default=[])
    parser.add_argument('--assemble_opts', action=OptSplitAction, help='assemble options (there are none yet)', default=[])
    parser.add_argument('--merge_opts', action=OptSplitAction, help='merge options', default=[])
    parser.add_argument('--clean_opts', action=OptSplitAction, help='clean options', default=[])
    parser.add_argument('--fixstart_opts', action=OptSplitAction, help='fixstart options', default=[])
    parser.add_argument('assembly', help='Name of original assembly', metavar='assembly.fasta')
    parser.add_argument('reads', help='Name of corrected reads FASTA file', metavar='reads.fasta')
    parser.add_argument('outdir', help='Name of output directory (must not already exist)', metavar='output directory')
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

    bam = '01.map_corrected_reads.bam'
    filtered_reads_prefix = '02.filtered_reads'
    filtered_reads =  filtered_reads_prefix + '.fasta'
    assembly_dir = '03.Assembly'
    reassembly = os.path.join(assembly_dir, 'scaffolds.fasta')
    merge_prefix = '04.merge'
    merged_fasta = merge_prefix + '.fasta'
    clean_prefix = '05.clean'
    clean_fasta = clean_prefix + '.fasta'
    fixstart_prefix = '06.final'
    fixstart_fasta = fixstart_prefix + '.fasta'

    options.mapreads_opts.extend([original_assembly, original_reads, bam])
    options.mapreads_opts = add_threads(options.mapreads_opts, options)
    options.mapreads_opts = add_verbose(options.mapreads_opts, options)
    print_message('{:_^79}'.format(' Running mapreads '), options)
    print_message('mapreads options:' + ' '.join(options.mapreads_opts), options)
    circlator.tasks.mapreads.run(args=options.mapreads_opts)


    options.bam2reads_opts.extend([bam, filtered_reads_prefix])
    print_message('{:_^79}'.format(' Running bam2reads '), options)
    print_message('bam2reads options:' + ' '.join(options.bam2reads_opts), options)
    circlator.tasks.bam2reads.run(args=options.bam2reads_opts)


    options.assemble_opts.extend([filtered_reads, assembly_dir])
    options.assemble_opts = add_threads(options.assemble_opts, options)
    options.assemble_opts = add_verbose(options.assemble_opts, options)
    print_message('{:_^79}'.format(' Running assemble '), options)
    print_message('assemble options:' + ' '.join(options.assemble_opts), options)
    circlator.tasks.assemble.run(args=options.assemble_opts)


    if not options.no_pair_merge:
        options.merge_opts.extend(['--reads', filtered_reads])
    options.merge_opts.extend([original_assembly, reassembly, merge_prefix])
    options.merge_opts = add_verbose(options.merge_opts, options)
    options.merge_opts = add_threads(options.merge_opts, options)
    print_message('{:_^79}'.format(' Running merge '), options)
    print_message('merge options:' + ' '.join(options.merge_opts), options)
    circlator.tasks.merge.run(args=options.merge_opts)

    merge_log = merge_prefix + '.log'
    contigs_to_keep = []
    contigs_to_not_fix_start = []
    with open(merge_log) as f:
        for line in f:
            name, circularised = line.rstrip().split()
            if circularised == '1': 
                contigs_to_keep.append(name)
            else:
                contigs_to_not_fix_start.append(name)

    keep_file = clean_prefix + '.contigs_to_keep'
    with open(keep_file, 'w') as f:
        print('\n'.join(contigs_to_keep), file=f)
        
    not_fix_start_file = fixstart_prefix + '.contigs_to_not_change'
    with open(not_fix_start_file, 'w') as f:
        print('\n'.join(contigs_to_not_fix_start), file=f)

    options.clean_opts.extend(['--keep', keep_file, merged_fasta, clean_prefix])
    print_message('{:_^79}'.format(' Running clean '), options)
    print_message('clean options:' + ' '.join(options.clean_opts), options)
    circlator.tasks.clean.run(args=options.clean_opts)


    options.fixstart_opts.extend(['--ignore', not_fix_start_file, clean_fasta, fixstart_prefix])
    print_message('{:_^79}'.format(' Running fixstart '), options)
    print_message('fixstart options:' + ' '.join(options.fixstart_opts), options)
    circlator.tasks.fixstart.run(args=options.fixstart_opts)

