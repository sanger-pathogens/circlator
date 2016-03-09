import argparse
import os
import shutil
import sys
import subprocess
import circlator

def run():
    parser = argparse.ArgumentParser(
        description = 'Run Circlator on a small test dataset',
        usage = 'ariba test [options] <outdir>')
    parser.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1, metavar='INT')
    parser.add_argument('outdir', help='Name of output directory')
    options = parser.parse_args()

    print('Running Circlator on test data...')

    try:
        os.mkdir(options.outdir)
        os.chdir(options.outdir)
    except:
        print('Error making output directory "', options.outdir, '". Cannot continue.', sep='', file=sys.stderr)
        sys.exit(1)

    print('Made output directory. Copying test data files into it:')

    modules_dir = os.path.dirname(os.path.abspath(circlator.__file__))
    test_data_dir = os.path.join(modules_dir, 'data')

    for filename in ['test_contigs.fa', 'test_reads.fq.gz']:
        shutil.copy(os.path.join(test_data_dir, filename), filename)
        print('    copied', filename)


    cmd = ' '.join([
        sys.argv[0],
        'all',
        '--threads', str(options.threads),
        '--verbose',
        '--assemble_spades_use_first',
        'test_contigs.fa',
        'test_reads.fq.gz',
        'OUT',
    ])

    print('\nRunning Circlator with:', cmd, '', sep='\n')

    return_code = subprocess.call(cmd, shell=True)

    if return_code != 0:
        print('\nSomething went wrong. See above for error message(s). Return code was', return_code)
        sys.exit(1)

    print('-' * 79)
    print('Finished run on test data OK')
