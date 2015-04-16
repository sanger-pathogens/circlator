import argparse
import sys
import os
import shutil
import tempfile
from bio_assembly_refinement import contig_cleanup

def run(args=None):
    parser = argparse.ArgumentParser(
        description = 'Clean contigs',
        usage = 'circlator clean [options] <in.fasta> <outprefix>')
    parser.add_argument('--min_length', help='Minimum contig length to keep [%(default)s]', default=2000)
    parser.add_argument('--keep', help='File of contigs names to keep in output file')
    parser.add_argument('--verbose', action='store_true', help='Be verbose')
    parser.add_argument('fasta_in', help='Name of input FASTA file', metavar='in.fasta')
    parser.add_argument('outprefix', help='Prefix of output files')
    options = parser.parse_args(args)

    options.fasta_in = os.path.abspath(options.fasta_in)
    options.outprefix = os.path.abspath(options.outprefix)
    if options.keep:
        options.keep = os.path.abspath(options.keep)

    original_dir = os.getcwd()
    tmpdir = tempfile.mkdtemp(prefix=options.outprefix + '.tmp.circlator.clean.', dir=os.getcwd())
    os.chdir(tmpdir)

    ccleaner = contig_cleanup.ContigCleanup(
        options.fasta_in,
        cutoff_contig_length=options.min_length,
        ids=options.keep,
    )
    ccleaner.run()

    # move output files and delete tmp dir
    os.rename(ccleaner.output_file, options.outprefix + '.fasta')
    os.rename(ccleaner.summary_file, options.outprefix + '.summary')
    os.chdir(original_dir)
    shutil.rmtree(tmpdir)
     

