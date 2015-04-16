import argparse
import shutil
import os
import sys
import tempfile
import circlator
from bio_assembly_refinement import contig_break_finder 

def run(args=None):
    parser = argparse.ArgumentParser(
        description = 'Change start point of each sequence in assembly',
        usage = 'circlator fixstart [options] <assembly.fasta> <outprefix>')
    parser.add_argument('--verbose', action='store_true', help='Be verbose')
    parser.add_argument('--genes_fa', help='Absolute path to FASTA file of genes to search for to use as start point', metavar='filename')
    parser.add_argument('--ignore', help='Absolute path to file of IDs of contigs to not change', metavar='filename')
    parser.add_argument('assembly_fa', help='Name of input FASTA file', metavar='assembly.fasta')
    parser.add_argument('outprefix', help='Prefix of output files')
    options = parser.parse_args(args)

    if options.genes_fa is None:
        d = os.path.abspath(os.path.dirname(circlator.__file__))
        options.genes_fa = os.path.join(d, 'data', 'dnaA.fasta')
        assert os.path.exists(options.genes_fa)
    else:
        options.genes_fa = os.path.abspath(options.genes_fa)

    options.assembly_fa = os.path.abspath(options.assembly_fa)
    options.outprefix = os.path.abspath(options.outprefix)

    original_dir = os.getcwd()
    tmpdir = tempfile.mkdtemp(prefix=options.outprefix + '.tmp.circlator.fixstart.', dir=os.getcwd())
    os.chdir(tmpdir)

    break_finder = contig_break_finder.ContigBreakFinder(
        fasta_file=options.assembly_fa,
        gene_file=options.genes_fa,
        rename=False,
        avoid=options.ignore,
    )
    break_finder.run()

    # move output files and delete tmp dir
    os.rename(break_finder.output_file, options.outprefix + '.fasta')
    os.rename(break_finder.summary_file, options.outprefix + '.summary')
    os.chdir(original_dir)
    shutil.rmtree(tmpdir)

