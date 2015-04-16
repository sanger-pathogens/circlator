import argparse
import shutil
import os
import sys
import tempfile
from bio_assembly_refinement import contig_break_finder 

def run(args=None):
    parser = argparse.ArgumentParser(
        description = 'Change start point of each circular sequence',
        usage = 'circlator assemble [options] <genes.fasta> <assembly.fasta> <outprefix>')
    parser.add_argument('--verbose', action='store_true', help='Be verbose')
    parser.add_argument('genes_fa', help='FASTA file of genes to search for to use as start point', metavar='genes.fasta')
    parser.add_argument('assembly_fa', help='Name of input FASTA file', metavar='assembly.fasta')
    parser.add_argument('outprefix', help='Prefix of output files')
    options = parser.parse_args(args)

    options.genes_fa = os.path.abspath(options.genes_fa)
    options.assembly_fa = os.path.abspath(options.assembly_fa)
    options.outprefix = os.path.abspath(options.outprefix)

    original_dir = os.getcwd()
    tmpdir = tempfile.mkdtemp(prefix=options.outprefix + '.tmp.circlator.fixstart.', dir=os.getcwd())
    os.chdir(tmpdir)

    break_finder = contig_break_finder.ContigBreakFinder(
        fasta_file=options.assembly_fa,
        gene_file=options.genes_fa,
    )
    break_finder.run()

    # move output files and delete tmp dir
    os.rename(break_finder.output_file, options.outprefix + '.fasta')
    os.rename(break_finder.summary_file, options.outprefix + '.summary')
    os.chdir(original_dir)
    shutil.rmtree(tmpdir)

