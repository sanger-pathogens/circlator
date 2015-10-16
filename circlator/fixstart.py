import os
import shutil
import tempfile
import circlator
from bio_assembly_refinement import contig_break_finder

class Error (Exception): pass

class StartFixer:
    def __init__(self,
        assembly_fa,
        outprefix,
        min_percent_identity=70,
        genes_fa=None,
        ignore=None,
    ):
        if genes_fa is None:
            d = os.path.abspath(os.path.dirname(circlator.__file__))
            self.genes_fa = os.path.join(d, 'data', 'dnaA.fasta')
            assert os.path.exists(self.genes_fa)
        else:
            self.genes_fa = os.path.abspath(genes_fa)

        self.min_percent_identity = min_percent_identity
        self.assembly_fa = os.path.abspath(assembly_fa)
        self.outprefix = os.path.abspath(outprefix)

        if ignore is None:
            self.ignore = None
        else:
            self.ignore = os.path.abspath(ignore)


    def run(self):
        original_dir = os.getcwd()
        tmpdir = tempfile.mkdtemp(prefix=self.outprefix + '.tmp.circlator.fixstart.', dir=os.getcwd())
        os.chdir(tmpdir)

        break_finder = contig_break_finder.ContigBreakFinder(
            fasta_file=self.assembly_fa,
            gene_file=self.genes_fa,
            hit_percent_id=self.min_percent_identity,
            skip=self.ignore,
            rename=False,
            summary_prefix='[fixstart]',
        )
        break_finder.run()

        # move output files and delete tmp dir
        assert os.path.exists(break_finder.output_file)
        assert os.path.exists(break_finder.summary_file)
        os.rename(break_finder.output_file, self.outprefix + '.fasta')
        os.rename(break_finder.summary_file, self.outprefix + '.log')
        os.chdir(original_dir)
        shutil.rmtree(tmpdir)
