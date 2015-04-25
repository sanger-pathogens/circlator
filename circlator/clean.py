import os
import shutil
import tempfile
from bio_assembly_refinement import contig_cleanup

class Error (Exception): pass

class Cleaner:
    def __init__(self,
        infile,
        outprefix,
        min_length=2000,
        keepfile=None,
    ):

        self.infile = os.path.abspath(infile)
        self.outprefix = os.path.abspath(outprefix)
        self.min_length = min_length
        if keepfile is not None:
            self.keepfile = os.path.abspath(keepfile)
        else:
            self.keepfile = None


    def run(self):
        original_dir = os.getcwd()
        tmpdir = tempfile.mkdtemp(prefix=self.outprefix + '.tmp.circlator.clean.', dir=os.getcwd())
        os.chdir(tmpdir)

        ccleaner = contig_cleanup.ContigCleanup(
            self.infile,
            cutoff_contig_length=self.min_length,
            skip=self.keepfile
        )
        ccleaner.run()

        # move output files and delete tmp dir
        assert os.path.exists(ccleaner.output_file)
        assert os.path.exists(ccleaner.summary_file)
        os.rename(ccleaner.output_file, self.outprefix + '.fasta')
        os.rename(ccleaner.summary_file, self.outprefix + '.log')
        os.chdir(original_dir)
        shutil.rmtree(tmpdir)
