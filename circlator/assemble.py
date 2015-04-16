import os
import pyfastaq
from circlator import common, external_progs

class Error (Exception): pass

class Assembler:
    def __init__(self,
      reads,
      outdir,
      threads=1,
      spades_kmer=127,
      verbose=False,
    ):
        self.outdir = os.path.abspath(outdir)
        try:
            os.mkdir(self.outdir)
        except:
            raise Error('Error mkdir ' + self.outdir)

        self.reads = os.path.abspath(reads)
        if not os.path.exists(self.reads):
            raise Error('Reads file not found:' + self.reads)
 
        self.verbose = verbose
        self.threads = threads
        self.spades = external_progs.make_and_check_prog('spades', verbose=self.verbose)
        self.spades_kmer = spades_kmer
        self.assembler = 'spades'


    def run_spades(self):
        cmd = ' '.join([
            self.spades.exe(),
            '-s', self.reads,
            '-k', str(self.spades_kmer),
            '--careful', 
            '--only-assembler',
            '-t', str(self.threads),
            '-o', self.outdir,
        ])

        common.syscall(cmd, verbose=self.verbose)


    def run(self):
        if self.assembler == 'spades':
            self.run_spades()
        else:
            raise Error('Unknown assembler: "' + self.assembler + '". cannot continue')
         
