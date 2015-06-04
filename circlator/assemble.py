import os
import sys
import shutil
import pyfastaq
from circlator import common, external_progs

class Error (Exception): pass

class Assembler:
    def __init__(self,
      reads,
      outdir,
      threads=1,
      spades_kmers=None,
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
        self.spades_kmers = self._build_spades_kmers(spades_kmers)
        self.assembler = 'spades'


    def _build_spades_kmers(self, kmers):
        if kmers is None:
            return [127,121,111,101,95,91,85,81,75,71]
        elif type(kmers) == str:
            try:
                kmer_list = [int(k) for k in kmers.split(',')]
            except:
                raise Error('Error getting list of kmers from:' + str(kmers))
            return kmer_list
        elif type(kmers) == list:
            return kmers
        else:
            raise Error('Error getting list of kmers from:' + str(kmers))


    def run_spades_once(self, kmer, outdir):
        cmd = ' '.join([
            self.spades.exe(),
            '-s', self.reads,
            '-k', str(kmer),
            '--careful',
            '--only-assembler',
            '-t', str(self.threads),
            '-o', outdir,
        ])

        return common.syscall(cmd, verbose=self.verbose, allow_fail=True)


    def run_spades(self):
        '''Runs spades, starting with biggest kmer. Takes the first result where SPAdes returns success.
           We do this because sometimes a large kmer makes SPAdes crash because the coverage is too low.
           Even if you give SPAdes a list of kmers, if one of the kmers doesn't work, then it stops.
           So need to run separately on each kmer'''
        spades_errors = []

        for k in self.spades_kmers:
            ok, errs = self.run_spades_once(k, self.outdir)
            if ok:
                return

            shutil.rmtree(self.outdir)
            spades_errors.append(errs)

        print('Error running SPAdes. Output from all runs follows...', file=sys.stderr)
        print('\n\n______________________________________________\n\n'.join(spades_errors), file=sys.stderr)
        print('\n\nError running SPAdes, most likely due to low read coverage.', file=sys.stderr)
        print('Please see above for the output from SPAdes. Cannot continue', file=sys.stderr)
        sys.exit(1)


    def run(self):
        if self.assembler == 'spades':
            self.run_spades()
        else:
            raise Error('Unknown assembler: "' + self.assembler + '". cannot continue')
