import os
import tempfile
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
      spades_use_first_success=False,
      spades_no_careful=False,
    ):
        self.outdir = os.path.abspath(outdir)
        self.reads = os.path.abspath(reads)
        if not os.path.exists(self.reads):
            raise Error('Reads file not found:' + self.reads)

        self.verbose = verbose
        self.threads = threads
        self.spades = external_progs.make_and_check_prog('spades', verbose=self.verbose)
        self.spades_kmers = self._build_spades_kmers(spades_kmers)
        self.spades_use_first_success = spades_use_first_success
        self.spades_no_careful = spades_no_careful
        self.samtools = external_progs.make_and_check_prog('samtools', verbose=self.verbose)
        self.assembler = 'spades'


    def _build_spades_kmers(self, kmers):
        if kmers is None:
            return [127,117,107,97,87,77]
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


    def run_spades_once(self, kmer, outdir,spades_no_careful):

        if spades_no_careful:
            careful = ''
        else:
            careful = '--careful'

        cmd = ' '.join([
            self.spades.exe(),
            '-s', self.reads,
            '-k', str(kmer),
            careful,
            '--only-assembler',
            '-t', str(self.threads),
            '-o', outdir,
        ])

        return common.syscall(cmd, verbose=self.verbose, allow_fail=True)


    def run_spades(self, stop_at_first_success=False,spades_no_careful=False):
        '''Runs spades on all kmers. Each a separate run because SPAdes dies if any kmer does
           not work. Chooses the 'best' assembly to be the one with the biggest N50'''
        n50 = {}
        kmer_to_dir = {}

        for k in self.spades_kmers:
            tmpdir = tempfile.mkdtemp(prefix=self.outdir + '.tmp.spades.' + str(k) + '.', dir=os.getcwd())
            kmer_to_dir[k] = tmpdir
            ok, errs = self.run_spades_once(k, tmpdir,spades_no_careful)
            if ok:
                contigs_fasta = os.path.join(tmpdir, 'contigs.fasta')
                contigs_fai = contigs_fasta + '.fai'
                common.syscall(self.samtools.exe() + ' faidx ' + contigs_fasta, verbose=self.verbose)
                stats = pyfastaq.tasks.stats_from_fai(contigs_fai)
                if stats['N50'] != 0:
                    n50[k] = stats['N50']

                    if stop_at_first_success:
                        break

        if len(n50) > 0:
            if self.verbose:
                print('[assemble]\tkmer\tN50')
                for k in sorted(n50):
                    print('[assemble]', k, n50[k], sep='\t')

            best_k = None

            for k in sorted(n50):
                if best_k is None or n50[k] >= n50[best_k]:
                    best_k = k

            assert best_k is not None

            for k, directory in kmer_to_dir.items():
                if k == best_k:
                    if self.verbose:
                        print('[assemble] using assembly with kmer', k)
                    os.rename(directory, self.outdir)
                else:
                    shutil.rmtree(directory)
        else:
            raise Error('Error running SPAdes. Output directories are:\n  ' + '\n  '.join(kmer_to_dir.values()) + '\nThe reason why should be in the spades.log file in each directory.')

    def run(self):
        if self.assembler == 'spades':
            self.run_spades(stop_at_first_success=self.spades_use_first_success,spades_no_careful=self.spades_no_careful)
        else:
            raise Error('Unknown assembler: "' + self.assembler + '". cannot continue')
