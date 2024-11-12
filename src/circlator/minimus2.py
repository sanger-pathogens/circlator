import os
import sys
import shutil
import tempfile
import pyfastaq
import logging
from circlator import common

class Error (Exception): pass

class Merger:
    def __init__(self,
      assembly_fasta,
      outprefix,
      merge=True,
    ):
        self.input_assembly = os.path.abspath(assembly_fasta)
        self.outprefix = os.path.abspath(outprefix)
        self.run_pre_merge = merge

        if not os.path.exists(self.input_assembly):
            raise Error('Input assembly "' + self.input_assembly + '" not found. Cannot continue')

        # need to have toAmos and minimus2 in user's path. We are
        # not using external_progs.py for this because toAmos and minimus2
        # don't output version information. We'll just check we can find
        # them and assume the version is OK
        for program in ['toAmos', 'minimus2']:
            self.check_path(program)

        self.log_file = self.outprefix + '.log'


    def check_path(self, program):
        if shutil.which(program) is None:
            print(program + ' not found in path! Cannot continue')
            raise Error(program + ' not found in path! Cannot continue')


    def _write_split_contig(self, contig, outfile):
        half_coord = int(len(contig) / 2)
        seq1 = contig.subseq(0, half_coord)
        seq2 = contig.subseq(half_coord, len(contig))
        seq1.id += '.1'
        seq2.id += '.2'
        f = pyfastaq.utils.open_file_write(outfile)
        print(seq1, file=f)
        print(seq2, file=f)
        pyfastaq.utils.close(f)


    def _run_minimus2(self, infile, outdir):
        amos_afg_prefix = os.path.join(outdir, 'minimus2')
        amos_afg = amos_afg_prefix + '.afg'
        cmd = 'toAmos -s ' + infile + ' -o ' + amos_afg
        common.syscall(cmd)
        cmd = 'minimus2 ' + amos_afg_prefix
        return common.syscall(cmd, allow_fail=True)


    def _run_minimus2_on_one_contig(self, contig):
        tmpdir = tempfile.mkdtemp(prefix=self.outprefix + 'tmp.minimus2.', dir=os.getcwd())
        split_fa = os.path.join(tmpdir, 'in.fasta')
        self._write_split_contig(contig, split_fa)
        finished_ok, errors = self._run_minimus2(split_fa, tmpdir)

        if not finished_ok:
            logging.info('Minimus2 returned non-zero error code. Keeping original contig ' + contig.id)
            shutil.rmtree(tmpdir)
            return contig, False

        singletons = {}
        circularised = {}
        pyfastaq.tasks.file_to_dict(os.path.join(tmpdir, 'minimus2.singletons.seq'), singletons)
        pyfastaq.tasks.file_to_dict(os.path.join(tmpdir, 'minimus2.fasta'), circularised)
        contig_was_circularised = False

        if len(singletons) == 0 and len(circularised) == 1:
            new_contig = list(circularised.values())[0]
            new_contig.id = contig.id + '.circularised'
            contig_was_circularised = True
        elif len(singletons) == 2:
            new_contig = contig
        else:
            print('Error running minimus2 on contig', contig.id, '- there were', len(singletons), 'singletons and', len(circularised), 'circularised contigs. Cannot continue', file=sys.stderr)
            print('minimus2 was run in this directory:', tmpdir, file=sys.stderr)
            raise Error('Error running minimus2. Cannot continue')

        shutil.rmtree(tmpdir)
        return new_contig, contig_was_circularised


    def _run_minimus2_on_each_contig(self, infile, outfile):
        seq_reader = pyfastaq.sequences.file_reader(infile)
        f_out = pyfastaq.utils.open_file_write(outfile)
        circularised_contigs = []

        for contig in seq_reader:
            logging.info('Circularising: ' + contig.id + ' ...')
            new_contig, circularised = self._run_minimus2_on_one_contig(contig)
            print(new_contig, file=f_out)
            if circularised:
                circularised_contigs.append(contig.id)
            logging.info(contig.id + ' circularised: ' + str(circularised))

        pyfastaq.utils.close(f_out)
        logging.info('Circularised contigs:\t' + '\t'.join([x for x in circularised_contigs]))


    def _run_initial_merge(self, infile, outdir):
        try:
            os.mkdir(outdir)
        except:
            raise Error('Error making directory ' + outdir)

        original_contigs = {}
        pyfastaq.tasks.file_to_dict(infile, original_contigs)
        finished_ok, errors = self._run_minimus2(infile, outdir)
        logging.info('First minimus2 run OK: ' +  str(finished_ok))

        if finished_ok:
            merge_fasta = self.outprefix + '.merge.fasta'
            merged_contigs = {}
            singleton_contigs = {}
            pyfastaq.tasks.file_to_dict(os.path.join(outdir, 'minimus2.fasta'), merged_contigs)
            pyfastaq.tasks.file_to_dict(os.path.join(outdir, 'minimus2.singletons.seq'), singleton_contigs)

            contigs_that_were_merged = [x for x in original_contigs if x not in singleton_contigs]
            logging.info('Number of original contigs: ' + str(len(original_contigs)))
            logging.info('Number of singleton contigs: ' + str(len(singleton_contigs)))
            logging.info('Number of contigs that were merged: ' +  str(len(contigs_that_were_merged)))
            logging.info('Merged contigs: ' + '\t'.join(merged_contigs))

            f_out = pyfastaq.utils.open_file_write(merge_fasta)
            for name, contig in merged_contigs.items():
                print(contig, file=f_out)
            for name, contig in singleton_contigs.items():
                print(contig, file=f_out)
            pyfastaq.utils.close(f_out)

            logging.info('Number of contigs to try to circularise: ' + str(len(merged_contigs) + len(singleton_contigs)))
            return merge_fasta
        else:
            print('Error running merge. Going to use original contigs. Minimus2 errors follow...\n', file=sys.stderr)
            print(errors, file=sys.stderr)
            print('... end of Minimus2 output.\n\n', file=sys.stderr)
            logging.info('Number of original contigs: ' + str(len(original_contigs)))
            logging.info('Number of singleton contigs: NA')
            logging.info('Number of contigs that were merged: NA')
            logging.info('Merged contigs: NA')
            return infile


    def run(self):
        logging.basicConfig(
            filename=self.log_file,
            format='%(message)s',
            filemode='w',
            level=logging.INFO
        )

        if self.run_pre_merge:
            merge_dir = self.outprefix + '.merge'
            merge_fasta = self._run_initial_merge(self.input_assembly, merge_dir)
        else:
            merge_fasta = self.input_assembly

        final_fasta = self.outprefix + '.circularise.fasta'
        self._run_minimus2_on_each_contig(merge_fasta, final_fasta)

