import unittest
import filecmp
import os
import shutil
import pyfastaq
from circlator import minimus2

modules_dir = os.path.dirname(os.path.abspath(minimus2.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

# It's likely that some people will not want to run the minimus2 pipeline at
# all, and just use the circlator pipeline. So we check if amos is installed
# and only run the tests if it is installed.
try:
    merger = minimus2.Merger(
        os.path.join(data_dir, 'minimus2_test.dummy_input_assembly.fa'),
        'test.minimus2.out'
    )
    can_run_tests = True
except:
    can_run_tests = False


@unittest.skipUnless(can_run_tests, 'Amos is not installed')
class TestMinimus2(unittest.TestCase):
    def setUp(self):
        self.merger = minimus2.Merger(
            os.path.join(data_dir, 'minimus2_test.dummy_input_assembly.fa'),
            'test.minimus2.out'
        )


    def test_write_split_contig(self):
        '''test _write_split_contig'''
        contig = pyfastaq.sequences.Fasta('seq', 'ACGT')
        tmp_file = 'tmp.minimus2_test.write_split_contig.out'
        self.merger._write_split_contig(contig, tmp_file)
        expected = os.path.join(data_dir, 'minimus2_test.write_split_contig.fa')
        self.assertTrue(filecmp.cmp(tmp_file, expected, shallow=False))
        os.unlink(tmp_file)


    def test_run_minimus2(self):
        '''test _run_minimus2'''
        input_assembly = os.path.join(data_dir, 'minimus2_test_run_minimus2.in.fa')
        outdir = 'tmp.test_run_minimus2'
        os.mkdir(outdir)
        self.merger._run_minimus2(input_assembly, outdir)
        self.assertTrue(os.path.exists(os.path.join(outdir, 'minimus2.fasta')))
        shutil.rmtree(outdir)


    def test_run_minimus2_on_one_contig_will_not_work(self):
        '''test_run_minimus2_on_one_contig that won't circularise'''
        contig = pyfastaq.sequences.Fasta('seq', 'ACGT')
        new_contig, circularised = self.merger._run_minimus2_on_one_contig(contig)
        self.assertEqual(contig, new_contig)
        self.assertFalse(circularised)


    def test_run_minimus2_on_one_contig_will_work(self):
        '''test_run_minimus2_on_one_contig that will circularise'''
        contigs = {}
        contigs_file = os.path.join(data_dir, 'minimus2_test_run_minimus2_on_one_contig.in.fa')
        pyfastaq.tasks.file_to_dict(contigs_file, contigs)
        assert len(contigs) == 1
        contig = contigs['contig1']
        new_contig, circularised = self.merger._run_minimus2_on_one_contig(contig)
        self.assertTrue(circularised)
        expected_contigs = {}
        expected_contigs_file = os.path.join(data_dir, 'minimus2_test_run_minimus2_on_one_contig.expected.fa')
        pyfastaq.tasks.file_to_dict(expected_contigs_file, expected_contigs)
        assert len(expected_contigs) == 1
        expected_contig = expected_contigs['contig1.circularised']
        self.assertEqual(new_contig, expected_contig)


    def test_run_minimus2_on_each_contig(self):
        '''test _run_minimus2_on_each_contig'''
        input_contigs = os.path.join(data_dir, 'minimus2_test_run_minimus2_on_each_contig.in.fa')
        expected_contigs = os.path.join(data_dir, 'minimus2_test_run_minimus2_on_each_contig.expected.fa')
        tmpfile = 'tmp.test.run_minimus2_on_each_contig.out.fa'
        self.merger._run_minimus2_on_each_contig(input_contigs, tmpfile)
        self.assertTrue(filecmp.cmp(tmpfile, expected_contigs, shallow=False))
        os.unlink(tmpfile)


    def test_run_initial_merge(self):
        '''test _run_initial_merge'''
        input_contigs = os.path.join(data_dir, 'minimus2_test_run_initial_merge.in.fa')
        expected_contigs = os.path.join(data_dir, 'minimus2_test_run_initial_merge.expected.fa')
        tmpdir = 'tmp.test.run_initial_merge.out'
        got_file = self.merger._run_initial_merge(input_contigs, tmpdir)
        self.assertTrue(filecmp.cmp(expected_contigs, got_file, shallow=False))
        os.unlink(got_file)
        shutil.rmtree(tmpdir)

