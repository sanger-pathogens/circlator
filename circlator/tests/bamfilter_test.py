import unittest
import filecmp
import os
import pyfastaq
from circlator import bamfilter

modules_dir = os.path.dirname(os.path.abspath(bamfilter.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestBamfilter(unittest.TestCase):
    def test_get_ref_lengths(self):
        '''test _get_ref_lengths'''
        b = bamfilter.BamFilter(os.path.join(data_dir, 'bamfilter_test_get_ref_lengths.bam'), 'out')
        expected = {
            'ref1': 41,
            'ref2': 42,
            'ref3': 43,
        }
        self.assertEqual(expected, b._get_ref_lengths())


    def test_get_contigs_to_use(self):
        '''test _get_contigs_to_use'''
        b = bamfilter.BamFilter(os.path.join(data_dir, 'bamfilter_test_get_contigs_to_use.bam'), 'out')
        test_file = os.path.join(data_dir, 'bamfilter_test_get_contigs_to_use.infile')
        self.assertEqual(b._get_contigs_to_use(test_file), {'contig42', 'contig4444244'})
        self.assertEqual(b._get_contigs_to_use(None), set())
        self.assertEqual(b._get_contigs_to_use({'42', '43'}), {'42', '43'})


    def test_check_contigs_to_use(self):
        '''test _check_contigs_to_use'''
        input_bam = os.path.join(data_dir, 'bamfilter_test_check_contigs_to_use.bam')
        b = bamfilter.BamFilter(input_bam, 'out')
        ref_lengths = b._get_ref_lengths()
        self.assertTrue(b._check_contigs_to_use(ref_lengths))

        b = bamfilter.BamFilter(input_bam, 'out', contigs_to_use={'1'})
        self.assertTrue(b._check_contigs_to_use(ref_lengths))

        b = bamfilter.BamFilter(input_bam, 'out', contigs_to_use={'1', '2'})
        self.assertTrue(b._check_contigs_to_use(ref_lengths))

        with self.assertRaises(bamfilter.Error):
            b = bamfilter.BamFilter(input_bam, 'out', contigs_to_use={'42'})
            self.assertTrue(b._check_contigs_to_use(ref_lengths))


    def test_all_reads_from_contig(self):
        '''test _all_reads_from_contig'''
        b = bamfilter.BamFilter(os.path.join(data_dir, 'bamfilter_test_all_reads_from_contig.bam'), 'out')
        tmp = 'tmp.test_all_reads_from_contig.out.fa'
        f = pyfastaq.utils.open_file_write(tmp)
        expected = os.path.join(data_dir, 'bamfilter_test_all_reads_from_contig.reads.fa')
        b._all_reads_from_contig('1', f)
        pyfastaq.utils.close(f)
        self.assertTrue(filecmp.cmp(expected, tmp, shallow=False))
        os.unlink(tmp)


    def test_get_all_unmapped_reads(self):
        '''test _get_all_unmapped_reads'''
        b = bamfilter.BamFilter(os.path.join(data_dir, 'bamfilter_test_get_all_unmapped_reads.bam'), 'out')
        expected = os.path.join(data_dir, 'bamfilter_test_get_all_unmapped_reads.reads.fa')
        tmp = 'tmp.test_get_all_unmapped_reads.out.fa'
        f = pyfastaq.utils.open_file_write(tmp)
        b._get_all_unmapped_reads(f)
        pyfastaq.utils.close(f)
        self.assertTrue(filecmp.cmp(expected, tmp, shallow=False))
        os.unlink(tmp)


    def test_break_reads(self):
        '''test _break_reads'''
        b = bamfilter.BamFilter(os.path.join(data_dir, 'bamfilter_test_break_reads.bam'), 'out')
        expected = os.path.join(data_dir, 'bamfilter_test_break_reads.broken_reads.fa')
        tmp = 'tmp.test_break_reads.out.fa'
        f = pyfastaq.utils.open_file_write(tmp)
        b._break_reads('contig1', 390, f, min_read_length=5)
        pyfastaq.utils.close(f)
        self.assertTrue(filecmp.cmp(expected, tmp, shallow=False))
        os.unlink(tmp)


    def test_exclude_region(self):
        '''test _exclude_region'''
        b = bamfilter.BamFilter(os.path.join(data_dir, 'bamfilter_test_exclude_region.bam'), 'out')
        expected = os.path.join(data_dir, 'bamfilter_test_exclude_region.reads.fa')
        tmp = 'tmp.test_exclude_reads.out.fa'
        f = pyfastaq.utils.open_file_write(tmp)
        b._exclude_region('1', 500, 700, f)
        pyfastaq.utils.close(f)
        self.assertTrue(filecmp.cmp(expected, tmp, shallow=False))
        os.unlink(tmp)


    def test_get_region_start(self):
        '''test _get_region start'''
        b = bamfilter.BamFilter(os.path.join(data_dir, 'bamfilter_test_get_region_start.bam'), 'out')
        expected = os.path.join(data_dir, 'bamfilter_test_get_region_start.reads.fa')
        tmp = 'tmp.test_get_region.out.fa'
        f = pyfastaq.utils.open_file_write(tmp)
        b._get_region('1', 0, 64, f, min_length=20)
        pyfastaq.utils.close(f)
        self.assertTrue(filecmp.cmp(expected, tmp, shallow=False))
        os.unlink(tmp)


    def test_get_region_end(self):
        '''test _get_region end'''
        b = bamfilter.BamFilter(os.path.join(data_dir, 'bamfilter_test_get_region_end.bam'), 'out')
        expected = os.path.join(data_dir, 'bamfilter_test_get_region_end.reads.fa')
        tmp = 'tmp.test_get_region.out.fa'
        f = pyfastaq.utils.open_file_write(tmp)
        b._get_region('2', 379, 499, f, min_length=20)
        pyfastaq.utils.close(f)
        self.assertTrue(filecmp.cmp(expected, tmp, shallow=False))
        os.unlink(tmp)


    def test_run_keep_unmapped_no_quals(self):
        '''test run keep unmapped bam has no quality scores'''
        outprefix = 'tmp.bamfilter_run'
        b = bamfilter.BamFilter(
            os.path.join(data_dir, 'bamfilter_test_run_no_qual.bam'),
            outprefix,
            length_cutoff=600,
            min_read_length=100,
            contigs_to_use={'contig1', 'contig3', 'contig4'}
        )
        b.run()
        expected = os.path.join(data_dir, 'bamfilter_test_run_keep_unmapped.out.reads.fa')
        self.assertTrue(filecmp.cmp(expected, outprefix + '.fasta', shallow=False))
        os.unlink(outprefix + '.fasta')
        os.unlink(outprefix + '.log')

        b = bamfilter.BamFilter(
            os.path.join(data_dir, 'bamfilter_test_run_no_qual.bam'),
            outprefix,
            fastq_out=True,
            length_cutoff=600,
            min_read_length=100,
            contigs_to_use={'contig1', 'contig3', 'contig4'}
        )
        b.run()
        expected = os.path.join(data_dir, 'bamfilter_test_run_keep_unmapped.out.reads.fa')
        self.assertTrue(filecmp.cmp(expected, outprefix + '.fastq', shallow=False))
        os.unlink(outprefix + '.fastq')
        os.unlink(outprefix + '.log')


    def test_run_keep_unmapped_with_quals(self):
        '''test run keep unmapped bam has quality scores'''
        outprefix = 'tmp.bamfilter_run'
        b = bamfilter.BamFilter(
            os.path.join(data_dir, 'bamfilter_test_run_with_qual.bam'),
            outprefix,
            fastq_out=False,
            length_cutoff=600,
            min_read_length=100,
            contigs_to_use={'contig1', 'contig3', 'contig4'}
        )
        b.run()
        expected = os.path.join(data_dir, 'bamfilter_test_run_keep_unmapped.out.reads.fa')
        self.assertTrue(filecmp.cmp(expected, outprefix + '.fasta', shallow=False))
        os.unlink(outprefix + '.fasta')
        os.unlink(outprefix + '.log')

        b = bamfilter.BamFilter(
            os.path.join(data_dir, 'bamfilter_test_run_with_qual.bam'),
            outprefix,
            fastq_out=True,
            length_cutoff=600,
            min_read_length=100,
            contigs_to_use={'contig1', 'contig3', 'contig4'}
        )
        b.run()
        expected = os.path.join(data_dir, 'bamfilter_test_run_keep_unmapped.out.reads.fq')
        self.assertTrue(filecmp.cmp(expected, outprefix + '.fastq', shallow=False))
        os.unlink(outprefix + '.fastq')
        os.unlink(outprefix + '.log')


    def test_run_discard_unmapped_no_quals(self):
        '''test run keep unmapped bam has no quality scores'''
        outprefix = 'tmp.bamfilter_run'
        b = bamfilter.BamFilter(
            os.path.join(data_dir, 'bamfilter_test_run_no_qual.bam'),
            outprefix,
            length_cutoff=600,
            min_read_length=100,
            contigs_to_use={'contig1', 'contig3', 'contig4'},
            discard_unmapped=True
        )
        b.run()
        expected = os.path.join(data_dir, 'bamfilter_test_run_discard_unmapped.out.reads.fa')
        self.assertTrue(filecmp.cmp(expected, outprefix + '.fasta', shallow=False))
        os.unlink(outprefix + '.fasta')
        os.unlink(outprefix + '.log')

