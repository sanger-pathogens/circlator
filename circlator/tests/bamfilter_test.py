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
        b = bamfilter.BamFilter(os.path.join(data_dir, 'bamfilter_test_break_reads.bam'), 'out', min_read_length=5)
        expected = os.path.join(data_dir, 'bamfilter_test_break_reads.broken_reads.fa')
        tmp = 'tmp.test_break_reads.out.fa'
        f = pyfastaq.utils.open_file_write(tmp) 
        b._break_reads('contig1', 390, f)
        pyfastaq.utils.close(f)
        self.assertTrue(filecmp.cmp(expected, tmp))
        os.unlink(tmp)


    def test_exclude_region(self):
        '''test _exclude_region'''
        b = bamfilter.BamFilter(os.path.join(data_dir, 'bamfilter_test_exclude_region.bam'), 'out')
        expected = os.path.join(data_dir, 'bamfilter_test_exclude_region.reads.fa')
        tmp = 'tmp.test_exclude_reads.out.fa'
        f = pyfastaq.utils.open_file_write(tmp) 
        b._exclude_region('1', 500, 700, f)
        pyfastaq.utils.close(f)
        self.assertTrue(filecmp.cmp(expected, tmp))
        os.unlink(tmp)

