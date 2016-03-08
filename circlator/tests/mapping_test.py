import unittest
import filecmp
import os
import pyfastaq
from circlator import mapping

modules_dir = os.path.dirname(os.path.abspath(mapping.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestMapping(unittest.TestCase):
    def test_bwa_index(self):
        '''test bwa_index'''
        ref = os.path.join(data_dir, 'mapping_test_bwa_index.fa')
        outprefix = 'tmp.mapping_test.bwa_index'
        mapping.bwa_index(ref, outprefix)
        expected_files = [outprefix + '.' + x for x in ['amb', 'ann', 'bwt', 'pac', 'sa']]

        for filename in expected_files:
            self.assertTrue(os.path.exists(filename))
            os.unlink(filename)


    def test_bwa_index_clean(self):
        '''test bwa_index_clean'''
        files_prefix = 'tmp.test_bwa_index_clean'
        test_files = [files_prefix + '.' + x for x in ['amb', 'ann', 'bwt', 'pac', 'sa']]
        for filename in test_files:
            with open(filename, 'w') as f:
                pass

            self.assertTrue(os.path.exists(filename))

        mapping.bwa_index_clean(files_prefix)

        for filename in test_files:
            self.assertFalse(os.path.exists(filename))


    def test_bwa_mem(self):
        '''test bwa_mem'''
        # FIXME
        pass


    def test_aligned_read_to_read(self):
        '''test aligned_read_to_read'''
        # FIXME
        pass
