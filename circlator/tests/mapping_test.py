import unittest
import filecmp
import os
import pyfastaq 
from circlator import mapping

modules_dir = os.path.dirname(os.path.abspath(mapping.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestMapping(unittest.TestCase):
    def test_bwa_mem(self):
        '''test bwa_mem'''
        # FIXME
        pass


    def test_aligned_read_to_read(self):
        '''test aligned_read_to_read'''
        # FIXME
        pass
