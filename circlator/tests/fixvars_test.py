import unittest
import filecmp
import os
from circlator import fixvars

modules_dir = os.path.dirname(os.path.abspath(fixvars.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestVariantFixer(unittest.TestCase):
    def test_make_vcf(self):
        '''test _make_vcf'''
        pass


    def test_get_variants_from_vcf(self):
        '''test _get_variants_from_vcf'''
        pass


    def test_fix_variants(self):
        '''test _fix_variants'''
        pass


    def test_fix_variant(self):
        '''test _fix_variant'''
        pass
