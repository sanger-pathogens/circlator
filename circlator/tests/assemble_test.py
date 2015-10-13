import unittest
import filecmp
import os
import shutil
import pyfastaq
from circlator import assemble

modules_dir = os.path.dirname(os.path.abspath(assemble.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestAssemble(unittest.TestCase):
    def setUp(self):
        self.tmp_assemble_dir = 'tmp.assemble_test'
        self.assembler = assemble.Assembler(
            os.path.join(data_dir, 'assemble_test.dummy_reads.fa'),
            self.tmp_assemble_dir
        )


    def test_build_spades_kmers(self):
        '''test _build_spades_kmers'''
        got = self.assembler._build_spades_kmers('127')
        self.assertEqual([127], got)

        got = self.assembler._build_spades_kmers('127,101,111')
        self.assertEqual([127, 101, 111], got)

        got = self.assembler._build_spades_kmers([55, 131, 75])
        self.assertEqual([55, 131, 75], got)

        with self.assertRaises(assemble.Error):
            self.assembler._build_spades_kmers(42)

        with self.assertRaises(assemble.Error):
            self.assembler._build_spades_kmers('shrubbery')

        with self.assertRaises(assemble.Error):
            self.assembler._build_spades_kmers('41,spam')

