import unittest
import filecmp
import os
import shutil
import pyfastaq
from circlator import assemble

modules_dir = os.path.dirname(os.path.abspath(assemble.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestAssembleSpades(unittest.TestCase):
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


    def test_make_spades_command(self):
        '''test _make_spades_command'''
        cmd_start = ' '.join([
            self.assembler.spades.exe(),
            '-s', os.path.join(data_dir, 'assemble_test.dummy_reads.fa'),
        ])

        self.assertEqual(cmd_start + ' -o out -t 1 -k 41 --careful --only-assembler', self.assembler._make_spades_command(41, 'out'))
        self.assertEqual(cmd_start + ' -o out -t 1 -k 43 --careful --only-assembler', self.assembler._make_spades_command(43, 'out'))
        self.assertEqual(cmd_start + ' -o out2 -t 1 -k 41 --careful --only-assembler', self.assembler._make_spades_command(41, 'out2'))

        self.assembler.careful = False
        self.assertEqual(cmd_start + ' -o out -t 1 -k 41 --only-assembler', self.assembler._make_spades_command(41, 'out'))
        self.assembler.careful = True

        self.assembler.only_assembler = False
        self.assertEqual(cmd_start + ' -o out -t 1 -k 41 --careful', self.assembler._make_spades_command(41, 'out'))
        self.assembler.only_assembler = True

        self.assembler.threads = 2
        self.assertEqual(cmd_start + ' -o out -t 2 -k 41 --careful --only-assembler', self.assembler._make_spades_command(41, 'out'))


class TestAssembleCanu(unittest.TestCase):
    def setUp(self):
        self.tmp_assemble_dir = 'tmp.assemble_test'


    def test_make_canu_command(self):
        '''test _make_canu_command'''
        tmp_assemble_dir = 'tmp.assemble_test'
        assembler = assemble.Assembler(
            os.path.join(data_dir, 'assemble_test.dummy_reads.fa'),
            tmp_assemble_dir,
            assembler='canu'
        )

        cmd_start = ' '.join([
            assembler.canu.exe(),
            '-useGrid=false',
            'gnuplotTested=true',
            '-assemble',
            'genomeSize=0.1m',
        ])

        reads = os.path.join(data_dir, 'assemble_test.dummy_reads.fa')
        self.assertEqual(cmd_start + ' -d out -p outname -pacbio-corrected ' + reads, assembler._make_canu_command('out', 'outname'))
        self.assertEqual(cmd_start + ' -d out2 -p outname2 -pacbio-corrected ' + reads, assembler._make_canu_command('out2', 'outname2'))
        assembler.data_type = 'pacbio-raw'
        self.assertEqual(cmd_start + ' -d out -p outname -pacbio-raw ' + reads, assembler._make_canu_command('out', 'outname'))


    def test_rename_canu_contigs(self):
        '''test _rename_canu_contigs'''
        infile = os.path.join(data_dir, 'assemble_test_rename_canu_contigs.in.fa')
        tmpfile = 'tmp.assemble_test_rename_canu_contigs.out.fa'
        expected = os.path.join(data_dir, 'assemble_test_rename_canu_contigs.expect.fa')
        assemble.Assembler._rename_canu_contigs(infile, tmpfile)
        self.assertTrue(filecmp.cmp(expected, tmpfile, shallow=False))
        os.unlink(tmpfile)
