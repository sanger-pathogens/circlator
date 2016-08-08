import unittest
import shutil
import filecmp
import copy
import os
import pyfastaq
from circlator import start_fixer

modules_dir = os.path.dirname(os.path.abspath(start_fixer.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestStartFixer(unittest.TestCase):
    def test_init(self):
        '''Test init'''
        contigs_fa = os.path.join(data_dir, 'start_fixer_init_contigs.fa')
        dna_da = os.path.join(data_dir, 'start_fixer_init_dnaA.fa')

        with self.assertRaises(start_fixer.Error):
            sfixer = start_fixer.StartFixer('notafile', 'outprefix')

        sfixer = start_fixer.StartFixer(contigs_fa, 'outprefix')
        expected_contigs = {
            'contig1': pyfastaq.sequences.Fasta('contig1', 'ACGT'),
            'contig2': pyfastaq.sequences.Fasta('contig2', 'AAAA'),
        }
        self.assertEqual(expected_contigs, sfixer.input_assembly)
        self.assertEqual(set(), sfixer.ignore)

        with self.assertRaises(start_fixer.Error):
            sfixer = start_fixer.StartFixer(contigs_fa, 'outprefix', genes_fa='notafile')

        ignore_file = os.path.join(data_dir, 'start_fixer_init_ignore_ids')
        sfixer = start_fixer.StartFixer(contigs_fa, 'outprefix', ignore=ignore_file)
        self.assertEqual({'ignore_me1', 'ignore_me2'}, sfixer.ignore)


    def test_max_length_from_fasta_file(self):
        '''Test _max_length_from_fasta_file'''
        infile = os.path.join(data_dir, 'start_fixer_max_length_from_fasta_file.fa')
        self.assertEqual(11, start_fixer.StartFixer._max_length_from_fasta_file(infile))


    def test_write_fasta_plus_circularized_ends(self):
        '''Test _write_fasta_plus_circularized_ends'''
        infile = os.path.join(data_dir, 'start_fixer_write_fasta_plus_circularized_ends.in.fa')
        expected = os.path.join(data_dir, 'start_fixer_write_fasta_plus_circularized_ends.out.fa')
        expected_ignore = os.path.join(data_dir, 'start_fixer_write_fasta_plus_circularized_ends.out.ignore.fa')
        tmp_out = 'tmp.test_write_fasta_plus_circularized_ends.fa'
        contigs = {}
        pyfastaq.tasks.file_to_dict(infile, contigs)
        start_fixer.StartFixer._write_fasta_plus_circularized_ends(contigs, tmp_out, 5)
        self.assertTrue(filecmp.cmp(expected, tmp_out, shallow=False))
        start_fixer.StartFixer._write_fasta_plus_circularized_ends(contigs, tmp_out, 5, ignore={'seq1', 'seq4'})
        self.assertTrue(filecmp.cmp(expected_ignore, tmp_out, shallow=False))
        os.unlink(tmp_out)
