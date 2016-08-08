import unittest
import shutil
import filecmp
import copy
import os
import sys
import pyfastaq
import pymummer
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


    def test_find_circular_using_promer(self):
        '''Test _find_circular_using_promer'''
        contigs_infile = os.path.join(data_dir, 'start_fixer_find_circular_using_promer.contigs.fa')
        ref_genes_fa = os.path.join(data_dir, 'start_fixer_find_circular_using_promer.refs.fa')
        tmp_outprefix = 'tmp.start_fixer_find_circular_using_promer'
        contigs_dict = {}
        pyfastaq.tasks.file_to_dict(contigs_infile, contigs_dict)
        end_extend = start_fixer.StartFixer._max_length_from_fasta_file(ref_genes_fa)

        expected = {
            'ctg_DNAA_ECOLI': pymummer.alignment.Alignment('1021\t2424\t1\t1404\t1404\t1404\t100.00\t.\t.\t3404\t1404\t1\t1\tctg_DNAA_ECOLI\tsp|P03004|DNAA_ECOLI\n'),
            'ctg_DNAA1_CHLPN': pymummer.alignment.Alignment('685\t2067\t1\t1383\t1383\t1383\t100.00\t.\t.\t2808\t1383\t1\t1\tctg_DNAA1_CHLPN__ends\tsp|Q9Z8M9|DNAA1_CHLPN\n'),
            'ctg_DNAA_ECOLI_2': pymummer.alignment.Alignment('901\t2304\t1\t1404\t1404\t1404\t100.00\t.\t.\t3044\t1404\t1\t1\tctg_DNAA_ECOLI_2\tsp|P03004|DNAA_ECOLI_2\n'),
}

        got = start_fixer.StartFixer._find_circular_using_promer(tmp_outprefix, ref_genes_fa, contigs_dict, 70, end_extend, sys.stdout)
        self.assertEqual(expected, got)
        os.unlink(tmp_outprefix + '.contigs_with_ends.fa')
        os.unlink(tmp_outprefix + '.promer')


    def test_find_circular_using_prodigal(self):
        '''test _find_circular_using_prodigal'''
        contigs_infile = os.path.join(data_dir, 'start_fixer_find_circular_using_prodigal.ctgs.fa')
        outprefix = 'tmp.test_find_circular_using_prodigal'
        contigs_dict = {}
        pyfastaq.tasks.file_to_dict(contigs_infile, contigs_dict)
        circular_from_promer = {'ctg3'}
        got = start_fixer.StartFixer._find_circular_using_prodigal(outprefix, contigs_dict, circular_from_promer, sys.stdout)
        expected = {
            'ctg1': 'ctg1\tProdigal_v2.60\tCDS\t2429\t2689\t3.3\t-\t0\tID=1_2;partial=00;start_type=ATG;rbs_motif=TAA;rbs_spacer=15bp;gc_cont=0.517;conf=67.99;score=3.28;cscore=-0.50;sscore=3.78;rscore=3.23;uscore=-3.07;tscore=4.12;'
        }
        self.assertEqual(expected, got)
        os.unlink(outprefix + '.for_prodigal.fa')
        os.unlink(outprefix + '.prodigal.gff')

