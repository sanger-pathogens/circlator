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


    def test_rename_contigs(self):
        '''Test _rename_contigs'''
        contigs_in = {
            'ctg1': pyfastaq.sequences.Fasta('ctg1', 'ACGT'),
            'ctg2 foo bar': pyfastaq.sequences.Fasta('ctg2 foo bar', 'AAA'),
        }
        expected_contigs = {
            'ctg1': pyfastaq.sequences.Fasta('ctg1', 'ACGT'),
            'ctg2': pyfastaq.sequences.Fasta('ctg2', 'AAA'),
        }
        expected_names = {'ctg1': 'ctg1', 'ctg2': 'ctg2 foo bar'}
        got_contigs, got_names = start_fixer.StartFixer._rename_contigs(contigs_in)
        self.assertEqual(expected_names, got_names)
        self.assertEqual(expected_contigs, got_contigs)
        contigs_in['ctg2 abc'] = pyfastaq.sequences.Fasta('ctg2 abc', 'AAA')
        with self.assertRaises(start_fixer.Error):
            start_fixer.StartFixer._rename_contigs(contigs_in)


    def test_write_renamed_contigs(self):
        '''Test _write_renamed_contigs'''
        contigs_dict = {
            'ctg1': pyfastaq.sequences.Fasta('ctg1', 'ACGT'),
            'ctg2': pyfastaq.sequences.Fasta('ctg2', 'AAA'),
        }
        rename_dict = {'ctg1': 'ctg1', 'ctg2': 'ctg2 foo bar'}
        tmp_out = 'tmp.test_write_renamed_contigs.fa'
        start_fixer.StartFixer._write_renamed_contigs(contigs_dict, rename_dict, tmp_out)
        expected = os.path.join(data_dir, 'start_fixer_write_renamed_contigs.fa')
        self.assertTrue(filecmp.cmp(expected, tmp_out, shallow=False))
        os.unlink(tmp_out)


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
        got = start_fixer.StartFixer._write_fasta_plus_circularized_ends(contigs, tmp_out, 5)
        self.assertEqual(6, got)
        self.assertTrue(filecmp.cmp(expected, tmp_out, shallow=False))
        got = start_fixer.StartFixer._write_fasta_plus_circularized_ends(contigs, tmp_out, 5, ignore={'seq1', 'seq4'})
        self.assertEqual(3, got)
        self.assertTrue(filecmp.cmp(expected_ignore, tmp_out, shallow=False))
        os.unlink(tmp_out)
        got = start_fixer.StartFixer._write_fasta_plus_circularized_ends(contigs, tmp_out, 5, ignore={'seq1', 'seq2', 'seq3', 'seq4'})
        self.assertEqual(0, got)
        self.assertFalse(os.path.exists(tmp_out))


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
        got = start_fixer.StartFixer._find_circular_using_promer(tmp_outprefix, ref_genes_fa, contigs_dict, 70, end_extend, sys.stdout, ignore={x for x in contigs_dict})
        self.assertEqual({}, got)


    def test_find_circular_using_prodigal(self):
        '''test _find_circular_using_prodigal'''
        contigs_infile = os.path.join(data_dir, 'start_fixer_find_circular_using_prodigal.ctgs.fa')
        outprefix = 'tmp.test_find_circular_using_prodigal'
        contigs_dict = {}
        pyfastaq.tasks.file_to_dict(contigs_infile, contigs_dict)
        circular_from_promer = {'ctg3': 'foo'}
        got = start_fixer.StartFixer._find_circular_using_prodigal(outprefix, contigs_dict, circular_from_promer, sys.stdout)
        # We should have a dictionary with one key == "ctg1" and a value like this:
        # 'ctg1\tProdigal_v2.60\tCDS\t841\t2244\t103.8\t+\t0\tID=1_1;partial=00;start_type=ATG;rbs_motif=TAAAA;rbs_spacer=4bp;gc_cont=0.486;conf=100.00;score=103.84;cscore=95.60;sscore=8.24;rscore=4.24;uscore=-0.13;tscore=4.12;'
        # ... but this depends on the version of prodigal, so is not predictable
        self.assertEqual({'ctg1'}, set(got.keys()))
        got_fields = got['ctg1'].split('\t')
        self.assertEqual(got_fields[0], 'ctg1')
        self.assertTrue(got_fields[1].startswith('Prodigal'))
        os.unlink(outprefix + '.for_prodigal.fa')
        os.unlink(outprefix + '.prodigal.gff')
        got = start_fixer.StartFixer._find_circular_using_prodigal(outprefix, contigs_dict, circular_from_promer, sys.stdout, ignore={'ctg1', 'ctg2'})
        self.assertEqual({}, got)


    def test_rearrange_contigs(self):
        '''test _rearrange_contigs'''
        # try all the combinations of promer and prodigal and ignore. Also whether or not to
        # reverse complement and DNAA spanning across start/end of contig
        # promer_ctg1: has sp|P03004|DNAA_ECOLI starting at position 2000 on fwd strand
        # promer_ctg2: has sp|P03004|DNAA_ECOLI_2 starting at position 3403 on reverse strand (=revcomp of promer_ctg1)
        # promer_ctg3: has sp|P03004|DNAA_ECOLI starts at 4745, split over end/start of contig
        # promer_ctg4: has sp|P03004|DNAA_ECOLI starts at 644, split over start/end of contig (=revcomp of promer_ctg3)
        # ignore_ctg: same as promer_ctg1, but in the ignore set, so should be unchanged
        # prodigal_ctg1: has a two genes, nearest to middle starts at 4529
        # prodigal_ctg2: same as ctg1, but reverse complement, gene starts at 4862 on rev strand
        contigs_infile = os.path.join(data_dir, 'start_fixer_rearrange_contigs.in.fa')
        ref_genes_fa = os.path.join(data_dir, 'start_fixer_rearrange_contigs.refs.fa')
        tmp_outprefix = 'tmp.test_rearrange_contigs'
        tmp_log = tmp_outprefix + '.log'
        contigs_dict = {}
        pyfastaq.tasks.file_to_dict(contigs_infile, contigs_dict)
        to_ignore = {'ignore_ctg'}
        end_extend = start_fixer.StartFixer._max_length_from_fasta_file(ref_genes_fa)
        circ_with_promer = start_fixer.StartFixer._find_circular_using_promer(tmp_outprefix, ref_genes_fa, contigs_dict, 70, end_extend, sys.stdout, ignore=to_ignore)
        circ_with_prodigal = start_fixer.StartFixer._find_circular_using_prodigal(tmp_outprefix, contigs_dict, circ_with_promer, sys.stdout, ignore=to_ignore)
        start_fixer.StartFixer._rearrange_contigs(contigs_dict, circ_with_promer, circ_with_prodigal, to_ignore, end_extend, tmp_log)
        expected_log = os.path.join(data_dir, 'start_fixer_rearrange_contigs.expect.log')
        self.assertTrue(filecmp.cmp(expected_log, tmp_log, shallow=False))
        expected_dict = {}
        expected_fa = os.path.join(data_dir, 'start_fixer_rearrange_contigs.expect.fa')
        pyfastaq.tasks.file_to_dict(expected_fa, expected_dict)
        self.assertEqual(expected_dict, contigs_dict)
        os.unlink(tmp_outprefix + '.contigs_with_ends.fa')
        os.unlink(tmp_outprefix + '.for_prodigal.fa')
        os.unlink(tmp_outprefix + '.prodigal.gff')
        os.unlink(tmp_outprefix + '.promer')
        os.unlink(tmp_log)


    def test_run_ignore_all(self):
        '''test run when ignoring all sequences'''
        tmp_prefix = 'tmp.start_fixer.test_run_when_ignoring_all'
        input_fa = os.path.join(data_dir, 'start_fixer_run_ignore_all.fa')
        to_ignore = os.path.join(data_dir, 'start_fixer_run_ignore_all.to_ignore')
        sfixer = start_fixer.StartFixer(input_fa, tmp_prefix, ignore=to_ignore)
        sfixer.run()
        expected_fa = os.path.join(data_dir, 'start_fixer_run_ignore_all.expect.fa')
        self.assertTrue(filecmp.cmp(expected_fa, tmp_prefix + '.fasta'))
        for suffix in ['detailed.log', 'fasta', 'log', 'promer.contigs_with_ends.fa', 'promer.promer']:
            try:
                os.unlink(tmp_prefix + '.' + suffix)
            except:
                pass


    def test_run_none_for_prodigal(self):
        '''test run when everthing sorted by promer'''
        tmp_prefix = 'tmp.start_fixer.test_run_none_for_prodigal'
        input_ctg = os.path.join(data_dir, 'start_fixer_run_none_for_prodigal.ctg.fa')
        input_ref = os.path.join(data_dir, 'start_fixer_run_none_for_prodigal.ref.fa')
        sfixer = start_fixer.StartFixer(input_ctg, tmp_prefix, genes_fa=input_ref)
        sfixer.run()
        expected_fa = os.path.join(data_dir, 'start_fixer_run_none_for_prodigal.expect.fa')
        self.assertTrue(filecmp.cmp(expected_fa, tmp_prefix + '.fasta'))
        for suffix in ['detailed.log', 'fasta', 'log', 'promer.contigs_with_ends.fa', 'promer.promer']:
            try:
                os.unlink(tmp_prefix + '.' + suffix)
            except:
                pass


    def test_run_bit_of_everything(self):
        '''test run with promer, prodigal, and ignoring a sequence'''
        tmp_prefix = 'tmp.start_fixer.test_run_bit_of_everything'
        input_ctg = os.path.join(data_dir, 'start_fixer_run_bit_of_everything.ctg.fa')
        input_ref = os.path.join(data_dir, 'start_fixer_run_bit_of_everything.ref.fa')
        input_ignore = os.path.join(data_dir, 'start_fixer_run_bit_of_everything.ignore')
        sfixer = start_fixer.StartFixer(input_ctg, tmp_prefix, genes_fa=input_ref, ignore=input_ignore)
        sfixer.run()
        expected_fa = os.path.join(data_dir, 'start_fixer_run_bit_of_everything.expect.fa')
        self.assertTrue(filecmp.cmp(expected_fa, tmp_prefix + '.fasta'))

        for suffix in ['detailed.log', 'fasta', 'log', 'promer.contigs_with_ends.fa', 'promer.promer', 'prodigal.for_prodigal.fa', 'prodigal.prodigal.gff']:
            try:
                os.unlink(tmp_prefix + '.' + suffix)
            except:
                pass
