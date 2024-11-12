import unittest
import shutil
import filecmp
import copy
import os
import pyfastaq
from circlator import assembly

modules_dir = os.path.dirname(os.path.abspath(assembly.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestAssembly(unittest.TestCase):
    def test_init_file_not_dir(self):
        '''Test _init_file_not_dir'''
        contigs_fasta = os.path.join(data_dir, 'assembly_test_init_contigs.fasta')
        a = assembly.Assembly(contigs_fasta, 'spades')
        self.assertEqual(a.contigs_fasta, contigs_fasta)
        self.assertIsNone(a.contigs_fastg)
        self.assertIsNone(a.assembly_graph_fastg)
        self.assertIsNone(a.assembler_dir)


    def test_init_pre_3_6_1_ok(self):
        '''Test _init_pre_3_6_1_ok'''
        test_dir = os.path.join(data_dir, 'assembly_test_init_spades_pre_3_6_1_ok')
        a = assembly.Assembly(test_dir, 'spades')
        self.assertEqual(a.contigs_fasta, os.path.join(test_dir, 'contigs.fasta'))
        self.assertEqual(a.contigs_fastg, os.path.join(test_dir, 'contigs.fastg'))
        self.assertEqual(a.assembler_dir, test_dir)
        self.assertIsNone(a.assembly_graph_fastg)
        self.assertIsNone(a.contigs_paths)


    def test_init_post_3_6_1_ok(self):
        '''Test _init_post_3_6_1_ok'''
        test_dir = os.path.join(data_dir, 'assembly_test_init_spades_post_3_6_1_ok')
        a = assembly.Assembly(test_dir, 'spades')
        self.assertEqual(a.assembler_dir, test_dir)
        self.assertEqual(a.contigs_fasta, os.path.join(test_dir, 'contigs.fasta'))
        self.assertEqual(a.assembly_graph_fastg, os.path.join(test_dir, 'assembly_graph.fastg'))
        self.assertEqual(a.contigs_paths, os.path.join(test_dir, 'contigs.paths'))
        self.assertIsNone(a.contigs_fastg)


    def test_init_raises_error(self):
        '''Test _init_raises_error'''
        test_dirs = [
            'assembly_test_init_no_fasta',
            'assembly_test_init_spades_pre_3_6_1_no_fastg',
            'assembly_test_init_spades_post_3_6_1_no_contigs_paths',
            'assembly_test_init_spades_post_3_6_1_no_assembly_graph',
            'assembly_test_init_spades_post_3_6_1_no_assembly_graph_or_contig_paths',
        ]

        for test_dir in test_dirs:
            with self.assertRaises(assembly.Error):
                a = assembly.Assembly(os.path.join(data_dir, test_dir), 'spades')


    def test_get_contigs(self):
        '''Test get_contigs'''
        a = assembly.Assembly(os.path.join(data_dir, 'assembly_test_get_contigs.fasta'), 'spades')
        expected = {
            'contig1': pyfastaq.sequences.Fasta('contig1', 'ACGT'),
            'contig2': pyfastaq.sequences.Fasta('contig2', 'AAAA'),
        }
        self.assertEqual(expected, a.get_contigs())


    def test_circular_contigs_from_spades_before_3_6_1(self):
        '''Test _circular_contigs_from_spades_before_3_6_1'''
        fastg_file = os.path.join(data_dir, 'assembly_test_circular_contigs_from_spades_before_3_6_1.fastg')
        got = assembly.Assembly._circular_contigs_from_spades_before_3_6_1(fastg_file)
        expected = {'NODE_1_length_5_cov_42.42_ID_1'}
        self.assertEqual(expected, got)


    def test_spades_contigs_paths_to_dict(self):
        '''Test _spades_contigs_paths_to_dict'''
        contigs_paths = os.path.join(data_dir, 'assembly_test_spades_contigs_paths_to_dict.contigs.paths')
        expected = {
            "NODE_1_length_86545_cov_16.2519_ID_2261": "1-,2+",
            "NODE_1_length_86545_cov_16.2519_ID_2261'": "2-,1+",
            "NODE_2_length_22106_cov_99.9797_ID_2263": "4+,10-,11+,10+,9-,5-,3+,5-",
            "NODE_2_length_22106_cov_99.9797_ID_2263'": "5+,3-,5+,9+,10-,11-,10+,4-",
            "NODE_3_length_21556_cov_27.4797_ID_2265": "6+",
            "NODE_3_length_21556_cov_27.4797_ID_2265'": "6-",
            "NODE_4_length_16761_cov_36.4471_ID_2267": "7+",
            "NODE_4_length_16761_cov_36.4471_ID_2267'": "7-",
            "NODE_5_length_8134_cov_20.6513_ID_2269": "14+",
            "NODE_5_length_8134_cov_20.6513_ID_2269'": "14-",
            "NODE_6_length_6712_cov_6.68474_ID_2271": "15+",
            "NODE_6_length_6712_cov_6.68474_ID_2271'": "15-",
            "NODE_7_length_4806_cov_3.33747_ID_2273": "13+",
            "NODE_7_length_4806_cov_3.33747_ID_2273'": "13-",
            "NODE_8_length_4566_cov_19.8709_ID_2275": "12+",
            "NODE_8_length_4566_cov_19.8709_ID_2275'": "12-",
            "NODE_9_length_3068_cov_67.4016_ID_2277": "8+",
            "NODE_9_length_3068_cov_67.4016_ID_2277'": "8-",
        }
        got = assembly.Assembly._spades_contigs_paths_to_dict(contigs_paths)
        self.assertEqual(expected, got)


    def test_circular_edges_to_edge_numbers_dict(self):
        '''Test _circular_edges_to_edge_numbers_dict'''
        edges = {'NODE_1_length_86545_cov_16.2519_ID_2261', 'NODE_5_length_8134_cov_20.6513_ID_2269'}
        expected = {
            '1': 'NODE_1_length_86545_cov_16.2519_ID_2261',
            '5': 'NODE_5_length_8134_cov_20.6513_ID_2269'
        }
        got = assembly.Assembly._circular_edges_to_edge_numbers_dict(edges)
        self.assertEqual(expected, got)


    def test_circular_contigs_from_spades_after_3_6_1(self):
        '''Test _circular_contigs_from_spades_after_3_6_1'''
        contigs_paths = os.path.join(data_dir, 'assembly_test_circular_contigs_from_spades_after_3_6_1.contigs.paths')
        assembly_graph_fastg = os.path.join(data_dir, 'assembly_test_circular_contigs_from_spades_after_3_6_1.assembly_graph.fastg')
        got = assembly.Assembly._circular_contigs_from_spades_after_3_6_1(assembly_graph_fastg, contigs_paths)
        expected = {
            'NODE_5_length_8134_cov_20.6513_ID_2269',
            'NODE_6_length_6712_cov_6.68474_ID_2271',
            'NODE_7_length_4806_cov_3.33747_ID_2273',
            'NODE_8_length_4566_cov_19.8709_ID_2275',
        }
        self.assertEqual(expected, got)


    def test_circular_contigs_from_canu_gfa(self):
        '''Test _circular_contigs_from_canu_gfa'''
        infile = os.path.join(data_dir, 'assembly_test_circular_contigs_from_canu.gfa')
        got = assembly.Assembly._circular_contigs_from_canu_gfa(infile)
        expected = {'tig00000042'}
        self.assertEqual(expected, got)


    def test_circular_contigs_spades_pre_3_6_1(self):
        '''Test circular_contigs with spades pre 3.6.1'''
        a = assembly.Assembly(os.path.join(data_dir, 'assembly_test_circular_contigs_spades_pre_3_6_1'), 'spades')
        got = a.circular_contigs()
        expected = {'NODE_1_length_5_cov_42.42_ID_1'}
        self.assertEqual(expected, got)


    def test_circular_contigs_spades_post_3_6_1(self):
        '''Test circular_contigs with spades post 3.6.1'''
        a = assembly.Assembly(os.path.join(data_dir, 'assembly_test_circular_contigs_spades_post_3_6_1'), 'spades')
        got = a.circular_contigs()
        expected = {
            'NODE_5_length_8134_cov_20.6513_ID_2269',
            'NODE_6_length_6712_cov_6.68474_ID_2271',
            'NODE_7_length_4806_cov_3.33747_ID_2273',
            'NODE_8_length_4566_cov_19.8709_ID_2275',
        }
        self.assertEqual(expected, got)


    def test_circular_contigs_just_fasta(self):
        '''Test circular_contigs when input is just fasta'''
        a = assembly.Assembly(os.path.join(data_dir, 'assembly_test_circular_contigs_only_contigs.fasta'), 'spades')
        got = a.circular_contigs()
        expected = set()
        self.assertEqual(expected, got)
