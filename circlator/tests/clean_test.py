import unittest
import filecmp
import os
import pymummer
from circlator import clean

modules_dir = os.path.dirname(os.path.abspath(clean.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestClean(unittest.TestCase):
    def test_get_contigs_to_keep(self):
        '''test _get_contigs_to_keep'''
        cleaner = clean.Cleaner('infile', 'outprefix')
        self.assertEqual(cleaner._get_contigs_to_keep(None), set())
        test_file = os.path.join(data_dir, 'clean_test_contigs_to_keep')
        self.assertEqual(cleaner._get_contigs_to_keep(test_file), {'contig42', 'contig64738'})


    def test_remove_small_contigs(self):
        '''test _remove_small_contigs'''
        cleaner = clean.Cleaner('infile', 'outprefix', min_contig_length=2)
        infile = os.path.join(data_dir, 'clean_test_remove_small_contigs.in.fa')
        expected_out = os.path.join(data_dir, 'clean_test_remove_small_contigs.out.fa')
        expected_removed = {'contig1'}
        expected_all = {'contig1', 'contig2', 'contig3', 'contig4'}
        got_out = 'tmp.test_remove_small_contigs.out.fa'
        got_all, got_removed = cleaner._remove_small_contigs(infile, got_out, {'contig4'})
        self.assertTrue(filecmp.cmp(expected_out, got_out, shallow=False))
        self.assertEqual(expected_removed, got_removed)
        self.assertEqual(expected_all, got_all)
        os.unlink(got_out)


    def test_load_nucmer_hits(self):
        '''test _load_nucmer_hits'''
        cleaner = clean.Cleaner('infile', 'outprefix')
        expected_lengths = {'contig1': 1000, 'contig2': 1002, 'contig3': 960, 'contig4': 540}
        hits = [
            '\t'.join(['1', '410', '1', '410', '410', '410', '100.00', '1000', '960', '1', '1', 'contig1', 'contig3']),
            '\t'.join(['481', '960', '481', '960', '480', '480', '100.00', '1002', '960', '1', '1', 'contig2', 'contig3']),
            '\t'.join(['481', '780', '241', '540', '300', '300', '100.00', '1002', '540', '1', '1', 'contig2', 'contig4']),
            '\t'.join(['1', '410', '1', '410', '410', '410', '100.00', '960', '1000', '1', '1', 'contig3', 'contig1']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]

        expected_hits = {
            'contig3': [hits[0], hits[1]],
            'contig4': [hits[2]],
            'contig1': [hits[3]],
        }

        coords_file = os.path.join(data_dir, 'clean_test_load_nucmer_hits.coords')
        got_lengths, got_hits = cleaner._load_nucmer_hits(coords_file)
        self.assertEqual(expected_lengths, got_lengths)
        self.assertEqual(expected_hits, got_hits)


    def test_contains(self):
        '''test _contains'''
        cleaner = clean.Cleaner('infile', 'outprefix')
        cleaner.contigs_to_keep = {'contig4'}
        hits = [
            '\t'.join(['1', '95', '1', '95', '95', '95', '100.00', '100', '100', '1', '1', 'contig1', 'contig2']),
            '\t'.join(['1', '95', '1', '95', '95', '95', '94.00', '100', '100', '1', '1', 'contig1', 'contig2']),
            '\t'.join(['2', '95', '1', '94', '94', '94', '100.00', '100', '100', '1', '1', 'contig1', 'contig2']),
            '\t'.join(['1', '95', '1', '95', '95', '95', '100.00', '100', '100', '1', '1', 'contig3', 'contig4']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]

        expected = [
            True,
            False,
            False,
            False,
        ]

        assert len(expected) == len(hits)

        for i in range(len(hits)):
            self.assertEqual(expected[i], cleaner._contains(hits[i]))


    def test_containing_contigs(self):
        '''test _containing_contigs'''
        cleaner = clean.Cleaner('infile', 'outprefix')
        hits = [
            '\t'.join(['1', '95', '1', '95', '95', '95', '100.00', '100', '100', '1', '1', 'contig1', 'contig2']),
            '\t'.join(['1', '95', '1', '95', '95', '95', '94.00', '100', '100', '1', '1', 'contig3', 'contig2']),
            '\t'.join(['2', '95', '1', '94', '94', '94', '100.00', '100', '100', '1', '1', 'contig4', 'contig2']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]

        self.assertEqual(cleaner._containing_contigs(hits), {'contig1'})


    def test_get_containing_contigs(self):
        '''test _get_containing_contigs'''
        cleaner = clean.Cleaner('infile', 'outprefix')
        hits = [
            '\t'.join(['1', '95', '1', '95', '95', '95', '100.00', '100', '100', '1', '1', 'contig1', 'contig2']),
            '\t'.join(['1', '95', '1', '95', '95', '95', '94.00', '100', '100', '1', '1', 'contig3', 'contig2']),
            '\t'.join(['2', '95', '1', '94', '94', '94', '100.00', '100', '100', '1', '1', 'contig4', 'contig2']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        dict_in = {'contig2': hits}
        expected = {'contig2': {'contig1'}}
        self.assertEqual(cleaner._get_containing_contigs(dict_in), expected)


    def test_get_all_containing(self):
        '''test _get_all_containing'''
        cleaner = clean.Cleaner('infile', 'outprefix')
        dict_in = {
            'a': {'b'},
            'b': {'c'},
            'c': {'d'},
            'e': {'f'},
            'f': {'e', 'g'},
        }

        name_and_expected = [
            ('a', {'b','c','d'}),
            ('b', {'c', 'd'}),
            ('c', {'d'}),
            ('e', {'f', 'g'}),
            ('f', {'e', 'g'})
        ]

        for name, expected in name_and_expected:
            self.assertEqual(expected, cleaner._get_all_containing(dict_in, name))


    def test_expand_containing_using_transitivity(self):
        '''test _expand_containing_using_transitivity'''
        cleaner = clean.Cleaner('infile', 'outprefix')
        dict_in = {
            'a': {'b'},
            'b': {'c'},
            'c': {'d'},
            'e': {'f', 'g'},
            'h': {'e'},
            'i': {'j'},
            'k': {'j'},
        }

        expected = {
            'a': {'b', 'c', 'd'},
            'b': {'c', 'd'},
            'c': {'d'},
            'e': {'f', 'g'},
            'h': {'e', 'f', 'g'},
            'i': {'j'},
            'k': {'j'},
        }

        got = cleaner._expand_containing_using_transitivity(dict_in)
        self.assertEqual(expected, got)
		
		

    def test_infinite_recursion(self):
        '''test _infinite_recursion'''
        cleaner = clean.Cleaner('infile', 'outprefix')
        dict_in = {
            'a': {'b'},
            'b': {'c'},
            'c': {'a', 'b'},
        }

        expected = {'a': {'b', 'c'}, 'b': {'a', 'c'}, 'c': {'a', 'b'}}

        got = cleaner._expand_containing_using_transitivity(dict_in)
        self.assertEqual(expected, got)

    def test_collapse_list_of_sets(self):
        '''test _collapse_list_of_sets'''
        cleaner = clean.Cleaner('infile', 'outprefix')
        tests = [
            ([{1}], [{1}]),
            ([{1,2}, {1}], [{1,2}]),
            ([{1,2,3}, {4,5,6}], [{1,2,3}, {4,5,6}]),
            ([{1,2,3}, {4,5,6}, {4,5,6}], [{1,2,3}, {4,5,6}]),
            ([{1,2,3}, {4,5,6}, {1,5,6}], [{1,2,3,4,5,6}])
        ]

        for list_in, expected in tests:
            self.assertEqual(cleaner._collapse_list_of_sets(list_in), expected)


    def test_get_identical_contigs(self):
        '''test _get_identical_contigs'''
        cleaner = clean.Cleaner('infile', 'outprefix')
        tests = [
            ({1: {2}, 2: {1}}, [{1,2}]),
            ({1: {2}, 2: {1}, 3: {1}}, [{1,2}]),
            ({1: {2,3}, 2: {1}, 3: {1}}, [{1,2,3}]),
            ({1: {2,3}, 2: {1}, 3: {1}, 4: {5}, 5: {6}}, [{1,2,3}]),
            ({1: {2,3}, 2: {1}, 3: {1}, 4: {5}, 5: {6}, 6: {4}}, [{1,2,3}])
        ]

        for dict_in, expected in tests:
            self.assertEqual(cleaner._get_identical_contigs(dict_in), expected)


    def test_longest_contig(self):
        '''test _longest_contig'''
        cleaner = clean.Cleaner('infile', 'outprefix')
        lengths = {'a': 42, 'b': 1, 'c': 123}
        names = {'a', 'b', 'c'}
        self.assertEqual('c', cleaner._longest_contig(names, lengths))


    def test_remove_identical_contigs(self):
        '''test _remove_identical_contigs'''
        cleaner = clean.Cleaner('infile', 'outprefix')
        dict_in = {
            'a': {'b'},
            'c': {'d'},
            'd': {'c'},
            'e': {'f', 'g'},
            'f': {'e'},
            'h': {'i', 'j'},
            'i': {'j', 'h'},
            'j': {'i', 'h'},
            'k': {'l'},
            'l': {'m'},
            'n': {'m'},
            'p': {'q', 'r', 's'},
            'q': {'r', 's'},
            'r': {'s'},
        }

        lengths = {
            'a': 10,
            'b': 20,
            'c': 40,
            'd': 42,
            'e': 10,
            'f': 11,
            'g': 10,
            'h': 20,
            'i': 21,
            'j': 20,
            'k': 20,
            'l': 20,
            'm': 20,
            'n': 20,
            'p': 20,
            'q': 20,
            'r': 20,
        }

        expected_contained = {
            'a': {'b'},
            'f': {'g'},
            'k': {'l'},
            'l': {'m'},
            'n': {'m'},
            'p': {'q', 'r', 's'},
            'q': {'r', 's'},
            'r': {'s'},
        }

        expected_replaced = {
            'c': 'd',
            'e': 'f',
            'h': 'i',
            'j': 'i'
        }

        got_contained, got_replaced = cleaner._remove_identical_contigs(dict_in, lengths)
        self.assertEqual(expected_contained, got_contained)
        self.assertEqual(expected_replaced, got_replaced)


    def test_clean_contigs(self):
        '''test _clean_contigs'''
        cleaner = clean.Cleaner('infile', 'outprefix')
        infile = os.path.join(data_dir, 'clean_test_clean_contigs.in.fa')
        expected_outfile = os.path.join(data_dir, 'clean_test_clean_contigs.out.fa')
        outfile = 'tmp.test_clean_contigs.out.fa'
        containing_contigs = {'contig1'}
        replaced_contigs = {'contig2'}
        cleaner._clean_contigs(infile, outfile, containing_contigs, replaced_contigs)
        self.assertTrue(filecmp.cmp(expected_outfile, outfile, shallow=False))
        os.unlink(outfile)


    def test_write_log(self):
        '''test _write_log'''
        cleaner = clean.Cleaner('infile', 'outprefix')
        cleaner.contigs_to_keep = {'2'}
        outfile = 'tmp.test_clean_write_log'
        all_contigs = {'1', '2', '3', '4', '5', '6', '7', '8', '9'}
        small_removed = '3'
        containing_contigs = {'4': {'5', '6'}}
        replaced_contigs = {'7': '8'}
        cleaner._write_log(outfile, '[clean]', all_contigs, small_removed, containing_contigs, replaced_contigs)
        expected = os.path.join(data_dir, 'clean_test_write_log.log')
        self.assertTrue(filecmp.cmp(expected, outfile, shallow=False))
        os.unlink(outfile)
