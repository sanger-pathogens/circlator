import unittest
import shutil
import filecmp
import copy
import os
import pymummer
import pyfastaq
from circlator import merge

modules_dir = os.path.dirname(os.path.abspath(merge.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestMerge(unittest.TestCase):
    def setUp(self):
        self.merger = merge.Merger(
            os.path.join(data_dir, 'merge_test_original.fa'),
            os.path.join(data_dir, 'merge_test_reassembly.fa'),
            'tmp.merge_test'
        )


    def test_load_nucmer_hits(self):
        '''test _load_nucmer_hits'''
        got = self.merger._load_nucmer_hits(os.path.join(data_dir, 'merge_test_load_nucmer_hits.coords'))
        lines = [
            '\t'.join(['61', '500', '61', '500', '440', '440', '100.00', '500', '500', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['10', '50', '11', '52', '51', '52', '99.42', '500', '500', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['1', '500', '1', '499', '500', '499', '99.40', '500', '499', '1', '1', 'ref2', 'qry2', '[IDENTITY]'])
        ]
        expected = {
            'ref1': [
                pymummer.alignment.Alignment(lines[0]),
                pymummer.alignment.Alignment(lines[1]),
            ],
            'ref2': [
                pymummer.alignment.Alignment(lines[2]),
            ]
        }
        self.assertEqual(got, expected)


    def test_hits_hashed_by_query(self):
        '''test _hits_hashed_by_query'''
        hits = [
            '\t'.join(['61', '500', '61', '500', '440', '440', '100.00', '500', '500', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['10', '50', '11', '52', '51', '52', '99.42', '500', '500', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['1', '500', '1', '499', '500', '499', '99.40', '500', '499', '1', '1', 'ref2', 'qry2', '[IDENTITY]'])
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        expected = {
            'qry1': [hits[0], hits[1]],
            'qry2': [hits[2]]
        }
        got = self.merger._hits_hashed_by_query(hits)
        self.assertEqual(got, expected)


    def test_get_longest_hit_by_ref_length(self):
        '''test _get_longest_hit_by_ref_length'''
        hits = [
            '\t'.join(['42', '542', '43', '543', '500', '502', '100.00', '500', '500', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['42', '542', '43', '543', '501', '500', '100.00', '500', '500', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['42', '542', '43', '543', '499', '500', '100.00', '500', '500', '1', '1', 'ref1', 'qry1']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        expected = hits[1]
        got = self.merger._get_longest_hit_by_ref_length(hits)
        self.assertEqual(expected, got)


    def test_is_at_ref_start(self):
        '''test _is_at_ref_start'''
        hits_at_start = [
            '\t'.join(['42', '542', '43', '543', '501', '500', '100.00', '500', '500', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['15000', '15542', '1', '543', '542', '500', '100.00', '500', '500', '1', '1', 'ref1', 'qry1']),
        ]
        hits_at_start = [pymummer.alignment.Alignment(x) for x in hits_at_start]

        hits_not_at_start = [
            '\t'.join(['19001', '20000', '1', '1000', '1000', '1000', '100.00', '20000', '500', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['15001', '15542', '1', '543', '542', '542', '100.00', '20000', '500', '1', '1', 'ref1', 'qry1']),
        ]
        hits_not_at_start = [pymummer.alignment.Alignment(x) for x in hits_not_at_start]

        for hit in hits_at_start:
            self.assertTrue(self.merger._is_at_ref_start(hit))

        for hit in hits_not_at_start:
            self.assertFalse(self.merger._is_at_ref_start(hit))


    def test_is_at_ref_end(self):
        '''test _is_at_ref_end'''
        hits_at_end = [
            '\t'.join(['19000', '20000', '1', '1000', '1000', '1000', '100.00', '20000', '1000', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['3999', '5001', '1', '1000', '1000', '1000', '100.00', '20000', '1000', '1', '1', 'ref1', 'qry1']),
        ]
        hits_at_end = [pymummer.alignment.Alignment(x) for x in hits_at_end]

        hits_not_at_end = [
            '\t'.join(['1000', '2000', '1', '1000', '1000', '1000', '100.00', '20000', '1000', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['3999', '5000', '1', '1000', '1000', '1000', '100.00', '20000', '1000', '1', '1', 'ref1', 'qry1']),
        ]
        hits_not_at_end = [pymummer.alignment.Alignment(x) for x in hits_not_at_end]

        for hit in hits_at_end:
            self.assertTrue(self.merger._is_at_ref_end(hit))

        for hit in hits_not_at_end:
            self.assertFalse(self.merger._is_at_ref_end(hit))


    def test_is_at_qry_start(self):
        '''test _is_at_qry_start'''
        hits_at_start = [
            '\t'.join(['42', '1042', '1', '1000', '1000', '1000', '100.00', '2000', '4242', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['1', '1000', '1000', '2000', '1000', '1000', '100.00', '2000', '4242', '1', '1', 'ref1', 'qry1']),
        ]
        hits_at_start = [pymummer.alignment.Alignment(x) for x in hits_at_start]

        hits_not_at_start = [
            '\t'.join(['1', '1000', '1001', '2000', '1000', '1000', '100.00', '2000', '4242', '1', '1', 'ref1', 'qry1']),
        ]
        hits_not_at_start = [pymummer.alignment.Alignment(x) for x in hits_not_at_start]

        for hit in hits_at_start:
            self.assertTrue(self.merger._is_at_qry_start(hit))

        for hit in hits_not_at_start:
            self.assertFalse(self.merger._is_at_qry_start(hit))


    def test_is_at_qry_end(self):
        '''test _is_at_qry_end'''
        hits_at_end = [
            '\t'.join(['42', '1042', '3242', '4242', '1000', '1000', '100.00', '2000', '4242', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['42', '1042', '2224', '3243', '1000', '1000', '100.00', '2000', '4242', '1', '1', 'ref1', 'qry1']),
        ]
        hits_at_end = [pymummer.alignment.Alignment(x) for x in hits_at_end]

        hits_not_at_end = [
            '\t'.join(['1', '1000', '2224', '3242', '1000', '1000', '100.00', '2000', '4242', '1', '1', 'ref1', 'qry1']),
        ]
        hits_not_at_end = [pymummer.alignment.Alignment(x) for x in hits_not_at_end]

        for hit in hits_at_end:
            self.assertTrue(self.merger._is_at_qry_end(hit))

        for hit in hits_not_at_end:
            self.assertFalse(self.merger._is_at_qry_end(hit))


    def test_get_hit_nearest_ref_start(self):
        '''test _get_hit_nearest_ref_start'''
        hits = [
            '\t'.join(['42', '1042', '1', '1000', '1001', '1001', '100.00', '20000', '4242', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['1', '1000', '1000', '2000', '1000', '1000', '100.00', '20000', '4242', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['16000', '18000', '2000', '4000', '2000', '2000', '100.00', '20000', '4242', '1', '1', 'ref1', 'qry1']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        expected = hits[1]
        self.assertEqual(expected, self.merger._get_hit_nearest_ref_start(hits))


    def test_get_hit_nearest_ref_end(self):
        '''test _get_hit_nearest_ref_end'''
        hits = [
            '\t'.join(['42', '1042', '1', '1000', '1001', '1001', '100.00', '20000', '4242', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['16000', '18000', '2000', '4000', '2000', '2000', '100.00', '20000', '4242', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['1', '1000', '1000', '2000', '1000', '1000', '100.00', '20000', '4242', '1', '1', 'ref1', 'qry1']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        expected = hits[1]
        self.assertEqual(expected, self.merger._get_hit_nearest_ref_end(hits))


    def test_get_longest_hit_at_ref_start(self):
        '''test _get_longest_hit_at_ref_start'''
        hits = [
            '\t'.join(['42', '1042', '1', '1000', '1001', '1001', '100.00', '20000', '4242', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['1', '1000', '1000', '2000', '1000', '1000', '100.00', '20000', '4242', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['16000', '18000', '2000', '4000', '2000', '2000', '100.00', '20000', '4242', '1', '1', 'ref1', 'qry1']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        expected = hits[0]
        self.assertEqual(expected, self.merger._get_longest_hit_at_ref_start(hits))


    def test_get_longest_hit_at_ref_end(self):
        '''test _get_longest_hit_at_ref_end'''
        hits = [
            '\t'.join(['15000', '18001', '1', '3000', '3001', '3001', '100.00', '20000', '4242', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['15001', '18001', '1', '3000', '3000', '3000', '100.00', '20000', '4242', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['1', '4000', '1000', '4000', '4000', '4000', '100.00', '20000', '4242', '1', '1', 'ref1', 'qry1']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        expected = hits[0]
        self.assertEqual(expected, self.merger._get_longest_hit_at_ref_end(hits))


    def test_get_longest_hit_at_qry_start(self):
        '''test _get_longest_hit_at_qry_start'''
        hits = [
            '\t'.join(['42', '1042', '1', '1000', '1001', '1001', '100.00', '20000', '4242', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['42', '1042', '1', '999', '1000', '1000', '100.00', '20000', '4242', '1', '1', 'ref2', 'qry1']),
            '\t'.join(['16000', '18000', '2000', '4000', '2000', '2000', '100.00', '20000', '4242', '1', '1', 'ref1', 'qry1']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        expected = hits[0]
        self.assertEqual(expected, self.merger._get_longest_hit_at_qry_start(hits))


    def test_get_longest_hit_at_qry_end(self):
        '''test _get_longest_hit_at_qry_end'''
        hits = [
            '\t'.join(['1', '1000', '40000', '42000', '2000', '2000', '100.00', '20000', '42000', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['1', '1001', '39999', '42000', '2001', '2001', '100.00', '20000', '42000', '1', '1', 'ref2', 'qry1']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        expected = hits[1]
        self.assertEqual(expected, self.merger._get_longest_hit_at_qry_end(hits))


    def test_hits_have_same_query(self):
        '''test _hits_have_same_query'''
        hits = [
            '\t'.join(['42', '43', '42', '43', '2', '2', '100.00', '424242', '4242', '1', '1', 'ref42', 'qry42']),
            '\t'.join(['42', '43', '42', '43', '2', '2', '100.00', '424242', '4242', '1', '1', 'ref43', 'qry42']),
            '\t'.join(['42', '43', '42', '43', '2', '2', '100.00', '424242', '4242', '1', '1', 'ref44', 'qry43']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        self.assertTrue(self.merger._hits_have_same_query(hits[0], hits[1]))
        self.assertFalse(self.merger._hits_have_same_query(hits[0], hits[2]))


    def test_hits_have_same_reference(self):
        '''test _hits_have_same_reference'''
        hits = [
            '\t'.join(['42', '43', '42', '43', '2', '2', '100.00', '424242', '4242', '1', '1', 'ref42', 'qry42']),
            '\t'.join(['42', '43', '42', '43', '2', '2', '100.00', '424242', '4242', '1', '1', 'ref42', 'qry43']),
            '\t'.join(['42', '43', '42', '43', '2', '2', '100.00', '424242', '4242', '1', '1', 'ref43', 'qry44']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        self.assertTrue(self.merger._hits_have_same_reference(hits[0], hits[1]))
        self.assertFalse(self.merger._hits_have_same_reference(hits[0], hits[2]))


    def test_min_qry_hit_length(self):
        '''test _min_qry_hit_length'''
        hits = [
            '\t'.join(['42', '43', '42', '43', '2', '2', '100.00', '424242', '4242', '1', '1', 'ref42', 'qry42']),
            '\t'.join(['42', '44', '42', '44', '3', '3', '100.00', '424242', '4242', '1', '1', 'ref43', 'qry42']),
            '\t'.join(['42', '45', '42', '45', '4', '4', '100.00', '424242', '4242', '1', '1', 'ref44', 'qry43']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        self.assertEqual(2, self.merger._min_qry_hit_length(hits))


    def test_has_qry_hit_longer_than(self):
        '''test _has_qry_hit_longer_than'''
        hits = [
            '\t'.join(['42', '43', '42', '43', '2', '2', '100.00', '424242', '4242', '1', '1', 'ref42', 'qry42']),
            '\t'.join(['42', '44', '42', '44', '3', '3', '100.00', '424242', '4242', '1', '1', 'ref43', 'qry42']),
            '\t'.join(['42', '45', '42', '45', '4', '4', '100.00', '424242', '4242', '1', '1', 'ref44', 'qry43']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        self.assertTrue(self.merger._has_qry_hit_longer_than(hits, 3))
        self.assertTrue(self.merger._has_qry_hit_longer_than(hits, 3, hits_to_exclude={hits[1]}))
        self.assertFalse(self.merger._has_qry_hit_longer_than(hits, 3, hits_to_exclude={hits[2]}))
        self.assertFalse(self.merger._has_qry_hit_longer_than(hits, 4))


    def test_can_circularise(self):
        '''test _can_circularise'''
        start_hits = [
            '\t'.join(['1', '1000', '1', '1000', '1000', '1000', '100.00', '420000', '4000', '1', '1', 'ref42', 'qry1']),
            '\t'.join(['1', '1000', '3000', '4000', '1000', '1000', '100.00', '420000', '4000', '1', '1', 'ref42', 'qry1']),
            '\t'.join(['1000', '1', '1', '1000', '1000', '1000', '100.00', '420000', '4000', '1', '1', 'ref42', 'qry1']),
            '\t'.join(['1000', '1', '3000', '4000', '1000', '1000', '100.00', '420000', '4000', '1', '1', 'ref42', 'qry1']),
        ]
        end_hits = [
            '\t'.join(['419000', '420000', '1', '1000', '1000', '1000', '100.00', '420000', '4000', '1', '1', 'ref42', 'qry1']),
            '\t'.join(['419000', '420000', '3000', '4000', '1000', '1000', '100.00', '420000', '4000', '1', '1', 'ref42', 'qry1']),
            '\t'.join(['420000', '419000', '1', '1000', '1000', '1000', '100.00', '420000', '4000', '1', '1', 'ref42', 'qry1']),
            '\t'.join(['420000', '419000', '3000', '4000', '1000', '1000', '100.00', '420000', '4000', '1', '1', 'ref42', 'qry1']),
        ]
        start_hits = [pymummer.alignment.Alignment(x) for x in start_hits]
        end_hits = [pymummer.alignment.Alignment(x) for x in end_hits]

        good_pairs = {(1, 0), (2, 3)}

        for i in range(len(start_hits)):
            for j in range(len(end_hits)):
                expected = (i, j) in good_pairs
                got = self.merger._can_circularise(start_hits[i], end_hits[j])
                self.assertEqual(expected, got)


    def test_get_possible_circular_ref_contigs(self):
        '''test _get_possible_circular_ref_contigs'''
        hits = [
            '\t'.join(['1', '1000', '3000', '4000', '1000', '1000', '100.00', '420000', '4000', '1', '1', 'ref42', 'qry1']),
            '\t'.join(['419000', '420000', '1', '1000', '1000', '1000', '100.00', '420000', '4000', '1', '1', 'ref42', 'qry1']),
            '\t'.join(['420000', '419000', '3000', '4000', '1000', '1000', '100.00', '420000', '4000', '1', '1', 'ref2', 'qry2']),
            '\t'.join(['1', '1000', '3000', '4000', '1000', '1000', '100.00', '420000', '4000', '1', '1', 'ref43', 'qry1']),
            '\t'.join(['419000', '420000', '1', '1000', '1000', '1000', '100.00', '420000', '4000', '1', '1', 'ref43', 'qry1']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        all_hits = {
            'ref42': [hits[0], hits[1]],
            'ref2': [hits[2]],
            'ref43': [hits[3], hits[4]],
        }

        expected = {
            'ref42': (hits[0], hits[1]),
            'ref43': (hits[3], hits[4]),
        }

        got = self.merger._get_possible_circular_ref_contigs(all_hits)
        self.assertEqual(expected, got)


    def test_get_possible_circular_ref_contigs_same_best_hits(self):
        '''test _get_possible_circular_ref_contigs same best hits'''
        hits = [
            '\t'.join(['1', '37978', '37978', '1', '37978', '37978', '99.98', '2177527', '68234', '1', '1', 'ref', 'qry']),
            '\t'.join(['2137113', '2177527', '68234', '27820', '40415', '40415', '99.95', '2177527', '68234', '1', '1', 'ref', 'qry']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        all_hits = {'ref': hits}
        expected = {'ref': (hits[0], hits[1])}
        got = self.merger._get_possible_circular_ref_contigs(all_hits)
        self.assertEqual(expected, got)


    def test_remove_keys_from_dict_with_nonunique_values(self):
        '''test _remove_keys_from_dict_with_nonunique_values'''
        d = {1:1, 2:1, 3:42, 42:43}
        expected = {3:42, 42:43}
        got = self.merger._remove_keys_from_dict_with_nonunique_values(d)
        self.assertEqual(expected, got)


    def test_make_circularised_contig_1(self):
        '''test _make_circularised_contig 1'''
        hits = [
            '\t'.join(['5', '360', '334', '689', '356', '356', '100.00', '817', '689', '1', '1', 'ref', 'reassembly']),
            '\t'.join(['481', '813', '1', '333', '333', '333', '100.00', '817', '689', '1', '1', 'ref', 'reassembly'])
        ]
        start_hit, end_hit = [pymummer.alignment.Alignment(x) for x in hits]
        ref_fasta = os.path.join(data_dir, 'merge_test_hits_to_new_seq.1.ref.fa')
        reassembly_fasta = os.path.join(data_dir, 'merge_test_hits_to_new_seq.1.reassembly.fa')
        expected_fasta = os.path.join(data_dir, 'merge_test_hits_to_new_seq.1.expected.fa')
        self.merger.original_contigs = {}
        self.merger.reassembly_contigs = {}
        pyfastaq.tasks.file_to_dict(ref_fasta, self.merger.original_contigs)
        pyfastaq.tasks.file_to_dict(reassembly_fasta, self.merger.reassembly_contigs)
        got = self.merger._make_circularised_contig(start_hit, end_hit)
        expected_ref_contigs = {}
        pyfastaq.tasks.file_to_dict(expected_fasta, expected_ref_contigs)
        self.assertEqual(list(expected_ref_contigs.values())[0].seq, got.seq)


    def test_make_circularised_contig_2(self):
        '''test _make_circularised_contig 2'''
        # same as test_make_circularised_contig_1 but reassembly is reverse complemented
        hits = [
            '\t'.join(['5', '360', '356', '1', '356', '356', '100.00', '817', '689', '1', '-1', 'ref', 'reassembly']),
            '\t'.join(['481', '813', '689', '357', '333', '333', '100.00', '817', '689', '1', '-1', 'ref', 'reassembly']),
        ]
        start_hit, end_hit = [pymummer.alignment.Alignment(x) for x in hits]
        ref_fasta = os.path.join(data_dir, 'merge_test_hits_to_new_seq.2.ref.fa')
        reassembly_fasta = os.path.join(data_dir, 'merge_test_hits_to_new_seq.2.reassembly.fa')
        expected_fasta = os.path.join(data_dir, 'merge_test_hits_to_new_seq.2.expected.fa')
        self.merger.original_contigs = {}
        self.merger.reassembly_contigs = {}
        pyfastaq.tasks.file_to_dict(ref_fasta, self.merger.original_contigs)
        pyfastaq.tasks.file_to_dict(reassembly_fasta, self.merger.reassembly_contigs)
        got = self.merger._make_circularised_contig(start_hit, end_hit)
        expected_ref_contigs = {}
        pyfastaq.tasks.file_to_dict(expected_fasta, expected_ref_contigs)
        self.assertEqual(list(expected_ref_contigs.values())[0].seq, got.seq)


    def test_orientation_ok_to_bridge_contigs(self):
        '''test _orientation_ok_to_bridge_contigs'''
        hits = [
            '\t'.join(['15000', '20000', '1', '5000', '5000', '5000', '100.00', '20000', '10000', '1', '-1', 'ref1', 'qry']),
            '\t'.join(['1', '5000', '5000', '10000', '5000', '5000', '100.00', '20000', '10000', '1', '-1', 'ref2', 'qry']),
            '\t'.join(['15000', '20000', '5000', '1', '5000', '5000', '100.00', '20000', '10000', '1', '-1', 'ref1', 'qry']),
            '\t'.join(['5000', '1', '5000', '10000', '5000', '5000', '100.00', '20000', '10000', '1', '-1', 'ref2', 'qry']),
            '\t'.join(['1', '5000', '5000', '1', '5000', '5000', '100.00', '20000', '10000', '1', '-1', 'ref2', 'qry']),
            '\t'.join(['15000', '20000', '10000', '5000', '5000', '5000', '100.00', '20000', '10000', '1', '-1', 'ref1', 'qry']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        self.assertTrue(self.merger._orientation_ok_to_bridge_contigs(hits[0], hits[1]))
        self.assertFalse(self.merger._orientation_ok_to_bridge_contigs(hits[0], hits[0]))
        self.assertFalse(self.merger._orientation_ok_to_bridge_contigs(hits[0], hits[2]))
        self.assertFalse(self.merger._orientation_ok_to_bridge_contigs(hits[0], hits[3]))
        self.assertTrue(self.merger._orientation_ok_to_bridge_contigs(hits[4], hits[5]))


    def test_get_possible_query_bridging_contigs(self):
        '''test _get_possible_query_bridging_contigs'''
        hits = [
            '\t'.join(['15000', '20000', '1', '5000', '5000', '5000', '100.00', '20000', '10000', '1', '-1', 'ref1', 'qry1']),
            '\t'.join(['1', '5000', '5000', '10000', '5000', '5000', '100.00', '20000', '10000', '1', '-1', 'ref2', 'qry1']),
            '\t'.join(['1', '5000', '5000', '1', '5000', '5000', '100.00', '20000', '10000', '1', '-1', 'ref2', 'qry2']),
            '\t'.join(['15000', '20000', '10000', '5000', '5000', '5000', '100.00', '20000', '10000', '1', '-1', 'ref1', 'qry2']),
            '\t'.join(['15500', '20000', '10000', '5500', '4500', '4500', '100.00', '20000', '10000', '1', '-1', 'ref3', 'qry2']),
            '\t'.join(['15500', '20000', '10000', '5500', '4500', '4500', '100.00', '20000', '10000', '1', '-1', 'ref4', 'qry3']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        input_hits = {
            'qry1': [hits[0], hits[1]],
            'qry2': [hits[2], hits[3], hits[4]],
            'qry3': [hits[5]]
        }
        got = self.merger._get_possible_query_bridging_contigs(input_hits)
        expected = {
            'qry1': (hits[0], hits[1]),
            'qry2': (hits[2], hits[3])
        }
        self.assertEqual(got, expected)


    def test_filter_bridging_contigs(self):
        '''test _filter_bridging_contigs'''
        hits = [
            '\t'.join(['15000', '20000', '1', '5000', '5000', '5000', '100.00', '20000', '10000', '1', '-1', 'ref1', 'qry1']),
            '\t'.join(['1', '5000', '5000', '10000', '5000', '5000', '100.00', '20000', '10000', '1', '-1', 'ref2', 'qry1']),
            '\t'.join(['15000', '20000', '1', '5000', '5000', '5000', '100.00', '20000', '10000', '1', '-1', 'ref3', 'qry2']),
            '\t'.join(['1', '5000', '5000', '10000', '5000', '5000', '100.00', '20000', '10000', '1', '-1', 'ref4', 'qry2']),
            '\t'.join(['15000', '20000', '1', '5000', '5000', '5000', '100.00', '20000', '10000', '1', '-1', 'ref3', 'qry3']),
            '\t'.join(['1', '5000', '5000', '10000', '5000', '5000', '100.00', '20000', '10000', '1', '-1', 'ref5', 'qry3']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        input_hits = {
            'qry1': (hits[0], hits[1]),
            'qry2': (hits[2], hits[3]),
            'qry3': (hits[4], hits[5]),
        }
        got = self.merger._filter_bridging_contigs(input_hits)
        expected = {'qry1': (hits[0], hits[1])}
        self.assertEqual(got, expected)


    def test_merge_bridged_contig_pair1(self):
        '''test _merge_bridged_contig_pair1'''
        # simple case: ref contigs 1 and 2 do not overlap each other and are joined by a reassembly contig
        hits = [
            '\t'.join(['721', '999', '5', '283', '279', '279', '100.00', '1000', '753', '1', '1', 'ref1', 'reassembly']),
            '\t'.join(['1', '420', '324', '743', '420', '420', '100.00', '1000', '753', '1', '1', 'ref2', 'reassembly']),
        ]
        start_hit, end_hit = [pymummer.alignment.Alignment(x) for x in hits]

        ref_fasta = os.path.join(data_dir, 'merge_test_merge_bridged_contig_pair.test1.ref.fa')
        reassembly_fasta = os.path.join(data_dir, 'merge_test_merge_bridged_contig_pair.test1.reassembly.fa')
        expected_fasta = os.path.join(data_dir, 'merge_test_merge_bridged_contig_pair.test1.expected.fa')
        ref_contigs = {}
        pyfastaq.tasks.file_to_dict(ref_fasta, ref_contigs)
        reassembly_contigs = {}
        pyfastaq.tasks.file_to_dict(reassembly_fasta, reassembly_contigs)
        self.merger._merge_bridged_contig_pair(start_hit, end_hit, ref_contigs, reassembly_contigs)
        self.assertEqual(len(ref_contigs), 1)
        expected_ref_contigs = {}
        pyfastaq.tasks.file_to_dict(expected_fasta, expected_ref_contigs)
        self.assertEqual(list(expected_ref_contigs.values())[0].seq, list(ref_contigs.values())[0].seq)


    def test_merge_bridged_contig_pair2(self):
        '''test _merge_bridged_contig_pair2'''
        # simple case: ref contigs 1 and 2 do not overlap each other and are joined by a reassembly contig
        hits = [
            '\t'.join(['1', '420', '430', '11', '420', '420', '100.00', '1000', '753', '1', '-1', 'ref2', 'reassembly']),
            '\t'.join(['721', '999', '749', '471', '279', '279', '100.00', '1000', '753', '1', '-1', 'ref1', 'reassembly']),
        ]
        start_hit, end_hit = [pymummer.alignment.Alignment(x) for x in hits]

        ref_fasta = os.path.join(data_dir, 'merge_test_merge_bridged_contig_pair.test2.ref.fa')
        reassembly_fasta = os.path.join(data_dir, 'merge_test_merge_bridged_contig_pair.test2.reassembly.fa')
        expected_fasta = os.path.join(data_dir, 'merge_test_merge_bridged_contig_pair.test2.expected.fa')
        ref_contigs = {}
        pyfastaq.tasks.file_to_dict(ref_fasta, ref_contigs)
        reassembly_contigs = {}
        pyfastaq.tasks.file_to_dict(reassembly_fasta, reassembly_contigs)
        self.merger._merge_bridged_contig_pair(start_hit, end_hit, ref_contigs, reassembly_contigs)
        self.assertEqual(len(ref_contigs), 1)
        expected_ref_contigs = {}
        pyfastaq.tasks.file_to_dict(expected_fasta, expected_ref_contigs)
        self.assertEqual(list(expected_ref_contigs.values())[0].seq, list(ref_contigs.values())[0].seq)


    def test_merge_bridged_contig_pair3(self):
        '''test _merge_bridged_contig_pair3'''
        # simple case: ref contigs 1 and 2 do not overlap each other and are joined by a reassembly contig
        hits = [
            '\t'.join(['541', '960', '1', '420', '420', '420', '100.00', '960', '840', '1', '1', 'ref1', 'reassembly']),
            '\t'.join(['1', '480', '361', '840', '480', '480', '100.00', '1060', '840', '1', '1', 'ref2', 'reassembly']),
        ]
        start_hit, end_hit = [pymummer.alignment.Alignment(x) for x in hits]

        ref_fasta = os.path.join(data_dir, 'merge_test_merge_bridged_contig_pair.test3.ref.fa')
        reassembly_fasta = os.path.join(data_dir, 'merge_test_merge_bridged_contig_pair.test3.reassembly.fa')
        expected_fasta = os.path.join(data_dir, 'merge_test_merge_bridged_contig_pair.test3.expected.fa')
        ref_contigs = {}
        pyfastaq.tasks.file_to_dict(ref_fasta, ref_contigs)
        reassembly_contigs = {}
        pyfastaq.tasks.file_to_dict(reassembly_fasta, reassembly_contigs)
        self.merger._merge_bridged_contig_pair(start_hit, end_hit, ref_contigs, reassembly_contigs)
        self.assertEqual(len(ref_contigs), 1)
        expected_ref_contigs = {}
        pyfastaq.tasks.file_to_dict(expected_fasta, expected_ref_contigs)
        self.assertEqual(list(expected_ref_contigs.values())[0].seq, list(ref_contigs.values())[0].seq)


    def test_merge_bridged_contig_pair4(self):
        '''test _merge_bridged_contig_pair4'''
        # simple case: ref contigs 1 and 2 do not overlap each other and are joined by a reassembly contig
        hits = [
            '\t'.join(['1', '480', '480', '1', '480', '480', '100.00', '1060', '840', '1', '-1', 'ref2', 'reassembly']),
            '\t'.join(['541', '960', '840', '421', '420', '420', '100.00', '960', '840', '1', '-1', 'ref1', 'reassembly']),
        ]
        start_hit, end_hit = [pymummer.alignment.Alignment(x) for x in hits]

        ref_fasta = os.path.join(data_dir, 'merge_test_merge_bridged_contig_pair.test4.ref.fa')
        reassembly_fasta = os.path.join(data_dir, 'merge_test_merge_bridged_contig_pair.test4.reassembly.fa')
        expected_fasta = os.path.join(data_dir, 'merge_test_merge_bridged_contig_pair.test4.expected.fa')
        ref_contigs = {}
        pyfastaq.tasks.file_to_dict(ref_fasta, ref_contigs)
        reassembly_contigs = {}
        pyfastaq.tasks.file_to_dict(reassembly_fasta, reassembly_contigs)
        self.merger._merge_bridged_contig_pair(start_hit, end_hit, ref_contigs, reassembly_contigs)
        self.assertEqual(len(ref_contigs), 1)
        expected_ref_contigs = {}
        pyfastaq.tasks.file_to_dict(expected_fasta, expected_ref_contigs)
        self.assertEqual(list(expected_ref_contigs.values())[0].seq, list(ref_contigs.values())[0].seq)


    def test_index_fasta(self):
        '''test _index_fasta'''
        fasta_file = os.path.join(data_dir, 'merge_test_index_fasta.fa')
        expected_fai = os.path.join(data_dir, 'merge_test_index_fasta.fa.fai')
        test_fa = 'tmp.test_index_fasta.fa'
        test_fai = test_fa + '.fai'
        shutil.copyfile(fasta_file, test_fa)
        self.merger._index_fasta(test_fa)
        self.assertTrue(os.path.exists(test_fai))
        self.assertTrue(filecmp.cmp(expected_fai, test_fai, shallow=False))
        os.unlink(test_fa)
        os.unlink(test_fai)


    def test_write_act_files(self):
        '''test _write_act_files'''
        original_ref_fasta = os.path.join(data_dir, 'merge_test_write_act_files.ref.fa')
        original_qry_fasta = os.path.join(data_dir, 'merge_test_write_act_files.qry.fa')
        original_coords_file = os.path.join(data_dir, 'merge_test_write_act_files.coords')
        ref_fasta = 'tmp.test_write_act_files.ref.fa'
        qry_fasta = 'tmp.test_write_act_files.qry.fa'
        coords_file = 'tmp.test_write_act_files.coords'
        shutil.copyfile(original_ref_fasta, ref_fasta)
        shutil.copyfile(original_qry_fasta, qry_fasta)
        shutil.copyfile(original_coords_file, coords_file)

        expected_bash_script = os.path.join(data_dir, 'merge_test_write_act_files.expected.sh')
        expected_crunch_file = os.path.join(data_dir, 'merge_test_write_act_files.expected.crunch')
        outprefix = 'tmp.test_write_act_files.out'
        got_bash_script = outprefix + '.start_act.sh'
        got_crunch_file = outprefix + '.crunch'

        self.merger._write_act_files(ref_fasta, qry_fasta, coords_file, outprefix)
        self.assertTrue(filecmp.cmp(expected_bash_script, got_bash_script, shallow=False))
        self.assertTrue(filecmp.cmp(expected_crunch_file, got_crunch_file, shallow=False))
        for f in [ref_fasta, ref_fasta + '.fai', qry_fasta, qry_fasta + '.fai', coords_file, got_bash_script, got_crunch_file]:
            os.unlink(f)


    def test_contigs_dict_to_file(self):
        '''test _contigs_dict_to_file'''
        d = {
            '3': pyfastaq.sequences.Fasta('3', 'A'),
            '2': pyfastaq.sequences.Fasta('2', 'AAC'),
            '1': pyfastaq.sequences.Fasta('1', 'ACGTA'),
        }
        tmpfile = 'tmp.test_contigs_dict_to_file.fa'
        self.merger._contigs_dict_to_file(d, tmpfile)
        self.assertTrue(filecmp.cmp(tmpfile, os.path.join(data_dir, 'merge_test_contigs_dict_to_file.fa'), shallow=False))
        os.unlink(tmpfile)


    def test_make_new_contig_from_nucmer_and_spades_not_hit(self):
        '''test _make_new_contig_from_nucmer_and_spades no hit'''
        circular = set()
        hits = [
            '\t'.join(['1', '42', '2', '43', '42', '42', '100.0', '1000', '2000', '1', '-1', 'ref', 'reassembly'])
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        got = self.merger._make_new_contig_from_nucmer_and_spades('contig_name', hits, circular)
        self.assertEqual(got, (None, None))


    def test_make_new_contig_from_nucmer_and_spades_with_hit(self):
        '''test _make_new_contig_from_nucmer_and_spades with hit'''
        circular = {'spades_node'}
        spades_node = pyfastaq.sequences.Fasta('spades_node', 'ACGTACGTACG')
        expected = pyfastaq.sequences.Fasta('contig_name', spades_node.seq), 'spades_node'
        self.merger.reassembly_contigs = {'spades_node': spades_node}
        hits = [
            '\t'.join(['1', '10', '1', '11', '11', '11', '100.0', '30', '11', '1', '-1', 'original', 'spades_node']),
            '\t'.join(['21', '30', '1', '11', '11', '11', '100.0', '30', '11', '1', '-1', 'original', 'spades_node']),
            '\t'.join(['11', '20', '1', '11', '11', '11', '100.0', '30', '11', '1', '-1', 'original', 'spades_node']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        got = self.merger._make_new_contig_from_nucmer_and_spades('contig_name', hits, circular)
        self.assertEqual(got, expected)


    def test_get_spades_circular_nodes(self):
        fastg = os.path.join(data_dir, 'merge_test_get_spades_circular_nodes.fastg')
        got = self.merger._get_spades_circular_nodes(fastg)
        expected = set(['NODE_1_length_5_cov_42.42_ID_1'])
        self.assertEqual(expected, got)

