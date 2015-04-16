import unittest
import filecmp
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


    def test_hits_to_new_seq1(self):
        '''test _hits_to_new_seq'''
        hits = [
            '\t'.join(['5', '360', '334', '689', '356', '356', '100.00', '817', '689', '1', '1', 'ref', 'reassembly']),
            '\t'.join(['481', '813', '1', '333', '333', '333', '100.00', '817', '689', '1', '1', 'ref', 'reassembly'])
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        self.merger.qry_end_tolerance = self.merger.ref_end_tolerance = 100
        ref_fasta = os.path.join(data_dir, 'merge_test_hits_to_new_seq.1.ref.fa')
        reassembly_fasta = os.path.join(data_dir, 'merge_test_hits_to_new_seq.1.reassembly.fa')
        expected_fasta = os.path.join(data_dir, 'merge_test_hits_to_new_seq.1.expected.fa')
        self.merger.original_contigs = {}
        self.merger.reassembly_contigs = {}
        pyfastaq.tasks.file_to_dict(ref_fasta, self.merger.original_contigs)
        pyfastaq.tasks.file_to_dict(reassembly_fasta, self.merger.reassembly_contigs)
        got = self.merger._hits_to_new_seq(hits)
        expected_ref_contigs = {}
        pyfastaq.tasks.file_to_dict(expected_fasta, expected_ref_contigs)
        self.assertEqual(list(expected_ref_contigs.values())[0].seq, got.seq)
        

    def test_hits_to_new_seq2(self):
        '''test _hits_to_new_seq 2'''
        # same as test_hits_to_new_seq1 but reassembly is reverse complemented
        hits = [
            '\t'.join(['5', '360', '356', '1', '356', '356', '100.00', '817', '689', '1', '-1', 'ref', 'reassembly']),
            '\t'.join(['481', '813', '689', '357', '333', '333', '100.00', '817', '689', '1', '-1', 'ref', 'reassembly']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        self.merger.qry_end_tolerance = self.merger.ref_end_tolerance = 100
        ref_fasta = os.path.join(data_dir, 'merge_test_hits_to_new_seq.2.ref.fa')
        reassembly_fasta = os.path.join(data_dir, 'merge_test_hits_to_new_seq.2.reassembly.fa')
        expected_fasta = os.path.join(data_dir, 'merge_test_hits_to_new_seq.2.expected.fa')
        self.merger.original_contigs = {}
        self.merger.reassembly_contigs = {}
        pyfastaq.tasks.file_to_dict(ref_fasta, self.merger.original_contigs)
        pyfastaq.tasks.file_to_dict(reassembly_fasta, self.merger.reassembly_contigs)
        got = self.merger._hits_to_new_seq(hits)
        expected_ref_contigs = {}
        pyfastaq.tasks.file_to_dict(expected_fasta, expected_ref_contigs)
        self.assertEqual(list(expected_ref_contigs.values())[0].seq, got.seq)


    def test_get_longest_hit(self):
        '''test _get_longest_hit'''
        hits = [
            '\t'.join(['61', '500', '61', '500', '440', '440', '100.00', '500', '500', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['10', '50', '11', '52', '51', '52', '99.42', '500', '500', '1', '1', 'ref1', 'qry1']),
            '\t'.join(['1', '500', '1', '499', '500', '499', '99.40', '500', '499', '1', '1', 'ref2', 'qry2', '[IDENTITY]'])
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        self.assertEqual(self.merger._get_longest_hit(hits), hits[2])


    def test_make_new_contig(self):
        '''test _make_new_contig'''
        # FIXME?
        pass


    def test_nucmer_hits_to_potential_join(self):
        '''test _nucmer_hits_to_potential_join'''
        hits = [
            '\t'.join(['721', '999', '5', '283', '279', '279', '100.00', '1000', '753', '1', '1', 'ref1', 'reassembly']),
            '\t'.join(['1', '420', '324', '743', '420', '420', '100.00', '1000', '753', '1', '1', 'ref2', 'reassembly']),
            '\t'.join(['200', '420', '324', '543', '420', '420', '100.00', '1000', '753', '1', '1', 'ref3', 'reassembly']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]
        self.merger.nucmer_min_length = 100
        self.merger.qry_end_tolerance = self.merger.ref_end_tolerance = 30
        
        ref_fasta = os.path.join(data_dir, 'merge_test_nucmer_hits_to_potential_join.ref.fa')
        reassembly_fasta = os.path.join(data_dir, 'merge_test_nucmer_hits_to_potential_join.reassembly.fa')
        ref_contigs = {}
        pyfastaq.tasks.file_to_dict(ref_fasta, ref_contigs)
        reassembly_contigs = {}
        pyfastaq.tasks.file_to_dict(reassembly_fasta, reassembly_contigs)
        got = self.merger._nucmer_hits_to_potential_join(hits, ref_contigs, reassembly_contigs)
        expected = hits[0], hits[1]
        self.assertEqual(got, expected)



    def test_merge_pair1(self):
        '''test _merge_pair 1'''
        # simple case: ref contigs 1 and 2 do not overlap each other and are joined by a reassembly contig
        hits = [
            '\t'.join(['721', '999', '5', '283', '279', '279', '100.00', '1000', '753', '1', '1', 'ref1', 'reassembly']),
            '\t'.join(['1', '420', '324', '743', '420', '420', '100.00', '1000', '753', '1', '1', 'ref2', 'reassembly']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]

        ref_fasta = os.path.join(data_dir, 'merge_test_merge_pair.test1.ref.fa')
        reassembly_fasta = os.path.join(data_dir, 'merge_test_merge_pair.test1.reassembly.fa')
        expected_fasta = os.path.join(data_dir, 'merge_test_merge_pair.test1.expected.fa')
        ref_contigs = {}
        pyfastaq.tasks.file_to_dict(ref_fasta, ref_contigs)
        reassembly_contigs = {}
        pyfastaq.tasks.file_to_dict(reassembly_fasta, reassembly_contigs)
        self.merger._merge_pair(hits, ref_contigs, reassembly_contigs)
        self.assertEqual(len(reassembly_contigs), 0)
        self.assertEqual(len(ref_contigs), 1)
        expected_ref_contigs = {}
        pyfastaq.tasks.file_to_dict(expected_fasta, expected_ref_contigs)
        self.assertEqual(list(expected_ref_contigs.values())[0].seq, list(ref_contigs.values())[0].seq)


    def test_merge_pair2(self):
        '''test _merge_pair 2'''
        # same as test_merge_pair1, except reassembly is reverse complemented
        hits = [
            '\t'.join(['1', '420', '430', '11', '420', '420', '100.00', '1000', '753', '1', '-1', 'ref2', 'reassembly']),
            '\t'.join(['721', '999', '749', '471', '279', '279', '100.00', '1000', '753', '1', '-1', 'ref1', 'reassembly']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]

        ref_fasta = os.path.join(data_dir, 'merge_test_merge_pair.test2.ref.fa')
        reassembly_fasta = os.path.join(data_dir, 'merge_test_merge_pair.test2.reassembly.fa')
        expected_fasta = os.path.join(data_dir, 'merge_test_merge_pair.test2.expected.fa')
        ref_contigs = {}
        pyfastaq.tasks.file_to_dict(ref_fasta, ref_contigs)
        reassembly_contigs = {}
        pyfastaq.tasks.file_to_dict(reassembly_fasta, reassembly_contigs)
        self.merger._merge_pair(hits, ref_contigs, reassembly_contigs)
        self.assertEqual(len(reassembly_contigs), 0)
        self.assertEqual(len(ref_contigs), 1)
        expected_ref_contigs = {}
        pyfastaq.tasks.file_to_dict(expected_fasta, expected_ref_contigs)
        self.assertEqual(list(expected_ref_contigs.values())[0].seq, list(ref_contigs.values())[0].seq)


    def test_merge_pair3(self):
        '''test _merge_pair 3'''
        #  the two reference contigs overlap
        hits = [
            '\t'.join(['541', '960', '1', '420', '420', '420', '100.00', '960', '840', '1', '1', 'ref1', 'reassembly']),
            '\t'.join(['1', '480', '361', '840', '480', '480', '100.00', '1060', '840', '1', '1', 'ref2', 'reassembly']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]

        ref_fasta = os.path.join(data_dir, 'merge_test_merge_pair.test3.ref.fa')
        reassembly_fasta = os.path.join(data_dir, 'merge_test_merge_pair.test3.reassembly.fa')
        expected_fasta = os.path.join(data_dir, 'merge_test_merge_pair.test3.expected.fa')
        ref_contigs = {}
        pyfastaq.tasks.file_to_dict(ref_fasta, ref_contigs)
        reassembly_contigs = {}
        pyfastaq.tasks.file_to_dict(reassembly_fasta, reassembly_contigs)
        self.merger._merge_pair(hits, ref_contigs, reassembly_contigs)
        self.assertEqual(len(reassembly_contigs), 0)
        self.assertEqual(len(ref_contigs), 1)
        expected_ref_contigs = {}
        pyfastaq.tasks.file_to_dict(expected_fasta, expected_ref_contigs)
        self.assertEqual(list(expected_ref_contigs.values())[0].seq, list(ref_contigs.values())[0].seq)


    def test_merge_pair4(self):
        '''test _merge_pair 4'''
        #  the two reference contigs overlap
        hits = [
            '\t'.join(['1', '480', '480', '1', '480', '480', '100.00', '1060', '840', '1', '-1', 'ref2', 'reassembly']),
            '\t'.join(['541', '960', '840', '421', '420', '420', '100.00', '960', '840', '1', '-1', 'ref1', 'reassembly']),
        ]
        hits = [pymummer.alignment.Alignment(x) for x in hits]

        ref_fasta = os.path.join(data_dir, 'merge_test_merge_pair.test4.ref.fa')
        reassembly_fasta = os.path.join(data_dir, 'merge_test_merge_pair.test4.reassembly.fa')
        expected_fasta = os.path.join(data_dir, 'merge_test_merge_pair.test4.expected.fa')
        ref_contigs = {}
        pyfastaq.tasks.file_to_dict(ref_fasta, ref_contigs)
        reassembly_contigs = {}
        pyfastaq.tasks.file_to_dict(reassembly_fasta, reassembly_contigs)
        self.merger._merge_pair(hits, ref_contigs, reassembly_contigs)
        self.assertEqual(len(reassembly_contigs), 0)
        self.assertEqual(len(ref_contigs), 1)
        expected_ref_contigs = {}
        pyfastaq.tasks.file_to_dict(expected_fasta, expected_ref_contigs)
        self.assertEqual(list(expected_ref_contigs.values())[0].seq, list(ref_contigs.values())[0].seq)


    def test_merge_contig_pairs(self):
        '''test _merge_contig_pairs'''
        # FIXME?
        pass


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
