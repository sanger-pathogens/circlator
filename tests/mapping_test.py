import unittest
import copy
import filecmp
import os
import pysam
import pyfastaq
from circlator import mapping

modules_dir = os.path.dirname(os.path.abspath(mapping.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestMapping(unittest.TestCase):
    def test_bwa_index(self):
        '''test bwa_index'''
        ref = os.path.join(data_dir, 'mapping_test_bwa_index.fa')
        outprefix = 'tmp.mapping_test.bwa_index'
        mapping.bwa_index(ref, outprefix)
        expected_files = [outprefix + '.' + x for x in ['amb', 'ann', 'bwt', 'pac', 'sa']]

        for filename in expected_files:
            self.assertTrue(os.path.exists(filename))
            os.unlink(filename)


    def test_bwa_index_clean(self):
        '''test bwa_index_clean'''
        files_prefix = 'tmp.test_bwa_index_clean'
        test_files = [files_prefix + '.' + x for x in ['amb', 'ann', 'bwt', 'pac', 'sa']]
        for filename in test_files:
            with open(filename, 'w') as f:
                pass

            self.assertTrue(os.path.exists(filename))

        mapping.bwa_index_clean(files_prefix)

        for filename in test_files:
            self.assertFalse(os.path.exists(filename))


    def test_bwa_mem(self):
        '''test bwa_mem'''
        ref = os.path.join(data_dir, 'mapping_test_bwa_mem.ref.fa')
        reads = os.path.join(data_dir, 'mapping_test_bwa_mem.reads.fq')
        outfile = 'tmp.mapping_test_bwa_mem.bam'
        mapping.bwa_mem(ref, reads, outfile)
        self.assertTrue(os.path.exists(outfile))
        self.assertTrue(os.path.exists(outfile + '.bai'))

        expected_reads = [
            '1:2:49:172',
            '1:1:113:235',
            '1:3:205:315',
            '2:5:27:142',
            '2:4:101:224',
            '2:6:214:330',
        ]

        sam_reader = pysam.Samfile(outfile, "rb")
        got_reads = []

        for read in sam_reader.fetch():
            got_reads.append(read.qname)
            self.assertFalse(read.is_unmapped)

        self.assertEqual(expected_reads, got_reads)
        os.unlink(outfile)
        os.unlink(outfile + '.bai')


    def test_aligned_read_to_read(self):
        '''test aligned_read_to_read'''
        infile = os.path.join(data_dir, 'mapping_test_aligned_read_to_read.bam')
        sam_reader = pysam.Samfile(infile, "rb")
        aln1, aln2 = [x for x in sam_reader.fetch()]
        read1_fq = pyfastaq.sequences.Fastq('read1', 'TGTGTAACACTCCACCTCTGGTTCCCAGAGTTCGGTATCCGGCCGATACTTGAGGATAGC', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIGHFEDCBA')
        read1_fa = pyfastaq.sequences.Fasta('read1', 'TGTGTAACACTCCACCTCTGGTTCCCAGAGTTCGGTATCCGGCCGATACTTGAGGATAGC')
        self.assertEqual(read1_fq, mapping.aligned_read_to_read(aln1))
        self.assertEqual(read1_fq, mapping.aligned_read_to_read(aln1, revcomp=False))
        self.assertEqual(read1_fa, mapping.aligned_read_to_read(aln1, ignore_quality=True))

        read2 = pyfastaq.sequences.Fastq('read2', 'GATCGTCACGAAAGAACCAAGCCGGATCGTGGGAGGGGTACAACTCAGGTGAATTAACGT', 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGFEDC')
        read2_rev = copy.copy(read2)
        read2_rev.revcomp()
        self.assertEqual(read2, mapping.aligned_read_to_read(aln2))
        self.assertEqual(read2_rev, mapping.aligned_read_to_read(aln2, revcomp=False))

