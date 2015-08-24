import unittest
import filecmp
import os
import pyfastaq
from circlator import dnaa

modules_dir = os.path.dirname(os.path.abspath(dnaa.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestDnaa(unittest.TestCase):
    def test_header_to_genus_species(self):
        '''test _header_to_genus_species'''
        downloader = dnaa.UniprotDownloader()
        good_header = '>sp|P03004|DNAA_ECOLI Chromosomal replication initiator protein DnaA OS=Escherichia coli (strain K12) GN=dnaA PE=1 SV=2'
        bad_header = '>sp|P03004|DNAA_ECOLI Chromosomal replication initiator protein DnaA coli (strain K12) GN=dnaA PE=1 SV=2',

        self.assertEqual(downloader._header_to_genus_species(good_header), ('Escherichia', 'coli'))
        self.assertEqual(downloader._header_to_genus_species(bad_header), None)


    def test_get_uniprot_url(self):
        '''test _get_uniprot_url'''
        tests = [
            ('dnaa', 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=no&query=dnaa&force=no&format=fasta'),
            ('spam eggs', 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=no&query=spam+eggs&force=no&format=fasta'),
            ('  blue yellow   green  ', 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=no&query=blue+yellow+green&force=no&format=fasta'),
        ]

        for search_term, expected in tests:
            downloader = dnaa.UniprotDownloader(uniprot_search=search_term)
            self.assertEqual(downloader._get_uniprot_url(), expected)


    def test_header_matches_regex(self):
        '''test _header_matches_regex'''
        tests = [
            ('dnaa', 'foo dnaa spam eggs', True, True),
            ('Dnaa', 'foo dnaa spam eggs', True, True),
            ('Dnaa', 'foo dnaa spam eggs', False, False),
            ('dnaa', 'foo DNAA spam eggs', True, True),
            ('dnaa', 'foo spam eggs', True, False),
            ('SPAM', 'foo spam eggs', True, True),
        ]

        for regex, sequence_name, ignorecase, expected in tests:
            downloader = dnaa.UniprotDownloader(header_regex=regex, header_regex_ignorecase=ignorecase)
            seq = pyfastaq.sequences.Fasta(sequence_name, 'ACGT')
            self.assertEqual(downloader._header_matches_regex(seq), expected)


    def test_check_sequence(self):
        '''test _check_sequence'''
        downloader = dnaa.UniprotDownloader(min_gene_length=3, max_gene_length=6)
        seen = set()
        tests = [
            (pyfastaq.sequences.Fasta('dnaa OS=genus1 species1 foo', 'MABC'), (True, None)),
            (pyfastaq.sequences.Fasta('Dnaa OS=genus1 species2 foo', 'MABC'), (True, None)),
            (pyfastaq.sequences.Fasta('DnaA OS=genus2 species1 foo', 'MABC'), (True, None)),
            (pyfastaq.sequences.Fasta('dnaA OS=genus2 species2 foo', 'MABC'), (True, None)),
            (pyfastaq.sequences.Fasta('DNAA OS=genus2 species3 foo', 'MABC'), (True, None)),
            (pyfastaq.sequences.Fasta('dnaa OS=genus3 species1 foo', 'MA'), (False, 'Too long or short')),
            (pyfastaq.sequences.Fasta('dnaa OS=genus3 species2 foo', 'MABCDEF'), (False, 'Too long or short')),
            (pyfastaq.sequences.Fasta('dnaa OS=genus3 species3 foo', 'XABC'), (False, 'Does not start with M')),
            (pyfastaq.sequences.Fasta('x OS=genus3 species4 foo', 'MABC'), (False, 'No match to regex in name')),
            (pyfastaq.sequences.Fasta('dnaa foo', 'MABC'), (False, 'Error getting genus species')),
            (pyfastaq.sequences.Fasta('dnaa OS=genus1 species1 foo', 'MABC'), (False, "Duplicate genus species ('genus1', 'species1')")),
        ]

        for fa, expected in tests:
            self.assertEqual(downloader._check_sequence(fa, seen), expected)


    def test_append_stop_to_seq(self):
        '''test _append_stop_to_seq'''
        downloader = dnaa.UniprotDownloader()
        fa = pyfastaq.sequences.Fasta('x', 'A')
        downloader._append_stop_to_seq(fa)
        self.assertEqual(fa.seq, 'A*')
        downloader._append_stop_to_seq(fa)
        self.assertEqual(fa.seq, 'A*')


    def test_reverse_translate(self):
        '''test _reverse_translate'''
        downloader = dnaa.UniprotDownloader()
        aa = pyfastaq.sequences.Fasta('x', 'W')
        nucl = pyfastaq.sequences.Fasta('x', 'TGG')
        downloader._reverse_translate(aa)
        self.assertEqual(aa, nucl)
