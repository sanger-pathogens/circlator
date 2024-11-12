import unittest
from circlator import external_progs

class TestExternalProgs(unittest.TestCase):
    
    def test_canu_version(self):
            '''Test canu version'''
            self.assertEqual('1.6', self.check_regex_version_extraction('canu', """
Canu 1.6
            """ ))
    
    def test_spades_version(self):
            '''Test spades version'''
            self.assertEqual('3.11.0', self.check_regex_version_extraction('spades', """
SPAdes v3.11.0
            """ ))
            self.assertEqual('3.7.1', self.check_regex_version_extraction('spades', """
SPAdes v3.7.1
            """ ))
            self.assertEqual('3.5.0', self.check_regex_version_extraction('spades', """
SPAdes genome assembler v.3.5.0
            """ ))
    
    def test_prodigal_version(self):
            '''Test prodigal version'''
            self.assertEqual('2.60', self.check_regex_version_extraction('prodigal', """

Prodigal V2.60: October, 2011

            """ ))
    
    def test_bwa_version(self):
        '''Test bwa version'''
        self.assertEqual('0.7.10', self.check_regex_version_extraction('bwa', """

Program: bwa (alignment via Burrows-Wheeler transformation)
Version: 0.7.10-r789
Contact: Heng Li <lh3@sanger.ac.uk>

        """ ))
        self.assertEqual('0.7.12', self.check_regex_version_extraction('bwa', """

Program: bwa (alignment via Burrows-Wheeler transformation)
Version: 0.7.12-r1039
Contact: Heng Li <lh3@sanger.ac.uk>
         """ ))

    def test_nucmer_version(self):    
        '''Test nucmer version'''
        self.assertEqual('3.1', self.check_regex_version_extraction('nucmer', """
nucmer
NUCmer (NUCleotide MUMmer) version 3.1
        """ ))
        self.assertEqual('4.0.0', self.check_regex_version_extraction('nucmer', """
4.0.0beta1
        """ ))
    
    def test_samtools_version(self):
        '''Test samtools version'''
        self.assertEqual('1.6', self.check_regex_version_extraction('samtools', """

Program: samtools (Tools for alignments in the SAM format)
Version: 1.6 (using htslib 1.6)
""" ))

    def test_samtools_original_version(self):
        '''Test samtools original version'''
        self.assertEqual('0.1.19', self.check_regex_version_extraction('samtools', """
Program: samtools (Tools for alignments in the SAM format)
Version: 0.1.19-44428cd

Usage:   samtools <command> [options]""" ))
         
    def check_regex_version_extraction(self, prog,  raw_version_output ):
        cmd, regex = external_progs.prog_to_version_cmd[prog]
        raw_output_lines = raw_version_output.splitlines()
        for line in raw_output_lines:
            hits = regex.search(line)
            if hits:
                return str(hits.group(1))
        return None
        