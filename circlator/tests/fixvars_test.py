import unittest
import filecmp
import os
from circlator import fixvars

modules_dir = os.path.dirname(os.path.abspath(fixvars.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


def load_key_info_from_vcf(infile):
    '''The header info in vcfs made in the tests is not predictable. This gets the snp/indel info that we need for later'''
    


class TestVariantFixer(unittest.TestCase):
    def test_make_vcf(self):
        '''test _make_vcf'''
        fa = os.path.join(data_dir, 'variant_fixer_test_make_vcf.ref.fa')
        bam = os.path.join(data_dir, 'variant_fixer_test_make_vcf.bam')
        tmp_vcf = 'tmp.test_make_vcf.vcf'
        vf = fixvars.VariantFixer(fa, bam, 'x')
        vf._make_vcf(tmp_vcf)
        with open(tmp_vcf) as f:
            got = [line.rstrip() for line in f if not line.startswith('#')]
        expected = [
            '1\t85\t.\tT\tA\t123\tPASS\tDP=3;VDB=0.0221621;SGB=-0.511536;MQ0F=0;DPR=0,3;AC=2;AN=2;DP4=0,0,3,0;MQ=60\tGT:PL:DV\t1/1:151,9,0:3',
            '1\t723\t.\tAGGG\tAGG\t52\tPASS\tINDEL;IDV=3;IMF=1;DP=3;VDB=0.0221621;SGB=-0.511536;MQ0F=0;DPR=0,3;AC=2;AN=2;DP4=0,0,3,0;MQ=60\tGT:PL:DV\t1/1:80,9,0:3'
        ]
        self.assertEqual(expected, got)
        os.unlink(tmp_vcf)


    def test_get_variants_from_vcf(self):
        '''test _get_variants_from_vcf'''
        fa = os.path.join(data_dir, 'variant_fixer_test_get_variants_from_vcf.ref.fa')
        bam = os.path.join(data_dir, 'variant_fixer_test_get_variants_from_vcf.bam')
        vcf = os.path.join(data_dir, 'variant_fixer_test_get_variants_from_vcf.vcf')
        vf = fixvars.VariantFixer(fa, bam, 'x')
        got_snps, got_indels = vf._get_variants_from_vcf(vcf)
        expected_snps = {
            '1': [(84, 'T', 'A'), (41, 'A', 'G')],
            '2': [(9, 'C', 'A')]
        }
        expected_indels = {'1': [(722, 'AGGG', 'AGG')]}
        self.assertEqual(got_snps, expected_snps)
        self.assertEqual(got_indels, expected_indels)
         
    
    def test_remove_overlapping_indels(self):
        '''test _remove_overlapping_indels'''
        fa = os.path.join(data_dir, 'variant_fixer_test_remove_overalpping_indels.ref.fa')
        bam = os.path.join(data_dir, 'variant_fixer_test_remove_overalpping_indels.bam')
        vf = fixvars.VariantFixer(fa, bam, 'x')
        indels = {
            'seq1': [(10, 'TG', 'T'), (1, 'GT', 'ACGT'), (0, 'AG', 'AGT')],
        }
        got_indels = vf._remove_overlapping_indels(indels)
        expected_indels = {
            'seq1': [(10, 'TG', 'T'), (0, 'AG', 'AGT')]
        }


    def test_fix_variants(self):
        '''test _fix_variants'''
        fa = os.path.join(data_dir, 'variant_fixer_test_fix_variants.in.fa')
        bam = os.path.join(data_dir, 'variant_fixer_test_fix_variants.bam')
        expected = os.path.join(data_dir, 'variant_fixer_test_fix_variants.out.fa')
        tmp_out = 'tmp.test_fix_variants.fa'
        vf = fixvars.VariantFixer(fa, bam, 'x')
        snps = {
            'seq1': [(60, 'T', 'A'), (0, 'T', 'G')],
            'seq2': [(3, 'A', 'T')]
        }
        indels = {
            'seq1': [(3, 'ATC', 'ATCC')]
        }
        vf._fix_variants(snps, indels, fa, tmp_out)
        self.assertTrue(filecmp.cmp(expected, tmp_out, shallow=False))
        os.unlink(tmp_out)


    def test_fix_variant(self):
        '''test _fix_variant'''
        fa = os.path.join(data_dir, 'variant_fixer_test_fix_variant.ref.fa')
        bam = os.path.join(data_dir, 'variant_fixer_test_fix_variant.bam')
        vf = fixvars.VariantFixer(fa, bam, 'x')
        seq = list('AGCAT')
        seq = vf._fix_variant((0, 'A', 'G'), seq)
        self.assertEqual(seq, list('GGCAT'))
        seq = vf._fix_variant((4, 'T', 'C'), seq)
        self.assertEqual(seq, list('GGCAC'))
        seq = vf._fix_variant((2, 'CA', 'CAG'), seq)
        self.assertEqual(seq, list('GGCAGC'))
        seq = vf._fix_variant((5, 'C', 'CT'), seq)
        self.assertEqual(seq, list('GGCAGCT'))

