import os

class Error (Exception): pass

class VariantFixer:
    def __init__(self,
        fasta_in,
        bam,
        outprefix,
        verbose=False
    ):
        self.fasta_in = os.path.abspath(fasta_in)
        self.bam_in = os.path.abspath(bam_in)
        self.outprefix = os.path.abspath(outprefix)
        self.verbose=verbose


    def _make_vcf(self, outfile):
        pass


    def _get_variants_from_vcf(self, infile):
        pass


    def _fix_variants(self, variants, infile, outfile):
        pass


    def _fix_variant(self, variant, sequence):
        pass


    def run(self):
        vcf_file = self.outprefix + '.vcf'
        self._make_vcf(vcf_file)
        variants = self._get_variants_from_vcf(vcf_file)
        self._fix_variants(variants, self.fasta_in, self.outprefix + '.fasta')

