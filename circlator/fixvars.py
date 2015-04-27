import os
from circlator import external_progs, common


class Error (Exception): pass

class VariantFixer:
    def __init__(self,
        fasta_in,
        bam_in,
        outprefix,
        verbose=False
    ):
        self.fasta_in = os.path.abspath(fasta_in)
        self.bam_in = os.path.abspath(bam_in)
        self.outprefix = os.path.abspath(outprefix)
        self.verbose=verbose
        self.samtools = external_progs.make_and_check_prog('samtools', verbose=verbose)
        self.bcftools = external_progs.make_and_check_prog('bcftools', verbose=verbose)

        for f in [self.fasta_in, self.bam_in]:
            if not os.path.exists(f):
                raise Error('The following file was not found. Cannot continue.\n' + f)

    def _make_vcf(self, outfile):
        cmd = ' '.join([
            self.samtools.exe(), 'mpileup',
            '-A -B -u',
            '-t INFO/DPR,DV',
            '-f', self.fasta_in,
            '-v', self.bam_in,
            '|',
            self.bcftools.exe(), 'call -m -v',
            '|',
            self.bcftools.exe(), 'filter',
            '-i "MIN(DV/DP)>=0.5 & DP>2"',
            '-o', outfile
        ])
        common.syscall(cmd, verbose=self.verbose)


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

