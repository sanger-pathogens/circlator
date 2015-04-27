import os
import pyfastaq
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
        snps = {}
        indels = {}
        f = pyfastaq.utils.open_file_read(infile)
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.rstrip().split('\t')
            name = fields[0]
            pos = int(fields[1]) - 1
            ref = fields[3]
            alt = fields[4]

            if len(ref) == len(alt):
                d = snps
            else:
                d = indels

            if name not in d:
                d[name] = []

            d[name].append((pos, ref, alt))

        for d in snps, indels:
            for name in d:
                d[name].sort()

        pyfastaq.utils.close(f)
        return snps, indels


    def _fix_variants(self, variants, infile, outfile):
        pass


    def _fix_variant(self, variant, sequence):
        pos, ref, alt = variant
        assert pos < len(sequence)
        assert ''.join(sequence[pos:pos+len(ref)]) == ref
        sequence = sequence[:pos] + list(alt) + sequence[pos + len(ref):]
        return sequence


    def run(self):
        vcf_file = self.outprefix + '.vcf'
        self._make_vcf(vcf_file)
        variants = self._get_variants_from_vcf(vcf_file)
        self._fix_variants(variants, self.fasta_in, self.outprefix + '.fasta')

