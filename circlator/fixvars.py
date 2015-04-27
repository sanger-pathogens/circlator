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
                d[name].sort(reverse=True)

        pyfastaq.utils.close(f)
        return snps, indels


    def _remove_overlapping_indels(self, indels):
        for name in indels:
            to_delete = set()
            previous_end = None
            for i in range(len(indels[name]) - 1, -1, -1):
                pos, ref, alt = indels[name][i]
                if previous_end is None:
                    previous_end = pos
                else:
                    if pos <= previous_end:
                        to_delete.add(i)
                    else:
                        previous_end = pos + max(len(ref), len(alt))

            indels[name] = [indels[name][i] for i in range(len(indels[name])) if i not in to_delete]

        return indels


    def _remove_indels_overlapping_snps(self, indels, snps):
        for name in indels:
            if name in snps:
                snp_positions = set([x[0] for x in snps[name]])
                to_remove = set()
                for i in range(len(indels[name])):
                    pos, ref, alt = indels[name][i]
                    positions = set(range(pos, pos + max(len(ref), len(alt))))
                    if len(positions.intersection(snp_positions)):
                        to_remove.add(i)

                indels[name] = [indels[name][i] for i in range(len(indels[name])) if i not in to_remove]

        return indels


    def _fix_variants(self, snps, indels, infile, outfile):
        seq_reader = pyfastaq.sequences.file_reader(infile)
        f_out = pyfastaq.utils.open_file_write(outfile)
        for seq in seq_reader:
            if seq.id in snps or seq.id in indels:
                bases = list(seq.seq)
                for d in snps, indels:
                    if seq.id in d:
                        for variant in d[seq.id]:
                            bases = self._fix_variant(variant, bases)

                seq.seq = ''.join(bases)

            print(seq, file=f_out)

        pyfastaq.utils.close(f_out)


    def _fix_variant(self, variant, sequence):
        pos, ref, alt = variant
        assert pos < len(sequence)
        expected = ''.join(sequence[pos:pos+len(ref)])
        if expected != ref:
            raise Error('Expected sequence "' + expected + '" but got "' + ref + '" at position ' + str(pos + 1))
           
        assert ''.join(sequence[pos:pos+len(ref)]) == ref
        sequence = sequence[:pos] + list(alt) + sequence[pos + len(ref):]
        return sequence


    def run(self):
        vcf_file = self.outprefix + '.vcf'
        self._make_vcf(vcf_file)
        snps, indels = self._get_variants_from_vcf(vcf_file)
        indels = self._remove_overlapping_indels(indels)
        indels = self._remove_indels_overlapping_snps(indels, snps)
        self._fix_variants(snps, indels, self.fasta_in, self.outprefix + '.fasta')

