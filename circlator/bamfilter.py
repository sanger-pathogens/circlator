import os
import pysam
import pyfastaq
from circlator import common, mapping

class Error (Exception): pass


class BamFilter:
    def __init__(
             self,
             bam,
             outprefix,
             length_cutoff=100000,
             min_read_length=250,
             contigs_to_use=None,
             discard_unmapped=False,
             log_prefix='[bamfilter]',
             verbose=False,
    ):
        self.bam = os.path.abspath(bam)
        if not os.path.exists(self.bam):
            raise Error('File not found:' + self.bam)

        self.length_cutoff = length_cutoff
        self.reads_fa = os.path.abspath(outprefix + '.fasta')
        self.log = os.path.abspath(outprefix + '.log')
        self.log_prefix = log_prefix
        self.contigs_to_use = self._get_contigs_to_use(contigs_to_use)
        self.discard_unmapped = discard_unmapped
        self.min_read_length = min_read_length
        self.verbose = verbose


    def _get_ref_lengths(self):
        '''Gets the length of each reference sequence from the header of the bam. Returns dict name => length'''
        sam_reader = pysam.Samfile(self.bam, "rb")
        return dict(zip(sam_reader.references, sam_reader.lengths))


    def _get_contigs_to_use(self, contigs_to_use):
        '''If contigs_to_use is a set, returns that set. If it's None, returns an empty set.
        Otherwise, assumes it's a file name, and gets names from the file'''
        if type(contigs_to_use) == set:
            return contigs_to_use
        elif contigs_to_use is None:
            return set()
        else:
            f = pyfastaq.utils.open_file_read(contigs_to_use)
            contigs_to_use = set([line.rstrip() for line in f])
            pyfastaq.utils.close(f)
            return contigs_to_use


    def _check_contigs_to_use(self, ref_dict):
        '''Checks that the set of contigs to use are all in the reference
        fasta lengths dict made by self._get_ref_lengths()'''
        if self.contigs_to_use is None:
            return True

        for contig in self.contigs_to_use:
            if contig not in ref_dict:
                raise Error('Requested to use contig "' + contig + '", but not found in input BAM file "' + self.bam + '"')

        return True


    def _all_reads_from_contig(self, contig, fout):
        '''Gets all reads from contig called "contig" and writes to fout'''
        sam_reader = pysam.Samfile(self.bam, "rb")
        for read in sam_reader.fetch(contig):
            print(mapping.aligned_read_to_read(read, ignore_quality=True), file=fout)


    def _get_all_unmapped_reads(self, fout):
        '''Writes all unmapped reads to fout'''
        sam_reader = pysam.Samfile(self.bam, "rb")
        for read in sam_reader.fetch(until_eof=True):
            if read.is_unmapped:
                print(mapping.aligned_read_to_read(read, ignore_quality=True), file=fout)


    def _break_reads(self, contig, position, fout, min_read_length=250):
        '''Get all reads from contig, but breaks them all at given position (0-based) in the reference. Writes to fout. Currently pproximate where it breaks (ignores indels in the alignment)'''
        sam_reader = pysam.Samfile(self.bam, "rb")
        for read in sam_reader.fetch(contig):
            seqs = []
            if read.pos < position < read.reference_end - 1:
                split_point = position - read.pos
                if split_point - 1 >= min_read_length:
                    sequence = mapping.aligned_read_to_read(read, revcomp=False, ignore_quality=True).subseq(0, split_point)
                    sequence.id += '.left'
                    seqs.append(sequence)
                if read.query_length - split_point >= min_read_length:
                    sequence = mapping.aligned_read_to_read(read, revcomp=False, ignore_quality=True).subseq(split_point, read.query_length)
                    sequence.id += '.right'
                    seqs.append(sequence)
            else:
                seqs.append(mapping.aligned_read_to_read(read, revcomp=False, ignore_quality=True))

            for seq in seqs:
                if read.is_reverse:
                    seq.revcomp()
                print(seq, file=fout)


    def _exclude_region(self, contig, start, end, fout):
        '''Writes reads not mapping to the given region of contig, start and end as per python convention'''
        sam_reader = pysam.Samfile(self.bam, "rb")
        exclude_interval = pyfastaq.intervals.Interval(start, end - 1)
        for read in sam_reader.fetch(contig):
            read_interval = pyfastaq.intervals.Interval(read.pos, read.reference_end - 1)
            if not read_interval.intersects(exclude_interval):
                print(mapping.aligned_read_to_read(read, ignore_quality=True), file=fout)


    def _get_region(self, contig, start, end, fout, min_length=250):
        '''Writes reads mapping to given region of contig, trimming part of read not in the region'''
        sam_reader = pysam.Samfile(self.bam, "rb")
        trimming_end = (start == 0)
        for read in sam_reader.fetch(contig, start, end):
            read_interval = pyfastaq.intervals.Interval(read.pos, read.reference_end - 1)
            seq = mapping.aligned_read_to_read(read, ignore_quality=True, revcomp=False)

            if trimming_end:
                bases_off_start = 0
                bases_off_end = max(0, read.reference_end - 1 - end)
                seq.seq = seq.seq[:read.query_alignment_end - bases_off_end]
            else:
                bases_off_start = max(0, start - read.pos + 1)
                seq.seq = seq.seq[bases_off_start  + read.query_alignment_start:]

            if read.is_reverse:
                seq.revcomp()

            if len(seq) >= min_length:
                print(seq, file=fout)


    def run(self):
        ref_lengths = self._get_ref_lengths()
        assert len(ref_lengths) > 0
        f_log = pyfastaq.utils.open_file_write(self.log)
        f_fa = pyfastaq.utils.open_file_write(self.reads_fa)
        print(self.log_prefix, '#contig', 'length', 'reads_kept', sep='\t', file=f_log)
        if self.verbose:
            print('Getting reads from BAM file', self.bam, flush=True)

        for contig in sorted(ref_lengths):
            if len(self.contigs_to_use) > 0 and contig not in self.contigs_to_use:
                print(self.log_prefix, contig, ref_lengths[contig], 'skipping', sep='\t', file=f_log)
                if self.verbose:
                    print('Skipping contig', contig, flush=True)
                continue

            if ref_lengths[contig] <= self.length_cutoff:
                self._all_reads_from_contig(contig, f_fa)
                print(self.log_prefix, contig, ref_lengths[contig], 'all', sep='\t', file=f_log)
                if self.verbose:
                    print('Getting all reads that map to contig', contig, flush=True)
            else:
                end_bases_keep = int(0.5 * self.length_cutoff)
                start = end_bases_keep - 1
                end = max(end_bases_keep - 1, ref_lengths[contig] - end_bases_keep)
                self._get_region(contig, 0, start, f_fa, min_length=self.min_read_length)
                self._get_region(contig, end, ref_lengths[contig], f_fa, min_length=self.min_read_length)
                coords_string = '1-' + str(start + 1) +  ';' + str(end + 1) + '-' + str(ref_lengths[contig])
                print(
                    self.log_prefix,
                    contig,
                    ref_lengths[contig],
                    coords_string,
                    sep='\t',
                    file=f_log
                )
                if self.verbose:
                    print('Getting all reads that map to ends of contig ', contig, '. Coords: ', coords_string, sep='', flush=True)

        if not self.discard_unmapped:
            if self.verbose:
                print('Getting all unmapped reads', sep='', flush=True)
            self._get_all_unmapped_reads(f_fa)

        pyfastaq.utils.close(f_fa)
        pyfastaq.utils.close(f_log)

        if self.verbose:
            print('Finished getting reads.')
            print('Log file:', self.log)
            print('Reads file:', self.reads_fa, flush=True)
