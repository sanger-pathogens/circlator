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
             min_length_short=5000,
             min_length_long=100000,
             min_read_length=250,
             min_remove_end_distance=45000,
    ):
        self.bam = os.path.abspath(bam)
        if not os.path.exists(self.bam):
            raise Error('File not found:' + self.bam)

        self.min_read_length = min_read_length
        self.min_length_short = min_length_short
        self.min_length_long = min_length_long
        self.min_remove_end_distance = min_remove_end_distance
        assert self.min_length_short < self.min_length_long

        self.reads_fa = os.path.abspath(outprefix + '.fasta')
        self.log = os.path.abspath(outprefix + '.log')



    def _get_ref_lengths(self):
        '''Gets the length of each reference sequence from the header of the bam. Returns dict name => length'''
        sam_reader = pysam.Samfile(self.bam, "rb")
        return dict(zip(sam_reader.references, sam_reader.lengths))


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


    def _break_reads(self, contig, position, fout):
        '''Get all reads from contig, but breaks them all at given position (0-based) in the reference. Writes to fout. Currently pproximate where it breaks (ignores indels in the alignment)'''
        sam_reader = pysam.Samfile(self.bam, "rb")
        for read in sam_reader.fetch(contig):
            seqs = []
            if read.pos < position < read.reference_end - 1:
                split_point = position - read.pos
                if split_point - 1 >= self.min_read_length:
                    sequence = mapping.aligned_read_to_read(read, revcomp=False, ignore_quality=True).subseq(0, split_point)
                    sequence.id += '.left'
                    seqs.append(sequence)
                if read.query_length - split_point >= self.min_read_length:
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


    def _choose_region(self, contig, contig_length):
        return pyfastaq.intervals.Interval(self.min_remove_end_distance, contig_length - self.min_remove_end_distance)
         

    def run(self):
        ref_lengths = self._get_ref_lengths()
        assert len(ref_lengths) > 0
        f_log = pyfastaq.utils.open_file_write(self.log)
        f_fq = pyfastaq.utils.open_file_write(self.reads_fa)

        for contig in sorted(ref_lengths):
            if ref_lengths[contig] < self.min_length_short:
                self._all_reads_from_contig(contig, f_fq)
                print(contig, ref_lengths[contig], 'unchanged', sep='\t', file=f_log)
            elif self.min_length_short <= ref_lengths[contig] < self.min_length_long:
                break_pos =  int(ref_lengths[contig] / 2)
                print(contig, ref_lengths[contig], 'break at position ' + str(break_pos + 1), sep='\t', file=f_log)
                self._break_reads(contig, break_pos, f_fq)
            else:
                region_to_remove = self._choose_region(contig, ref_lengths[contig])
                if region_to_remove is not None:
                    print(contig, ref_lengths[contig], 'remove region ' + str(region_to_remove.start + 1) +  '-' + str(region_to_remove.end + 1), sep='\t', file=f_log)
                    self._exclude_region(contig, region_to_remove.start, region_to_remove.end, f_fq)
                else:
                    break_pos =  int(ref_lengths[contig] / 2)
                    print(contig, ref_lengths[contig], 'could not find region to remove - break at position ' + str(break_pos + 1), sep='\t', file=f_log)
                    self._break_reads(contig, int(ref_lengths[contig] / 2), f_fq)

        self._get_all_unmapped_reads(f_fq)
        pyfastaq.utils.close(f_fq)
        pyfastaq.utils.close(f_log)
