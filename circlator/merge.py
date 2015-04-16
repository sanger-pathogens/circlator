import os
import copy
import shutil
import pymummer
import pyfastaq
import circlator

class Error (Exception): pass

class Merger:
    def __init__(
          self,
          original_assembly,
          reassembly,
          outprefix,
          reads=None,
          nucmer_min_id=99,
          nucmer_min_length=5000,
          nucmer_breaklen=500,
          ref_end_tolerance=20000,
          qry_end_tolerance=20000,
          verbose=False,
          threads=1,
    ):
        for f in [original_assembly, reassembly]:
            if not os.path.exists(f):
                raise Error('File not found:' + f)

        self.original_fasta = original_assembly
        self.reassembly_fasta = reassembly
        self.reads = reads
        self.outprefix = outprefix
        self.nucmer_min_id = nucmer_min_id
        self.nucmer_min_length = nucmer_min_length
        self.nucmer_breaklen = nucmer_breaklen
        self.ref_end_tolerance = ref_end_tolerance
        self.qry_end_tolerance = qry_end_tolerance
        self.verbose = verbose
        self.threads = threads
        self.original_contigs = {}
        self.reassembly_contigs = {}
        pyfastaq.tasks.file_to_dict(self.original_fasta, self.original_contigs)
        pyfastaq.tasks.file_to_dict(self.reassembly_fasta, self.reassembly_contigs)


    def _run_nucmer(self, ref, qry, outfile):
        '''Run nucmer of new assembly vs original assembly'''
        n = pymummer.nucmer.Runner(
            ref,
            qry,
            outfile,
            min_id=self.nucmer_min_id,
            min_length=self.nucmer_min_length,
            maxmatch=True,
            breaklen=self.nucmer_breaklen,
            verbose=self.verbose
        )
        n.run()


    def _load_nucmer_hits(self, infile):
        '''Returns dict ref name => list of nucmer hits from infile'''
        hits = {}
        file_reader = pymummer.coords_file.reader(infile)
        for al in file_reader:
            if al.ref_name not in hits:
                hits[al.ref_name] = []
            hits[al.ref_name].append(al)
        return hits


    def _hits_hashed_by_query(self, hits):
        '''Input: list of nucmer hits. Output: dictionary, keys are query names, values are lists of hits'''
        d = {}
        for hit in hits:
            if hit.qry_name not in d:
                d[hit.qry_name] = []
            d[hit.qry_name].append(hit)
        return d


    def _hits_to_new_seq(self, hits):
        '''Input hits = list of nucmer hits, all with same query and ref names. Tries to make a new circularised reference contig using the hits'''
        if len(hits) < 2:
            return None

        ref_hits_at_start = []
        ref_hits_at_end = []
        ref_name = hits[0].ref_name
        ref_length = len(self.original_contigs[ref_name])
        qry_name = hits[0].qry_name
        qry_length = len(self.reassembly_contigs[qry_name])
        ref_start_interval = pyfastaq.intervals.Interval(0, self.ref_end_tolerance - 1)
        ref_end_interval = pyfastaq.intervals.Interval(ref_length - self.ref_end_tolerance, ref_length - 1)
        qry_start_interval = pyfastaq.intervals.Interval(0, self.qry_end_tolerance - 1)
        qry_end_interval = pyfastaq.intervals.Interval(qry_length - self.qry_end_tolerance, qry_length - 1)
     
        hits_at_ref_start = [
            x for x in hits if \
                x.ref_coords().intersects(ref_start_interval) and (
                    ( x.qry_coords().intersects(qry_end_interval) and x.on_same_strand() ) or \
                    ( x.qry_coords().intersects(qry_start_interval) and not x.on_same_strand() )
                )
        ]

        hits_at_ref_end = [
            x for x in hits if \
                x.ref_coords().intersects(ref_end_interval) and (
                    ( x.qry_coords().intersects(qry_start_interval) and x.on_same_strand() ) or \
                    ( x.qry_coords().intersects(qry_end_interval) and not x.on_same_strand() )
                )
        ]

        if len(hits_at_ref_start) == 0 or len(hits_at_ref_end) == 0:
            return None

        ref_start_hit = self._get_longest_hit(hits_at_ref_start)
        ref_end_hit = self._get_longest_hit(hits_at_ref_end)

        if ref_start_hit.on_same_strand() != ref_end_hit.on_same_strand():
            return None

        ref_start_coords = ref_start_hit.ref_coords()
        ref_end_coords = ref_end_hit.ref_coords()

        if ref_start_coords.intersects(ref_end_coords):
            new_ctg = copy.copy(self.reassembly_contigs[qry_name])
            new_ctg.id = ref_name
            return new_ctg

        if ref_start_hit.on_same_strand():
            qry_start_coords = ref_end_hit.qry_coords()
            qry_end_coords = ref_start_hit.qry_coords()
            bases = self.original_contigs[ref_name][ref_start_coords.end+1:ref_end_coords.start] + \
                    self.reassembly_contigs[qry_name][qry_start_coords.start:qry_end_coords.end+1]
            return pyfastaq.sequences.Fasta(ref_name, bases)
        else:
            qry_start_coords = ref_start_hit.qry_coords()
            qry_end_coords = ref_end_hit.qry_coords()
            tmp_seq = pyfastaq.sequences.Fasta('x', self.reassembly_contigs[qry_name][qry_start_coords.start:qry_end_coords.end+1])
            tmp_seq.revcomp()
            return pyfastaq.sequences.Fasta(ref_name, self.original_contigs[ref_name][ref_start_coords.end+1:ref_end_coords.start] + tmp_seq.seq)

        return None


    def _get_longest_hit(self, hits):
        '''Returns the longest hit from a list of nucmer hits'''
        max_length = -1
        max_index = -1
        for i in range(len(hits)):
            if hits[i].hit_length_ref > max_length:
                max_index = i
                max_length = hits[i].hit_length_ref
        assert max_length != -1 and max_index != -1
        return hits[i]


    def _make_new_contig(self, original_contig, hits):
        '''Makes a new contig from the contig with name original_contig, using a list of nucmer hits all with original_contig as the ref_name'''
        hits_by_query = self._hits_hashed_by_query(hits)
        for qry_name, l in hits_by_query.items():
            new_contig = self._hits_to_new_seq(l)
            if new_contig is not None:
                return new_contig
            else:
                continue
 
        return None


    def _nucmer_hits_to_potential_join(self, hits, genome_contigs, reassembly_contigs):
        '''Given a list of numcer hits, all to the same query, returns a pair of nucmer hits that could be used to join two of the genome contigs'''
        if len(hits) < 2:
            return None

        reassembly_contig = reassembly_contigs[hits[0].qry_name]
        reassembly_contig_length = len(reassembly_contig)
        if reassembly_contig_length < self.nucmer_min_length:
            return None
        
        qry_start_interval = pyfastaq.intervals.Interval(0, self.qry_end_tolerance)
        qry_end_interval = pyfastaq.intervals.Interval(max(0, reassembly_contig_length - self.qry_end_tolerance), reassembly_contig_length)
        hits_at_start = []
        hits_at_end = []
        for hit in hits:
            ref_len = len(genome_contigs[hit.ref_name])
            ref_start_interval = pyfastaq.intervals.Interval(0, self.ref_end_tolerance)
            ref_end_interval = pyfastaq.intervals.Interval(max(0, ref_len - self.ref_end_tolerance), ref_len)

            if hit.qry_coords().intersects(qry_start_interval) and (
                   (hit.ref_coords().intersects(ref_end_interval) and hit.on_same_strand()) or \
                   (hit.ref_coords().intersects(ref_start_interval) and not hit.on_same_strand())
            ):
                hits_at_start.append(hit)

            if hit.qry_coords().intersects(qry_end_interval) and (
                   (hit.ref_coords().intersects(ref_start_interval) and hit.on_same_strand()) or \
                   (hit.ref_coords().intersects(ref_end_interval) and not hit.on_same_strand())
            ):
                hits_at_end.append(hit)

        if len(hits_at_start) == len(hits_at_end) == 1 and hits_at_start[0].ref_name != hits_at_end[0].ref_name:
            return hits_at_start[0], hits_at_end[0]
        

    def _merge_pair(self, hits, ref_contigs, reassembly_contigs):
        '''Merges two reference contigs together that are bridged by a reassembly contig. Hits between the contigs are in the list "hits" - there should be exactly two of them'''
        assert len(hits) == 2
        assert hits[0].qry_name == hits[1].qry_name
        assert hits[0].ref_name != hits[1].ref_name
        start_hit, end_hit = hits
        bridging_contig = reassembly_contigs.pop(start_hit.qry_name)
        start_contig = ref_contigs.pop(start_hit.ref_name)
        end_contig = ref_contigs.pop(end_hit.ref_name)
        bridge_seq = bridging_contig[start_hit.qry_coords().start:end_hit.qry_coords().end + 1]
        if start_hit.on_same_strand():
            start_seq = start_contig[0:start_hit.ref_coords().start]
        else:
            tmp_seq = start_contig.subseq(start_hit.ref_coords().end + 1, len(start_contig))
            tmp_seq.revcomp()
            start_seq = tmp_seq.seq

        if end_hit.on_same_strand():
            end_seq = end_contig[end_hit.ref_coords().end + 1:]
        else:
            tmp_seq = end_contig.subseq(0, end_hit.ref_coords().start)
            tmp_seq.revcomp()
            end_seq = tmp_seq.seq

        new_contig = pyfastaq.sequences.Fasta(start_contig.id + '.' + end_contig.id, start_seq + bridge_seq + end_seq)
        ref_contigs[new_contig.id] = new_contig


    def _merge_contig_pairs(self, outprefix):
        '''Iteratively merges contig pairs using ovelapping contigs from reassembly, until no more can be merged'''
        if self.reads is None:
            return self.original_fasta, self.reassembly_fasta

        nucmer_coords = outprefix + '.tmp.coords'
        genome_fasta = outprefix + '.merged.fa'
        reassembly_fasta = outprefix + '.reassembly.fa'
        bam = outprefix + '.bam'
        reads_prefix = outprefix + '.reads'
        self._contigs_dict_to_file(self.original_contigs, genome_fasta)
        self._contigs_dict_to_file(self.reassembly_contigs, reassembly_fasta)
        made_join = True

        while made_join:
            made_join = False
            self._run_nucmer(genome_fasta, reassembly_fasta, nucmer_coords)
            nucmer_hits_by_ref = self._load_nucmer_hits(nucmer_coords)
            all_hits = []
            for l in nucmer_hits_by_ref.values():
                all_hits.extend(l)
            nucmer_hits_by_qry = self._hits_hashed_by_query(all_hits)
            potential_join_pairs = {}

            for reassemble_contig in nucmer_hits_by_qry:
                hits = self._nucmer_hits_to_potential_join(
                    nucmer_hits_by_qry[reassemble_contig],
                    self.original_contigs,
                    self.reassembly_contigs)

                if hits is not None:
                    assert hits[0].qry_name == hits[1].qry_name
                    key = hits[0].qry_name
                    if key not in potential_join_pairs:
                        potential_join_pairs[key] = []
                    potential_join_pairs[key].append(hits)

            potential_join_pairs = {x:potential_join_pairs[x][0] for x in potential_join_pairs if len(potential_join_pairs[x]) == 1}

            if len(potential_join_pairs):
                ref_seq_counts = {}
                for qry in potential_join_pairs:
                    for hit in potential_join_pairs[qry]:
                        ref_seq_counts[hit.ref_name] = ref_seq_counts.get(hit.ref_name, 0) + 1

                for qry in potential_join_pairs:
                    hits = potential_join_pairs[qry]
                    assert len(hits) == 2
                    if ref_seq_counts[hits[0].ref_name] == ref_seq_counts[hits[1].ref_name] == 1:
                        self._merge_pair(hits, self.original_contigs, self.reassembly_contigs)
                        self._contigs_dict_to_file(self.original_contigs, genome_fasta)
                        if os.path.exists(reads_prefix + '.fasta'):
                           reads_to_map = reads_prefix + '.fasta'
                        else:
                           reads_to_map = self.reads
                        circlator.mapping.bwa_mem(
                          genome_fasta,
                          reads_to_map,
                          bam,
                          threads=self.threads,
                          verbose=self.verbose,
                        )
                        bam_filter = circlator.bamfilter.BamFilter(bam, reads_prefix)
                        bam_filter.run()
                        assembler_dir = outprefix + '.assembly'
                        a = circlator.assemble.Assembler(reads_prefix + '.fasta', assembler_dir, threads=self.threads, verbose=self.verbose)
                        a.run()
                        os.rename(os.path.join(assembler_dir, 'scaffolds.fasta'), reassembly_fasta)
                        shutil.rmtree(assembler_dir)
                        pyfastaq.tasks.file_to_dict(reassembly_fasta, self.reassembly_contigs)
                        made_join = True

        os.unlink(nucmer_coords)
        return genome_fasta, reassembly_fasta


    def _contigs_dict_to_file(self, contigs, fname):
        '''Writes dictionary of contigs to file'''
        f = pyfastaq.utils.open_file_write(fname)
        for contig in sorted(contigs, key=lambda x:len(contigs[x]), reverse=True):
            print(contigs[contig], file=f)
        pyfastaq.utils.close(f)



    def run(self):
        out_fasta = os.path.abspath(self.outprefix + '.fasta')
        out_log = os.path.abspath(self.outprefix + '.log')
        merge_pairs_dir = os.path.abspath(self.outprefix + '.merge_pairs')

        self.original_fasta, self.reassembly_fasta = self._merge_contig_pairs(merge_pairs_dir)

        nucmer_coords = os.path.abspath(self.outprefix + '.coords')
        self._run_nucmer(self.original_fasta, self.reassembly_fasta, nucmer_coords)
        nucmer_hits = self._load_nucmer_hits(nucmer_coords)
        fasta_fh = pyfastaq.utils.open_file_write(out_fasta)
        log_fh = pyfastaq.utils.open_file_write(out_log)

        for original_name, original_contig in sorted(self.original_contigs.items()):
            if original_name in nucmer_hits:
                new_contig = self._make_new_contig(original_name, nucmer_hits[original_name])
                circularised = 0

                if new_contig is not None:
                    original_contig = new_contig
                    circularised = 1
                    
                print(original_contig, file=fasta_fh)
                print(original_name, circularised, sep='\t', file=log_fh)
            else:
                print(original_contig, file=fasta_fh)
                print(original_name, 0, sep='\t', file=log_fh)
        
        pyfastaq.utils.close(fasta_fh)
        pyfastaq.utils.close(log_fh)

