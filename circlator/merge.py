import os
import copy
import shutil
import collections
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
          nucmer_diagdiff=25,
          nucmer_min_id=95,
          nucmer_min_length=500,
          nucmer_min_length_for_merges=4000,
          nucmer_breaklen=500,
          min_spades_circular_percent=95,
          ref_end_tolerance=15000,
          qry_end_tolerance=1000,
          verbose=False,
          threads=1,
          log_prefix='merge',
    ):
        for f in [original_assembly, reassembly]:
            if not os.path.exists(f):
                raise Error('File not found:' + f)

        self.original_fasta = original_assembly
        self.reassembly_fasta = reassembly
        self.reads = reads
        self.outprefix = outprefix
        self.nucmer_diagdiff = nucmer_diagdiff
        self.nucmer_min_id = nucmer_min_id
        self.nucmer_min_length = nucmer_min_length
        self.nucmer_min_length_for_merges = nucmer_min_length_for_merges
        self.nucmer_breaklen = nucmer_breaklen
        self.min_spades_circular_percent = min_spades_circular_percent
        self.ref_end_tolerance = ref_end_tolerance
        self.qry_end_tolerance = qry_end_tolerance
        self.verbose = verbose
        self.threads = threads
        self.log_prefix = log_prefix
        self.merges = []
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
            diagdiff=self.nucmer_diagdiff,
            maxmatch=True,
            breaklen=self.nucmer_breaklen,
            simplify=False,
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


    def _get_longest_hit_by_ref_length(self, nucmer_hits):
        '''Input: list of nucmer hits. Returns the longest hit, taking hit length on the reference'''
        if len(nucmer_hits) == 0:
            return None

        max_length = None
        longest_hit = None
        for hit in nucmer_hits:
            if max_length is None or hit.hit_length_ref > max_length:
                max_length = hit.hit_length_ref
                longest_hit = copy.copy(hit)

        assert longest_hit is not None
        return longest_hit


    def _is_at_ref_start(self, nucmer_hit):
        '''Returns True iff the hit is "close enough" to the start of the reference sequence'''
        hit_coords = nucmer_hit.ref_coords()
        return hit_coords.start < self.ref_end_tolerance


    def _is_at_ref_end(self, nucmer_hit):
        '''Returns True iff the hit is "close enough" to the end of the reference sequence'''
        hit_coords = nucmer_hit.ref_coords()
        return hit_coords.end >= nucmer_hit.ref_length - self.ref_end_tolerance


    def _is_at_qry_start(self, nucmer_hit):
        '''Returns True iff the hit is "close enough" to the start of the query sequence'''
        hit_coords = nucmer_hit.qry_coords()
        return hit_coords.start < self.qry_end_tolerance


    def _is_at_qry_end(self, nucmer_hit):
        '''Returns True iff the hit is "close enough" to the end of the query sequence'''
        hit_coords = nucmer_hit.qry_coords()
        return hit_coords.end >= nucmer_hit.qry_length - self.qry_end_tolerance


    def _get_hit_nearest_ref_start(self, hits):
        '''Returns the hit nearest to the start of the ref sequence from the input list of hits'''
        nearest_to_start = hits[0]
        for hit in hits[1:]:
            if hit.ref_coords().start < nearest_to_start.ref_coords().start:
                nearest_to_start = hit
        return nearest_to_start


    def _get_hit_nearest_ref_end(self, hits):
        '''Returns the hit nearest to the end of the ref sequence from the input list of hits'''
        nearest_to_end = hits[0]
        for hit in hits[1:]:
            if hit.ref_coords().end > nearest_to_end.ref_coords().end:
                nearest_to_end = hit
        return nearest_to_end


    def _get_longest_hit_at_ref_start(self, nucmer_hits, hits_to_exclude=None):
        '''Input: list of nucmer hits to the same reference. Returns the longest hit to the start of the reference, or None if there is no such hit'''
        if hits_to_exclude is None:
            hits_to_exclude = set()
        hits_at_start = [hit for hit in nucmer_hits if self._is_at_ref_start(hit) and hit not in hits_to_exclude]
        return self._get_longest_hit_by_ref_length(hits_at_start)


    def _get_longest_hit_at_ref_end(self, nucmer_hits, hits_to_exclude=None):
        '''Input: list of nucmer hits to the same reference. Returns the longest hit to the end of the reference, or None if there is no such hit'''
        if hits_to_exclude is None:
            hits_to_exclude = set()
        hits_at_end = [hit for hit in nucmer_hits if self._is_at_ref_end(hit) and hit not in hits_to_exclude]
        return self._get_longest_hit_by_ref_length(hits_at_end)


    def _get_longest_hit_at_qry_start(self, nucmer_hits):
        '''Input: list of nucmer hits to the same query. Returns the longest hit to the start of the query, or None if there is no such hit'''
        hits_at_start = [hit for hit in nucmer_hits if self._is_at_qry_start(hit)]
        return self._get_longest_hit_by_ref_length(hits_at_start)


    def _get_longest_hit_at_qry_end(self, nucmer_hits):
        '''Input: list of nucmer hits to the same query. Returns the longest hit to the end of the query, or None if there is no such hit'''
        hits_at_end = [hit for hit in nucmer_hits if self._is_at_qry_end(hit)]
        return self._get_longest_hit_by_ref_length(hits_at_end)


    def _hits_have_same_query(self, nucmer_hit1, nucmer_hit2):
        '''Returns True iff the two nucmer hits are to the same query sequence'''
        return nucmer_hit1.qry_name == nucmer_hit2.qry_name


    def _hits_have_same_reference(self, nucmer_hit1, nucmer_hit2):
        '''Returns True iff the two nucmer hits are to the same reference sequence'''
        return nucmer_hit1.ref_name == nucmer_hit2.ref_name


    def _min_qry_hit_length(self, nucmer_hits):
        '''Returns the minimum query hit length from list of nucmer hits'''
        return min([hit.hit_length_qry for hit in nucmer_hits])


    def _has_qry_hit_longer_than(self, nucmer_hits, min_length, hits_to_exclude=None):
        '''Returns True iff list of nucmer_hits has a hit longer than min_length, not counting the hits in hits_to_exclude'''
        if hits_to_exclude is None:
            to_exclude = set()
        else:
            to_exclude = hits_to_exclude
        long_hits = [hit.hit_length_qry for hit in nucmer_hits if hit not in to_exclude and hit.hit_length_qry > min_length]
        return len(long_hits) > 0


    def _can_circularise(self, start_hit, end_hit):
        '''Returns true iff the two hits can be used to circularise the reference sequence of the hits'''
        if not(self._is_at_ref_start(start_hit) or self._is_at_ref_end(end_hit)):
            return False

        if self._is_at_qry_end(start_hit) \
          and self._is_at_qry_start(end_hit) \
          and start_hit.on_same_strand() \
          and end_hit.on_same_strand():
            return True

        if self._is_at_qry_start(start_hit) \
          and self._is_at_qry_end(end_hit) \
          and (not start_hit.on_same_strand()) \
          and (not end_hit.on_same_strand()):
            return True

        return False


    def _get_possible_circular_ref_contigs(self, nucmer_hits):
        '''Returns a dict ref name => tuple(hit at start, hit at end) for each ref sequence in the hash nucmer_hits (each value is a list of nucmer hits)'''
        maybe_circular = {}
        all_nucmer_hits = []
        for l in nucmer_hits.values():
            all_nucmer_hits.extend(l)
        nucmer_hits_by_qry = self._hits_hashed_by_query(all_nucmer_hits)

        for ref_name, list_of_hits in nucmer_hits.items():
            longest_start_hit = self._get_longest_hit_at_ref_start(list_of_hits)
            longest_end_hit = self._get_longest_hit_at_ref_end(list_of_hits)
            if longest_start_hit == longest_end_hit:
                second_longest_start_hit = self._get_longest_hit_at_ref_start(list_of_hits, hits_to_exclude={longest_start_hit})
                second_longest_end_hit = self._get_longest_hit_at_ref_end(list_of_hits, hits_to_exclude={longest_end_hit})
                if second_longest_start_hit is not None:
                    longest_start_hit = self._get_hit_nearest_ref_start([longest_start_hit, second_longest_start_hit])
                if second_longest_end_hit is not None:
                    longest_end_hit = self._get_hit_nearest_ref_end([longest_end_hit, second_longest_end_hit])

            if (
              longest_start_hit is not None
              and longest_end_hit is not None
              and longest_start_hit != longest_end_hit
              and self._hits_have_same_query(longest_start_hit, longest_end_hit)
            ):
                shortest_hit_length = self._min_qry_hit_length([longest_start_hit, longest_end_hit])
                has_longer_hit = self._has_qry_hit_longer_than(
                    nucmer_hits_by_qry[longest_start_hit.qry_name],
                    shortest_hit_length,
                    hits_to_exclude={longest_start_hit, longest_end_hit}
                )

                if (not has_longer_hit) and self._can_circularise(longest_start_hit, longest_end_hit):
                    maybe_circular[ref_name] = (longest_start_hit, longest_end_hit)

        return maybe_circular


    def _remove_keys_from_dict_with_nonunique_values(self, d):
        '''Returns a new dictionary, with keys from input dict removed if their value was not unique'''
        value_counts = collections.Counter(d.values())
        return {key:d[key] for key in d if value_counts[d[key]] == 1}


    def _make_circularised_contig(self, ref_start_hit, ref_end_hit):
        '''Given a nucmer ref_start_hit and ref_end_hit, returns a new contig. Assumes that these hits can be used to circularise the reference contig of the hits using the query contig'''
        assert ref_start_hit.ref_name == ref_end_hit.ref_name
        assert ref_start_hit.qry_name == ref_end_hit.qry_name


        qry_name = ref_start_hit.qry_name
        ref_name = ref_start_hit.ref_name
        ref_start_coords = ref_start_hit.ref_coords()
        ref_end_coords = ref_end_hit.ref_coords()

        if self.verbose:
            print('[' + self.log_prefix + '] Using the following two nucmer hits to circularise', ref_name)
            print('[' + self.log_prefix + ']    ', ref_start_hit)
            print('[' + self.log_prefix + ']    ', ref_end_hit)

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


    def _circularise_contigs(self, nucmer_hits):
        long_nucmer_hits = {}
        for name, hits in nucmer_hits.items():
            long_nucmer_hits[name] = [x for x in hits if x.hit_length_qry >= self.nucmer_min_length_for_merges]
        to_circularise_with_nucmer = self._get_possible_circular_ref_contigs(long_nucmer_hits)
        if self.verbose:
            for name, hits in to_circularise_with_nucmer.items():
                print('[' + self.log_prefix + ']', name, 'maybe circular:')
                print('[' + self.log_prefix + ']    ', hits[0])
                print('[' + self.log_prefix + ']    ', hits[1])

        to_circularise_with_nucmer = self._remove_keys_from_dict_with_nonunique_values(to_circularise_with_nucmer)
        used_spades_contigs = set()
        reassembly_fastg = self.reassembly_fasta[:-1] + 'g'
        called_as_circular_by_spades = self._get_spades_circular_nodes(reassembly_fastg)
        fate_of_contigs = {
            x: set() for x in [
                'repetitive_deleted',
                'circl_using_nucmer',
                'circl_using_spades',
            ]
        }

        for ref_name in self.original_contigs:
            if ref_name in nucmer_hits:
                new_contig, spades_contig = self._make_new_contig_from_nucmer_and_spades(ref_name, nucmer_hits[ref_name], called_as_circular_by_spades)

                if new_contig is None:
                    if ref_name in to_circularise_with_nucmer:
                        start_hit, end_hit = to_circularise_with_nucmer[ref_name]
                        assert start_hit.ref_name == end_hit.ref_name == ref_name
                        new_contig = self._make_circularised_contig(start_hit, end_hit)
                        fate_of_contigs['circl_using_nucmer'].add(ref_name)
                else:
                    assert new_contig.id == ref_name
                    assert spades_contig is not None
                    if spades_contig in used_spades_contigs:
                        if self.verbose:
                            print('[' + self.log_prefix + ']    ', ref_name, 'is circular, but duplicate sequence, so deleting it')
                        fate_of_contigs['repetitive_deleted'].add(ref_name)
                    else:
                        self.original_contigs[new_contig.id] = new_contig
                        fate_of_contigs['circl_using_spades'].add(ref_name)
                        used_spades_contigs.add(spades_contig)
            else:
                new_contig = None


            if new_contig is not None:
                assert new_contig.id == ref_name
                self.original_contigs[new_contig.id] = new_contig


        out_fasta = os.path.abspath(self.outprefix + '.fasta')
        fasta_fh = pyfastaq.utils.open_file_write(out_fasta)
        out_log = os.path.abspath(self.outprefix + '.circularise.log')
        log_fh = pyfastaq.utils.open_file_write(out_log)

        print(
            '[' + self.log_prefix + ' circularised]',
            '#Contig',
            'repetitive_deleted',
            'circl_using_nucmer',
            'circl_using_spades',
            'circularised',
            sep='\t', file=log_fh
        )

        for ref_name, contig in self.original_contigs.items():
            circl_using_nucmer = 1 if ref_name in fate_of_contigs['circl_using_nucmer'] else 0
            circl_using_spades = 1 if ref_name in fate_of_contigs['circl_using_spades'] else 0
            repetitive_deleted = 1 if ref_name in fate_of_contigs['repetitive_deleted'] else 0
            assert circl_using_nucmer * circl_using_spades == 0
            print(
                '[' + self.log_prefix + ' circularised]',
                ref_name,
                repetitive_deleted,
                circl_using_nucmer,
                circl_using_spades,
                circl_using_nucmer + circl_using_spades,
                sep='\t', file=log_fh
            )
            if repetitive_deleted == 0:
                print(contig, file=fasta_fh)

        for name in fate_of_contigs['repetitive_deleted']:
            del self.original_contigs[ref_name]

        pyfastaq.utils.close(log_fh)
        pyfastaq.utils.close(fasta_fh)


    def _orientation_ok_to_bridge_contigs(self, start_hit, end_hit):
        '''Returns True iff the orientation of the hits means that the query contig of both hits can bridge the reference contigs of the hits'''
        assert start_hit.qry_name == end_hit.qry_name
        if start_hit.ref_name == end_hit.ref_name:
            return False

        if (
            (self._is_at_ref_end(start_hit) and start_hit.on_same_strand())
            or (self._is_at_ref_start(start_hit) and not start_hit.on_same_strand())
        ):
            start_hit_ok = True
        else:
            start_hit_ok = False

        if (
            (self._is_at_ref_start(end_hit) and end_hit.on_same_strand())
            or (self._is_at_ref_end(end_hit) and not end_hit.on_same_strand())
        ):
            end_hit_ok = True
        else:
            end_hit_ok = False

        return start_hit_ok and end_hit_ok


    def _get_possible_query_bridging_contigs(self, nucmer_hits):
        '''Input is dict qry_name -> list of nucmer hits to that qry. Returns dict qry_name -> tuple(start hit, end hit)'''
        bridges = {}

        for qry_name, hits_to_qry in nucmer_hits.items():
            if len(hits_to_qry) < 2:
                continue

            longest_start_hit = self._get_longest_hit_at_qry_start(hits_to_qry)
            longest_end_hit = self._get_longest_hit_at_qry_end(hits_to_qry)

            if (
                None in (longest_start_hit, longest_end_hit)
              or longest_start_hit.ref_name == longest_end_hit.ref_name
              or self._hits_have_same_reference(longest_start_hit, longest_end_hit)
            ):
                continue

            shortest_hit_length = self._min_qry_hit_length([longest_start_hit, longest_end_hit])
            has_longer_hit = self._has_qry_hit_longer_than(
                hits_to_qry,
                shortest_hit_length,
                hits_to_exclude={longest_start_hit, longest_end_hit}
            )

            if self._orientation_ok_to_bridge_contigs(longest_start_hit, longest_end_hit) and not has_longer_hit:
                bridges[qry_name] = (longest_start_hit, longest_end_hit)

        return bridges


    def _filter_bridging_contigs(self, bridges):
        '''Input is dict qry_name -> tuple(start hit, end hit) made by _get_possible_query_bridging_contigs. Removes key/values where the value(==ref contig) has a hit at the start, or at the end, to more than one key(==qry contig)'''
        ref_hits_start = {}
        ref_hits_end = {}
        for qry_name, (start_hit, end_hit) in bridges.items():
            for hit in start_hit, end_hit:
                assert self._is_at_ref_start(hit) or self._is_at_ref_end(hit)
                if self._is_at_ref_start(hit):
                    ref_hits_start[hit.ref_name] = ref_hits_start.get(hit.ref_name, 0) + 1
                elif self._is_at_ref_end(hit):
                    ref_hits_end[hit.ref_name] = ref_hits_end.get(hit.ref_name, 0) + 1

        qry_names_to_remove = set()

        for qry_name, (start_hit, end_hit) in bridges.items():
            for hit in start_hit, end_hit:
                remove = False
                if (
                  (self._is_at_ref_start(hit) and ref_hits_start.get(hit.ref_name, 0) > 1)
                  or (self._is_at_ref_end(hit) and ref_hits_end.get(hit.ref_name, 0) > 1)
                ):
                    qry_names_to_remove.add(qry_name)

        for name in qry_names_to_remove:
            del bridges[name]

        return bridges


    def _merge_bridged_contig_pair(self, start_hit, end_hit, ref_contigs, qry_contigs):
        assert start_hit.qry_name == end_hit.qry_name
        assert start_hit.ref_name != end_hit.ref_name
        assert self._orientation_ok_to_bridge_contigs(start_hit, end_hit)
        bridging_contig = qry_contigs[start_hit.qry_name]
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

        if self.verbose:
            print('[' + self.log_prefix + '] Using the following two nucmer hits to merge contigs', start_hit.ref_name, end_hit.ref_name)
            print('[' + self.log_prefix + ']    ', start_hit)
            print('[' + self.log_prefix + ']    ', end_hit)
        new_id = start_contig.id + '.' + end_contig.id
        new_contig = pyfastaq.sequences.Fasta(new_id, start_seq + bridge_seq + end_seq)
        self.merges.append([new_id, start_contig.id, end_contig.id])
        ref_contigs[new_contig.id] = new_contig


    def _merge_all_bridged_contigs(self, nucmer_hits, ref_contigs, qry_contigs):
        '''Input is dict of nucmer_hits. Makes any possible contig merges.
           Returns True iff any merges were made'''
        if len(nucmer_hits) == 0:
            return
        all_nucmer_hits = []
        for l in nucmer_hits.values():
            all_nucmer_hits.extend(l)
        nucmer_hits_by_qry = self._hits_hashed_by_query(all_nucmer_hits)
        bridges = self._get_possible_query_bridging_contigs(nucmer_hits_by_qry)
        bridges = self._filter_bridging_contigs(bridges)
        merged = set()
        made_a_join = False

        for qry_name, (start_hit, end_hit) in bridges.items():
            if start_hit.ref_name in merged or end_hit.ref_name in merged:
                continue
            self._merge_bridged_contig_pair(start_hit, end_hit, ref_contigs, qry_contigs)
            merged.add(start_hit.ref_name)
            merged.add(end_hit.ref_name)
            made_a_join = True

        return made_a_join


    def _iterative_bridged_contig_pair_merge(self, outprefix):
        '''Iteratively merges contig pairs using bridging contigs from reassembly, until no more can be merged'''
        if self.reads is None:
            return self.original_fasta, self.reassembly_fasta, None

        genome_fasta = self.original_fasta
        reassembly_fasta = self.reassembly_fasta
        nucmer_coords = outprefix + '.iter.0.coords'
        reads_to_map = self.reads
        made_a_join = True
        iteration = 0

        while made_a_join:
            iteration += 1
            self._run_nucmer(genome_fasta, reassembly_fasta, nucmer_coords)
            nucmer_hits_by_ref = self._load_nucmer_hits(nucmer_coords)
            made_a_join = self._merge_all_bridged_contigs(nucmer_hits_by_ref, self.original_contigs, self.reassembly_contigs)

            if made_a_join:
                if self.verbose:
                    print('[' + self.log_prefix + '] Contig merges were made. Running new read mapping and assembly')

                nucmer_coords = outprefix + '.iter.' + str(iteration) + '.coords'
                genome_fasta = outprefix + '.iter.' + str(iteration) + '.merged.fasta'
                reassembly_fasta = outprefix + '.iter.' + str(iteration) + '.reassembly.fasta'
                reassembly_fastg = outprefix + '.iter.' + str(iteration) + '.reassembly.fastg'
                self._contigs_dict_to_file(self.original_contigs, genome_fasta)
                self._contigs_dict_to_file(self.reassembly_contigs, reassembly_fasta)
                bam = outprefix + '.iter.' + str(iteration) + '.bam'

                circlator.mapping.bwa_mem(
                  genome_fasta,
                  reads_to_map,
                  bam,
                  threads=self.threads,
                  verbose=self.verbose,
                )

                reads_prefix = outprefix + '.iter.' + str(iteration) + '.reads'
                reads_to_map =  reads_prefix + '.fasta'
                bam_filter = circlator.bamfilter.BamFilter(bam, reads_prefix)
                bam_filter.run()
                assembler_dir = outprefix + '.iter.' + str(iteration) + '.assembly'
                a = circlator.assemble.Assembler(
                    reads_prefix + '.fasta',
                    assembler_dir,
                    threads=self.threads,
                    verbose=self.verbose
                )
                a.run()
                os.rename(os.path.join(assembler_dir, 'contigs.fasta'), reassembly_fasta)
                os.rename(os.path.join(assembler_dir, 'contigs.fastg'), reassembly_fastg)
                shutil.rmtree(assembler_dir)
                pyfastaq.tasks.file_to_dict(reassembly_fasta, self.reassembly_contigs)
            elif iteration == 1 and self.verbose:
                print('[' + self.log_prefix + '] No contig merges were made')

        return genome_fasta, reassembly_fasta, nucmer_coords


    def _contigs_dict_to_file(self, contigs, fname):
        '''Writes dictionary of contigs to file'''
        f = pyfastaq.utils.open_file_write(fname)
        for contig in sorted(contigs, key=lambda x:len(contigs[x]), reverse=True):
            print(contigs[contig], file=f)
        pyfastaq.utils.close(f)


    def _get_spades_circular_nodes(self, fastg):
        '''Returns set of names of nodes in SPAdes fastg file that are circular. Names will match those in spades fasta file'''
        seq_reader = pyfastaq.sequences.file_reader(fastg)
        names = set([x.id.rstrip(';') for x in seq_reader if ':' in x.id])
        found_fwd = set()
        found_rev = set()
        for name in names:
            l = name.split(':')
            if len(l) != 2:
                continue
            if l[0] == l[1]:
                if l[0][-1] == "'":
                    found_rev.add(l[0][:-1])
                else:
                    found_fwd.add(l[0])

        return found_fwd.intersection(found_rev)


    def _make_new_contig_from_nucmer_and_spades(self, original_contig, hits, circular_spades):
        '''Tries to make new circularised contig from contig called original_contig. hits = list of nucmer hits, all with ref=original contg. circular_spades=set of query contig names that spades says are circular'''
        hits_to_circular_contigs = [x for x in hits if x.qry_name in circular_spades]
        if len(hits_to_circular_contigs) == 0:
            return None, None

        if self.verbose:
            print('[' + self.log_prefix + '] Hits to spades circular contig', hits[0].qry_name )
            for hit in hits_to_circular_contigs:
                print('[' + self.log_prefix + ']    ', hit)

        for hit in hits_to_circular_contigs:
            percent_query_covered = 100 * (hit.hit_length_qry / hit.qry_length)
            if self.verbose:
                print('[' + self.log_prefix + ']     percent query covered:', percent_query_covered, '...', hit)

            if self.min_spades_circular_percent <= percent_query_covered:
                # the spades contig hit is long enough, but now check that
                # the input contig is covered by hits from this spades contig
                hit_intervals = [x.ref_coords() for x in hits_to_circular_contigs if x.qry_name == hit.qry_name]

                if len(hit_intervals) > 0:
                    pyfastaq.intervals.merge_overlapping_in_list(hit_intervals)
                    percent_covered = 100 * pyfastaq.intervals.length_sum_from_list(hit_intervals) / hit.ref_length
                    if self.verbose:
                        print('[' + self.log_prefix + ']       reference bases covered by spades contig:', ', '.join([str(x) for x in hit_intervals]))
                        print('[' + self.log_prefix + ']      ...which is', percent_covered, 'percent of', hit.ref_length, 'bases')
                    if self.min_spades_circular_percent <= 100 * pyfastaq.intervals.length_sum_from_list(hit_intervals) / hit.ref_length:
                        if self.verbose:
                            print('[' + self.log_prefix + ']    ', hit.ref_name, 'is circular')
                        return pyfastaq.sequences.Fasta(original_contig, self.reassembly_contigs[hit.qry_name].seq), hit.qry_name

        return None, None


    def _write_merge_log(self, filename):
        if len(self.merges):
            log_fh = pyfastaq.utils.open_file_write(filename)
            print('[' + self.log_prefix + ' contig_merge]', '#new_name', 'previous_contig1', 'previous_contig2', sep='\t', file=log_fh)
            for l in self.merges:
                print('[' + self.log_prefix + ' contig_merge]', '\t'.join(l), sep='\t', file=log_fh)
            pyfastaq.utils.close(log_fh)


    def run(self):
        out_log = os.path.abspath(self.outprefix + '.merge.log')
        self.original_fasta, self.reassembly_fasta, nucmer_coords_file = self._iterative_bridged_contig_pair_merge(self.outprefix + '.merge')
        self._write_merge_log(self.outprefix + '.merge.log')
        nucmer_coords = os.path.abspath(self.outprefix + '.coords')

        if nucmer_coords_file is None:
            self._run_nucmer(self.original_fasta, self.reassembly_fasta, nucmer_coords)
        else:
            os.symlink(nucmer_coords_file, nucmer_coords)

        nucmer_hits = self._load_nucmer_hits(nucmer_coords)
        self._circularise_contigs(nucmer_hits)
