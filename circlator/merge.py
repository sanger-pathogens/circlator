import os
import sys
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
          spades_kmers=None,
          spades_use_first_success=False,
          spades_careful=True,
          spades_only_assembler=True,
          ref_end_tolerance=15000,
          qry_end_tolerance=1000,
          verbose=False,
          threads=1,
          log_prefix='merge',
    ):
        if not os.path.exists(original_assembly):
            raise Error('File not found:' + original_assembly)

        self.original_fasta = original_assembly
        self.reassembly = circlator.assembly.Assembly(reassembly)
        self.reads = reads
        self.outprefix = outprefix
        self.nucmer_diagdiff = nucmer_diagdiff
        self.nucmer_min_id = nucmer_min_id
        self.nucmer_min_length = nucmer_min_length
        self.nucmer_min_length_for_merges = nucmer_min_length_for_merges
        self.nucmer_breaklen = nucmer_breaklen
        self.min_spades_circular_percent = min_spades_circular_percent
        self.spades_kmers = spades_kmers
        self.spades_use_first_success = spades_use_first_success
        self.spades_careful = spades_careful
        self.spades_only_assembler = spades_only_assembler
        self.ref_end_tolerance = ref_end_tolerance
        self.qry_end_tolerance = qry_end_tolerance
        self.verbose = verbose
        self.threads = threads
        self.log_prefix = log_prefix
        self.merges = []
        self.original_contigs = {}
        self.reassembly_contigs = self.reassembly.get_contigs()
        pyfastaq.tasks.file_to_dict(self.original_fasta, self.original_contigs)


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
            simplify=True,
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


    def _get_possible_circular_ref_contigs(self, nucmer_hits, log_fh=None, log_outprefix=None):
        '''Returns a dict ref name => tuple(hit at start, hit at end) for each ref sequence in the hash nucmer_hits (each value is a list of nucmer hits)'''
        writing_log_file = None not in [log_fh, log_outprefix]
        maybe_circular = {}
        all_nucmer_hits = []
        for l in nucmer_hits.values():
            all_nucmer_hits.extend(l)
        nucmer_hits_by_qry = self._hits_hashed_by_query(all_nucmer_hits)

        for ref_name, list_of_hits in nucmer_hits.items():
            if writing_log_file:
                print(log_outprefix, ref_name, 'Checking ' + str(len(list_of_hits)) + ' nucmer hits', sep='\t', file=log_fh)

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
                if writing_log_file:
                    print(log_outprefix, ref_name, 'potential pair of nucmer hits for circularization:', sep='\t', file=log_fh)
                    print(log_outprefix, ref_name, '', longest_start_hit, sep='\t', file=log_fh)
                    print(log_outprefix, ref_name, '', longest_end_hit, sep='\t', file=log_fh)

                shortest_hit_length = self._min_qry_hit_length([longest_start_hit, longest_end_hit])
                has_longer_hit = self._has_qry_hit_longer_than(
                    nucmer_hits_by_qry[longest_start_hit.qry_name],
                    shortest_hit_length,
                    hits_to_exclude={longest_start_hit, longest_end_hit}
                )

                if writing_log_file and has_longer_hit:
                    print(log_outprefix, ref_name, 'cannot use this pair because longer match was found', sep='\t', file=log_fh)

                can_circularise = self._can_circularise(longest_start_hit, longest_end_hit)

                if writing_log_file and not can_circularise:
                    print(log_outprefix, ref_name, 'cannot use this pair because positions/orientations of matches no good', sep='\t', file=log_fh)

                if (not has_longer_hit) and can_circularise:
                    print(log_outprefix, ref_name, 'can use this pair of hits', sep='\t', file=log_fh)
                    maybe_circular[ref_name] = (longest_start_hit, longest_end_hit)

        return maybe_circular


    def _remove_keys_from_dict_with_nonunique_values(self, d, log_fh=None, log_outprefix=None):
        '''Returns a new dictionary, with keys from input dict removed if their value was not unique'''
        value_counts = collections.Counter(d.values())
        new_d = {}
        writing_log_file = None not in [log_fh, log_outprefix]

        for key in d:
            if value_counts[d[key]] == 1:
                new_d[key] = d[key]
            elif writing_log_file:
                print(log_outprefix, 'Reject because non-unique:', d[key], sep='\t', file=log_fh)

        return new_d


    def _make_circularised_contig(self, ref_start_hit, ref_end_hit):
        '''Given a nucmer ref_start_hit and ref_end_hit, returns a new contig. Assumes that these hits can be used to circularise the reference contig of the hits using the query contig'''
        assert ref_start_hit.ref_name == ref_end_hit.ref_name
        assert ref_start_hit.qry_name == ref_end_hit.qry_name


        qry_name = ref_start_hit.qry_name
        ref_name = ref_start_hit.ref_name
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


    def _circularise_contigs(self, nucmer_hits):
        log_fh = pyfastaq.utils.open_file_write(self.outprefix + '.circularise_details.log')
        log_outprefix = '[merge circularise_details]'

        long_nucmer_hits = {}
        for name, hits in nucmer_hits.items():
            long_nucmer_hits[name] = [x for x in hits if x.hit_length_qry >= self.nucmer_min_length_for_merges]
        to_circularise_with_nucmer = self._get_possible_circular_ref_contigs(long_nucmer_hits, log_fh=log_fh, log_outprefix=log_outprefix)
        to_circularise_with_nucmer = self._remove_keys_from_dict_with_nonunique_values(to_circularise_with_nucmer, log_fh=log_fh, log_outprefix=log_outprefix)
        used_spades_contigs = set()
        called_as_circular_by_spades = self.reassembly.circular_contigs()

        if len(called_as_circular_by_spades):
            circular_string = ','.join(sorted(called_as_circular_by_spades))
        else:
            circular_string = 'None'
        print(log_outprefix, file=log_fh)
        print(log_outprefix, '\tSPAdes reassembly contigs that are circular: ', circular_string, sep='', file=log_fh)

        fate_of_contigs = {
            x: set() for x in [
                'repetitive_deleted',
                'circl_using_nucmer',
                'circl_using_spades',
            ]
        }

        for ref_name in self.original_contigs:
            print(log_outprefix, sep='\t', file=log_fh)

            if ref_name in nucmer_hits:
                print(log_outprefix, ref_name, 'Trying to circularize. Has nucmer hits to check...', sep='\t', file=log_fh)
                new_contig, spades_contig = self._make_new_contig_from_nucmer_and_spades(ref_name, nucmer_hits[ref_name], called_as_circular_by_spades, log_fh=log_fh, log_outprefix=log_outprefix)

                if new_contig is None:
                    print(log_outprefix, ref_name, 'Could not circularize using matches to SPAdes circular contigs', sep='\t', file=log_fh)

                    if ref_name in to_circularise_with_nucmer:
                        start_hit, end_hit = to_circularise_with_nucmer[ref_name]
                        assert start_hit.ref_name == end_hit.ref_name == ref_name
                        print(log_outprefix, ref_name, 'Circularizing using this pair of nucmer matches to SPAdes contig:', sep='\t', file=log_fh)
                        print(log_outprefix, '\t', ref_name, '\t\t', start_hit, sep='', file=log_fh)
                        print(log_outprefix, '\t', ref_name, '\t\t', end_hit, sep='', file=log_fh)
                        new_contig = self._make_circularised_contig(start_hit, end_hit)
                        fate_of_contigs['circl_using_nucmer'].add(ref_name)
                    else:
                        print(log_outprefix, ref_name, 'Cannot circularize: no suitable nucmer hits', sep='\t', file=log_fh)
                else:
                    assert new_contig.id == ref_name
                    assert spades_contig is not None
                    if spades_contig in used_spades_contigs:
                        print(log_outprefix, ref_name, 'Is circular, but duplicate sequence, so deleting it', sep='\t', file=log_fh)
                        fate_of_contigs['repetitive_deleted'].add(ref_name)
                    else:
                        self.original_contigs[new_contig.id] = new_contig
                        fate_of_contigs['circl_using_spades'].add(ref_name)
                        print(log_outprefix, ref_name, 'Circularized using matches to SPAdes circular contigs', sep='\t', file=log_fh)
                        used_spades_contigs.add(spades_contig)
            else:
                print(log_outprefix, ref_name, 'Cannot circularize: no nucmer hits', sep='\t', file=log_fh)
                new_contig = None


            if new_contig is None:
                print(log_outprefix, ref_name, 'Circularized: no', sep='\t', file=log_fh)
            else:
                assert new_contig.id == ref_name
                self.original_contigs[new_contig.id] = new_contig
                print(log_outprefix, ref_name, 'Circularized: yes', sep='\t', file=log_fh)

        pyfastaq.utils.close(log_fh)
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


    def _get_possible_query_bridging_contigs(self, nucmer_hits, log_fh=None, log_outprefix=None):
        '''Input is dict qry_name -> list of nucmer hits to that qry. Returns dict qry_name -> tuple(start hit, end hit)'''
        bridges = {}
        writing_log_file = None not in [log_fh, log_outprefix]

        for qry_name, hits_to_qry in nucmer_hits.items():
            if len(hits_to_qry) < 2:
                continue

            if writing_log_file:
                print(log_outprefix, '\t', qry_name, ': checking nucmer matches', sep='', file=log_fh)

            longest_start_hit = self._get_longest_hit_at_qry_start(hits_to_qry)
            longest_end_hit = self._get_longest_hit_at_qry_end(hits_to_qry)

            if (
                None in (longest_start_hit, longest_end_hit)
              or longest_start_hit.ref_name == longest_end_hit.ref_name
              or self._hits_have_same_reference(longest_start_hit, longest_end_hit)
            ):
                if writing_log_file:
                    print(log_outprefix, '\t', qry_name, ': no potential pairs of hits to merge contigs', sep='', file=log_fh)
                continue

            if writing_log_file:
                print(log_outprefix, '\t', qry_name, ': potential pair of hits to merge contigs...', sep='', file=log_fh)
                print(log_outprefix, '\t', qry_name, ':    ', longest_start_hit, sep='', file=log_fh)
                print(log_outprefix, '\t', qry_name, ':    ', longest_end_hit, sep='', file=log_fh)

            shortest_hit_length = self._min_qry_hit_length([longest_start_hit, longest_end_hit])
            has_longer_hit = self._has_qry_hit_longer_than(
                hits_to_qry,
                shortest_hit_length,
                hits_to_exclude={longest_start_hit, longest_end_hit}
            )

            if has_longer_hit and writing_log_file:
                print(log_outprefix, '\t', qry_name, ':    rejected - there is a longer hit to elsewhere', sep='', file=log_fh)

            orientation_ok = self._orientation_ok_to_bridge_contigs(longest_start_hit, longest_end_hit)

            if writing_log_file and not orientation_ok:
                print(log_outprefix, '\t', qry_name, ':    rejected - orientation/distance from ends not correct to make a merge', sep='', file=log_fh)

            if orientation_ok and not has_longer_hit:
                if writing_log_file:
                    print(log_outprefix, '\t', qry_name, ':    might be used - no longer hits elsewhere and orientation/distance to ends OK', sep='', file=log_fh)
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


    def _merge_bridged_contig_pair(self, start_hit, end_hit, ref_contigs, qry_contigs, log_fh=None, log_outprefix=None):
        writing_log_file = None not in [log_fh, log_outprefix]
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

        if writing_log_file:
            print(log_outprefix, '\tMerging contigs ', start_hit.ref_name, ' and ', end_hit.ref_name, ' using these two hits:', sep='', file=log_fh)
            print(log_outprefix, start_hit, sep='\t\t', file=log_fh)
            print(log_outprefix, end_hit, sep='\t\t', file=log_fh)

        new_id = start_contig.id + '.' + end_contig.id
        new_contig = pyfastaq.sequences.Fasta(new_id, start_seq + bridge_seq + end_seq)
        self.merges.append([new_id, start_contig.id, end_contig.id])
        ref_contigs[new_contig.id] = new_contig


    def _merge_all_bridged_contigs(self, nucmer_hits, ref_contigs, qry_contigs, log_fh=None, log_outprefix=None):
        '''Input is dict of nucmer_hits. Makes any possible contig merges.
           Returns True iff any merges were made'''
        writing_log_file = None not in [log_fh, log_outprefix]

        if len(nucmer_hits) == 0:
            if writing_log_file:
                print(log_outprefix, 'No nucmer hits, so will not make any merges', sep='\t', file=log_fh)
            return

        all_nucmer_hits = []
        for l in nucmer_hits.values():
            all_nucmer_hits.extend(l)
        nucmer_hits_by_qry = self._hits_hashed_by_query(all_nucmer_hits)
        bridges = self._get_possible_query_bridging_contigs(nucmer_hits_by_qry, log_fh=log_fh, log_outprefix=log_outprefix)
        if writing_log_file:
            print(log_outprefix, '\tPotential contigs to use for merging: ', ' '.join(sorted(bridges.keys())), sep='', file=log_fh)
        bridges = self._filter_bridging_contigs(bridges)
        if writing_log_file:
            print(log_outprefix, '\tContigs to use for merging after uniqueness filtering: ', ' '.join(sorted(bridges.keys())), sep='', file=log_fh)

        merged = set()
        made_a_join = False

        for qry_name, (start_hit, end_hit) in bridges.items():
            if start_hit.ref_name in merged or end_hit.ref_name in merged:
                continue
            self._merge_bridged_contig_pair(start_hit, end_hit, ref_contigs, qry_contigs, log_fh=log_fh, log_outprefix=log_outprefix)
            merged.add(start_hit.ref_name)
            merged.add(end_hit.ref_name)
            made_a_join = True

        if writing_log_file:
            print(log_outprefix, '\tMade at least one contig join: ', made_a_join, sep='', file=log_fh)

        return made_a_join


    def _index_fasta(self, infile):
        fai = infile + '.fai'
        if not os.path.exists(fai):
            samtools = circlator.external_progs.make_and_check_prog('samtools')
            circlator.common.syscall(samtools.exe() + ' faidx ' + infile, verbose=self.verbose)


    def _write_act_files(self, ref_fasta, qry_fasta, coords_file, outprefix):
        '''Writes crunch file and shell script to start up ACT, showing comparison of ref and qry'''
        if self.verbose:
            print('Making ACT files from', ref_fasta, qry_fasta, coords_file)
        ref_fasta = os.path.relpath(ref_fasta)
        qry_fasta = os.path.relpath(qry_fasta)
        coords_file = os.path.relpath(coords_file)
        outprefix = os.path.relpath(outprefix)
        self._index_fasta(ref_fasta)
        self._index_fasta(qry_fasta)
        crunch_file = outprefix + '.crunch'
        pymummer.coords_file.convert_to_msp_crunch(
            coords_file,
            crunch_file,
            ref_fai=ref_fasta + '.fai',
            qry_fai=qry_fasta + '.fai'
        )

        bash_script = outprefix + '.start_act.sh'
        with open(bash_script, 'w') as f:
            print('#!/usr/bin/env bash', file=f)
            print('act', ref_fasta, crunch_file, qry_fasta, file=f)

        pyfastaq.utils.syscall('chmod +x ' + bash_script)


    def _iterative_bridged_contig_pair_merge(self, outprefix):
        '''Iteratively merges contig pairs using bridging contigs from reassembly, until no more can be merged'''
        if self.reads is None:
            if self.verbose:
                print('Skipping iterative contig merging because no reads given (see --reads option)')
            return self.original_fasta, None, None

        log_file = outprefix + '.iterations.log'
        log_fh = pyfastaq.utils.open_file_write(log_file)
        genome_fasta = self.original_fasta
        nucmer_coords = outprefix + '.iter.1.coords'
        reads_to_map = self.reads
        act_prefix = None
        made_a_join = True
        iteration = 1

        while made_a_join:
            this_log_prefix = '[' + self.log_prefix + ' iterative_merge ' + str(iteration) + ']'
            print(this_log_prefix, '\tUsing nucmer matches from ', nucmer_coords, sep='', file=log_fh)
            self._run_nucmer(genome_fasta, self.reassembly.contigs_fasta, nucmer_coords)
            act_prefix = outprefix + '.iter.' + str(iteration)
            print(this_log_prefix, '\tYou can view the nucmer matches with ACT using: ./', act_prefix, '.start_act.sh', sep='', file=log_fh)
            self._write_act_files(genome_fasta, self.reassembly.contigs_fasta, nucmer_coords, act_prefix)
            nucmer_hits_by_ref = self._load_nucmer_hits(nucmer_coords)
            made_a_join = self._merge_all_bridged_contigs(nucmer_hits_by_ref, self.original_contigs, self.reassembly_contigs, log_fh, this_log_prefix)
            iteration += 1

            if made_a_join:
                print(this_log_prefix, '\tMade at least one merge. Remapping reads and reassembling',sep='', file=log_fh)
                nucmer_coords = outprefix + '.iter.' + str(iteration) + '.coords'
                genome_fasta = outprefix + '.iter.' + str(iteration) + '.merged.fasta'
                self._contigs_dict_to_file(self.original_contigs, genome_fasta)
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
                    careful=self.spades_careful,
                    only_assembler=self.spades_only_assembler,
                    verbose=self.verbose,
                    spades_kmers=self.spades_kmers,
                    spades_use_first_success=self.spades_use_first_success,
                )
                a.run()
                self.reassembly = circlator.assembly.Assembly(assembler_dir)
                self.reassembly_contigs = self.reassembly.get_contigs()
            elif iteration <= 2:
                print(this_log_prefix, '\tNo contig merges were made',sep='', file=log_fh)

        pyfastaq.utils.close(log_fh)
        return genome_fasta, nucmer_coords, act_prefix + '.start_act.sh'


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


    def _make_new_contig_from_nucmer_and_spades(self, original_contig, hits, circular_spades, log_fh=None, log_outprefix=None):
        '''Tries to make new circularised contig from contig called original_contig. hits = list of nucmer hits, all with ref=original contg. circular_spades=set of query contig names that spades says are circular'''
        writing_log_file = None not in [log_fh, log_outprefix]
        hits_to_circular_contigs = [x for x in hits if x.qry_name in circular_spades]
        if len(hits_to_circular_contigs) == 0:
            if writing_log_file:
                print(log_outprefix, original_contig, 'No matches to SPAdes circular contigs', sep='\t', file=log_fh)
            return None, None

        for hit in hits_to_circular_contigs:
            print(log_outprefix, original_contig, 'Checking hit:', hit, sep='\t', file=log_fh)
            percent_query_covered = 100 * (hit.hit_length_qry / hit.qry_length)

            if self.min_spades_circular_percent <= percent_query_covered:
                print(log_outprefix, '\t', original_contig, '\t\tHit is long enough. Percent of contig covered by hit is ', percent_query_covered, sep='', file=log_fh)
                # the spades contig hit is long enough, but now check that
                # the input contig is covered by hits from this spades contig
                hit_intervals = [x.ref_coords() for x in hits_to_circular_contigs if x.qry_name == hit.qry_name]

                if len(hit_intervals) > 0:
                    pyfastaq.intervals.merge_overlapping_in_list(hit_intervals)
                    percent_covered = 100 * pyfastaq.intervals.length_sum_from_list(hit_intervals) / hit.ref_length

                    if writing_log_file:
                        print(log_outprefix, '\t', original_contig, '\t\treference bases covered by spades contig:', ', '.join([str(x) for x in hit_intervals]), sep='', file=log_fh)
                        print(log_outprefix, '\t', original_contig, '\t\t   ... which is ', percent_covered, ' percent of ', hit.ref_length, ' bases', sep='', file=log_fh)

                    if self.min_spades_circular_percent <= percent_covered:
                        if writing_log_file:
                            print(log_outprefix, original_contig, '\tUsing hit to call as circular (enough bases covered)', sep='\t', file=log_fh)

                        return pyfastaq.sequences.Fasta(original_contig, self.reassembly_contigs[hit.qry_name].seq), hit.qry_name
                    elif writing_log_file:
                        print(log_outprefix, original_contig, '\tNot using hit to call as circular (not enough bases covered)', sep='\t', file=log_fh)
            else:
                print(log_outprefix, original_contig, '\tNot using hit to call as circular (hit too short)', sep='\t', file=log_fh)

        if writing_log_file:
            print(log_outprefix, original_contig, 'No suitable matches to SPAdes circular contigs', sep='\t', file=log_fh)

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
        self.original_fasta, nucmer_coords_file, act_script = self._iterative_bridged_contig_pair_merge(self.outprefix + '.merge')
        self._write_merge_log(self.outprefix + '.merge.log')
        nucmer_circularise_coords = os.path.abspath(self.outprefix + '.circularise.coords')

        if nucmer_coords_file is None:
            self._run_nucmer(self.original_fasta, self.reassembly.contigs_fasta, nucmer_circularise_coords)
            self._write_act_files(self.original_fasta, self.reassembly.contigs_fasta, nucmer_circularise_coords, self.outprefix + '.circularise')
        else:
            os.symlink(nucmer_coords_file, nucmer_circularise_coords)
            os.symlink(act_script, self.outprefix + '.circularise.start_act.sh')

        nucmer_hits = self._load_nucmer_hits(nucmer_circularise_coords)
        self._circularise_contigs(nucmer_hits)
