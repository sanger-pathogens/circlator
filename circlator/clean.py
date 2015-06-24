import os
import shutil
import tempfile
import pymummer
import pyfastaq

class Error (Exception): pass

class Cleaner:
    def __init__(self,
        infile,
        outprefix,
        min_contig_length=2000,
        min_contig_percent_match=95,
        nucmer_diagdiff=25,
        nucmer_min_id=95,
        nucmer_min_length=500,
        nucmer_breaklen=500,
        keepfile=None,
        verbose=False
    ):
        self.infile = os.path.abspath(infile)
        self.outprefix = os.path.abspath(outprefix)
        self.min_contig_length = min_contig_length
        self.min_contig_percent_match = min_contig_percent_match
        self.contigs_to_keep = self._get_contigs_to_keep(keepfile)
        self.nucmer_diagdiff = nucmer_diagdiff
        self.nucmer_min_id = nucmer_min_id
        self.nucmer_min_length = nucmer_min_length
        self.nucmer_breaklen = nucmer_breaklen
        self.verbose = verbose


    def _get_contigs_to_keep(self, filename):
        '''Returns a set of names from file called filename. If filename is None, returns an empty set'''
        if filename is None:
            return set()

        with open(filename) as f:
            return {line.rstrip() for line in f}


    def _remove_small_contigs(self, infile, outfile, keep=None):
        '''Writes a new file with small contigs removed.
           Returns lists of all names and names of removed contigs'''
        removed = set()
        all_names = set()
        if keep is None:
            keep = set()

        file_reader = pyfastaq.sequences.file_reader(infile)
        fout = pyfastaq.utils.open_file_write(outfile)

        for seq in file_reader:
            all_names.add(seq.id)
            if len(seq) >= self.min_contig_length or seq.id in keep:
                print(seq, file=fout)
            else:
                removed.add(seq.id)

        pyfastaq.utils.close(fout)
        return all_names, removed


    def _run_nucmer(self, infile, outfile):
        '''Run nucmer of assembly against itself'''
        n = pymummer.nucmer.Runner(
            infile,
            infile,
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
        '''Returns two dictionaries:
           1) name=>contig length.
           2) Second is dictionary of nucmer hits (ignoring self matches).
              contig name => list of hits'''
        hits = {}
        lengths = {}
        file_reader = pymummer.coords_file.reader(infile)
        for al in file_reader:
            if al.qry_name == al.ref_name:
                continue
            elif al.qry_name not in hits:
                hits[al.qry_name] = []
            hits[al.qry_name].append(al)
            lengths[al.qry_name] = al.qry_length
            lengths[al.ref_name] = al.ref_length
        return lengths, hits


    def _contains(self, hit):
        '''Returns True iff (the query contig is contained in the reference contig and
           the query contig is not flagged to be kept)'''
        return (
            hit.qry_name not in self.contigs_to_keep
            and hit.qry_name != hit.ref_name
            and (100 * hit.hit_length_qry / hit.qry_length >= self.min_contig_percent_match)
            and hit.percent_identity >= self.nucmer_min_id
        )


    def _containing_contigs(self, hits):
        '''Given a list of hits, all with same query,
           returns a set of the contigs containing that query'''
        return {hit.ref_name for hit in hits if self._contains(hit)}


    def _get_containing_contigs(self, hits_dict):
        '''Given dictionary of nucmer hits (made by self._load_nucmer_hits()), returns a dictionary.
           key=contig name. Value = set of contigs that contain the key.'''
        containing = {}

        for qry_name in hits_dict:
            d = self._containing_contigs(hits_dict[qry_name])
            if len(d):
                containing[qry_name] = d

        return containing


    def _get_all_containing(self, containing_contigs, name, exclude=None):
        '''containing_contigs is a dict:
             key=contig name. Value = set of contigs that contain the key.
           Returns alls contigs called "name" that contain that contig'''
        contains_name = set()
        if name in containing_contigs:
            for containing_contig in containing_contigs[name]:
                # if we have a contains b and b contains a, then this stops infinite recursion
                if containing_contig==exclude:
                    continue
                contains_name.add(containing_contig)
                new_names = self._get_all_containing(containing_contigs, containing_contig, exclude=name)
                new_names.discard(name)
                contains_name.update(new_names)
        return contains_name


    def _expand_containing_using_transitivity(self, containing_contigs):
        '''This uses a contined in b, and b contained in c to force a contained in c.
           Just in case a contained in c wasn't already found by nucmer'''
        for name in containing_contigs:
            containing_contigs[name] = self._get_all_containing(containing_contigs, name)
        return containing_contigs


    def _collapse_list_of_sets(self, sets):
        '''Input is a list of sets. Merges any intersecting sets in the list'''
        found = True
        while found:
            found = False
            to_intersect = None
            for i in range(len(sets)):
                for j in range(len(sets)):
                    if i == j:
                        continue
                    elif sets[i].intersection(sets[j]):
                        to_intersect = i, j
                        break

                if to_intersect is not None:
                    break

            if to_intersect is not None:
                found = True
                sets[i].update(sets[j])
                sets.pop(j)

        return sets


    def _get_identical_contigs(self, hits_dict):
        '''Input is a dict:
             key=contig name. Value = set of contigs that contain the key.
           Returns a list of sets of contigs that are equivalent'''
        equivalent_contigs = []

        for qry_name, containing in hits_dict.items():
            equivalent = set()
            for containing_name in containing:
                if containing_name in hits_dict and qry_name in hits_dict[containing_name]:
                    equivalent.add(containing_name)
                    equivalent.add(qry_name)

            if len(equivalent):
                equivalent_contigs.append(equivalent)
                equivalent_contigs = self._collapse_list_of_sets(equivalent_contigs)

        return equivalent_contigs


    def _longest_contig(self, contig_set, contig_lengths):
        '''Returns the name of the longest contig, from the set of names contig_set. contig_lengths
           is expected to be a dictionary of contig name => length.'''
        longest_name = None
        max_length = -1
        for name in contig_set:
            if contig_lengths[name] > max_length:
                longest_name = name
                max_length = contig_lengths[name]
        assert max_length != -1
        assert longest_name is not None
        return longest_name


    def _remove_identical_contigs(self, containing_contigs, contig_lengths):
        '''Input is dictionary of containing contigs made by self._expand_containing_using_transitivity().
           Removes redundant identical contigs, leaving one representative (the longest) of
           each set of identical contigs.
           Returns new version of dictionary, and a dictionary of contig name => contig it was replaced with'''
        identical_contigs = self._get_identical_contigs(containing_contigs)
        to_replace = {}  # contig name => name to replace it with
        for contig_set in identical_contigs:
            longest_contig = self._longest_contig(contig_set, contig_lengths)
            for name in contig_set - {longest_contig}:
                assert name not in to_replace
                to_replace[name] = longest_contig

        for name, replace_with in to_replace.items():
            if replace_with not in containing_contigs:
                containing_contigs[replace_with] = set()

            if name in containing_contigs:
                containing_contigs[replace_with].update(containing_contigs[name])
                del containing_contigs[name]

        to_delete = set()

        for name, names_set in containing_contigs.items():
            assert name not in to_replace
            new_set = {to_replace.get(x, x) for x in names_set}
            new_set.discard(name)
            if len(new_set) > 0:
                containing_contigs[name] = new_set
            else:
                to_delete.add(name)

        for name in to_delete:
            del containing_contigs[name]

        return containing_contigs, to_replace


    def _clean_contigs(self, infile, outfile, containing_contigs, replaced_contigs):
        file_reader = pyfastaq.sequences.file_reader(infile)
        fout = pyfastaq.utils.open_file_write(outfile)

        for seq in file_reader:
            if seq.id not in replaced_contigs and seq.id not in containing_contigs:
                print(seq, file=fout)

        pyfastaq.utils.close(fout)


    def _write_log(self, outfile, prefix, all_contigs, small_removed, containing_contigs, replaced_contigs):
        f = pyfastaq.utils.open_file_write(outfile)

        for name in sorted(all_contigs):
            if name in self.contigs_to_keep:
                print(prefix, name, 'user_kept', sep='\t', file=f)
            elif name in small_removed:
                print(prefix, name, 'small_removed', sep='\t', file=f)
            elif name in containing_contigs:
                print(prefix, name, 'contained in', '\t'.join(sorted(containing_contigs[name])), sep='\t', file=f)
            elif name in replaced_contigs:
                print(prefix, name, 'replaced with', replaced_contigs[name], sep='\t', file=f)
            else:
                print(prefix, name, 'kept', sep='\t', file=f)

        pyfastaq.utils.close(f)


    def run(self):
        removed_small_file = self.outprefix + '.remove_small.fa'
        names_all, names_small = self._remove_small_contigs(self.infile, removed_small_file, keep=self.contigs_to_keep)
        nucmer_coords_file = self.outprefix + '.coords'
        self._run_nucmer(removed_small_file, nucmer_coords_file)
        contig_lengths, nucmer_hits = self._load_nucmer_hits(nucmer_coords_file)
        containing_contigs = self._get_containing_contigs(nucmer_hits)
        if self.verbose and len(containing_contigs) > 0:
            print('\nContig\tContained in')
            for x in containing_contigs:
                print(x, containing_contigs[x])
            print()

        containing_contigs = self._expand_containing_using_transitivity(containing_contigs)
        containing_contigs, replaced_contigs = self._remove_identical_contigs(containing_contigs, contig_lengths)
        self._clean_contigs(removed_small_file, self.outprefix + '.fasta', containing_contigs, replaced_contigs)
        self._write_log(self.outprefix + '.log', '[clean]', names_all, names_small, containing_contigs, replaced_contigs)

