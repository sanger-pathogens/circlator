import os
import pyfastaq

class Error (Exception): pass

class Assembly:
    def __init__(self, path):
        '''path can be a directory or a filename. If directory, assumes the name of a SPAdes
           output directory. If a file, assumes it is a fasta file of contigs'''
        if not os.path.exists(path):
            raise Error('Input path to Assembly.__init__ not found: ' + path)
        elif os.path.isdir(path):
            self.spades_dir = os.path.abspath(path)
        else:
            self.contigs_fasta = os.path.abspath(path)
            self.spades_dir = None

        self._set_filenames()


    @staticmethod
    def _file_exists(filename):
        if os.path.exists(filename):
            return filename
        else:
            return None


    def _set_filenames(self):
        self.contigs_fastg = None
        self.contigs_paths = None
        self.assembly_graph_fastg = None

        if self.spades_dir is None:
            return

        contigs_fasta = os.path.join(self.spades_dir, 'contigs.fasta')
        self.contigs_fasta = self._file_exists(os.path.join(contigs_fasta))

        if self.contigs_fasta is None:
            raise Error('Error finding SPAdes contigs file: ' + contigs_fasta)

        self.contigs_fastg = self._file_exists(os.path.join(self.spades_dir, 'contigs.fastg'))
        self.contigs_paths = self._file_exists(os.path.join(self.spades_dir, 'contigs.paths'))
        self.assembly_graph_fastg = self._file_exists(os.path.join(self.spades_dir, 'assembly_graph.fastg'))

        if None == self.contigs_fastg == self.contigs_paths == self.assembly_graph_fastg or \
           ( self.contigs_fastg is None and (None in {self.contigs_paths, self.assembly_graph_fastg}) ) or \
           ( self.contigs_fastg is not None and (self.contigs_paths is not None or self.assembly_graph_fastg is not None) ):
            error_message = '\n'.join([
                 'Error finding SPAdes graph files in the directory ' + self.spades_dir,
                 'Expected either:',
                 '    contigs.fastg (SPAdes <3.6.1)',
                 'or:',
                 '    contigs.paths and assembly_graph.fastg (SPAdes >3.6.1)'
            ])
            raise Error(error_message)


    def get_contigs(self):
        '''Returns a dictionary of contig_name -> pyfastaq.Sequences.Fasta object'''
        contigs = {}
        pyfastaq.tasks.file_to_dict(self.contigs_fasta, contigs)
        return contigs


    @classmethod
    def _circular_contigs_from_spades_before_3_6_1(cls, fastg_file):
        seq_reader = pyfastaq.sequences.file_reader(fastg_file)
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


    @classmethod
    def _spades_contigs_paths_to_dict(cls, filename):
        d = {}
        node = None

        with open(filename) as f:
            for line in f:
                if node is None:
                    if not line.startswith('NODE_'):
                        raise Error('Error loading info from SPAdes contigs path file ' + filename)
                    node = line.rstrip()
                else:
                    if line.startswith('NODE_') or node in d:
                        raise Error('Error loading info from SPAdes contigs path file ' + filename)
                    d[node] = line.rstrip()
                    node = None

        return d


    @classmethod
    def _circular_edges_to_edge_numbers_dict(cls, circular_edges):
        '''Edge names are like this: EDGE_1_length_44766_cov_16.555.
           This returns dict of edge number -> edge name.
           eg 1 -> EDGE_1_length_44766_cov_16.555.'''
        return {x.split('_')[1]: x for x in circular_edges}


    @staticmethod
    def _circular_contigs_from_spades_after_3_6_1(assembly_graph_fastg, contigs_paths):
        circular_graph_edges = Assembly._circular_contigs_from_spades_before_3_6_1(assembly_graph_fastg)
        circular_edge_dict = Assembly._circular_edges_to_edge_numbers_dict(circular_graph_edges)
        paths_dict = Assembly._spades_contigs_paths_to_dict(contigs_paths)
        circular_nodes = set()

        for node in paths_dict:
            if node.endswith("'"):
                continue

            edges = paths_dict[node].split(',')
            rev_node = node + "'"

            if len(edges) != 1 or rev_node not in paths_dict:
                continue

            rev_edges = paths_dict[rev_node].split(',')

            if len(rev_edges) != 1:
                continue

            edge = list(edges)[0][:-1]
            edge_strand = list(edges)[0][-1]
            rev_edge = list(rev_edges)[0][:-1]
            rev_edge_strand = list(rev_edges)[0][-1]
            if {'-', '+'} == {edge_strand, rev_edge_strand} and edge == rev_edge and edge in circular_edge_dict:
                circular_nodes.add(node)

        return circular_nodes


    def circular_contigs(self):
        '''Returns a set of the contig names that are circular'''
        if self.contigs_fastg is not None:
            return self._circular_contigs_from_spades_before_3_6_1(self.contigs_fastg)
        elif None not in [self.contigs_paths, self.assembly_graph_fastg]:
            return self._circular_contigs_from_spades_after_3_6_1(self.assembly_graph_fastg, self.contigs_paths)
        else:
            return set()
