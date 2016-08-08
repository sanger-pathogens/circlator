import os
import shutil
import tempfile
import pyfastaq
import pymummer
import circlator

class Error (Exception): pass

class StartFixer:
    def __init__(self,
        input_assembly_fa,
        outprefix,
        min_percent_identity=70,
        genes_fa=None,
        ignore=None,
        verbose=False,
    ):
        if genes_fa is None:
            d = os.path.abspath(os.path.dirname(circlator.__file__))
            self.genes_fa = os.path.join(d, 'data', 'dnaA.fasta')
            assert os.path.exists(self.genes_fa)
        else:
            if not os.path.exists(genes_fa):
                raise Error('Genes file "' + genes_fa + '" not found. Cannot continue')
            self.genes_fa = os.path.abspath(genes_fa)

        if not os.path.exists(input_assembly_fa):
            raise Error('Error! File not found: ' + input_assembly_fa)
        self.input_assembly = {}
        pyfastaq.tasks.file_to_dict(input_assembly_fa, self.input_assembly)

        self.min_percent_identity = min_percent_identity
        self.outprefix = os.path.abspath(outprefix)

        if ignore is None:
            self.ignore = set()
        else:
            with open(ignore) as f:
                self.ignore = {x.rstrip().split()[0] for x in f}


    @classmethod
    def _max_length_from_fasta_file(cls, infile):
        freader = pyfastaq.sequences.file_reader(infile)
        return max([len(x) for x in freader])


    @classmethod
    def _write_fasta_plus_circularized_ends(cls, contigs, outfile, end_length, ignore=None):
        if ignore is None:
            ignore = set()

        f = pyfastaq.utils.open_file_write(outfile)
        used_names = set(contigs.keys())

        for contig_name, contig in sorted(contigs.items()):
            if contig_name in ignore:
                continue

            print(contig, file=f)

            if len(contig) >= 2 * end_length:
                start_coord = end_length - 1
                end_coord = len(contig) - end_length
                new_name = contig.id + '__ends'
                assert new_name not in used_names
                used_names.add(new_name)
                new_contig = pyfastaq.sequences.Fasta(new_name, contig[end_coord:] + contig[:start_coord + 1])
                print(new_contig, file=f)

        pyfastaq.utils.close(f)

