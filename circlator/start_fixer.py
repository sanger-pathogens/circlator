import os
import shutil
import tempfile
import pyfastaq
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



    @staticmethod
    def _max_length_from_fasta_file(infile):
        freader = pyfastaq.sequences.file_reader(infile)
        return max([len(x) for x in freader])

