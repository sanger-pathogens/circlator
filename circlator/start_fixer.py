import os
import shutil
import copy
import tempfile
from distutils.version import LooseVersion
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

        self.contigs_fa_with_ends = self.outprefix + '.contigs_plus_ends.fa'


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


    @classmethod
    def _find_circular_using_promer(cls, outprefix, ref_genes_fa, contigs_dict, min_percent_id, end_length, log_fh, ignore=None):
        if ignore is None:
            ignore = set()
        promer_out = outprefix + '.promer'
        contigs_with_ends = outprefix + '.contigs_with_ends.fa'
        StartFixer._write_fasta_plus_circularized_ends(contigs_dict, contigs_with_ends, end_length, ignore=ignore)

        prunner = pymummer.nucmer.Runner(
            contigs_with_ends,
            ref_genes_fa,
            promer_out,
            min_id=min_percent_id,
            promer=True,
            verbose=True,
        )
        prunner.run()

        circularized = {} # original_contig_name -> promer match
        file_reader = pymummer.coords_file.reader(promer_out)

        for match in file_reader:
            if match.ref_name.endswith('__ends'):
                original_contig_name = match.ref_name[:-len('__ends')]
            else:
                original_contig_name = match.ref_name

            if match.hit_length_qry == match.qry_length and match.percent_identity >= min_percent_id:
                if original_contig_name not in circularized or circularized[original_contig_name].percent_identity < match.percent_identity:
                    circularized[original_contig_name] = copy.copy(match)

        if len(circularized):
            print('Using the following promer matches to circularize contigs:', file=log_fh)
            for contig_name in circularized:
                print(circularized[contig_name], file=log_fh)
        else:
            print('No suitable promer matches found', file=log_fh)

        return circularized


    @classmethod
    def _find_circular_using_prodigal(cls, outprefix, contigs_dict, circular_from_promer, log_fh):
        if len(contigs_dict) == len(circular_from_promer):
            return {}

        prodigal = circlator.external_progs.make_and_check_prog('prodigal')
        total_contig_length = sum([len(contigs_dict[x]) for x in contigs_dict if x not in circular_from_promer])
        prodigal_input_file = outprefix + '.for_prodigal.fa'
        prodigal_output_file = outprefix + '.prodigal.gff'
        with open(prodigal_input_file, 'w') as f:
            for seq_name in sorted(contigs_dict):
                if seq_name not in circular_from_promer:
                    print(contigs_dict[seq_name], file=f)

        if (total_contig_length < 20000):
            # prodigal needs -p meta option for sequences less than 20000
            # annoyingly newer version of prodigal has different -p option!
            if LooseVersion(prodigal.version()) >= LooseVersion('3.0'):
                p_option = "-p anon"
            else:
                p_option = "-p meta"

            cmd = ' '.join([
              'prodigal',
              '-i', prodigal_input_file,
              '-o', prodigal_output_file,
              '-f gff -c -q',
              p_option
            ])

        circlator.common.syscall(cmd)
        circularized = {}
        best_dist = {}

        with open(prodigal_output_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                contig, x, y, start, end, *z = line.rstrip().split('\t')
                new_inter = pyfastaq.intervals.Interval(int(start), int(end))
                new_dist = new_inter.distance_to_point(int(len(contigs_dict[contig])))

                if contig in circularized:
                    assert contig in best_dist
                    if new_dist < best_dist[contig]:
                        best_dist[contig] = new_dist
                        circularized[contig] = line.rstrip()
                else:
                    best_dist[contig] = new_dist
                    circularized[contig] = line.rstrip()

        if len(circularized):
            print('Using the following prodigal predictions to circularize contigs:', file=log_fh)
            for contig_name in circularized:
                print(circularized[contig_name], file=log_fh)
        else:
            print('No prodigal matches found', file=log_fh)

        return circularized
