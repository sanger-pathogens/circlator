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
        promer_mincluster=None,
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
        self.promer_mincluster = promer_mincluster
        self.outprefix = os.path.abspath(outprefix)
        self.verbose = verbose

        if ignore is None:
            self.ignore = set()
        else:
            with open(ignore) as f:
                self.ignore = {x.rstrip().split()[0] for x in f}


    @classmethod
    def _rename_contigs(cls, contigs_dict):
        rename_dict = {}
        new_contigs_dict = {}
        for old_name, contig in contigs_dict.items():
            new_name = old_name.split()[0]
            if new_name in rename_dict:
                raise Error('Non-unique name "' + new_name + '" in input contigs. Cannot continue')
            rename_dict[new_name] = old_name
            new_contigs_dict[new_name] = pyfastaq.sequences.Fasta(new_name, contig.seq)

        return new_contigs_dict, rename_dict


    @classmethod
    def _write_renamed_contigs(cls, contigs_dict, rename_dict, outfile):
        f = pyfastaq.utils.open_file_write(outfile)

        for old_name, contig in sorted(contigs_dict.items()):
            contig.id = rename_dict[old_name]
            print(contig, file=f)

        pyfastaq.utils.close(f)


    @classmethod
    def _max_length_from_fasta_file(cls, infile):
        freader = pyfastaq.sequences.file_reader(infile)
        return max([len(x) for x in freader])


    @classmethod
    def _write_fasta_plus_circularized_ends(cls, contigs, outfile, end_length, ignore=None):
        if ignore is None:
            ignore = set()

        if ignore.issuperset(set(contigs.keys())):
            return 0

        f = pyfastaq.utils.open_file_write(outfile)
        used_names = set(contigs.keys())
        seqs_written = 0

        for contig_name, contig in sorted(contigs.items()):
            if contig_name in ignore:
                continue

            print(contig, file=f)
            seqs_written += 1

            if len(contig) >= 2 * end_length:
                start_coord = end_length - 1
                end_coord = len(contig) - end_length
                new_name = contig.id + '__ends'
                assert new_name not in used_names
                used_names.add(new_name)
                new_contig = pyfastaq.sequences.Fasta(new_name, contig[end_coord:] + contig[:start_coord + 1])
                print(new_contig, file=f)
                seqs_written += 1

        pyfastaq.utils.close(f)
        return seqs_written


    @classmethod
    def _find_circular_using_promer(cls, outprefix, ref_genes_fa, contigs_dict, min_percent_id, end_length, log_fh, ignore=None, promer_mincluster=None):
        if ignore is None:
            ignore = set()
        promer_out = outprefix + '.promer'
        contigs_with_ends = outprefix + '.contigs_with_ends.fa'
        sequences_written = StartFixer._write_fasta_plus_circularized_ends(contigs_dict, contigs_with_ends, end_length, ignore=ignore)
        if sequences_written == 0:
            print('No sequences to check with promer. Skipping', file=log_fh)
            return {}

        prunner = pymummer.nucmer.Runner(
            contigs_with_ends,
            ref_genes_fa,
            promer_out,
            min_id=min_percent_id,
            promer=True,
            verbose=False,
            maxmatch=True,
            mincluster=promer_mincluster,
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
    def _find_circular_using_prodigal(cls, outprefix, contigs_dict, circular_from_promer, log_fh, ignore=None):
        if ignore is None:
            ignore = set()

        all_to_ignore = copy.copy(ignore)
        all_to_ignore.update(set(circular_from_promer.keys()))
        if all_to_ignore.issuperset(set(contigs_dict.keys())):
            print('No sequences left for which to look for genes using prodigal. Skipping', file=log_fh)
            return {}

        prodigal = circlator.external_progs.make_and_check_prog('prodigal')
        total_contig_length = 0
        prodigal_input_file = outprefix + '.for_prodigal.fa'
        prodigal_output_file = outprefix + '.prodigal.gff'
        with open(prodigal_input_file, 'w') as f:
            for seq_name in sorted(contigs_dict):
                if seq_name not in circular_from_promer and seq_name not in ignore:
                    total_contig_length += len(contigs_dict[seq_name])
                    print(contigs_dict[seq_name], file=f)


        if (total_contig_length < 20000):
            # prodigal needs -p meta option for sequences less than 20000
            # annoyingly newer version of prodigal has different -p option!
            if LooseVersion(prodigal.version) >= LooseVersion('3.0'):
                p_option = "-p anon"
            else:
                p_option = "-p meta"
        else:
            p_option = ''

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
                new_dist = new_inter.distance_to_point(int(len(contigs_dict[contig])/2))

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


    @classmethod
    def _rearrange_contigs(cls, contigs_dict, circular_from_promer, circular_from_prodigal, to_ignore, end_length, log_filename):
        log_fh = pyfastaq.utils.open_file_write(log_filename)
        log_prefix = '[fixstart]'
        print(log_prefix, 'id', 'break_point', 'gene_name', 'gene_reversed', 'new_name', 'skipped', sep='\t', file=log_fh)
        for contig_name, contig in sorted(contigs_dict.items()):
            original_length = len(contig)
            if contig.id in circular_from_promer:
                match = circular_from_promer[contig.id]
                if contig.id == match.ref_name:
                    if match.on_same_strand():
                        position = min(match.ref_start, match.ref_end)
                        contig.seq = contig.seq[position:] + contig.seq[:position]
                        print(log_prefix, contig.id, position + 1, match.qry_name, 'no', '-', '-', sep='\t', file=log_fh)
                    else:
                        position = max(match.ref_start, match.ref_end) + 1
                        contig.seq = contig.seq[position:] + contig.seq[:position]
                        contig.revcomp()
                        print(log_prefix, contig.id, position, match.qry_name, 'yes', '-', '-', sep='\t', file=log_fh)
                else:
                    assert contig.id + '__ends' == match.ref_name
                    if match.on_same_strand():
                        start_in_ends_contig = min(match.ref_start, match.ref_end)
                    else:
                        start_in_ends_contig = max(match.ref_start, match.ref_end) + 1

                    if start_in_ends_contig < end_length:
                        position = len(contig) - end_length + start_in_ends_contig
                    else:
                        position = start_in_ends_contig - end_length

                    if match.on_same_strand():
                        contig.seq = contig.seq[position:] + contig.seq[:position]
                        print(log_prefix, contig.id, position + 1, match.qry_name, 'no', '-', '-', sep='\t', file=log_fh)
                    else:
                        contig.seq = contig.seq[position:] + contig.seq[:position]
                        contig.revcomp()
                        print(log_prefix, contig.id, position, match.qry_name, 'yes', '-', '-', sep='\t', file=log_fh)
            elif contig.id in circular_from_prodigal:
                gff_contig_name, x, y, start, end, score, strand, *z = circular_from_prodigal[contig.id].rstrip().split('\t')
                assert contig.id == gff_contig_name

                if strand == '+':
                    position = int(start) - 1
                    contig.seq = contig.seq[position:] + contig.seq[:position]
                    print(log_prefix, contig.id, position + 1, 'prodigal', 'no', '-', '-', sep='\t', file=log_fh)
                else:
                    assert strand == '-'
                    position = int(end)
                    contig.seq = contig.seq[position:] + contig.seq[:position]
                    contig.revcomp()
                    print(log_prefix, contig.id, position + 1, 'prodigal', 'yes', '-', '-', sep='\t', file=log_fh)
            elif contig.id in to_ignore:
                print(log_prefix, contig.id, '-', '-', '-', '-', 'skipped', sep='\t', file=log_fh)
            else:
                print(log_prefix, contig.id, '-', '-', '-', '-', '-', sep='\t', file=log_fh)

            assert original_length == len(contig)

        pyfastaq.utils.close(log_fh)


    def run(self):
        renamed_contigs, contig_rename_dict = StartFixer._rename_contigs(self.input_assembly)
        if self.verbose:
            print('[fixstart] loaded input contigs', flush=True)

        end_extend = StartFixer._max_length_from_fasta_file(self.genes_fa)
        log_fh = pyfastaq.utils.open_file_write(self.outprefix + '.detailed.log')

        if self.verbose:
            print('[fixstart] Running promer to look for reference genes', flush=True)
        circ_with_promer = StartFixer._find_circular_using_promer(
            self.outprefix + '.promer',
            self.genes_fa,
            renamed_contigs,
            self.min_percent_identity,
            end_extend,
            log_fh,
            ignore=self.ignore,
            promer_mincluster=self.promer_mincluster,
        )

        if self.verbose:
            print('[fixstart] Running prodigal on sequences that had no promer match', flush=True)
        circ_with_prodigal = StartFixer._find_circular_using_prodigal(
            self.outprefix + '.prodigal',
            renamed_contigs,
            circ_with_promer,
            log_fh,
            ignore=self.ignore)

        if self.verbose:
            print('[fixstart] Fixing start positions of contigs', flush=True)
        StartFixer._rearrange_contigs(
            renamed_contigs,
            circ_with_promer,
            circ_with_prodigal,
            self.ignore,
            end_extend,
            self.outprefix + '.log'
        )

        pyfastaq.utils.close(log_fh)

        if self.verbose:
            print('[fixstart] Writing final FASTA file of contigs', flush=True)
        StartFixer._write_renamed_contigs(
            renamed_contigs,
            contig_rename_dict,
            self.outprefix + '.fasta'
        )
