import re
import pyfastaq


class Error (Exception): pass


genetic_code = 11
aa_to_dna = {}
for nucleotide, aa in pyfastaq.genetic_codes.codes[genetic_code].items():
    aa_to_dna[aa] = nucleotide


class UniprotDownloader:
    def __init__(self, min_gene_length=333, max_gene_length=500, uniprot_search='dnaa', header_regex='dnaa', header_regex_ignorecase=True):
        self.min_gene_length = min_gene_length
        self.max_gene_length = max_gene_length
        self.uniprot_search = re.sub('\s+', '+', uniprot_search.strip())
        if header_regex_ignorecase:
            self.header_regex = re.compile(header_regex, re.IGNORECASE)
        else:
            self.header_regex = re.compile(header_regex)


    def _get_uniprot_url(self):
        return 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=no&query=' + self.uniprot_search + '&force=no&format=fasta'


    def _download_from_uniprot(self, outfile):
        pyfastaq.utils.syscall('wget -O ' + outfile + ' "' + self._get_uniprot_url() + '"')


    def _header_to_genus_species(self, header):
        if 'OS=' not in header:
            return None

        pos = header.find('OS=')
        assert pos != -1
        genus_species = header[pos+3:].split(' ')[:2]
        return tuple(genus_species)


    def _header_matches_regex(self, sequence):
        return self.header_regex.search(sequence.id) is not None


    def _check_sequence(self, sequence, seen_genus_species):
        if not(self.min_gene_length <= len(sequence) <= self.max_gene_length):
            return False, 'Too long or short'
        elif not sequence[0] == 'M':
            return False, 'Does not start with M'
        elif not self._header_matches_regex(sequence):
            return False, 'No match to regex in name'
        else:
            genus_species = self._header_to_genus_species(sequence.id)
            if genus_species is None:
                return False, 'Error getting genus species'
            elif genus_species in seen_genus_species:
                return False, 'Duplicate genus species ' + str(genus_species)
            else:
                seen_genus_species.add(genus_species)
                return True, None


    def _append_stop_to_seq(self, sequence):
        if not sequence[-1] == '*':
            sequence.seq += '*'


    def _reverse_translate(self, sequence):
        amino_acids = list(sequence.seq)
        try:
            nucleotides = [aa_to_dna[aa] for aa in amino_acids]
        except:
            raise Error('Error reverse translate')

        sequence.seq = ''.join(nucleotides)


    def run(self, outprefix):
        aa_file = outprefix + '.aa.fa'
        self._download_from_uniprot(aa_file)
        log_fh = pyfastaq.utils.open_file_write(outprefix + '.log')
        fasta_fh = pyfastaq.utils.open_file_write(outprefix + '.nucleotides.fa')
        seen_genus_species = set()

        for seq in pyfastaq.sequences.file_reader(aa_file):
            ok, error = self._check_sequence(seq, seen_genus_species)

            if ok:
                self._append_stop_to_seq(seq)
                try:
                    self._reverse_translate(seq)
                except:
                    print('Failed reverse translation', seq.id, sep='\t', file=log_fh)
                    continue

                print(seq, file=fasta_fh)
            else:
                print(error, seq.id, sep='\t', file=log_fh)

        pyfastaq.utils.close(log_fh)
        pyfastaq.utils.close(fasta_fh)
