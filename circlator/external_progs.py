import re
from circlator import program

class Error (Exception): pass

prog_to_env_var = {
    'spades': 'CIRCLATOR_SPADES',
}


prog_to_version_cmd = {
    'bwa': ('', re.compile('^Version: ([0-9\.]+)')),
    'nucmer': ('--version', re.compile('^NUCmer \(NUCleotide MUMmer\) version ([0-9\.]+)')),
    'samtools': ('', re.compile('^Version: ([0-9\.]+)')),
    'spades': ('', re.compile('^SPAdes genome assembler v.([0-9\.]+)')),
}


min_versions = {
    'bwa': '0.7.12',
    'nucmer': '3.1',
    'samtools': '0.1.19',
    'spades': '3.5.0',
}


prog_name_to_default = {
    'bwa': 'bwa',
    'nucmer': 'nucmer',
    'spades': 'spades.py',
    'samtools': 'samtools',
}

def make_and_check_prog(name, verbose=False):
    p = program.Program(
        prog_name_to_default[name], 
        prog_to_version_cmd[name][0],
        prog_to_version_cmd[name][1],
        environment_var=prog_to_env_var.get(name, None)
    )

    if not p.version_at_least(min_versions[name]):
        raise Error('Version of ' + name + ' too low. I found ' + p.version() + ', but must be at least ' + min_versions[name])

    if verbose:
        print('Using', name, 'version', p.version())

    return p


def check_all_progs(verbose=False):
    for prog in prog_name_to_default:
        make_and_check_prog(prog, verbose=verbose)
