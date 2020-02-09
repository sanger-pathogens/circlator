import re
import sys
from circlator import program, common
from circlator import __version__ as circlator_version
import shutil
import pyfastaq

class Error (Exception): pass

prog_to_env_var = {
    'samtools': 'CIRCLATOR_SAMTOOLS',
    'spades': 'CIRCLATOR_SPADES',
    'canu': 'CIRCLATOR_CANU',
    'racon': 'CIRCLATOR_RACON',
}

prog_to_version_cmd = {
    'bwa': ('', re.compile(r'^Version: ([0-9\.]+)')),
    'nucmer': ('--version', re.compile(r'([0-9\.]+)')),
    'prodigal': ('-v', re.compile(r'^Prodigal V([0-9\.]+):')),
    'samtools': ('', re.compile(r'Version: (\d+\.\d+[\.\d]*)')),
    'spades': ('-v', re.compile(r'v.?([0-9][0-9\.]+)')),
    'canu': ('-version', re.compile(r'^Canu \D*([\d][\d\.]+)')),
    'minimap2': ('-V', re.compile(r'([0-9\.]+)')),
    'miniasm': ('-V', re.compile(r'([0-9\.]+)')),
    'racon': ('--version', re.compile(r'v.?([0-9][0-9\.]+)')),
}

min_versions = {
    'bwa': '0.7.12',
    'nucmer': '3.1',
    'prodigal': '2.6',
    'samtools': '0.1.19',
    'spades': '3.11.1', # this is the first version to support python3
    'canu': '0.0',
    'minimap2': '2.10',
    'miniasm': '0.2',
    'racon': '1.2.1',
}


bad_versions = {
    'spades': '3.6.1'
}


prog_name_to_default = {
    'bwa': 'bwa',
    'nucmer': 'nucmer',
    'prodigal': 'prodigal',
    'spades': 'spades.py',
    'samtools': 'samtools',
    'canu': 'canu',
    'minimap2':'minimap2',
    'miniasm': 'miniasm',
    'racon': 'racon'
}

not_required = {'spades', 'canu','racon'}

def handle_error(message, raise_error=True):
    if raise_error:
        raise Error(message + '\nCannot continue')
    else:
        print(message)


def make_and_check_prog(name, verbose=False, raise_error=True, filehandle=None, debug=False, required=False):
    p = program.Program(
        prog_name_to_default[name],
        prog_to_version_cmd[name][0],
        prog_to_version_cmd[name][1],
        environment_var=prog_to_env_var.get(name, None),
        debug=debug
    )

    if not p.in_path():
        if required:
            die = True
        else:
            die = raise_error and (name not in not_required)
        handle_error("WARNING: Didn't find " + name + " in path. Looked for:" + p.path, raise_error=die)
        return p

    version = p.version

    if version is None:
        handle_error('Found ' + name + " but couldn't get version.", raise_error=raise_error)
        return p

    if not p.version_at_least(min_versions[name]):
        handle_error('Version of ' + name + ' too low. I found ' + p.version + ', but must be at least ' + min_versions[name] + '. Found here:\n' + p.from_which, raise_error=raise_error)
        return p

    if name == 'spades' and p.version == bad_versions['spades']:
        handle_error('ERROR! SPAdes version ' + bad_versions['spades'] + ' is incompatible with Circlator. Please use SPAdes 3.7.1', raise_error=raise_error)
        return p

    if name == 'spades' and not p.version.startswith('3.7.'):
        print('WARNING: SPAdes version', p.version, 'is being used. It will work, but better results are usually obtained from Circlator using SPAdes version 3.7.1. Although 3.7.1 is not the latest version, we recommend it for Circlator.', file=sys.stderr)

    if verbose:
        print(name, p.version, p.from_which, sep='\t')

    if filehandle:
        print(name, p.version, p.from_which, sep='\t', file=filehandle)

    return p


def check_all_progs(verbose=False, raise_error=False, filehandle=None, debug=False, assembler=None):
    for prog in sorted(prog_name_to_default):
        if debug:
            print('__________ checking', prog, '____________', flush=True)
        make_and_check_prog(prog, verbose=verbose, raise_error=raise_error, filehandle=filehandle, debug=debug, required=prog==assembler)
