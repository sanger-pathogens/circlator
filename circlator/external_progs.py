import re
import sys
from circlator import program, common
import shutil
import pyfastaq

class Error (Exception): pass


# Travis is using python3.4, and actually "python" in travis means
# python3.4, not python2. SPAdes throws an error about not being
# compatible with python3.4.
# This means we need to explicitly run SPAdes with python2.
class Spades(program.Program):
    def exe(self):
        return 'python2 ' + shutil.which(self.path)


prog_to_env_var = {
    'samtools': 'CIRCLATOR_SAMTOOLS',
    'spades': 'CIRCLATOR_SPADES',
}


prog_to_version_cmd = {
    'bwa': ('', re.compile('^Version: ([0-9\.]+)')),
    'nucmer': ('--version', re.compile('^NUCmer \(NUCleotide MUMmer\) version ([0-9\.]+)')),
    'prodigal': ('-v', re.compile('^Prodigal V([0-9\.]+):')),
    'samtools': ('', re.compile('^Version: ([0-9\.]+)')),
    'spades': ('', re.compile('^SPAdes genome assembler v.([0-9\.]+)')),
}


min_versions = {
    'bwa': '0.7.12',
    'nucmer': '3.1',
    'prodigal': '2.6',
    'samtools': '0.1.19',
    'spades': '3.5.0',
}


max_versions = {
    'spades': '3.6.0'
}

prog_name_to_default = {
    'bwa': 'bwa',
    'nucmer': 'nucmer',
    'prodigal': 'prodigal',
    'spades': 'spades.py',
    'samtools': 'samtools',
}


prog_builders = {
    'spades': Spades
}


def handle_error(message, raise_error=True):
    if raise_error:
        raise Error(message + ' Cannot continue')
    else:
        print(message)


def make_and_check_prog(name, verbose=False, raise_error=True, filehandle=None):
    builder = prog_builders.get(name, program.Program)
    p = builder(
        prog_name_to_default[name],
        prog_to_version_cmd[name][0],
        prog_to_version_cmd[name][1],
        environment_var=prog_to_env_var.get(name, None)
    )

    if not p.in_path():
        handle_error("Didn't find " + name + " in path. Looked for:" + p.path, raise_error=raise_error)
        return p

    version = p.version()

    if version is None:
        handle_error('Found ' + name + " but couldn't get version.", raise_error=raise_error)
        return p

    if not p.version_at_least(min_versions[name]):
        handle_error('Version of ' + name + ' too low. I found ' + p.version() + ', but must be at least ' + min_versions[name] + '. Found here:\n' + p.full_path, raise_error=raise_error)
        return p

    if name in max_versions and not p.version_at_most(max_versions[name]):
        handle_error('Version of ' + name + ' too high. I found ' + p.version() + ', but must be at most ' + max_versions[name] + '. Found here:\n' + p.full_path, raise_error=raise_error)
        return p

    if verbose:
        print('Using', name, 'version', p.version(), 'as', p.full_path)

    if filehandle:
        print('Using', name, 'version', p.version(), 'as', p.full_path, file=filehandle)

    return p


def check_all_progs(verbose=False, raise_error=False, filehandle=None):
    if filehandle is not None:
        print(sys.argv[0], 'version', common.version, file=filehandle)

    for prog in sorted(prog_name_to_default):
        make_and_check_prog(prog, verbose=verbose, raise_error=raise_error, filehandle=filehandle)

