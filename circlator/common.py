import sys
import os
import subprocess

class Error (Exception): pass

def syscall(cmd, allow_fail=False, verbose=False):
    if verbose:
        print('syscall:', cmd, flush=True)
    try:
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:
        errors = error.output.decode()
        if allow_fail:
            return False, errors
        else:
            print('The following command failed with exit code', error.returncode, file=sys.stderr)
            print(cmd, file=sys.stderr)
            print('\nThe output was:\n', file=sys.stderr)
            print(errors, file=sys.stderr, flush=True)
            sys.exit(1)

    return True, None


def syscall_get_stdout(cmd):
    try:
        out = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE).communicate()[0].decode('utf-8').rstrip()
        return out.split('\n')
    except:
        raise Error('Error in system call. I tried to run:\n' + str(cmd))


def decode(x):
    try:
        s = x.decode()
    except:
        return x
    return s


def check_files_exist(filenames):
    '''Dies if any files in the list of filenames does not exist'''
    files_not_found = [x for x in filenames if not os.path.exists(x)]
    if len(files_not_found):
        for filename in files_not_found:
            print('File not found: "', filename, '"', sep='', file=sys.stderr)
        raise Error('File(s) not found. Cannot continue')
