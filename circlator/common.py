import sys
import subprocess

version = '0.2.0'

def syscall(cmd, allow_fail=False, verbose=False):
    if verbose:
        print('syscall:', cmd, flush=True)
    try:
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:
        errors = error.output.decode()
        print('The following command failed with exit code', error.returncode, file=sys.stderr)
        print(cmd, file=sys.stderr)
        print('\nThe output was:\n', file=sys.stderr)
        print(errors, file=sys.stderr, flush=True)

        if allow_fail:
            return False, errors
        else:
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
