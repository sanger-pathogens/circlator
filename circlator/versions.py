import sys
import openpyxl
import pyfastaq
import pymummer
import pysam
from circlator import external_progs, __version__
from circlator import __version__ as circlator_version


def get_all_versions(filehandle, raise_error=True, debug=False):
    if filehandle is not None:
        print('Circlator version:', circlator_version, file=filehandle)
        print('\nExternal dependencies:', file=filehandle)

    external_progs.check_all_progs(verbose=False, raise_error=raise_error, filehandle=filehandle, debug=debug)

    if filehandle is not None:
        print('\nPython version:', file=filehandle)
        print(sys.version, file=filehandle)
        print('\nPython dependencies:', file=filehandle)

    found_bad_module = False

    for module in ['openpyxl', 'pyfastaq', 'pymummer', 'pysam']:
        try:
            version = eval(module + '.__version__')
            path = eval(module + '.__file__')
        except:
            version = 'NOT_FOUND'
            path = 'NOT_FOUND'
            found_bad_module = True

        if filehandle is not None:
            print(module + '\t' + version + '\t' + path, file=filehandle)

    if raise_error and found_bad_module:
        print('Some dependencies not satisfied. Cannot continue. Try running: circlator progcheck', file=sys.stderr)
        sys.exit(1)

