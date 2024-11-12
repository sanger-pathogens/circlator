import argparse
import sys
from circlator import versions

def run():
    parser = argparse.ArgumentParser(
        description = 'Checks all dependencies are found and are correct versions',
        usage = 'circlator progcheck'
    )
    parser.add_argument('--debug', action='store_true', help='Debug mode with very verbose output')
    options = parser.parse_args()
    versions.get_all_versions(sys.stdout, raise_error=False, debug=options.debug)
