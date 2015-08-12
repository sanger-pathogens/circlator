import argparse
import sys
from circlator import external_progs

def run():
    parser = argparse.ArgumentParser(
        description = 'Checks all dependencies are found and are correct versions',
        usage = 'circlator progcheck')
    options = parser.parse_args()

    external_progs.check_all_progs(verbose=True, raise_error=False)
