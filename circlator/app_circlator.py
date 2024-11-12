#!/usr/bin/env python3

import argparse
import sys

tasks = {
    'all': 'Run mapreads, bam2reads, assemble, merge, clean, fixstart',
    'progcheck': 'Checks dependencies are installed',
    'mapreads': 'Map reads to assembly',
    'bam2reads': 'Make reads from mapping to be reassembled',
    'assemble': 'Run assembly using reads from bam2reads',
    'merge': 'Merge original assembly and new assembly made by assemble',
    'clean': 'Remove small and completely contained contigs from assembly',
    'fixstart': 'Change start position of circular sequences',
    'minimus2': 'Run the minimus2 based circularisation pipeline',
    'get_dnaa': 'Download file of dnaA (or other of user\'s choice) genes',
    'test': 'Run Circlator on a small test set',
    'version': 'Print version and exit',
}


ordered_tasks = [
    'all',
    'mapreads',
    'bam2reads',
    'assemble',
    'merge',
    'clean',
    'fixstart',
    'minimus2',
    'get_dnaa',
    'progcheck',
    'test',
    'version'
]


def print_usage_and_exit():
    print('Usage: circlator <command> [options] <required arguments>', file=sys.stderr)
    print('\nTo get minimal usage for a command use:\ncirclator command', file=sys.stderr)
    print('\nTo get full help for a command use one of:\ncirclator command -h\ncirclator command --help\n', file=sys.stderr)
    print('\nAvailable commands:\n', file=sys.stderr)
    max_task_length = max([len(x) for x in list(tasks.keys())]) + 1
    for task in ordered_tasks:
        print('{{0: <{}}}'.format(max_task_length).format(task), tasks[task], sep=' ', file=sys.stderr)
    sys.exit(0)

if len(sys.argv) == 1 or sys.argv[1] in ['-h', '-help', '--help']:
    print_usage_and_exit()

task = sys.argv.pop(1)

if task not in tasks:
    print('Task "' + task + '" not recognised. Cannot continue.\n', file=sys.stderr)
    print_usage_and_exit()

exec('import circlator.tasks.' + task)
exec('circlator.tasks.' + task + '.run()')
