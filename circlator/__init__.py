from pkg_resources import get_distribution

try:
    __version__ = get_distribution('circlator').version
except:
    __version__ = 'local'


__all__ = [
    'assemble',
    'assembly',
    'bamfilter',
    'clean',
    'common',
    'dnaa',
    'external_progs',
    'mapping',
    'merge',
    'minimus2',
    'program',
    'start_fixer',
    'tasks',
    'versions',
]

from circlator import *

