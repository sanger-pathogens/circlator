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
    'fixstart',
    'mapping',
    'merge',
    'minimus2',
    'program',
    'tasks',
    'versions',
]

from circlator import *

