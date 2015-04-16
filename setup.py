import os
import shutil
import sys
import glob
from setuptools import setup, find_packages


setup(
    name='circlator',
    version='0.2.0',
    description='circlator: a tool to circularise bacterial genome assemblies',
    packages = find_packages(),
    package_data={'circlator': ['data/*']},
    author='Martin Hunt, Nishadi De Silva',
    author_email='path-help@sanger.ac.uk',
    url='https://github.com/sanger-pathogens/circlator',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    install_requires=[
        'nose >= 1.3',
        'openpyxl',
        'pyfastaq >= 3.3.0',
        'pysam >= 0.8.1',
        'pymummer>=0.2.0',
        'bio_assembly_refinement>=0.1.0',
    ],
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)
