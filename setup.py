import os
import shutil
import sys
import glob
from setuptools import setup, find_packages


setup(
    name='circlator',
    version='1.4.0',
    description='circlator: a tool to circularise genome assemblies',
    packages = find_packages(),
    package_data={'circlator': ['data/*']},
    author='Martin Hunt, Nishadi De Silva',
    author_email='path-help@sanger.ac.uk',
    url='https://github.com/sanger-pathogens/circlator',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    tests_require=['nose >= 1.3'],
    install_requires=[
        'openpyxl',
        'pyfastaq >= 3.12.1',
        'pysam >= 0.8.1',
        'pymummer>=0.9.0',
    ],
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)
