# Circlator
A tool to circularize genome assemblies.

PLEASE NOTE: we currently do not have the resources to provide support for Circlator, so please do not expect a reply if you flag any issue. Apologies for the inconvenience.

[![Unmaintained](http://unmaintained.tech/badge.svg)](http://unmaintained.tech/)  
[![Build Status](https://travis-ci.org/sanger-pathogens/circlator.svg?branch=master)](https://travis-ci.org/sanger-pathogens/circlator)   
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/circlator/blob/master/LICENSE)   
[![status](https://img.shields.io/badge/GenomeBiology-10.1186-brightgreen.svg)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0849-0)  
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/circlator/README.html)  
[![Container ready](https://img.shields.io/badge/container-ready-brightgreen.svg)](https://quay.io/repository/biocontainers/circlator)  
[![Docker Build Status](https://img.shields.io/docker/build/sangerpathogens/circlator.svg)](https://hub.docker.com/r/sangerpathogens/circlator)  
[![Docker Pulls](https://img.shields.io/docker/pulls/sangerpathogens/circlator.svg)](https://hub.docker.com/r/sangerpathogens/circlator)  
[![codecov](https://codecov.io/gh/sanger-pathogens/circlator/branch/master/graph/badge.svg)](https://codecov.io/gh/sanger-pathogens/circlator)

## Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
  * [Usage](#usage)
  * [License](#license)
  * [Feedback/Issues](#feedbackissues)
  * [Citation](#citation)
  * [Further Information](#further-information)

## Introduction
Circlator will attempt to identify each circular sequence and output a linearised version of it. It does this by assembling all reads that map to contig ends and comparing the resulting contigs with the input assembly.

## Installation
For how to install Circlator, please see the [Circlator website](http://sanger-pathogens.github.io/circlator/).

## Usage
```
Usage: circlator <command> [options] <required arguments>

To get minimal usage for a command use:
circlator command

To get full help for a command use one of:
circlator command -h
circlator command --help


Available commands:

all        Run mapreads, bam2reads, assemble, merge, clean, fixstart
mapreads   Map reads to assembly
bam2reads  Make reads from mapping to be reassembled
assemble   Run assembly using reads from bam2reads
merge      Merge original assembly and new assembly made by assemble
clean      Remove small and completely contained contigs from assembly
fixstart   Change start position of circular sequences
minimus2   Run the minimus2 based circularisation pipeline
get_dnaa   Download file of dnaA (or other of user's choice) genes
progcheck  Checks dependencies are installed
test       Run Circlator on a small test set
version    Print version and exit
```
For more information on how to use Circlator, please see the [Circlator wiki page](https://github.com/sanger-pathogens/circlator/wiki).

## License
Circlator is free software, licensed under [GPLv3](https://github.com/sanger-pathogens/circlator/blob/master/LICENSE).

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/sanger-pathogens/circlator/issues) or email circlator-help@sanger.ac.uk.

PLEASE NOTE: we currently do not have the resources to provide support for Circlator, so please do not expect a reply if you flag any issue. Apologies for the inconvenience.

## Citation
If you use this software please cite:

__Circlator: automated circularization of genome assemblies using long sequencing reads__   
Hunt et al, Genome Biology 2015 Dec 29;16(1):294. doi: [10.1186/s13059-015-0849-0](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0849-0)

__[BWA](http://arxiv.org/abs/1303.3997):__ Li, H et al. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997.  

__[MUMmer](http://genomebiology.com/content/5/2/R12):__ Kurtz, S. et al. Versatile and open software for comparing large genomes. Genome Biol. 5, R12 (2004).  

__[Prodigal](http://www.biomedcentral.com/1471-2105/11/119):__ Hyatt, D. et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11, 119 (2010).  

__[SAMtools](http://bioinformatics.oxfordjournals.org/content/25/16/2078.abstract):__ Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078–9 (2009).  

__[SPAdes](http://online.liebertpub.com/doi/abs/10.1089/cmb.2012.0021)__: Bankevich, A. et al. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. J. Comput. Biol. 19, 455–77 (2012)  

## Further Information
For more information on this software see:
* the [Circlator website](http://sanger-pathogens.github.io/circlator/)
* the [Circlator wiki page](https://github.com/sanger-pathogens/circlator/wiki)