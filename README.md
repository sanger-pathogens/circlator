circlator
=========

A tool to circularise bacterial genome assemblies.

For how to use Circlator, please see the [Circlator wiki page] [circlator wiki].


Installation
------------

circlator has the following dependencies, which need to be installed:
  * [BWA] [BWA] version >= 0.7.12
  * [prodigal] [prodigal] >= 2.6
  * [samtools] [samtools] (version 0.1.9 or 1.2)
  * [MUMmer] [mummer] version >= 3.23
  * [SPAdes] [spades] version >= 3.5.0


Once the dependencies are installed, install circlator using pip3:

    pip3 install circlator

Alternatively, you can download the latest release from this github repository,
or clone the repository. Then run the tests:

    python3 setup.py test

If the tests all pass, install:

    python3 setup.py install


Usage
-----

Please read the [Circlator wiki page] [circlator wiki] for usage instructions.


  [BWA]: http://bio-bwa.sourceforge.net/
  [circlator wiki]: https://github.com/sanger-pathogens/circlator/wiki
  [mummer]: http://mummer.sourceforge.net/
  [prodigal]: http://prodigal.ornl.gov/
  [samtools]: http://www.htslib.org/
  [spades]: http://bioinf.spbau.ru/spades

