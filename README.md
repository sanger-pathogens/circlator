circlator
=========

A tool to circularise bacterial genome assemblies.

For how to use Circlator, please see the [Circlator wiki page] [circlator wiki].


Installation
------------

circlator has the following dependencies, which need to be installed:
  * [BWA] [BWA] version >= 0.7.12
  * [prodigal] [prodigal] >= 2.6
  * [SAMtools] [samtools] (version 0.1.9 or 1.2)
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


References
----------

[BWA] [BWA paper]: Li, H et al. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997.

[MUMmer] [mummer paper]: Kurtz, S. et al. Versatile and open software for comparing large genomes. Genome Biol. 5, R12 (2004).

[Prodigal] [prodigal paper]: Hyatt, D. et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11, 119 (2010).

[SAMtools] [samtools paper]: Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078–9 (2009).

[SPAdes] [spades paper]: Bankevich, A. et al. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. J. Comput. Biol. 19, 455–77 (2012).


  [BWA]: http://bio-bwa.sourceforge.net/
  [BWA paper]: http://arxiv.org/abs/1303.3997
  [circlator wiki]: https://github.com/sanger-pathogens/circlator/wiki
  [mummer]: http://mummer.sourceforge.net/
  [mummer paper]: http://genomebiology.com/content/5/2/R12
  [prodigal]: http://prodigal.ornl.gov/
  [prodigal paper]: http://www.biomedcentral.com/1471-2105/11/119
  [samtools]: http://www.htslib.org/
  [samtools paper]: http://bioinformatics.oxfordjournals.org/content/25/16/2078.abstract
  [spades]: http://bioinf.spbau.ru/spades
  [spades paper]: http://online.liebertpub.com/doi/abs/10.1089/cmb.2012.0021

