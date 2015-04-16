import os
import pysam
import pyfastaq
from circlator import common, external_progs

class Error (Exception): pass


index_extensions = [
        'amb',
        'ann',
        'bwt',
        'pac',
        'sa'
]


def bwa_index(infile, outprefix=None, bwa='bwa', verbose=False):
    if outprefix is None: 
        outprefix = infile

    missing = [not os.path.exists(outprefix + '.' + x) for x in index_extensions]
    if True not in missing:
        return

    cmd = ' '.join([
        bwa,  'index',
        '-p', outprefix,
        infile
    ])
    common.syscall(cmd, verbose=verbose)


def bwa_index_clean(prefix):
    for e in index_extensions:
        try:
            os.unlink(prefix + '.' + e)
        except:
            pass


def bwa_mem(
      ref,
      reads,
      outfile,
      threads=1,
      bwa_options = '-x pacbio',
      verbose=False,
      index=None
    ):

    samtools = external_progs.make_and_check_prog('samtools', verbose=verbose)
    bwa = external_progs.make_and_check_prog('bwa', verbose=verbose)
    unsorted_bam = outfile + '.tmp.unsorted.bam'
    tmp_index = outfile + '.tmp.bwa_index'
    bwa_index(ref, outprefix=tmp_index, verbose=verbose, bwa=bwa.exe())

    cmd = ' '.join([
        bwa.exe(), 'mem',
        bwa_options,
        '-t', str(threads),
        tmp_index,
        reads, 
        '|',
        samtools.exe(), 'view',
        '-F 0x0800',
        '-T', ref,
        '-b',
        '-o', unsorted_bam,
        '-',
    ])

    common.syscall(cmd, verbose=verbose)
    bwa_index_clean(tmp_index)
    threads = min(4, threads)
    thread_mem = int(500 / threads)
    
    cmd = ' '.join([
        samtools.exe(), 'sort',
        '-@', str(threads),
        '-m', str(thread_mem) + 'M',
        unsorted_bam,
        outfile[:-4]
    ])
 
    common.syscall(cmd, verbose=verbose)
    os.unlink(unsorted_bam)

    cmd = samtools.exe() + ' index ' + outfile
    common.syscall(cmd, verbose=verbose)


def aligned_read_to_read(read, revcomp=True, qual=None, ignore_quality=False):
    '''Returns Fasta or Fastq sequence from pysam aligned read'''
    if read.qual is None or ignore_quality:
        if qual is None or ignore_quality:
            seq = pyfastaq.sequences.Fasta(read.qname, common.decode(read.seq))
        else:
            seq = pyfastaq.sequences.Fastq(read.qname, common.decode(read.seq), qual * read.query_length)
    else:
        if qual is None:
            seq = pyfastaq.sequences.Fastq(read.qname, common.decode(read.seq), common.decode(read.qual))
        else:
            seq = pyfastaq.sequences.Fastq(read.qname, common.decode(read.seq), qual * read.query_length)

    if read.is_reverse and revcomp:
        seq.revcomp()

    return seq
