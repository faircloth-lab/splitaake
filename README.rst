Purpose
========

`splitaake` is a reasonably easy method for demultiplexing Illumina reads using either Hamming distance or Levenshtein (edit) distance sequence tags.  A the moment, Hamming distance is the default, which is similar to the approach Illumina uses in their own software.  `splitaake` differs from the Illumina software in that it can demultiplex sequence tags of many lengths, rather than just the "standard" TruSeq index length of 6 nucleotides.

Notes
======

`splitaake` is under development.  This means that it may break, be a pain to get running, etc.  I'll improve as I have time.  Please feel free to suggest contributions/changes/additions.

Design
=======

`splitaake` is currently designed as a single-core application meaning that it does not parallelize the process of demultiplexing.  After a number of tests, I've found that for most files, a single core approach is reasonably fast (and sometimes faster) than multi-core options, particularly when you wish to work with `gzipped` fastq files.

I'm still testing additional ways of demultiplexing data in parallel.  Hopefully more on this front soon...

Dependencies
=============

* seqtools_ ("working" branch)
* jellyfish_ (at the moment - quite fast Hamming implementation)

Running
========

- generate a config file mapping indexes to filenames.  This file is named
  `map.conf`, as used below::

    TruSeq1:ATCACGATCT
    TruSeq2:CGATGTATCT
    TruSeq3:TTAGGCATCT
    TruSeq4:TGACCAATCT
    TruSeq5:ACAGTGATCT
    TruSeq6:GCCAATATCT
    TruSeq7:CAGATCATCT
    TruSeq8:ACTTGAATCT
    TruSeq9:GATCAGATCT
    TruSeq10:TAGCTTATCT
    TruSeq11:GGCTACATCT
    TruSeq12:CTTGTAATCT

- run splitaake::

    python splitaake.py L007_R1.fastq.gz L007_R2.fastq.gz L007_R3.fastq.gz map.conf --section taxa

- this will identify your reads and create a directory `dmux` containing your
  reads in interleaved, fastq, gzip files, like so::

    dmux/
        TruSeq1.fastq.gz
	TruSeq2.fastq.gz
	TruSeq3.fastq.gz
	...

.. _jellyfish: https://github.com/sunlightlabs/jellyfish
.. _seqtools: https://github.com/faircloth-lab/seqtools
