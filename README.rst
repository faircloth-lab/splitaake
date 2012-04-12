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

