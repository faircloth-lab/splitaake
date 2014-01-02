#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 02 January 2014 10:00 PST (-0800)
"""


from __future__ import absolute_import


def readfq(f):
    """From https://github.com/lh3/readfq
    This is a generator function
    """
    # this is a buffer keeping the last unprocessed line
    last = None
    while True:
        # the first record or a record following a fastq
        if not last:
            # search for the start of the next record
            for l in f:
                # fasta/q header line
                if l[0] in '>@':
                    # save this line
                    last = l[:-1]
                    break
        if not last:
            break
        #pdb.set_trace()
        name, seqs, last = last[1:], [], None
        # read the sequence
        for l in f:
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        # this is a fasta record
        if not last or last[0] != '+':
            # yield a fasta record
            yield name, ''.join(seqs), None
            if not last:
                break
        # this is a fastq record
        else:
            seq, leng, seqs = ''.join(seqs), 0, []
            # read the quality
            for l in f:
                seqs.append(l[:-1])
                leng += len(l) - 1
                # have read enough quality
                if leng >= len(seq):
                    last = None
                    # yield a fastq record
                    yield name, seq, ''.join(seqs)
                    break
            # reach EOF before reading enough quality
            if last:
                # yield a fasta record instead
                yield name, seq, None
                break


def writefq(f, seq):
    f.write("@{0}\n{1}\n+\n{2}\n".format(*seq))
