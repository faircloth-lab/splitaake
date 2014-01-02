#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 02 January 2014 09:38 PST (-0800)

Description: A quality-aware fastq file demultiplexer for Illumina (Casava 1.8+)
data that demultiplexes Levenshtein/edit- distance sequence tags rather than
using  Hamming distance (but can use Hamming distance, too).  Works with gzip
and uncompressed outputs from Casava in single or multiple files.

Requires separate R1 (Read 1), R2 (Index read), R3 (Read 2) data as output from
Casava when choosing not to demultiplex, e.g.:

configureBclToFastq.pl \
    --input-dir /your/path/to/run_folder/Data/Intensities/BaseCalls \
    --output-dir /where/you/want/the/output \
    --sample-sheet /my/path/to/SampleSheet.csv \
    --use-bases-mask Y*,Y7,Y*

When sample sheet looks like:

FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject
D109LACXX,8,nodmux,,,No demultiplexing,N,D109LACXX,BCF,nodmux
"""


import os
import gzip
import glob
import numpy
import ConfigParser
from itertools import izip
from splitaake.fastq import readfq, writefq
from collections import defaultdict
from splitaake.core import motd

from jellyfish import levenshtein_distance as levenshtein
from jellyfish import hamming_distance as hamming

import pdb


class Tags:
    """Tagging object to hold tag and tag mapping information"""
    def __init__(self, conf, section, no_correct=False):
        config = ConfigParser.ConfigParser()
        config.read(conf)
        tags = config.items(section)
        if no_correct:
            nc = config.items('no correct')
            self.__get_no_correct_set(nc)
        else:
            self.no_correct = []
        self.__name_to_tag(tags)
        self.__tag_to_name()
        self.__get_min_length()
        self.seqs = self.seq_dict.keys()
        self.count = len(self.name_dict)

    def __name_to_tag(self, tags):
        self.name_dict = dict(tags)

    def __tag_to_name(self):
        self.seq_dict = {k: v for v, k in self.name_dict.iteritems()}

    def __get_min_length(self):
        l = set([len(s) for s in self.seq_dict.keys()])
        assert len(l) == 1
        self.length = l.pop()

    def __get_no_correct_set(self, nc):
        if nc:
            self.no_correct = set([v for k, v in nc])
        else:
            self.no_correct = []

    def __open_zip_file(self, pth):
        return gzip.open(pth, 'wb')

    def create_zip_files(self, pth):
        if not os.path.exists(pth):
            os.makedirs(pth)
        self.files = defaultdict(dict)
        for name, seq in self.name_dict.iteritems():
            for i in [1, 2]:
                r = "{0}_R{1}.fastq.gz".format(name, i)
                self.files[seq][i] = self.__open_zip_file(
                    os.path.join(pth, r)
                )
        for name in ['unknown', 'lowqual']:
            for i in [1, 2, 3]:
                r = "{0}-R{1}.fastq.gz".format(name, i)
                self.files[name][i] = self.__open_zip_file(
                    os.path.join(pth, r)
                )

    def close_zip_files(self):
        for f in self.files.values():
            for i in [1, 2]:
                f[i].close()


def change_read_num(sequence, read, index):
    """update the reads with preferred header info giving index"""
    si = sequence[0].split(' ')
    sil = list(si[1])
    if read == 2:
        sil[0] = '2'
    new_si = "{0}{1}".format(''.join(sil), index)
    new_head = ' '.join([si[0], new_si])
    sl = list(sequence)
    sl[0] = new_head
    #pdb.set_trace()
    return tuple(sl)


def get_quality(string, min_accept, min_mean_accept, ascii_min=33,
                quality_min=0):
    """return True if the min quality >= 10 and mean quality >= 20"""
    qa = numpy.fromstring(string, dtype='uint8')
    sanger = qa - ascii_min + quality_min
    if min(sanger) >= min_accept and numpy.mean(sanger) >= min_mean_accept:
        return True
    else:
        return False


def filtered(read):
    """Ensure we're skipping any reads that are filtered"""
    try:
        if read[0].split(' ')[1].split(':')[1] == 'N':
            return False
        else:
            return True
    except IndexError as e:
        if e.message == 'list index out of range':
            # raise warning in log
            return False


def get_index(index, length):
    """trim an index sequence down if the sequence used
    is shorter than the sequence converted from BCL"""
    if length is not None:
        idx = index[1][:length]
        idx_qual = index[2][:length]
    else:
        idx = index[1]
        idx_qual = index[2]
    return idx, idx_qual


def read_fastq(file):
    if os.path.splitext(file)[1] == ".gz":
        infile = gzip.open(file, "rU")
    elif os.path.splitext(file)[1] == ".fastq":
        infile = open(file, "rU")
    else:
        raise IOError("The files don't appear to be *.fastq or *.fastq.gz")
    return readfq(infile)


def write_to_r1_and_r2(tags, match, r1, r2):
    writefq(tags.files[match][1], r1)
    writefq(tags.files[match][2], r2)


def write_to_unk(tags, fset, r1, i, r2):
    writefq(tags.files[fset][1], r1)
    writefq(tags.files[fset][2], i)
    writefq(tags.files[fset][3], r2)


def main(args, parser):
    """main loop"""
    if not args.quiet:
        motd()
    # setup tags object
    tags = Tags(args.tagmap, args.section, args.no_correct)
    # create output files
    tags.create_zip_files(args.output)
    if not args.hamming:
        # vectorize the levenshtein function so we only call once per read
        distance = numpy.vectorize(levenshtein)
    else:
        # vectorize the hamming function so we only call once per read
        distance = numpy.vectorize(hamming)
    read = 0
    for f in glob.glob(os.path.join(args.reads, '*_R1_*')):
        read1 = read_fastq(f)
        # get basename of R1 file
        read1_basename = os.path.basename(f)
        index = read_fastq(os.path.join(
            args.reads,
            read1_basename.replace('R1', 'R2')
        ))
        read2 = read_fastq(os.path.join(
            args.reads,
            read1_basename.replace('R1', 'R3')
        ))
        # read all of our files into fastq iterators
        for r1, i, r2 in izip(read1, index, read2):
            if not args.quiet:
                if read % 100000 == 0:
                    print "{:,}".format(read)
            if not filtered(r1):
                # see if we need to trim the index
                idx, idx_qual = get_index(i, args.tag_length)
                # get index sequence differences from tags
                dist = distance(tags.seqs, idx)
                # get quality values
                good_quality = get_quality(
                    idx_qual,
                    args.min_qual,
                    args.min_mean_qual
                )
                # find tags with 0 distance from other tags
                positions = numpy.where(dist == 0)[0]
                # if not a perfect match, check distance 1 matches
                if (positions.size == 0 and good_quality and
                   idx not in tags.no_correct):
                    positions = numpy.where(dist == 1)[0]
                if positions.size == 1 and good_quality:
                    # assert headers match
                    try:
                        assert r1[0].split(' ')[0] == r2[0].split(' ')[0]
                    except:
                        raise IOError("Mismatch between R1 and R2 headers")
                    # get tag for match
                    match = tags.seqs[positions[0]]
                    # change header to add tag
                    r1 = change_read_num(r1, 1, match)
                    r2 = change_read_num(r2, 2, match)
                    # write to output
                    write_to_r1_and_r2(tags, match, r1, r2)
                # put low quality tags into their own file
                elif ((positions is not None) and (len(positions) == 1) and
                      not good_quality):
                    write_to_unk(tags, 'lowqual', r1, i, r2)
                # put everything else into unknown
                else:
                    write_to_unk(tags, 'unknown', r1, i, r2)
            # if for some reason there are reads not passing filter,
            # put those into unknown, too.
            else:
                write_to_unk(tags, 'unknown', r1, i, r2)
            read += 1
        # close the input read files
        for f in [read1, index, read2]:
            f.close()
    tags.close_zip_files()