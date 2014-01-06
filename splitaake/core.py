#!/usr/bin/env python
# encoding: utf-8
"""
File: core.py
Author: Brant Faircloth

Created by Brant Faircloth on 16 June 2012 15:06 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: core functions for demuxipy

"""


import os
import sys
import gzip
import numpy
import errno
import string
import itertools

from jellyfish import levenshtein_distance as levenshtein

from splitaake import pairwise2

import pdb


class AlignScore:
    def __init__(self, allowed_errors):
        self.tag = None
        self.seq = None
        self.score = 0
        self.start = None
        self.end = None
        self.offset = 0
        self.matches = 0
        self.errors = None
        self.allowed_errors = allowed_errors

    def set(self, tag, seq_match, seq_span, score, start, end, match, errors, offset):
        self.tag = tag
        self.seq_match = seq_match
        self.seq_span = seq_span
        self.score = score
        self.start = start
        self.end = end
        self.matches = match
        self.errors = errors
        self.offset = offset


class Demux:
    '''Trimming, tag, and sequence data for individual reads'''
    def __init__(self, identifier):
        # super(Params, self).__init__()
        #assert isinstance(sequence, FastaSequence), \
        #    'The Record class must be instantiated with a FastaSequence object'
        # a biopython sequence object
        self.name = identifier
        self.r1 = None
        self.r2 = None
        self.individual = None
        self.drop = False
        self.r1site = SeqSearchResult('site')
        self.r2site = SeqSearchResult('site')
        self.r1tag = SeqSearchResult('tag')
        self.r2tag = SeqSearchResult('tag')
        self.r1overrun = SeqSearchResult('overrun')
        self.r2overrun = SeqSearchResult('overrun')

    def __str__(self):
        return "{0}({1})".format(self.__class__, self.name)

    def __repr__(self):
        return "<{0} instance at {1}>".format(self.__class__, hex(id(self)))


class SeqSearchResult:
    def __init__(self, typ):
        self.type = typ
        self.name = None
        self.seq = None
        self.tag = None
        self.start = None
        self.end = None
        self.match = False
        self.match_type = None

    def __str__(self):
        return "{0}({1})".format(self.__class__, self.__dict__['name'])

    def __repr__(self):
        return "<{0} instance at {1}>".format(self.__class__, hex(id(self)))

    def reset(self):
        self.name = None
        self.seq = None
        self.tag = None
        self.start = None
        self.end = None
        self.match = False
        self.match_type = None


def motd():
    """Print a welcome message """
    motd = """
    ###############################################################
    #           _____                          _                  #
    #          /  ___|     | (_) |            | |                 #
    #          \ `--. _ __ | |_| |_ __ _  __ _| | _____           #
    #           `--. \ '_ \| | | __/ _` |/ _` | |/ / _ \          #
    #          /\__/ / |_) | | | || (_| | (_| |   <  __/          #
    #          \____/| .__/|_|_|\__\__,_|\__,_|_|\_\___|          #
    #                | |                                          #
    #                |_|                                          #
    #                                                             #
    # Demultiplexing of:                                          #
    #                                                             #
    #     * Illumina Reads                                        #
    #     * Illumina sequenced, combinatorially-tagged amplicons  #
    #     * Illumina sequence 2RAD data                           #
    #                                                             #
    #                                                             #
    # Copyright (c) 2009-2013 Brant C. Faircloth                  #
    # Available under 3-clause BSD license                        #
    #                                                             #
    # Ecology and Evolutionary Biology                            #
    # 621 Charles E. Young Drive                                  #
    # University of California, Los Angeles, 90095, USA           #
    ###############################################################\n
    """
    print motd


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def complement(s):
    complement = string.maketrans(
        'acgtrymkbdhvACGTRYMKBDHV',
        'tgcayrkmvhdbTGCAYRKMVHDB'
    )
    return s.translate(complement)


def reverse_complement(s):
    return complement(s)[::-1]


def merge_fastq(a, b):
    for i, j in itertools.izip(a, b):
        yield i, j


def split_and_check_reads(pair):
    assert len(pair) == 2
    r1, r2 = pair
    identifiers = [i.identifier.split(' ')[0] for i in [r1, r2]]
    # make sure fastq headers are the same
    assert identifiers[0] == identifiers[1]
    return r1, r2, identifiers[0]


def split_every(n, iterable):
    i = iter(iterable)
    piece = list(itertools.islice(i, n))
    while piece:
        yield piece
        piece = list(itertools.islice(i, n))


def matches(tag, seq_span, tag_span, allowed_errors):
    """Determine the gap/error counts for a particular match"""
    # deal with case where tag match might be perfect, but extremely gappy,
    # e.g. ACGTCGTGCGGA-------------------------ATC
    if tag_span.count('-') > allowed_errors or seq_span.count('-') > allowed_errors:
        return 0, 100
    else:
        seq_array = numpy.array(list(seq_span))
        tag_array = numpy.array(list(tag_span))
        matches = sum(seq_array == tag_array)
        error = sum(seq_array != tag_array) + (tag.tag_len - len(tag_span.replace('-', '')))
        return matches, error


def align(seq, tags, allowed_errors):
    """Alignment method for aligning tags with their respective
    sequences.  Only called when regular expression matching patterns fail.
    Inspired by http://github.com/chapmanb/bcbb/tree/master"""
    score = AlignScore(allowed_errors)
    for tag in tags:
        try:
            seq_match, tag_match, aln_score, start, end = pairwise2.align.localms(
                    seq,
                    tag.string,
                    5.0,
                    -4.0,
                    -9.0,
                    -0.5,
                    one_alignment_only=True
                )[0]
            # the offset used below is a mechanism to deal with trimmng the trailing
            # lowercase base when the tag sequence is ATCGac or similar.  We want to
            # compare the sequence (in upper()) to the tag, but we want to trim the
            # sequence + the spacer from the resulting read.  This is not a problem
            # for tags like acATCG, since we're trimming from the end of the tag forward.
            if tag.string[0].islower():
                offset = 0
                try:
                    start = max([k for k, v in enumerate(tag_match) if v.islower()]) + 1
                except ValueError, e:
                    if e.message == 'max() arg is an empty sequence':
                        start = 0
                    else:
                        raise ValueError(e)
                end = start + tag.tag_len
            elif tag.string[-1].islower():
                offset = 0
                try:
                    positions = [k for k, v in enumerate(tag_match) if v.islower()]
                    end = min(positions)
                    offset = (max(positions) - end) + 1
                except ValueError, e:
                    if e.message == 'max() arg is an empty sequence':
                        start = 0
                    else:
                        raise ValueError(e)
                start = end - tag.tag_len
            else:
                offset = 0
            if start is None:
                start = 0
            seq_span = seq_match[start:end]
            tag_span = tag_match[start:end]
            match, errors = matches(tag, seq_span, tag_span, allowed_errors)
            if start > 0:
                for pos, base in enumerate(tag_match[:start]):
                    if seq_match[pos] == base or seq_match[pos] == '-':
                        start_ok = True
                    else:
                        start_ok = False
            else:
                start_ok = True
            #print "{}\t\t{}, matches = {}, errors = {}, good = {}".format(seq_match, seq_span, match, errors, start_ok)
            #print "{}\t\t{}".format(tag_match, tag_span)
            if start_ok and match >= (tag.tag_len - allowed_errors) and (match > score.matches) and (errors <= score.allowed_errors):
                score.set(tag, seq_match, seq_span, score, start, end, match, errors, offset)
        except IndexError:
            pass
    if score.matches is not 0:
        return score
    else:
        return None


def trim_tags_from_reads(read, start, end):
    return read.slice(end, len(read.sequence), False)


def get_align_match_position(seq_match_span, start, stop):
    # slice faster than ''.startswith()
    if seq_match_span[0] == '-':
        start = start + seq_match_span.count('-')
    else:
        stop = stop - seq_match_span.count('-')
    return start, stop


def find_which_reads_have_tags(params, r1, r2, dmux, p=False):
    pair = [r1, r2]
    for index, read in enumerate(pair):
        dmux.r1tag = find_tag(read.sequence, params.tags, 'r1', 'regex', dmux.r1tag)
        if dmux.r1tag.match:
            dmux.r1 = read
            #dmux.r1tag = forward
            trim_tags_from_reads(dmux.r1, dmux.r1tag.start, dmux.r1tag.end + dmux.r1tag.offset)
            if p:
                print dmux.r1.identifier
                print "regex forward: tag={0} seq={1} type={2} sequence={3}".format(dmux.r1tag.tag, dmux.r1tag.seq, dmux.r1tag.match_type, dmux.r1.sequence[:20])
            #pdb.set_trace()
            # remove the read from further consideration
            pair.pop(index)
            break
    for index, read in enumerate(pair):
        dmux.r2tag = find_tag(read.sequence, params.tags, 'r2', 'regex', dmux.r2tag)
        if dmux.r2tag.match:
            dmux.r2 = read
            #dmux.r2tag = reverse
            trim_tags_from_reads(dmux.r2, dmux.r2tag.start, dmux.r2tag.end + dmux.r2tag.offset)
            if p:
                print dmux.r2.identifier
                print "regex reverse: tag={0} seq={1} type={2} sequence={3}".format(dmux.r2tag.tag, dmux.r2tag.seq, dmux.r2tag.match_type, dmux.r2.sequence[:20])
            pair.pop(index)
            break
    if not dmux.r1tag.match:
        for index, read in enumerate(pair):
            dmux.r1tag = find_tag(read.sequence, params.tags, 'r1', 'fuzzy', dmux.r1tag)
            if dmux.r1tag.match:
                dmux.r1 = read
                #dmux.r1tag = forward
                trim_tags_from_reads(dmux.r1, dmux.r1tag.start, dmux.r1tag.end + dmux.r1tag.offset)
                if p:
                    print dmux.r1.identifier
                    print "fuzzy forward: tag={0} seq={1} type={2} sequence={3}".format(dmux.r1tag.tag, dmux.r1tag.seq, dmux.r1tag.match_type, dmux.r1.sequence[:20])
                #pdb.set_trace()
                # remove the read from further consideration
                pair.pop(index)
                break
    if not dmux.r2tag.match:
        for index, read in enumerate(pair):
            dmux.r2tag = find_tag(read.sequence, params.tags, 'r2', 'fuzzy', dmux.r2tag)
            if dmux.r2tag.match:
                dmux.r2 = read
                #dmux.r2tag = reverse
                trim_tags_from_reads(dmux.r2, dmux.r2tag.start, dmux.r2tag.end + dmux.r2tag.offset)
                if p:
                    print dmux.r2.identifier
                    print "fuzzy reverse: tag={0} seq={1} type={2} sequence={3}".format(dmux.r2tag.tag, dmux.r2tag.seq, dmux.r2tag.match_type, dmux.r2.sequence[:20])
                break
    if p:
        pdb.set_trace()
    return dmux


def find_tag(seq, tags, read, match_type, result):
    """Matching methods for left linker - regex first, followed by fuzzy (SW)
    alignment, if the option is passed"""
    #result = SeqSearchResult('tag')
    #pdb.set_trace()
    if read == 'r1':
        tag_group = tags.r1
    elif read == 'r2':
        tag_group = tags.r2
    if match_type == 'regex':
        for tag in tag_group:
            match = tag.regex.search(seq)
            if match is not None:
                result.match = True
                # by default, this is true
                assert match.groups()[0] == tag.tag
                result.match_type = 'regex'
                result.start, result.end = match.start(), match.end()
                result.tag = tag.string
                result.seq = seq[result.start:result.end]
                result.name = tags.name_lookup(result.tag)
                result.offset = 0
                break
    elif match_type == 'fuzzy':
        match = align(seq[:tags.max_tag_length], tag_group, tags.errors)
        if match:
            result.match = True
            result.match_type = 'fuzzy'
            result.tag = match.tag.tag
            result.seq = match.seq_span
            result.name = tags.name_lookup(match.tag.string)
            result.start, result.end = get_align_match_position(match.seq_span, match.start, match.end)
            result.offset = match.offset
    # check ALL resulting values for correct levenshtein distance or
    # reset match parameters to None/False
    if result.match:
        if not levenshtein(result.seq.upper(), result.tag.upper()) <= tags.errors:
            result.reset()
    return result


def find_which_reads_have_sites(params, dmux):
    find_site(dmux.r1.sequence, params.site, 'r1', dmux.r1site)
    if dmux.r1site.match:
        trim_tags_from_reads(dmux.r1, dmux.r1site.start, dmux.r1site.end)
    find_site(dmux.r2.sequence, params.site, 'r2', dmux.r2site)
    if dmux.r2site.match:
        trim_tags_from_reads(dmux.r2, dmux.r2site.start, dmux.r2site.end)
    return dmux


def find_site(seq, sites, read, result):
    if read == 'r1':
        site_group = sites.r1
    elif read == 'r2':
        site_group = sites.r2
    for site in site_group:
        match = site.regex.search(seq)
        if match is not None:
            result.match = True
            # by default, this is true
            assert match.groups()[0] == site.string
            result.match_type = 'regex'
            result.start, result.end = match.start(), match.end()
            result.tag = site.string
            result.seq = seq[result.start:result.end]
            break
    if match is None and sites.fuzzy:
        match = align(seq[:sites.max_site_length], site_group, sites.errors)
        if match:
            result.match = True
            result.match_type = 'fuzzy'
            result.tag = match.tag.string
            result.seq = match.seq_span
            result.start, result.end = get_align_match_position(match.seq_span, match.start, match.end)
    return result


def find_which_reads_have_overruns(params, dmux):
    p = params.overruns["{},{}".format(
            dmux.r1tag.name,
            dmux.r2tag.name
            )
        ]
    find_overrun(dmux.r1.sequence, p, 'r1', dmux.r1overrun)
    find_overrun(dmux.r2.sequence, p, 'r2', dmux.r2overrun)
    # if we overrun on one read, we should always overrun on the opposite read,
    # so only trim when this is True for both
    if dmux.r1overrun.match == True and dmux.r2overrun.match == True:
        dmux.r1.slice(0, dmux.r1overrun.start, False)
        dmux.r2.slice(0, dmux.r2overrun.start, False)
    return dmux

def find_overrun(seq, sites, read, result):
    if read == 'r1':
        site_group = sites['r1']
    elif read == 'r2':
        site_group = sites['r2']
    for site in site_group:
        match = site.regex.search(seq)
        if match is not None:
            result.match = True
            # by default, this is true
            assert match.groups()[0] == site.string
            result.match_type = 'regex'
            result.start, result.end = match.start(), match.end()
            result.tag = site.string
            result.seq = seq[result.start:result.end]
            break
    '''
    if match is None and sites.fuzzy:
        match = align(seq[:sites.max_site_length], site_group, sites.errors)
        if match:
            result.match = True
            result.match_type = 'fuzzy'
            result.tag = match.tag.string
            result.seq = match.seq_span
            result.start, result.end = get_align_match_position(match.seq_span, match.start, match.end)
    '''
    return result


def get_sequence_count(input, kind):
    """Determine the number of sequence reads in the input"""
    if kind == 'fasta':
        return sum([1 for line in open(input, 'rU') if line.startswith('>')])
    elif kind == 'fastq' and input.endswith('gz'):
        return sum([1 for line in gzip.open(input, 'rb')]) / 4
    else:
        return sum([1 for line in open(input, 'rU')]) / 4

'''
def split_fasta_reads_into_groups(reads, num_reads, num_procs):
    job_size = num_reads / num_procs
    print "Parsing reads into groups of {} reads".format(job_size)
    i = iter(reads)
    chunk = list(itertools.islice(i, job_size))
    while chunk:
        yield chunk
        chunk = list(itertools.islice(i, job_size))
'''


def progress(count, interval, big_interval):
    """give a rudimentary indication of progress"""
    if count % big_interval == 0:
        sys.stdout.write('%')
        sys.stdout.flush()
    elif count % interval == 0:
        sys.stdout.write('.')
        sys.stdout.flush()
