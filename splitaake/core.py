#!/usr/bin/env python
# encoding: utf-8
"""
File: core.py
Author: Brant Faircloth

Created by Brant Faircloth on 16 June 2012 15:06 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: core functions for demuxipy

"""


#import os
import sys
#import re
import gzip
#import time
import numpy
#import string
#import cPickle
#import sqlite3
import argparse
import itertools
#import ConfigParser

#from multiprocessing import Process, Queue, JoinableQueue, cpu_count

#from seqtools.sequence.fastq import FastqReader
#from seqtools.sequence.fasta import FastaQualityReader
from seqtools.sequence.transform import DNA_reverse_complement
from jellyfish import levenshtein_distance as levenshtein

#from demuxipy import db
from splitaake import pairwise2

from splitaake.lib2 import FullPaths, ListQueue, Tagged, Parameters

import pdb


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
    #       * Reads                                               #
    #       * Hierarchically tagged amplicons                     #
    #       * Combinatorially-tagged amplicons                    #
    #       * 2RAD data                                           #
    #                                                             #
    # from Illumina sequencing.                                   #
    #                                                             #
    #                                                             #
    # Copyright (c) 2009-2012 Brant C. Faircloth                  #
    # Available under 3-clause BSD license                        #
    #                                                             #
    # Ecology and Evolutionary Biology                            #
    # 621 Charles E. Young Drive                                  #
    # University of California, Los Angeles, 90095, USA           #
    ###############################################################\n
    """
    print motd


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


def get_align_match_position(seq_match_span, start, stop):
    # slice faster than ''.startswith()
    if seq_match_span[0] == '-':
        start = start + seq_match_span.count('-')
    else:
        stop = stop - seq_match_span.count('-')
    return start, stop


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


def trim_one(tagged, regexes, strings, buff, length, fuzzy, errors, trim=0):
    """Remove the MID tag from the sequence read"""
    #if sequence.id == 'MID_No_Error_ATACGACGTA':
    #    pdb.set_trace()
    #pdb.set_trace()
    mid = find_left_tag(tagged.read.sequence,
                regexes,
                strings,
                buff,
                length,
                fuzzy,
                errors
            )
    if mid:
        target, match_type, match = mid[0], mid[1], mid[4]
        #tagged.mid, tagged.m_type, tagged.seq_match = mid[0],mid[1],mid[4]
        tagged.read = tagged.read.slice(mid[3] + trim, len(tagged.read.sequence), False)
        #tagged.mid_name = params.sequence_tags.reverse_mid_lookup[tagged.mid]
        return tagged, target, match_type, match
    else:
        return tagged, None, None, None


def trim_two(tagged, fregex, fstring, rregex, rstring, buff,
        length, fuzzy, errors, trim = 0, revcomp = True):
    """Use regular expression and (optionally) fuzzy string matching
    to locate and trim linkers from sequences"""

    left = find_left_tag(tagged.read.sequence,
                fregex,
                fstring,
                buff,
                length,
                fuzzy,
                errors
            )
    
    right = find_right_tag(tagged.read.sequence,
                rregex,
                rstring,
                buff,
                length,
                fuzzy,
                errors,
                tagged,
                revcomp
            )

    # we can have 5 types of matches - tags on left and right sides,
    # tags on left side only, tags on right side only, mismatching tags
    # and no tags at all.  We take care of matching position (i,e. we
    # want only matches at the ends) by regex and slicing in the
    # search methods above.  If for some reason you want to check
    # for concatemers, then turn that function on.

    if left is not None \
            and right is not None \
            and left[0] == right[0]:
        # trim the read
        tagged.read = tagged.read.slice(left[3], right[2], False)
        # left and right are identical so largely pass back the left
        # info... except for m_type which can be a combination
        target, match = left[0], left[4]
        match_type = "{}-{}-both".format(left[1], right[1])

    elif left is not None \
            and right is not None\
            and left[0] != right[0]:
        # these are no good.  check for within gaps
        # to make sure it's not a spurious match
        #pdb.set_trace()
        tagged.trimmed = None
        target, match = None, None
        match_type = "tag-mismatch"

    elif right is None and left and left[2] <= buff:
        tagged.read = tagged.read.slice(left[3], len(tagged.read.sequence), False)
        target, match = left[0], left[4]
        match_type = "{}-left".format(left[1])

    elif left is None and right and right[2] >= (len(tagged.read.sequence) - (len(right[0]) + buff)):
        tagged.read = tagged.read.slice(0, right[2], False)
        target, match = right[0], right[4]
        match_type = "{}-right".format(right[1])

    else:
        target, match_type, match = None, None, None

    return tagged, target, match_type, match


def concat_check(tagged, params):
    """Check screened sequence for the presence of concatemers by scanning 
    for all possible tags - after the 5' and 3' tags have been removed"""
    s = tagged.read.sequence
    m_type = None
    #pdb.set_trace()
    for tag in params.sequence_tags.all_tags[str(tagged.outer_seq)]['regex']:
        match = tag.search(s)
        if match:
            tagged.concat_seq = tag.pattern
            tagged.concat_type = "regex-concat"
            tagged.concat_match = tagged.read.sequence[match.start():match.end()]
            break
    if match is None and params.concat_fuzzy:
        match = align(s,
                params.sequence_tags.all_tags[str(tagged.outer_seq)]['string'], 
                params.concat_allowed_errors
            )
        if match:
            tagged.concat_seq = match[0]
            tagged.concat_type = "fuzzy-concat"
            tagged.concat_match = match[3]
    return tagged


def progress(count, interval, big_interval):
    """give a rudimentary indication of progress"""
    if count % big_interval == 0:
        sys.stdout.write('%')
        sys.stdout.flush()
    elif count % interval == 0:
        sys.stdout.write('.')
        sys.stdout.flush()


def get_args():
    """get arguments (config file location)"""
    parser = argparse.ArgumentParser(description = "splitaake.py:  sequence " + \
        "demultiplexing for hierarchically-tagged samples")
    parser.add_argument('config', help="The input configuration file",
            action=FullPaths)
    return parser.parse_args()


def get_sequence_count(input, kind):
    """Determine the number of sequence reads in the input"""
    if kind == 'fasta':
        return sum([1 for line in open(input, 'rU') if line.startswith('>')])
    elif kind == 'fastq' and input.endswith('gz'):
        return sum([1 for line in gzip.open(input, 'rb')]) / 4
    else:
        return sum([1 for line in open(input, 'rU')]) / 4


def split_fasta_reads_into_groups(reads, num_reads, num_procs):
    job_size = num_reads / num_procs
    print "Parsing reads into groups of {} reads".format(job_size)
    i = iter(reads)
    chunk = list(itertools.islice(i, job_size))
    while chunk:
        yield chunk
        chunk = list(itertools.islice(i, job_size))


def merge_fastq(a, b):
    for i, j in itertools.izip(a, b):
        yield i, j
