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

#from demuxipy import db
from splitaake import pairwise2

from splitaake.lib import FullPaths, ListQueue, Tagged, Parameters

import pdb


def motd():
    """Print a welcome message """
    motd = """
    ###############################################################
    #                       demuxi.py                             #
    #                                                             #
    # Demultiplexing of hierarchically tagged, massively parallel #
    # DNA sequence reads                                          #
    #                                                             #
    #                                                             #
    # Copyright (c) 2009-2011 Brant C. Faircloth                  #
    #                                                             #
    #                                                             #
    # Ecology and Evolutionary Biology                            #
    # 621 Charles E. Young Drive                                  #
    # University of California, Los Angeles, 90095, USA           #
    ###############################################################\n
    """
    print motd


def matches(tag, seq_match_span, tag_match_span, allowed_errors):
    """Determine the gap/error counts for a particular match"""
    # deal with case where tag match might be perfect, but extremely gappy, 
    # e.g. ACGTCGTGCGGA-------------------------ATC
    if tag_match_span.count('-') > allowed_errors or \
        seq_match_span.count('-') > allowed_errors:
        return 0, 0
    else:
        seq_array   = numpy.array(list(seq_match_span))
        tag_array   = numpy.array(list(tag_match_span))
        matches     = sum(seq_array == tag_array)
        error       = sum(seq_array != tag_array) + (len(tag) - \
            len(tag_match_span.replace('-','')))
        return matches, error


def align(seq, tags, allowed_errors):
    """Alignment method for aligning tags with their respective
    sequences.  Only called when regular expression matching patterns fail.
    Inspired by http://github.com/chapmanb/bcbb/tree/master"""
    high_score = {'tag':None, 'seq_match':None, 'mid_match':None, 'score':None, 
        'start':None, 'end':None, 'matches':None, 'errors':allowed_errors}
    for tag in tags:
        #pdb.set_trace()
        try:
            seq_match, tag_match, score, start, end = pairwise2.align.localms(seq, 
            tag, 5.0, -4.0, -9.0, -0.5, one_alignment_only=True)[0]
            seq_match_span  = seq_match[start:end]
            tag_match_span  = tag_match[start:end]
            match, errors   = matches(tag, seq_match_span, tag_match_span, allowed_errors)
            if match >= len(tag)-allowed_errors and match > high_score['matches'] \
                and errors <= high_score['errors']:
                high_score['tag'] = tag
                high_score['seq_match'] = seq_match
                high_score['tag_match'] = tag_match
                high_score['score'] = score
                high_score['start'] = start
                high_score['end'] = end
                high_score['matches'] = match
                high_score['seq_match_span'] = seq_match_span
                high_score['errors'] = errors
        except IndexError:
            pass
    if high_score['matches']:
        return high_score['tag'], high_score['matches'], \
        high_score['seq_match'], high_score['seq_match_span'], \
        high_score['start'], high_score['end']
    else:
        return None


def get_align_match_position(seq_match_span, start, stop):
    # slice faster than ''.startswith()
    if seq_match_span[0] == '-':
        start = start + seq_match_span.count('-')
    else:
        stop = stop - seq_match_span.count('-')
    return start, stop


def find_left_tag(s, tag_regexes, tag_strings, max_gap_char, tag_len, fuzzy, errors):
    """Matching methods for left linker - regex first, followed by fuzzy (SW)
    alignment, if the option is passed"""
    for regex in tag_regexes:
        match = regex.search(s)
        if match is not None:
            m_type = 'regex'
            start, stop = match.start(), match.end()
            # by default, this is true
            tag_matched = regex.pattern.split('}')[1]
            seq_matched = s[start:stop]
            break
    if match is None and fuzzy:
        match = align(s[:max_gap_char + tag_len], tag_strings, errors)
        # we can trim w/o regex
        if match:
            m_type = 'fuzzy'
            tag_matched = match[0]
            seq_matched = match[3]
            start, stop = get_align_match_position(match[3], match[4], match[5])
    if match:
        return tag_matched, m_type, start, stop, seq_matched
    else:
        return None


def find_right_tag(s, tag_regexes, tag_strings, max_gap_char, tag_len,
        fuzzy, errors, tagged, revcomp = True):
    """Matching methods for right linker - regex first, followed by fuzzy (SW)
    alignment, if the option is passed"""
    #if 'MID15_NoError_SimpleX1_NoError_F_NEQ_R' in tagged.read.identifier:
    #    pdb.set_trace()
    for regex in tag_regexes:
        match = regex.search(s)
        if match is not None:
            m_type = 'regex'
            start, stop = match.start(), match.end()
            # by default, this is true
            tag_matched = regex.pattern.split('[')[0]
            seq_matched = s[start:stop]
            break
    if match is None and fuzzy:
        match = align(s[-(tag_len + max_gap_char):], tag_strings, errors)
        # we can trim w/o regex
        if match:
            # correct match_position
            start_of_slice = len(s) - (tag_len + max_gap_char)
            m_type = 'fuzzy'
            tag_matched = match[0]
            seq_matched = match[3]
            start, stop = get_align_match_position(match[3], match[4], match[5])
            start, stop = start + start_of_slice, stop + start_of_slice
    if match and revcomp:
        return DNA_reverse_complement(tag_matched), m_type, start, stop, DNA_reverse_complement(seq_matched)
    elif match and not revcomp:
        return tag_matched, m_type, start, stop, seq_matched
    else:
        return None


def trim_one(tagged, regexes, strings, buff, length, fuzzy, errors, trim = 0):
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
    parser = argparse.ArgumentParser(description = "demuxi.py:  sequence " + \
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
